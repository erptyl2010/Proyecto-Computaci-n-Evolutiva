import numpy as np
from pymoo.indicators.hv import HV
from pymoo.util.nds.non_dominated_sorting import NonDominatedSorting
from ce_peptidos.datos import concentration, is_micromolar, similarity_groups, peptide_reason
from ce_peptidos.modelos import features
from ce_peptidos.evolucion import (BASE, domain, nondominated, hypervolume, TwoChangesRepair,
                                  CategoricalMutation, UniformCrossover, run_search)

def test_concentration_censoring_and_units():
    for raw in ["<10", ">200", "10-20", "10 ± 2", "0", "-1", "nan", None, "inf"]:
        assert concentration(raw) is None
    assert concentration(" 6.25 ") == 6.25
    assert concentration("1e2") == 100
    assert is_micromolar("µM") and is_micromolar("μM")
    assert not is_micromolar("mg/L")

def test_chemical_annotation_and_features():
    p = {"sequence": BASE, "complexity": {"name": "Monomer"}}
    assert peptide_reason(p) == ""
    assert peptide_reason(p | {"cTerminus": {"name": "AMD"}}) == ""
    assert peptide_reason(p | {"cTerminus": {"name": "OTHER"}}) == "terminal_no_representado"
    assert peptide_reason(p | {"sequence": BASE.lower()}) == "secuencia_no_canonica_o_longitud"
    missing = features(BASE)
    amide = features(BASE, ct="AMD")
    assert len(missing) == 27 and len(features(BASE, ordered=True)) == 42
    assert missing["c_annotation_missing"] == 1 and amide["c_annotation_missing"] == 0
    assert missing != amide

def test_transitive_groups():
    g = similarity_groups(["AAAA", "AACC", "CCCC", "WWWW", "AAAAA"])
    assert g[0] == g[1] == g[2]
    assert len(set(g)) == 3

def test_domain_and_categorical_variation():
    X = domain()
    assert X.shape == (76210, 21)
    assert len(np.unique(X, axis=0)) == 76210
    assert np.max(np.sum(X != X[0], axis=1)) == 2
    np.random.seed(123)
    parents = np.array([X[100:200], X[1000:1100]])
    crossed = UniformCrossover(.9)._do(None, parents).reshape(-1, 21)
    mutated = CategoricalMutation(1.)._do(None, crossed)
    assert np.all(mutated != crossed)
    repaired = TwoChangesRepair()._do(None, mutated)
    assert np.all((repaired >= 0) & (repaired < 20))
    assert np.max(np.sum(repaired != X[0], axis=1)) <= 2

def test_pareto_ties_and_hv_independent_reference():
    F = np.array([[1, 4], [2, 2], [2, 2], [3, 3], [4, 1], [2, 4]])
    assert nondominated(F).tolist() == [True, True, True, False, True, False]
    assert hypervolume(F, np.array([5, 5])) == 11.
    rng = np.random.default_rng(14)
    F = rng.integers(0, 20, (300, 2)).astype(float)
    indices = NonDominatedSorting().do(F, only_non_dominated_front=True)
    assert set(np.flatnonzero(nondominated(F))) == set(indices)
    assert np.isclose(hypervolume(F, np.array([21, 21])), HV(ref_point=np.array([21, 21]))(F))

class ToyLandscape:
    """Prueba de software: no es un resultado del estudio de péptidos."""
    def __init__(self):
        from ce_peptidos.evolucion import decode
        self.base = "AAA"
        self.X = domain(self.base)
        self.lookup = {tuple(x): i for i, x in enumerate(self.X)}
        self.sequences = [decode(x) for x in self.X]
        self.F = self.predict(np.arange(len(self.X)))
        self.reference = self.F.max(axis=0) + 1
        self.hv_star = hypervolume(self.F, self.reference)

    def predict(self, indices):
        X = self.X[indices]
        return np.column_stack((X.sum(axis=1) / 20, ((X - 10)**2).sum(axis=1) / 100))

def test_exact_unique_budget_pairing_and_repeatability():
    landscape = ToyLandscape()
    a, ha, ea = run_search(landscape, "NSGA-II", 42, population=20, mutation=.3, budget=113)
    b, hb, eb = run_search(landscape, "NSGA-II", 42, population=20, mutation=.3, budget=113)
    _, _, er = run_search(landscape, "Aleatoria", 42, population=20, budget=113)
    assert a["unique_evaluations"] == b["unique_evaluations"] == 113
    assert ea.sequence.nunique() == 113
    assert ea.equals(eb)
    assert ha.drop(columns="elapsed_seconds").equals(hb.drop(columns="elapsed_seconds"))
    assert ea.sequence.iloc[:20].tolist() == er.sequence.iloc[:20].tolist()
