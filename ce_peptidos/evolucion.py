"""NSGA-II categórico y referencia exhaustiva del paisaje PREDICHO."""
from itertools import combinations, product
from time import perf_counter
import numpy as np
import pandas as pd
from pymoo.algorithms.moo.nsga2 import NSGA2
from pymoo.core.problem import Problem
from pymoo.core.crossover import Crossover
from pymoo.core.mutation import Mutation
from pymoo.core.repair import Repair
from pymoo.core.duplicate import DefaultDuplicateElimination
from .datos import AA
from .modelos import matrix

BASE = "GFMKYIKPLIPHAVKAIKKLI"

def encode(sequence):
    return np.array([AA.index(a) for a in sequence], dtype=np.int16)

def decode(x):
    return "".join(AA[int(a)] for a in x)

def domain(base=BASE):
    b = encode(base)
    rows = [b.copy()]
    for k in (1, 2):
        for positions in combinations(range(len(b)), k):
            options = [[a for a in range(20) if a != b[p]] for p in positions]
            for alleles in product(*options):
                x = b.copy()
                x[list(positions)] = alleles
                rows.append(x)
    return np.array(rows)

def nondominated(F):
    """Minimización 2D, conservando TODAS las filas empatadas no dominadas."""
    F = np.asarray(F, dtype=float)
    keep = np.zeros(len(F), dtype=bool)
    order = np.lexsort((F[:, 1], F[:, 0]))
    best_y, best_x = np.inf, np.inf
    for i in order:
        x, y = F[i]
        if y < best_y:
            keep[i], best_x, best_y = True, x, y
        elif y == best_y and x == best_x:
            keep[i] = True
    return keep

def hypervolume(F, reference):
    F = np.asarray(F, dtype=float).reshape(-1, 2)
    F = F[np.all(F < reference, axis=1)]
    if not len(F):
        return 0.
    front = np.unique(F[nondominated(F)], axis=0)
    previous = float(reference[1])
    area = 0.
    for x, y in front:
        area += (reference[0] - x) * (previous - y)
        previous = y
    return float(area)

class UniformCrossover(Crossover):
    def __init__(self, probability):
        super().__init__(2, 2, prob=probability)

    def _do(self, problem, X, **kwargs):
        mask = np.random.random(X.shape[1:]) < .5
        return np.array([np.where(mask, X[0], X[1]), np.where(mask, X[1], X[0])])

class CategoricalMutation(Mutation):
    def __init__(self, per_gene):
        super().__init__(prob=1.)
        self.per_gene = per_gene

    def _do(self, problem, X, **kwargs):
        mask = np.random.random(X.shape) < self.per_gene
        # Sumar entre 1 y 19 módulo 20 escoge uniformemente un aminoácido distinto.
        alternatives = (X + np.random.randint(1, 20, X.shape)) % 20
        return np.where(mask, alternatives, X).astype(np.int16)

class TwoChangesRepair(Repair):
    def __init__(self, base=BASE):
        super().__init__()
        self.base = encode(base)

    def _do(self, problem, X, **kwargs):
        X = X.astype(np.int16, copy=True)
        for row in X:
            changed = np.flatnonzero(row != self.base)
            if len(changed) > 2:
                restore = np.random.choice(changed, len(changed) - 2, replace=False)
                row[restore] = self.base[restore]
        return X

class GlobalDuplicates(DefaultDuplicateElimination):
    def __init__(self, problem):
        super().__init__()
        self.problem = problem

    def _do(self, pop, other, is_duplicate):
        is_duplicate = super()._do(pop, other, is_duplicate)
        for i, row in enumerate(pop.get("X")):
            if tuple(row) in self.problem.seen:
                is_duplicate[i] = True
        return is_duplicate

class Landscape:
    def __init__(self, models, base=BASE):
        start = perf_counter()
        self.base = base
        self.X = domain(base)
        self.sequences = [decode(x) for x in self.X]
        self.lookup = {tuple(x): i for i, x in enumerate(self.X)}
        self.models = models
        frame = pd.DataFrame({"sequence": self.sequences, "n_terminus": "SIN_ANOTACION",
                              "c_terminus": "SIN_ANOTACION"})
        self.features = {key: matrix(frame, model["ordered"]).to_numpy()
                         for key, model in models.items()}
        self.feature_seconds = perf_counter() - start
        start = perf_counter()
        self.F = self.predict(np.arange(len(self.X)))
        self.exhaustive_seconds = perf_counter() - start
        self.reference = self.F.max(axis=0) + .15
        self.front_mask = nondominated(self.F)
        self.hv_star = hypervolume(self.F, self.reference)

    def predict(self, indices):
        mic = self.models["MIC"]["model"].predict(self.features["MIC"][indices])
        hc = self.models["HC50"]["model"].predict(self.features["HC50"][indices])
        return np.column_stack((mic, -hc))

class PeptideProblem(Problem):
    def __init__(self, landscape):
        super().__init__(n_var=len(landscape.base), n_obj=2, xl=0, xu=19, vtype=int)
        self.landscape = landscape
        self.seen, self.indices, self.values = set(), [], []

    def _evaluate(self, X, out, *args, **kwargs):
        keys = [tuple(x) for x in X]
        if len(set(keys)) != len(keys) or any(key in self.seen for key in keys):
            raise RuntimeError("Evaluación duplicada: el presupuesto debe contar secuencias únicas.")
        indices = [self.landscape.lookup[key] for key in keys]
        # Los algoritmos consultan el modelo; nunca reciben los valores del barrido exhaustivo.
        F = self.landscape.predict(indices)
        self.seen.update(keys)
        self.indices.extend(indices)
        self.values.extend(F.tolist())
        out["F"] = F

def checkpoint(problem, current, generation, start):
    F = np.array(problem.values)
    return {"generation_or_batch": generation, "n_evaluations": len(F),
            "mean_log10_MIC": float(current[:, 0].mean()),
            "best_log10_MIC": float(current[:, 0].min()),
            "mean_log10_HC50": float(-current[:, 1].mean()),
            "best_log10_HC50": float(-current[:, 1].min()),
            "mean_MIC_uM": float(np.power(10., current[:, 0]).mean()),
            "best_MIC_uM": float(np.power(10., current[:, 0].min())),
            "mean_HC50_uM": float(np.power(10., -current[:, 1]).mean()),
            "best_HC50_uM": float(np.power(10., -current[:, 1].min())),
            "archive_best_log10_MIC": float(F[:, 0].min()),
            "archive_best_log10_HC50": float(-F[:, 1].min()),
            "hv_percent": 100 * hypervolume(F, problem.landscape.reference) / problem.landscape.hv_star,
            "elapsed_seconds": perf_counter() - start}

def run_search(landscape, method, seed, population=60, mutation=1/21, crossover=.9, budget=2100):
    if budget > len(landscape.X) or population > budget:
        raise ValueError("Presupuesto/población fuera del dominio.")
    problem = PeptideProblem(landscape)
    random_ids = np.random.default_rng(seed).choice(len(landscape.X), budget, replace=False)
    initial = landscape.X[random_ids[:population]].copy()
    start = perf_counter()
    history = []
    if method == "NSGA-II":
        algorithm = NSGA2(pop_size=population, sampling=initial,
                          crossover=UniformCrossover(crossover), mutation=CategoricalMutation(mutation),
                          repair=TwoChangesRepair(landscape.base), eliminate_duplicates=GlobalDuplicates(problem))
        algorithm.setup(problem, termination=("n_eval", budget), seed=seed, verbose=False)
        while algorithm.has_next():
            algorithm.n_offsprings = min(population, budget - len(problem.indices))
            previous = len(problem.indices)
            algorithm.next()
            if len(problem.indices) == previous:
                raise RuntimeError("La búsqueda no generó nuevas secuencias antes de agotar el presupuesto.")
            history.append(checkpoint(problem, algorithm.pop.get("F"), len(history)+1, start))
        current_meaning = "Poblacion superviviente; ambos extremos pueden corresponder a individuos diferentes."
    elif method == "Aleatoria":
        for offset in range(0, budget, population):
            out = {}
            problem._evaluate(landscape.X[random_ids[offset:offset+population]], out)
            history.append(checkpoint(problem, out["F"], len(history)+1, start))
        current_meaning = "Lote aleatorio recien evaluado; no es una poblacion evolutiva."
    else:
        raise ValueError(method)
    seconds = perf_counter() - start
    if len(problem.indices) != budget or len(problem.seen) != budget:
        raise RuntimeError("Presupuesto incompleto o duplicado.")
    F = np.array(problem.values)
    if not np.allclose(F, landscape.F[problem.indices], rtol=0, atol=1e-12):
        raise RuntimeError("La referencia exhaustiva y el predictor consultado discrepan.")
    front = nondominated(F)
    summary = {"method": method, "seed": seed, "population": population, "mutation": mutation,
               "crossover": crossover, "budget": budget, "unique_evaluations": len(problem.seen),
               "hv_percent": history[-1]["hv_percent"], "seconds": seconds,
               "n_front_sequences": int(front.sum()), "n_front_objective_pairs": len(np.unique(F[front], axis=0)),
               "current_fitness_meaning": current_meaning}
    archive = pd.DataFrame({"evaluation": np.arange(1, len(F)+1), "domain_index": problem.indices,
                            "sequence": [landscape.sequences[i] for i in problem.indices],
                            "pred_log10_MIC": F[:, 0], "pred_log10_HC50": -F[:, 1],
                            "pred_MIC_uM": 10**F[:, 0], "pred_HC50_uM": 10**(-F[:, 1]),
                            "nondominated_in_run": front})
    return summary, pd.DataFrame(history), archive
