"""Auditoría de resultados reales; no vuelve a consultar la web."""
from pathlib import Path
import argparse
import hashlib
import json
import sys
import re
import unicodedata
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import joblib
import numpy as np
import pandas as pd
from sklearn.metrics import mean_absolute_error, r2_score
from pymoo.indicators.hv import HV
from ce_peptidos.modelos import matrix
from ce_peptidos.datos import similarity_groups
from ce_peptidos.evolucion import BASE, Landscape, run_search, nondominated

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--repeat-one", action="store_true", help="Repite de verdad una corrida final de 2 100 secuencias.")
    args = parser.parse_args()
    report = {}
    manifest = json.loads((ROOT / "data/raw/manifest.json").read_text(encoding="utf-8"))
    assert manifest["status"] == "complete" and not manifest["failed_records"]
    for rel, h in manifest["sha256"].items():
        assert hashlib.sha256((ROOT / "data/raw" / rel).read_bytes()).hexdigest() == h, rel
    report["raw_files_verified"] = len(manifest["sha256"])
    partitions = pd.read_csv(ROOT / "results/models/particiones.csv")
    metrics = pd.read_csv(ROOT / "results/models/validacion_reservada.csv")
    models = joblib.load(ROOT / "results/models/modelos_congelados.joblib")
    references = json.loads((ROOT / "data/processed/referencias_dbaasp.json").read_text(encoding="utf-8"))
    for objective in ["MIC", "HC50"]:
        df = pd.read_csv(ROOT / f"data/processed/dataset_{objective}.csv", keep_default_na=False)
        part = partitions[partitions.objective == objective]
        groups = similarity_groups(df.sequence.tolist())
        merged = df.merge(part, on="chemical_id", validate="one_to_one")
        assert merged.sequence.tolist() == df.sequence.tolist()
        assert (groups == merged.group.to_numpy()).all()
        assert not set(merged[merged.split == "train"].group) & set(merged[merged.split == "holdout"].group)
        test = merged[merged.split == "holdout"]
        pred = models[objective]["model"].predict(matrix(test, models[objective]["ordered"]).to_numpy())
        y = np.log10(test.value_uM)
        metric = metrics[metrics.objective == objective].iloc[0]
        assert np.isclose(mean_absolute_error(y, pred), metric.mae_log10)
        assert np.isclose(r2_score(y, pred), metric.r2_log10)
        pubs = {}
        for split, subgroup in merged.groupby("split"):
            pubs[split] = set()
            for source_ids in subgroup.article_ids_all:
                for article_id in json.loads(source_ids):
                    article = references.get(str(article_id), {})
                    title = re.sub(r"[^\w]+", "", unicodedata.normalize("NFKC", article.get("title", "")).casefold())
                    if title:
                        pubs[split].add((title, str(article.get("year", ""))))
        assert len(pubs['train'] & pubs['holdout']) == metric.shared_publications_train_holdout
    report["reserved_partitions_and_metrics"] = "verified"
    opt = ROOT / "results/optimization"
    reference = json.loads((opt / "referencia.json").read_text(encoding="utf-8"))
    oracle = pd.read_csv(opt / "dominio_exhaustivo.csv")
    assert len(oracle) == oracle.sequence.nunique() == 76210
    assert oracle.sequence.map(lambda s: len(s) == 21 and sum(a != b for a, b in zip(s, BASE)) <= 2).all()
    F = oracle[["pred_log10_MIC", "pred_log10_HC50"]].to_numpy() * [1, -1]
    assert (nondominated(F) == oracle.nondominated).all()
    hv = HV(ref_point=np.array(reference["reference_point"]))
    assert np.isclose(hv(F[oracle.nondominated]), reference["hv_exhaustive"])
    report["exhaustive_domain_and_hv"] = "verified"
    comparison = pd.read_csv(opt / "comparacion_final.csv")
    sensitivity = pd.read_csv(opt / "sensibilidad.csv")
    all_runs = pd.concat([comparison, sensitivity], ignore_index=True).drop_duplicates("run_id")
    for row in all_runs.itertuples():
        folder = opt / "runs" / row.run_id
        details = json.loads((folder / "complete.json").read_text(encoding="utf-8"))
        assert details['spec']['models_sha256'] == hashlib.sha256((ROOT / "results/models/modelos_congelados.joblib").read_bytes()).hexdigest()
        assert details['spec']['source_sha256'] == hashlib.sha256((ROOT / "ce_peptidos/evolucion.py").read_bytes()).hexdigest()
        for filename, h in details["sha256"].items():
            assert hashlib.sha256((folder / filename).read_bytes()).hexdigest() == h
        evaluations = pd.read_csv(folder / "evaluaciones.csv")
        history = pd.read_csv(folder / "historia.csv")
        assert len(evaluations) == evaluations.sequence.nunique() == 2100
        f = evaluations[["pred_log10_MIC", "pred_log10_HC50"]].to_numpy() * [1, -1]
        assert np.allclose(f, F[evaluations.domain_index], rtol=0, atol=1e-12)
        assert np.isclose(hv(f[nondominated(f)]) / reference["hv_exhaustive"] * 100, row.hv_percent)
        assert np.isclose(history.iloc[-1].hv_percent, row.hv_percent)
        assert (np.diff(history.hv_percent) >= -1e-10).all()
    for seed, g in comparison.groupby("seed"):
        first = [pd.read_csv(opt / "runs" / r.run_id / "evaluaciones.csv").sequence.iloc[:int(r.population)].tolist() for r in g.itertuples()]
        assert len(first) == 2 and first[0] == first[1]
    report["unique_runs_verified"] = len(all_runs)
    report["final_pairs_initial_population"] = "identical_in_each_pair"
    candidates = pd.read_csv(opt / "candidatos_encontrados.csv")
    found = set()
    for row in comparison[comparison.method == "NSGA-II"].itertuples():
        found.update(pd.read_csv(opt / "runs" / row.run_id / "evaluaciones.csv").sequence)
    assert set(candidates.sequence) <= found
    report["candidate_provenance"] = "observed_in_final_nsga_runs"
    if args.repeat_one:
        row = comparison[comparison.method == "NSGA-II"].iloc[0]
        landscape = Landscape(models)
        summary, history, archive = run_search(landscape, "NSGA-II", int(row.seed), population=int(row.population),
                                               mutation=float(row.mutation), crossover=float(row.crossover))
        saved = pd.read_csv(opt / "runs" / row.run_id / "evaluaciones.csv")
        assert archive.sequence.tolist() == saved.sequence.tolist()
        assert np.allclose(archive.pred_log10_MIC, saved.pred_log10_MIC, rtol=0, atol=1e-12)
        assert np.allclose(archive.pred_log10_HC50, saved.pred_log10_HC50, rtol=0, atol=1e-12)
        assert np.isclose(summary["hv_percent"], row.hv_percent)
        saved_history = pd.read_csv(opt / "runs" / row.run_id / "historia.csv")
        assert np.allclose(history.drop(columns="elapsed_seconds"), saved_history.drop(columns="elapsed_seconds"), rtol=0, atol=1e-12)
        report["repeated_real_run"] = {"seed": int(row.seed), "exact_sequence_order": True, "prediction_atol": 1e-12,
                                       "note": "Se excluye el tiempo del requisito de igualdad."}
    (ROOT / "results/auditoria_resultados.json").write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")
    print(json.dumps(report, indent=2, ensure_ascii=False))

if __name__ == "__main__":
    main()
