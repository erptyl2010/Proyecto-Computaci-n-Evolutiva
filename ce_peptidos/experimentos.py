from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd
from scipy.stats import wilcoxon
from .evolucion import Landscape, run_search, BASE

def experiments(root, models):
    out = root / "results" / "optimization"
    out.mkdir(parents=True, exist_ok=True)
    landscape = Landscape(models)
    oracle = pd.DataFrame({"sequence": landscape.sequences,
                           "pred_log10_MIC": landscape.F[:, 0], "pred_log10_HC50": -landscape.F[:, 1],
                           "pred_MIC_uM": 10**landscape.F[:, 0], "pred_HC50_uM": 10**(-landscape.F[:, 1]),
                           "n_changes": np.sum(landscape.X != landscape.X[0], axis=1),
                           "nondominated": landscape.front_mask})
    oracle.to_csv(out / "dominio_exhaustivo.csv", index=False)
    oracle[oracle.nondominated].to_csv(out / "frente_exhaustivo.csv", index=False)
    description = {"base": BASE, "n_sequences": len(oracle), "n_front_sequences": int(oracle.nondominated.sum()),
                   "n_front_objective_pairs": len(np.unique(landscape.F[landscape.front_mask], axis=0)),
                   "reference_point": landscape.reference.tolist(), "hv_exhaustive": landscape.hv_star,
                   "feature_preparation_seconds": landscape.feature_seconds,
                   "exhaustive_prediction_seconds": landscape.exhaustive_seconds,
                   "oracle_scope": "Optimo respecto a los dos modelos congelados, no respecto a resultados biologicos.",
                   "terminal_context": "SIN_ANOTACION en ambos extremos; no equivale a demostrar grupos terminales libres.",
                   "budget": 2100, "base_prediction": oracle.iloc[0].to_dict()}
    (out / "referencia.json").write_text(json.dumps(description, indent=2, ensure_ascii=False), encoding="utf-8")
    source_hash = hashlib.sha256((root / "ce_peptidos" / "evolucion.py").read_bytes()).hexdigest()
    model_hash = hashlib.sha256((root / "results" / "models" / "modelos_congelados.joblib").read_bytes()).hexdigest()
    def execute(method, seed, config):
        spec = {"method": method, "seed": seed, **config, "budget": 2100,
                "source_sha256": source_hash, "models_sha256": model_hash}
        digest = hashlib.sha256(json.dumps(spec, sort_keys=True).encode()).hexdigest()[:20]
        folder = out / "runs" / digest
        if (folder / "complete.json").exists():
            saved = json.loads((folder / "complete.json").read_text(encoding="utf-8"))
            if saved["spec"] != spec:
                raise RuntimeError("Clave de caché inconsistente.")
            for file, expected in saved["sha256"].items():
                if hashlib.sha256((folder / file).read_bytes()).hexdigest() != expected:
                    raise RuntimeError("Salida de una corrida modificada.")
            return saved["summary"] | {"run_id": digest}
        folder.mkdir(parents=True, exist_ok=True)
        summary, history, archive = run_search(landscape, method, seed, **config)
        history.to_csv(folder / "historia.csv", index=False)
        archive.to_csv(folder / "evaluaciones.csv", index=False)
        manifest = {"spec": spec, "summary": summary,
                    "sha256": {file: hashlib.sha256((folder / file).read_bytes()).hexdigest()
                               for file in ["historia.csv", "evaluaciones.csv"]}}
        (folder / "complete.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8")
        print(f"{method} seed={seed} pop={config['population']} mut={config['mutation']:.5f} "
              f"cross={config['crossover']:.2f}: HV={summary['hv_percent']:.4f}% ({summary['seconds']:.2f}s)", flush=True)
        return summary | {"run_id": digest}
    config = {"population": 60, "mutation": 1/21, "crossover": .9}
    sensitivity = []
    selections = []
    # Protocolo fijado antes de observar resultados: OAT secuencial, cinco semillas por opción.
    for parameter, values in [("population", [30, 60, 100]), ("mutation", [.5/21, 1/21, 2/21]),
                              ("crossover", [.5, .7, .9])]:
        stage = []
        for value in values:
            trial = config | {parameter: value}
            for seed in range(201, 206):
                row = execute("NSGA-II", seed, trial) | {"sensitivity_parameter": parameter, "sensitivity_value": value}
                sensitivity.append(row)
                stage.append(row)
        table = pd.DataFrame(stage).groupby("sensitivity_value").hv_percent.agg(["mean", "std"]).reset_index()
        selected = table.sort_values(["mean", "std", "sensitivity_value"], ascending=[False, True, True]).iloc[0]
        config[parameter] = int(selected.sensitivity_value) if parameter == "population" else float(selected.sensitivity_value)
        selections.append({"parameter": parameter, "selected": config[parameter], "reason": "Mayor media de HV; desempate por menor desviacion y menor valor."})
    pd.DataFrame(sensitivity).to_csv(out / "sensibilidad.csv", index=False)
    final = []
    for seed in range(401, 411):
        # Alternar orden reduce el sesgo temporal sistemático de medir siempre un método primero.
        for method in (["NSGA-II", "Aleatoria"] if seed % 2 else ["Aleatoria", "NSGA-II"]):
            final.append(execute(method, seed, config))
    table = pd.DataFrame(final)
    table.to_csv(out / "comparacion_final.csv", index=False)
    paired = table.pivot(index="seed", columns="method", values="hv_percent")
    delta = (paired["NSGA-II"] - paired["Aleatoria"]).to_numpy()
    p = 1. if np.all(delta == 0) else float(wilcoxon(delta, alternative="greater", method="auto").pvalue)
    rng = np.random.default_rng(20261017)
    interval = np.quantile(rng.choice(delta, size=(10000, len(delta)), replace=True).mean(axis=1), [.025, .975])
    stats = {"selected_configuration": config, "selection_steps": selections,
             "tuning_seeds": list(range(201, 206)), "final_seeds": list(range(401, 411)),
             "sensitivity_design": "Un factor a la vez, secuencial; no estima interacciones ni garantiza optimo global de parametros.",
             "paired_mean_difference_percentage_points": float(delta.mean()),
             "paired_bootstrap_95_percent_interval": interval.tolist(),
             "wilcoxon_alternative": "greater: NSGA-II > aleatoria; hipotesis fijada antes de ejecutar",
             "wilcoxon_pvalue": p, "n_paired_seeds": len(delta),
             "scope": "Variabilidad de busqueda en un unico paisaje predicho congelado; no incluye incertidumbre de entrenamiento o medicion.",
             "per_method": {method: {"mean_hv_percent": float(g.hv_percent.mean()),
                                       "std_hv_percent_ddof1": float(g.hv_percent.std()),
                                       "mean_seconds": float(g.seconds.mean())} for method, g in table.groupby("method")}}
    (out / "estadisticas.json").write_text(json.dumps(stats, indent=2, ensure_ascii=False), encoding="utf-8")
    return stats
