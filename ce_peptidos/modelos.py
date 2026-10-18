from pathlib import Path
import json
import hashlib
import re
import unicodedata

import joblib
import numpy as np
import pandas as pd
from sklearn.dummy import DummyRegressor
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import mean_absolute_error, r2_score
from sklearn.model_selection import GroupShuffleSplit, GroupKFold

from .datos import AA, similarity_groups

HYDRO = dict(zip(AA, [1.8, 2.5, -3.5, -3.5, 2.8, -0.4, -3.2, 4.5, -3.9, 3.8,
                       1.9, -3.5, -1.6, -3.5, -4.5, -0.8, -0.7, 4.2, -0.9, -1.3]))
CHARGE = {aa: {"K": 1., "R": 1., "H": .1, "D": -1., "E": -1.}.get(aa, 0.) for aa in AA}

def features(sequence, nt="SIN_ANOTACION", ct="SIN_ANOTACION", ordered=False):
    if len(sequence) < 3 or any(a not in AA for a in sequence):
        raise ValueError("Secuencia no compatible con los descriptores.")
    n = len(sequence)
    f = {"length": n} | {f"fraction_{a}": sequence.count(a) / n for a in AA}
    f.update({"charge_approx": sum(CHARGE[a] for a in sequence),
              "hydrophobic_fraction": sum(a in "AILMFWVY" for a in sequence) / n,
              "n_acetylated": int(nt == "ACT"), "c_amidated": int(ct == "AMD"),
              "n_annotation_missing": int(nt == "SIN_ANOTACION"),
              "c_annotation_missing": int(ct == "SIN_ANOTACION")})
    if ordered:
        for label, scale in [("hydro", HYDRO), ("charge", CHARGE),
                             ("aromatic", {a: float(a in "FWY") for a in AA})]:
            x = np.array([scale[a] for a in sequence])
            f.update({f"{label}_n3": x[:3].mean(), f"{label}_c3": x[-3:].mean(),
                      f"{label}_position": np.mean(x * np.linspace(-1, 1, n)),
                      f"{label}_lag1": np.mean(x[:-1] * x[1:]),
                      f"{label}_lag2": np.mean(x[:-2] * x[2:])})
    return f

def matrix(df, ordered):
    return pd.DataFrame([features(r.sequence, r.n_terminus, r.c_terminus, ordered)
                         for r in df.itertuples()])

def train(root, datasets, seed=20261017):
    out = root / "results" / "models"
    out.mkdir(parents=True, exist_ok=True)
    models, cv_rows, metrics, predictions, membership = {}, [], [], [], []
    references = json.loads((root / "data/processed/referencias_dbaasp.json").read_text(encoding="utf-8"))
    def publication_keys(ids):
        keys = set()
        for article_id in ids:
            article = references.get(str(article_id), {})
            title = re.sub(r"[^\w]+", "", unicodedata.normalize("NFKC", article.get("title", "")).casefold())
            if title:
                keys.add((title, str(article.get("year", ""))))
        return keys
    for objective, df in datasets.items():
        if len(df) < 40:
            raise RuntimeError(f"{objective}: menos de 40 entidades; revisar el alcance antes de entrenar.")
        groups = similarity_groups(df.sequence.tolist())
        y = np.log10(df.value_uM.to_numpy())
        splitter = GroupShuffleSplit(n_splits=1, test_size=.2, random_state=seed)
        train_idx, test_idx = next(splitter.split(df, y, groups))
        if len(test_idx) < 2 or len(np.unique(groups[train_idx])) < 5:
            raise RuntimeError("Particion insuficiente para validacion agrupada de cinco pliegues.")
        representations = {key: matrix(df, key == "ordered") for key in ["composition", "ordered"]}
        for label, X in representations.items():
            splits = GroupKFold(n_splits=5).split(X.iloc[train_idx], y[train_idx], groups[train_idx])
            for fold, (a, b) in enumerate(splits):
                ia, ib = train_idx[a], train_idx[b]
                model = RandomForestRegressor(n_estimators=120, min_samples_leaf=3, n_jobs=1,
                                              random_state=seed + fold)
                model.fit(X.iloc[ia].to_numpy(), y[ia])
                cv_rows.append({"objective": objective, "representation": label, "fold": fold + 1,
                                "mae_log10": mean_absolute_error(y[ib], model.predict(X.iloc[ib].to_numpy())),
                                "n_train": len(ia), "n_validation": len(ib)})
        current = pd.DataFrame([r for r in cv_rows if r["objective"] == objective])
        selected = current.groupby("representation").mae_log10.mean().sort_values(kind="stable").index[0]
        X = representations[selected]
        model = RandomForestRegressor(n_estimators=120, min_samples_leaf=3, n_jobs=1, random_state=seed)
        model.fit(X.iloc[train_idx].to_numpy(), y[train_idx])
        baseline = DummyRegressor(strategy="mean").fit(X.iloc[train_idx].to_numpy(), y[train_idx])
        pred = model.predict(X.iloc[test_idx].to_numpy())
        dummy = baseline.predict(X.iloc[test_idx].to_numpy())
        articles_train = {v for value in df.iloc[train_idx].article_ids_all for v in json.loads(value)}
        articles_test = {v for value in df.iloc[test_idx].article_ids_all for v in json.loads(value)}
        models[objective] = {"model": model, "ordered": selected == "ordered", "features": X.columns.tolist(),
                             "n_fit": len(train_idx), "n_holdout": len(test_idx), "seed": seed}
        metrics.append({"objective": objective, "representation": selected, "n_features": X.shape[1],
                        "n_entities": len(df), "n_train": len(train_idx), "n_holdout": len(test_idx),
                        "n_groups": len(set(groups)), "mae_log10": mean_absolute_error(y[test_idx], pred),
                        "r2_log10": r2_score(y[test_idx], pred),
                        "mae_baseline_log10": mean_absolute_error(y[test_idx], dummy),
                        "shared_article_ids_train_holdout": len(articles_train & articles_test),
                        "shared_publications_train_holdout": len(publication_keys(articles_train) & publication_keys(articles_test))})
        for i, predicted, base_pred in zip(test_idx, pred, dummy):
            predictions.append({"objective": objective, "chemical_id": df.iloc[i].chemical_id,
                                "experimental_log10": y[i], "predicted_log10": predicted,
                                "baseline_log10": base_pred})
        train_set = set(train_idx)
        for i, row in df.iterrows():
            membership.append({"objective": objective, "chemical_id": row.chemical_id,
                               "group": int(groups[i]), "split": "train" if i in train_set else "holdout"})
    pd.DataFrame(cv_rows).to_csv(out / "seleccion_representacion_cv.csv", index=False)
    pd.DataFrame(metrics).to_csv(out / "validacion_reservada.csv", index=False)
    pd.DataFrame(predictions).to_csv(out / "predicciones_reservadas.csv", index=False)
    pd.DataFrame(membership).to_csv(out / "particiones.csv", index=False)
    joblib.dump(models, out / "modelos_congelados.joblib")
    protocol = {"seed": seed, "outer": "GroupShuffleSplit 20% de grupos, una particion reservada",
                "inner": "GroupKFold 5 pliegues solo en entrenamiento; menor MAE elige representacion",
                "grouping": "Componentes transitivas, igual longitud, Hamming <=2; no elimina toda homologia ni separa por publicacion",
                "scope": "Prueba interna reservada de la misma cohorte; no validacion biologica de variantes ni validacion externa",
                "publication_audit": "Coincidencia de titulo normalizado y anio entre todos los articulos de las fichas; los IDs internos distintos no garantizan publicaciones distintas. No verifica la atribucion de cada medicion.",
                "fit_after_selection": "Solo entrenamiento. No se reentrena con los casos reservados antes de optimizar.",
                "model_sha256": hashlib.sha256((out / "modelos_congelados.joblib").read_bytes()).hexdigest()}
    (out / "protocolo.json").write_text(json.dumps(protocol, ensure_ascii=False, indent=2), encoding="utf-8")
    return models, metrics
