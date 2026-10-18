from collections import Counter
from pathlib import Path
import hashlib
import json
import math
import re
import unicodedata

import numpy as np
import pandas as pd

AA = "ACDEFGHIKLMNPQRSTVWY"
TARGET = "Escherichia coli ATCC 25922"

def name(value):
    return (value or {}).get("name", "") if isinstance(value, (dict, type(None))) else str(value)

def concentration(value):
    text = str(value).strip()
    if not re.fullmatch(r"\+?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?", text):
        return None
    result = float(text)
    return result if math.isfinite(result) and result > 0 else None

def is_micromolar(unit):
    return unicodedata.normalize("NFKC", unit.strip()) in {"μM", "uM"}

def peptide_reason(record):
    seq = str(record.get("sequence", "")).strip()
    if not re.fullmatch("[" + AA + "]{10,30}", seq):
        return "secuencia_no_canonica_o_longitud"
    if name(record.get("complexity")).lower() != "monomer":
        return "no_monomero"
    if any(record.get(k) for k in ["unusualAminoAcids", "intrachainBonds", "interchainBonds", "coordinationBonds"]):
        return "quimica_no_representada"
    nt = name(record.get("nTerminus") or record.get("NTerminus"))
    ct = name(record.get("cTerminus") or record.get("CTerminus"))
    if nt not in {"", "ACT"} or ct not in {"", "AMD"}:
        return "terminal_no_representado"
    return ""

def build_datasets(root):
    raw = root / "data" / "raw"
    out = root / "data" / "processed"
    out.mkdir(parents=True, exist_ok=True)
    manifest = json.loads((raw / "manifest.json").read_text(encoding="utf-8"))
    if manifest["status"] != "complete":
        raise RuntimeError("La descarga no esta completa; no se entrena sobre una cohorte parcial.")
    for rel, expected in manifest["sha256"].items():
        if hashlib.sha256((raw / rel).read_bytes()).hexdigest() != expected:
            raise RuntimeError(f"El archivo original cambio: {rel}")
    ids = json.loads((raw / "cohort_ids.json").read_text())
    audit, accepted, articles = [], [], {}
    for pid in ids:
        p = json.loads((raw / "peptides" / f"{pid}.json").read_text(encoding="utf-8"))
        reason = peptide_reason(p)
        nt = name(p.get("nTerminus") or p.get("NTerminus"))
        ct = name(p.get("cTerminus") or p.get("CTerminus"))
        article_ids = sorted(int(a["id"]) for a in p.get("articles", []))
        for article in p.get("articles", []):
            articles[int(article["id"])] = article
        for objective, key in [("MIC", "targetActivities"), ("HC50", "hemoliticCytotoxicActivities")]:
            for a in p.get(key, []):
                target = name(a.get("targetSpecies")) if objective == "MIC" else name(a.get("targetCell"))
                measure = a.get("activityMeasureValue", "") if objective == "MIC" else a.get("activityMeasureForLysisValue", "")
                expected_target = TARGET if objective == "MIC" else "Human erythrocytes"
                expected_measure = "MIC" if objective == "MIC" else "50% Hemolysis"
                value = concentration(a.get("concentration"))
                why = reason
                if not why and target != expected_target:
                    why = "otro_objetivo_biologico"
                if not why and measure != expected_measure:
                    why = "otra_medida"
                if not why and not is_micromolar(name(a.get("unit"))):
                    why = "unidad_no_micromolar"
                if not why and value is None:
                    why = "concentracion_no_puntual_positiva"
                item = {"peptide_id": pid, "measurement_id": a.get("id"),
                        "objective": objective, "sequence": p.get("sequence"),
                        "n_terminus": nt or "SIN_ANOTACION", "c_terminus": ct or "SIN_ANOTACION",
                        "target": target, "measure": measure, "concentration_raw": a.get("concentration"),
                        "unit_raw": name(a.get("unit")), "value_uM": value,
                        "medium": name(a.get("medium")), "cfu": a.get("cfu"), "ph": a.get("ph"),
                        "ionic_strength": a.get("ionicStrength"), "reference_raw": a.get("reference"),
                        "article_ids_all": json.dumps(article_ids),
                        "source_url": f"https://dbaasp.org/peptides/{pid}",
                        "raw_sha256": manifest["sha256"][f"peptides/{pid}.json"],
                        "included": not bool(why), "exclusion_reason": why}
                audit.append(item)
                if not why:
                    accepted.append(item)
    audit_df = pd.DataFrame(audit)
    audit_df.to_csv(out / "auditoria_mediciones.csv", index=False)
    pd.DataFrame(accepted).to_csv(out / "mediciones_aceptadas.csv", index=False)
    (out / "referencias_dbaasp.json").write_text(json.dumps(articles, ensure_ascii=False, indent=2), encoding="utf-8")
    datasets = {}
    accepted_df = pd.DataFrame(accepted)
    for objective in ["MIC", "HC50"]:
        subset = accepted_df[accepted_df.objective == objective]
        rows = []
        for (sequence, nt, ct), group in subset.groupby(["sequence", "n_terminus", "c_terminus"], sort=True):
            identity = f"{sequence}|{nt}|{ct}"
            pubs = sorted({a for value in group.article_ids_all for a in json.loads(value)})
            rows.append({"chemical_id": hashlib.sha256(identity.encode()).hexdigest()[:16],
                         "sequence": sequence, "n_terminus": nt, "c_terminus": ct,
                         "value_uM": float(group.value_uM.median()), "n_measurements": len(group),
                         "peptide_ids": json.dumps(sorted(set(group.peptide_id.astype(int)))),
                         "article_ids_all": json.dumps(pubs),
                         "n_conditions": len(group[["medium", "cfu", "ph", "ionic_strength"]].drop_duplicates())})
        datasets[objective] = pd.DataFrame(rows)
        datasets[objective].to_csv(out / f"dataset_{objective}.csv", index=False)
    summary = {"cohort_api_ids": len(ids), "measurements_audited": len(audit),
               "accepted_measurements": {k: sum(r["objective"] == k for r in accepted) for k in datasets},
               "chemical_entities": {k: len(v) for k, v in datasets.items()},
               "exclusion_reasons": dict(Counter(r["exclusion_reason"] for r in audit if not r["included"])),
               "aggregation": "Mediana de mediciones puntuales por secuencia y anotaciones terminales, combinando IDs equivalentes. Las condiciones no se armonizan.",
               "scope": manifest["scope"],
               "reference_mapping": "article_ids_all incluye todos los articulos de la ficha. reference_raw conserva la referencia de la medicion; no se asume una correspondencia no verificada."}
    (out / "resumen_datos.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    return datasets, summary

def similarity_groups(sequences, changes=2):
    """Componentes transitivas de Hamming <=2, solamente a igual longitud."""
    parent = np.arange(len(sequences))
    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    buckets = {}
    for i, seq in enumerate(sequences):
        buckets.setdefault(len(seq), []).append(i)
    for indices in buckets.values():
        encoded = np.array([list(sequences[i]) for i in indices])
        for i in range(len(indices)):
            neighbours = np.flatnonzero(np.count_nonzero(encoded[i + 1:] != encoded[i], axis=1) <= changes)
            for relative in neighbours:
                a, b = find(indices[i]), find(indices[i + 1 + relative])
                parent[b] = a
    return np.array([find(i) for i in range(len(sequences))])
