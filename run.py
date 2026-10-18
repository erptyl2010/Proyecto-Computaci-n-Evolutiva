"""Punto de entrada: python run.py --stage all (datos crudos ya descargados)."""
import argparse
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import platform
import sys
from datetime import datetime, timezone
ROOT = Path(__file__).resolve().parent
os.environ.setdefault("MPLBACKEND", "Agg")
os.environ.setdefault("MPLCONFIGDIR", str(ROOT / ".tmp" / "matplotlib"))
import joblib
import pandas as pd

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--stage", choices=["all", "data", "train", "experiment", "report"], default="all")
    args = parser.parse_args()
    if args.stage in ["all", "data"]:
        from ce_peptidos.datos import build_datasets
        datasets, summary = build_datasets(ROOT)
        print(json.dumps(summary, ensure_ascii=False, indent=2), flush=True)
    if args.stage in ["all", "train"]:
        from ce_peptidos.modelos import train
        datasets = {k: pd.read_csv(ROOT / "data" / "processed" / f"dataset_{k}.csv", keep_default_na=False) for k in ["MIC", "HC50"]}
        models, metrics = train(ROOT, datasets)
        print(pd.DataFrame(metrics).to_string(index=False), flush=True)
    if args.stage in ["all", "experiment"]:
        from ce_peptidos.experimentos import experiments
        models = joblib.load(ROOT / "results" / "models" / "modelos_congelados.joblib")
        print(json.dumps(experiments(ROOT, models), ensure_ascii=False, indent=2), flush=True)
    if args.stage in ["all", "report"]:
        from ce_peptidos.informe import report
        report(ROOT)
    (ROOT / "results").mkdir(exist_ok=True)
    files = [p for p in ROOT.rglob("*") if p.is_file() and
             p.relative_to(ROOT).parts[0] in {"ce_peptidos", "scripts", "tests", "results", "docs"} and
             "__pycache__" not in p.parts and p.name != "execution_manifest.json"]
    files.extend([ROOT / "run.py", ROOT / "requirements.txt"])
    manifest = {"finished_utc": datetime.now(timezone.utc).isoformat(), "stage": args.stage,
                "python": sys.version, "platform": platform.platform(),
                "processor": platform.processor(), "logical_cpus": os.cpu_count(),
                "packages": {p: importlib.metadata.version(p) for p in ["numpy", "pandas", "scipy", "scikit-learn", "pymoo", "matplotlib"]},
                "sha256": {str(p.relative_to(ROOT)).replace("\\", "/"): hashlib.sha256(p.read_bytes()).hexdigest() for p in files}}
    (ROOT / "results" / "execution_manifest.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8")

if __name__ == "__main__":
    main()
