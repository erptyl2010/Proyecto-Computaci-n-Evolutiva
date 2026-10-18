"""Paquete portable con manifiesto SHA-256; excluye entornos y temporales."""
from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json
import os
from zipfile import ZipFile, ZIP_DEFLATED

ROOT = Path(__file__).resolve().parents[1]
EXCLUDE = {".venv", ".tmp", ".pytest_cache", "__pycache__", ".git"}

def main():
    files = []
    for folder, dirs, names in os.walk(ROOT):
        dirs[:] = [d for d in dirs if d not in EXCLUDE]
        for name in names:
            path = Path(folder) / name
            if name != "manifest_entrega.json" and path.suffix not in {".pyc", ".part"}:
                files.append(path)
    files.sort()
    manifest = {"created_utc": datetime.now(timezone.utc).isoformat(),
                "description": "Reconstruccion ejecutada con instantanea publica, codigo, modelos, resultados y documentos.",
                "sha256": {p.relative_to(ROOT).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest() for p in files}}
    mpath = ROOT / "manifest_entrega.json"
    mpath.write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")
    target = ROOT.parent / "Entrega_Proyecto_CE_Peptidos.zip"
    with ZipFile(target, "w", ZIP_DEFLATED, compresslevel=6) as archive:
        for path in files + [mpath]:
            archive.write(path, ROOT.name + "/" + path.relative_to(ROOT).as_posix())
    with ZipFile(target) as archive:
        assert archive.testzip() is None
        for rel, expected in manifest["sha256"].items():
            assert hashlib.sha256(archive.read(ROOT.name + "/" + rel)).hexdigest() == expected, rel
    print(json.dumps({"zip": str(target), "files": len(files)+1,
                      "bytes": target.stat().st_size, "sha256": hashlib.sha256(target.read_bytes()).hexdigest()}, indent=2))

if __name__ == "__main__":
    main()
