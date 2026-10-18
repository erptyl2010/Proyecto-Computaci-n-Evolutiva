"""Descarga una cohorte pública completa y conserva los JSON originales.

Solo biblioteca estándar. Ejecutar de nuevo reanuda descargas; no reemplaza
registros ya conservados. No usa los resultados históricos como datos nuevos.
"""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
import hashlib
import json
import time
import urllib.parse
import urllib.request

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
API = "https://dbaasp.org"
FILTERS = {"sequenceLength.value": "10-30", "complexity.value": "Monomer",
           "unusualAminoAcid.value": "without_modification",
           "targetSpecies.value": "Escherichia coli ATCC 25922"}

def get(url):
    for attempt in range(4):
        try:
            request = urllib.request.Request(url, headers={"Accept": "application/json",
                "User-Agent": "CE-academic-reproducibility/1.0 (DBAASP public API)"})
            with urllib.request.urlopen(request, timeout=45) as response:
                data = response.read()
            json.loads(data)
            return data
        except Exception:
            if attempt == 3:
                raise
            time.sleep(2 ** attempt)

def cached(url, path):
    if not path.exists():
        data = get(url)
        temp = path.with_suffix(".part")
        temp.write_bytes(data)
        temp.replace(path)
    return json.loads(path.read_bytes())

def main():
    RAW.mkdir(parents=True, exist_ok=True)
    pages = RAW / "pages"
    peptides = RAW / "peptides"
    pages.mkdir(exist_ok=True)
    peptides.mkdir(exist_ok=True)
    metadata_path = RAW / "manifest.json"
    manifest = json.loads(metadata_path.read_text(encoding="utf-8")) if metadata_path.exists() else {
        "started_utc": datetime.now(timezone.utc).isoformat(), "filters": FILTERS,
        "source": API, "status": "in_progress",
        "scope": "Todos los IDs devueltos por el filtro declarado. HC50 se extrae de esta misma cohorte; no reproduce el conjunto historico de 181/153 registros."}
    if manifest["filters"] != FILTERS:
        raise RuntimeError("El cache pertenece a otros filtros; use otra carpeta para otra cohorte.")
    metadata_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")
    cached(API + "/v3/api-docs", RAW / "api_schema.json")
    false_filters = FILTERS | {"targetSpecies.value": "Organismo_Inexistente_XYZ_987654",
                               "limit": 1, "offset": 0}
    control = cached(API + "/peptides?" + urllib.parse.urlencode(false_filters), RAW / "negative_control.json")
    if int(control["totalCount"]) != 0:
        raise RuntimeError("El control negativo no valida el filtro de especie.")
    offset = 0
    records = {}
    total = None
    while total is None or offset < total:
        url = API + "/peptides?" + urllib.parse.urlencode(FILTERS | {"offset": offset, "limit": 100})
        page = cached(url, pages / f"page_{offset:06d}.json")
        current_total = int(page["totalCount"])
        if total is None:
            total = current_total
        elif current_total != total:
            raise RuntimeError("Cambio el total durante la paginacion; no se considera una descarga completa.")
        data = page["data"]
        if not data:
            raise RuntimeError("Pagina vacia antes de completar el total.")
        for item in data:
            records[int(item["id"])] = item
        offset += len(data)
        print(f"Listado: {offset}/{total}", flush=True)
        time.sleep(0.1)
    if len(records) != total:
        raise RuntimeError("IDs duplicados o faltantes en la paginacion.")
    ids = sorted(records)
    (RAW / "cohort_ids.json").write_text(json.dumps(ids), encoding="utf-8")
    manifest.update({"total_api": total, "ids": len(ids), "status": "downloading_records"})
    metadata_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")

    def download(pid):
        path = peptides / f"{pid}.json"
        payload = cached(f"{API}/peptides/{pid}", path)
        if int(payload["id"]) != pid:
            raise RuntimeError(f"Respuesta no corresponde al ID {pid}")
        time.sleep(0.15)
        return pid

    failures = []
    with ThreadPoolExecutor(max_workers=4) as pool:
        tasks = {pool.submit(download, pid): pid for pid in ids}
        for count, task in enumerate(as_completed(tasks), 1):
            try:
                task.result()
            except Exception as error:
                failures.append({"id": tasks[task], "error": str(error)})
            if count % 100 == 0 or count == len(ids):
                print(f"Registros: {count}/{len(ids)}; fallos: {len(failures)}", flush=True)
    hashes = {str(path.relative_to(RAW)).replace("\\", "/"): hashlib.sha256(path.read_bytes()).hexdigest()
              for path in sorted(RAW.rglob("*.json")) if path != metadata_path}
    manifest.update({"status": "complete" if not failures else "incomplete",
                     "finished_utc": datetime.now(timezone.utc).isoformat(),
                     "failed_records": failures, "sha256": hashes})
    metadata_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")
    print("Estado:", manifest["status"], flush=True)
    if failures:
        raise SystemExit(1)

if __name__ == "__main__":
    main()
