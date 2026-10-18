"""Crea y ejecuta un cuaderno de lectura sin inventar salidas de entrenamiento.

Las celdas Python se ejecutan secuencialmente en este proceso; el cuaderno declara
ese motor y usa solamente Python estándar, pandas y archivos ya producidos.
"""
from pathlib import Path
import contextlib
import io
import json
import nbformat as nbf

ROOT = Path(__file__).resolve().parents[1]

def main():
    cells = []
    def md(s):
        cells.append(nbf.v4.new_markdown_cell(s))
    def code(s):
        cells.append(nbf.v4.new_code_cell(s))
    md("""# Proyecto de Computación Evolutiva · lectura verificable

Este cuaderno recorre los resultados de la **reconstrucción nueva**, no las salidas históricas.
Las celdas de lectura se ejecutaron secuencialmente al generar este archivo.
El entrenamiento y la optimización se ejecutan mediante `run.py`; no se afirma que este cuaderno los haya repetido.

Para reconstruir todo, ejecutar desde la raíz del proyecto:
```python
import subprocess, sys
subprocess.run([sys.executable, 'run.py', '--stage', 'all'], check=True)
```
La instantánea de DBAASP debe estar completa. Instalar antes las dependencias de `requirements.txt`.
""")
    code("""from pathlib import Path
import json
import pandas as pd
ROOT = Path.cwd()
if not (ROOT / 'run.py').exists():
    ROOT = ROOT.parent
assert (ROOT / 'run.py').exists(), 'Abrir desde la raíz del proyecto o desde notebooks/'
print('Raíz:', ROOT)
def read_json(relative):
    return json.loads((ROOT / relative).read_text(encoding='utf-8'))
""")
    md("## 1. Problema y procedencia\nSe minimiza log10 MIC y se maximiza log10 HC50 predichas. El dominio tiene 21 posiciones y hasta dos cambios. Las concentraciones generadas por el modelo no son nuevas mediciones.")
    code("""raw = read_json('data/raw/manifest.json')
print(json.dumps({k: raw[k] for k in ['status', 'started_utc', 'finished_utc', 'filters', 'ids', 'scope']}, ensure_ascii=False, indent=2))
print(json.dumps(read_json('data/processed/resumen_datos.json'), ensure_ascii=False, indent=2))
""")
    md("## 2. Modelo del proceso\nLa representación se selecciona dentro del entrenamiento. La reserva es interna, agrupada por similitud; no equivale a validación externa ni evita toda dependencia por publicaciones.")
    code("""print(pd.read_csv(ROOT / 'results/models/validacion_reservada.csv').to_string(index=False))
print(json.dumps(read_json('results/models/protocolo.json'), ensure_ascii=False, indent=2))
""")
    md("![Validación interna](../results/figures/validacion_predictores.png)")
    md("## 3. Sensibilidad y comparación\nNSGA-II tiene cruce uniforme, mutación categórica y reparación. Se cuentan secuencias únicas. El ajuste de parámetros y las corridas finales usan semillas diferentes.")
    code("""s = pd.read_csv(ROOT / 'results/optimization/sensibilidad.csv')
print(s.groupby(['sensitivity_parameter','sensitivity_value']).agg(HV_media=('hv_percent','mean'), HV_DE=('hv_percent','std'), segundos=('seconds','mean')).to_string())
final = pd.read_csv(ROOT / 'results/optimization/comparacion_final.csv')
print(final[['method','seed','unique_evaluations','hv_percent','seconds']].to_string(index=False))
print(json.dumps(read_json('results/optimization/estadisticas.json'), ensure_ascii=False, indent=2))
""")
    md("![Sensibilidad](../results/figures/sensibilidad.png)\n\n![Comparación](../results/figures/comparacion_calidad_tiempo.png)")
    md("## 4. Convergencia\nEl promedio y el mejor valor se registran por objetivo en cada población superviviente. Los extremos pueden corresponder a individuos distintos. El hipervolumen usa todo el archivo acumulado.")
    code("""h = pd.read_csv(ROOT / 'results/optimization/historia_comparacion.csv')
example = h[(h.method == 'NSGA-II') & (h.seed == 401)]
print(example[['generation_or_batch','n_evaluations','mean_log10_MIC','best_log10_MIC','mean_log10_HC50','best_log10_HC50','hv_percent']].to_string(index=False))
""")
    md("![Aptitud por generación](../results/figures/aptitud_promedio_y_mejor.png)\n\n![Hipervolumen](../results/figures/convergencia_hipervolumen.png)")
    md("## 5. Candidatos y comprobaciones\nLa procedencia indica en qué semillas finales NSGA-II evaluó cada candidato. El cociente HC50/MIC sigue siendo una razón de predicciones.")
    code("""c = pd.read_csv(ROOT / 'results/optimization/candidatos_encontrados.csv')
print(c[['criterio','sequence','mutations_1based','pred_MIC_uM','pred_HC50_uM','seeds_found']].to_string(index=False))
print(json.dumps(read_json('results/auditoria_resultados.json'), ensure_ascii=False, indent=2))
""")
    md("![Frentes](../results/figures/frentes_pareto.png)\n\nLa interpretación y bibliografía verificable se encuentran en [el informe](../docs/INFORME_RESULTADOS.md). No se han realizado ensayos biológicos ni una nueva ejecución de ColabFold.")
    env = {}
    import os
    old = Path.cwd()
    os.chdir(ROOT)
    try:
        count = 0
        for cell in cells:
            if cell.cell_type == "code":
                count += 1
                buf = io.StringIO()
                with contextlib.redirect_stdout(buf):
                    exec(compile(cell.source, f"cuaderno_celda_{count}", "exec"), env)
                cell.execution_count = count
                cell.outputs = [nbf.v4.new_output("stream", name="stdout", text=buf.getvalue())]
    finally:
        os.chdir(old)
    notebook = nbf.v4.new_notebook(cells=cells, metadata={
        "kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
        "language_info": {"name": "python", "version": "3.12"},
        "execution_note": "Celdas de lectura ejecutadas en orden mediante Python exec en scripts/crear_cuaderno.py. Entrenamiento y optimizacion ejecutados por run.py, no repetidos en este cuaderno."})
    nbf.validate(notebook)
    folder = ROOT / "notebooks"
    folder.mkdir(exist_ok=True)
    nbf.write(notebook, folder / "PROYECTO_CE_REPRODUCIBLE.ipynb")
    print(f"Cuaderno verificado: {count} celdas de lectura ejecutadas.")

if __name__ == "__main__":
    main()
