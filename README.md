# Proyecto de Computación Evolutiva: péptidos y NSGA-II

Esta es una reconstrucción reproducible del proyecto, separada del informe y los notebooks originales. Estudia la búsqueda de variantes con dos propiedades **predichas**. No contiene validación experimental de los candidatos ni una nueva ejecución de ColabFold.

## Equipo y repositorio

Ing. MIGUEL HERRERA, WILLY ROBERTO; Ing. PARIHUANA TRAVEZAÑO, OMAR EDUARDO; Ing. TORRES CRISTOBAL, WILBER.

Repositorio: https://github.com/erptyl2010/Proyecto-Computaci-n-Evolutiva

## Qué leer primero

- [Informe académico en PDF, 13 páginas](output/pdf/Informe_CE_Peptidos.pdf).
- [Informe ampliado y editable en Markdown](docs/INFORME_RESULTADOS.md).
- [Cuaderno ejecutado de lectura y comprobación](notebooks/PROYECTO_CE_REPRODUCIBLE.ipynb).
- [Fase I para el parcial](docs/FASE_I.md).
- [Correspondencia con el enunciado y el sílabo](docs/TRAZABILIDAD_CURSO.md).
- [Auditoría final contra el enunciado](docs/AUDITORIA_FINAL_CIERRE.md).
- [Guion de defensa de 30 minutos](docs/GUION_DEFENSA.md).
- `results/auditoria_resultados.json`: verificaciones sobre los resultados reales.

## Reproducir en Windows

Desde esta carpeta, con Python instalado:

```powershell
python -m venv .venv
.\.venv\Scripts\python.exe -m pip install -r requirements-lock.txt
.\.venv\Scripts\python.exe -m pytest -q
.\.venv\Scripts\python.exe run.py --stage all
.\.venv\Scripts\python.exe scripts\verificar_resultados.py --repeat-one
```

La instantánea de datos originales está incluida: `run.py` no necesita descargar de nuevo. Si se parte únicamente del código, ejecutar primero `python scripts/descargar_dbaasp.py`, con acceso a Internet. La API puede cambiar, por lo que una descarga futura no garantiza la misma cohorte. No se reemplazan los JSON existentes; para otro estudio se debe crear una copia separada sin caché.

En Linux/Colab, instalar `requirements.txt` en un entorno propio y ejecutar `python run.py --stage all` desde la raíz del proyecto. El bloqueo completo se obtuvo en el entorno local indicado en `results/execution_manifest.json`; la ejecución en Colab no se ha certificado. No se requiere GPU. Los modelos `.joblib` deben provenir de esta entrega o de una fuente de confianza.

Las etapas disponibles son `data`, `train`, `experiment`, `report` y `all`. Las tres primeras tienen dependencias: datos completos → datos procesados → modelos → experimentos. La etapa `report` regenera la narración y las figuras a partir de los CSV existentes. El programa falla si faltan insumos; no rellena resultados.

El cuaderno muestra salidas de celdas de lectura ejecutadas con Python por `scripts/crear_cuaderno.py`; no se atribuye a ese cuaderno una nueva ejecución del entrenamiento. Para regenerar el PDF, instalar adicionalmente `requirements-pdf.txt` y ejecutar `python scripts/crear_pdf.py` después de la auditoría. El cuaderno se regenera con `python scripts/crear_cuaderno.py`.

## Estructura y trazabilidad

| Carpeta | Contenido |
|---|---|
| `data/raw` | Esquema API, páginas de consulta, fichas originales, IDs, fechas, filtros y SHA-256 |
| `data/processed` | Todas las mediciones auditadas, inclusiones/exclusiones, conjuntos por objetivo, referencias de fichas |
| `ce_peptidos` | Procesamiento, modelos, operadores, experimentos e informe |
| `results/models` | Selección de representación, particiones, reserva interna, modelos congelados |
| `results/optimization` | Dominio completo, frentes, sensibilidad, comparación, estadística y candidatos |
| `results/optimization/runs` | Especificación, semilla, evaluaciones únicas e historia de cada corrida |
| `results/figures` | Seis figuras generadas exclusivamente desde estos resultados |
| `tests` | Comprobaciones de datos, operadores, métrica y reproducibilidad del software |
| `docs` | Fase I, informe, correspondencia con el curso y guion |

Las corridas se identifican por una huella de configuración, semilla, código del buscador y modelos. Una configuración idéntica reutiliza resultados cuyos archivos se verifican por SHA-256. No se contabiliza esa reutilización como una réplica independiente. Para cronometrar de nuevo todas las corridas, guardar primero los resultados existentes y usar otra carpeta de proyecto sin `results/optimization/runs`.

## Decisiones metodológicas

La cohorte incluye todas las fichas devueltas por el filtro declarado de DBAASP frente a *E. coli* ATCC 25922; HC50 se extrae de esa misma cohorte. Solo se aceptan concentraciones puntuales positivas en µM. Los extremos sin anotación se marcan expresamente, sin tratarlos como una modificación química conocida. La mediana combina condiciones no armonizadas: se conserva esa limitación.

Los modelos se eligen con validación agrupada dentro del entrenamiento y se evalúan una vez en una reserva interna. No se reajustan con la reserva antes de la optimización. El agrupamiento por Hamming y longitud no elimina toda homología ni la dependencia por publicaciones. El objetivo del algoritmo evolutivo es el paisaje de esos modelos congelados.

NSGA-II y el muestreo uniforme sin reemplazo consultan 2 100 secuencias únicas por corrida. Comparten la población inicial y se comparan en diez semillas finales distintas de las de sensibilidad. El frente exhaustivo de 76 210 secuencias es una referencia del modelo. El hipervolumen relativo no mide eficacia biológica ni exactitud del predictor.

## Entrega académica

El desarrollo cubre las cuatro fases en código y análisis, con un PDF de 13 páginas revisado visualmente. Incluye integrantes y enlace al repositorio proporcionados por el equipo. El guion prepara 30 minutos de sustentación. Falta ensayar y realizar la exposición, entregar en el medio indicado por el docente y confirmar la discrepancia de semana 15/16. No se presupone que el parcial de semana 08 ya haya sido entregado.

Las referencias académicas y de software están al final del informe. Las mediciones se respaldan en las fichas originales y sus referencias bibliográficas, no en valores inventados. Los documentos del curso y los archivos históricos permanecen en las carpetas originales del usuario.

## Auditoría de cierre

Una segunda ejecución desde código e instantánea cruda, sin modelos ni corridas previos, reprodujo 122 tablas y las 55 corridas en la misma máquina y entorno. Se compararon evaluaciones, historias y métricas con tolerancia absoluta 1e-12; los tiempos se excluyeron de la igualdad. Los modelos coincidieron por SHA-256. Ver `results/auditoria_cierre.json` y `results/reproduccion_desde_cero.log`. Esto no certifica portabilidad a otro sistema ni validez biológica.

Para recalcular los controles descriptivos: `python scripts/auditoria_cierre.py`. Para verificar una reproducción independiente previamente ejecutada: `python scripts/auditoria_cierre.py --reproduction-root RUTA_DE_LA_COPIA`. El modo sin esa opción informa explícitamente que no verificó una reproducción limpia. Guardar la evidencia de cierre antes de volver a generar sus archivos.

El PDF de cierre requiere `auditoria_cierre.json` con reproducción limpia verificada. Para repetir esa comprobación, copiar a otra carpeta únicamente `ce_peptidos`, `scripts`, `run.py`, `requirements.txt` y `data/raw`; ejecutar allí `run.py --stage all` con las mismas versiones y pasar su ruta a `--reproduction-root` desde esta carpeta. El ZIP incluye el resultado auditado, sin el entorno virtual ni la copia temporal. `manifest_entrega.json` identifica los archivos de la entrega; `execution_manifest.json` registra el entorno y los archivos en el momento de su etapa de ejecución.
