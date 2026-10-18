# Correspondencia con el curso

Fuentes: `../Curso/Enunciado_Proyecto_Final_Computacion_Evolutiva.docx` y `../Curso/SILABO 402 - COMPUTACIÓN EVOLUTIVA_firma.pdf`, tomando como raíz `Proyecto_CE_Peptidos/`. No se anticipa una calificación.

| Requisito | Implementación/evidencia del proyecto | Observación |
|---|---|---|
| Fase I: problema y formulación | `docs/FASE_I.md`; sección I del informe | 21 variables categóricas, dos objetivos y restricción Hamming |
| Modelo del proceso | `ce_peptidos/datos.py`, `modelos.py`; `results/models/` | Predictores empíricos con prueba interna reservada |
| Paradigma y representación | `ce_peptidos/evolucion.py` | NSGA-II sobre secuencias, una representación por secuencia |
| Operadores y selección | Cruce uniforme, mutación categórica, reparación; pymoo 0.6.1.5 | Se explica el sesgo de reparación |
| Implementación funcional | `run.py`, módulos y pruebas | Se ejecuta por etapas sin depender de variables de Colab |
| Métrica de evaluación | Hipervolumen, frente exhaustivo y pruebas independientes | Diferente del vector de aptitud y de la exactitud del predictor |
| Al menos una sensibilidad | Población, mutación y cruce; `sensibilidad.csv` | Diseño secuencial de un factor; no factorial completo |
| Promedio y mejor aptitud | `results/figures/aptitud_promedio_y_mejor.png` | Por objetivo, en población superviviente; extremos pueden ser individuos diferentes |
| Convergencia y eficiencia | HV por evaluaciones, tiempo y umbrales | Se explican costos excluidos y corridas que no alcanzan cada umbral |
| Solución e interpretación | `candidatos_encontrados.csv` | Candidatos realmente encontrados; predicciones no equivalen a mediciones |
| Comparación de referencia | Búsqueda aleatoria y barrido de 76 210 secuencias | Igual presupuesto e inicio compartido por pareja |
| Conclusiones y futuro | Sección IV del informe | Diferencia desempeño del optimizador y validez del evaluador |
| Reproducibilidad | JSON originales, SHA-256, particiones, modelos, semillas, versiones, CSV | Reproducción desde la instantánea; el sitio público puede cambiar |
| Informe de hasta unas 15 páginas | `output/pdf/Informe_CE_Peptidos.pdf` y fuente ampliada en Markdown | PDF de 13 páginas; integrantes y URL incorporados |
| Link de repositorio en informe | URL indicada por el equipo e incorporada al PDF | https://github.com/erptyl2010/Proyecto-Computaci-n-Evolutiva |
| Exposición oral | `docs/GUION_DEFENSA.md` | Guion de 30 minutos y preguntas preparado; ensayo y exposición pendientes |

La fase I corresponde al parcial de semana 08. Hay una discrepancia entre semana 15 del enunciado y semana 16 del sílabo para el proyecto final; debe confirmarse con el docente, sin deducir una fecha calendario.

La rúbrica del sílabo revisado asigna 4 puntos a formulación, 4 a diseño, 4 a experimentación, 3 a reproducibilidad/documentación/responsabilidad, 2 a informe y 3 a defensa. La tabla ubica evidencias, no certifica cumplimiento de la entrega administrativa ni promete una nota.

La auditoría final detallada está en `docs/AUDITORIA_FINAL_CIERRE.md`. La reproducción limpia comparó 122 tablas y las 55 corridas. La discusión del presupuesto, costo completo del ajuste y recuperación del frente se amplió en la página 9 del PDF.
