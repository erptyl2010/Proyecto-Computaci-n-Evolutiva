# Auditoría final del proyecto de Computación Evolutiva

Fecha local: 17 de octubre de 2026. Versión auditada: reconstrucción ejecutada con la instantánea pública incluida, no los notebooks históricos.

Equipo: Ing. MIGUEL HERRERA, WILLY ROBERTO; Ing. PARIHUANA TRAVEZAÑO, OMAR EDUARDO; Ing. TORRES CRISTOBAL, WILBER.

Repositorio indicado por el equipo: https://github.com/erptyl2010/Proyecto-Computaci-n-Evolutiva

## Dictamen

**El desarrollo técnico cubre las cuatro fases del enunciado y es defendible como un estudio de Computación Evolutiva sobre propiedades predichas.** La revisión de cierre reforzó la discusión del presupuesto, la sensibilidad en calidad y tiempo, el costo total del ajuste y la diferencia entre hipervolumen y recuperación exacta del frente. No se identificó un requisito técnico obligatorio que exija cambiar de tema o añadir otro algoritmo para cumplir el enunciado.

Esto no equivale a certificar péptidos eficaces o seguros, ni garantiza una nota. El punto científico más débil es el modelo de HC50, cuyo R² reservado es aproximadamente 0,059. Se puede defender cómo buscó NSGA-II en ese paisaje, pero no convertir ese resultado en una validación biológica. Si el equipo mantiene esta distinción, la conclusión es honesta y coherente con lo ejecutado.

El informe queda en 13 páginas, con los integrantes y el enlace proporcionados. El guion se ajusta a 30 minutos. El ensayo, la exposición y la entrega en el medio indicado por el docente todavía deben realizarse; no pueden darse por cumplidos mediante una auditoría documental. La publicación de esta versión se comprueba por separado con Git.

## Fuentes que fijan los requisitos

Se leyó el enunciado completo y las siete páginas del sílabo local. Sus huellas permiten identificar exactamente qué versiones se revisaron:

| Documento del curso | SHA-256 |
|---|---|
| Enunciado_Proyecto_Final_Computacion_Evolutiva.docx | `2412d6fe21adfee0a9e33d0def779ba37a96b32d245c1ee668946dcd46178e7f` |
| SILABO 402 - COMPUTACIÓN EVOLUTIVA_firma.pdf | `1ea33f868ec29781093b6decf48947e820beacfa687b1c96c180de7c1feade8e` |

Los documentos originales permanecen en la carpeta `Curso` del espacio local. Esta entrega no presupone autorización para redistribuir esos documentos docentes. La correspondencia con la rúbrica es una evaluación razonada, no una calificación oficial.

## Revisión de las cuatro fases

| Requisito del enunciado | Evidencia concreta | Dictamen y alcance |
|---|---|---|
| I. Problema real y pregunta de optimización | Informe, pp. 2–4; `docs/FASE_I.md` | Cubierto. Diseño de variantes con dos propiedades de interés predichas. |
| I. Variables, objetivos y restricciones | 21 categorías, minimizar log MIC y menos log HC50, Hamming ≤2 | Cubierto. Dominio exacto de 76 210 secuencias; restricciones biológicas adicionales no se inventan. |
| I. Modelo del proceso explícito | `datos.py`, `modelos.py`, modelos y reserva | Cubierto. Dos modelos empíricos construidos por el equipo; no simulación mecanística. |
| II. Paradigma y justificación | Informe, p. 5; NSGA-II, Deb et al. | Cubierto. Multiobjetivo, elitismo y preservación de diversidad; no se exige implementar todo el sílabo. |
| II. Codificación, aptitud, operadores y selección | `evolucion.py`, pruebas y tabla del informe | Cubierto. Codificación categórica, cruce uniforme, mutación, reparación, torneo y supervivencia. |
| II. Biblioteca y justificación | pymoo 0.6.1.5, versión fijada y código verificado | Cubierto. Permite operadores propios y supervivencia multiobjetivo; se distingue documentación web posterior. |
| III. Implementación funcional | `run.py`, módulos y pruebas | Cubierto. Reejecución limpia verificada desde datos crudos incluidos. |
| III. Métrica clara y separada de aptitud | Vector de predicciones frente a HV externo | Cubierto. Referencia y escala comunes, denominador exhaustivo y comprobación independiente. |
| III. Al menos una sensibilidad | Población, mutación y cruce; cinco semillas por opción | Cubierto. Diseño secuencial; no se atribuyen interacciones ni óptimo global de parámetros. |
| IV. Promedio y mejor aptitud por generación | Figura 5 y `historia.csv` por corrida | Cubierto. Ambos objetivos en la población; sus mejores valores pueden ser individuos diferentes. |
| IV. Convergencia, presupuesto y eficiencia | Figuras 3–5; nueva discusión p. 9 | Cubierto con matices. Se cuantifican umbrales, mejoras tardías, estancamiento y costos excluidos. |
| IV. Efecto de sensibilidad sobre calidad y tiempo | Figura 2, tabla y discusión p. 9 | Cubierto. La elección prioriza HV, no tiempo; diferencias pequeñas siguen siendo descriptivas. |
| IV. Solución e interpretación práctica | Tres candidatos con semillas de procedencia; p. 10 | Cubierto como propuestas computacionales, sin atribuir efectos causales ni novedad demostrada. |
| IV. Comparación y pertinencia de métrica | Aleatoria con igual presupuesto e inicio; frente exhaustivo | Cubierto. HV compara conjuntos; su interpretación no es biológica ni conteo de soluciones. |
| IV. Logros, límites y futuro viable | Discusión, p. 11 | Cubierto. Se distinguen mejoras algorítmicas de validación del evaluador. |
| Entrega: informe de hasta unas 15 páginas y enlace | PDF de 13 páginas y URL real proporcionada por el equipo | Documento preparado. El contenido publicado debe corresponder con esta versión. |
| Entrega: exposición oral | `docs/GUION_DEFENSA.md`, recorrido de 30 minutos y preguntas | Preparada como guion; falta ensayo y realización por el equipo. |

La rúbrica distribuye 4 puntos a formulación, 4 a diseño, 4 a experimentación, 3 a reproducibilidad/documentación/responsabilidad, 2 a informe y 3 a defensa. Las evidencias anteriores permiten sustentar esos apartados; no corresponde asignarse automáticamente sus puntajes.

## Qué se comprobó, además de leer el informe

- Integridad previa: 4 577 archivos del manifiesto de entrega anterior coincidieron con sus SHA-256 antes de editar.
- Auditoría existente: 4 355 archivos crudos, particiones y métricas reservadas, dominio y HV, 55 corridas diferentes, presupuestos de 2 100 secuencias únicas y procedencia de candidatos.
- Reproducción completa sin caché: se creó una copia con código e instantánea cruda, sin datos procesados, modelos ni corridas. Se ejecutaron procesamiento, entrenamiento, sensibilidad, comparación y generación del informe.
- Comparación de esa reproducción: **122 tablas coincidieron**, incluidas las evaluaciones e historias de las 55 corridas. Tolerancia numérica absoluta 10⁻¹², relativa cero; los tiempos se excluyeron de la igualdad. La huella binaria de los modelos también coincidió.
- Recálculo separado de diferencia emparejada, Wilcoxon, selección de parámetros, recuperación del frente y costo de las corridas distintas.
- Comprobación posterior de HV con tres márgenes de referencia. Se mantuvo la ventaja media de NSGA-II sin reajustar el algoritmo ni sustituir la métrica primaria.
- Consulta de presencia exacta de los tres candidatos en los conjuntos procesados: ninguno aparece allí. Esto no certifica novedad en toda la base ni en la literatura.

Los resultados están en `results/auditoria_cierre.json`, `auditoria_convergencia.csv` y `auditoria_referencia_hv.csv`. El procedimiento está en `scripts/auditoria_cierre.py`; el registro de la nueva ejecución, en `results/reproduccion_desde_cero.log`. Se usaron la misma máquina y entorno: no se certificó reproducción en Colab, Linux u otro equipo. La revisión visual del PDF se documenta en `results/qa_documento.json`.

## Hallazgos que deben quedar claros al sustentar

**1. Calidad alta no significa frente completo en cada corrida.** NSGA-II obtiene 98,915% de HV medio frente a 87,313% de aleatoria. Cada corrida NSGA-II recupera entre dos y cinco de las seis secuencias exactas. La unión de diez corridas recupera las seis, con 21 000 consultas sumadas entre corridas. No se debe presentar esa unión como resultado de una sola corrida de 2 100 consultas.

**2. El presupuesto no certifica convergencia.** Ocho de diez semillas llegan a 99% de HV. Las semillas 401, 407 y 409 mejoran en el último lote; la 410 gana 5,983 puntos en las últimas 600 evaluaciones. La 406 se estanca en un valor subóptimo. La formulación correcta es “buena aproximación en la mayoría de las corridas”, sin garantía de completitud.

**3. La eficiencia tiene un costo previo.** La sensibilidad contiene 45 filas, pero 35 corridas distintas, que suman 42,93 s de búsqueda. Las 20 corridas finales agregan 19,87 s. Los descriptores se calcularon para todo el dominio. Consultar 2,76% del dominio por corrida no demuestra un ahorro total de 97,24%.

**4. La enumeración es viable y rápida en este caso.** Su predicción por lotes tomó 0,84 s; no se demuestra ventaja temporal de NSGA-II sobre enumeración. El dominio pequeño sirve como experimento controlado de búsqueda con presupuesto, no como prueba de escalabilidad o necesidad práctica de una heurística.

**5. La generalización biológica sigue abierta.** HC50 tiene capacidad explicativa limitada. Existen condiciones no armonizadas, exclusión de censura y publicaciones compartidas entre particiones: 133 en MIC y 9 en HC50, por título y año entre referencias de fichas. No se verificó una correspondencia artículo-medición para todas ellas. No hay validación externa ni ensayos de candidatos. El cociente HC50/MIC predicho no demuestra ventana terapéutica.

Estos límites restringen las conclusiones científicas; no son requisitos adicionales inventados del curso. El enunciado no obliga a realizar ensayos de laboratorio ni a ejecutar ColabFold. El material estructural histórico no se usa como validación de las variantes actuales.

## Cierre documental y siguientes pasos del equipo

Se completaron autoría y enlace, se fortaleció el análisis de fase IV y se preparó un guion de 30 minutos. La entrega reproducible incluye código, instantánea cruda, modelos, tablas, figuras, cuaderno e informe. El informe ampliado y el PDF corresponden al mismo experimento.

Antes de la presentación, el equipo necesita ensayar el guion, comprobar que todos entienden las limitaciones y realizar la entrega en el medio del curso. También debe confirmar el calendario: el enunciado indica semana 15 y el sílabo semana 16. No se deduce de ello una fecha calendario. Disponer del documento de fase I no demuestra que el parcial de semana 08 haya sido presentado: ese cumplimiento administrativo solo puede confirmarlo el equipo o el docente.

No hace falta añadir PSO, algoritmos genéticos paralelos u otros métodos solo porque aparecen en el sílabo. Añadirlos ampliaría el estudio, pero no repara ningún requisito obligatorio faltante identificado. La prioridad para una buena defensa es explicar con precisión las decisiones y las cifras que ya se pueden verificar.

## Referencias y alcance de las fuentes

- Deb et al. (2002), NSGA-II: https://doi.org/10.1109/4235.996017 . Respalda el algoritmo, no los resultados de este proyecto.
- Pirtskhalava et al. (2021), DBAASP v3: https://doi.org/10.1093/nar/gkaa991 . Respalda el recurso de datos, no que las variantes propuestas funcionen.
- Walsh et al. (2021), DOME: https://doi.org/10.1038/s41592-021-01205-4 . Respalda buenas prácticas de validación; no convierte esta reserva en validación externa.
- pymoo, NSGA-II: https://pymoo.org/algorithms/moo/nsga2.html . Se distingue la documentación actual de la versión 0.6.1.5 utilizada.
- pymoo, indicadores: https://pymoo.org/misc/indicators.html . Respalda la definición del hipervolumen; los valores concretos provienen de los CSV y su recálculo.

La bibliografía completa está en el informe. Las afirmaciones del experimento se verifican en los archivos entregados, no se atribuyen a autores externos.
