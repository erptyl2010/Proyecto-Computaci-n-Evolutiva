# Guion de sustentación: 30 minutos

Equipo: Ing. MIGUEL HERRERA, WILLY ROBERTO; Ing. PARIHUANA TRAVEZAÑO, OMAR EDUARDO; Ing. TORRES CRISTOBAL, WILBER.

Repositorio: https://github.com/erptyl2010/Proyecto-Computaci-n-Evolutiva

El equipo indicó una duración de 30 minutos. La distribución siguiente es una propuesta de ensayo; no presupone que el docente incluya las preguntas dentro de ese tiempo. Se plantean 25 minutos de exposición y cinco de margen o preguntas. Si las preguntas van fuera, usar ese margen para la demostración. El reparto por integrante es una sugerencia organizativa.

## Recorrido de la exposición

| Tiempo | Responsable sugerido | Qué explicar y mostrar |
|---|---|---|
| 0–2 min | Willy | Pregunta, objetivo del curso y alcance computacional. Abrir con el resumen del informe. |
| 2–5 min | Willy | Qué significan MIC y HC50, variables, dos objetivos y restricción de dos cambios. Mostrar la fórmula de 76 210 posibilidades. |
| 5–8 min | Willy | Datos DBAASP, filtros, unidades, agrupación y reserva. Mostrar la tabla de datos y validación. |
| 8–10 min | Omar | Separar predictor y optimizador. Explicar la debilidad de HC50 antes de presentar resultados. |
| 10–14 min | Omar | Cromosoma, selección, cruce, mutación, reparación, duplicados y archivo. Mostrar una corrida real en el código o cuaderno. |
| 14–17 min | Omar | Presupuesto justo, búsqueda aleatoria, hipervolumen, semillas de ajuste y finales. Mostrar sensibilidad. |
| 17–20 min | Wilber | Comparación final: calidad y tiempo, diferencia emparejada y alcance del intervalo. |
| 20–23 min | Wilber | Promedio/mejor por generación, presupuesto insuficiente para recuperar todo en cada corrida, candidatos encontrados. |
| 23–25 min | Wilber | Qué se concluye, qué no queda demostrado y siguientes mejoras. Mostrar auditoría y enlace del repositorio. |
| 25–30 min | Todos | Preguntas o demostración breve; evitar ejecutar otra descarga o un entrenamiento completo durante la defensa. |

## Apertura que pueden adaptar a su forma de hablar

“Queríamos estudiar cómo un algoritmo evolutivo busca soluciones cuando mejorar una propiedad puede perjudicar otra. Trabajamos con secuencias de péptidos y dos propiedades estimadas a partir de datos públicos. Nuestra pregunta no es si ya fabricamos un péptido mejor, sino qué tan bien NSGA-II aprovecha un número limitado de evaluaciones frente a una búsqueda aleatoria. Para comprobarlo usamos un dominio pequeño cuyo frente de referencia podemos calcular por completo.”

## Las ideas que cada integrante debe comprender

**MIC** es la concentración mínima inhibitoria registrada frente a la cepa elegida. En este proyecto se desea un valor menor. **HC50** es la concentración registrada para 50% de hemólisis en eritrocitos humanos; se desea un valor mayor. HC50 no mide toda la toxicidad ni establece por sí sola seguridad clínica. El informe cita DBAASP para el origen y significado de las mediciones.

**Random Forest y NSGA-II cumplen funciones distintas.** Los bosques aprenden a estimar dos propiedades desde los descriptores. NSGA-II consulta esas estimaciones y busca compromisos. Mejorar la búsqueda no elimina el error del modelo. Aquí el R² reservado de HC50 es aproximadamente 0,059: hay poca capacidad para explicar su variación. Esta es una limitación central, no un detalle para esconder.

**Dominancia y frente.** Una solución domina a otra si no empeora ningún objetivo y mejora al menos uno. El frente reúne soluciones no dominadas dentro del conjunto considerado. El frente observado de una corrida puede contener puntos dominados por secuencias que esa corrida nunca vio. No confundir su tamaño con el número de soluciones del frente exhaustivo recuperadas.

**Por qué el signo negativo.** La biblioteca trabaja con minimización. Por eso minimizamos `(log10 MIC, -log10 HC50)`. Así, aumentar HC50 mejora el segundo objetivo. El logaritmo ayuda a trabajar en una escala relativa de concentraciones; conserva el orden, pero cambia el valor geométrico del hipervolumen.

**Por qué la codificación es categórica.** Cada posición es uno de 20 símbolos. Que un símbolo tenga índice mayor no significa que sea físicamente mayor o más parecido a otro. El cruce intercambia posiciones; la mutación elige otro aminoácido. La reparación vuelve a la referencia las posiciones excedentes hasta dejar como máximo dos sustituciones.

**Qué se compara.** Diez parejas de semillas, 2 100 secuencias distintas por método y el mismo inicio de 60 secuencias en cada pareja. El presupuesto es por secuencia, con dos predicciones por secuencia. Las semillas finales 401–410 son distintas de las de sensibilidad 201–205.

**Qué significa 98,915%.** Es el HV medio de NSGA-II como porcentaje del HV exhaustivo, en un espacio de predicciones y con una referencia declarada. Aleatoria obtuvo 87,313%. La diferencia es 11,602 puntos porcentuales. No es porcentaje de exactitud, de péptidos eficaces ni de soluciones del frente recuperadas. La definición del indicador se respalda en https://pymoo.org/misc/indicators.html .

## Preguntas difíciles y respuestas defendibles

**1. ¿Por qué NSGA-II?**

Porque se buscan dos objetivos sin imponer un peso único y se quiere conservar un conjunto de compromisos. El elitismo y la diversidad por hacinamiento son adecuados para ese propósito. Es una decisión justificada, no una prueba de que sea el mejor algoritmo posible. El fundamento es Deb et al. (2002), https://doi.org/10.1109/4235.996017 .

**2. ¿Por qué no enumerar todo?**

En este dominio sí es viable. Lo enumeramos para conocer la referencia exacta del modelo y medir la calidad lograda con presupuesto reducido. No afirmamos que NSGA-II sea necesario ni más rápido aquí: la predicción exhaustiva en lote tomó 0,84 s y la búsqueda NSGA-II promedió 1,23 s, con operaciones cronometradas diferentes. La contribución es evaluar el comportamiento del algoritmo; escalabilidad requeriría otro experimento.

**3. ¿Ahorramos realmente 97% del trabajo?**

No puede decirse así. Una corrida consulta el 2,7555% de las secuencias, pero se precomputaron descriptores para todo el dominio y hubo entrenamiento, enumeración y ajuste. La afirmación correcta es sobre consultas únicas durante cada búsqueda, no sobre ahorro integral de cómputo.

**4. ¿Cuántas generaciones usamos y bastaron?**

35 registros incluyendo la población inicial: 34 generaciones de descendencia con población 60. Ocho de diez corridas alcanzan 99% de HV; ninguna recupera sola las seis soluciones exactas. Tres aún mejoran en el último lote. El presupuesto funciona bien para aproximar, sin garantizar completitud ni convergencia.

**5. ¿Entonces por qué la figura muestra todo el frente?**

Porque muestra la unión de diez corridas NSGA-II. Esa unión recupera las seis secuencias, con 21 000 consultas sumadas entre corridas. Una corrida individual encuentra entre dos y cinco. Presentar la unión como resultado de solo 2 100 consultas sería incorrecto.

**6. ¿Hubo convergencia prematura?**

La semilla 406 se estanca en 97,363% desde 1 440 evaluaciones; es una meseta subóptima observada. No aislamos si se debe a pérdida de diversidad, operadores u otra causa. Tampoco demostramos que aumentar generaciones garantice llegar al frente completo.

**7. ¿Cómo elegimos parámetros?**

Variamos población, mutación y cruce, un factor a la vez de forma secuencial y con cinco semillas por opción. Elegimos mayor HV medio. Resultado: población 60, mutación 0,5/21 por posición, cruce 0,7. La tabla tiene 45 filas pero solo 35 corridas distintas por reutilización de configuraciones iguales. No estima interacciones ni asegura un óptimo global de parámetros.

**8. ¿Cuánto costó el ajuste?**

35 corridas distintas y 42,93 s sumados de búsqueda. La comparación final agrega 20 corridas y 19,87 s. Son tiempos de estas ejecuciones; excluyen descarga, preparación, entrenamiento y escritura. El ajuste también forma parte del costo del estudio.

**9. ¿Qué garantiza la prueba estadística?**

La diferencia emparejada tiene IC bootstrap de 95% [10,029; 12,949] puntos de HV y Wilcoxon unilateral p=0,000976563. Describe variación entre semillas en un único paisaje congelado. No mide incertidumbre biológica, probabilidad de que la hipótesis sea cierta ni superioridad universal.

**10. ¿El conjunto reservado es externo?**

No. Es una reserva interna por grupos de secuencias conectadas mediante Hamming ≤2 y misma longitud. Persisten publicaciones compartidas y condiciones heterogéneas. No debemos decir que validamos otros laboratorios o familias independientes. DOME respalda la importancia de separar selección y evaluación: https://doi.org/10.1038/s41592-021-01205-4 .

**11. ¿Los tres candidatos son nuevos?**

No tienen coincidencia exacta en los dos conjuntos procesados de esta cohorte. Eso no prueba novedad en toda DBAASP ni en la literatura. Son candidatos computacionales para este experimento; no se realizaron ensayos de ellos.

**12. ¿Qué significa que HC50 sea mayor?**

Dentro del modelo, se requiere mayor concentración para alcanzar el nivel de hemólisis registrado. Es una dirección deseada para ese ensayo específico. No autoriza a llamar al candidato “seguro”, “no tóxico” o “terapéutico”.

**13. ¿Por qué se excluyeron valores como >100?**

Porque no son concentraciones puntuales: sabemos una desigualdad, no el valor exacto. El protocolo actual los excluye de forma trazable. Esa decisión puede sesgar los datos; tratar explícitamente la censura sería una mejora futura, no algo ya realizado.

**14. ¿Qué aporta ColabFold?**

El anexo previo es exploratorio y no valida las variantes de la reconstrucción. No se reejecutó. No lo presentamos como prueba de actividad, hemólisis o estabilidad.

**15. ¿Qué demostramos al repetir todo desde cero?**

Que, en la misma máquina y entorno, la instantánea cruda y el código producen los mismos datos, modelos, evaluaciones y métricas dentro de la tolerancia declarada. No es una repetición independiente en otro laboratorio o sistema. Los tiempos pueden variar.

**16. ¿Qué mejoraríamos primero?**

La validez del evaluador: armonizar condiciones, separar por publicación/familia, tratar censura e incertidumbre. En el estudio evolutivo: varias referencias y radios, interacción de parámetros y otros operadores. Para afirmar mejoras biológicas se necesitaría validación independiente de los candidatos.

## Demostración breve y comprobable

1. Abrir `notebooks/PROYECTO_CE_REPRODUCIBLE.ipynb` o el informe con las tablas reales.
2. Mostrar `results/optimization/comparacion_final.csv` y ubicar una semilla.
3. Abrir su carpeta `results/optimization/runs/<run_id>/`: enseñar `evaluaciones.csv`, `historia.csv` y `complete.json`.
4. Relacionar esa corrida con la gráfica de convergencia. Mostrar que el presupuesto final es 2 100 y la historia termina en generación 35.
5. Mostrar `results/auditoria_cierre.json` y `docs/AUDITORIA_FINAL_CIERRE.md`. Explicar cómo se verificaron las cifras sin pretender leer el JSON completo ante el jurado.

## Cierre sugerido

“En las condiciones evaluadas, NSGA-II aprovechó mejor las consultas que el muestreo aleatorio para aproximar el frente del modelo. Lo pudimos comprobar porque controlamos presupuesto, semillas y referencia. El aprendizaje principal es que una buena búsqueda y un buen predictor son problemas distintos: resolvimos y auditamos el primero en este caso, mientras que la validez biológica del segundo limita lo que podemos afirmar sobre los candidatos.”

Antes de exponer, cada integrante debería responder sin leer las preguntas 2, 3, 4, 5, 9 y 10. Esas preguntas comprueban que el equipo entiende el experimento y sus límites.
