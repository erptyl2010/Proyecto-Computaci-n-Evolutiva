# Fase I: definición y modelado del problema

**Título propuesto:** Evaluación de NSGA-II para la búsqueda multiobjetivo de variantes de péptidos mediante propiedades predichas.

**Pregunta:** ¿Con qué calidad aproxima NSGA-II el frente de compromiso entre MIC y HC50 predichas, frente a una búsqueda aleatoria con igual número de secuencias evaluadas?

Queremos estudiar un problema real de diseño con objetivos que pueden entrar en conflicto. Una MIC menor es la dirección deseada para la actividad frente a la cepa elegida; una HC50 mayor es la dirección deseada respecto de hemólisis. La tarea del curso consiste en formular y evaluar la búsqueda. La evidencia que obtendremos es computacional: no se realizarán ensayos de las variantes.

## Variables, objetivos y restricciones

La referencia es `GFMKYIKPLIPHAVKAIKKLI`, de 21 posiciones. Su ficha DBAASP 2514 está conservada en `data/raw/peptides/2514.json`. Las 21 variables son categóricas, cada una con un aminoácido de `ACDEFGHIKLMNPQRSTVWY`. Se mantiene longitud 21 y se permiten hasta dos sustituciones. Por ello:

`Ω = {x ∈ A²¹ : d_H(x,s₀) ≤ 2}`

`|Ω| = 1 + 21×19 + 210×19² = 76 210`.

El vector a minimizar es `F(x) = (m_MIC(x), −m_HC50(x))`. Los regresores `m` devuelven log10 de la concentración expresada en µM. La solución es un conjunto de compromisos no dominados; no se asignan pesos arbitrarios para colapsarlo en una sola propiedad. Se mantiene el contexto terminal sin anotación de la referencia. Esto no equivale a demostrar que los grupos terminales sean libres.

## Modelo del proceso

Se construyen dos modelos supervisados Random Forest con mediciones públicas de DBAASP. MIC corresponde exclusivamente a *Escherichia coli* ATCC 25922; HC50, a eritrocitos humanos. La cohorte se define antes del entrenamiento y las exclusiones conservan un motivo auditable. Se distinguen secuencia y anotaciones terminales al agrupar mediciones.

El flujo es: API y registros originales → filtrado trazable → componentes de secuencias similares → reserva interna → selección de representación solo en entrenamiento → modelos congelados → búsqueda → comparación y análisis. El modelo del proceso es empírico; no es un simulador mecanístico de membranas.

## Método de resolución y evaluación

Se usa NSGA-II con cruce uniforme, mutación categórica, reparación a un máximo de dos cambios, torneo y supervivencia elitista. Se eliminan secuencias repetidas durante toda la corrida. La hipótesis es que su selección permite aprovechar mejor 2 100 consultas que una muestra uniforme sin reemplazo.

La aptitud es el vector de predicciones. La calidad del conjunto se evalúa mediante hipervolumen relativo a la referencia exhaustiva del mismo paisaje. Se examinan sensibilidad, promedio y mejor aptitud por generación, tiempo y variabilidad entre semillas. Se reservan semillas finales distintas de las empleadas en ajuste de parámetros.

## Alcance y límites anticipados

El barrido completo es viable por el tamaño del dominio. Su propósito es validar la métrica y estudiar presupuestos menores; no demuestra que el método escale automáticamente. El error del modelo, la heterogeneidad de mediciones y la exclusión de valores censurados limitan la interpretación. La actividad antimicrobiana, hemólisis y seguridad reales de las variantes seguirán sin estar establecidas.

## Referencias de apoyo

- Pirtskhalava et al. DBAASP v3 (2021). https://doi.org/10.1093/nar/gkaa991
- Deb et al. NSGA-II (2002). https://doi.org/10.1109/4235.996017
- Breiman. Random Forests (2001). https://doi.org/10.1023/A:1010933404324
- Walsh et al. DOME (2021). https://doi.org/10.1038/s41592-021-01205-4

Este documento responde a la fase I del enunciado local: problema, formulación y modelo del proceso. Los resultados de la reconstrucción se documentan por separado en `INFORME_RESULTADOS.md`.
