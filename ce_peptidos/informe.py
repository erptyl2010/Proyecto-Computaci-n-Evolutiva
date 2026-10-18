"""Figuras y narración alimentadas exclusivamente por los CSV de esta ejecución."""
import json
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from .evolucion import BASE, nondominated

def md_table(df, decimals=4):
    def fmt(v):
        if isinstance(v, (float, np.floating)):
            return f"{v:.{decimals}f}"
        return str(v).replace("|", "/")
    return "| " + " | ".join(map(str, df.columns)) + " |\n| " + " | ".join(["---"]*len(df.columns)) + " |\n" + "\n".join("| " + " | ".join(fmt(v) for v in row) + " |" for row in df.itertuples(index=False, name=None))

def report(root):
    out = root / "results" / "optimization"
    figures = root / "results" / "figures"
    docs = root / "docs"
    figures.mkdir(exist_ok=True)
    docs.mkdir(exist_ok=True)
    read_json = lambda p: json.loads(p.read_text(encoding="utf-8"))
    stats = read_json(out / "estadisticas.json")
    ref = read_json(out / "referencia.json")
    data = read_json(root / "data" / "processed" / "resumen_datos.json")
    raw = read_json(root / "data" / "raw" / "manifest.json")
    metrics = pd.read_csv(root / "results" / "models" / "validacion_reservada.csv")
    comparison = pd.read_csv(out / "comparacion_final.csv")
    sensitivity = pd.read_csv(out / "sensibilidad.csv")
    exhaustive = pd.read_csv(out / "dominio_exhaustivo.csv")
    histories, archives = [], []
    for r in comparison.itertuples():
        folder = out / "runs" / r.run_id
        h = pd.read_csv(folder / "historia.csv")
        histories.append(h.assign(method=r.method, seed=r.seed))
        if r.method == "NSGA-II":
            a = pd.read_csv(folder / "evaluaciones.csv")
            archives.append(a.assign(seed=r.seed, run_id=r.run_id))
    history = pd.concat(histories, ignore_index=True)
    history.to_csv(out / "historia_comparacion.csv", index=False)
    found = pd.concat(archives, ignore_index=True)
    provenance = found.groupby("sequence").seed.apply(lambda x: ",".join(map(str, sorted(set(x)))))
    found = found.drop_duplicates("sequence").copy()
    found["seeds_found"] = found.sequence.map(provenance)
    F = found[["pred_log10_MIC", "pred_log10_HC50"]].to_numpy() * [1, -1]
    front = found[nondominated(F)].copy()
    front["pred_ratio_HC50_MIC"] = front.pred_HC50_uM / front.pred_MIC_uM
    front["mutations_1based"] = front.sequence.map(lambda s: ";".join(f"{b}{i+1}{a}" for i, (a,b) in enumerate(zip(s,BASE)) if a != b) or "Original")
    front.to_csv(out / "frente_encontrado_nsga.csv", index=False)
    selections = [("Menor MIC predicha", front.pred_MIC_uM.idxmin()),
                  ("Mayor HC50 predicha", front.pred_HC50_uM.idxmax()),
                  ("Mayor cociente predicho", front.pred_ratio_HC50_MIC.idxmax())]
    candidates = pd.DataFrame([{"criterio": label, **front.loc[index].to_dict()} for label, index in selections])
    candidates.to_csv(out / "candidatos_encontrados.csv", index=False)
    plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False,
                         "figure.dpi": 130, "savefig.dpi": 170})
    colors = {"NSGA-II": "#087e8b", "Aleatoria": "#d16b32"}
    fig, ax = plt.subplots(figsize=(8.3, 4.3), layout="constrained")
    for method, g in history.groupby("method"):
        s = g.groupby("n_evaluations").hv_percent.agg(["mean", "std"])
        ax.plot(s.index, s["mean"], label=method, color=colors[method])
        ax.fill_between(s.index, s["mean"]-s["std"], s["mean"]+s["std"], alpha=.15, color=colors[method])
    ax.axhline(100, color="gray", linestyle=":", label="Referencia exhaustiva")
    ax.set(xlabel="Secuencias únicas evaluadas", ylabel="Hipervolumen relativo (%)", title="Calidad de la búsqueda · media ± desviación, 10 semillas")
    ax.legend(loc="lower right")
    fig.savefig(figures / "convergencia_hipervolumen.png")
    plt.close(fig)
    fig, axes = plt.subplots(1, 2, figsize=(10.4, 4), layout="constrained")
    evolutionary = history[history.method == "NSGA-II"]
    for ax, objective, direction in zip(axes, ["MIC", "HC50"], ["menor es mejor", "mayor es mejor"]):
        for statistic, label, color in [("mean", "Promedio poblacional", "#087e8b"), ("best", "Mejor de la población", "#d16b32")]:
            g = evolutionary.groupby("generation_or_batch")[f"{statistic}_log10_{objective}"].agg(["mean", "std"])
            ax.plot(g.index, g["mean"], label=label, color=color)
            ax.fill_between(g.index, g["mean"]-g["std"], g["mean"]+g["std"], color=color, alpha=.12)
        ax.set(xlabel="Generación (1 = población inicial)", ylabel=f"log10 {objective} predicha [µM]", title=f"{objective}: {direction}")
        ax.legend(fontsize=8)
    fig.savefig(figures / "aptitud_promedio_y_mejor.png")
    plt.close(fig)
    fig, axes = plt.subplots(1, 2, figsize=(9.5, 4), layout="constrained")
    for ax, field, title in zip(axes, ["hv_percent", "seconds"], ["Hipervolumen relativo (%)", "Tiempo de búsqueda (s)"]):
        values = [comparison[comparison.method == m][field] for m in ["NSGA-II", "Aleatoria"]]
        ax.boxplot(values, tick_labels=["NSGA-II", "Aleatoria"], showfliers=False)
        for i, v in enumerate(values, 1):
            ax.scatter(np.linspace(i-.09, i+.09, len(v)), v, s=22, color=list(colors.values())[i-1])
        ax.set_ylabel(title)
    fig.suptitle("Comparación final · 2 100 secuencias únicas por corrida")
    fig.savefig(figures / "comparacion_calidad_tiempo.png")
    plt.close(fig)
    fig, axes = plt.subplots(2, 3, figsize=(11.8, 6), layout="constrained")
    for col, (parameter, label) in enumerate([("population", "Población"), ("mutation", "Mutación por posición"), ("crossover", "Probabilidad de cruce")]):
        g = sensitivity[sensitivity.sensitivity_parameter == parameter]
        for row, field, ylabel in [(0, "hv_percent", "HV relativo (%)"), (1, "seconds", "Tiempo (s)")]:
            s = g.groupby("sensitivity_value")[field].agg(["mean", "std"])
            ax = axes[row, col]
            ax.errorbar(np.arange(len(s)), s["mean"], yerr=s["std"], fmt="o-", color="#087e8b", capsize=4)
            ax.set_xticks(np.arange(len(s)), [f"{v:.4g}" for v in s.index])
            ax.set(xlabel=label, ylabel=ylabel)
    fig.suptitle("Sensibilidad secuencial · media ± desviación, cinco semillas por opción")
    fig.savefig(figures / "sensibilidad.png")
    plt.close(fig)
    fig, ax = plt.subplots(figsize=(7.5, 4.8), layout="constrained")
    ax.scatter(exhaustive.pred_log10_MIC, exhaustive.pred_log10_HC50, s=3, alpha=.08, color="#728791", rasterized=True, label="Dominio completo predicho")
    pf = exhaustive[exhaustive.nondominated]
    ax.scatter(pf.pred_log10_MIC, pf.pred_log10_HC50, s=25, facecolors="none", edgecolors="#d16b32", label="Frente exhaustivo")
    ax.scatter(front.pred_log10_MIC, front.pred_log10_HC50, s=12, color="#087e8b", label="Frente de unión de 10 corridas NSGA-II")
    ax.scatter([ref["base_prediction"]["pred_log10_MIC"]], [ref["base_prediction"]["pred_log10_HC50"]], marker="*", s=110, color="black", label="Secuencia inicial")
    ax.set(xlabel="log10 MIC predicha [µM] · minimizar", ylabel="log10 HC50 predicha [µM] · maximizar", title="Compromisos entre las dos propiedades predichas")
    ax.legend(fontsize=8)
    fig.savefig(figures / "frentes_pareto.png")
    plt.close(fig)
    pred = pd.read_csv(root / "results" / "models" / "predicciones_reservadas.csv")
    fig, axes = plt.subplots(1, 2, figsize=(9, 4), layout="constrained")
    for ax, (objective, g) in zip(axes, pred.groupby("objective", sort=False)):
        ax.scatter(g.experimental_log10, g.predicted_log10, color="#087e8b", s=13, alpha=.6)
        limits = [min(g.experimental_log10.min(), g.predicted_log10.min()), max(g.experimental_log10.max(), g.predicted_log10.max())]
        ax.plot(limits, limits, "--", color="gray")
        ax.set(xlabel="log10 concentración observada [µM]", ylabel="log10 concentración predicha [µM]", title=f"{objective} · prueba interna reservada")
    fig.savefig(figures / "validacion_predictores.png")
    plt.close(fig)
    hit_rows = []
    for (method, seed), g in history.groupby(["method", "seed"]):
        for threshold in [99., 99.9]:
            reached = g[g.hv_percent >= threshold]
            hit_rows.append({"method": method, "seed": seed, "threshold_percent": threshold,
                             "reached": len(reached)>0, "first_recorded_evaluations": None if reached.empty else int(reached.iloc[0].n_evaluations)})
    hits = pd.DataFrame(hit_rows)
    hits.to_csv(out / "evaluaciones_hasta_umbral.csv", index=False)
    hit_summary = hits.groupby(["method", "threshold_percent"]).agg(successes=("reached", "sum"),
                    median_evaluations_among_successes=("first_recorded_evaluations", "median")).reset_index()
    summary_table = comparison.groupby("method").agg(HV_media=("hv_percent", "mean"), HV_DE=("hv_percent", "std"),
                                                     segundos_media=("seconds", "mean")).reset_index()
    sensitivity_summary = sensitivity.groupby(["sensitivity_parameter", "sensitivity_value"]).agg(
        HV_media=("hv_percent", "mean"), HV_DE=("hv_percent", "std"), segundos_media=("seconds", "mean")).reset_index()
    def sensitivity_value(parameter, value, column):
        selected = sensitivity_summary[(sensitivity_summary.sensitivity_parameter == parameter) &
            np.isclose(sensitivity_summary.sensitivity_value, value)]
        return float(selected.iloc[0][column])
    tuning_cost = float(sensitivity.drop_duplicates('run_id').seconds.sum())
    final_cost = float(comparison.seconds.sum())
    delta = stats["paired_mean_difference_percentage_points"]
    ci = stats["paired_bootstrap_95_percent_interval"]
    conclusion = ("NSGA-II alcanzó mayor hipervolumen medio que la búsqueda aleatoria" if delta > 0 else
                  "NSGA-II no alcanzó mayor hipervolumen medio que la búsqueda aleatoria")
    significance = ("La prueba unilateral favorece esa hipótesis al umbral descriptivo de 0,05." if stats["wilcoxon_pvalue"] < .05 else
                    "La prueba unilateral no permite sostener superioridad al umbral de 0,05.")
    model_notes = []
    for r in metrics.itertuples():
        better = "menor" if r.mae_log10 < r.mae_baseline_log10 else "mayor o igual"
        model_notes.append(f"En {r.objective}, el MAE del bosque es {better} que el del predictor constante de media del entrenamiento. Un error absoluto medio de {r.mae_log10:.3f} en log10 equivale a un factor geométrico medio de discrepancia de {10**r.mae_log10:.2f}; no es un intervalo de confianza individual.")
        if r.r2_log10 < .1:
            model_notes.append(f"El R² reservado de {r.objective} ({r.r2_log10:.3f}) muestra una capacidad limitada para explicar la variación observada, aun cuando mejore el MAE del referente constante.")
    config = stats["selected_configuration"]
    full_text = f'''# Búsqueda evolutiva multiobjetivo de variantes de péptidos mediante propiedades predichas

Informe técnico de desarrollo · Computación Evolutiva · Maestría en Inteligencia Artificial

Versión reconstruida con datos públicos. Las tablas y figuras se generan desde los archivos de esta ejecución, sin trasladar las cifras históricas de los notebooks originales. Este documento es la versión ampliada en Markdown. El informe académico se maquetó en 13 páginas en `output/pdf/Informe_CE_Peptidos.pdf`; Integrantes: Ing. MIGUEL HERRERA, WILLY ROBERTO; Ing. PARIHUANA TRAVEZAÑO, OMAR EDUARDO; Ing. TORRES CRISTOBAL, WILBER. Repositorio: https://github.com/erptyl2010/Proyecto-Computaci-n-Evolutiva .

## Resumen

Estudiamos si NSGA-II aprovecha un presupuesto reducido para aproximar el frente de compromiso entre menor MIC y mayor HC50 predichas. La referencia es un dominio de {ref['n_sequences']:,} secuencias de 21 aminoácidos, con hasta dos sustituciones. Reconstruimos una cohorte de {data['cohort_api_ids']:,} fichas de DBAASP y entrenamos dos Random Forest con una partición interna reservada. En diez semillas finales, cada método evaluó 2 100 secuencias únicas. {conclusion}, con una diferencia emparejada de {delta:.4f} puntos porcentuales. Esta conclusión describe la búsqueda sobre modelos congelados. No demuestra actividad, seguridad ni superioridad biológica de los candidatos.

## I. Definición y modelado del problema

### Motivación y pregunta

El interés está en el compromiso entre actividad antimicrobiana y hemólisis: una MIC menor y una HC50 mayor son las direcciones deseadas para estos dos indicadores específicos. HC50 no resume todas las formas de toxicidad. DBAASP recoge secuencias, modificaciones y mediciones biológicas con sus referencias, lo que permite construir un problema verificable, aunque heterogéneo [1]. Boone et al. ofrecen un antecedente de integración de algoritmos genéticos y aprendizaje automático para péptidos; su sistema y organismo objetivo son distintos a los de este trabajo [2]. No presentamos esta combinación general como una innovación inédita.

La pregunta es: **con 2 100 consultas únicas a dos predictores congelados, ¿qué calidad de frente alcanza NSGA-II frente a una muestra uniforme sin reemplazo?** La hipótesis algorítmica es que NSGA-II obtiene mayor hipervolumen medio. El alcance es una evaluación computacional de presupuesto limitado sobre un solo paisaje predicho.

### Formulación

La referencia es `s₀ = {BASE}`. Cada variable `xᵢ` representa un aminoácido del alfabeto de 20 símbolos canónicos, para `i = 1,…,21`. Se impone `d_H(x,s₀) ≤ 2`, sin inserciones ni eliminaciones. Los extremos permanecen en el contexto de anotación ausente usado por los modelos; no se infiere su estado químico desde un campo nulo. La secuencia está presente en la ficha DBAASP 2514, conservada como JSON en esta reconstrucción.

El problema minimiza el vector `F(x) = (m_MIC(x), −m_HC50(x))`, donde cada `m` predice `log10(concentración/µM)`. Es una optimización multiobjetivo: puede haber soluciones de compromiso y no se presupone un único ganador. Las concentraciones mostradas son `10^m`, no mediciones nuevas. El cociente predicho HC50/MIC se usa únicamente para describir una elección posible del frente; no sustituye los dos objetivos.

La cardinalidad exacta es `1 + 21×19 + C(21,2)×19² = 76 210`. No hay restricciones implementadas de síntesis, estabilidad, agregación o toxicidad general. Esta delimitación permite una referencia exhaustiva del modelo y hace posible cuantificar cuánto se obtiene con el {2100/76210*100:.4f}% del dominio por corrida. No demuestra escalabilidad a secuencias o radios mayores.

### Modelo del proceso y datos

El modelo del proceso se construye con dos regresores Random Forest [3], uno por propiedad. Los datos proceden de la API oficial de DBAASP [1,4]. La instantánea se inició en `{raw['started_utc']}` y terminó en `{raw['finished_utc']}`. El filtro solicita monómeros de longitud 10–30, sin aminoácidos inusuales y con registro frente a *E. coli* ATCC 25922. Se descargan todas las fichas devueltas; HC50 se extrae de esa misma cohorte. Esto limita su representatividad frente a un muestreo independiente de todos los datos de hemólisis.

Se conservan únicamente MIC de esa cepa y 50% de hemólisis en eritrocitos humanos, expresadas explícitamente en µM, con concentración puntual, positiva y finita. Se excluyen rangos, desigualdades, valores ±, otras unidades y química no representada. Se usa `concentration`, no el campo `activity`. Las anotaciones terminales ACT y AMD se diferencian de la ausencia de anotación. Se agrupan mediciones de una misma secuencia y anotaciones terminales mediante la mediana. Esto reduce redundancia, pero no armoniza medios, inóculo ni demás condiciones.

Tras procesar {data['measurements_audited']:,} mediciones de las fichas, se aceptaron {data['accepted_measurements']['MIC']:,} mediciones MIC y {data['accepted_measurements']['HC50']:,} HC50, correspondientes a {data['chemical_entities']['MIC']:,} y {data['chemical_entities']['HC50']:,} entidades con secuencia y anotaciones terminales. Cada fila conserva ID, concentración original, condiciones disponibles, referencia original, URL y SHA-256. Los artículos de cada ficha se conservan sin atribuir una medición a un artículo concreto cuando esa correspondencia no se ha verificado.

Se comparan una representación de composición de 27 atributos y otra de 42 que agrega descriptores de orden. Ambas incluyen longitud, frecuencias, carga aproximada, fracción hidrofóbica e indicadores terminales. El orden se resume mediante segmentos terminales, posición y productos a desfases 1 y 2 de hidropatía, carga y aromaticidad. La escala de hidropatía se basa en Kyte–Doolittle [5]. La carga es una aproximación definida en el código, no un cálculo de ionización dependiente del pH; los productos a desfase no son autocorrelaciones normalizadas.

Las secuencias de igual longitud conectadas por Hamming ≤ 2 forman componentes transitivas. Se reserva el 20% de grupos con semilla 20261017; no necesariamente el 20% de filas. En entrenamiento, cinco pliegues agrupados seleccionan la representación por MAE medio. Cada bosque tiene 120 árboles, mínimo de tres observaciones por hoja y un hilo. El modelo final se ajusta solo con entrenamiento. El conjunto reservado se evalúa después de elegir la representación, conforme al principio de separar selección y evaluación [6]. La agrupación no elimina toda homología y no separa por publicación.

{md_table(metrics[[c for c in ['objective','representation','n_features','n_train','n_holdout','mae_log10','r2_log10','mae_baseline_log10','shared_publications_train_holdout'] if c in metrics.columns]])}

{' '.join(model_notes)} Se auditan publicaciones compartidas por título normalizado y año entre las referencias de las fichas. Los IDs internos de artículos no son suficientes para identificar publicaciones únicas. La coincidencia, cuando existe, advierte dependencia residual posible, sin demostrar que todas esas referencias correspondan a las mediciones incluidas. Esta prueba es interna y no una validación externa. Los valores reservados no se usan para elegir candidatos del optimizador.

![Predicción y observación en la reserva interna](../results/figures/validacion_predictores.png)

## II. Diseño del algoritmo evolutivo

NSGA-II se elige por su selección elitista y tratamiento explícito de objetivos en conflicto [7]. La implementación emplea pymoo 0.6.1.5 [8], elegida por permitir integrar operadores categóricos y reparación propios con su supervivencia multiobjetivo. Se fijan las versiones y se contrasta el hipervolumen con un cálculo independiente. El cromosoma tiene 21 enteros categóricos con correspondencia uno a uno con la secuencia. El cruce uniforme intercambia posiciones entre dos progenitores. La mutación sustituye cada posición elegida por uno de los otros 19 aminoácidos con igual probabilidad. Si quedan más de dos cambios, se restauran posiciones al azar a la referencia hasta cumplir la restricción. Esta reparación introduce un sesgo que debe reconocerse; se estudia la variante implementada, no cualquier posible NSGA-II.

La selección parental de esta versión de pymoo usa torneo binario por dominancia y distancia de hacinamiento; la supervivencia usa ordenamiento no dominado y hacinamiento. No se aplican distancias continuas al índice alfabético de los aminoácidos. Se impide evaluar nuevamente una secuencia ya consultada durante la corrida, además de eliminar duplicados dentro del lote. La población inicial de cada comparación es compartida entre métodos. El archivo de soluciones conserva todo lo evaluado y calcula el frente observado, aunque algunos puntos ya no pertenezcan a la población superviviente.

## III. Implementación y experimentación

La aptitud es el vector de dos predicciones. La métrica externa es el hipervolumen en `(log10 MIC, −log10 HC50)` respecto de un punto común: `{ref['reference_point']}`. Se fija tomando el máximo de cada objetivo en el dominio más 0,15. Se informa `100 × HV(archivo)/HV(frente exhaustivo)`. El punto usa conocimiento del dominio como referencia de evaluación; los operadores no consultan los valores exhaustivos. El hipervolumen mide el área dominada respecto de la referencia y se maximiza [9]. Es pertinente para comparar conjuntos de compromisos bajo una misma escala; no es precisión predictiva, probabilidad de éxito ni fracción de secuencias del frente recuperadas. El logaritmo conserva la dominancia pero modifica la geometría del área.

El barrido exhaustivo encuentra {ref['n_front_sequences']:,} secuencias no dominadas que representan {ref['n_front_objective_pairs']:,} pares objetivos diferentes. Se conservan los empates, evitando equiparar secuencias y pares objetivos. La preparación común de descriptores tomó {ref['feature_preparation_seconds']:.2f} s y la predicción exhaustiva {ref['exhaustive_prediction_seconds']:.2f} s en esta ejecución. Ambos métodos consultan realmente los modelos en sus lotes, usando matrices de descriptores preparadas. Los tiempos comparados son de búsqueda y seguimiento; excluyen descarga, entrenamiento, preparación y escritura a disco.

### Sensibilidad

Se estudiaron población 30/60/100, mutación por posición 0,5/21, 1/21 y 2/21, y cruce 0,5/0,7/0,9. El diseño varía un factor a la vez de forma secuencial, con semillas 201–205 y 2 100 evaluaciones por corrida. Se elige mayor HV medio, luego menor desviación si hay empate y finalmente menor valor. Las configuraciones repetidas reutilizan exactamente sus resultados verificados por hash. No son réplicas nuevas. La configuración final fue población {config['population']}, mutación {config['mutation']:.6f} por posición y cruce {config['crossover']:.2f}.

{md_table(sensitivity_summary)}

![Calidad y tiempo en sensibilidad](../results/figures/sensibilidad.png)

Cinco semillas por opción ofrecen una exploración inicial, no una caracterización exhaustiva. El ajuste secuencial no estima interacciones. Todo ese costo de ajuste es adicional al presupuesto de una corrida final, y no debe ocultarse al describir la eficiencia.

### Comparación final y convergencia

Las semillas 401–410 no participan en la sensibilidad. La referencia aleatoria muestrea sin reemplazo de manera uniforme todo el dominio. Cada pareja comparte las primeras {config['population']} secuencias. Se alterna el orden temporal de los métodos entre semillas. Ambos terminan en exactamente 2 100 consultas únicas, incluso si la última generación es parcial.

{md_table(summary_table)}

La diferencia media emparejada NSGA-II menos aleatoria es {delta:.4f} puntos porcentuales. El intervalo percentil bootstrap del 95% sobre diez diferencias, con 10 000 remuestreos, es [{ci[0]:.4f}, {ci[1]:.4f}]. Wilcoxon unilateral, con hipótesis NSGA-II > aleatoria fijada en el protocolo, produce `p = {stats['wilcoxon_pvalue']:.9f}`. {significance} Las semillas miden variación de búsqueda dentro de este paisaje; el intervalo no incorpora incertidumbre del conjunto de datos ni del entrenamiento.

![Calidad y tiempo de los dos métodos](../results/figures/comparacion_calidad_tiempo.png)

![Convergencia de hipervolumen](../results/figures/convergencia_hipervolumen.png)

![Promedio y mejor aptitud por objetivo](../results/figures/aptitud_promedio_y_mejor.png)

Las bandas son una desviación estándar entre diez semillas. Las curvas de aptitud describen la población superviviente de NSGA-II. Cada mejor objetivo puede corresponder a una secuencia distinta. El hipervolumen usa el archivo acumulado y por ello su progreso no implica mejora monótona de cada promedio poblacional. Una curva plana indica estancamiento observado; por sí sola no prueba óptimo ni convergencia prematura. Se mide además cuándo se alcanzan umbrales de la referencia, en los puntos de control registrados:

{md_table(hit_summary, 1)}

Cada denominador es diez corridas. La mediana se calcula solo entre las que alcanzan el umbral; `nan` indica ninguna. Los otros casos quedan censurados a 2 100. No se confunde la mediana condicional con el rendimiento de todas las corridas.

## IV. Análisis, conclusiones y límites

### Presupuesto, convergencia y costo

Con población 60, las historias registran 35 generaciones incluyendo la inicial. Ocho de diez corridas finales alcanzan 99% de HV; ninguna recupera las seis secuencias del frente exacto por sí sola. La unión de diez corridas sí lo recupera, con 21 000 consultas sumadas entre corridas y posibles repeticiones entre ellas. Esto respalda una buena aproximación con el presupuesto, sin garantizar recuperación completa.

Las semillas 401, 407 y 409 mejoran aún en el último lote. La 410 gana 5,983 puntos porcentuales entre las evaluaciones 1 500 y 2 100. La 406 se estanca desde 1 440 en 97,363%: existe estancamiento subóptimo observado, pero no se aisló su causa ni se demostró que más generaciones aseguren completitud. No corresponde declarar convergencia garantizada.

En la sensibilidad de población, 60 obtuvo {sensitivity_value('population',60,'HV_media'):.3f}% de HV y {sensitivity_value('population',60,'segundos_media'):.3f} s; 30 obtuvo {sensitivity_value('population',30,'HV_media'):.3f}% y {sensitivity_value('population',30,'segundos_media'):.3f} s; 100 obtuvo {sensitivity_value('population',100,'HV_media'):.3f}% y {sensitivity_value('population',100,'segundos_media'):.3f} s. Se eligió la calidad media, no el menor tiempo. A presupuesto fijo, las poblaciones 30/60/100 registran 70/35/21 generaciones; el número de lotes puede explicar parte del costo, sin constituir una relación causal aislada. La mutación 0,5/21 mejoró a 99,254% frente a 98,320% de 1/21, con tiempo de {sensitivity_value('mutation',0.5/21,'segundos_media'):.3f} frente a {sensitivity_value('mutation',1/21,'segundos_media'):.3f} s. En la etapa posterior, cruce 0,7 obtuvo 99,693%, frente a 99,539% y 99,254% con 0,5 y 0,9. Cinco semillas por opción no bastan para afirmar significación de diferencias pequeñas.

Las 45 filas de sensibilidad representan 35 corridas distintas y {tuning_cost:.2f} s sumados de búsqueda. Las 20 corridas finales agregan {final_cost:.2f} s. En total son 115 500 evaluaciones de candidatos y 231 000 predicciones escalares entre corridas. No son 115 500 secuencias globalmente distintas y no se incluyen entrenamiento ni barrido exhaustivo. Los descriptores se prepararon para todo el dominio: evaluar 2,76% de secuencias por corrida no equivale a realizar solo 2,76% del trabajo total.

La predicción exhaustiva en lote tomó {ref['exhaustive_prediction_seconds']:.2f} s, por debajo del tiempo de búsqueda medio de NSGA-II, pero los cronómetros no incluyen operaciones equivalentes. No se sostiene ventaja temporal frente a enumeración ni necesidad práctica del algoritmo en este dominio pequeño. En el cierre, un reanálisis descriptivo posterior con márgenes de referencia HV 0,05/0,15/0,30 mantuvo la ventaja media: NSGA-II 98,840/98,915/99,011% frente a aleatoria 86,372/87,313/88,500%. No se cambió el criterio primario ni se reajustó el algoritmo. Los cálculos están en `scripts/auditoria_cierre.py` y sus salidas en `results/auditoria_cierre.json`.

### Soluciones y significado práctico

Los candidatos siguientes proceden de la unión de las secuencias realmente evaluadas en las diez corridas finales de NSGA-II. `seeds_found` permite rastrear su hallazgo. No se atribuyen al algoritmo puntos conocidos únicamente por el barrido exhaustivo. Dos criterios pueden seleccionar la misma secuencia. En la referencia, MIC y HC50 predichas fueron {ref['base_prediction']['pred_MIC_uM']:.4f} y {ref['base_prediction']['pred_HC50_uM']:.4f} µM, respectivamente.

{md_table(candidates[['criterio','sequence','mutations_1based','pred_MIC_uM','pred_HC50_uM','pred_ratio_HC50_MIC','seeds_found']])}

![Paisaje y frentes predichos](../results/figures/frentes_pareto.png)

Estas propuestas muestran cómo cambian las predicciones al modificar la secuencia dentro de la región estudiada. No permiten identificar posiciones causalmente responsables ni afirmar una mejora experimental. Incluso un optimizador exacto puede explotar errores del predictor. Si la representación elegida ignora parte del orden, algunos cambios distintos pueden resultar indistinguibles para el modelo. El cociente de dos predicciones agrega la incertidumbre de ambas y no demuestra una ventana terapéutica.

### Conclusión del experimento

{conclusion} en las condiciones declaradas. La diferencia observada es {delta:.4f} puntos porcentuales, y los tiempos de cada método están publicados por separado. La contribución para Computación Evolutiva es la formulación discreta, la adaptación categórica, el control del presupuesto, la sensibilidad y el análisis frente a una referencia exacta del paisaje predicho. No se sostiene que NSGA-II sea universalmente superior ni más rápido por definición.

La limitación central es la validez del evaluador biológico. La reserva interna ayuda a describir su error, pero la agregación entre condiciones, los valores censurados excluidos, las anotaciones faltantes, la selección de la cohorte y las publicaciones compartidas limitan la generalización. El radio de dos cambios tampoco sustituye un análisis formal de dominio de aplicabilidad. Antes de avanzar como investigación biológica se necesitarían datos mejor armonizados, evaluación por publicaciones/familias, tratamiento explícito de censura y validación de candidatos independientes. Como continuidad algorítmica, conviene evaluar varias referencias y radios, incluir el costo total del ajuste y contrastar operadores o reparación alternativos.

El anexo ColabFold del material anterior permanece como evidencia histórica exploratoria. No se reejecutó ni se empleó para validar las variantes de esta versión. Sus métricas estructurales no sustituyen los objetivos ni las mediciones biológicas.

### Correspondencia con los entregables

Las cuatro fases están desarrolladas en este manuscrito y en el PDF de 13 páginas. El código funcional, fuentes crudas, particiones, modelos, semillas, evaluaciones y figuras están en el proyecto. `docs/FASE_I.md` permite preparar el parcial y `docs/TRAZABILIDAD_CURSO.md` ubica cada requisito. Se incorporan los tres integrantes y el repositorio indicado por el equipo. `docs/GUION_DEFENSA.md` organiza la exposición de 30 minutos. Su ensayo y presentación siguen siendo tareas del equipo; la preparación documental no demuestra una exposición ya realizada. El calendario debe confirmarse con el docente: el enunciado indica semana 15 y el sílabo semana 16 para el proyecto final.

## Referencias

1. Pirtskhalava, M., et al. (2021). DBAASP v3: database of antimicrobial/cytotoxic activity and structure of peptides as a resource for development of new therapeutics. *Nucleic Acids Research*, 49(D1), D288–D297. https://doi.org/10.1093/nar/gkaa991
2. Boone, K., et al. (2021). Combining genetic algorithm with machine learning strategies for designing potent antimicrobial peptides. *BMC Bioinformatics*. https://doi.org/10.1186/s12859-021-04156-x
3. Breiman, L. (2001). Random Forests. *Machine Learning*, 45, 5–32. https://doi.org/10.1023/A:1010933404324
4. DBAASP. API oficial y esquema consultados para esta reconstrucción: https://dbaasp.org/api y https://dbaasp.org/v3/api-docs . El esquema y la fecha se conservan en `data/raw/`.
5. Kyte, J., y Doolittle, R. F. (1982). A simple method for displaying the hydropathic character of a protein. *Journal of Molecular Biology*, 157(1), 105–132. https://doi.org/10.1016/0022-2836(82)90515-0
6. Walsh, I., et al. (2021). DOME: recommendations for supervised machine learning validation in biology. *Nature Methods*, 18, 1122–1127. https://doi.org/10.1038/s41592-021-01205-4
7. Deb, K., Pratap, A., Agarwal, S., y Meyarivan, T. (2002). A fast and elitist multiobjective genetic algorithm: NSGA-II. *IEEE Transactions on Evolutionary Computation*, 6(2), 182–197. https://doi.org/10.1109/4235.996017
8. pymoo. Documentación oficial de NSGA-II: https://pymoo.org/algorithms/moo/nsga2.html . La documentación web consultada muestra 0.6.2; esta ejecución fija 0.6.1.5 y verifica sus componentes en el código instalado.

9. pymoo. Performance Indicator: Hypervolume. https://pymoo.org/misc/indicators.html . Consulta: 17 de octubre de 2026.

Fuentes académicas contrastadas en sus páginas primarias. Los valores del experimento se respaldan con archivos de ejecución, no con atribuciones a artículos externos. Los documentos locales que fijan los requisitos son `Curso/Enunciado_Proyecto_Final_Computacion_Evolutiva.docx` y `Curso/SILABO 402 - COMPUTACIÓN EVOLUTIVA_firma.pdf`.
'''
    (docs / "INFORME_RESULTADOS.md").write_text(full_text, encoding="utf-8")
    print("Informe y seis figuras generados desde resultados verificados.", flush=True)
