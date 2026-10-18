"""Informe académico de 13 páginas, generado desde resultados reales.

Dependencias documentales opcionales: reportlab y pypdf. No modifica modelos.
"""
from pathlib import Path
import csv
import json
from datetime import datetime, timezone, timedelta
from html import escape
import reportlab
from reportlab.pdfgen import canvas
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle, Image, PageBreak
from reportlab.lib import colors
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.enums import TA_JUSTIFY
from reportlab.lib.pagesizes import A4
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from pypdf import PdfReader

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output/pdf/Informe_CE_Peptidos.pdf"
OPT = ROOT / "results/optimization"
FIG = ROOT / "results/figures"
def js(path):
    return json.loads(path.read_text(encoding="utf-8"))
def rows(path):
    with path.open(encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f))
def num(v, places=3):
    return f"{float(v):.{places}f}"

def main():
    stats = js(OPT / "estadisticas.json")
    ref = js(OPT / "referencia.json")
    data = js(ROOT / "data/processed/resumen_datos.json")
    raw = js(ROOT / "data/raw/manifest.json")
    audit = js(ROOT / "results/auditoria_resultados.json")
    closure = js(ROOT / "results/auditoria_cierre.json")
    metrics = rows(ROOT / "results/models/validacion_reservada.csv")
    candidates = rows(OPT / "candidatos_encontrados.csv")
    sensitivity = rows(OPT / "sensibilidad.csv")
    histories = rows(OPT / "historia_comparacion.csv")
    environment = js(ROOT / "results/execution_manifest.json")
    arial = Path("C:/Windows/Fonts/arial.ttf")
    font_dir = Path(reportlab.__file__).parent / "fonts"
    regular = arial if arial.exists() else font_dir / "Vera.ttf"
    bold = Path("C:/Windows/Fonts/arialbd.ttf") if arial.exists() else font_dir / "VeraBd.ttf"
    pdfmetrics.registerFont(TTFont("Report", str(regular)))
    pdfmetrics.registerFont(TTFont("ReportBold", str(bold)))
    pdfmetrics.registerFontFamily("Report", normal="Report", bold="ReportBold", italic="Report", boldItalic="ReportBold")
    styles = getSampleStyleSheet()
    styles.add(ParagraphStyle(name="BodyCE", fontName="Report", fontSize=10.2, leading=14.2, spaceAfter=8, alignment=TA_JUSTIFY))
    styles.add(ParagraphStyle(name="SmallCE", fontName="Report", fontSize=8.4, leading=11.2, spaceAfter=6))
    styles.add(ParagraphStyle(name="TitleCE", fontName="ReportBold", fontSize=22, leading=27, textColor=colors.HexColor("#113b50"), spaceAfter=17))
    styles.add(ParagraphStyle(name="H1CE", fontName="ReportBold", fontSize=16, leading=20, textColor=colors.HexColor("#113b50"), spaceAfter=13))
    styles.add(ParagraphStyle(name="H2CE", fontName="ReportBold", fontSize=11.5, leading=15, textColor=colors.HexColor("#087e8b"), spaceBefore=9, spaceAfter=7))
    styles.add(ParagraphStyle(name="CellCE", fontName="Report", fontSize=8.7, leading=11.5))
    story = []
    width = A4[0] - 100
    def p(text, style="BodyCE"):
        story.append(Paragraph(text, styles[style]))
    def heading(text):
        p(text, "H1CE")
    def sub(text):
        p(text, "H2CE")
    def newpage():
        story.append(PageBreak())
    def table(headers, data_rows, widths=None):
        def cell(x):
            return Paragraph(escape(str(x)), styles["CellCE"])
        t = Table([[cell(x) for x in headers]] + [[cell(x) for x in row] for row in data_rows],
                  colWidths=widths, repeatRows=1, hAlign="LEFT")
        t.setStyle(TableStyle([("BACKGROUND", (0,0), (-1,0), colors.HexColor("#e8f1f3")),
                              ("VALIGN", (0,0), (-1,-1), "TOP"),
                              ("LINEBELOW", (0,0), (-1,0), .7, colors.HexColor("#087e8b")),
                              ("LINEBELOW", (0,1), (-1,-1), .3, colors.HexColor("#d5dce0")),
                              ("TOPPADDING", (0,0), (-1,-1), 7), ("BOTTOMPADDING", (0,0), (-1,-1), 7)]))
        story.append(t)
        story.append(Spacer(1, 10))
    def fig(filename, caption, scale=1.):
        image = Image(str(FIG / filename))
        natural_width, natural_height = image.imageWidth, image.imageHeight
        image.drawWidth = width * scale
        image.drawHeight = natural_height * image.drawWidth / natural_width
        story.append(image)
        p(caption, "SmallCE")
    ns = stats["per_method"]["NSGA-II"]
    ra = stats["per_method"]["Aleatoria"]
    delta = stats["paired_mean_difference_percentage_points"]
    ci = stats["paired_bootstrap_95_percent_interval"]
    config = stats["selected_configuration"]

    # 1
    p("MAESTRÍA EN INTELIGENCIA ARTIFICIAL · COMPUTACIÓN EVOLUTIVA", "SmallCE")
    p("Búsqueda evolutiva multiobjetivo de variantes de péptidos", "TitleCE")
    p("Evaluación de NSGA-II sobre propiedades predichas", "H1CE")
    p("Informe técnico · 17 de octubre de 2026<br/>Ing. MIGUEL HERRERA, WILLY ROBERTO<br/>Ing. PARIHUANA TRAVEZAÑO, OMAR EDUARDO<br/>Ing. TORRES CRISTOBAL, WILBER", "SmallCE")
    p("Código y evidencia: <link href='https://github.com/erptyl2010/Proyecto-Computaci-n-Evolutiva' color='#087e8b'>github.com/erptyl2010/Proyecto-Computaci-n-Evolutiva</link>", "SmallCE")
    sub("Resumen")
    p(f"Se evalúa si NSGA-II aprovecha un presupuesto reducido para aproximar el compromiso entre menor MIC y mayor HC50 predichas. Se reconstruyó una cohorte de {data['cohort_api_ids']} fichas de DBAASP y se entrenaron dos Random Forest con reserva interna. El dominio contiene 76 210 secuencias de 21 aminoácidos, con hasta dos sustituciones respecto de una referencia. En diez semillas finales, cada método consultó 2 100 secuencias únicas.")
    table(["Indicador final", "NSGA-II", "Aleatoria"], [
        ["Hipervolumen relativo medio (%)", num(ns['mean_hv_percent']), num(ra['mean_hv_percent'])],
        ["Desviación estándar (puntos porcentuales)", num(ns['std_hv_percent_ddof1']), num(ra['std_hv_percent_ddof1'])],
        ["Tiempo medio de búsqueda (s)", num(ns['mean_seconds']), num(ra['mean_seconds'])]], [255,120,120])
    p(f"La diferencia emparejada fue {delta:.3f} puntos porcentuales. La conclusión se limita a este paisaje predicho y a la implementación evaluada. No se realizaron ensayos de los candidatos; no se demuestra actividad, seguridad ni superioridad biológica.")
    sub("Contribución al curso")
    p("Se desarrollan las cuatro fases: formulación del problema, diseño de operadores categóricos, implementación y experimentación, y análisis de resultados. Se incluyen sensibilidad, promedio y mejor aptitud por generación, comparación con un método simple y referencia exhaustiva del modelo.")
    p("Palabras clave: NSGA-II; optimización multiobjetivo; péptidos; modelos sustitutos; reproducibilidad.", "SmallCE")

    # 2
    newpage(); heading("I. Problema y formulación")
    sub("Pregunta de investigación")
    p("¿Con qué calidad aproxima NSGA-II el frente de compromiso entre MIC y HC50 predichas, frente a una búsqueda aleatoria con el mismo número de secuencias evaluadas? La hipótesis algorítmica es que la selección evolutiva alcanza mayor hipervolumen medio con 2 100 consultas únicas.")
    p("DBAASP ofrece secuencias, modificaciones y mediciones biológicas documentadas [1]. Una MIC menor es la dirección deseada frente a la cepa elegida; una HC50 mayor es la dirección deseada respecto de hemólisis. HC50 no representa todas las formas de toxicidad. Boone et al. constituyen un antecedente de integración de algoritmos genéticos y aprendizaje automático para péptidos, con otro organismo y otra formulación [2]. No se reclama novedad para esa combinación general.")
    sub("Variables y restricciones")
    p("La secuencia inicial es <b>GFMKYIKPLIPHAVKAIKKLI</b>, registrada en la ficha DBAASP 2514. Cada variable representa uno de los 20 aminoácidos canónicos en una de las 21 posiciones. El cromosoma tiene una correspondencia uno a uno con la secuencia.")
    table(["Elemento", "Definición"], [
        ["Variables", "x = (x1, ..., x21), con cada posición en ACDEFGHIKLMNPQRSTVWY."],
        ["Dominio", "Longitud fija 21 y distancia de Hamming a la referencia menor o igual que 2."],
        ["Objetivos a minimizar", "F(x) = (m_MIC(x), -m_HC50(x)). Cada m predice log10(concentración/µM)."],
        ["Contexto terminal", "Ambos extremos sin anotación en las entradas del modelo; un campo nulo no demuestra grupos libres."],
        ["No implementado", "Restricciones de estabilidad, síntesis, agregación o toxicidad general."]], [110,385])
    sub("Cardinalidad y solución")
    p("El espacio tiene <b>1 + 21×19 + C(21,2)×19² = 76 210</b> secuencias. La solución es un conjunto no dominado, no necesariamente una única secuencia. Las concentraciones propuestas son 10 elevado a la predicción logarítmica; no son nuevas mediciones.")
    p("El barrido exhaustivo permite evaluar un presupuesto del 2,7555% del dominio por corrida. Es una referencia exacta de los modelos, no de las propiedades reales. Esta experiencia no demuestra escalabilidad a radios o longitudes mayores.")

    # 3
    newpage(); heading("I. Datos y trazabilidad")
    lima = timezone(timedelta(hours=-5))
    started = datetime.fromisoformat(raw['started_utc']).astimezone(lima)
    finished = datetime.fromisoformat(raw['finished_utc']).astimezone(lima)
    p(f"La descarga de la API oficial [4] se realizó el {started:%d/%m/%Y}, entre las {started:%H:%M:%S} y las {finished:%H:%M:%S} (hora de Lima, UTC-5). Se conservaron las páginas, IDs, esquema, filtros y fichas, con SHA-256. La instantánea contiene {data['cohort_api_ids']} fichas; una descarga futura puede diferir porque la base cambia.")
    p("El filtro solicita monómeros de 10 a 30 residuos, sin aminoácidos inusuales y con registro frente a Escherichia coli ATCC 25922. HC50 se extrae de esa misma cohorte, lo que limita su representatividad frente al conjunto completo de datos de hemólisis.")
    table(["Paso", "MIC", "HC50"], [
        ["Objetivo biológico", "E. coli ATCC 25922", "Eritrocitos humanos"],
        ["Medida aceptada", "MIC", "50% Hemolysis"],
        ["Mediciones puntuales incluidas", data['accepted_measurements']['MIC'], data['accepted_measurements']['HC50']],
        ["Entidades por secuencia y terminales", data['chemical_entities']['MIC'], data['chemical_entities']['HC50']]], [245,125,125])
    p("Se emplea la concentración reportada y su unidad explícita, no el campo distinto llamado activity. Solo se aceptan concentraciones puntuales positivas en µM. Se distinguen ACT, AMD y ausencia de anotación. Se excluyen otras unidades, valores censurados o rangos y química que los descriptores no representan.")
    sub("Exclusiones auditadas")
    names = {"otro_objetivo_biologico":"Otro objetivo biológico", "unidad_no_micromolar":"Unidad distinta de µM", "otra_medida":"Otra medida", "quimica_no_representada":"Química no representada", "concentracion_no_puntual_positiva":"Concentración no puntual/positiva", "terminal_no_representado":"Terminal no representado"}
    table(["Primer motivo de exclusión aplicable", "Mediciones"], [[names.get(k,k),v] for k,v in data['exclusion_reasons'].items()], [375,120])
    p("Los conteos no son causas independientes: cada fila recibe el primer motivo aplicable. Se toma la mediana por secuencia y anotaciones terminales, combinando IDs equivalentes. Las condiciones no se armonizan. Se preservan las referencias de las fichas y reference_raw, sin inventar correspondencias entre artículos y mediciones.", "SmallCE")

    # 4
    newpage(); heading("I. Modelo del proceso y validación")
    p("Se ajusta un Random Forest por propiedad [3]. Se comparan 27 atributos de composición y anotación con 42 que agregan resúmenes de orden. Estos incluyen hidropatía de Kyte-Doolittle [5], carga aproximada y aromaticidad en extremos, posición y productos a desfases. La carga no es un cálculo de ionización dependiente del pH.")
    p("Se agrupan componentes transitivas de secuencias de igual longitud con Hamming ≤ 2. Se reserva el 20% de grupos con semilla 20261017. Cinco pliegues agrupados dentro del entrenamiento seleccionan representación por MAE. Los bosques usan 120 árboles, mínimo tres observaciones por hoja y un hilo; no se reajustan con la reserva [6].")
    table(["Objetivo", "Ajuste / reserva", "MAE log10", "R²", "MAE constante"], [[r['objective'],f"{r['n_train']} / {r['n_holdout']}",num(r['mae_log10']),num(r['r2_log10']),num(r['mae_baseline_log10'])] for r in metrics], [70,110,95,70,150])
    p("Ambos objetivos seleccionaron la representación de 42 atributos. El referente constante predice la media del entrenamiento. La reserva es interna; no demuestra generalización a laboratorios, publicaciones o familias independientes.")
    fig("validacion_predictores.png", "Figura 1. Predicciones frente a observaciones en los casos reservados. La línea representa igualdad.")
    shares = "; ".join(f"{r['objective']}: {r.get('shared_publications_train_holdout','pendiente')}" for r in metrics)
    p(f"El R² de HC50 es {float(metrics[1]['r2_log10']):.3f}: la capacidad para explicar su variación es limitada. Publicaciones compartidas por título normalizado y año: {shares}. IDs internos distintos no garantizan artículos distintos. El agrupamiento tampoco elimina toda homología.", "SmallCE")

    # 5
    newpage(); heading("II. Diseño del algoritmo evolutivo")
    p("NSGA-II incorpora elitismo, clasificación no dominada y preservación de diversidad [7]. Se usa pymoo 0.6.1.5, verificando su código instalado; la página oficial consultada muestra una versión posterior [8]. El tratamiento categórico evita interpretar el índice alfabético de un aminoácido como una distancia continua.")
    p("Se eligió pymoo porque permite integrar cruce, mutación, reparación y control de duplicados propios con su supervivencia multiobjetivo. La versión se fija y la métrica se contrasta independientemente.", "SmallCE")
    table(["Componente", "Implementación"], [
        ["Inicialización", "Muestra uniforme sin reemplazo del dominio; misma población inicial en cada pareja de métodos."],
        ["Codificación", "21 enteros categóricos: una única representación por secuencia."],
        ["Cruce", "Uniforme de dos progenitores, con intercambio por posición y dos descendientes."],
        ["Mutación", "Por posición seleccionada, escoger uniformemente uno de los otros 19 aminoácidos."],
        ["Reparación", "Restaurar al azar cambios excedentes a la referencia hasta tener como máximo dos."],
        ["Selección parental", "Torneo binario por dominancia y distancia de hacinamiento, según la versión instalada."],
        ["Supervivencia", "Unión de progenitores y descendientes; frentes y hacinamiento, con elitismo."],
        ["Duplicados", "Eliminar secuencias del lote y cualquiera ya evaluada durante la misma corrida."],
        ["Parada", "Exactamente 2 100 consultas únicas; ajustar el tamaño del último lote si es necesario."]], [115,380])
    sub("Flujo de una corrida")
    p("1. Crear población válida y evaluar los dos modelos congelados.<br/>2. Seleccionar padres, cruzar, mutar y reparar.<br/>3. Eliminar duplicados antes de consultar los modelos.<br/>4. Evaluar y aplicar supervivencia elitista.<br/>5. Actualizar archivo, promedios, mejores valores e hipervolumen.<br/>6. Repetir hasta consumir el presupuesto.")
    p("El archivo conserva todas las soluciones consultadas, incluso las que dejan la población. La reparación y la exclusión global de duplicados alteran la dinámica del algoritmo; los resultados corresponden a esta variante implementada.")

    # 6
    newpage(); heading("III. Protocolo y sensibilidad")
    p("La aptitud es el vector predicho; la métrica externa es el hipervolumen (HV) del archivo: área dominada respecto de una referencia, donde más es mejor [9]. Se usa F = (log10 MIC, -log10 HC50) y una referencia común: máximo de cada componente del dominio más 0,15. La métrica relativa es 100×HV(archivo)/HV(exhaustivo). Los operadores consultan los modelos, no las respuestas exhaustivas.")
    p("La sensibilidad varía un factor a la vez, secuencialmente, con semillas 201-205. Se escoge la mayor media de HV, luego menor desviación y menor valor en caso de empate. Las configuraciones repetidas reutilizan archivos verificados; no constituyen réplicas nuevas. Este diseño no estima interacciones.")
    groups = {}
    for r in sensitivity:
        groups.setdefault((r['sensitivity_parameter'],r['sensitivity_value']),[]).append(float(r['hv_percent']))
    import statistics
    labels = {"population":"Población", "mutation":"Mutación/posición", "crossover":"Cruce"}
    table(["Factor", "Valor", "HV medio (%)", "DE"], [[labels[k[0]],num(k[1],5) if k[0]=='mutation' else num(k[1],2),num(statistics.mean(v)),num(statistics.stdev(v))] for k,v in groups.items()], [180,100,125,90])
    fig("sensibilidad.png", "Figura 2. Calidad y tiempo en sensibilidad; media ± una desviación entre cinco semillas.", .96)
    p(f"Configuración elegida: población {config['population']}, mutación {config['mutation']:.6f} por posición y cruce {config['crossover']:.2f}. El costo del ajuste es adicional al de una corrida final.", "SmallCE")

    # 7
    newpage(); heading("III. Comparación final")
    p("Las diez semillas finales, 401-410, no participaron en la sensibilidad. La referencia aleatoria muestrea uniformemente sin reemplazo. Cada pareja comparte las primeras 60 secuencias. Se alterna el orden temporal de los métodos para reducir sesgos de medición. Cada corrida evalúa 2 100 secuencias diferentes.")
    table(["Medida", "Resultado"], [
        ["Diferencia media NSGA-II - aleatoria", f"{delta:.4f} puntos porcentuales de HV"],
        ["IC bootstrap percentil 95%", f"[{ci[0]:.4f}, {ci[1]:.4f}] puntos porcentuales"],
        ["Remuestreo", "10 000 muestras de las diez diferencias; semilla 20261017"],
        ["Wilcoxon unilateral (NSGA-II > aleatoria)", f"p = {stats['wilcoxon_pvalue']:.9f}"],
        ["Alcance de la inferencia", "Variación de búsqueda en un único paisaje predicho y congelado"]], [250,245])
    fig("comparacion_calidad_tiempo.png", "Figura 3. Distribución entre diez semillas finales. Los puntos representan corridas individuales.")
    p("La ventaja de calidad observada no implica menor tiempo. La búsqueda aleatoria resultó más rápida en estas mediciones. El tiempo incluye consultas de los modelos, operadores y seguimiento; excluye descarga, entrenamiento, descriptores comunes y escritura a disco.")
    p(f"La preparación común de descriptores tomó {ref['feature_preparation_seconds']:.2f} s; las predicciones exhaustivas, {ref['exhaustive_prediction_seconds']:.2f} s. El barrido halló {ref['n_front_sequences']} secuencias no dominadas y {ref['n_front_objective_pairs']} pares objetivos distintos. No deben confundirse ambos conteos.")

    # 8
    newpage(); heading("III. Convergencia y aptitud")
    fig("convergencia_hipervolumen.png", "Figura 4. Progreso del HV acumulado; media ± una desviación estándar, diez semillas.")
    fig("aptitud_promedio_y_mejor.png", "Figura 5. Promedio y mejor valor de cada objetivo en la población superviviente, por generación.")
    p("Los mejores valores de MIC y HC50 pueden corresponder a individuos distintos. El HV mide el archivo acumulado; sus mejoras no obligan a que cada promedio poblacional mejore de forma monótona. Las bandas muestran dispersión, no intervalos de confianza.")
    thresholds = rows(OPT / "evaluaciones_hasta_umbral.csv")
    summaries = []
    for method in ["NSGA-II", "Aleatoria"]:
        for threshold in [99.,99.9]:
            values = [float(r['first_recorded_evaluations']) for r in thresholds if r['method']==method and float(r['threshold_percent'])==threshold and r['reached']=='True']
            summaries.append([method,f"{threshold}%",f"{len(values)}/10",num(statistics.median(values),0) if values else "No alcanzado"])
    table(["Método", "Umbral", "Éxitos", "Mediana de evaluaciones"], summaries, [120,80,70,225])
    p("Medianas solo entre éxitos; los demás casos quedan censurados a 2 100. Se detecta el primer punto de control registrado, no una evaluación intermedia. Una meseta por sí sola no demuestra optimalidad ni convergencia prematura.", "SmallCE")

    # 9: cierre del análisis exigido en la fase IV
    newpage(); heading("IV. Presupuesto y eficiencia")
    sub("¿Fueron suficientes 2 100 evaluaciones?")
    p("Con población 60 se registran 35 generaciones, incluida la inicial. Ocho de diez corridas alcanzan 99% de HV, pero ninguna recupera por sí sola las seis secuencias del frente exacto: encuentran entre dos y cinco. La unión de diez corridas sí recupera las seis; consume 21 000 consultas sumadas entre corridas, con posibles repeticiones entre ellas.")
    p("El presupuesto logra una aproximación alta en la mayoría de semillas, pero no garantiza recuperar todo el frente. Las semillas 401, 407 y 409 aún mejoran en el último lote; la 410 gana 5,983 puntos de HV entre las evaluaciones 1 500 y 2 100. La 406 se estanca desde 1 440 en 97,363%. Hay estancamiento subóptimo observado, sin evidencia suficiente para atribuir su causa ni afirmar que más iteraciones asegurarían el óptimo.")
    sub("Qué aporta la sensibilidad")
    p("En la etapa de población, 60 consigue 98,320% de HV medio y 1,052 s; 30 consigue 97,840% y 2,111 s; 100 obtiene 97,693% y 0,841 s. Se elige por calidad, no por menor tiempo. Cambiar población con presupuesto fijo cambia también el número de generaciones (70, 35 y 21). El mayor número de lotes es una explicación plausible del costo observado, no una causa aislada experimentalmente.")
    p("Con población 60, mutación 0,5/21 mejora la media a 99,254% frente a 98,320% con 1/21; su tiempo aumenta de 1,052 a 1,160 s. El cruce 0,7 alcanza después 99,693% y 1,193 s, frente a 99,539% con 0,5 y 99,254% con 0,9. Son comparaciones descriptivas de cinco semillas, sin estimar interacciones ni probar diferencias pequeñas.")
    sub("Costo y alcance de la métrica")
    cost = closure['cost_accounting']
    p(f"La sensibilidad suma {cost['distinct_tuning_runs']} corridas distintas y {cost['tuning_search_seconds_sum']:.2f} s de búsqueda. Las 20 corridas finales suman {cost['final_search_seconds_sum']:.2f} s. En total son 115 500 evaluaciones de candidatos entre corridas, excluyendo entrenamiento y barrido exhaustivo. Se prepararon descriptores para todo el dominio: consultar el 2,76% por corrida no significa realizar solo el 2,76% del trabajo total.")
    p("La predicción exhaustiva en lote tardó 0,84 s; esta experiencia no demuestra una necesidad práctica de NSGA-II ni una ventaja de velocidad sobre enumeración. Los cronómetros no incluyen operaciones equivalentes. HV es pertinente para comparar conjuntos de compromisos bajo la misma escala y referencia, pero no mide la fracción de soluciones recuperadas ni éxito biológico.")
    p("El logaritmo conserva el orden de dominancia, pero cambia la geometría de HV. Una comprobación posterior con márgenes 0,05/0,15/0,30 mantiene la ventaja media de NSGA-II (98,840/98,915/99,011%) frente a aleatoria (86,372/87,313/88,500%). Es un análisis descriptivo adicional; no se cambió el criterio primario ni se reajustó el algoritmo.", "SmallCE")

    # 10
    newpage(); heading("IV. Soluciones e interpretación")
    p("Estos candidatos pertenecen a la unión de las secuencias realmente evaluadas en las diez corridas finales de NSGA-II. No se atribuyen a la búsqueda soluciones obtenidas solo por enumeración. Los tres criterios son elecciones explícitas y pueden seleccionar una misma secuencia.")
    for r in candidates:
        sub(r['criterio'])
        p(f"<b>{r['sequence']}</b><br/>Cambios (posición desde 1): {escape(r['mutations_1based'])}<br/>MIC predicha: {num(r['pred_MIC_uM'])} µM; HC50 predicha: {num(r['pred_HC50_uM'])} µM; cociente: {num(r['pred_ratio_HC50_MIC'])}.", "SmallCE")
        p("Encontrada en semillas: " + r['seeds_found'] + ".", "SmallCE")
    fig("frentes_pareto.png", "Figura 6. Paisaje predicho, frente exhaustivo y frente de la unión de diez corridas finales.", .88)
    p(f"En la referencia, MIC y HC50 predichas fueron {ref['base_prediction']['pred_MIC_uM']:.3f} y {ref['base_prediction']['pred_HC50_uM']:.3f} µM. Las diferencias anteriores se refieren a estos predictores. No identifican efectos causales ni prueban una ventana terapéutica.", "SmallCE")

    # 11
    newpage(); heading("IV. Discusión y conclusiones")
    sub("Qué queda respaldado")
    p(f"NSGA-II obtuvo mayor HV medio que la búsqueda aleatoria en esta comparación: {ns['mean_hv_percent']:.3f}% frente a {ra['mean_hv_percent']:.3f}%, con un costo temporal superior. El hallazgo está condicionado a una secuencia de referencia, dos modelos, operadores concretos y el presupuesto declarado. No demuestra superioridad universal.")
    p("La contribución para el curso es una cadena verificable: formulación discreta, selección de operadores, control de secuencias únicas, sensibilidad, curvas de promedio y mejor aptitud, y comparación con una referencia conocida. La enumeración permite medir calidad con presupuesto reducido, aunque aquí sea viable evaluar todo.")
    sub("Qué limita la interpretación")
    p("<b>Evaluador:</b> el error reservado, en especial la capacidad limitada de HC50, restringe la credibilidad biológica del paisaje. El optimizador puede aprovechar errores del predictor. Una solución óptima para el modelo no es automáticamente óptima en laboratorio.")
    p("<b>Datos:</b> la mediana combina condiciones no armonizadas; la exclusión de censura puede sesgar la distribución. La cohorte de HC50 está condicionada por la búsqueda frente a E. coli. Anotaciones faltantes y coincidencias bibliográficas requieren atención. La prueba es interna y el criterio Hamming no elimina toda homología.")
    p("<b>Experimento:</b> cinco semillas por opción dan una sensibilidad inicial; el diseño no estima interacciones. El intervalo estadístico final refleja la aleatoriedad de búsqueda, sin incluir variación por datos o entrenamiento. No se estableció un dominio formal de aplicabilidad para las variantes.")
    sub("Continuación razonable")
    p("Primero, armonizar condiciones y evaluar separación por publicación o familia, además de tratar explícitamente valores censurados. Luego, estudiar incertidumbre de las predicciones y varias referencias o radios de búsqueda. Las mejoras de operadores deben compararse incluyendo el costo de ajuste. Una afirmación biológica requiere validación independiente de los candidatos.")
    p("El anexo ColabFold anterior se conserva como material histórico exploratorio. No se reejecutó y no valida las variantes nuevas. El proyecto actual puede defenderse como estudio de Computación Evolutiva si mantiene explícitos estos límites.")

    # 12
    newpage(); heading("Reproducibilidad y entrega del curso")
    p(f"La auditoría comprobó {audit['raw_files_verified']} archivos crudos y {audit['unique_runs_verified']} corridas distintas, además de presupuestos, métricas y procedencia. Se repitió una corrida real y se verificó el mismo orden de secuencias y predicciones hasta tolerancia 10<super>-12</super>. Los tiempos no se exigen idénticos.")
    clean = closure['clean_reproduction']
    if clean['status'] != 'verified':
        raise RuntimeError('Falta verificar la reproducción limpia para cerrar este informe.')
    p(f"En el cierre se volvió a procesar la instantánea, entrenar y ejecutar todas las corridas en una copia sin resultados previos. Coincidieron los modelos y {clean['tables_compared']} tablas, incluidas las evaluaciones e historias de las 55 corridas, con tolerancia absoluta 10<super>-12</super>. Se usaron la misma máquina y entorno; no se certifica portabilidad a otro sistema.")
    table(["Evidencia", "Archivo o carpeta"], [
        ["Origen, filtros y huellas", "data/raw/manifest.json y peptides/"],
        ["Mediciones y exclusiones", "data/processed/auditoria_mediciones.csv"],
        ["Modelos, reserva y selección", "results/models/"],
        ["Corridas, semillas e historias", "results/optimization/runs/"],
        ["Candidatos y su procedencia", "results/optimization/candidatos_encontrados.csv"],
        ["Auditoría de resultados reales", "results/auditoria_resultados.json"],
        ["Comandos de reproducción", "README.md y run.py"],
        ["Código y pruebas del núcleo", "ce_peptidos/ y tests/test_core.py"]], [215,280])
    p("El bloqueo de paquetes está en requirements-lock.txt. El entorno medido usa Python " + environment['python'].split()[0] + ". Las versiones principales son: " + "; ".join(f"{k} {v}" for k,v in environment['packages'].items()) + ".", "SmallCE")
    sub("Correspondencia y pendientes formales")
    p("Fase I: definición, formulación y modelo empírico. Fase II: paradigma, codificación, selección y operadores. Fase III: ejecución, sensibilidad y métrica. Fase IV: convergencia, comparación, candidatos, límites y trabajo futuro. El detalle se encuentra en docs/TRAZABILIDAD_CURSO.md.")
    p("El PDF desarrolla las cuatro fases en trece páginas. La portada identifica a los integrantes e incluye el repositorio indicado por el equipo. El guion organiza una sustentación de 30 minutos; su preparación no acredita que ya se haya realizado. El calendario debe confirmarse con el docente: semana 15 en el enunciado y semana 16 en el sílabo.")
    p("Las referencias bibliográficas sustentan métodos y contexto. Las cifras nuevas se sustentan en los archivos de ejecución, sin adjudicarlas a los artículos ni mezclarlas con los resultados históricos.")

    # 13
    newpage(); heading("Referencias")
    bibliography = [
        ("1", "Pirtskhalava, M., et al. (2021). DBAASP v3: database of antimicrobial/cytotoxic activity and structure of peptides as a resource for development of new therapeutics. Nucleic Acids Research, 49(D1), D288-D297.", "https://doi.org/10.1093/nar/gkaa991"),
        ("2", "Boone, K., Wisdom, C., Camarda, K., et al. (2021). Combining genetic algorithm with machine learning strategies for designing potent antimicrobial peptides. BMC Bioinformatics, 22, 239.", "https://doi.org/10.1186/s12859-021-04156-x"),
        ("3", "Breiman, L. (2001). Random Forests. Machine Learning, 45, 5-32.", "https://doi.org/10.1023/A:1010933404324"),
        ("4", "DBAASP. API oficial y esquema OpenAPI. Consulta e instantánea documentadas el 17 de octubre de 2026. El esquema original se incluye en data/raw/api_schema.json.", "https://dbaasp.org/api"),
        ("5", "Kyte, J., y Doolittle, R. F. (1982). A simple method for displaying the hydropathic character of a protein. Journal of Molecular Biology, 157(1), 105-132.", "https://doi.org/10.1016/0022-2836(82)90515-0"),
        ("6", "Walsh, I., et al. (2021). DOME: recommendations for supervised machine learning validation in biology. Nature Methods, 18, 1122-1127.", "https://doi.org/10.1038/s41592-021-01205-4"),
        ("7", "Deb, K., Pratap, A., Agarwal, S., y Meyarivan, T. (2002). A fast and elitist multiobjective genetic algorithm: NSGA-II. IEEE Transactions on Evolutionary Computation, 6(2), 182-197.", "https://doi.org/10.1109/4235.996017"),
        ("8", "pymoo. Documentación oficial de NSGA-II. La página consultada corresponde a 0.6.2; el experimento fija y revisa el código de 0.6.1.5.", "https://pymoo.org/algorithms/moo/nsga2.html"),
        ("9", "pymoo. Performance Indicator: Hypervolume. Definición de área dominada y punto de referencia. Consulta: 17 de octubre de 2026.", "https://pymoo.org/misc/indicators.html")]
    for index, text, link in bibliography:
        p(f"<b>[{index}]</b> {escape(text)}<br/><link href='{escape(link, quote=True)}' color='#087e8b'>{escape(link)}</link>", "SmallCE")
        story.append(Spacer(1,6))
    sub("Documentos que fijan los requisitos")
    p("Universidad Nacional de Ingeniería. MAE-IA: Computación Evolutiva 2026-2. Enunciado del proyecto final y sílabo proporcionados por el estudiante en la carpeta Curso. Se revisaron sus fases, entregables, rúbrica y calendario. No se ha inventado una fecha de entrega ni una URL de repositorio.", "SmallCE")
    sub("Disponibilidad")
    p("Los JSON originales, datos procesados, código, modelos y resultados acompañan este informe en Proyecto_CE_Peptidos. Las fuentes históricas permanecen aparte. README.md explica la reproducción y la función de cada archivo.", "SmallCE")
    OUT.parent.mkdir(parents=True, exist_ok=True)
    def footer(c, doc):
        c.setStrokeColor(colors.HexColor("#c8d9de"))
        c.line(50, 43, A4[0]-50, 43)
        c.setFont("Report", 8)
        c.setFillColor(colors.HexColor("#51636b"))
        c.drawString(50, 29, "Computación Evolutiva · Estudio sobre propiedades predichas")
        c.drawRightString(A4[0]-50, 29, str(doc.page))
    doc = SimpleDocTemplate(str(OUT), pagesize=A4, rightMargin=50, leftMargin=50, topMargin=46, bottomMargin=57,
                            title="Búsqueda evolutiva multiobjetivo de variantes de péptidos", author="MIGUEL HERRERA, WILLY ROBERTO; PARIHUANA TRAVEZAÑO, OMAR EDUARDO; TORRES CRISTOBAL, WILBER")
    doc.build(story, onFirstPage=footer, onLaterPages=footer)
    reader = PdfReader(str(OUT))
    print("Páginas:", len(reader.pages))
    if len(reader.pages) != 13:
        raise RuntimeError("Revisar paginación: el diseño prevé 13 páginas.")
    print(OUT)

if __name__ == "__main__":
    main()
