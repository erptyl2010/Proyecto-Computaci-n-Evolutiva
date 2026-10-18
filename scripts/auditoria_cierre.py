"""Comprobaciones adicionales de cierre sobre resultados existentes.

No cambia modelos, candidatos ni configuración. --reproduction-root compara
una ejecución independiente iniciada sin procesados, modelos o caché de corridas.
Los tiempos se excluyen de la igualdad porque dependen de la máquina y su carga.
"""
from pathlib import Path
from datetime import datetime, timezone
import argparse
import hashlib
import json
import numpy as np
import pandas as pd
from pandas.testing import assert_frame_equal
from pymoo.indicators.hv import HV
from scipy.stats import wilcoxon

ROOT = Path(__file__).resolve().parents[1]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def require(condition, message):
    if not condition:
        raise ValueError(message)


def compare_table(a, b):
    left, right = pd.read_csv(a, keep_default_na=False), pd.read_csv(b, keep_default_na=False)
    columns = [c for c in left.columns if c not in {'seconds', 'elapsed_seconds'}]
    assert_frame_equal(left[columns], right[columns], check_exact=False, rtol=0, atol=1e-12)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--reproduction-root', type=Path)
    args = parser.parse_args()
    opt = ROOT / 'results/optimization'
    final = pd.read_csv(opt / 'comparacion_final.csv')
    tuning = pd.read_csv(opt / 'sensibilidad.csv')
    stats = json.loads((opt / 'estadisticas.json').read_text(encoding='utf-8'))
    oracle = pd.read_csv(opt / 'dominio_exhaustivo.csv')
    F = oracle[['pred_log10_MIC', 'pred_log10_HC50']].to_numpy() * [1, -1]
    exact = set(oracle.loc[oracle.nondominated, 'sequence'])
    report = {'created_utc': datetime.now(timezone.utc).isoformat(),
              'scope': 'Auditoría del experimento computacional; no validación biológica.',
              'source_tables_sha256': {p.name: sha(p) for p in [
                  opt/'comparacion_final.csv', opt/'sensibilidad.csv', opt/'dominio_exhaustivo.csv']}}
    paired = final.pivot(index='seed', columns='method', values='hv_percent')
    differences = paired['NSGA-II'].to_numpy() - paired['Aleatoria'].to_numpy()
    require(np.isclose(differences.mean(), stats['paired_mean_difference_percentage_points'], rtol=0, atol=1e-12), 'Diferencia media incoherente')
    pvalue = float(wilcoxon(differences, alternative='greater').pvalue)
    require(pvalue == stats['wilcoxon_pvalue'], 'Wilcoxon incoherente')
    report['paired_difference_recomputed'] = float(differences.mean())
    report['wilcoxon_recomputed'] = pvalue
    summary = tuning.groupby(['sensitivity_parameter', 'sensitivity_value']).agg(
        mean_hv=('hv_percent', 'mean'), sd_hv=('hv_percent', 'std'), mean_seconds=('seconds', 'mean')).reset_index()
    for step in stats['selection_steps']:
        choices = summary[summary.sensitivity_parameter == step['parameter']]
        selected = choices.sort_values(['mean_hv', 'sd_hv', 'sensitivity_value'], ascending=[False, True, True]).iloc[0]
        require(np.isclose(selected.sensitivity_value, step['selected']), 'Selección de sensibilidad incoherente')
    report['sensitivity_summary'] = summary.to_dict('records')
    unique_tuning = tuning.drop_duplicates('run_id')
    all_runs = pd.concat([unique_tuning, final], ignore_index=True).drop_duplicates('run_id')
    report['cost_accounting'] = {
        'tuning_table_rows': len(tuning), 'distinct_tuning_runs': len(unique_tuning),
        'tuning_search_seconds_sum': float(unique_tuning.seconds.sum()),
        'final_runs': len(final), 'final_search_seconds_sum': float(final.seconds.sum()),
        'distinct_runs_total': len(all_runs), 'candidate_evaluations_across_runs': int(all_runs.unique_evaluations.sum()),
        'scalar_model_predictions_across_runs': int(2 * all_runs.unique_evaluations.sum()),
        'note': 'Los conteos suman corridas: una secuencia puede reaparecer en otras semillas. Excluyen el barrido exhaustivo, entrenamiento y preparación.'}
    convergence, hv_rows, union = [], [], set()
    for row in final.itertuples():
        folder = opt / 'runs' / row.run_id
        archive = pd.read_csv(folder/'evaluaciones.csv')
        history = pd.read_csv(folder/'historia.csv')
        found = exact & set(archive.sequence)
        if row.method == 'NSGA-II':
            union.update(found)
            previous = history[history.n_evaluations <= 1500].iloc[-1]
            improvements = history[history.hv_percent.diff() > 1e-10]
            convergence.append({'seed': int(row.seed), 'generations_including_initial': len(history),
                'final_hv_percent': float(history.iloc[-1].hv_percent),
                'hv_at_1500': float(previous.hv_percent),
                'last_600_evaluations_gain_pp': float(history.iloc[-1].hv_percent - previous.hv_percent),
                'last_recorded_improvement_evaluations': int(improvements.iloc[-1].n_evaluations),
                'exact_front_sequences_found': len(found)})
        af = archive[['pred_log10_MIC', 'pred_log10_HC50']].to_numpy() * [1, -1]
        # Reanálisis descriptivo posterior: no se reajusta ni selecciona el margen.
        for margin in [.05, .15, .30]:
            indicator = HV(ref_point=F.max(axis=0) + margin)
            relative = float(100 * indicator(af) / indicator(F[oracle.nondominated]))
            hv_rows.append({'margin': margin, 'method': row.method, 'seed': int(row.seed), 'hv_percent': relative})
            if margin == .15:
                require(np.isclose(relative, row.hv_percent, rtol=0, atol=1e-10), 'HV primario incoherente')
    pd.DataFrame(convergence).to_csv(ROOT/'results/auditoria_convergencia.csv', index=False)
    hv = pd.DataFrame(hv_rows)
    hv.to_csv(ROOT/'results/auditoria_referencia_hv.csv', index=False)
    report['convergence'] = convergence
    report['exact_front_recovery'] = {'reference_sequences': len(exact), 'union_ten_final_nsga_runs': len(union),
        'individual_complete_runs': sum(x['exact_front_sequences_found'] == len(exact) for x in convergence)}
    report['posthoc_hv_reference_check'] = hv.groupby(['margin', 'method']).hv_percent.mean().reset_index().to_dict('records')
    candidates = pd.read_csv(opt/'candidatos_encontrados.csv').sequence.unique()
    familiarity = []
    for objective in ['MIC', 'HC50']:
        data = pd.read_csv(ROOT/f'data/processed/dataset_{objective}.csv', keep_default_na=False)
        for seq in candidates:
            same_length = data[data.sequence.str.len() == len(seq)].sequence
            minimum = min(sum(a != b for a, b in zip(seq, s)) for s in same_length)
            familiarity.append({'objective': objective, 'sequence': seq,
                'exact_sequence_in_processed_cohort': bool((data.sequence == seq).any()),
                'minimum_hamming_to_same_length_in_cohort': minimum})
    report['candidate_presence_in_processed_cohort'] = familiarity
    report['novelty_scope'] = 'La búsqueda de coincidencias solo cubre esta cohorte procesada; no certifica novedad ni dominio de aplicabilidad.'
    if args.reproduction_root:
        other = args.reproduction_root.resolve()
        source_files = [ROOT/'run.py', *sorted((ROOT/'ce_peptidos').glob('*.py'))]
        # El generador de narración puede revisarse sin alterar la reproducción numérica.
        source_files = [p for p in source_files if p.name != 'informe.py']
        for p in source_files:
            require(sha(p) == sha(other/p.relative_to(ROOT)), f'Código numérico distinto: {p.name}')
        checks = []
        for directory in ['data/processed', 'results/models']:
            for path in sorted((ROOT/directory).glob('*.csv')):
                rel = path.relative_to(ROOT)
                compare_table(path, other/rel)
                checks.append(rel.as_posix())
        for name in ['dominio_exhaustivo.csv', 'comparacion_final.csv', 'sensibilidad.csv', 'candidatos_encontrados.csv']:
            rel = Path('results/optimization')/name
            compare_table(ROOT/rel, other/rel)
            checks.append(rel.as_posix())
        for run in all_runs.run_id:
            for name in ['evaluaciones.csv', 'historia.csv']:
                rel = Path('results/optimization/runs')/run/name
                compare_table(ROOT/rel, other/rel)
                checks.append(rel.as_posix())
        model = Path('results/models/modelos_congelados.joblib')
        require(sha(ROOT/model) == sha(other/model), 'La huella del modelo no coincide')
        report['clean_reproduction'] = {
            'status': 'verified', 'distinct_runs_compared': len(all_runs), 'tables_compared': len(checks),
            'absolute_tolerance': 1e-12, 'relative_tolerance': 0,
            'excluded_time_columns': ['seconds', 'elapsed_seconds'], 'model_sha256': sha(ROOT/model),
            'numerical_sources_sha256': {p.relative_to(ROOT).as_posix(): sha(p) for p in source_files},
            'checked_tables': checks,
            'environment': 'Misma máquina y entorno virtual; ejecución desde copia limpia de código e instantánea cruda. No prueba portabilidad a otro sistema.'}
    else:
        report['clean_reproduction'] = {'status': 'not_requested'}
    path = ROOT/'results/auditoria_cierre.json'
    path.write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding='utf-8')
    print(json.dumps(report, indent=2, ensure_ascii=False))


if __name__ == '__main__':
    main()
