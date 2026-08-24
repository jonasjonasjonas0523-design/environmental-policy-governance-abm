"""Gate C4 per V1_2_CODE_CHECKS.md: module, constant, and stream isolation.
Underwrites reuse of all archived S/T results."""
import numpy as np, json, sys, os, inspect
import pandas as pd
sys.path.insert(0, os.path.dirname(__file__))
import model as m11
import model_v12 as m12

OUT = {}
# ---- static: function-level byte identity except the declared diff ----
shared = ['logistic', 'seed_data', 'batch', '_quad', 'F_of_tau', 'solve_tau_star',
          'M_of', '_base_out', 'run_standard', 'run_charge', '_finish', 'export_config', 'run_all']
static = {}
for fn in shared:
    static[fn] = inspect.getsource(getattr(m11, fn)) == inspect.getsource(getattr(m12, fn))
changed_ok = inspect.getsource(m11.run_voluntary) != inspect.getsource(m12.run_voluntary)
cfg_diff = {k for k in set(m11.CONFIG) ^ set(m12.CONFIG)} | \
           {k for k in set(m11.CONFIG) & set(m12.CONFIG) if m11.CONFIG[k] != m12.CONFIG[k]}
tau_frozen = 0.09326697147369385
tau_cfg = json.load(open('results/verification/tau_star.json'))['tau_star']
OUT['static'] = dict(shared_functions_byte_identical=static,
                     all_shared_identical=bool(all(static.values())),
                     run_voluntary_changed=bool(changed_ok),
                     config_diff_keys=sorted(cfg_diff),
                     config_diff_ok=(cfg_diff == {'nu', 'nu0', 'eta'}),
                     tau_star_config_byte_equal=(repr(tau_cfg) == repr(tau_frozen)),
                     f_S=repr(m12.CONFIG['m_S'] * tau_frozen),
                     f_S_expected='0.1865339429473877',
                     f_S_ok=(m12.CONFIG['m_S'] * tau_frozen == 0.1865339429473877))

# ---- dynamic: one principal S and one principal T condition, 3 seeds, vs archived v1.1 ----
raw = pd.read_parquet('results/principal/principal_raw_results.parquet')
seeds3 = [100000, 100001, 100002]
b = m12.batch(seeds3)
KR, HR = m12.CONFIG['ref']['k'], m12.CONFIG['ref']['h']
# pick two archived cells: reference (0.5,0.5) for S and (0.9,0.9) for T (both on the 41-grid)
pp = np.linspace(0, 1, 41)
cells = {'S': (float(pp[20]), float(pp[20])), 'T': (float(pp[36]), float(pp[36]))}
dyn = {}
metrics = ['A','C_ab','C_social','nM','nA','nI','T_charges','T_fines','viol','under','part','CE']
key = {'T_charges':'T_ch','T_fines':'T_fi'}
for instr, (p, c) in cells.items():
    run = m12.run_standard if instr == 'S' else m12.run_charge
    o = run(b, dict(p=p, c=c, k=KR, h=HR), tau_frozen)
    arch = raw[(raw.instrument == instr) & (np.isclose(raw.p, p)) & (np.isclose(raw.c, c)) &
               (raw.seed.isin(seeds3))].sort_values('seed')
    ok_all = True; det = {}
    for mt in metrics:
        new = np.asarray(o[key.get(mt, mt)], dtype=np.float32)
        old = arch[mt].to_numpy(dtype=np.float32)
        same = bool(np.array_equal(new, old, equal_nan=True))
        det[mt] = same; ok_all &= same
    dyn[f'{instr}@p={p:.3f},c={c:.3f}'] = dict(bitwise_identical=ok_all, per_metric=det)
OUT['dynamic'] = dyn

# ---- stream audit: RNG consumption order and per-instrument stream usage ----
sd_src = inspect.getsource(m12.seed_data)
order_ok = sd_src.find("normal(cfg['mu_g']") < sd_src.find("normal(cfg['mu_E']") < \
           sd_src.find("U_mon=") < sd_src.find("U_imp=") < sd_src.find("U_imp2=")
use = {}
for fn in ('run_standard', 'run_charge', 'run_voluntary'):
    s = inspect.getsource(getattr(m12, fn))
    use[fn] = dict(U_mon='U_mon' in s, U_imp="['U_imp']" in s, U_imp2='U_imp2' in s)
stream_ok = (order_ok and use['run_standard'] == dict(U_mon=True, U_imp=True, U_imp2=False)
             and use['run_charge'] == dict(U_mon=True, U_imp=True, U_imp2=True)
             and use['run_voluntary'] == dict(U_mon=True, U_imp=False, U_imp2=False))
OUT['stream_audit'] = dict(seed_data_draw_order='gamma, E0, U_mon[T,N], U_imp[T,N], U_imp2[T,N] (fixed)',
                           order_ok=bool(order_ok), per_instrument_usage=use,
                           note='streams pre-drawn per seed; instruments consume identical shared arrays (CRN); '
                                'v1.2 voluntary consumes ONLY the monitoring stream (verification draw), no new draws',
                           ok=bool(stream_ok))

OUT['status'] = 'PASS' if (OUT['static']['all_shared_identical'] and OUT['static']['config_diff_ok']
                           and OUT['static']['tau_star_config_byte_equal'] and OUT['static']['f_S_ok']
                           and all(v['bitwise_identical'] for v in dyn.values()) and stream_ok) else 'FAIL'
json.dump(OUT, open('results_v12/gate_C4.json', 'w'), indent=1)
print('C4:', OUT['status'])
print(' static shared identical:', OUT['static']['all_shared_identical'],
      '| config diff:', OUT['static']['config_diff_keys'], '| tau byte-equal:', OUT['static']['tau_star_config_byte_equal'])
for k, v in dyn.items(): print(' dynamic', k, '->', v['bitwise_identical'])
print(' streams ok:', stream_ok)
