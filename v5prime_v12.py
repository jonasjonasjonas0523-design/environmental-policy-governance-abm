"""V5-prime per V1_2_VALIDATION_MEMO (executed BEFORE any principal v1.2 voluntary
output is inspected). Scored: Delta_grad 95 pct simultaneous band strictly above zero
at the Dutch central placement (max-t over replicate seeds, B=2000; validation seeds
90000..90152, principal-disjoint). 8 one-at-a-time placement variants reported (sign
flips reportable, never reclassifying the central verdict). Unscored: A_V feasibility
with the a-priori note. Theta-ablation reported as MECHANISM DIAGNOSTIC only."""
import numpy as np, json, os, sys
import pandas as pd
sys.path.insert(0, os.path.dirname(__file__))
from model_v12 import CONFIG, batch
from utils_v12 import delta_grad, boot_stats, run_voluntary_instr

TAU = json.load(open('results/verification/tau_star.json'))['tau_star']
os.makedirs('results_v12/validation', exist_ok=True)
n = json.load(open('results/pilot/replication_choice.json'))['chosen_n']
seeds = list(range(90000, 90000 + n)); b = batch(seeds)

CEN = dict(p=0.75, c=0.85, k=0.05, h=0.95)
BANDS = dict(p=(0.60, 0.90), c=(0.70, 0.95), k=(0.02, 0.12), h=(0.90, 0.98))
conds = [('central', CEN)]
for dim in 'pckh':
    for tag, val in (('lo', BANDS[dim][0]), ('hi', BANDS[dim][1])):
        g = dict(CEN); g[dim] = val
        conds.append((f'{dim}_{tag}', g))

res = {name: delta_grad(b, g, TAU) for name, g in conds}
# Band family: per-condition. The SCORED family (V1_2_VALIDATION_MEMO) contains only
# the central Delta_grad, so its 95 pct simultaneous band is the one-member max-t band;
# each of the 8 reported variants carries its own analogous band. This construction
# reproduces the archived v12_v5prime_report numbers exactly (9/9 at tol 5e-6); the
# joint-9 family was tested and rejected by the reproduction gate (crit 2.40 vs ~1.97).
mu = np.empty(9); se = np.empty(9); lo = np.empty(9); hi = np.empty(9); crits = np.empty(9)
for i, (name, _) in enumerate(conds):
    m1, s1, c1, l1, h1 = boot_stats(res[name]['delta'][:, None], B=CONFIG['B'], seed=11)
    mu[i], se[i], crits[i], lo[i], hi[i] = m1[0], s1[0], c1, l1[0], h1[0]
crit = crits[0]  # scored (central) critical value

rows = []
for i, (name, g) in enumerate(conds):
    r = res[name]
    rows.append(dict(condition=name, p=g['p'], c=g['c'], k=g['k'], h=g['h'],
                     delta_grad=float(mu[i]), lo=float(lo[i]), hi=float(hi[i]),
                     intensity_top=float(r['itop'].mean()), intensity_bot=float(r['ibot'].mean()),
                     A_V=float(r['A'].mean()), part=float(r['part'].mean()),
                     sign=int(np.sign(mu[i])), sign_flip=bool(mu[i] < 0)))
df = pd.DataFrame(rows)
scored_pass = bool(lo[0] > 0)

# unscored feasibility + realized quantile
nui = CONFIG['nu0'] * b['E0'] ** CONFIG['eta']
realized_q = float((0.75 * nui >= 1.0).mean())
# theta-ablation (mechanism diagnostic only)
cfg_ab = dict(CONFIG); cfg_ab['theta0'] = -700.0; cfg_ab['lam'] = 0.0
o_ab = run_voluntary_instr(b, CEN, TAU, cfg_ab)
A_ab = float(o_ab['A'].mean())

df.to_csv('results_v12/validation/v5prime_results.csv', index=False)
summary = dict(scored_verdict='PASS' if scored_pass else 'FAIL',
               central=dict(delta=float(mu[0]), lo=float(lo[0]), hi=float(hi[0])),
               band_family=('per-condition bootstrap max-t; scored family = central only '
                            '(one-member simultaneous band); variants reported with own bands'),
               crit_central=float(crit), crit_by_condition={c[0]: float(v) for c, v in zip(conds, crits)},
               B=CONFIG['B'], seeds=f'90000..{90000+n-1}', n=n,
               sign_flips=[r['condition'] for r in rows if r['sign_flip']],
               unscored_feasibility=dict(A_V_central=float(res['central']['A'].mean()),
                                         quantile_frozen=0.0855, quantile_realized=realized_q),
               theta_ablation_diagnostic=dict(A_V=float(res['central']['A'].mean()), A_V_theta0=A_ab))
json.dump(summary, open('results_v12/validation/v5prime_summary.json', 'w'), indent=1)

with open('results_v12/validation/v5prime_report.md', 'w') as f:
    f.write('# V5-prime (v1.2) held-out validation\n\n')
    f.write(f'SCORED verdict at the Dutch central point: **{summary["scored_verdict"]}** '
            f'(Delta_grad = {mu[0]:.4f}, 95 pct band [{lo[0]:.4f}, {hi[0]:.4f}]).\n\n')
    f.write(df.to_markdown(index=False))
    f.write(f'\n\nUnscored feasibility check: A_V = {res["central"]["A"].mean():.1f} > 0 at the central '
            f'point, as near-certain a priori (P(0.75 nu_i >= 1) = 0.0855 frozen; realized '
            f'{realized_q:.4f}). No validation weight.\n\n')
    f.write(f'Theta-ablation (mechanism diagnostic only): A_V {res["central"]["A"].mean():.1f} -> '
            f'{A_ab:.1f} with the threat channel removed.\n')

print(f'V5-prime SCORED: {summary["scored_verdict"]}  Delta_grad={mu[0]:.6f}  band=[{lo[0]:.6f},{hi[0]:.6f}]  crit_central={crit:.4f}')
print(df[['condition','delta_grad','lo','hi','A_V','part','sign_flip']].to_string(index=False))
print(f'A_V central={res["central"]["A"].mean():.2f}  realized quantile={realized_q:.4f}  ablated A_V={A_ab:.1f}')

# consistency check vs the archived prior-session report
arch = dict(central=(0.118094, 0.11505, 0.121138), p_lo=(0.111112, 0.10812, 0.114104),
            p_hi=(0.123575, 0.120722, 0.126428), c_lo=(0.116135, 0.112945, 0.119326),
            c_hi=(0.119311, 0.116281, 0.122341), k_lo=(0.118499, 0.115468, 0.12153),
            k_hi=(0.117218, 0.11411, 0.120325), h_lo=(0.124283, 0.121211, 0.127356),
            h_hi=(0.116959, 0.114011, 0.119907))
ok = True
for i, (name, _) in enumerate(conds):
    d0, l0, h0 = arch[name]
    m = (abs(mu[i]-d0) < 5e-7) and (abs(lo[i]-l0) < 5e-6) and (abs(hi[i]-h0) < 5e-6)
    ok &= m
    if not m: print('MISMATCH', name, mu[i], lo[i], hi[i], 'vs', arch[name])
print('REPRODUCTION vs archived v12_v5prime_report:', 'EXACT MATCH' if ok else 'MISMATCH')
json.dump(dict(reproduces_archived=bool(ok)), open('results_v12/validation/v5prime_reproduction.json','w'))
