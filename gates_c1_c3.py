"""Gates C1-C3 per V1_2_CODE_CHECKS.md, executed on the UNMODIFIED v1.1 build
before any v1.2 model modification. Function-level / unit-trace only (no ABM run
beyond the specified traces)."""
import numpy as np, sys, os, json
sys.path.insert(0, os.path.dirname(__file__))
from model import CONFIG, M_of, logistic, seed_data, run_voluntary, batch

TAU = 0.09326697147369385
F_S = CONFIG['m_S'] * TAU
OUT = {}

# ---------------- C1: v1.1 robustness-corner activation check ----------------
# (a) one-period M, nu=0.75, theta forced 0.909, complier firm (m/K = 1),
#     h in {1.0, 0.95}, p swept 0..1 by 0.01. v1.1 rule: part iff K <= p*(nu*K + delta*h*theta*M)
#     with M/K = 1 -> switch at p >= 1/(nu + delta*h*theta).
nu, th, dl = 0.75, 0.909, CONFIG['delta']
res_a = {}
for h, expect in ((1.0, 0.6198), (0.95, 0.6368)):
    ps = np.round(np.arange(0, 1.0001, 0.01), 2)
    part = ps * (nu + dl * h * th * 1.0) >= 1.0          # v1.1 rule, ties -> participate (<= form)
    p_on = float(ps[part][0]) if part.any() else None
    ok = p_on is not None and abs(p_on - expect) <= 0.01 + 1e-12
    res_a[f'h={h}'] = dict(first_p_on=p_on, expected=expect, within_one_grid_step=bool(ok))
# (b) v1.1 baseline (nu=0.25, theta<=0.6900): False on full 5^4 probe grid, 1000 probe firms
rng = np.random.default_rng(12345)                        # probe-firm seed (implementation detail, logged)
gam = np.exp(rng.normal(CONFIG['mu_g'], CONFIG['sg'], 1000))
E0  = np.exp(rng.normal(CONFIG['mu_E'], CONFIG['sE'], 1000))
K   = 0.5 * gam * (0.2 * E0)**2
grid = [0.0, 0.25, 0.5, 0.75, 1.0]
theta_cap = 0.6900
any_part = False; n_pts = 0
for p in grid:
    for c in grid:
        for k in grid:
            M, _, _ = M_of(p, c, k, F_S, gam, E0)
            for h in grid:
                part = K <= p * (0.25 * K + dl * h * theta_cap * M)   # v1.1 rule at theta ceiling
                any_part |= bool(part.any()); n_pts += 1
OUT['C1'] = dict(activation=res_a, baseline_probe=dict(points=n_pts, firms=1000,
                 any_participation=bool(any_part), probe_seed=12345),
                 status='PASS' if (all(v['within_one_grid_step'] for v in res_a.values()) and not any_part) else 'FAIL')

# ---------------- C2: theta timing (code inspection refs + 2-period unit trace) --------
# model.py: theta initialized at line `theta = np.full(S, logistic(cfg['theta0']))` (theta_1, P2);
# inside the t-loop the decision uses the CURRENT theta (built from s_{t-1});
# the update `theta = logistic(theta0 + lam*s)` happens AFTER the decision of period t.
b1 = batch([31337]); o = run_voluntary(b1, dict(p=.8, c=.85, k=.05, h=.95), TAU,
                                       {**CONFIG, 'nu': 1.5})   # live channel so s_t is informative
tp, sp = o['theta_path'], o['s_path']                            # [T,S]
ok2 = (abs(float(tp[0, 0]) - logistic(CONFIG['theta0'])) < 1e-15 and
       abs(float(tp[1, 0]) - logistic(CONFIG['theta0'] + CONFIG['lam'] * float(sp[0, 0]))) < 1e-15)
OUT['C2'] = dict(theta1=float(tp[0, 0]), expected_theta1=logistic(CONFIG['theta0']),
                 theta2=float(tp[1, 0]), expected_theta2=logistic(CONFIG['theta0'] + CONFIG['lam'] * float(sp[0, 0])),
                 code_refs='model.py run_voluntary: init logistic(theta0) pre-loop; decision uses theta before update; update post-decision',
                 status='PASS' if ok2 else 'FAIL')

# ---------------- C3: exemption bookkeeping (3-period trace, stateless) ----------------
# Frozen behavior: exemption exists only inside the current period's belief term;
# no regime switch in-run; no persistent exemption state. Trace: force alternating
# verified/unverified achievement via crafted U_mon; verify decisions each period
# depend only on (theta_t, current-period comparison), never on past z.
cfg3 = {**CONFIG, 'nu': 1.5}
d = seed_data(777)
g3 = dict(p=.8, c=.85, k=.05, h=.95)
o3 = run_voluntary({k_: v[None] if v.ndim == 1 else v[None] for k_, v in d.items()}, g3, TAU, cfg3)
# independent 3-period re-derivation with NO exemption state:
gam3, E03 = d['gamma'], d['E0']
M3, _, _ = M_of(g3['p'], g3['c'], g3['k'], F_S, gam3, E03)
Ct3 = 0.5 * gam3 * (0.2 * E03)**2; rho3 = 1.5 * Ct3
Ast = CONFIG['A_star_frac'] * E03.sum()
th3 = logistic(CONFIG['theta0']); ok3 = True
for t in range(3):
    part = Ct3 <= g3['p'] * (rho3 + dl * g3['h'] * th3 * M3)      # depends ONLY on current theta
    z = part & (d['U_mon'][t] < g3['p'])
    s = max(0.0, (Ast - (0.2 * E03 * z).sum()) / Ast)
    # cross-check against the model's realized paths (statelessness: identical decisions)
    ok3 &= abs(float(o3['theta_path'][t, 0]) - th3) < 1e-15
    ok3 &= abs(float(o3['s_path'][t, 0]) - s) < 1e-12
    th3 = logistic(CONFIG['theta0'] + CONFIG['lam'] * s)
# static inspection: no exemption flag persists (no state variable carried across t in run_voluntary
# other than theta and accumulators; verified by reading the function source)
import inspect, re
src = inspect.getsource(run_voluntary)
code_only = '\n'.join(ln.split('#', 1)[0] for ln in src.splitlines())   # strip comments
no_flag = not re.search(r'\bexempt\w*\s*=|\bregime\w*\s*=', code_only, re.I)
OUT['C3'] = dict(three_period_trace_matches=bool(ok3), no_persistent_exemption_state=bool(no_flag),
                 status='PASS' if (ok3 and no_flag) else 'FAIL')

json.dump(OUT, open('results_v12/gates_C1_C3.json', 'w'), indent=1)
print(json.dumps({k: v['status'] for k, v in OUT.items()}, indent=0))
print('C1 activation:', OUT['C1']['activation'])
print('C1 baseline probe any participation:', OUT['C1']['baseline_probe']['any_participation'])
