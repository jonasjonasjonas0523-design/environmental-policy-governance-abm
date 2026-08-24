"""Shared v1.2 harness utilities. The instrumented voluntary runner executes the
model's OWN run_voluntary source (model_v12.py untouched) with per-period
participation masks captured; V17 established that these masks equal the frozen
closed form bitwise."""
import numpy as np, inspect, sys, os
sys.path.insert(0, os.path.dirname(__file__))
import model_v12 as M
from scipy.stats import norm

Z95 = norm.ppf(0.975)

def _make_instrumented():
    src = inspect.getsource(M.run_voluntary)
    src = src.replace("o = _base_out(S); theta_path=[]; s_path=[]; Rp=np.zeros(S)",
                      "o = _base_out(S); theta_path=[]; s_path=[]; Rp=np.zeros(S); pmasks=[]")
    src = src.replace("        z = part & (bat['U_mon'][:,t,:] < p)",
                      "        pmasks.append(part.copy())\n        z = part & (bat['U_mon'][:,t,:] < p)")
    src = src.replace("    o['theta_path']=np.array(theta_path)",
                      "    o['part_masks']=np.array(pmasks)\n    o['theta_path']=np.array(theta_path)")
    ns = dict(np=np, logistic=M.logistic, M_of=M.M_of, _base_out=M._base_out,
              _finish=M._finish, Ann=M.Ann, CONFIG=M.CONFIG)
    exec(src, ns)
    return ns['run_voluntary']

run_voluntary_instr = _make_instrumented()

def delta_grad(bat, g, tau, cfg=None):
    """Per-seed Delta_grad at governance g (V1_2_VALIDATION_MEMO Sec 2.2):
    intensity(Q) = cumulative voluntary abatement of quintile Q / (T x sum E0 in Q);
    quintiles by E0 within each seed's realized population. Returns per-seed arrays."""
    cfg = cfg or M.CONFIG
    o = run_voluntary_instr(bat, g, tau, cfg)
    E0 = bat['E0']; S, N = E0.shape; T = cfg['T']
    nper = o['part_masks'].sum(0)               # [S,N] participating periods per firm
    cumA = 0.2 * E0 * nper                      # per-firm cumulative voluntary abatement
    q = N // 5
    order = np.argsort(E0, axis=1)
    bot = order[:, :q]; top = order[:, -q:]
    r = np.arange(S)[:, None]
    itop = cumA[r, top].sum(1) / (T * E0[r, top].sum(1))
    ibot = cumA[r, bot].sum(1) / (T * E0[r, bot].sum(1))
    return dict(delta=itop - ibot, itop=itop, ibot=ibot,
                A=o['A'], part=o['part'], o=o)

def boot_stats(X, B=2000, seed=11):
    """Identical construction to v1.1 pipeline.boot_stats (max-t over the columns of X)."""
    S, K = X.shape; rng = np.random.default_rng(seed)
    W = rng.multinomial(S, [1 / S] * S, size=B) / S
    mb = W @ X
    mu = X.mean(0); se = mb.std(0, ddof=1)
    se_safe = np.where(se > 0, se, np.inf)
    t = np.abs((mb - mu) / se_safe).max(1)
    crit = float(np.quantile(t, 0.95)) if np.isfinite(t).any() else Z95
    lo = np.where(se > 0, mu - crit * se, mu); hi = np.where(se > 0, mu + crit * se, mu)
    return mu, se, crit, lo, hi
