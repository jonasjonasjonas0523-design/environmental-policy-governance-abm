"""Re-execute validation V1-V4 on the v1.2 build (V1_2_IMPLEMENTATION_SPECIFICATION
Sec 0 step 5). Expected bit-identical to the archived v1.1 rows; any change is a
verification failure. The v1.1 V5 row is preserved as the permanent record (not re-run:
its rule is the v1.1 decision rule, retained in Table 4 per AMENDMENTS_V1_2 Sec 3.3)."""
import numpy as np, json, os, sys
import pandas as pd
sys.path.insert(0, os.path.dirname(__file__))
from model_v12 import CONFIG, batch, run_standard, run_charge
from utils_v12 import boot_stats

CFG = CONFIG
TAU = json.load(open('results/verification/tau_star.json'))['tau_star']
os.makedirs('results_v12/validation', exist_ok=True)

n = json.load(open('results/pilot/replication_choice.json'))['chosen_n']
seeds = list(range(90000, 90000 + n)); b = batch(seeds)
gam, E0 = b['gamma'], b['E0']; f_S = CFG['m_S'] * TAU
V = CFG['val']; rows = []
# V1
c1 = V['V1']
o5 = run_standard(b, dict(p=c1['p0'], c=c1['c'], k=c1['k'], h=c1['h']), TAU)
o9 = run_standard(b, dict(p=c1['p1'], c=c1['c'], k=c1['k'], h=c1['h']), TAU)
X = np.stack([o9['A'] - o5['A'], o9['viol'] - o5['viol']], 1)
mu, se, crit, lo, hi = boot_stats(X)
p1 = bool(lo[0] > 0 and hi[1] < 0)
rows.append(dict(ID='V1', criterion='band(dA)>0 AND band(dviol)<0',
    result=f'dA={mu[0]:.1f} [{lo[0]:.1f},{hi[0]:.1f}]; dviol={mu[1]:.3f} [{lo[1]:.3f},{hi[1]:.3f}]',
    status='PASS' if p1 else 'FAIL', classification='fatal' if not p1 else ''))
# V2
def emis(p):
    q = p * c1['c'] * (1 - c1['k']) * f_S; a = np.minimum(q / gam, 0.2 * E0); return E0 - a
E5, E9 = emis(c1['p0']), emis(c1['p1'])
nc = E5 > 0.8 * E0 + 1e-12
m5 = np.array([E5[s][nc[s]].mean() for s in range(n)])
m9 = np.array([E9[s][nc[s]].mean() for s in range(n)])
eff = (m5 - m9) / m5 * 100 / (c1['p1'] - c1['p0'])
muE, seE, _, loE, hiE = boot_stats(eff[:, None])
p2 = bool(10 <= muE[0] <= 40)
rows.append(dict(ID='V2', criterion='effect in [10,40] pct per inspection; Z reported',
    result=f'{muE[0]:.1f} pct [{loE[0]:.1f},{hiE[0]:.1f}]; Z pending AUTHOR ACTION '
           f'(empirical coeff+SE from ref [35] tables); sim SE={seE[0]:.2f}',
    status='PASS' if p2 else 'FAIL', classification='reportable' if not p2 else ''))
# V3
e, d = V['V3']['eng'], V['V3']['dis']
oe = run_standard(b, dict(zip('pckh', e)), TAU); od = run_standard(b, dict(zip('pckh', d)), TAU)
X = np.stack([oe['A'] - od['A'], od['A']], 1)
mu, se, crit, lo, hi = boot_stats(X)
p3i = bool(lo[0] > 0); p3ii = bool(mu[1] < 0.25 * CFG['A_cum_nominal'])
rows.append(dict(ID='V3', criterion='band(A_eng-A_dis)>0; A_dis<0.25 A*_cum',
    result=f'diff={mu[0]:.1f} [{lo[0]:.1f},{hi[0]:.1f}]; A_dis={mu[1]:.1f} (<100: {p3ii})',
    status='PASS' if (p3i and p3ii) else ('PARTIAL' if p3i else 'FAIL'),
    classification=('' if (p3i and p3ii) else ('reportable(ii)' if p3i else 'fatal(i)'))))
# V4
c4 = V['V4']; o4 = run_charge(b, dict(p=c4['p'], c=c4['c'], k=c4['k'], h=c4['h']), TAU)
X = np.stack([o4['A'], o4['under']], 1); mu, se, crit, lo, hi = boot_stats(X)
p4 = bool(mu[0] >= 0.8 * CFG['A_cum_nominal'] and mu[1] < 0.05)
rows.append(dict(ID='V4', criterion='A>=0.80 A*_cum AND under-reporting<5 pct',
    result=f'A={mu[0]:.1f} [{lo[0]:.1f},{hi[0]:.1f}] (>=320); under={100*mu[1]:.2f} pct',
    status='PASS' if p4 else 'FAIL', classification='fatal' if not p4 else ''))

new = pd.DataFrame(rows)
arch = pd.read_csv('results/validation/validation_results.csv')
cmp = []
for _, r in new.iterrows():
    a = arch[arch.ID == r.ID].iloc[0]
    same = (str(a.criterion) == r.criterion and str(a.result) == r.result and str(a.status) == r.status)
    cmp.append(dict(ID=r.ID, bit_identical=bool(same)))
new.to_csv('results_v12/validation/validation_V1_V4_reexec.csv', index=False)
allok = all(c['bit_identical'] for c in cmp)
json.dump(dict(comparison=cmp, all_bit_identical=allok,
               note='v1.1 V5 row preserved as permanent record (not re-run); V5-prime executed separately'),
          open('results_v12/validation/V1_V4_bit_identity.json', 'w'), indent=1)
print(new[['ID', 'status']].to_string(index=False))
print('bit-identical to archive:', cmp, '->', 'PASS' if allok else 'FAIL')
