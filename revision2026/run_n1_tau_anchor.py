"""N1 (R2 main point 1 / AQ2-6): tau* anchoring statement + alternative anchoring run.
Frozen anchoring: full implementation integrity (p=c=1, k=0, truthful reporting).
Alternative: tau_ref solved so expected abatement at REFERENCE governance
(c,k)=(0.50,0.30), truthful reporting (p=0.50 is the frozen indifference point,
tie -> truthful), equals the common ambition: E[min(c(1-k)tau/gamma, E0)] = 0.2.
Charge arm re-run at tau_ref with f_T = m_T*tau_ref (reporting threshold preserved
at 0.5); standard and voluntary arms UNCHANGED (f_S stays 2*tau* frozen, so the
variant isolates the charge's anchoring). Per-plane family both sides."""
import numpy as np, json
from common import *
from model_v12 import _quad

gq,Eq = _quad(10**6, 424242)
EFF = 0.5*(1-0.30)  # c(1-k) at reference
def F(tau): return float(np.mean(np.minimum(EFF*tau/gq, Eq)) - 0.20)
lo_,hi_ = 1e-5, 2.0
while F(hi_) < 0: hi_ *= 2
for _ in range(200):
    mid=0.5*(lo_+hi_)
    if abs(F(mid))<1e-6: break
    if F(mid)<0: lo_=mid
    else: hi_=mid
tau_ref=mid
print(f'tau_ref = {tau_ref:.10f}  (= {tau_ref/TAU:.4f} x tau*; naive tau*/0.35 = {TAU/0.35:.6f})')

b = batch(seeds)
base = np.load('revision2026/ref_plane_out.npy')
st0 = stage(base); anc = anchor('reference')
s0_summary,_ ,Al0,_ = summarize('baseline_reference', st0, anc, base)

outv = plane_out(b, CONFIG['ref']['k'], CONFIG['ref']['h'], T_tau=tau_ref, progress='tau_ref')
stv = stage(outv)
sv_summary, labv, Alv, Clv = summarize('tau_ref_anchored', stv, anc, outv)

# charge delivery at the reference cell and across anchored cells
icr = np.where((np.abs(np.repeat(pp,Rg)-0.5)<1e-9)&(np.abs(np.tile(pp,Rg)-0.5)<1e-9))[0][0]
rep=dict(tau_ref=tau_ref, ratio=tau_ref/TAU,
    anchoring_frozen='full integrity p=c=1,k=0, truthful (quad seed 424242, n=1e6)',
    anchoring_alt='reference governance (c,k)=(0.5,0.30), truthful at the p=0.5 indifference point',
    f_T_rule='m_T x tau_ref (reporting threshold preserved at 0.5)', f_S='unchanged 2 x tau* (S and V arms untouched)',
    baseline=s0_summary, variant=sv_summary,
    sign_agreement=sign_agreement(stv['signs'], st0['signs']),
    T_delivery_ref_cell=dict(baseline=float(base[:,icr,1,0].mean()), variant=float(outv[:,icr,1,0].mean()),
                             ambition=400.0),
    anchored_A_leader_T=dict(baseline=s0_summary['A_leader_anchored']['T'], variant=sv_summary['A_leader_anchored']['T']),
    anchored_C_leader_S=dict(baseline=s0_summary['C_leader_anchored']['S'], variant=sv_summary['C_leader_anchored']['S']))
np.savez(f'{OUT}/n1_tau_ref_bands.npz', **stv, out_mean=outv.mean(0))
json.dump(rep, open(f'{OUT}/n1_tau_anchor.json','w'), indent=1, default=float)
print(json.dumps({k:rep[k] for k in ('tau_ref','ratio','sign_agreement','T_delivery_ref_cell','anchored_A_leader_T','anchored_C_leader_S')}, indent=1, default=float))
print('baseline labels:', s0_summary['labels']); print('variant labels:', sv_summary['labels'])
print('baseline A-leaders:', s0_summary['A_leader_anchored']); print('variant A-leaders:', sv_summary['A_leader_anchored'])
