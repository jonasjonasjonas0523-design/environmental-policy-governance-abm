"""N6 (R2 points 3/4 context): channel decomposition of the V5-prime size gradient at
the frozen Dutch placement (reputation-only and threat-only Delta_grad), quantifying
whether the scored sign was preordained. N7 (R2 point 9): tally of the preregistered
directional expectations E1-E4 on the archived/recomputed surfaces."""
import numpy as np, json, pandas as pd
from common import *
from utils_v12 import delta_grad

nv=153; vseeds=list(range(90000,90000+nv)); bv=batch(vseeds)
CEN=dict(p=0.75,c=0.85,k=0.05,h=0.95)
full=delta_grad(bv,CEN,TAU)
cfg_rep=dict(CONFIG); cfg_rep['theta0']=-700.0; cfg_rep['lam']=0.0     # threat off
cfg_thr=dict(CONFIG); cfg_thr['nu0']=0.0                               # reputation off
rep_only=delta_grad(bv,CEN,TAU,cfg_rep)
thr_only=delta_grad(bv,CEN,TAU,cfg_thr)
n6=dict(full=dict(delta=float(full['delta'].mean()), itop=float(full['itop'].mean()),
                  ibot=float(full['ibot'].mean()), A=float(full['A'].mean())),
        reputation_only=dict(delta=float(rep_only['delta'].mean()), itop=float(rep_only['itop'].mean()),
                  ibot=float(rep_only['ibot'].mean()), A=float(rep_only['A'].mean())),
        threat_only=dict(delta=float(thr_only['delta'].mean()), itop=float(thr_only['itop'].mean()),
                  ibot=float(thr_only['ibot'].mean()), A=float(thr_only['A'].mean())))
json.dump(n6, open(f'{OUT}/n6_grad_decomposition.json','w'), indent=1)
print('N6 Delta_grad: full', round(n6['full']['delta'],4),
      '| reputation-only', round(n6['reputation_only']['delta'],4),
      '| threat-only', round(n6['threat_only']['delta'],4))
print('   intensities top/bot: full', round(n6['full']['itop'],3), round(n6['full']['ibot'],3),
      '| rep-only', round(n6['reputation_only']['itop'],3), round(n6['reputation_only']['ibot'],3),
      '| thr-only', round(n6['threat_only']['itop'],3), round(n6['threat_only']['ibot'],3))

# ---- N7: E1-E4 ----
b=batch(seeds)
planes={}
for nm,k,h in (('reference',0.30,0.75),('k010',0.10,0.75),('k055',0.55,0.75),('h050',0.30,0.50),('h095',0.30,0.95)):
    f=f'{OUT}/n2_base_{nm}_out.npy' if nm in ('k055','h095') else f'{OUT}/n7_base_{nm}_out.npy'
    if nm=='reference': f='revision2026/ref_plane_out.npy'
    if not os.path.exists(f):
        o=plane_out(b,k,h,progress=nm); np.save(f,o)
    planes[nm]=np.load(f).mean(0)[:,:,0].reshape(Rg,Rg,3)   # mean A [p,c,instr]
def mono_viol(M, axis):
    d=np.diff(M,axis=axis)
    return float(-d.min()) if d.size else 0.0
tol=1e-9
e1=dict(S_p=mono_viol(planes['reference'][:,:,0],0), S_c=mono_viol(planes['reference'][:,:,0],1))
e2=dict(T_c=mono_viol(planes['reference'][:,:,1],1),
        T_k_010_vs_ref=float((planes['k010'][:,:,1]-planes['reference'][:,:,1]).min()),
        T_k_ref_vs_055=float((planes['reference'][:,:,1]-planes['k055'][:,:,1]).min()))
e3=dict(V_p=mono_viol(planes['reference'][:,:,2],0),
        V_h_095_vs_ref=float((planes['h095'][:,:,2]-planes['reference'][:,:,2]).min()),
        V_h_ref_vs_050=float((planes['reference'][:,:,2]-planes['h050'][:,:,2]).min()))
rA1=run_all(b,dict(p=.80,c=.85,k=.05,h=.95),TAU)
e4={j: float(rA1[j]['A'].mean()) for j in 'STV'}
n7=dict(E1=dict(met=bool(max(e1.values())<=tol), detail=e1),
        E2=dict(met=bool(e2['T_c']<=tol and e2['T_k_010_vs_ref']>=-tol and e2['T_k_ref_vs_055']>=-tol), detail=e2),
        E3=dict(met=bool(e3['V_p']<=tol and e3['V_h_095_vs_ref']>=-tol and e3['V_h_ref_vs_050']>=-tol), detail=e3),
        E4=dict(met=bool(min(e4.values())>200.0), detail=e4, rule='all three > 0.5 x 400 at A1 central'))
json.dump(n7, open(f'{OUT}/n7_expectations.json','w'), indent=1, default=float)
for k_,v_ in n7.items(): print(f'N7 {k_}: met={v_["met"]}  {v_["detail"]}')
print('N6+N7 done')
