"""N3 (R2 main point 8 / AQ2-8): correlated firm heterogeneity at additional values.
Frozen program ran corr(log gamma, log E0) = -0.3 (scale economies). Added here:
r = +0.3 and r = -0.6, reference plane, per-plane families, same construction as the
frozen corr03 variant (robust_v12.batch_corr; r=0 reduces bitwise to the frozen batch)."""
import numpy as np, json
from common import *
from robust_v12 import batch_corr

base = np.load('revision2026/ref_plane_out.npy'); st0=stage(base); anc=anchor('reference')
s0,_,_,_ = summarize('baseline',st0,anc,base)
rep=dict(baseline=s0)
z=np.load(f'{RES}/robustness/corr03_bands.npz')
rep['corr_m03_frozen']=dict(sign_agreement=float((z['signs']==st0['signs']).mean()),
                            crit_A=float(z['crit'][0]), crit_C=float(z['crit'][1]),
                            supported=int((z['signs']!=0).sum()))
for tag,r_ in (('corr_p03',+0.3),('corr_m06',-0.6)):
    f=f'{OUT}/n3_{tag}_out.npy'
    if not os.path.exists(f):
        bc = batch_corr(seeds, CONFIG, r_)
        ov = plane_out(bc, CONFIG['ref']['k'], CONFIG['ref']['h'], progress=tag); np.save(f,ov)
    ov=np.load(f); stv=stage(ov)
    sv,_,_,_=summarize(tag,stv,anc,ov)
    rep[tag]=dict(summary=sv, sign_agreement=sign_agreement(stv['signs'],st0['signs']))
    print(tag,'sign agreement', round(rep[tag]['sign_agreement'],4),
          'A-leaders', {k:round(v,3) for k,v in sv['A_leader_anchored'].items()},
          'labels', sv['labels'], flush=True)
json.dump(rep, open(f'{OUT}/n3_corr.json','w'), indent=1, default=float)
print('N3 done')
