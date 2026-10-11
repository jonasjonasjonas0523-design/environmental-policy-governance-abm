"""N2 (R2 main point 7 / AQ2-7): reputational-symmetry variant (rep_standard) on the
high-leakage (k=0.55) and long-horizon (h=0.95) slices, where the voluntary arm's
baseline leadership lives, plus the reference plane record from the frozen robustness
run. Baseline and variant both under per-plane families."""
import numpy as np, json, pandas as pd
from common import *

b = batch(seeds)
rep={}
for plane,k,h in (('k055',0.55,CONFIG['ref']['h']), ('h095',CONFIG['ref']['k'],0.95)):
    f_b=f'{OUT}/n2_base_{plane}_out.npy'; f_v=f'{OUT}/n2_repstd_{plane}_out.npy'
    if not os.path.exists(f_b):
        ob = plane_out(b,k,h,progress=f'base {plane}'); np.save(f_b,ob)
    ob=np.load(f_b)
    # reproduction check of the plane against the archived final summary (band-free means)
    arch=pd.read_csv(f'{RES}/principal/plane_{plane}_summary_final.csv')
    mA_arch=arch.pivot_table(index=['p','c'],columns='instrument',values='mean_A').sort_index()
    mA_mine=ob[:,:,:,0].mean(0)
    chk=float(np.abs(mA_arch.values - mA_mine.reshape(Rg*Rg,3)).max())
    if not os.path.exists(f_v):
        ov = plane_out(b,k,h,v=dict(rule='rep_standard'),progress=f'repstd {plane}'); np.save(f_v,ov)
    ov=np.load(f_v)
    st_b=stage(ob); st_v=stage(ov); anc=anchor(plane)
    sb,_,Alb,_=summarize(f'{plane}_base',st_b,anc,ob)
    sv,_,Alv,_=summarize(f'{plane}_repstd',st_v,anc,ov)
    # supported crossover existence on S-V|A and T-V|A (anchored, per-plane family)
    xo={}
    for tag,st in (('base',st_b),('repstd',st_v)):
        for pi,(pn,_,_) in enumerate(pairs):
            posm=(st['lo'][:,pi,0]>0)&anc; negm=(st['hi'][:,pi,0]<0)&anc
            xo[f'{tag}:{pn}|A']=bool(posm.any() and negm.any())
    rep[plane]=dict(plane_reproduction_max_meanA_diff=chk, baseline=sb, variant=sv,
        sign_agreement=sign_agreement(st_v['signs'],st_b['signs']),
        V_A_leader_anchored=dict(baseline=sb['A_leader_anchored']['V'], variant=sv['A_leader_anchored']['V']),
        supported_crossover_A=xo)
    print(plane, 'plane repro max|dmeanA|:', f'{chk:.2e}')
    print('  V anchored A-leader share: base', round(sb['A_leader_anchored']['V'],3),
          '-> repstd', round(sv['A_leader_anchored']['V'],3),
          '| T:', round(sb['A_leader_anchored']['T'],3),'->',round(sv['A_leader_anchored']['T'],3),
          '| sign agreement', round(rep[plane]['sign_agreement'],3))
    print('  crossovers A:', xo, flush=True)
json.dump(rep, open(f'{OUT}/n2_repstd_slices.json','w'), indent=1, default=float)
print('N2 done')
