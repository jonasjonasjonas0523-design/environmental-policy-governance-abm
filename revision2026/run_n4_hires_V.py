"""N4 (R2 main point 6 / AQ2-9): extend the archived 81x81 S/T boundary check to the
voluntary arm, so every H3-carrying pair (S-V|A, T-V|A) has a coarse-to-fine record.
Same construction as the archived check: principal seeds, pilot nearest-matched
crossing metric (tol = 0.5 x finer spacing), plus the labeled post-hoc shared-point
2SE supported-sign diagnostic at the 1681 shared grid points."""
import numpy as np, json, time
from common import *
from pipeline import all_crossings, cross_disp

R2g=81; pp2=np.linspace(0,1,R2g)
assert np.allclose(pp2[::2], pp)
b = batch(seeds)
f=f'{OUT}/n4_hires81_out.npy'
if not os.path.exists(f):
    out2=np.zeros((n,R2g*R2g,3,2)); t0=time.time()
    for ic in range(R2g*R2g):
        g=dict(p=pp2[ic//R2g], c=pp2[ic%R2g], k=CONFIG['ref']['k'], h=CONFIG['ref']['h'])
        r=run_all(b,g,TAU)
        for ji,j in enumerate('STV'):
            out2[:,ic,ji,0]=r[j]['A']; out2[:,ic,ji,1]=r[j]['C_social']
        if ic%1200==0: print(f'  81x81 cell {ic}/{R2g*R2g} ({time.time()-t0:.0f}s)', flush=True)
    np.save(f,out2); print(f'81x81 computed in {time.time()-t0:.0f}s')
out2=np.load(f)
out1=np.load('revision2026/ref_plane_out.npy')
D1=np.stack([out1[:,:,i,:]-out1[:,:,j,:] for _,i,j in pairs],2)  # [n,1681,3,2]
D2=np.stack([out2[:,:,i,:]-out2[:,:,j,:] for _,i,j in pairs],2)
tol=0.5/(R2g-1)
rep=dict(tol=tol)
for pi,(pn,_,_) in enumerate(pairs):
    for mi,m in enumerate(('A','C_social')):
        m1=D1[:,:,pi,mi].mean(0).reshape(Rg,Rg); s1=(D1[:,:,pi,mi].std(0,ddof=1)/np.sqrt(n)).reshape(Rg,Rg)
        m2=D2[:,:,pi,mi].mean(0).reshape(R2g,R2g); s2=(D2[:,:,pi,mi].std(0,ddof=1)/np.sqrt(n)).reshape(R2g,R2g)
        cr1=all_crossings(m1,s1,pp); cr2=all_crossings(m2,s2,pp2)
        disp,matched,mismatch=cross_disp(cr1,cr2)
        # shared-point pointwise 2SE sign diagnostic
        sh2=out2.reshape(n,R2g,R2g,3,2)[:, ::2, ::2].reshape(n,Rg*Rg,3,2)
        dD=np.stack([sh2[:,:,i,:]-sh2[:,:,j,:] for _,i,j in pairs],2)[:,:,pi,mi]
        d1=D1[:,:,pi,mi]
        def psign(X):
            mu=X.mean(0); se=X.std(0,ddof=1)/np.sqrt(n)
            return np.where(mu-2*se>0,1,np.where(mu+2*se<0,-1,0))
        g1,g2=psign(d1),psign(dD)
        both=(g1!=0)&(g2!=0)
        agree=float((g1[both]==g2[both]).mean()) if both.any() else 1.0
        contra=int(((g1*g2)==-1).sum())
        rep[f'{pn}|{m}']=dict(max_crossing_shift=disp, matched=matched, mismatched=mismatch,
            within_tol=bool(disp is not None and disp<tol),
            shared_sign_agreement=agree, sign_contradictions=contra,
            supported_41=int((g1!=0).sum()), supported_81=int((g2!=0).sum()))
        print(f'{pn}|{m}: shift={disp if disp is None else round(disp,5)} (tol {tol:.5f}) matched={matched} mm={mismatch} '
              f'| shared-sign agree={agree:.4f} contradictions={contra}', flush=True)
json.dump(rep, open(f'{OUT}/n4_hires_V.json','w'), indent=1, default=float)
print('N4 done')
