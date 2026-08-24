# v1_2_calibration_script.py -- FROZEN v1.2 calibration (audit-corrected, 2026-08-21)
# Deterministic mean-field population algebra. NOT an ABM simulation.
# Audit corrections embodied here:
#  (1) tau* is the frozen v1.1 constant (seed 424242, n_quad=1e6) imported below and
#      NEVER re-solved for use; an independent re-solve is printed as a consistency
#      diagnostic only (tolerance 5e-4).
#  (2) Calibration statistic = EVER-participated share over t=1..T (construct match
#      to 33/50 membership counts); the previous time-average annual statistic is
#      retired (construct mismatch; see V1_2_CALIBRATION_MEMO Sec 3).
#  (3) The ONLY Dutch-placement quantity computed anywhere is the static quantile
#      P(0.75*nu_i >= 1) ordered by the external audit for the V5p-i assessment.
#      No Dutch dynamics, no archetype evaluation, by construction (gate C5).
# Protocol: screening at 100k draws (first 100k of the seeded 400k), ladder from the
# central placement outward in 0.01 steps within the frozen band; first admissible
# rung is confirmed at 400k with the fine grid; confirmation governs adoption.
import numpy as np
TAU_STAR=0.09326697147369385   # v1.1 frozen. DO NOT RE-SOLVE FOR USE.
FS=2.0*TAU_STAR                 # 0.1865339429473877
SEED=20260821; NF=400_000; NS=100_000
delta,T=0.95,10; theta0,lam=-2.2,3.0
c_us,k_us,h_us=0.85,0.05,0.95
P_CENTRAL,P_LO,P_HI=0.35,0.20,0.50
TIER=0.075; S_ALL_T,S_TOP_T=0.15,0.64
TOL_ALL,TOL_TOP=0.03,0.08
def sig(x): return 1.0/(1.0+np.exp(-x))
r=np.random.default_rng(SEED)
gF=r.lognormal(-0.5,1.0,NF); eF=r.lognormal(-0.5,1.0,NF)
gS,eS=gF[:NS],eF[:NS]
topF=eF>=np.quantile(eF,1-TIER); topS=eS>=np.quantile(eS,1-TIER)
dh=delta*h_us
def ever_stats(nu0,eta,p,g,e0,top):
    q=p*c_us*(1-k_us)*FS; u=g*e0; v=np.minimum(5.0*q/u,1.0); mK=2.0*v-v*v
    nui=nu0*np.power(e0,eta)
    ever=np.zeros(g.shape[0],dtype=bool); th=sig(theta0); traj=[]
    for t in range(1,T+1):
        n=T-t; Ann=(1.0-dh**n)/(1.0-dh) if n>0 else 0.0
        P=p*(nui+dh*th*mK*Ann)>=1.0
        ever|=P; traj.append(float(P.mean()))
        w=e0[P].sum()/e0.sum() if P.any() else 0.0
        th=sig(theta0+lam*max(0.0,1.0-p*w))
    return float(ever.mean()),float(ever[top].mean()),traj
def logit(x):
    x=min(max(x,1e-6),1-1e-6); return np.log(x/(1-x))
def search(p,g,e0,top,stt=S_TOP_T,grid=(16,11)):
    best=(1e18,None)
    for eta in np.linspace(0.0,2.0,grid[1]):
        for nu0 in np.linspace(0.0,3.0,grid[0]):
            sa,st,_=ever_stats(nu0,eta,p,g,e0,top)
            d=(logit(sa)-logit(S_ALL_T))**2+(logit(st)-logit(stt))**2
            if d<best[0]: best=(d,(nu0,eta,sa,st))
    nu0,eta=best[1][0],best[1][1]
    for span in (0.12,0.03,0.008):
        cand=best
        for dn in np.linspace(-4*span,4*span,9):
            for de in np.linspace(-4*span,4*span,9):
                n2,e2=max(0.0,nu0+dn),min(2.0,max(0.0,eta+de))
                sa,st,_=ever_stats(n2,e2,p,g,e0,top)
                d=(logit(sa)-logit(S_ALL_T))**2+(logit(st)-logit(stt))**2
                if d<cand[0]: cand=(d,(n2,e2,sa,st))
        best=cand; nu0,eta=best[1][0],best[1][1]
    return best
def ladder():
    steps=[P_CENTRAL]; d=0.01
    while True:
        lo2,hi2=round(P_CENTRAL-d,2),round(P_CENTRAL+d,2); added=False
        if lo2>=P_LO-1e-9: steps.append(lo2); added=True
        if hi2<=P_HI+1e-9: steps.append(hi2); added=True
        if not added: break
        d+=0.01
    return steps
if __name__=="__main__":
    lo,hi=0.01,0.5
    for _ in range(60):
        m=0.5*(lo+hi)
        if np.minimum(m/gF,eF).mean()<0.2: lo=m
        else: hi=m
    print(f"[diag] independent tau* re-solve = {0.5*(lo+hi):.6f}; frozen = {TAU_STAR:.17f}; |diff| = {abs(0.5*(lo+hi)-TAU_STAR):.2e} (<= 5e-4 required)")
    print(f"[const] fS = {FS:.16f}; Ann(0.9025,8) = {(1-0.9025**8)/0.0975:.6f}; Ann(0.9025,9) = {(1-0.9025**9)/0.0975:.6f}")
    print(f"[diag] threat-ignition p (theta ceiling {sig(theta0+lam):.4f}, max Ann {(1-0.9025**8)/0.0975:.4f}): p >= {1.0/(0.9025*((1-0.9025**8)/0.0975)*sig(theta0+lam)):.4f}")
    adopted=None
    for p in ladder():
        d,(nu0,eta,sa,st)=search(p,gS,eS,topS)
        ok=(abs(sa-S_ALL_T)<=TOL_ALL) and (abs(st-S_TOP_T)<=TOL_TOP)
        print(f"screen p={p:.2f}: EVER_all={sa:.4f} EVER_top={st:.4f} admissible={ok}")
        if ok:
            dF,(nu0F,etaF,saF,stF)=search(p,gF,eF,topF,grid=(31,21))
            okF=(abs(saF-S_ALL_T)<=TOL_ALL) and (abs(stF-S_TOP_T)<=TOL_TOP)
            print(f"CONFIRM 400k p={p:.2f}: nu0={nu0F:.4f} eta={etaF:.4f} EVER_all={saF:.4f} EVER_top={stF:.4f} admissible={okF}")
            if okF: adopted=(p,nu0F,etaF,saF,stF); break
    assert adopted, "no admissible rung: calibration fails; return to decision stage"
    p_ad,nu0f,etaf,saf,stf=adopted
    sa_f,st_f,traj=ever_stats(nu0f,etaf,p_ad,gF,eF,topF)
    print(f"ADOPTED: p_US={p_ad:.2f} nu_0={nu0f:.4f} eta={etaf:.4f} EVER_all={sa_f:.4f} EVER_top={st_f:.4f}")
    print(f"path t=1..10: {[round(x,3) for x in traj]}")
    nui=nu0f*np.power(eF,etaf)
    print(f"nu_i: median={np.median(nui):.3f} p90={np.quantile(nui,0.90):.3f} p99={np.quantile(nui,0.99):.3f} share>1={float((nui>1).mean()):.4f}")
    topD=eF>=np.quantile(eF,1-0.10)
    dD,(nD,eD,sAD,sTD)=search(p_ad,gF,eF,topD,grid=(31,21))
    print(f"[unc decile] TIER=0.10 at p={p_ad:.2f}: nu0={nD:.4f} eta={eD:.4f} EVER_all={sAD:.4f} EVER_top10={sTD:.4f} (bound solution if eta=2.0)")
    for tag,pp in (("corner p=0.20",0.20),("corner p=0.50",0.50)):
        dd,(n2,e2,sA,sT)=search(pp,gS,eS,topS)
        med=np.median(n2*np.power(eS,e2))
        print(f"[unc] {tag}: nu0={n2:.4f} eta={e2:.4f} EVER_all={sA:.4f} EVER_top={sT:.4f} median_nu={med:.3f}")
    print(f"[V5p audit] P(0.75*nu_i >= 1) = {float((0.75*nui>=1.0).mean()):.4f}  (static quantile only; no Dutch dynamics computed)")
