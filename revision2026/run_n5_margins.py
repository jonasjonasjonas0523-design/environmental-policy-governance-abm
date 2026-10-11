"""N5 (R2 main point 2 + AQ2-2/AQ2-3): failure-margin analysis for the scored tests and
the V2 parameter-box pass fraction. Validation seeds 90000..90152 (frozen). The V2
statistic is the validate_v12 construction verbatim, parametrized over the predeclared
robustness ranges that enter it (m_S in [1,4]; sigma_g, sigma_E in [0.5,1.5], mean-1
recentered per the frozen F4 convention). Uniform weighting over the box, stated as a
convention. Deterministic population algebra; no new model claims."""
import numpy as np, json
from common import *
nv=153; vseeds=list(range(90000,90000+nv))
c1=CONFIG['val']['V1']

def vbatch(sg,sE):
    cfg=dict(CONFIG); cfg.update(sg=sg, mu_g=-sg*sg/2, sE=sE, mu_E=-sE*sE/2)
    return batch(vseeds,cfg), cfg

def v2_stat(mS, sg=1.0, sE=1.0, bcache={}):
    key=(sg,sE)
    if key not in bcache: bcache[key]=vbatch(sg,sE)
    bb,cfg=bcache[key]
    gam,E0=bb['gamma'],bb['E0']; f_S=mS*TAU
    def emis(p):
        q=p*c1['c']*(1-c1['k'])*f_S; a=np.minimum(q/gam,0.2*E0); return E0-a
    E5,E9=emis(c1['p0']),emis(c1['p1'])
    nc=E5>0.8*E0+1e-12
    m5=np.array([E5[s][nc[s]].mean() for s in range(nv)])
    m9=np.array([E9[s][nc[s]].mean() for s in range(nv)])
    return float(((m5-m9)/m5*100/(c1['p1']-c1['p0'])).mean())

base=v2_stat(2.0)
print(f'V2 baseline reproduction: {base:.6f} (frozen 10.026535)')
assert abs(base-10.026535)<5e-4
# m_S crossing points
lo_,hi_=1.0,2.0
for _ in range(60):
    mid=0.5*(lo_+hi_)
    if v2_stat(mid)<10: lo_=mid
    else: hi_=mid
mS_fail_low=0.5*(lo_+hi_)
hi2=4.0; v40=v2_stat(4.0)
rep=dict(V2=dict(baseline=base, band_frozen=[9.8542,10.1989],
    mS_at_lower_bound=mS_fail_low, pct_drop_in_mS_to_fail=100*(2.0-mS_fail_low)/2.0,
    stat_at_mS4=v40))
# box fraction
mSg=np.linspace(1,4,61); sgg=[0.5,0.75,1.0,1.25,1.5]; sEg=[0.5,0.75,1.0,1.25,1.5]
vals=np.array([[[v2_stat(m,sg,sE) for sE in sEg] for sg in sgg] for m in mSg])
inside=(vals>=10)&(vals<=40)
rep['V2']['box']=dict(fraction_inside=float(inside.mean()), vmin=float(vals.min()), vmax=float(vals.max()),
    frac_below_10=float((vals<10).mean()), frac_above_40=float((vals>40).mean()),
    grid='m_S 61x linspace(1,4); sigma_g, sigma_E in {0.5,0.75,1.0,1.25,1.5}; uniform convention')
print('V2 box: fraction inside [10,40] =', round(rep['V2']['box']['fraction_inside'],4),
      'range', round(vals.min(),2), '-', round(vals.max(),2),
      '| m_S lower-bound crossing at', round(mS_fail_low,4))

# ---- V3: engaged-minus-disengaged + level criterion ----
from utils_v12 import boot_stats
bb,_=vbatch(1.0,1.0)
e,d=CONFIG['val']['V3']['eng'],CONFIG['val']['V3']['dis']
oe=run_standard(bb,dict(zip('pckh',e)),TAU); od=run_standard(bb,dict(zip('pckh',d)),TAU)
X=np.stack([oe['A']-od['A'],od['A']],1); mu,se,crit,lo,hi=boot_stats(X)
def Adis_of_mS(mS):
    q=d[0]*d[1]*(1-d[2])*mS*TAU; a=np.minimum(q/bb['gamma'],0.2*bb['E0'])
    return float(CONFIG['T']*a.sum(1).mean())
lo_,hi_=4.0,40.0
while Adis_of_mS(hi_)<100: hi_*=1.5
for _ in range(60):
    mid=0.5*(lo_+hi_)
    if Adis_of_mS(mid)<100: lo_=mid
    else: hi_=mid
rep['V3']=dict(diff=float(mu[0]), diff_lo=float(lo[0]), diff_t=float(mu[0]/se[0]),
    A_dis=float(mu[1]), threshold=100.0, margin_ratio=float(100/mu[1]),
    mS_to_fail_level=0.5*(lo_+hi_),
    note='sign component monotone by construction under paired CRN seeds')
print('V3: diff t-stat', round(mu[0]/se[0],1), '| A_dis', round(mu[1],1),'vs 100; m_S to fail level criterion:', round(rep['V3']['mS_to_fail_level'],1))

# ---- V4: delivery and evasion margins over the population box ----
c4=CONFIG['val']['V4']
o4=run_charge(bb,dict(p=c4['p'],c=c4['c'],k=c4['k'],h=c4['h']),TAU)
Amin_box=1e9; worst=None
for sg in sgg:
    for sE in sEg:
        b2,_=vbatch(sg,sE)
        o=run_charge(b2,dict(p=c4['p'],c=c4['c'],k=c4['k'],h=c4['h']),TAU)
        a=float(o['A'].mean())
        if a<Amin_box: Amin_box=a; worst=(sg,sE)
rep['V4']=dict(A=float(o4['A'].mean()), threshold=0.8*400, under=float(o4['under'].mean()),
    A_min_over_box=Amin_box, box_worst=worst,
    evasion_bound='under-reporting at p=0.95 requires m_T < 1/19 = 0.053; predeclared floor 0.5',
    fails_in_box=bool(Amin_box<320))
print('V4: A =', round(rep['V4']['A'],1), 'vs 320; min over population box', round(Amin_box,1), 'at', worst)

# ---- V1 margins ----
o5=run_standard(bb,dict(p=c1['p0'],c=c1['c'],k=c1['k'],h=c1['h']),TAU)
o9=run_standard(bb,dict(p=c1['p1'],c=c1['c'],k=c1['k'],h=c1['h']),TAU)
X=np.stack([o9['A']-o5['A'],o9['viol']-o5['viol']],1); mu,se,crit,lo,hi=boot_stats(X)
rep['V1']=dict(dA=float(mu[0]), dA_t=float(mu[0]/se[0]), dviol=float(mu[1]), dviol_t=float(mu[1]/se[1]),
    note='both components monotone in p by construction; margins reported in t units')
print('V1: dA t', round(mu[0]/se[0],1), '| dviol t', round(mu[1]/se[1],1))
json.dump(rep, open(f'{OUT}/n5_margins.json','w'), indent=1, default=float)
print('N5 done')
