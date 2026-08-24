"""Archetype probabilistic readouts (manuscript Sec V-D protocol, VI-E reporting).
Draws come from the SAME seed-2026 joint copula stream that defines the frozen
anchor masks (per-archetype block: standard_normal(1e5 x 4) then choice(1e5,3000)),
taking the FIRST NDRAW draws of each archetype's block -- an ex-ante instantiation
documented for the audit. Per archetype and metric: P(instrument leads), margin
over runner-up conditional on each draw, probability of overlapping the nearest
crossover band (nearest plane in (k,h); nearest grid cell; transition cells of
FINAL joint bands restricted to pair|metric with a supported within-plane
crossover), placement-sensitivity statement, and local gradients of each pairwise
difference Delta_jl on (p,c,k,h) by least squares over the placement region.
"""
import numpy as np, json, os, sys, time
import pandas as pd
sys.path.insert(0, os.path.dirname(__file__))
from model_v12 import CONFIG, batch, run_all
from scipy.stats import beta as beta_dist, norm

RES='results_v12'; CFG=CONFIG
TAU = json.load(open(f'{RES}/verification/tau_star.json'))['tau_star']
os.makedirs(f'{RES}/archetypes', exist_ok=True)
n  = json.load(open('results/pilot/replication_choice.json'))['chosen_n']
seeds = list(range(100000,100000+n))
NDRAW = 500
pairs = [('S-T',0,1),('S-V',0,2),('T-V',1,2)]
PLANES = [('reference', CFG['ref']['k'], CFG['ref']['h']),
          ('k010', 0.10, CFG['ref']['h']), ('k055', 0.55, CFG['ref']['h']),
          ('h050', CFG['ref']['k'], 0.50), ('h095', CFG['ref']['k'], 0.95)]
Rg = 41; pp = np.linspace(0,1,Rg)

def archetype_draws():
    """Replicates anchor_masks' rng(2026) consumption; returns first NDRAW joint draws per archetype."""
    L = np.linalg.cholesky(np.array(CFG['copula_R']))
    rng = np.random.default_rng(2026)
    out={}
    for name,ab in CFG['archetypes'].items():
        z = rng.standard_normal((100000,4))@L.T; u = norm.cdf(z)
        draws = np.column_stack([beta_dist.ppf(u[:,i],*ab[dim])
                                 for i,dim in enumerate(('p','c','k','h'))])
        rng.choice(100000,3000,replace=False)   # consume, exactly as anchor_masks
        out[name]=draws[:NDRAW]
    return out

def transition_lookup():
    """plane -> dict[(pair,metric)] -> (transition_mask [41,41], supported flag)."""
    sup = pd.read_csv(f'{RES}/principal/crossover_summary_final_all_planes.csv')
    look={}
    for pl,_,_ in PLANES:
        bl = pd.read_parquet(f'{RES}/principal/plane_{pl}_bands_final_joint.parquet')
        d={}
        for pn,_,_ in pairs:
            for m in ('A','C_social'):
                r = bl[(bl['pair']==pn)&(bl['metric']==m)].sort_values('cell')
                tr = ((r['lo'].values<=0)&(r['hi'].values>=0)).reshape(Rg,Rg)
                s = bool(sup[(sup['plane']==pl)&(sup['pair']==pn)&(sup['outcome']==m)]['supported_crossover'].iloc[0])
                d[(pn,m)] = (tr, s)
        look[pl]=d
    return look

def main():
    b = batch(seeds)
    draws = archetype_draws()
    look = transition_lookup()
    pk = np.array([[k,h] for _,k,h in PLANES]); pnames=[p[0] for p in PLANES]
    rows=[]; grads=[]; overl=[]
    cache_f = f'{RES}/archetypes/draw_outcomes.npz'
    if os.path.exists(cache_f):
        z = np.load(cache_f, allow_pickle=True); OUT = z['OUT'].item()
    else:
        OUT={}
        for name,D in draws.items():
            t0=time.time()
            arr = np.zeros((NDRAW,3,2))
            for di,(p,c,k,h) in enumerate(D):
                r = run_all(b, dict(p=float(p),c=float(c),k=float(k),h=float(h)), TAU)
                for ji,j in enumerate('STV'):
                    arr[di,ji,0]=r[j]['A'].mean(); arr[di,ji,1]=r[j]['C_social'].mean()
            OUT[name]=arr
            print(f'  {name}: {NDRAW} draws evaluated in {time.time()-t0:.0f}s', flush=True)
        np.savez(cache_f, OUT=np.array(OUT,dtype=object))
    A1_flag=None
    for name,D in draws.items():
        arr = OUT[name]
        for mi,m in enumerate(('A','C_social')):
            vals = arr[:,:,mi]
            best = vals.argmax(1) if m=='A' else vals.argmin(1)
            srt  = np.sort(vals,1)
            marg = (srt[:,-1]-srt[:,-2]) if m=='A' else (srt[:,1]-srt[:,0])
            Pl = [(best==j).mean() for j in range(3)]
            modal = int(np.argmax(Pl))
            rows.append(dict(archetype=name, metric=m,
                P_lead_S=Pl[0], P_lead_T=Pl[1], P_lead_V=Pl[2],
                modal_leader='STV'[modal], P_displaced=1-Pl[modal],
                margin_mean=float(marg.mean()), margin_p05=float(np.quantile(marg,.05)),
                margin_p50=float(np.quantile(marg,.50)), margin_p95=float(np.quantile(marg,.95)),
                n_draws=NDRAW))
            if name=='A1' and m=='A' and Pl[2]>0.10:
                A1_flag=f'A1 voluntary flag: P(V leads on A) = {Pl[2]:.2f} under high-enforcement placement.'
        # crossover-band overlap: nearest plane in (k,h), nearest grid cell
        dpl = ((D[:,2:3]-pk[:,0][None])**2 + (D[:,3:4]-pk[:,1][None])**2).argmin(1)
        ip = np.clip(np.round(D[:,0]*(Rg-1)).astype(int),0,Rg-1)
        icx= np.clip(np.round(D[:,1]*(Rg-1)).astype(int),0,Rg-1)
        anyhit = np.zeros(NDRAW,bool)
        for pn,_,_ in pairs:
            for m in ('A','C_social'):
                hit = np.zeros(NDRAW,bool)
                for pli,pl in enumerate(pnames):
                    tr,s = look[pl][(pn,m)]
                    if not s: continue
                    selp = dpl==pli
                    hit[selp] = tr[ip[selp],icx[selp]]
                anyhit |= hit
                overl.append(dict(archetype=name,pair=pn,metric=m,P_overlap=float(hit.mean())))
        overl.append(dict(archetype=name,pair='ANY',metric='ANY',P_overlap=float(anyhit.mean())))
        # local gradients of Delta_jl on (p,c,k,h) over the placement region
        X = np.column_stack([np.ones(NDRAW), D])
        for pi,(pn,i,j) in enumerate(pairs):
            for mi,m in enumerate(('A','C_social')):
                y = arr[:,i,mi]-arr[:,j,mi]
                beta,_res,_rk,_sv = np.linalg.lstsq(X,y,rcond=None)
                grads.append(dict(archetype=name,pair=pn,metric=m,
                    dD_dp=beta[1],dD_dc=beta[2],dD_dk=beta[3],dD_dh=beta[4],
                    intercept=beta[0], mean_delta=float(y.mean())))
    rd = pd.DataFrame(rows); rd.to_csv(f'{RES}/archetypes/archetype_readouts.csv',index=False)
    pd.DataFrame(overl).to_csv(f'{RES}/archetypes/archetype_crossover_overlap.csv',index=False)
    pd.DataFrame(grads).to_csv(f'{RES}/archetypes/archetype_local_gradients.csv',index=False)
    with open(f'{RES}/archetypes/ARCHETYPE_READOUTS.md','w') as f:
        f.write('# Archetype probabilistic readouts (v1.2)\n\n'
                'Every statement below is conditional: "under governance conditions consistent with '
                'this placement evidence" (manuscript Sec V-D). Draws: first %d per archetype from the '
                'frozen seed-2026 joint copula stream; outcomes averaged over the %d principal seeds.\n\n'%(NDRAW,n))
        f.write(rd.to_markdown(index=False)+'\n\n')
        ov = pd.DataFrame(overl); f.write('## Crossover-band overlap (final joint bands, supported pairs only)\n\n')
        f.write(ov.pivot_table(index='archetype',columns=['pair','metric'],values='P_overlap').to_markdown()+'\n\n')
        f.write('## Local gradients of Delta_jl (least squares over placement draws)\n\n')
        f.write(pd.DataFrame(grads).to_markdown(index=False)+'\n\n')
        if A1_flag: f.write('**'+A1_flag+'**\n')
    print(rd.to_string(index=False))
    if A1_flag: print(A1_flag)
    print('archetype readouts written')

if __name__=='__main__':
    main()
