"""Global sensitivity (CODE_INPUT Sec 15; manuscript Sec V-C): variance-based
Sobol indices via Saltelli design, base N=1024 (power of 2, qmc.Sobol scrambled,
seeds logged), 10 CRN replicate seeds per point (150000..150009, disjoint from
pilot/validation/principal; outcomes = mean over replicates), outputs = A and
C_social per instrument plus each pairwise difference Delta_jl (12 channels).
Groups reported separately:
  governance (d=4): p [0.05,0.95], c [0.10,0.95], k [0.02,0.55], h [0.50,0.99]
     (Table 2 spans), non-governance parameters at baseline.
  non-governance (d=9): delta [0.90,0.97]; sigma_g, sigma_E [0.5,1.5] with
     mean-1 recentering mu = -sigma^2/2 (CRN preserved: Generator.normal keeps
     the underlying stream under (mu, sigma) changes); m_S [1,4]; m_T [0.5,2];
     nu-coordinate = nu_0 [0.05,0.75] with eta fixed at its v1.2 value (frozen
     nu-coordinate reading); lambda [1.5,4.5]; theta_0 [-2.94,-1.73];
     kappa-scale log-uniform [0.2,5] multiplying all three kappas.
     Governance fixed at the reference point; tau* held at its frozen value
     (instrument design constant, never re-solved).
Estimators: first-order Saltelli-2010/Jansen S_i = mean(f_B*(f_ABi-f_A))/V;
total Jansen T_i = mean((f_A-f_ABi)^2)/(2V); V = var([f_A;f_B]).
"""
import numpy as np, json, os, sys, time
import pandas as pd
sys.path.insert(0, os.path.dirname(__file__))
from model_v12 import CONFIG, batch, run_all
from scipy.stats import qmc

RES='results_v12'; CFG=CONFIG
TAU = json.load(open(f'{RES}/verification/tau_star.json'))['tau_star']
os.makedirs(f'{RES}/sobol', exist_ok=True)
NBASE = 1024
CRN_SEEDS = list(range(150000, 150010))
QSEED = dict(gov=31001, nongov=31002)
CHAN = ['A_S','A_T','A_V','C_S','C_T','C_V','D_ST_A','D_SV_A','D_TV_A','D_ST_C','D_SV_C','D_TV_C']

GOV = [('p',0.05,0.95,'lin'),('c',0.10,0.95,'lin'),('k',0.02,0.55,'lin'),('h',0.50,0.99,'lin')]
NONGOV = [('delta',0.90,0.97,'lin'),('sigma_g',0.5,1.5,'lin'),('sigma_E',0.5,1.5,'lin'),
          ('m_S',1.0,4.0,'lin'),('m_T',0.5,2.0,'lin'),('nu0',0.05,0.75,'lin'),
          ('lam',1.5,4.5,'lin'),('theta0',-2.94,-1.73,'lin'),('kappa_scale',0.2,5.0,'log')]

def scale(u, spec):
    x = np.empty_like(u)
    for i,(_,lo,hi,tp) in enumerate(spec):
        if tp=='log': x[:,i] = np.exp(np.log(lo)+u[:,i]*(np.log(hi)-np.log(lo)))
        else:         x[:,i] = lo+u[:,i]*(hi-lo)
    return x

def eval_gov(X, b):
    out = np.zeros((len(X), len(CHAN)))
    for r,(p,c,k,h) in enumerate(X):
        rr = run_all(b, dict(p=float(p),c=float(c),k=float(k),h=float(h)), TAU)
        v = [rr['S']['A'].mean(),rr['T']['A'].mean(),rr['V']['A'].mean(),
             rr['S']['C_social'].mean(),rr['T']['C_social'].mean(),rr['V']['C_social'].mean()]
        out[r,:6]=v
        out[r,6:9]  = v[0]-v[1], v[0]-v[2], v[1]-v[2]
        out[r,9:12] = v[3]-v[4], v[3]-v[5], v[4]-v[5]
    return out

def eval_nongov(X):
    g = dict(CFG['ref'])
    out = np.zeros((len(X), len(CHAN)))
    for r,row in enumerate(X):
        d = dict(zip([s[0] for s in NONGOV], row))
        cfg = dict(CFG)
        cfg['delta']=d['delta']; cfg['sg']=d['sigma_g']; cfg['sE']=d['sigma_E']
        cfg['mu_g']=-d['sigma_g']**2/2; cfg['mu_E']=-d['sigma_E']**2/2
        cfg['m_S']=d['m_S']; cfg['m_T']=d['m_T']
        cfg['nu0']=d['nu0']                      # eta stays at frozen v1.2 value
        cfg['lam']=d['lam']; cfg['theta0']=d['theta0']
        for kk in ('kappa_M','kappa_A','kappa_I'): cfg[kk]=CFG[kk]*d['kappa_scale']
        b = batch(CRN_SEEDS, cfg)                # CRN: same seeds, stream-preserving redraw
        rr = run_all(b, g, TAU, cfg)
        v = [rr['S']['A'].mean(),rr['T']['A'].mean(),rr['V']['A'].mean(),
             rr['S']['C_social'].mean(),rr['T']['C_social'].mean(),rr['V']['C_social'].mean()]
        out[r,:6]=v
        out[r,6:9]  = v[0]-v[1], v[0]-v[2], v[1]-v[2]
        out[r,9:12] = v[3]-v[4], v[3]-v[5], v[4]-v[5]
    return out

def saltelli_indices(fA, fB, fAB):
    """fA,fB [N,ch]; fAB [d,N,ch] -> (S1, ST) each [d,ch]"""
    V = np.var(np.concatenate([fA,fB],0), axis=0, ddof=1)
    V = np.where(V>0, V, np.inf)
    S1 = np.stack([(fB*(fAB[i]-fA)).mean(0)/V for i in range(fAB.shape[0])])
    ST = np.stack([((fA-fAB[i])**2).mean(0)/(2*V) for i in range(fAB.shape[0])])
    return S1, ST

def run_group(tag, spec, evalf, needs_batch):
    d = len(spec)
    f_np = f'{RES}/sobol/{tag}_f.npz'
    if os.path.exists(f_np):
        z=np.load(f_np); fA,fB,fAB = z['fA'],z['fB'],z['fAB']
    else:
        eng = qmc.Sobol(2*d, scramble=True, seed=QSEED[tag])
        U = eng.random_base2(int(np.log2(NBASE)))
        A = scale(U[:,:d], spec); Bm = scale(U[:,d:], spec)
        t0=time.time()
        if needs_batch:
            b = batch(CRN_SEEDS)
            fA = evalf(A,b); fB = evalf(Bm,b)
            fAB = np.stack([evalf(np.column_stack([A[:,:i],Bm[:,i:i+1],A[:,i+1:]]),b) for i in range(d)])
        else:
            fA = evalf(A); fB = evalf(Bm)
            fAB = np.stack([evalf(np.column_stack([A[:,:i],Bm[:,i:i+1],A[:,i+1:]])) for i in range(d)])
        np.savez_compressed(f_np, fA=fA, fB=fB, fAB=fAB)
        print(f'  {tag}: {NBASE*(d+2)} evaluations in {time.time()-t0:.0f}s', flush=True)
    S1, ST = saltelli_indices(fA,fB,fAB)
    rows=[]
    for i,(nm,*_) in enumerate(spec):
        for ci,ch in enumerate(CHAN):
            rows.append(dict(group=tag,param=nm,channel=ch,S1=S1[i,ci],ST=ST[i,ci]))
    return pd.DataFrame(rows)

def main():
    dfg = run_group('gov', GOV, eval_gov, needs_batch=True)
    dfn = run_group('nongov', NONGOV, eval_nongov, needs_batch=False)
    df = pd.concat([dfg,dfn],ignore_index=True)
    df.to_csv(f'{RES}/sobol/sobol_indices.csv',index=False)
    json.dump(dict(base=NBASE, crn_seeds=f'{CRN_SEEDS[0]}..{CRN_SEEDS[-1]}', qmc_seeds=QSEED,
                   estimators='S1 Saltelli-2010, ST Jansen', channels=CHAN,
                   gov_ranges={s[0]:[s[1],s[2]] for s in GOV},
                   nongov_ranges={s[0]:[s[1],s[2],s[3]] for s in NONGOV},
                   notes=['non-gov evaluated at reference governance; tau* frozen',
                          'sigma variants recentered mu=-sigma^2/2 (mean-1), CRN stream-preserving',
                          'nu-coordinate = nu_0 with eta fixed (v1.2 reading)']),
              open(f'{RES}/sobol/sobol_design.json','w'),indent=1)
    for tag in ('gov','nongov'):
        piv_T = df[df.group==tag].pivot(index='param',columns='channel',values='ST')[CHAN]
        print(f'--- {tag}: total-order (Jansen) ---'); print(piv_T.round(3).to_string())
    with open(f'{RES}/sobol/SOBOL_REPORT.md','w') as f:
        f.write('# Global sensitivity (Sobol/Saltelli, v1.2)\n\nDesign per CODE_INPUT Sec 15; see sobol_design.json.\n\n')
        for tag in ('gov','nongov'):
            f.write(f'## {tag} group: total-order (Jansen)\n\n')
            f.write(df[df.group==tag].pivot(index='param',columns='channel',values='ST')[CHAN].round(4).to_markdown()+'\n\n')
            f.write(f'## {tag} group: first-order (Saltelli-2010)\n\n')
            f.write(df[df.group==tag].pivot(index='param',columns='channel',values='S1')[CHAN].round(4).to_markdown()+'\n\n')
    print('sobol written')

if __name__=='__main__':
    main()
