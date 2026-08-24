"""Principal v1.2 stage. Constructions copied VERBATIM from the frozen v1.1
pipeline (outcome_mat, classify, anchor_masks, rng-77 multinomial band weights,
chunked DoD H1), bound to model_v12 and the v1.2 tau_star record (field-equal to
v1.1). Stages:
  A  reference plane (k=0.30,h=0.75), per-plane max-t family exactly as v1.1.
     Gates: (i) S/T raw rows bitwise-identical to the archived v1.1 parquet
     (operationalized reuse); (ii) staged H verdicts equal H_verdicts_v12.json.
  B  four secondary slices (k=0.10,k=0.55 at h=0.75; h=0.50,h=0.95 at k=0.30),
     same per-plane construction (staged records).
  C  paper-wide joint max-t family across ALL 5 planes' cells x 3 pairs per
     metric (single rng-77 weight draw), final per-plane bands/classifications/
     crossovers, final H1 (DoD over pooled anchored cells, memory-chunked),
     H2, H3; staged reference-plane record retained alongside.
Checkpoints per plane under results_v12/principal/ for cheap re-runs.
"""
import numpy as np, json, os, sys, time
import pandas as pd
sys.path.insert(0, os.path.dirname(__file__))
from model_v12 import CONFIG, batch, run_all
from scipy.stats import beta as beta_dist, norm, gaussian_kde

RES = 'results_v12'; CFG = CONFIG
TAU = json.load(open(f'{RES}/verification/tau_star.json'))['tau_star']
os.makedirs(f'{RES}/principal', exist_ok=True)
Rg = json.load(open('results/pilot/grid_resolution_choice.json'))['chosen_resolution']
n  = json.load(open('results/pilot/replication_choice.json'))['chosen_n']
seeds = list(range(100000, 100000 + n))
pp = np.linspace(0, 1, Rg)
pairs = [('S-T', 0, 1), ('S-V', 0, 2), ('T-V', 1, 2)]
B = CFG['B']

PLANES = [('reference', CFG['ref']['k'], CFG['ref']['h']),
          ('k010', 0.10, CFG['ref']['h']), ('k055', 0.55, CFG['ref']['h']),
          ('h050', CFG['ref']['k'], 0.50), ('h095', CFG['ref']['k'], 0.95)]

# ------------- verbatim v1.1 cores (model_v12-bound) -------------
def outcome_mat(bat, cells, want_ext=False):
    S = bat['gamma'].shape[0]
    out = np.zeros((S, len(cells), 3, 2), dtype=np.float64)
    ext = [] if want_ext else None
    for ic,(p,c,k,h) in enumerate(cells):
        r = run_all(bat, dict(p=p,c=c,k=k,h=h), TAU)
        for ji,j in enumerate('STV'):
            out[:,ic,ji,0]=r[j]['A']; out[:,ic,ji,1]=r[j]['C_social']
        if want_ext:
            for ji,j in enumerate('STV'):
                o=r[j]
                ext.append(dict(cell=ic,p=p,c=c,k=k,h=h,instr=j,
                    A=o['A'],C_ab=o['C_ab'],C_social=o['C_social'],nM=o['nM'],nA=o['nA'],
                    nI=o['nI'],T_charges=o['T_ch'],T_fines=o['T_fi'],viol=o['viol'],
                    under=o['under'],part=o['part'],CE=o['CE']))
    return out, ext

def anchor_masks(pp, KR, HR):
    L = np.linalg.cholesky(np.array(CFG['copula_R']))
    masks={}; cent={}
    rng = np.random.default_rng(2026)
    for name,ab in CFG['archetypes'].items():
        z = rng.standard_normal((100000,4))@L.T; u = norm.cdf(z)
        draws = np.column_stack([beta_dist.ppf(u[:,i],*ab[dim])
                                 for i,dim in enumerate(('p','c','k','h'))])
        kde = gaussian_kde(draws.T)
        sub = draws[rng.choice(100000,3000,replace=False)]
        thr = np.quantile(kde(sub.T),0.05)          # 95 percent HDR threshold
        P,C = np.meshgrid(pp,pp,indexing='ij')
        pts = np.column_stack([P.ravel(),C.ravel(),np.full(P.size,KR),np.full(P.size,HR)])
        masks[name] = (kde(pts.T)>=thr).reshape(len(pp),len(pp))
        cent[name] = (ab['p'][0]/sum(ab['p']), ab['c'][0]/sum(ab['c']))
    union = np.zeros_like(next(iter(masks.values())))
    for m in masks.values(): union |= m
    return union, masks, cent

def classify(loA,hiA,loC,hiC):
    """per cell, per ordered pair (first-minus-second)."""
    first_dom  = (loA>0)&(hiC<0); second_dom = (hiA<0)&(loC>0)
    both = ((loA>0)|(hiA<0)) & ((loC>0)|(hiC<0))
    trade = both & ~first_dom & ~second_dom
    a_only = ((loA>0)|(hiA<0)) & ~((loC>0)|(hiC<0))
    c_only = ~((loA>0)|(hiA<0)) & ((loC>0)|(hiC<0))
    unres = ~((loA>0)|(hiA<0)) & ~((loC>0)|(hiC<0))
    return first_dom,second_dom,trade,a_only,c_only,unres

def labels_and_leaders(lo, hi, ncell):
    lab = np.zeros(ncell,dtype=int); Aled = np.zeros(ncell,dtype=int); Cled = np.zeros(ncell,dtype=int)
    anyunres = np.zeros(ncell,dtype=bool)
    dom = {j:np.ones(ncell,dtype=bool) for j in range(3)}
    for pi,(pn,i,j) in enumerate(pairs):
        fd,sd_,tr,ao,co,un = classify(lo[:,pi,0],hi[:,pi,0],lo[:,pi,1],hi[:,pi,1])
        anyunres |= un
        dom[i] &= fd; dom[j] &= sd_
    for j in range(3):
        supA = np.ones(ncell,bool); supC = np.ones(ncell,bool)
        for pi,(pn,i,jj) in enumerate(pairs):
            if i==j: supA &= lo[:,pi,0]>0; supC &= hi[:,pi,1]<0
            elif jj==j: supA &= hi[:,pi,0]<0; supC &= lo[:,pi,1]>0
        Aled[supA]=j+1; Cled[supC]=j+1
    lab[dom[0]]=1; lab[dom[1]]=2; lab[dom[2]]=3
    lab[(lab==0)&anyunres]=0
    lab[(lab==0)&~anyunres]=4
    return lab, Aled, Cled

def crossover_rows(lo, hi, anc, cells):
    xr=[]; H3=False
    for pi,(pn,i,j) in enumerate(pairs):
        for mi,m in enumerate(('A','C_social')):
            posm=(lo[:,pi,mi]>0)&anc; negm=(hi[:,pi,mi]<0)&anc
            trans=(lo[:,pi,mi]<=0)&(hi[:,pi,mi]>=0)&anc
            sup=bool(posm.any() and negm.any())
            H3 |= sup
            tp=[cells[ic][:2] for ic in np.where(trans)[0]]
            xr.append(dict(pair=pn,outcome=m,n_pos=int(posm.sum()),n_neg=int(negm.sum()),
                n_transition=int(trans.sum()),supported_crossover=sup,
                transition_p_range=f'[{min(x[0] for x in tp):.2f},{max(x[0] for x in tp):.2f}]' if tp else '',
                transition_c_range=f'[{min(x[1] for x in tp):.2f},{max(x[1] for x in tp):.2f}]' if tp else ''))
    return xr, H3

def h1_dod(mb_anc, Dhat_anc, ckpt=None):
    """mb_anc [B,nA,3,2], Dhat_anc [nA,3,2]; returns dict pair|metric -> bool (verbatim chunked DoD).
    ckpt: optional path prefix for per-combo checkpoints (resumable)."""
    import time as _t
    H1res={}
    for pi,(pn,i,j) in enumerate(pairs):
        for mi,m in enumerate(('A','C_social')):
            cf = f'{ckpt}_{pn}_{m}.json' if ckpt else None
            if cf and os.path.exists(cf):
                H1res[f'{pn}|{m}'] = json.load(open(cf))['supported']; continue
            _t0=_t.time()
            dlt = mb_anc[:,:,pi,mi]-Dhat_anc[None,:,pi,mi]
            cov = dlt.T@dlt/B; var=np.diag(cov)
            sed = np.sqrt(np.maximum(var[:,None]+var[None,:]-2*cov,1e-30))
            np.fill_diagonal(sed,np.inf)
            ts=[]; nA_=dlt.shape[1]
            ch = max(2, int(2.0e8/(max(nA_,1)**2*8)))
            for s0 in range(0,B,ch):
                dd = dlt[s0:s0+ch,:,None]-dlt[s0:s0+ch,None,:]
                ts.append(np.abs(dd/sed[None]).reshape(dd.shape[0],-1).max(1))
            critD=float(np.quantile(np.concatenate(ts),0.95))
            mdod = np.abs(Dhat_anc[:,pi,mi][:,None]-Dhat_anc[:,pi,mi][None,:])
            sup = bool((mdod/sed > critD).any())
            H1res[f'{pn}|{m}'] = sup
            if cf: json.dump(dict(supported=sup, crit=critD), open(cf,'w'))
            print(f'    DoD {pn}|{m}: {sup} ({_t.time()-_t0:.0f}s)', flush=True)
    return H1res

def raw_frame(ext):
    recs=[]
    for e in ext:
        recs.append(pd.DataFrame({ 'seed':np.array(seeds,dtype=np.int64),
            'p':np.float32(e['p']),'c':np.float32(e['c']),'k':np.float32(e['k']),
            'h':np.float32(e['h']),'instrument':e['instr'],
            **{kk:np.asarray(e[kk],dtype=np.float32) for kk in
               ('A','C_ab','C_social','nM','nA','nI','T_charges','T_fines','viol','under','part','CE')}}))
    raw = pd.concat(recs,ignore_index=True)
    raw['A_norm']=raw['A']/np.float32(CFG['A_cum_nominal'])
    return raw

def plane_cells(k, h):
    return [(pp[ip],pp[icx],k,h) for ip in range(Rg) for icx in range(Rg)]

def compute_plane(name, k, h, b):
    f_out = f'{RES}/principal/plane_{name}_out.npy'
    f_raw = f'{RES}/principal/plane_{name}_raw.parquet'
    cells = plane_cells(k,h)
    if os.path.exists(f_out) and os.path.exists(f_raw):
        return np.load(f_out), cells
    t0=time.time()
    out, ext = outcome_mat(b, cells, want_ext=True)
    raw_frame(ext).to_parquet(f_raw, index=False)
    np.save(f_out, out)
    print(f'  plane {name}: computed in {time.time()-t0:.0f}s')
    return out, cells

def band_stage(out, W):
    D = np.stack([out[:,:,i,:]-out[:,:,j,:] for _,i,j in pairs],2)   # [S,cells,3,2]
    Dhat = D.mean(0)
    mb = (W@D.reshape(n,-1)).reshape(B,D.shape[1],3,2)
    se = mb.std(0,ddof=1)
    return D, Dhat, mb, se

def crit_from(mb, Dhat, se):
    crit=np.zeros(2)
    for mi in range(2):
        s = np.where(se[:,:,mi]>0, se[:,:,mi], np.inf)
        t = np.abs((mb[:,:,:,mi]-Dhat[None,:,:,mi])/s[None]).reshape(B,-1).max(1)
        crit[mi]=float(np.quantile(t,0.95))
    return crit

def bands_frame(cells, Dhat, se, lo, hi):
    return pd.DataFrame([dict(cell=ic,p=cells[ic][0],c=cells[ic][1],pair=pn,metric=m,
                delta=Dhat[ic,pi,mi],se=se[ic,pi,mi],lo=lo[ic,pi,mi],hi=hi[ic,pi,mi])
                for ic in range(len(cells)) for pi,(pn,_,_) in enumerate(pairs)
                for mi,m in enumerate(('A','C_social'))])

def summary_frame(cells, means, anc, lab, Aled, Cled):
    summ=[]
    for ic,(p,c,k,h) in enumerate(cells):
        for ji,j in enumerate('STV'):
            summ.append(dict(p=p,c=c,instrument=j,mean_A=means[ic,ji,0],
                mean_C_social=means[ic,ji,1],anchored=bool(anc[ic]),
                cell_label=['unresolved','S_dominant','T_dominant','V_dominant','no_dominance_mixed'][lab[ic]],
                A_leader=['none','S','T','V'][Aled[ic]],C_leader=['none','S','T','V'][Cled[ic]]))
    return pd.DataFrame(summ)

# ---------------------------------------------------------------
def main():
    b = batch(seeds)
    plane_data = {}
    for name,k,h in PLANES:
        out, cells = compute_plane(name,k,h,b)
        plane_data[name] = dict(out=out, cells=cells, k=k, h=h)

    # ---- Gate: S/T bitwise reuse + V liveness on the reference plane ----
    new = pd.read_parquet(f'{RES}/principal/plane_reference_raw.parquet')
    old = pd.read_parquet('results/principal/principal_raw_results.parquet')
    mcols = ['A','C_ab','C_social','nM','nA','nI','T_charges','T_fines','viol','under','part','CE','A_norm']
    stn, sto = new[new.instrument!='V'].reset_index(drop=True), old[old.instrument!='V'].reset_index(drop=True)
    assert len(stn)==len(sto) and (stn['instrument'].values==sto['instrument'].values).all()
    assert all((stn[c].values==sto[c].values).all() for c in ('seed','p','c','k','h'))
    def _beq(a, b):  # bitwise incl. NaN-position identity (CE holds NaNs by construction)
        return bool((( a==b ) | (np.isnan(a) & np.isnan(b))).all())
    st_bit = all(_beq(stn[c].values.astype(np.float64), sto[c].values.astype(np.float64)) for c in mcols)
    v_old_inert = bool((old[old.instrument=='V']['A'].values==0).all())
    v_new_live  = bool((new[new.instrument=='V']['A'].values>0).any())
    print(f'GATE S/T bitwise reuse: {"PASS" if st_bit else "FAIL"}; archived V inert={v_old_inert}, v1.2 V live={v_new_live}')
    assert st_bit and v_old_inert and v_new_live, 'STOP: raw reuse gate failed'
    # canonical filename for the reference plane (v1.1-compatible)
    new.to_parquet(f'{RES}/principal/principal_raw_results.parquet', index=False)

    # ---- anchors per plane ----
    anchors = {}
    for name,k,h in PLANES:
        f_anc = f'{RES}/principal/plane_{name}_anchor.npz'
        if os.path.exists(f_anc):
            z = np.load(f_anc, allow_pickle=True)
            anchors[name] = dict(union=z['union'], masks=z['masks'].item(), cent=z['cent'].item())
        else:
            union, masks, cent = anchor_masks(pp, k, h)
            np.savez(f_anc, union=union, masks=np.array(masks,dtype=object), cent=np.array(cent,dtype=object))
            anchors[name] = dict(union=union, masks=masks, cent=cent)
        print(f'  anchors {name}: anchored fraction {anchors[name]["union"].mean():.3f}')

    # ================= Stage A: reference plane, per-plane family =================
    rngA = np.random.default_rng(77)
    W_A = rngA.multinomial(n,[1/n]*n,size=B)/n
    ref = plane_data['reference']
    D,Dhat,mb,se = band_stage(ref['out'], W_A)
    # consistency: S-T delta/se must equal archived v1.1 bands (S,T unchanged, same W)
    ob = pd.read_parquet('results/principal/simultaneous_bands.parquet')
    ost = ob[(ob['pair']=='S-T')].sort_values(['cell','metric']).reset_index(drop=True)
    mine = bands_frame(ref['cells'], Dhat, se, Dhat, Dhat)  # lo/hi placeholders
    mst = mine[(mine['pair']=='S-T')].sort_values(['cell','metric']).reset_index(drop=True)
    st_delta = np.allclose(ost['delta'].values, mst['delta'].values, rtol=0, atol=0)
    st_se    = np.allclose(ost['se'].values,    mst['se'].values,    rtol=0, atol=0)
    print(f'S-T delta/se vs archived v1.1 bands: delta {"EXACT" if st_delta else "DIFF"}, se {"EXACT" if st_se else "DIFF"}')
    assert st_delta and st_se, 'STOP: S-T band inputs deviate from archived construction'
    critA = crit_from(mb, Dhat, se)
    lo = Dhat - critA[None,None,:]*se; hi = Dhat + critA[None,None,:]*se
    json.dump(dict(B=B,alpha=0.05,family='per-plane (staged): reference-plane cells x 3 pairwise surfaces per metric',
                   crit_A=critA[0],crit_C=critA[1],n=n,grid=Rg,seeds='100000..%d'%(100000+n-1)),
              open(f'{RES}/principal/max_t_summary_staged_reference.json','w'),indent=1)
    bands_frame(ref['cells'],Dhat,se,lo,hi).to_parquet(f'{RES}/principal/simultaneous_bands.parquet',index=False)
    anc = anchors['reference']['union'].reshape(-1)
    lab,Aled,Cled = labels_and_leaders(lo,hi,len(ref['cells']))
    means = ref['out'].mean(0)
    summary_frame(ref['cells'],means,anc,lab,Aled,Cled).to_csv(f'{RES}/principal/principal_summary.csv',index=False)
    xr,H3 = crossover_rows(lo,hi,anc,ref['cells'])
    pd.DataFrame(xr).to_csv(f'{RES}/principal/crossover_summary.csv',index=False)
    aidx = np.where(anc)[0]
    H1res = h1_dod(mb[:,aidx], Dhat[aidx])
    H2res = {pn: bool(((lo[aidx,pi,0]>0)|(hi[aidx,pi,0]<0)).any() or
                      ((lo[aidx,pi,1]>0)|(hi[aidx,pi,1]<0)).any()) for pi,(pn,_,_) in enumerate(pairs)}
    H1 = any(H1res.values()); H2 = any(H2res.values())
    staged = dict(H1_supported=H1,H1_detail=H1res,H2_supported=H2,H2_detail=H2res,
                  H3_supported=H3,H3_REJECTED=(not H3),
                  scope='reference plane k=0.30, h=0.75; anchor cells only; slices pending',
                  v5prime='external voluntary interpretation retained (V5-prime PASS)')
    json.dump(staged, open(f'{RES}/principal/H_verdicts_staged_reference_plane.json','w'),indent=1)
    arch = json.load(open('reference_v12_artifacts/H_verdicts_v12.json'))
    gate = (staged['H1_detail']==arch['H1_detail'] and staged['H2_detail']==arch['H2_detail']
            and staged['H1_supported']==arch['H1_supported'] and staged['H2_supported']==arch['H2_supported']
            and staged['H3_supported']==arch['H3_supported'])
    print(f'GATE staged reference-plane H verdicts vs archived H_verdicts_v12.json: {"MATCH" if gate else "MISMATCH"}')
    print('  staged:', json.dumps({k:staged[k] for k in ('H1_supported','H2_supported','H3_supported')}),
          'crit_A=%.4f crit_C=%.4f'%(critA[0],critA[1]),
          'anchored fraction %.2f'%anc.mean(), 'unresolved %.2f'%(lab==0).mean())
    assert gate, 'STOP: staged verdicts do not reproduce archived v1.2 reference-plane record'
    np.save(f'{RES}/principal/staged_reference_lo.npy', lo); np.save(f'{RES}/principal/staged_reference_hi.npy', hi)

    # ================= Stage C: paper-wide joint family =================
    rngC = np.random.default_rng(77)
    W = rngC.multinomial(n,[1/n]*n,size=B)/n
    tmax = {0: np.zeros(B), 1: np.zeros(B)}
    stash = {}
    for name,k,h in PLANES:
        D,Dh,mb_,se_ = band_stage(plane_data[name]['out'], W)
        for mi in range(2):
            s = np.where(se_[:,:,mi]>0, se_[:,:,mi], np.inf)
            t = np.abs((mb_[:,:,:,mi]-Dh[None,:,:,mi])/s[None]).reshape(B,-1).max(1)
            tmax[mi] = np.maximum(tmax[mi], t)
        aidx_ = np.where(anchors[name]['union'].reshape(-1))[0]
        stash[name] = dict(Dhat=Dh, se=se_, mb_anc=mb_[:,aidx_], Dhat_anc=Dh[aidx_], aidx=aidx_)
        del mb_
    critF = np.array([float(np.quantile(tmax[0],0.95)), float(np.quantile(tmax[1],0.95))])
    json.dump(dict(B=B,alpha=0.05,
                   family='paper-wide joint: ALL 5 planes cells x 3 pairwise surfaces per metric',
                   planes=[p[0] for p in PLANES],crit_A=critF[0],crit_C=critF[1],n=n,grid=Rg,
                   seeds='100000..%d'%(100000+n-1)),
              open(f'{RES}/principal/max_t_summary_final.json','w'),indent=1)
    print(f'paper-wide joint crit_A={critF[0]:.4f} crit_C={critF[1]:.4f} (staged ref: {critA[0]:.4f}/{critA[1]:.4f})')

    all_xr=[]; H3f=False; H2f={pn:False for pn,_,_ in pairs}
    mb_anc_all=[]; Dhat_anc_all=[]
    for name,k,h in PLANES:
        st = stash[name]; cells = plane_data[name]['cells']
        loF = st['Dhat'] - critF[None,None,:]*st['se']; hiF = st['Dhat'] + critF[None,None,:]*st['se']
        bands_frame(cells,st['Dhat'],st['se'],loF,hiF).to_parquet(
            f'{RES}/principal/plane_{name}_bands_final_joint.parquet',index=False)
        ancv = anchors[name]['union'].reshape(-1)
        labF,AledF,CledF = labels_and_leaders(loF,hiF,len(cells))
        summary_frame(cells,plane_data[name]['out'].mean(0),ancv,labF,AledF,CledF)\
            .to_csv(f'{RES}/principal/plane_{name}_summary_final.csv',index=False)
        xrF,h3p = crossover_rows(loF,hiF,ancv,cells)
        for r in xrF: r['plane']=name
        all_xr += xrF; H3f |= h3p
        aidx_ = st['aidx']
        for pi,(pn,_,_) in enumerate(pairs):
            H2f[pn] |= bool(((loF[aidx_,pi,0]>0)|(hiF[aidx_,pi,0]<0)).any() or
                            ((loF[aidx_,pi,1]>0)|(hiF[aidx_,pi,1]<0)).any())
        mb_anc_all.append(st['mb_anc']); Dhat_anc_all.append(st['Dhat_anc'])
        print(f'  final {name}: unresolved {(labF==0).mean():.2f}, anchored {ancv.mean():.2f}, '
              f'labels S/T/V/mixed = {[(labF==z).sum() for z in (1,2,3,4)]}')
    pd.DataFrame(all_xr).to_csv(f'{RES}/principal/crossover_summary_final_all_planes.csv',index=False)
    mb_anc = np.concatenate(mb_anc_all,axis=1); Dh_anc = np.concatenate(Dhat_anc_all,axis=0)
    print(f'  pooled anchored cells across planes: {Dh_anc.shape[0]}; running final DoD H1 (chunked)...')
    t0=time.time()
    H1resF = h1_dod(mb_anc, Dh_anc, ckpt=f'{RES}/principal/h1_final')
    print(f'  DoD H1 done in {time.time()-t0:.0f}s')
    H1F = any(H1resF.values()); H2F = any(H2f.values())
    json.dump(dict(H1_supported=H1F,H1_detail=H1resF,H2_supported=H2F,H2_detail=H2f,
                   H3_supported=H3f,H3_REJECTED=(not H3f),
                   H3_note='H3 evaluated within-plane (both supported signs inside the same plane anchored region); overall = any plane',
                   scope='final: all 5 planes (reference + k010/k055/h050/h095), paper-wide joint family, anchor cells only',
                   staged_record='H_verdicts_staged_reference_plane.json retained',
                   v5prime='external voluntary interpretation retained (V5-prime PASS)'),
              open(f'{RES}/principal/H_verdicts_final.json','w'),indent=1)
    print('FINAL H1',H1F,json.dumps(H1resF),'H2',H2F,json.dumps(H2f),'H3',H3f)

if __name__ == '__main__':
    main()
