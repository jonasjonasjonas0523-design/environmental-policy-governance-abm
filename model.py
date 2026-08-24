"""Core ABM per CODE_INPUT_SPECIFICATION.md (frozen 2026-08-20), Rev2 manuscript,
PRE_REGISTRATION_SPECIFICATION.md. Model version v1.1-prerun-amended.

PRE-RUN AMENDMENTS (author-authorized "correct the objectively wrong ones"; see
AMENDMENTS_PRE_RUN.md):
 A1. Voluntary participation comparison carries a theta-differential: a participant
     whose achievement is VERIFIED (prob p, same monitoring environment) is exempt
     from the threatened future regime. Decision rule: participate-and-achieve iff
        C_i(target) <= p * ( rho_i + delta*h*theta_t*M_i(g) ).
     Rationale: literal eq.(13) puts the identical delta*h*theta*M_i term in every
     option, so theta cancels from the decision, contradicting the abstract, II-C,
     ref [14], validation V5(iii), and ablation m4. Baseline numerics (nu, theta_0,
     lambda) are NOT changed.
 A2. Verification item 12 anchor restated: tau* must agree with an independent
     quadrature evaluation (tol 1e-6 on F/N); 0.2/e is an order-of-magnitude floor
     only (measured cap correction is reported, not assumed "slight").

PIN-DOWNS implemented as documented in MISSING_OR_CONFLICTING_INPUTS.md v2:
 P1 A* roles: s_t uses per-seed A* = 0.2*sum(E0); thresholds/epsilon_A/A_min use
    nominal A*_cum = 400.  P2 theta timing: firms act on published theta from
    s_{t-1}; theta_1 = logistic(theta_0).  P3 violation rate = share firm-periods
    with E > Ebar.  P4 under-reporting share = sum(u)/sum(E).
 P5 kappa_I counts implemented CONSEQUENCES (standard penalties; charge evasion
    recoveries); routine charge collection is administration (kappa_A assessment).

CRN: per replicate seed, one population draw and shared uniform stream arrays
U_mon[T,N], U_imp[T,N] (primary implementation), U_imp2[T,N] (charge-only
secondary implementation: evasion recovery), consumed identically across
instruments and reused across governance cells. Leakage is a deterministic (1-k)
scaling of implemented consequences (PRE_REG Sec 6 item 7).
"""
import numpy as np, yaml, os

CONFIG = dict(
    N=200, T=10, delta=0.95, A_star_frac=0.20, A_cum_nominal=400.0, A_min=20.0,
    m_S=2.0, m_T=1.0, nu=0.25, theta0=-2.2, lam=3.0,
    kappa_M=0.002, kappa_A=0.0005, kappa_I=0.001,
    ref=dict(p=0.50, c=0.50, k=0.30, h=0.75),
    k_slices=[0.10, 0.55], h_slices=[0.50, 0.95],
    mu_g=-0.5, sg=1.0, mu_E=-0.5, sE=1.0,
    alpha=0.05, eps_A=20.0, B=2000,
    seed_pilot_lo=1, seed_pilot_hi=10**4, seed_principal_lo=10**5,
    archetypes=dict(
        A1=dict(p=(14.4,3.6), c=(17,3),  k=(2,38),      h=(57,3)),
        A2=dict(p=(11,9),     c=(12,8),  k=(6.25,18.75),h=(25.5,4.5)),
        A3=dict(p=(8,12),     c=(10,10), k=(7.5,17.5),  h=(18.6,11.4)),
        A4=dict(p=(3.96,18.04),c=(4.5,10.5),k=(8.36,13.64),h=(20,5)),
        A5=dict(p=(7.7,14.3), c=(9,11),  k=(6.6,15.4),  h=(46.5,3.5))),
    copula_R=[[1,.5,-.4,0],[.5,1,-.5,0],[-.4,-.5,1,0],[0,0,0,1]],
    val=dict(V1=dict(c=.85,k=.05,h=.95,p0=.50,p1=.90),
             V3=dict(eng=(.40,.50,.30,.62), dis=(.15,.25,.40,.80)),
             V4=dict(p=.95,c=.90,k=.03,h=.95),
             V5=dict(A1=(.80,.85,.05,.95), A4=(.18,.30,.38,.80))),
)

def export_config(path):
    with open(path,'w') as f: yaml.safe_dump(CONFIG, f, sort_keys=False)

def logistic(x): return 1.0/(1.0+np.exp(-x))

# ---------------- population and CRN streams ----------------
def seed_data(seed, cfg=CONFIG):
    rng = np.random.default_rng(int(seed))
    N,T = cfg['N'], cfg['T']
    gamma = np.exp(rng.normal(cfg['mu_g'], cfg['sg'], N))
    E0    = np.exp(rng.normal(cfg['mu_E'], cfg['sE'], N))
    return dict(gamma=gamma, E0=E0,
                U_mon=rng.random((T,N)), U_imp=rng.random((T,N)),
                U_imp2=rng.random((T,N)))

def batch(seeds, cfg=CONFIG):
    ds = [seed_data(s,cfg) for s in seeds]
    return {k: np.stack([d[k] for d in ds]) for k in ds[0]}   # [S,...]

# ---------------- tau* ----------------
def _quad(nq, qseed, cfg=CONFIG):
    r = np.random.default_rng(qseed)
    return (np.exp(r.normal(cfg['mu_g'],cfg['sg'],nq)),
            np.exp(r.normal(cfg['mu_E'],cfg['sE'],nq)))

def F_of_tau(tau, gq, Eq, cfg=CONFIG):
    return float(np.mean(np.minimum(tau/gq, Eq)) - cfg['A_star_frac'])

def solve_tau_star(nq=10**6, qseed=424242, tol=1e-6, cfg=CONFIG):
    gq,Eq = _quad(nq,qseed,cfg); lo,hi = 1e-5, 1.0
    while F_of_tau(hi,gq,Eq,cfg) < 0: hi *= 2
    for _ in range(200):
        mid = 0.5*(lo+hi)
        if abs(F_of_tau(mid,gq,Eq,cfg)) < tol: return mid
        if F_of_tau(mid,gq,Eq,cfg) < 0: lo = mid
        else: hi = mid
    return 0.5*(lo+hi)

# ---------------- derived M_i ----------------
def M_of(p,c,k,f_S,gamma,E0):
    q = p*c*(1.0-k)*f_S
    aS = np.minimum(q/gamma, 0.2*E0)
    return 0.5*gamma*aS**2 + q*np.maximum(0.2*E0-aS,0.0), aS, q

# ---------------- instruments (batched over seeds S, firms N) ----------------
def _base_out(S): 
    z = lambda: np.zeros(S)
    return dict(A=z(),C_ab=z(),nM=z(),nA=z(),nI=z(),T_ch=z(),T_fi=z(),
                viol=z(),under=z(),part=z(),C_social=z())

def run_standard(bat, g, tau_s, cfg=CONFIG):
    p,c,k = g['p'],g['c'],g['k']; f_S = cfg['m_S']*tau_s
    gam,E0 = bat['gamma'],bat['E0']; S,N = gam.shape; T = cfg['T']
    M,a,q = M_of(p,c,k,f_S,gam,E0)     # standard optimum a == a_S
    v = np.maximum(0.2*E0 - a, 0.0); E = E0 - a
    insp = bat['U_mon'] < p                      # [S,T,N]
    violm = (v > 1e-15)[:,None,:]
    impl = insp & violm & (bat['U_imp'] < c)
    o = _base_out(S)
    o['A']    = T*a.sum(1)
    o['C_ab'] = T*(0.5*gam*a**2).sum(1)
    o['nM']   = insp.sum((1,2)).astype(float)
    o['nA']   = float(N*T)*np.ones(S)
    o['nI']   = impl.sum((1,2)).astype(float)
    o['T_fi'] = ((1-k)*f_S*v[:,None,:]*impl).sum((1,2))
    o['viol'] = (v>1e-15).mean(1)                # P3: constant over ticks under the standard
    o['under']= np.zeros(S)
    _finish(o,cfg); return o

def run_charge(bat, g, tau_s, cfg=CONFIG):
    p,c,k = g['p'],g['c'],g['k']; tau=tau_s; f_T=cfg['m_T']*tau_s
    gam,E0 = bat['gamma'],bat['E0']; S,N = gam.shape; T = cfg['T']
    thr = tau/(tau+f_T) if (tau+f_T) > 0 else 0.0
    truthful = (p >= thr - 1e-12)                # tie -> truthful; = p>=0.5 at baseline
    price = c*(1-k)*tau if truthful else p*c*(1-k)*(tau+f_T)
    a = np.minimum(price/gam, E0); E = E0-a
    Etil = E if truthful else np.zeros_like(E); u = E-Etil
    insp = bat['U_mon'] < p
    coll = bat['U_imp'] < c
    det  = insp & (u > 1e-15)[:,None,:]
    rec  = det & (bat['U_imp2'] < c)
    o = _base_out(S)
    o['A']=T*a.sum(1); o['C_ab']=T*(0.5*gam*a**2).sum(1)
    o['nM']=insp.sum((1,2)).astype(float); o['nA']=float(N*T)*np.ones(S)
    o['nI']=rec.sum((1,2)).astype(float)
    o['T_ch']=((1-k)*tau*Etil[:,None,:]*coll).sum((1,2)) + ((1-k)*tau*u[:,None,:]*rec).sum((1,2))
    o['T_fi']=((1-k)*f_T*u[:,None,:]*rec).sum((1,2))
    denom = np.maximum(T*E.sum(1),1e-12)
    o['under']=T*u.sum(1)/denom
    o['viol']=np.zeros(S)
    _finish(o,cfg); return o

def run_voluntary(bat, g, tau_s, cfg=CONFIG):
    p,c,k,h = g['p'],g['c'],g['k'],g['h']; f_S=cfg['m_S']*tau_s
    gam,E0 = bat['gamma'],bat['E0']; S,N = gam.shape; T=cfg['T']; dl=cfg['delta']
    M,_,_ = M_of(p,c,k,f_S,gam,E0)
    Ct  = 0.5*gam*(0.2*E0)**2                    # cost of hitting the target
    rho = cfg['nu']*Ct
    Astar_seed = cfg['A_star_frac']*E0.sum(1)    # per-seed A* per period (P1)
    theta = np.full(S, logistic(cfg['theta0']))  # theta_1 (P2)
    o = _base_out(S); theta_path=[]; s_path=[]; Rp=np.zeros(S)
    for t in range(T):
        threat = dl*h*theta[:,None]*M
        part = Ct <= p*(rho + threat)            # AMENDMENT A1 (verification-gated exemption)
        a = 0.2*E0*part
        z = part & (bat['U_mon'][:,t,:] < p)     # verified achievement
        Ahat = (0.2*E0*z).sum(1)
        s = np.maximum((Astar_seed-Ahat)/Astar_seed, 0.0)
        theta_path.append(theta.copy()); s_path.append(s.copy())
        theta = logistic(cfg['theta0']+cfg['lam']*s)
        Rp += (rho*z).sum(1)                     # paid only on d=1 and verified z=1
        o['A']    += a.sum(1)
        o['C_ab'] += (0.5*gam*a**2).sum(1)
        o['nM']   += z.sum(1)
        o['nA']   += part.sum(1)
        o['part'] += part.mean(1)/T
    o['theta_path']=np.array(theta_path)         # [T,S]
    o['s_path']=np.array(s_path); o['R_paid']=Rp
    _finish(o,cfg); return o

def _finish(o,cfg):
    o['C_social'] = o['C_ab'] + cfg['kappa_M']*o['nM'] + cfg['kappa_A']*o['nA'] + cfg['kappa_I']*o['nI']
    o['T_fiscal'] = o['T_ch'] + o['T_fi']
    A = o['A']
    o['CE'] = np.where(A >= cfg['A_min'], o['C_social']/np.maximum(A,1e-12), np.nan)

RUN = dict(S=run_standard, T=run_charge, V=run_voluntary)

def run_all(bat, g, tau_s, cfg=CONFIG):
    return {j: RUN[j](bat,g,tau_s,cfg) for j in ('S','T','V')}

if __name__ == '__main__':
    os.makedirs('config',exist_ok=True); export_config('config/baseline.yaml')
    ts = solve_tau_star(10**5)
    b = batch([1,2,3]); out = run_all(b, dict(p=.5,c=.5,k=.3,h=.75), ts)
    print('smoke ok; tau*(1e5 quad) =', round(ts,5),
          {j: round(float(out[j]['A'].mean()),2) for j in out})
