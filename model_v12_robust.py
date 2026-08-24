"""Robustness-variant instrument implementations for v1.2 (families F1/F2/F4),
PRE_REG Sec 8 + AMENDMENTS_V1_2 Sec 3.2 + PARAMETER_PROVENANCE ranges.
The freeze names each variant and its ranges; operational functional forms were
left open and are instantiated HERE, ex ante, as the documented record (mirrors
the earlier session's settled list; see FINAL_COMPUTATIONAL_AUDIT):

F1 behavioral
  cara r in {1,10}: per-instrument certainty-equivalent decisions.
    S: choose a in [0,req] minimizing per-period CE(a) = C(a) +
       (1/r) ln(1-pc+pc e^{r(1-k)f_S v}), v=req-a (40-step bisection on the CE
       first-order condition; CE convex in a for the quadratic form).
    T: branch CE with collection Bernoulli(c) (truthful) or recovery
       Bernoulli(pc) (evading); a per branch by the same bisection on the
       branch's exposure; branch with lower CE chosen.
    V: certain cost K_i vs risky gain B = rho_i + threat received w.p. p:
       CE_gain = -(1/r) ln(1-p+p e^{-rB}); participate iff CE_gain > K_i
       (recovers p*B > K_i as r -> 0; ties to non-participation).
  satisficing (compliance-target rule): binary {0, target} expected-cost choice.
    S: a in {0, req}: full-target cost C(req) vs expected fine q*req at a=0.
    T: a in {0, 0.2 E0} (common-ambition target) under the charge's expected
       payment at the firm's reporting branch. V: already binary (unchanged).
  normative additive term omega = 0.25*tau* per unit abatement: effective
    abatement price q -> q+omega (S), price -> price+omega (T); V: K_i reduced
    by omega*target in the decision only (realized K_i cost unchanged).
  cost powers phi in {1.5, 3}: C(a) = gamma a^phi / phi (phi=2 is baseline up to
    the 1/2 factor, which the phi-form reproduces exactly at phi=2);
    interior a = (price/gamma)^{1/(phi-1)}, capped; K and threat m recomputed
    under the same form (coherence rule below).
  piecewise-linear cost: marginal cost s1=0.05 gamma E0 for a<0.1E0, s2=0.15
    gamma E0 beyond; corner optima a* = 1{price>s1} min(cap,0.1E0) +
    1{price>s2} max(cap-0.1E0,0).
  evasion off: identical to ablation m1; outputs REUSED (documented, same seeds).
  convex evasion penalty phi_c = f_T: interior under-reporting
    u* = max(0, [tau(1-p) - p f_T]/(2 p phi_c)) capped at E (the recorded
    ex-ante formula); abatement at the truthful margin; charges on E-u*,
    recovery on u* as baseline.

F2 instrument design
  threshold rule variants (standard): uniform absolute cut Ebar = 0.8 (the
    distribution-mean normalizer E[E0]=1), req_i = max(E0_i - 0.8, 0);
    cost-conditioned req_i = 0.2 sum(E0) * w_i/sum(w), w_i = 1/gamma_i,
    capped at E0_i (single cap, no redistribution; per-seed).
  target rule variants (voluntary): uniform absolute target 0.2; cost-
    conditioned target = the cost-conditioned req_i above.
  rep_standard: private cost adds p*rho_i*1{v>0}; firm compares interior-optimum
    total vs full compliance; rho private, outside C_social.
  theta forms: probit theta = Phi(a0 + b0 s) endpoint-matched to the logistic at
    s=0,1; piecewise-linear theta = theta(0) + [theta(1)-theta(0)] s.
  collective preemption: exemption is collective, individual threat term drops;
    participate iff p*nu_i > 1 (free-riding); theta path still tracked.
  within-tick theta: theta_t solves theta = logistic(theta0 + lam s(theta))
    within the tick, s from the SAME period's verified achievement at that
    theta; 60-step bisection on [0,1], smallest fixed point.
  one-period M: Ann -> 1 (v1.1 convention, v1.2 nu). nu members set (nu0,eta)
    or cap: top-decile (0.2920,2.0000); uniform (eta=0) and visibility (eta=1)
    with nu_0 refit by the frozen single-moment rule (driver); tail-cap
    min(nu_i,2).

Coherence rule (documented): any variant altering the standard's cost structure
or threshold recomputes the voluntary threat value m_i(g) under the SAME
variant (the threatened regime is the variant standard); variants touching only
T or only V leave m_i at the frozen baseline M_of. K_i (real cost of the
voluntary target) always uses the variant cost form and variant target.
Realized event accounting (draw streams, fines, charges, kappa counts,
C_social) is unchanged in structure throughout.
"""
import numpy as np
from model_v12 import CONFIG, M_of, logistic, Ann, _base_out, _finish

PHI = None  # scipy imported lazily for probit

def _norm_cdf(x):
    from scipy.stats import norm
    return norm.cdf(x)

# ---------- cost primitives ----------
def cost_of(a, gam, E0, v):
    f = v.get('cost','quad')
    if f=='quad':  return 0.5*gam*a**2
    if f=='power': return gam*np.power(np.maximum(a,0.0), v['phi'])/v['phi']
    if f=='pwl':
        s1,s2,bp = 0.05*gam*E0, 0.15*gam*E0, 0.1*E0
        a1 = np.minimum(a,bp); a2 = np.maximum(a-bp,0.0)
        return s1*a1 + s2*a2
    raise ValueError(f)

def interior_a(price, gam, E0, cap, v):
    f = v.get('cost','quad')
    if f=='quad':  a = price/gam
    elif f=='power': a = np.power(np.maximum(price,0.0)/gam, 1.0/(v['phi']-1.0))
    elif f=='pwl':
        s1,s2,bp = 0.05*gam*E0, 0.15*gam*E0, 0.1*E0
        a = np.where(price>s2, np.maximum(cap,0.0),
            np.where(price>s1, np.minimum(bp,cap), 0.0))
        return np.minimum(a, cap)
    else: raise ValueError(f)
    return np.minimum(a, cap)

def cara_a(exposure_rate, prob, gam, E0, cap, r, v):
    """min_a C(a) + (1/r) ln(1-prob+prob e^{r*rate*(cap-a)}); bisection on dCE/da=0.
    Marginal benefit of a: rate*prob*e^{rz}/(1-prob+prob*e^{rz}), z=r*rate*(cap-a),
    decreasing in a; marginal cost C'(a) increasing -> unique interior root."""
    lo = np.zeros_like(cap); hi = np.maximum(cap,0.0)
    for _ in range(40):
        mid = 0.5*(lo+hi)
        z = np.clip(r*exposure_rate*(cap-mid), 0.0, 500.0)
        ez = np.exp(z)
        mb = exposure_rate*prob*ez/(1.0-prob+prob*ez)
        mc = (gam*mid if v.get('cost','quad')=='quad' else
              gam*np.power(np.maximum(mid,1e-300), v['phi']-1.0) if v.get('cost')=='power'
              else np.where(mid<0.1*E0, 0.05*gam*E0, 0.15*gam*E0))
        lo = np.where(mc<mb, mid, lo); hi = np.where(mc<mb, hi, mid)
    return 0.5*(lo+hi)

def cara_ce(a, exposure_rate, prob, gam, E0, r, cap, v):
    z = np.clip(r*exposure_rate*(cap-a), 0.0, 500.0)
    return cost_of(a,gam,E0,v) + (1.0/r)*np.log(1.0-prob+prob*np.exp(z))

# ---------- required abatement under the standard ----------
def std_req(gam, E0, v):
    tr = v.get('threshold','base')
    if tr=='base':    return 0.2*E0
    if tr=='uniform': return np.maximum(E0-0.8, 0.0)
    if tr=='cost':
        w = 1.0/gam
        req = 0.2*E0.sum(1, keepdims=True)*w/w.sum(1, keepdims=True)
        return np.minimum(req, E0)
    raise ValueError(tr)

# ---------- standard ----------
def run_standard_R(bat, g, tau_s, cfg, v):
    p,c,k = g['p'],g['c'],g['k']; f_S = cfg['m_S']*tau_s
    gam,E0 = bat['gamma'],bat['E0']; S,N = gam.shape; T = cfg['T']
    req = std_req(gam,E0,v)
    q = p*c*(1.0-k)*f_S + v.get('omega',0.0)
    if v.get('rule')=='cara':
        a = cara_a((1.0-k)*f_S, p*c, gam, E0, req, v['r'], v)
    elif v.get('rule')=='satisficing':
        full = cost_of(req,gam,E0,v); zero = q*req
        a = np.where(full<=zero, req, 0.0)
    elif v.get('rule')=='rep_standard':
        ai = interior_a(q,gam,E0,req,v)
        rho = (cfg['nu0']*np.power(E0,cfg['eta'])) * (0.02*gam*E0**2)
        tot_i = cost_of(ai,gam,E0,v) + q*np.maximum(req-ai,0.0) + p*rho*((req-ai)>1e-15)
        tot_f = cost_of(req,gam,E0,v)
        a = np.where(tot_f<=tot_i, req, ai)
    else:
        a = interior_a(q,gam,E0,req,v)
    vio = np.maximum(req-a, 0.0)
    insp = bat['U_mon'] < p
    violm = (vio>1e-15)[:,None,:]
    impl = insp & violm & (bat['U_imp'] < c)
    o = _base_out(S)
    o['A'] = T*a.sum(1); o['C_ab'] = T*cost_of(a,gam,E0,v).sum(1)
    o['nM'] = insp.sum((1,2)).astype(float); o['nA'] = float(N*T)*np.ones(S)
    o['nI'] = impl.sum((1,2)).astype(float)
    o['T_fi'] = ((1-k)*f_S*vio[:,None,:]*impl).sum((1,2))
    o['viol'] = (vio>1e-15).mean(1); o['under'] = np.zeros(S)
    _finish(o,cfg)
    # variant one-period standard cost (threat value), same variant machinery
    q0 = p*c*(1.0-k)*f_S
    a0 = interior_a(q0,gam,E0,req,v)
    o['m_threat'] = cost_of(a0,gam,E0,v) + q0*np.maximum(req-a0,0.0)
    return o

# ---------- charge ----------
def run_charge_R(bat, g, tau_s, cfg, v):
    p,c,k = g['p'],g['c'],g['k']; tau = tau_s; f_T = cfg['m_T']*tau_s
    gam,E0 = bat['gamma'],bat['E0']; S,N = gam.shape; T = cfg['T']
    omg = v.get('omega',0.0)
    if v.get('convex_evasion'):
        phi_c = f_T
        price = c*(1.0-k)*tau + omg
        a = interior_a(price,gam,E0,E0,v); E = E0-a
        if p > 0:
            ustar = max((tau*(1.0-p)-p*f_T)/(2.0*p*phi_c), 0.0)
            u = np.minimum(np.full_like(E, ustar), E)
        else:
            u = E.copy()
        Etil = E-u
    elif v.get('rule')=='cara':
        r = v['r']
        aT = cara_a((1.0-k)*tau, c, gam, E0, E0, r, v)
        ceT = cara_ce(aT,(1.0-k)*tau,c,gam,E0,r,E0,v)
        aE = cara_a((1.0-k)*(tau+f_T), p*c, gam, E0, E0, r, v)
        ceE = cara_ce(aE,(1.0-k)*(tau+f_T),p*c,gam,E0,r,E0,v)
        tru = ceT<=ceE
        a = np.where(tru,aT,aE); E = E0-a; Etil = np.where(tru,E,0.0); u = E-Etil
    else:
        thr = tau/(tau+f_T) if (tau+f_T)>0 else 0.0
        truthful = (p >= thr-1e-12)
        price = (c*(1.0-k)*tau if truthful else p*c*(1.0-k)*(tau+f_T)) + omg
        if v.get('rule')=='satisficing':
            tgt = 0.2*E0
            eff = c*(1.0-k)*tau if truthful else p*c*(1.0-k)*(tau+f_T)
            pay0, payt = eff*E0, eff*(E0-tgt)
            a = np.where(cost_of(tgt,gam,E0,v)+payt <= pay0, tgt, 0.0)
        else:
            a = interior_a(price,gam,E0,E0,v)
        E = E0-a; Etil = E if truthful else np.zeros_like(E); u = E-Etil
    insp = bat['U_mon'] < p; coll = bat['U_imp'] < c
    det = insp & (u>1e-15)[:,None,:]; rec = det & (bat['U_imp2'] < c)
    o = _base_out(S)
    o['A'] = T*a.sum(1); o['C_ab'] = T*cost_of(a,gam,E0,v).sum(1)
    o['nM'] = insp.sum((1,2)).astype(float); o['nA'] = float(N*T)*np.ones(S)
    o['nI'] = rec.sum((1,2)).astype(float)
    o['T_ch'] = ((1-k)*tau*Etil[:,None,:]*coll).sum((1,2)) + ((1-k)*tau*u[:,None,:]*rec).sum((1,2))
    o['T_fi'] = ((1-k)*f_T*u[:,None,:]*rec).sum((1,2))
    o['under'] = T*u.sum(1)/np.maximum(T*(E0-a).sum(1),1e-12)
    o['viol'] = np.zeros(S)
    _finish(o,cfg); return o

# ---------- voluntary ----------
def vol_target(gam, E0, v):
    tg = v.get('target','base')
    if tg=='base':    return 0.2*E0
    if tg=='uniform': return 0.2*np.ones_like(E0)
    if tg=='cost':
        w = 1.0/gam
        return np.minimum(0.2*E0.sum(1,keepdims=True)*w/w.sum(1,keepdims=True), E0)
    raise ValueError(tg)

def theta_of(s, cfg, v):
    tf = v.get('theta_form','logistic')
    if tf=='logistic': return logistic(cfg['theta0']+cfg['lam']*s)
    t0 = logistic(cfg['theta0']); t1 = logistic(cfg['theta0']+cfg['lam'])
    if tf=='linear': return t0 + (t1-t0)*s
    if tf=='probit':
        from scipy.stats import norm
        a0 = norm.ppf(t0); b0 = norm.ppf(t1)-a0
        return _norm_cdf(a0+b0*s)
    raise ValueError(tf)

def run_voluntary_R(bat, g, tau_s, cfg, v, m_threat=None):
    p,c,k,h = g['p'],g['c'],g['k'],g['h']; f_S = cfg['m_S']*tau_s
    gam,E0 = bat['gamma'],bat['E0']; S,N = gam.shape; T = cfg['T']; dl = cfg['delta']
    tgt = vol_target(gam,E0,v)
    K = cost_of(tgt,gam,E0,v)
    nu0 = v.get('nu0', cfg['nu0']); eta = v.get('eta', cfg['eta'])
    nu = nu0*np.power(E0, eta)
    if v.get('nu_cap') is not None: nu = np.minimum(nu, v['nu_cap'])
    rho = nu*K
    m = m_threat if m_threat is not None else M_of(p,c,k,f_S,gam,E0)[0]
    mK = np.where(K>0, m/K, 0.0)
    Kdec = K - v.get('omega',0.0)*tgt
    Astar_seed = cfg['A_star_frac']*E0.sum(1)
    o = _base_out(S); theta_path=[]; s_path=[]; Rp=np.zeros(S)
    ann1 = v.get('one_period_M', False)
    theta = np.full(S, logistic(cfg['theta0']) if v.get('theta_form','logistic')=='logistic'
                    else theta_of(0.0,cfg,v))

    def part_mask(th_vec, t):
        AnnT = 1.0 if ann1 else Ann(dl*h, T-(t+1))
        if v.get('collective'):
            return p*nu > 1.0
        if v.get('rule')=='cara':
            r = v['r']
            Bmon = rho + dl*h*th_vec[:,None]*m*AnnT
            ce = -(1.0/r)*np.log(1.0-p+p*np.exp(np.clip(-r*Bmon,-500.0,0.0)))
            return ce > Kdec
        if v.get('omega',0.0) > 0.0:
            return p*(rho + dl*h*th_vec[:,None]*m*AnnT) > Kdec
        # exact frozen v1.2 expression (ties to non-participation)
        return p*(nu + dl*h*th_vec[:,None]*mK*AnnT) > 1.0

    for t in range(T):
        if v.get('within_tick'):
            # theta -> s(theta) nonincreasing, theta_of nondecreasing => unique fixed point
            lov = np.zeros(S); hiv = np.ones(S)
            for _ in range(25):
                mid = 0.5*(lov+hiv)
                part = part_mask(mid, t)
                z = part & (bat['U_mon'][:,t,:] < p)
                s = np.maximum((Astar_seed-(tgt*z).sum(1))/Astar_seed, 0.0)
                f = theta_of(s,cfg,v) - mid
                lov = np.where(f>0, mid, lov); hiv = np.where(f>0, hiv, mid)
            theta = 0.5*(lov+hiv)
        part = part_mask(theta, t)
        a = tgt*part
        z = part & (bat['U_mon'][:,t,:] < p)
        Ahat = (tgt*z).sum(1)
        s = np.maximum((Astar_seed-Ahat)/Astar_seed, 0.0)
        theta_path.append(theta.copy()); s_path.append(s.copy())
        if not v.get('within_tick'):
            theta = theta_of(s,cfg,v)
        Rp += (rho*z).sum(1)
        o['A'] += a.sum(1); o['C_ab'] += cost_of(a,gam,E0,v).sum(1)
        o['nM'] += z.sum(1); o['nA'] += part.sum(1); o['part'] += part.mean(1)/T
    o['theta_path']=np.array(theta_path); o['s_path']=np.array(s_path); o['R_paid']=Rp
    _finish(o,cfg); return o

STD_COUPLED = ('threshold','cost','rule','omega')   # keys that alter the standard

def run_all_R(bat, g, tau_s, cfg=CONFIG, v=None):
    v = v or {}
    oS = run_standard_R(bat,g,tau_s,cfg,v)
    oT = run_charge_R(bat,g,tau_s,cfg,v)
    mth = oS.pop('m_threat')
    couple = any(kk in v for kk in STD_COUPLED)
    oV = run_voluntary_R(bat,g,tau_s,cfg,v, m_threat=(mth if couple else None))
    return dict(S=oS,T=oT,V=oV)
