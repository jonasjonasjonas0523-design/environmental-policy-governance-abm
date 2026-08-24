"""Ablation variants of the v1.2 instrument functions (PRE_REG Sec 7). The frozen
model_v12.py is NOT modified; the two mechanism switches that cannot be expressed
through CONFIG are implemented here as minimal, explicitly documented diffs:

  m1 'no_evasion'            (Etilde = E): run_charge forces the truthful branch
                             regardless of p vs thr; u = 0 by construction. Other
                             instruments have no evasion channel.
  m2 'no_collection_friction' (c=1, k=0 in COLLECTION, evasion untouched): every
                             collection/implementation role of (c, k) is replaced
                             by (1, 0) in all three instruments -- the standard's
                             expected-penalty optimum M_of(p,1,0,...), its
                             implementation draw and fine leakage; the charge's
                             expected price, routine collection, recovery draw and
                             leakage; the voluntary threat value M_of(p,1,0,...).
                             The evasion RULE (thr comparison, u accounting,
                             detection by inspection) is untouched.

Diffs are line-level minimal against model_v12.run_standard / run_charge /
run_voluntary; everything else is imported unchanged.
"""
import numpy as np
from model_v12 import (CONFIG, M_of, logistic, Ann, _base_out, _finish)

def run_standard_abl(bat, g, tau_s, cfg=CONFIG, no_collection_friction=False, **_):
    p,c,k = g['p'],g['c'],g['k']; f_S = cfg['m_S']*tau_s
    if no_collection_friction: c,k = 1.0, 0.0                       # m2
    gam,E0 = bat['gamma'],bat['E0']; S,N = gam.shape; T = cfg['T']
    M,a,q = M_of(p,c,k,f_S,gam,E0)
    v = np.maximum(0.2*E0 - a, 0.0); E = E0 - a
    insp = bat['U_mon'] < p
    violm = (v > 1e-15)[:,None,:]
    impl = insp & violm & (bat['U_imp'] < c)
    o = _base_out(S)
    o['A']    = T*a.sum(1)
    o['C_ab'] = T*(0.5*gam*a**2).sum(1)
    o['nM']   = insp.sum((1,2)).astype(float)
    o['nA']   = float(N*T)*np.ones(S)
    o['nI']   = impl.sum((1,2)).astype(float)
    o['T_fi'] = ((1-k)*f_S*v[:,None,:]*impl).sum((1,2))
    o['viol'] = (v>1e-15).mean(1)
    o['under']= np.zeros(S)
    _finish(o,cfg); return o

def run_charge_abl(bat, g, tau_s, cfg=CONFIG, no_evasion=False, no_collection_friction=False, **_):
    p,c,k = g['p'],g['c'],g['k']; tau=tau_s; f_T=cfg['m_T']*tau_s
    cc,kc = (1.0,0.0) if no_collection_friction else (c,k)          # m2 collection roles
    gam,E0 = bat['gamma'],bat['E0']; S,N = gam.shape; T = cfg['T']
    thr = tau/(tau+f_T) if (tau+f_T) > 0 else 0.0
    truthful = True if no_evasion else (p >= thr - 1e-12)           # m1
    price = cc*(1-kc)*tau if truthful else p*cc*(1-kc)*(tau+f_T)
    a = np.minimum(price/gam, E0); E = E0-a
    Etil = E if truthful else np.zeros_like(E); u = E-Etil
    insp = bat['U_mon'] < p
    coll = bat['U_imp'] < cc
    det  = insp & (u > 1e-15)[:,None,:]
    rec  = det & (bat['U_imp2'] < cc)
    o = _base_out(S)
    o['A']=T*a.sum(1); o['C_ab']=T*(0.5*gam*a**2).sum(1)
    o['nM']=insp.sum((1,2)).astype(float); o['nA']=float(N*T)*np.ones(S)
    o['nI']=rec.sum((1,2)).astype(float)
    o['T_ch']=((1-kc)*tau*Etil[:,None,:]*coll).sum((1,2)) + ((1-kc)*tau*u[:,None,:]*rec).sum((1,2))
    o['T_fi']=((1-kc)*f_T*u[:,None,:]*rec).sum((1,2))
    denom = np.maximum(T*E.sum(1),1e-12)
    o['under']=T*u.sum(1)/denom
    o['viol']=np.zeros(S)
    _finish(o,cfg); return o

def run_voluntary_abl(bat, g, tau_s, cfg=CONFIG, no_collection_friction=False, **_):
    p,c,k,h = g['p'],g['c'],g['k'],g['h']; f_S=cfg['m_S']*tau_s
    if no_collection_friction: c,k = 1.0, 0.0                       # m2 threat value
    gam,E0 = bat['gamma'],bat['E0']; S,N = gam.shape; T=cfg['T']; dl=cfg['delta']
    M,_,_ = M_of(p,c,k,f_S,gam,E0)
    Ct  = 0.5*gam*(0.2*E0)**2
    nui = cfg['nu0']*E0**cfg['eta']
    rho = nui*Ct
    mK  = M/Ct
    Astar_seed = cfg['A_star_frac']*E0.sum(1)
    theta = np.full(S, logistic(cfg['theta0']))
    o = _base_out(S); theta_path=[]; s_path=[]; Rp=np.zeros(S)
    for t in range(T):
        AnnV = Ann(dl*h, T-(t+1))
        part = p*(nui + dl*h*theta[:,None]*mK*AnnV) > 1.0
        a = 0.2*E0*part
        z = part & (bat['U_mon'][:,t,:] < p)
        Ahat = (0.2*E0*z).sum(1)
        s = np.maximum((Astar_seed-Ahat)/Astar_seed, 0.0)
        theta_path.append(theta.copy()); s_path.append(s.copy())
        theta = logistic(cfg['theta0']+cfg['lam']*s)
        Rp += (rho*z).sum(1)
        o['A']    += a.sum(1)
        o['C_ab'] += (0.5*gam*a**2).sum(1)
        o['nM']   += z.sum(1)
        o['nA']   += part.sum(1)
        o['part'] += part.mean(1)/T
    o['theta_path']=np.array(theta_path)
    o['s_path']=np.array(s_path); o['R_paid']=Rp
    _finish(o,cfg); return o

def run_all_abl(bat, g, tau_s, cfg=CONFIG, **abl):
    return dict(S=run_standard_abl(bat,g,tau_s,cfg,**abl),
                T=run_charge_abl(bat,g,tau_s,cfg,**abl),
                V=run_voluntary_abl(bat,g,tau_s,cfg,**abl))
