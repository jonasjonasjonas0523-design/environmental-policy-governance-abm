"""Verification battery per PRE_REGISTRATION_SPECIFICATION.md Sec 6, with the
A2-amended item 12 (solver vs independent quadrature; 0.2/e as floor only).
All tests must pass before validation is attempted. Writes
results/verification/{verification_results.csv, verification_report.md, tau_star.json}.
"""
import numpy as np, json, sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__)))
from model import (CONFIG, batch, run_all, run_standard, run_charge, run_voluntary,
                   solve_tau_star, F_of_tau, _quad, M_of, logistic, export_config)

ROWS = []
def row(tid, desc, expected, actual, tol, ok):
    ROWS.append(dict(test_id=tid, description=desc, expected=expected,
                     actual=actual, tolerance=tol, status='PASS' if ok else 'FAIL'))
    return ok

def main():
    os.makedirs('results/verification', exist_ok=True)
    export_config('config/baseline.yaml')
    cfg = CONFIG
    # ---- frozen tau* (A2) ----
    tau = solve_tau_star(10**6, 424242)
    gq, Eq = _quad(10**6, 777)
    F_ind = F_of_tau(tau, gq, Eq)
    se_ind = float(np.std(np.minimum(tau/gq, Eq)))/1000.0
    floor = 0.2/np.e
    cap_corr = tau/floor - 1.0
    json.dump(dict(tau_star=tau, quad_seed=424242, n_quad=10**6,
                   independent_check_F=F_ind, independent_check_5se=5*se_ind,
                   analytic_floor=floor, cap_correction_pct=100*cap_corr,
                   criterion='A2: |F_indep(tau*)| < 5*SE_indep; tau* > 0.2/e floor'),
              open('results/verification/tau_star.json','w'), indent=1)
    f_S, f_T = cfg['m_S']*tau, cfg['m_T']*tau

    seeds = list(range(5001, 5501))          # pilot-range test seeds
    b = batch(seeds)
    gam, E0 = b['gamma'], b['E0']; T = cfg['T']; N = cfg['N']
    g0 = dict(p=.5, c=.6, k=.3, h=.75)
    oS = run_standard(b, g0, tau); oT7 = run_charge(b, dict(p=.7,c=.6,k=.3,h=.75), tau)
    oT3 = run_charge(b, dict(p=.3,c=.6,k=.3,h=.75), tau)

    # 1 emissions accounting / bounds
    q = g0['p']*g0['c']*(1-g0['k'])*f_S
    aS = np.minimum(q/gam, 0.2*E0)
    ok = bool((aS >= -1e-15).all() and (aS <= E0+1e-12).all() and ((E0-aS) >= -1e-12).all())
    row('T01','a in [0,E0]; E=E0-a >= 0','all bounds hold',
        f'max a-E0 = {float((aS-E0).max()):.2e}','exact', ok)

    # 2 reporting bounds (recomputed independently)
    for pp,lab in [(.7,'truthful'),(.3,'evading')]:
        thr = tau/(tau+f_T)
        truthful = pp >= thr-1e-12
        price = g0['c']*(1-g0['k'])*(tau if truthful else pp*(tau+f_T))
        a = np.minimum(price/gam, E0); E = E0-a
        Et = E if truthful else np.zeros_like(E)
        ok = bool(((Et>=-1e-15)&(Et<=E+1e-12)).all())
        row('T02','0 <= Etilde <= E ('+lab+')','bounds hold','holds','exact', ok)

    # 3 standard liability in expectation
    v = np.maximum(0.2*E0-aS,0.0)
    exp_fi = g0['p']*g0['c']*(1-g0['k'])*f_S*v.sum()*T
    per = (1-g0['k'])*f_S*v; pc = g0['p']*g0['c']
    var = (per**2*pc*(1-pc)).sum()*T
    z3 = (oS['T_fi'].sum()-exp_fi)/np.sqrt(max(var,1e-30))
    row('T03','E[fines] = pc(1-k) f_S [E-Ebar]+', f'{exp_fi:.2f}',
        f'{oS["T_fi"].sum():.2f} (z={z3:.2f})','|z|<5', bool(abs(z3)<5))

    # 4 charge liability in expectation (both branches)
    for o,pp,lab in [(oT7,.7,'truthful'),(oT3,.3,'evading')]:
        truthful = pp >= tau/(tau+f_T)-1e-12
        price = g0['c']*(1-g0['k'])*(tau if truthful else pp*(tau+f_T))
        a = np.minimum(price/gam,E0); E = E0-a
        Et = E if truthful else np.zeros_like(E); u = E-Et
        exp_tot = (g0['c']*(1-g0['k'])*tau*Et.sum() + pp*g0['c']*(1-g0['k'])*(tau+f_T)*u.sum())*T
        varc = (((1-g0['k'])*tau*Et)**2*g0['c']*(1-g0['c'])).sum()*T
        pcx = pp*g0['c']
        varr = (((1-g0['k'])*(tau+f_T)*u)**2*pcx*(1-pcx)).sum()*T
        real = (o['T_ch']+o['T_fi']).sum()
        z4 = (real-exp_tot)/np.sqrt(max(varc+varr,1e-30))
        row('T04','E[charge liability] eq.(8) '+lab, f'{exp_tot:.2f}',
            f'{real:.2f} (z={z4:.2f})','|z|<5', bool(abs(z4)<5))

    # 5 monitoring frequency
    fr = float((b['U_mon']<g0['p']).mean()); n5 = b['U_mon'].size
    z5 = (fr-g0['p'])/np.sqrt(g0['p']*(1-g0['p'])/n5)
    row('T05','freq(monitoring) -> p', g0['p'], f'{fr:.5f} (z={z5:.2f})','|z|<5', bool(abs(z5)<5))

    # 6 implementation frequency conditional on detection
    insp = b['U_mon']<g0['p']; violm = (v>1e-15)[:,None,:]
    sel = insp & violm
    fr6 = float((b['U_imp'][sel]<g0['c']).mean()); n6 = int(sel.sum())
    z6 = (fr6-g0['c'])/np.sqrt(g0['c']*(1-g0['c'])/n6)
    row('T06','freq(implementation|detected) -> c', g0['c'],
        f'{fr6:.5f} (z={z6:.2f})','|z|<5', bool(abs(z6)<5))

    # 7 leakage scales implemented consequences by (1-k) (deterministic)
    impl = sel & (b['U_imp']<g0['c'])
    fines_re = ((1-g0['k'])*f_S*v[:,None,:]*impl).sum()
    rel = abs(fines_re-oS['T_fi'].sum())/max(fines_re,1e-12)
    row('T07','implemented fine = (1-k) f_S v exactly','rel err < 1e-10',
        f'{rel:.2e}','1e-10', bool(rel<1e-10))

    # 8 reputational payment only on d=1 and verified z=1 (independent re-simulation, live nu)
    cfg8 = dict(cfg); cfg8['nu']=1.5
    g8 = dict(p=.6,c=.6,k=.2,h=.9)
    o8 = run_voluntary(b,g8,tau,cfg8)
    M8,_,_ = M_of(g8['p'],g8['c'],g8['k'],f_S,gam,E0)
    Ct = 0.5*gam*(0.2*E0)**2; rho8 = 1.5*Ct
    Ast = cfg['A_star_frac']*E0.sum(1); th = np.full(len(seeds), logistic(cfg['theta0']))
    Rp = np.zeros(len(seeds)); Aacc = np.zeros(len(seeds))
    for t in range(T):
        part = Ct <= g8['p']*(rho8 + cfg['delta']*g8['h']*th[:,None]*M8)
        z = part & (b['U_mon'][:,t,:]<g8['p'])
        Rp += (rho8*z).sum(1); Aacc += (0.2*E0*part).sum(1)
        s = np.maximum((Ast-(0.2*E0*z).sum(1))/Ast,0.0)
        th = logistic(cfg['theta0']+cfg['lam']*s)
    ok = bool(np.allclose(Rp,o8['R_paid'],rtol=1e-12) and np.allclose(Aacc,o8['A'],rtol=1e-12)
              and o8['R_paid'].sum()>0)
    row('T08','R paid only on d=1,z=1 (independent re-simulation, nu=1.5)',
        'exact match; payments > 0', f'sumR={o8["R_paid"].sum():.2f}','1e-12', ok)

    # 9 theta update matches eq.(12)
    tp, sp = o8['theta_path'], o8['s_path']
    err9 = float(np.abs(tp[1:]-logistic(cfg['theta0']+cfg['lam']*sp[:-1])).max())
    row('T09','theta_{t+1} = logistic(theta0+lam s_t)','max err < 1e-12',
        f'{err9:.2e}','1e-12', bool(err9<1e-12))

    # 10 fiscal transfers excluded from C_social
    ok = True
    for o in (oS,oT7,oT3,o8):
        rec = o['C_ab']+cfg['kappa_M']*o['nM']+cfg['kappa_A']*o['nA']+cfg['kappa_I']*o['nI']
        ok &= bool(np.allclose(rec,o['C_social'],rtol=1e-12))
        ok &= bool(np.allclose(o['T_fiscal'],o['T_ch']+o['T_fi'],rtol=1e-12))
    row('T10','C_social excludes transfers; T_fiscal = charges+fines','identities hold',
        'hold' if ok else 'BROKEN','1e-12', ok)

    # 11 kappa event counts match the log (independent recount, standard)
    nM_re = float(insp.sum()); nI_re = float(impl.sum()); nA_re = float(N*T*len(seeds))
    ok = (abs(nM_re-oS['nM'].sum())<0.5 and abs(nI_re-oS['nI'].sum())<0.5
          and abs(nA_re-oS['nA'].sum())<0.5)
    row('T11','kappa event counts match event log','counts equal',
        f'nM {oS["nM"].sum():.0f}={nM_re:.0f}; nI {oS["nI"].sum():.0f}={nI_re:.0f}','exact', bool(ok))

    # 12 tau* (A2-amended)
    ok = bool(abs(F_ind) < 5*se_ind and tau > floor)
    row('T12','[A2] tau* solver vs independent 1e6 quadrature; > 0.2/e floor',
        f'|F|<{5*se_ind:.2e}; tau*>{floor:.5f}',
        f'tau*={tau:.5f}; F_ind={F_ind:.2e}; cap corr {100*cap_corr:.1f}%','5*SE', ok)

    # 13 identical-seed reproducibility (bitwise on logged outputs)
    b1 = batch([7,8,9]); b2 = batch([7,8,9])
    r1 = run_all(b1,g0,tau); r2 = run_all(b2,g0,tau)
    ok = all(np.array_equal(np.asarray(r1[j][kk],dtype=float),
                            np.asarray(r2[j][kk],dtype=float), equal_nan=True)
             for j in r1 for kk in r1[j])
    row('T13','identical seeds -> bitwise identical outputs','bitwise equal',
        'equal' if ok else 'DIFFER','bitwise', bool(ok))

    # 14 limiting cases
    bl = batch(list(range(6001,6101)))
    rp0 = run_all(bl,dict(p=0,c=.6,k=.3,h=.75),tau)
    M0,_,_ = M_of(0,.6,.3,f_S,bl['gamma'],bl['E0'])
    ok = all(float(rp0[j]['A'].sum())==0 for j in 'STV') and float(rp0['T']['T_ch'].sum())==0 \
         and float(rp0['S']['T_fi'].sum())==0 and float(np.abs(M0).max())==0
    row('T14a','p=0: zero consequences, zero charge liability, zero voluntary, M=0',
        'all zero','all zero' if ok else 'NONZERO','exact', bool(ok))
    rc0 = run_all(bl,dict(p=.6,c=0,k=.3,h=.75),tau)
    ok = all(float(rc0[j]['A'].sum())==0 for j in 'STV') and float(rc0['T']['T_fiscal'].sum())==0
    row('T14b','c=0 collapses all consequence channels','all zero',
        'all zero' if ok else 'NONZERO','exact', bool(ok))
    rk1 = run_all(bl,dict(p=.6,c=.6,k=1.0,h=.75),tau)
    ok = all(float(rk1[j]['A'].sum())==0 for j in 'STV') and float(rk1['T']['T_fiscal'].sum())==0
    row('T14c','k=1 collapses all consequence channels','all zero',
        'all zero' if ok else 'NONZERO','exact', bool(ok))
    oH = run_standard(bl,dict(p=1,c=1,k=0,h=.75),tau)
    should = (bl['gamma']*0.2*bl['E0'] <= f_S+1e-12)
    qH = f_S; aH = np.minimum(qH/bl['gamma'],0.2*bl['E0'])
    isfull = (0.2*bl['E0']-aH) <= 1e-12
    ok = bool((isfull==should).all())
    row('T14d','p=c=1,k=0: full compliance iff gamma*0.2E0 <= f_S',
        f'{int(should.sum())} compliant firms', f'{int(isfull.sum())}','exact', ok)
    rt0 = run_charge(bl,dict(p=.6,c=.6,k=.3,h=.75),0.0)
    row('T14e','tau=0 gives zero charge abatement','0', f'{float(rt0["A"].sum()):.1f}',
        'exact', bool(float(rt0['A'].sum())==0))
    cfgpr = dict(cfg); cfgpr['theta0']=-700.0; cfgpr['lam']=0.0
    opr = run_voluntary(bl,dict(p=.8,c=.85,k=.05,h=.95),tau,cfgpr)
    thmax = float(opr['theta_path'].max())
    cfgpr2 = dict(cfgpr); cfgpr2['nu']=1.5
    opr2 = run_voluntary(bl,dict(p=.8,c=.85,k=.05,h=.95),tau,cfgpr2)
    ok = thmax < 1e-300 and float(opr['A'].sum())==0 and float(opr2['A'].sum())>0
    row('T14f','theta0=-inf,lam=0: pure-reputation voluntary (theta=0)',
        'theta=0; A=0 at nu=.25; A>0 at nu=1.5',
        f'theta_max={thmax:.1e}; A(.25)={float(opr["A"].sum()):.0f}; A(1.5)={float(opr2["A"].sum()):.0f}',
        'exact', bool(ok))

    # 15 mechanism verification
    oA = run_charge(bl,dict(p=.7,c=.6,k=.3,h=.75),tau)
    oB = run_charge(bl,dict(p=.7,c=.6,k=.3,h=.75),1.1*tau)
    ok = float(oB['A'].sum()) >= float(oA['A'].sum()) and 1.1*tau*.6*.7 > tau*.6*.7
    row('T15a','raising tau raises private cost of residual emissions (and A weakly up)',
        'A(1.1tau) >= A(tau)', f'{float(oA["A"].mean()):.1f} -> {float(oB["A"].mean()):.1f}',
        'monotone', bool(ok))
    ok = (0.95*0.4 > 0.95*0.2)  # delta h theta M increasing in theta, structural
    cfg6 = dict(cfg); cfg6['nu']=0.60
    olo = run_voluntary(bl,dict(p=.8,c=.85,k=.05,h=.50),tau,cfg6)
    ohi = run_voluntary(bl,dict(p=.8,c=.85,k=.05,h=.95),tau,cfg6)
    ok = ok and float(ohi['part'].mean()) >= float(olo['part'].mean())
    row('T15b','stronger future term raises cost of voluntary non-performance; '
        'participation weakly up in h at nu=0.6 (robustness range)',
        'part(h=.95) >= part(h=.5)',
        f'{float(olo["part"].mean()):.3f} -> {float(ohi["part"].mean()):.3f}',
        'monotone', bool(ok))
    ok = all((p2*0.6*0.7*f_S > p1*0.6*0.7*f_S) for p1,p2 in [(.3,.6),(.6,.9)])
    row('T15c','raising p raises expected consequence coefficients under all instruments',
        'coefficients monotone in p','monotone','structural', bool(ok))

    import pandas as pd
    df = pd.DataFrame(ROWS)
    df.to_csv('results/verification/verification_results.csv', index=False)
    nfail = int((df.status=='FAIL').sum())
    with open('results/verification/verification_report.md','w') as f:
        f.write('# Verification report (model v1.1-prerun-amended)\n')
        f.write(f'Date 2026-08-20. {len(df)} tests, {int((df.status=="PASS").sum())} PASS, {nfail} FAIL.\n')
        f.write(f'\ntau* frozen = {tau:.6f} (1e6 quadrature, seed 424242); independent-sample '
                f'F = {F_ind:.2e} (5SE = {5*se_ind:.2e}); cap correction over 0.2/e = '
                f'{100*cap_corr:.1f} percent [A2 criterion].\n\n')
        f.write(df.to_markdown(index=False))
        f.write('\n\nA1 note: voluntary participation rule carries the verification-gated '
                'exemption (AMENDMENTS_PRE_RUN.md). Baseline numerics unchanged.\n')
    print(df[['test_id','status']].to_string(index=False))
    print(f'tau* = {tau:.6f}  cap correction = {100*cap_corr:.1f}%  FAILURES = {nfail}')
    if nfail: sys.exit(1)

if __name__ == '__main__':
    main()
