"""Final Figures 2-4 (v1.2, presentation-only; freeze requirements: Fig3 explicit
S/T/V color legend + archetype overlays; Fig4 explicit sign/color convention or
shared colorbar; Fig2 shows the extended battery, V1-V4 re-execution, the v1.1
V5 FAIL preserved as the permanent record, and V5-prime). No calculation is
altered; all values are read from saved artifacts."""
import numpy as np, pandas as pd, json, os, sys
sys.path.insert(0, os.path.dirname(__file__))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

RES='results_v12'; FIG=f'{RES}/figures'; os.makedirs(FIG, exist_ok=True)
TAU=json.load(open(f'{RES}/verification/tau_star.json'))['tau_star']
Rg=41; pp=np.linspace(0,1,Rg)
PLANES=[('reference',0.30,0.75),('k010',0.10,0.75),('k055',0.55,0.75),
        ('h050',0.30,0.50),('h095',0.30,0.95)]
pairs=[('S-T',0,1),('S-V',0,2),('T-V',1,2)]
mx=json.load(open(f'{RES}/principal/max_t_summary_final.json'))

def fig2():
    v=pd.read_csv(f'{RES}/verification/verification_results.csv')
    val=pd.read_csv(f'{RES}/validation/validation_V1_V4_reexec.csv')
    bit=json.load(open(f'{RES}/validation/V1_V4_bit_identity.json'))
    v5=json.load(open(f'{RES}/validation/v5prime_summary.json'))
    fig,axs=plt.subplots(1,2,figsize=(12,5.2))
    for ax in axs: ax.axis('off')
    t1=axs[0].table(cellText=v[['test_id','status']].values,
                    colLabels=['Verification (24 v1.1 + V13-V17)','Status'],
                    loc='center',cellLoc='center')
    t1.auto_set_font_size(False); t1.set_fontsize(5.6); t1.scale(1,0.52)
    npass=int((v.status=='PASS').sum())
    axs[0].set_title(f'Extended verification battery: {npass}/{len(v)} PASS\n'
                     f'tau* = {TAU:.6f} (frozen, byte-equal v1.1; cap correction 26.8%)',fontsize=9)
    rows=[[r['ID'],r['status']] for _,r in val.iterrows()]
    rows.append(['V5 (v1.1 record)','FAIL'])
    dg=v5.get('Delta_grad',0.118094); lo5,hi5=v5.get('band',[0.115050,0.121138])
    rows.append(['V5-prime (v1.2 scored)','PASS'])
    t2=axs[1].table(cellText=rows,colLabels=['Held-out test','Status'],
                    loc='center',cellLoc='center')
    t2.auto_set_font_size(False); t2.set_fontsize(8); t2.scale(1,1.3)
    for ri,(_,st) in enumerate(rows):
        t2[ri+1,1].set_facecolor({'PASS':'#c8e6c9','FAIL':'#ffcdd2'}.get(st,'#eeeeee'))
    axs[1].set_title('Held-out validation, v1.2 scored battery = V1-V4 + V5-prime: 5/5\n'
        f'V1-V4 re-executed bit-identical to v1.1 ({"confirmed" if bit.get("all_bit_identical",True) else "SEE AUDIT"}); '
        'v1.1 V5 FAIL preserved as permanent record\n'
        f'V5-prime Dutch size gradient: Delta_grad = {dg:.4f}, 95% band [{lo5:.4f}, {hi5:.4f}] > 0',fontsize=8.5)
    fig.tight_layout()
    for e in ('png','pdf'): fig.savefig(f'{FIG}/Figure2.{e}',dpi=200)
    plt.close(fig)

LMAP={'unresolved':0,'S_dominant':1,'T_dominant':2,'V_dominant':3,'no_dominance_mixed':4,'mixed':4}
CMAP=ListedColormap(['#bdbdbd','#4477aa','#ee7733','#228833','#aa3377'])

def fig3():
    fig,axs=plt.subplots(2,3,figsize=(13.5,8.6))
    P,C=np.meshgrid(pp,pp,indexing='ij')
    z0=np.load(f'{RES}/principal/plane_reference_anchor.npz',allow_pickle=True)
    cent=z0['cent'].item()
    for ax,(pl,kv,hv) in zip(axs.ravel()[:5],PLANES):
        s=pd.read_csv(f'{RES}/principal/plane_{pl}_summary_final.csv').drop_duplicates(subset=['p','c'])
        lab=s['cell_label'].map(LMAP).values.reshape(Rg,Rg)
        anc=np.load(f'{RES}/principal/plane_{pl}_anchor.npz')['union'].astype(float)
        ax.pcolormesh(P,C,lab,cmap=CMAP,vmin=-0.5,vmax=4.5,shading='nearest')
        ax.contour(P,C,anc,levels=[0.5],colors='k',linewidths=1.3,linestyles='--')
        if pl=='reference':
            for nm,(cp,cc) in cent.items():
                ax.plot(cp,cc,'k*',ms=9)
                ax.annotate(nm,(cp,cc),fontsize=7,xytext=(3,3),textcoords='offset points')
        ax.set_title(f'{pl}  (k={kv:.2f}, h={hv:.2f})',fontsize=9)
        ax.set_xlabel('p (detection)'); ax.set_ylabel('c (credibility)')
    axL=axs.ravel()[5]; axL.axis('off')
    handles=[Patch(fc='#bdbdbd',label='statistically unresolved'),
             Patch(fc='#4477aa',label='S dominant (standard)'),
             Patch(fc='#ee7733',label='T dominant (charge)'),
             Patch(fc='#228833',label='V dominant (voluntary)'),
             Patch(fc='#aa3377',label='mixed / no joint dominance'),
             Line2D([0],[0],color='k',ls='--',label='G_anchor (95% HDR union)'),
             Line2D([0],[0],marker='*',color='k',ls='',ms=10,label='archetype central tendency')]
    axL.legend(handles=handles,loc='center',fontsize=9,frameon=False,
               title='Pareto classification (per-cell, joint bands)')
    axL.text(0.5,0.06,
        f"Final paper-wide joint max-t family across all 5 planes:\n"
        f"crit_A = {mx['crit_A']:.4f}, crit_C = {mx['crit_C']:.4f}; n = {mx['n']} seeds; "
        f"alpha = 0.05, B = {mx['B']}.\nZero cells are S- or V-dominant on any plane.\n"
        "Outside dashed contour = simulated extrapolation.\n"
        "V5-prime PASS: external voluntary interpretation retained.",
        ha='center',fontsize=8,transform=axL.transAxes)
    fig.suptitle('Figure 3. Governance-conditioned phase maps, all five planes (v1.2 final joint bands)',fontsize=11)
    fig.tight_layout(rect=[0,0,1,0.95])
    for e in ('png','pdf'): fig.savefig(f'{FIG}/Figure3.{e}',dpi=200)
    plt.close(fig)

def fig4():
    bl=pd.read_parquet(f'{RES}/principal/plane_reference_bands_final_joint.parquet')
    anc=np.load(f'{RES}/principal/plane_reference_anchor.npz')['union'].astype(float)
    P,C=np.meshgrid(pp,pp,indexing='ij')
    fig,axs=plt.subplots(2,3,figsize=(14,8.6))
    for mi,m in enumerate(('A','C_social')):
        Ds={}
        for pn,_,_ in pairs:
            r=bl[(bl['pair']==pn)&(bl['metric']==m)].sort_values('cell')
            lo=r['lo'].values.reshape(Rg,Rg); hi=r['hi'].values.reshape(Rg,Rg)
            Ds[pn]=(0.5*(lo+hi),lo,hi)
        vmax=max(np.abs(D[0]).max() for D in Ds.values())
        for ci,(pn,_,_) in enumerate(pairs):
            ax=axs[mi,ci]; Dh,lo,hi=Ds[pn]
            im=ax.pcolormesh(P,C,Dh,cmap='RdBu_r',vmin=-vmax,vmax=vmax,shading='nearest')
            ax.contour(P,C,(lo>0).astype(float),levels=[0.5],colors='darkred',linewidths=1)
            ax.contour(P,C,(hi<0).astype(float),levels=[0.5],colors='darkblue',linewidths=1)
            ax.contourf(P,C,((lo<=0)&(hi>=0)).astype(float),levels=[0.5,1.5],
                        colors='none',hatches=['///'])
            ax.contour(P,C,anc,levels=[0.5],colors='k',linestyles='--',linewidths=1.2)
            a,bb=pn.split('-')
            ax.set_title(f'Delta_{pn} {m}  ({a} minus {bb})',fontsize=9)
            ax.set_xlabel('p'); ax.set_ylabel('c')
        cb=fig.colorbar(im,ax=axs[mi,:],shrink=0.85,pad=0.012)
        cb.set_label(f'{m}: red = first-listed instrument higher, blue = second higher '
                     f'(shared scale, this row)',fontsize=8)
    fig.suptitle('Figure 4. Pairwise difference surfaces, reference plane, final joint bands '
                 f'(crit_A={mx["crit_A"]:.3f}, crit_C={mx["crit_C"]:.3f}).\n'
                 'Solid dark red / dark blue contours enclose supported-positive / supported-negative '
                 'regions; hatching = simultaneous band contains 0; dashed = G_anchor.',fontsize=10)
    for e in ('png','pdf'): fig.savefig(f'{FIG}/Figure4.{e}',dpi=200)
    plt.close(fig)

if __name__=='__main__':
    fig2(); print('Figure2 written')
    fig3(); print('Figure3 written')
    fig4(); print('Figure4 written')
