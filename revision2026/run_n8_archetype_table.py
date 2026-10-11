"""N8 (R2 point 11 / AQ3-4): assemble the main-text archetype table from the archived
readout files, plus reference-plane rep_standard leader shares from the frozen bands."""
import numpy as np, json, pandas as pd
from common import *
rd=pd.read_csv(f'{RES}/archetypes/archetype_readouts.csv')
ov=pd.read_csv(f'{RES}/archetypes/archetype_crossover_overlap.csv')
gr=pd.read_csv(f'{RES}/archetypes/archetype_local_gradients.csv')
rows=[]
for a in ('A1','A2','A3','A4','A5'):
    rA=rd[(rd.archetype==a)&(rd.metric=='A')].iloc[0]
    rC=rd[(rd.archetype==a)&(rd.metric=='C_social')].iloc[0]
    ovA=ov[(ov.archetype==a)&(ov.metric=='A')]['P_overlap'].max()
    ovC=ov[(ov.archetype==a)&(ov.metric=='C_social')]['P_overlap'].max()
    gA=gr[(gr.archetype==a)&(gr.metric=='A')]
    dom=gA[['dD_dp','dD_dc','dD_dk','dD_dh']].abs().max()
    rows.append(dict(archetype=a,
        A_leader=rA['modal_leader'], P_A=float(rA[f'P_lead_{rA["modal_leader"]}']),
        A_margin_p50=float(rA['margin_p50']), A_margin_p05=float(rA['margin_p05']),
        P_displaced_A=float(rA['P_displaced']),
        C_leader=rC['modal_leader'], P_C=float(rC[f'P_lead_{rC["modal_leader"]}']),
        C_margin_p50=float(rC['margin_p50']),
        P_overlap_A_max=float(ovA), P_overlap_C_max=float(ovC),
        steepest_gradient_dir=str(dom.idxmax()).replace('dD_d','')))
df=pd.DataFrame(rows); df.to_csv(f'{OUT}/n8_archetype_table.csv',index=False)
print(df.to_string(index=False))
z=np.load(f'{RES}/robustness/rep_standard_bands.npz')
anc=anchor('reference')
lab,Al,Cl=leaders(z['lo'],z['hi'])
print('reference-plane rep_standard anchored A-leaders:',
      {j: round(float((Al[anc]==ji+1).mean()),3) for ji,j in enumerate('STV')},
      'labels:', {nm:int((lab==zz).sum()) for zz,nm in enumerate(('unres','S','T','V','mixed'))})
