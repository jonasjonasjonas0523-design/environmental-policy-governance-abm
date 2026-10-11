# Revision-stage analyses (2026-10-08), IEEE Access resubmission Access-2026-42785

Status: POST-HOC, reviewer-requested, labeled. The frozen v1.2 baseline record is
untouched; every analysis below runs inside the frozen protocol (model_v12 code,
principal seeds 100000..100152 or validation seeds 90000..90152, rng-77 per-plane
joint max-t families, B=2000, pinned numpy 2.4.4/scipy 1.17.1/pandas 3.0.2).

## 0. Reproduction gates (all PASS before any new run)
- tau* re-solve (quad seed 424242, n=1e6): bit-exact 0.09326697147369385.
- V1-V4 re-execution: bit-identical to the archived record (all four PASS).
- V5-prime: EXACT MATCH to the archived report (9/9 conditions, tol 5e-6).
- Reference plane + staged rng-77 family: lo/hi arrays reproduce to max|diff| 1.4e-14,
  labels identical, crit (3.1101, 3.4758) identical.

## N1. tau anchoring (R2 #1, AQ2-6; his top item)
Frozen anchoring (now to be stated in the main text): full implementation integrity,
p=c=1, k=0, truthful reporting. Alternative anchoring executed: tau_ref solves
E[min(c(1-k) tau / gamma, E0)] = 0.2 at reference governance (c,k)=(0.50,0.30),
truthful at the p=0.50 indifference point; same quadrature protocol.
- tau_ref = 0.2664767158 = 2.8571 x tau* (naive tau*/0.35 = 0.266477; cap correction negligible at this rate).
- Charge arm re-run at tau_ref with f_T = m_T x tau_ref (reporting threshold preserved at 0.5);
  standard and voluntary arms untouched (f_S stays 2 x tau*), isolating the charge's anchoring.
- Reference plane, per-plane families both sides:
  - Charge delivery at the reference cell: 164.4/400 (41%) frozen -> 396.8/400 (99%) re-anchored.
  - Anchored abatement leadership T: 0.547 -> 0.940 (V: 0.268 -> 0.023).
  - Anchored cost leadership S: 0.780 -> 1.000.
  - T_dominant cells: 192 -> 0; mixed 1368 -> 1578; unresolved 121 -> 103. Sign agreement 0.829.
- Reading: the frozen full-integrity anchoring is CONSERVATIVE for the charge's abatement
  leadership (understates it everywhere); the charge's joint-dominance pockets are
  anchoring-dependent and vanish under reference anchoring because the charge's cost
  disadvantage grows with the higher rate. Both directions reported.

## N2. Reputational symmetry on the slices where V leads (R2 #7, AQ2-7)
rep_standard (frozen F2 instantiation: standard's private cost adds p*rho_i*1{v>0},
rho outside C_social) run on k055 and h095; per-plane families; baseline planes
reproduce archived final summaries to 5.7e-14.
- k055: anchored A-leader V 0.985 -> 0.000; S 0.000 -> 1.000. T_dom 98 -> 0. Sign agreement 0.478.
- h095: anchored A-leader V 0.408 -> 0.000, T 0.392 -> 0.000; S -> 1.000. T_dom 90 -> 0. Sign agreement 0.463.
- Reference (frozen robustness record, same construction): S -> 1.000 likewise.
- Reading: the voluntary arm's high-leakage and long-horizon abatement leadership does NOT
  survive granting the standard the same calibrated reputational channel; the reputational
  standard becomes the anchored abatement leader on all three planes examined, while joint
  dominance disappears (mixed/trade-off structure persists). Baseline leadership claims are
  conditional on the baseline's asymmetric reputational assignment; governance-conditioning
  itself persists.

## N3. Correlated heterogeneity (R2 #8, AQ2-8)
batch_corr construction (frozen; r=0 reduces bitwise to the frozen batch). Frozen program
had r = -0.3 (sign agreement 0.9416). Added r = +0.3 and r = -0.6, reference plane:
- +0.3: sign agreement 0.939; anchored A-leader T 0.516 (V 0.279); T_dom 157.
- -0.6: sign agreement 0.890; anchored A-leader T 0.478 (V 0.270; S 0.071); T_dom 144; unresolved 253.
- Reading: the charge's anchored abatement leadership and the map structure are robust to
  cost-size correlation in both directions; stronger negative correlation mainly grows the
  statistically unresolved set.

## N4. Grid-resolution record completed for the V pairs (R2 #6, AQ2-9)
81x81 all-arms reference plane (principal seeds), pilot nearest-matched crossing metric
(tol = 0.5 x finer spacing = 0.00625) + the labeled shared-point 2SE sign diagnostic:
- S-V|A shift 0.00125 (within tol; 3 mismatched lines), T-V|A 0.00125 (0 mm),
  S-V|C 0.00375, T-V|C 0.00375 (both within tol), S-T|A 0.00125;
  S-T|C 0.73375 (the known near-flat-line crossing-pair reconfiguration).
- Shared-point supported-sign agreement 1.0000 with 0 contradictions for ALL six pair|metric
  surfaces.
- Full ladder record for the paper: 6->11 0.0538, 11->21 0.0488, 21->41 0.0238 (pilot),
  41->81 0.001-0.004 for every H3-carrying boundary. Finding: crossovers and supported-sign
  regions are resolution-stable; the only moving object is crossing-pairings on near-flat
  S-T cost surfaces. No formal convergence claim is added.

## N5. Failure-margin analysis of the scored battery (R2 #2, AQ2-2, AQ2-3)
Constructions verbatim from validate_v12; validation seeds; parameter movements measured
against the predeclared robustness ranges (uniform-weighting convention over the box,
tau* held frozen per the F2/F4/Sobol convention).
- V1: dA t = 98.8, dviol t = -73.3; both components monotone in p by construction under
  paired seeds -> structural consistency check, no realistic failure room.
- V2: baseline statistic reproduces exactly (10.026535). Binding test: fails below
  m_S = 1.9812 (a 0.94 percent decrease from the frozen 2.0); over the predeclared box
  (m_S in [1,4], sigma_g and sigma_E in [0.5,1.5]) only 54.4 percent of parameter
  combinations land inside [10,40] (range 4.89-15.27; failures all below 10). The pass
  both constrains the model and was knife-edge at the lower bound.
- V3: sign component t = 188 (structural); level component A_dis = 21.9 vs threshold 100,
  failing only at m_S = 12.1, three times the predeclared ceiling -> near-unfailable.
- V4: A = 358.4 vs threshold 320 at baseline; evasion component structural at p=0.95
  (failure needs m_T < 0.053 vs floor 0.5). Real exposure within the predeclared
  population box: at (sigma_g, sigma_E) = (0.5, 1.5) delivery falls to 181.0 and the
  test FAILS -> partially failable (population structure), not failable via enforcement
  parameters.
- V5-prime: see N6; sign genuinely contested between channels; no magnitude criterion
  (conceded; post-hoc magnitude comparison to be added from the Dutch sources).
Tally for the reframed abstract: two tests with substantial failure exposure (V2, V5-prime),
one with partial exposure (V4, population corner), two structural consistency checks
(V1, V3).

## N6. V5-prime channel decomposition (R2 #3/#4 context)
At the frozen Dutch placement, validation seeds:
- Full model Delta_grad = 0.1181 (intensities top 0.167 / bottom 0.049).
- Reputation-only (theta off): +0.1276 (top 0.128 / bottom 0.000).
- Threat-only (nu0 = 0): -0.0478 NEGATIVE (top 0.032 / bottom 0.080).
Reading: the two channels push the gradient in OPPOSITE directions, exactly as the frozen
memo Sec 2.3 argued ex ante; a threat-dominated mix fails the sign test, so the sign was
not preordained by eta. The magnitude criterion remains absent and is conceded.

## N7. Preregistered directional expectations E1-E4 (R2 #9, AQ2-1)
Map-level check on the frozen/reproduced surfaces: E1 met (S monotone in p and c; zero
violations), E2 met (T monotone in c; T abatement ordered across k010/ref/k055 cellwise),
E3 met (V monotone in p; ordered across h050/ref/h095), E4 met (A1 central: S 236.6,
T 335.9, V 232.6, all above 200 = half the cumulative ambition). Tally: 4 of 4 met,
0 failed.

## N8. Main-text archetype table (R2 #11, AQ3-4)
Assembled from the archived readout files (500 frozen rng-2026 draws per archetype x 153
seeds): per archetype, modal A-leader and probability, median and p05 margin over the
runner-up, displacement probability, modal C-leader and probability, max crossover-band
overlap probability, and steepest local placement-sensitivity direction (c for all five).
A1 T 1.000; A2 T 0.924; A3 T 0.774; A4 V 0.834; A5 T 0.588/V 0.412; C-leader S at all
five (A4 near-tie, median margin 0.010). File: revision2026/out/n8_archetype_table.csv.
