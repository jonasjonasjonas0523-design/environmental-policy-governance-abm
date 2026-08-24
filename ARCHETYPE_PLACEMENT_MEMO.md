# ARCHETYPE_PLACEMENT_MEMO.md
Placement record for Table 3 (Rev2). Date: 2026-08-20. Archetypes are stylized governance configurations, not country predictions; the anchoring jurisdiction named for each is the evidentiary anchor, not a claim about that jurisdiction's future behavior. This qualification travels with every named use, per Section IV-C.

## 1. Decision: five archetypes retained (A5 kept)
The manuscript permitted four or five. Five are retained because the regional evidence supports a configuration that A1-A4 cannot span: a long policy horizon combined with only moderate detection and credibility. Regime durability and enforcement capacity are positively correlated in most global samples, and A1-A4 respect that correlation; Vietnam-type single-party durability decouples h from (p, c), and losing that decoupling would artificially shrink G_anchor exactly where the shadow-of-regulation mechanism (which loads on h) is most distinctive. Retention is therefore evidence-motivated, not row-filling. A5 retention is tested directly through the predeclared drop-A5 governance-measurement variant; the final manuscript does not rely on a companion study for this placement decision.

## 2. Representation
Each archetype r is a joint distribution P_r(p, c, k, h): Beta marginals (parameters below) coupled by a Gaussian copula with correlation matrix
rho(p,c) = +0.5, rho(p,k) = -0.4, rho(c,k) = -0.5, h independent of the others.
Rationale for the copula: within a regime, detection capacity and follow-through capacity co-move (both are functions of administrative capacity), and leakage co-occurs with weak detection and weak follow-through (the Gujarat evidence, ref [12], is a joint observation of weak effective detection with purchased waiver, not two independent facts). h is left independent within-archetype because the between-archetype placement already carries the horizon information and no within-jurisdiction evidence supports a tighter within-archetype coupling. Robustness (predeclared): independence copula; each rho entry shifted by plus or minus 0.2; all marginal means shifted by plus or minus 0.10.

Beta marginals are parameterized by (alpha, beta) fitted to the stated mean with a concentration chosen so the central 90 percent interval matches the placement interval below (calibrated numerically; the intervals, not the Beta family, are the empirical content).

## 3. Placement table (means, 90 percent intervals, Beta parameters)

### A1. High detection, high implementation credibility, low leakage, long horizon (high-enforcement anchor). Anchoring: Singapore; enforcement-pattern evidence from the strongest documented regimes (US majors).
- p: mean 0.80, 90% [0.63, 0.93], Beta(14.4, 3.6)
- c: mean 0.85, 90% [0.71, 0.95], Beta(17, 3)
- k: mean 0.05, 90% [0.01, 0.12], Beta(2, 38)
- h: mean 0.95, 90% [0.90, 0.99], Beta(57, 3)
Basis. Tier 1/2: near-annual inspection coverage of major sources with systematic escalation is documented for the strongest regimes (refs [35] institutional fact, [10]); Singapore's environmental statutes (Environmental Protection and Management Act lineage, NEA and PUB enforcement regimes) have been continuously operative since the 1970s with routinely published enforcement actions. Tier 3: Singapore sits in the top percentiles of WGI rule-of-law, government-effectiveness, and control-of-corruption rankings (ref [26]; used for ordering only). Final submission metadata: WGI ordering support uses the World Bank WGI 2025 Revision, 2024 country observations, accessed August 2026; this memo asserts orderings only and no WGI score is converted into a model probability.

### A2. High formal capacity with elevated leakage in implementation. Anchoring: Malaysia.
- p: mean 0.55, 90% [0.37, 0.72], Beta(11, 9)
- c: mean 0.60, 90% [0.42, 0.77], Beta(12, 8)
- k: mean 0.25, 90% [0.12, 0.40], Beta(6.25, 18.75)
- h: mean 0.85, 90% [0.73, 0.94], Beta(25.5, 4.5)
Basis. Tier 1: a mature statutory apparatus with instrument sophistication unusual for the region and period: the Environmental Quality (Prescribed Premises) (Crude Palm Oil) Regulations 1977 established licence-based discharge standards with an effluent-related licence fee schedule (ref [45]), and the DOE publishes annual Environmental Quality Reports with inspection and prosecution counts. High formal capacity and policy continuity (EQA 1974 continuously in force) support mid-high p and c and high h. Tier 3: WGI places Malaysia clearly above the rest of the anchor set except Singapore on government effectiveness, with control-of-corruption substantially lower than its capacity ranking (ref [26], ordering only); that wedge between capacity and integrity is what k = 0.25 encodes.

### A3. Moderate detection, credible sanction conditional on detection, shorter policy horizon. Anchoring: Philippines and Thailand (contested-politics exemplars).
- p: mean 0.40, 90% [0.23, 0.58], Beta(8, 12)
- c: mean 0.50, 90% [0.32, 0.68], Beta(10, 10)
- k: mean 0.30, 90% [0.16, 0.46], Beta(7.5, 17.5)
- h: mean 0.62, 90% [0.48, 0.76], Beta(18.6, 11.4)
Basis. Tier 1: the Philippine mining episode shows that detection and formal sanction do occur under political salience (audit of 41 mines, 23 closure orders and 5 suspensions, Feb 2017; DAO 2017-10 open-pit ban, Apr 2017) but that implementation and persistence fail: MICC recommended reversal by Oct 2017, EO 130 (Apr 2021) lifted the new-agreements moratorium, DAO 2021-40 (Dec 2021) repealed the ban (refs [46], [47]). Thailand: NCPO Orders 3/2016, 4/2016, 9/2016 under Section 44 altered environmental-assessment sequencing by decree, and coups in 2006 and 2014 bound the credible horizon of any environmental commitment (ref [48]). Transformation for h: documented 2-5 year lifetimes of headline measures imply annualized persistence 0.6-0.8 (CALIBRATION_MEMO Sec 1). Tier 3: WGI/V-Dem ordering places both below A2 on rule of law and control of corruption (refs [26], [27], ordering only).

### A4. Low routine detection, weak implementation, elevated leakage. Anchoring: Indonesia, 1990s-2000s routine-enforcement baseline (pre- and early-PROPER; decentralization era).
- p: mean 0.18, 90% [0.07, 0.32], Beta(3.96, 18.04)
- c: mean 0.30, 90% [0.13, 0.50], Beta(4.5, 10.5)
- k: mean 0.38, 90% [0.22, 0.55], Beta(8.36, 13.64)
- h: mean 0.80, 90% [0.66, 0.91], Beta(20, 5)
Basis. Tier 1: at the June 1995 PROPER pilot, about two thirds of the 187 rated factories were noncompliant (Black plus Red 65 percent) under standards that had formally existed for years, and compliance moved from about 35 to 51 percent only after disclosure (ref [34]); Black-rated firms rarely faced prosecution. Tier 2: Olken's 24 percent missing-expenditures estimate for Indonesian public works (ref [44]) independently anchors the leakage fraction of implemented public processes in the same environment; the Gujarat audit-market evidence (ref [12]) shows how corrupted verification produces the same joint signature (low effective p, high k) and justifies the negative p-k copula entry. h is placed high-moderate: the framework laws persisted even while enforcement was weak, so the horizon of the regime exceeded the credibility of its enforcement. NOTE the deliberate historical anchoring: contemporary Indonesian enforcement (post-PROPER revival, post-OSS reforms) is plausibly stronger; A4 is the weak-routine-enforcement configuration the regional record documents, not a statement about Indonesia today. This is exactly the stylized-configuration-not-country-prediction qualification of Section IV-C.

### A5. Durable single-party regime with moderate detection and credibility. Anchoring: Vietnam.
- p: mean 0.35, 90% [0.20, 0.52], Beta(7.7, 14.3)
- c: mean 0.45, 90% [0.28, 0.62], Beta(9, 11)
- k: mean 0.30, 90% [0.16, 0.46], Beta(6.6, 15.4)
- h: mean 0.93, 90% [0.86, 0.98], Beta(46.5, 3.5)
Basis. Tier 3 primary: V-Dem regime-durability information (single-party continuity since 1976) supports the high-h placement; WGI places Vietnam's regulatory quality and control of corruption in the same band as A3/A4 (refs [26], [27], ordering only). Tier 1 qualitative (AUTHOR TO SOURCE): the 2016 Formosa Ha Tinh marine pollution case, in which a large penalty was assessed and actually collected, indicates c is not negligible when the center engages; MONRE inspection capacity expanded after 2016. This memo flags that episode as knowledge-based context requiring a citable source before submission; it is not cited in the manuscript tables. The load-bearing content of A5 is the h-decoupling, which rests on refs [26], [27].

## 4. G_anchor
G_anchor = union over r of G_r^0.95, where G_r^0.95 is the 95 percent highest-density region of the joint P_r (copula plus Beta marginals). This is a union of joint regions, not a Cartesian product of the marginal spans in Table 2, exactly as eq. (20) requires: for example, the corner (p = 0.9, h = 0.55) lies inside the marginal spans but inside no archetype's joint region, and is therefore labelled simulated extrapolation. The executed analysis computes G_r^0.95 by Monte Carlo (10^5 draws per archetype, kernel HDR at the 95 percent level) and exports the indicator on the sweep grid; the definition is fixed here.

## 5. Empirically anchored vs extrapolated
Empirically anchored: every sweep cell inside G_anchor. Simulated extrapolation: everything else on [0,1]^4, prominently including the near-origin corner (p, c below 0.05), the high-p-low-c wedge beyond A2's region, and all cells with h below 0.45. H2 and H3 are evaluated over G_anchor only; extrapolated cells are computed, reported, and labelled.

## 6. Robustness (predeclared, from Section VIII governance-measurement family)
(i) independence copula; (ii) copula entries plus or minus 0.2; (iii) all marginal means plus or minus 0.10 with intervals re-centered; (iv) A4 re-anchored to a contemporary (stronger) Indonesian configuration, mean p 0.30, c 0.40, k 0.30; (v) drop A5 and recompute G_anchor. Archetype readouts and the H2/H3 verdicts are re-reported under each; features that do not survive are reported as fragile.

## 7. Honest limits
The marginals are placement intervals dressed as Beta distributions for samplability; the empirical content is the interval and the ordering, not the density shape. Nothing here uses a WGI or V-Dem score as a probability. The single most contestable placements are c for A3 (episodic, salience-driven follow-through is hard to average) and the A4 historical anchoring; both are inside the predeclared robustness families, and the readout protocol reports placement sensitivity per archetype as Section V-D requires.
