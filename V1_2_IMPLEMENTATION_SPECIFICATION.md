# V1_2_IMPLEMENTATION_SPECIFICATION.md
Complete handoff for implementing model v1.2. Audit-corrected regeneration, frozen 2026-08-21; supersedes the prior same-named file in full. Companions (all part of the freeze): AMENDMENTS_V1_2.md (incl. the audit-correction addendum), V1_2_CALIBRATION_MEMO.md, V1_2_VALIDATION_MEMO.md, V1_2_PARAMETER_CHANGES.csv, V1_2_LOCKED_PARAMETERS.md, V1_2_CODE_CHECKS.md, v1_2_calibration_script.py. The implementation workspace implements and executes; every scientific decision is already made here.

## 0. Order of operations (binding)
1. git tag v1.2-freeze containing this package unmodified; record timestamp and SHA-256 per file in CHANGELOG. No v1.2 seed before the tag.
2. Gates C1-C5 (V1_2_CODE_CHECKS). Any FAIL: stop, record, return.
3. Implement Sec 1 inside the voluntary module only, respecting V1_2_LOCKED_PARAMETERS (M_of() and tau* byte-identical).
4. Extended verification battery (incl. V13-V17). Any FAIL: return to implementation.
5. Re-execute validation V1-V4 (bit-identical expected; any change is a verification failure).
6. Execute V5-prime: the scored size-gradient test, its 8 placement variants, and the unscored feasibility check, BEFORE inspecting any principal v1.2 voluntary output.
7. Principal v1.2 voluntary arm (primary plane, slices, archetype draws) reusing archived populations and CRN streams; recompute all simultaneous bands (the joint max-t family spans all pairwise surfaces, so S-T bands are recomputed from archived S/T outputs plus new V outputs); regenerate three-instrument maps, H1-H3, readouts, ablations m3-m4 and implicated pairs, V-involving robustness incl. the new family members, V-involving Sobol entries, Figures 2-4.
8. Report per the manuscript impact register; no Results text beyond the predeclared structure.

## 1. The v1.2 voluntary specification (only these changes)
Per firm i, period t = 1..T, governance g, all else as v1.1:
- K_i = 0.02 gamma_i E0_i^2; nu_i = nu_0 E0_i^eta with nu_0 = 0.2460, eta = 1.9520 (normalizer: distribution mean E0 = 1); rho_i = nu_i K_i
- m_i(g) = M_of(g, i), the UNCHANGED shared one-period function (also used by run_standard(); do not modify it)
- Ann(x, n) = (1 - x^n)/(1 - x), Ann(x, 0) = 0; M_i(g, t) = M_of(g, i) x Ann(delta h, T - t), computed ONLY in the voluntary decision path
- theta_1 = sigma(theta_0); theta_{t+1} = sigma(theta_0 + lambda s_t), s_t = max(0, 1 - Ahat_t / A*) (unchanged)
- Decision (expected cost; ties to non-participation): participate and achieve iff p [ nu_i + delta h theta_t (m_i / K_i) Ann(delta h, T - t) ] >= 1; J branches per V1_2_CALIBRATION_MEMO Sec 2
- Realized accounting unchanged: K_i real resource cost; rho_i private, outside C_social; verification draw from the shared monitoring stream; belief term only, no regime switch, no persistent exemption state
- Final period: M = 0; participation iff p nu_i >= 1.

## 2. Frozen numerical values (complete)
tau* = 0.09326697147369385 (byte-equal v1.1 import; NEVER re-solved for use); f_S = 0.1865339429473877; f_T = tau*; nu_0 = 0.2460; eta = 1.9520; theta_0 = -2.2; lambda = 3.0; N = 200; T = 10; delta = 0.95; A*_frac = 0.20; A*_cum = 400; A_min = 20; kappa_M = 0.002; kappa_A = 0.0005; kappa_I = 0.001; populations LN(-0.5, 1.0) x LN(-0.5, 1.0) independent; governance reference (0.50, 0.50, 0.30, 0.75); slices k {0.10, 0.55}, h {0.50, 0.95}; alpha = 0.05; B = 2,000; archetypes, copula, G_anchor per the locked files. Calibration-uncertainty members: top-decile pair (0.2920, 2.0000); uniform-nu (refit rule); tail-cap min(nu_i, 2); one-period-M; collective-preemption; within-tick theta. Constants for unit tests: Ann(0.9025, 8) = 5.742291; Ann(0.9025, 9) = 6.182417.

## 3. Validation execution constants
V1-V4 unchanged. V5-prime: Dutch central (0.75, 0.85, 0.05, 0.95), bands p [0.60, 0.90], c [0.70, 0.95], k [0.02, 0.12], h [0.90, 0.98]; SCORED: Delta_grad (top-vs-bottom E0-quintile abatement intensity gap) 95 percent band strictly above zero, FATAL, one-shot; 8 one-at-a-time placement variants reported (sign flips reportable); UNSCORED: A_V feasibility check reported with the a-priori note (P(0.75 nu_i >= 1) = 0.0855); context items per V1_2_VALIDATION_MEMO Sec 4.

## 4. Reuse and regeneration
Reuse without rerun (conditional on C4 PASS): raw S and T seed-level outputs; S-T difference surfaces as point estimates; archetype masks; S/T-only Sobol outputs; Figure 1. Regenerate: step 7 of Sec 0. Budget: about one third of principal simulation plus V-involving robustness; band recomputation cheap.

## 5. Final audit status
Corrections resolved on the record at this regeneration: tau* restoration; calibration-statistic construct correction with full ladder log; S_top = 0.39 member removal; V5-prime restructure (single scored component); annuity constants. All in AMENDMENTS_V1_2 Sec 7. Superseded v1.0/v1.1 sections per the override map; manuscript sentences per the impact register. Equations (6)-(13) implemented verbatim with the adaptive-theta disclosure; no other contradiction found across the manuscript, CODE_INPUT_SPECIFICATION, CALIBRATION_MEMO, PRE_REGISTRATION_SPECIFICATION, DECISION_MEMO_VOLUNTARY_ARM, and this package.
