# V1_2_LOCKED_PARAMETERS.md
Locked-from-v1.1 register, audit-corrected regeneration. Frozen 2026-08-21. Everything below is byte-identical between the executed v1.1 and v1.2; deviation is a blocker. The ONLY changes are the two items at the bottom plus the declared bookkeeping in AMENDMENTS_V1_2 Sec 3.

| Item | Status | Locked value / form |
|---|---|---|
| tau* | LOCKED, byte-equal import from the v1.1 configuration | 0.09326697147369385 (v1.1 solve: seed 424242, n_quad = 10^6). Never re-solved for use. f_S = 0.1865339429473877; f_T = tau* |
| M_of(): the one-period minimized standard cost function | LOCKED, byte-identical (run_standard() depends on it) | as v1.1; the v1.2 annuity multiplication wraps its RETURN VALUE inside the voluntary decision path only |
| Standard instrument: eq. (7), Ebar_i = 0.8 E0_i, closed-form optimum, run_standard() | LOCKED | as v1.1 |
| Charge instrument: eq. (8), corner reporting solution, tau* solve method (for verification only) | LOCKED | as v1.1 |
| p, c, k, h definitions; Table 2 spans; reference point (0.50, 0.50, 0.30, 0.75); slices k {0.10, 0.55}, h {0.50, 0.95} | LOCKED | as v1.1 |
| Population distributions LN(-0.5, 1.0) x LN(-0.5, 1.0), independence; N = 200 | LOCKED | as v1.1 |
| T = 10; delta = 0.95; A* = 0.20/period; A_min = 20; allocation rules | LOCKED | as v1.1 |
| kappa_M = 0.002; kappa_A = 0.0005; kappa_I = 0.001; C_social and T_fiscal accounting | LOCKED | as v1.1 |
| theta_0 = -2.2; lambda = 3.0; logistic form; shortfall eq. (11); recognition via shared verification | LOCKED | as v1.1 |
| A1 structure (verification-gated exemption; branch comparison; ties to non-participation); adaptive theta timing (C2) | LOCKED | as v1.1 |
| Archetypes A1-A5; copula; G_anchor | LOCKED | as v1.1 |
| alpha = 0.05; B = 2,000 max-t; epsilon rules; replicate-count rule; classification tree; grid rule; CRN protocol; seed ranges | LOCKED | as v1.1 |
| Validation V1-V4 | LOCKED | as v1.1; re-executed under v1.2 as a gate, expected bit-identical |
| Verification battery v1.1 items; Sobol groups; ablation list m1-m5; outcome metrics; paired estimand | LOCKED | as v1.1 (additive items appended; nu-coordinates read the v1.2 construction) |

CHANGED (complete list):
1. M_i in the voluntary path: continuation value M_of(g, i) x Ann(delta h, T - t) per Amendment A2.
2. rho_i construction: nu_i = 0.2460 x E0_i^1.9520.
Declared bookkeeping: verification additions; robustness additions; Table 4 restructure (v1.1 V5 record + V5-prime single scored component); provenance rows superseded by V1_2_PARAMETER_CHANGES.csv.
Flag rule: any apparent need to touch a locked row stops work and returns for a decision. Nothing changes silently.
