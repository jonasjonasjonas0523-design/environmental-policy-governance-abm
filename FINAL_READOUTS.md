# FINAL READOUTS (v1.2, all planes, final paper-wide joint bands)

Scope: 41x41 p-x-c planes at reference (k=0.30,h=0.75) and slices k in {0.10,0.55}, h in {0.50,0.95}; n=153 principal seeds; paper-wide joint max-t family (crit_A=3.1597, crit_C=3.5874); anchored = G_anchor 95 percent HDR union.

## Per-plane classification (counts of 1681 cells)

| plane | anchored | T_dominant | mixed/no-dominance | unresolved | A-leader T (anchored) | C-leader S (anchored) |
|---|---|---|---|---|---|---|
| reference | 481 | 191 | 0 | 122 | 0.55 | 0.78 |
| k010 | 114 | 169 | 0 | 111 | 1.00 | 0.99 |
| k055 | 131 | 96 | 0 | 152 | 0.00 | 0.30 |
| h050 | 287 | 226 | 0 | 122 | 0.72 | 0.93 |
| h095 | 444 | 89 | 0 | 122 | 0.39 | 0.81 |

Zero cells carry S_dominant or V_dominant labels on any plane. H1 6/6 pair-metric combinations, H2 3/3 pairs, H3 SUPPORTED (within-plane; reference, h050, h095) per H_verdicts_final.json.

## Supported crossover bands (connected components, 4-connectivity)

Component = anchored transition cells adjacent to both supported signs; width = n_cells / max cell extent along the longer axis.

| plane     | pair   | metric   |   component |   n_cells |   p_min |   p_max |   c_min |   c_max |   extent_p_cells |   extent_c_cells |   width_cells |
|:----------|:-------|:---------|------------:|----------:|--------:|--------:|--------:|--------:|-----------------:|-----------------:|--------------:|
| reference | S-T    | C_social |           3 |         2 |   0.55  |   0.55  |   0.375 |   0.4   |                1 |                2 |       1       |
| reference | S-T    | C_social |           4 |         2 |   0.575 |   0.575 |   0.425 |   0.45  |                1 |                2 |       1       |
| reference | S-T    | C_social |           5 |         2 |   0.6   |   0.6   |   0.475 |   0.5   |                1 |                2 |       1       |
| reference | S-T    | C_social |           6 |         2 |   0.625 |   0.625 |   0.55  |   0.575 |                1 |                2 |       1       |
| reference | S-T    | C_social |           7 |         2 |   0.65  |   0.65  |   0.6   |   0.625 |                1 |                2 |       1       |
| reference | S-T    | C_social |           8 |         2 |   0.675 |   0.675 |   0.65  |   0.675 |                1 |                2 |       1       |
| reference | S-T    | C_social |           9 |         4 |   0.7   |   0.725 |   0.7   |   0.75  |                2 |                3 |       1.33333 |
| reference | S-V    | A        |           1 |       156 |   0.05  |   0.75  |   0.325 |   0.65  |               29 |               14 |       5.37931 |
| reference | T-V    | A        |           1 |        88 |   0.05  |   0.65  |   0.325 |   0.525 |               25 |                9 |       3.52    |
| h050      | S-T    | C_social |           1 |         2 |   0.55  |   0.55  |   0.375 |   0.4   |                1 |                2 |       1       |
| h050      | S-T    | C_social |           2 |         2 |   0.575 |   0.575 |   0.425 |   0.45  |                1 |                2 |       1       |
| h050      | S-T    | C_social |           3 |         2 |   0.6   |   0.6   |   0.475 |   0.5   |                1 |                2 |       1       |
| h050      | S-T    | C_social |           4 |         2 |   0.625 |   0.625 |   0.55  |   0.575 |                1 |                2 |       1       |
| h050      | S-T    | C_social |           5 |         2 |   0.65  |   0.65  |   0.6   |   0.625 |                1 |                2 |       1       |
| h050      | S-V    | A        |           1 |        99 |   0.175 |   0.65  |   0.425 |   0.55  |               20 |                6 |       4.95    |
| h050      | T-V    | A        |           1 |        45 |   0.175 |   0.5   |   0.325 |   0.5   |               14 |                8 |       3.21429 |
| h095      | S-T    | C_social |           3 |         1 |   0.525 |   0.525 |   0.35  |   0.35  |                1 |                1 |       1       |
| h095      | S-T    | C_social |           4 |         2 |   0.55  |   0.55  |   0.375 |   0.4   |                1 |                2 |       1       |
| h095      | S-T    | C_social |           5 |         2 |   0.575 |   0.575 |   0.425 |   0.45  |                1 |                2 |       1       |
| h095      | S-T    | C_social |           6 |         2 |   0.6   |   0.6   |   0.475 |   0.5   |                1 |                2 |       1       |
| h095      | S-T    | C_social |           7 |         2 |   0.625 |   0.625 |   0.55  |   0.575 |                1 |                2 |       1       |
| h095      | S-T    | C_social |           8 |         2 |   0.65  |   0.65  |   0.6   |   0.625 |                1 |                2 |       1       |
| h095      | S-T    | C_social |           9 |         2 |   0.675 |   0.675 |   0.65  |   0.675 |                1 |                2 |       1       |
| h095      | S-T    | C_social |          10 |         4 |   0.7   |   0.725 |   0.7   |   0.75  |                2 |                3 |       1.33333 |
| h095      | S-V    | A        |           2 |        68 |   0.125 |   0.425 |   0.4   |   0.725 |               13 |               14 |       4.85714 |
| h095      | T-V    | A        |           1 |        92 |   0.1   |   0.75  |   0.375 |   0.7   |               27 |               14 |       3.40741 |

## A_min sensitivity (report-only; PARAMETER_PROVENANCE row)

CE-defined cells (mean A >= A_min) per plane x instrument at A_min in {8, 20, 40} (= {0.02, 0.05, 0.10} x A*_cum). Classification never uses CE; identical labels at all three thresholds by construction.

| plane     |   ('A_min_20_defined', 'S') |   ('A_min_20_defined', 'T') |   ('A_min_20_defined', 'V') |   ('A_min_40_defined', 'S') |   ('A_min_40_defined', 'T') |   ('A_min_40_defined', 'V') |   ('A_min_8_defined', 'S') |   ('A_min_8_defined', 'T') |   ('A_min_8_defined', 'V') |
|:----------|----------------------------:|----------------------------:|----------------------------:|----------------------------:|----------------------------:|----------------------------:|---------------------------:|---------------------------:|---------------------------:|
| h050      |                        1426 |                        1410 |                        1558 |                        1263 |                        1255 |                        1435 |                       1542 |                       1542 |                       1599 |
| h095      |                        1426 |                        1410 |                        1558 |                        1263 |                        1255 |                        1435 |                       1542 |                       1542 |                       1599 |
| k010      |                        1460 |                        1464 |                        1558 |                        1326 |                        1326 |                        1435 |                       1559 |                       1559 |                       1599 |
| k055      |                        1338 |                        1326 |                        1558 |                        1120 |                        1099 |                        1435 |                       1499 |                       1487 |                       1599 |
| reference |                        1426 |                        1410 |                        1558 |                        1263 |                        1255 |                        1435 |                       1542 |                       1542 |                       1599 |

## Higher-resolution S/T boundary check (81x81, post-principal robustness only)

Max nearest-matched crossing shift 41->81: 0.7338 (tol 0.00625, within: False); A: shift 0.0012500000000000011, matched 43, mismatched 0; C_social: shift 0.73375, matched 51, mismatched 0. POST-PRINCIPAL ROBUSTNESS ONLY: the preregistered 41x41 grid was a fallback to the finest candidate and did NOT formally converge; this S/T-only 81x81 check does not replace the 41x41 result and no formal convergence is claimed.

## Diagnostics and robustness (pointers)

- Ablations m1-m5 + 6 implicated pairs: results_v12/ablations/ (identity gate PASS; m1 evasion moves only T-channels, m3 rho=0 only V-channels with a vanishing V-pair class, m2 collection friction moves all three).
- Robustness F1/F2/F4 (49 executed members incl. kappa post-processing and the m1 reuse): results_v12/robustness/ROBUSTNESS_REPORT.md; strict Sec-8 core nearly empty on abatement orderings (F2 leverage corners + F1 behavioral flips responsible), F4 population core large; uniform-nu and nu-visibility recorded infeasible-to-calibrate under the frozen single-moment rule (nu_refits.json).
- F3 governance measurement (9 variants): results_v12/f3/f3_report.md; H2 3/3 and H3 SUPPORTED under every variant.
- Sobol, archetype readouts, V5-prime: respective results_v12 subfolders.

Added diagnostic (post hoc, labeled): at the 1681 shared grid points the pointwise 2SE-supported Delta_ST signs agree between 41x41 and 81x81 at 1.0000 (A) and 1.0000 (C_social), with 0 and 0 sign contradictions respectively; the large crossing-metric shift therefore reflects crossing-set reconfiguration on near-flat lines (pairings appearing/disappearing against the 2SE rule), not displacement of the supported-sign regions that the classification uses.
