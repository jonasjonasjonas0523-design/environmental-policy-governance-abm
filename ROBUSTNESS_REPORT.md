# Robustness families F1/F2/F4: robust core on the reference plane

Variants evaluated: 50 (+9 F3 measurement variants reported separately in f3_report.md). Baseline = staged reference-plane classification. Robust core = supported sign preserved across baseline and ALL variants (PRE_REG Sec 8 intersection rule).

- S-T | A: baseline supported 1542 cells -> robust core 0 (0% retained); top breakers: {'satisficing': 1262, 'rep_standard': 1262, 'thr_cost': 1260, 'mS4': 1099, 'mT05': 613, 'power15': 536}
- S-T | C_social: baseline supported 1418 cells -> robust core 7 (0% retained); top breakers: {'mS4': 801, 'rep_standard': 801, 'thr_cost': 801, 'normative': 590, 'mT05': 552, 'power3': 510}
- S-V | A: baseline supported 1391 cells -> robust core 0 (0% retained); top breakers: {'rep_standard': 759, 'power3': 694, 'sE05': 684, 'cara_r10': 659, 'power15': 632, 'mS1': 632}
- S-V | C_social: baseline supported 1640 cells -> robust core 41 (2% retained); top breakers: {'rep_standard': 1599, 'tgt_cost': 1599, 'tgt_uniform': 1576, 'cara_r10': 1316, 'sE05': 1025, 'tail_cap': 712}
- T-V | A: baseline supported 1477 cells -> robust core 38 (3% retained); top breakers: {'power15': 817, 'sE15': 797, 'satisficing': 765, 'sg05': 706, 'power3': 600, 'cara_r10': 562}
- T-V | C_social: baseline supported 1640 cells -> robust core 41 (2% retained); top breakers: {'tgt_cost': 1599, 'tgt_uniform': 1395, 'cara_r10': 1278, 'sE05': 990, 'tail_cap': 713, 'cara_r1': 521}
- label map: baseline classified 1560 -> identical under all variants 118

| variant        | family   |   crit_A |   crit_C |   sign_agreement |   supported_cells |
|:---------------|:---------|---------:|---------:|-----------------:|------------------:|
| cara_r1        | F1       |  3.32285 |  3.52522 |         0.752627 |              9577 |
| cara_r10       | F1       |  3.44768 |  3.49967 |         0.491771 |              9781 |
| convex_evasion | F1       |  3.09991 |  3.28931 |         0.947254 |              9467 |
| evasion_off    | F1       |  3.09991 |  3.27822 |         0.944676 |              9467 |
| normative      | F1       |  2.9183  |  2.77917 |         0.751834 |              9799 |
| power15        | F1       |  3.14627 |  3.61007 |         0.717827 |              9087 |
| power3         | F1       |  2.96169 |  3.33692 |         0.730914 |              9702 |
| pwl_cost       | F1       |  3.49299 |  3.6228  |         0.831846 |              9209 |
| satisficing    | F1       |  3.70464 |  3.80664 |         0.681836 |              8285 |
| Astar010       | F2       |  3.10636 |  3.4758  |         0.977097 |              9132 |
| Astar030       | F2       |  3.11005 |  3.4758  |         0.976899 |              9080 |
| T20            | F2       |  3.11029 |  3.42939 |         0.971743 |              9113 |
| T5             | F2       |  3.10636 |  3.50484 |         0.97204  |              9137 |
| collective     | F2       |  3.1114  |  3.4758  |         0.963415 |              9154 |
| lam15          | F2       |  3.10639 |  3.4758  |         0.972834 |              9131 |
| lam45          | F2       |  3.11581 |  3.4758  |         0.945866 |              9066 |
| linear_theta   | F2       |  3.11012 |  3.4758  |         0.988796 |              9101 |
| mS1            | F2       |  2.96929 |  3.41087 |         0.816082 |              9393 |
| mS4            | F2       |  3.00287 |  3.33059 |         0.720603 |              9412 |
| mT05           | F2       |  2.98202 |  3.50006 |         0.845429 |              9163 |
| mT2            | F2       |  3.0729  |  3.42897 |         0.9651   |              9311 |
| one_period_M   | F2       |  3.10945 |  3.4758  |         0.970652 |              9127 |
| probit_theta   | F2       |  3.11001 |  3.4758  |         0.997422 |              9110 |
| rep_standard   | F2       |  3.89967 |  4.31034 |         0.501091 |              9582 |
| tail_cap       | F2       |  3.23661 |  3.60805 |         0.767896 |              9379 |
| tgt_cost       | F2       |  3.30042 |  3.51402 |         0.525977 |              9515 |
| tgt_uniform    | F2       |  3.2928  |  3.50487 |         0.547591 |              9051 |
| theta0_hi      | F2       |  3.11214 |  3.4758  |         0.969463 |              9082 |
| theta0_lo      | F2       |  3.10636 |  3.4758  |         0.976403 |              9124 |
| thr_cost       | F2       |  3.11288 |  3.70744 |         0.731112 |              7851 |
| thr_uniform    | F2       |  2.91838 |  3.43284 |         0.810331 |              9508 |
| top_decile     | F2       |  3.10522 |  3.47331 |         0.928614 |              9150 |
| within_tick    | F2       |  3.11018 |  3.4758  |         0.995439 |              9100 |
| N1000          | F4       |  3.16244 |  3.48764 |         0.960738 |              9504 |
| N50            | F4       |  3.09012 |  3.509   |         0.904025 |              8310 |
| corr03         | F4       |  3.10446 |  3.54312 |         0.941602 |              9090 |
| delta090       | F4       |  3.10636 |  3.4758  |         0.990284 |              9114 |
| delta097       | F4       |  3.1101  |  3.4758  |         0.994646 |              9110 |
| kappaA_d5      | F4       |  3.11006 |  3.4758  |         1        |              9108 |
| kappaA_x5      | F4       |  3.11006 |  3.4758  |         1        |              9108 |
| kappaI_d5      | F4       |  3.11006 |  3.4135  |         0.983938 |              9159 |
| kappaI_x5      | F4       |  3.11006 |  3.5693  |         0.9767   |              9079 |
| kappaM_d5      | F4       |  3.11006 |  3.4758  |         1        |              9108 |
| kappaM_x5      | F4       |  3.11006 |  3.4758  |         1        |              9108 |
| kappa_d5       | F4       |  3.11006 |  3.4135  |         0.983938 |              9159 |
| kappa_x5       | F4       |  3.11006 |  3.5693  |         0.9767   |              9079 |
| sE05           | F4       |  3.26596 |  3.6062  |         0.58606  |              8769 |
| sE15           | F4       |  3.01771 |  3.48926 |         0.77563  |              9536 |
| sg05           | F4       |  3.10188 |  3.61492 |         0.741523 |              8738 |
| sg15           | F4       |  3.07991 |  3.34654 |         0.78951  |              9511 |

## Diagnostic decomposition (NOT preregistered; interpretive aid to the strict core)
Strict core = sign supported identically in baseline and every variant (the Sec 8 rule, reported above).
Per-family cores isolate the responsible families; the contradiction-free set keeps baseline-supported
cells never OPPOSITELY supported by any variant (pure de-support from precision loss allowed).

| pair\|metric | baseline | strict core | F1-only | F2-only | F4-only | contradiction-free |
|---|---|---|---|---|---|---|
| S-T\|A | 1542 | 0 | 0 | 0 | 781 | 6 |
| S-T\|C_social | 1418 | 7 | 395 | 165 | 644 | 21 |
| S-V\|A | 1391 | 0 | 46 | 0 | 75 | 29 |
| S-V\|C_social | 1640 | 41 | 324 | 41 | 615 | 41 |
| T-V\|A | 1477 | 38 | 44 | 704 | 135 | 54 |
| T-V\|C_social | 1640 | 41 | 362 | 41 | 650 | 41 |

Reading: the S-T|A and S-V|A strict cores are emptied primarily by the F2 leverage corners
(m_S in {1,4} x tau*, m_T in {0.5,2} x tau*) which re-lever the instrument orderings by design, and by
structural behavioral variants (satisficing, CARA r=10, rep_standard). Where signs survive within-family
(F4 population core), fragility is attributable to instrument-design leverage, not to population or
statistical artifacts. The v1.2 voluntary-construction members (one_period_M, top_decile, tail_cap,
collective, within_tick) are individually reported in robustness_summary.csv.

EXCLUSIONS ON THE RECORD: uniform-nu (eta=0) and nu-visibility (eta=1) are infeasible to calibrate
under the frozen single-moment rule at the adopted rung (EVER_all(nu0=0)=0.4131 > 0.15 target,
monotone in nu0; mechanism = V1_2_CALIBRATION_MEMO Sec 4 ignition; see nu_refits.json).
