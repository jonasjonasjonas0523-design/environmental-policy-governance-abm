# F3 governance-measurement robustness (memo Sec 6 set + stress pair)

Anchors and archetype readouts recomputed per variant on the EXISTING final joint bands and the frozen seed-2026 draw stream (CRN across variants); H2/H3 re-verdicts per the memo ("readouts and the H2/H3 verdicts are re-reported under each").

| variant | H2 S-T/S-V/T-V | H3 | anchored frac (ref plane) |
|---|---|---|---|
| independence | 1/1/1 | SUPPORTED | 0.405 |
| cop_p02 | 1/1/1 | SUPPORTED | 0.240 |
| cop_m02 | 1/1/1 | SUPPORTED | 0.206 |
| means_p010 | 1/1/1 | SUPPORTED | 0.193 |
| means_m010 | 1/1/1 | SUPPORTED | 0.181 |
| cont_A4 | 1/1/1 | SUPPORTED | 0.271 |
| drop_A5 | 1/1/1 | SUPPORTED | 0.286 |
| stress_p20 | 1/1/1 | SUPPORTED | 0.208 |
| stress_m20 | 1/1/1 | SUPPORTED | 0.238 |

## Modal-leader stability (A metric)

| archetype   | baseline   | cont_A4   | cop_m02   | cop_p02   | drop_A5   | independence   | means_m010   | means_p010   | stress_m20   | stress_p20   |
|:------------|:-----------|:----------|:----------|:----------|:----------|:---------------|:-------------|:-------------|:-------------|:-------------|
| A1          | T          | T         | T         | T         | T         | T              | T            | T            | T            | T            |
| A2          | T          | T         | T         | T         | T         | T              | T            | T            | T            | T            |
| A3          | T          | T         | T         | T         | T         | T              | T            | T            | T            | T            |
| A4          | V          | T         | V         | V         | V         | V              | V            | V            | V            | V            |
| A5          | T          | T         | T         | T         | nan       | T              | V            | V            | V            | V            |

## Modal-leader stability (C_social)

| archetype   | baseline   | cont_A4   | cop_m02   | cop_p02   | drop_A5   | independence   | means_m010   | means_p010   | stress_m20   | stress_p20   |
|:------------|:-----------|:----------|:----------|:----------|:----------|:---------------|:-------------|:-------------|:-------------|:-------------|
| A1          | S          | S         | S         | S         | S         | S              | S            | T            | S            | S            |
| A2          | S          | S         | S         | S         | S         | S              | S            | S            | S            | S            |
| A3          | S          | S         | S         | S         | S         | S              | S            | S            | S            | S            |
| A4          | S          | S         | S         | S         | S         | S              | S            | S            | S            | S            |
| A5          | S          | S         | S         | S         | nan       | S              | S            | S            | S            | S            |
