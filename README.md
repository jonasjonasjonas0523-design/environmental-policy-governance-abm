# Source-code map

The repository preserves both the archived v1.1 implementation and the final v1.2 implementation because the v1.2 audit explicitly checks reuse/bit-identity of the standard and charge arms.

## Final v1.2 entry points

- `model_v12.py` — final v1.2 core model.
- `battery_v12.py` — extended v1.2 verification battery.
- `validate_v12.py` — v1.2 re-execution / identity checks for V1–V4.
- `v5prime_v12.py` — held-out V5-prime execution.
- `principal_v12.py` — principal five-plane runs and final joint bands.
- `analysis_v12.py` — crossover / high-resolution diagnostics.
- `archetypes_v12.py` — archetype readouts.
- `sobol_v12.py` — global sensitivity.
- `ablations_v12.py` — mechanism ablations.
- `robust_v12.py` — F1/F2/F4 robustness.
- `f3_v12.py` — F3 governance-measurement robustness.
- `figs_v12.py` — computational figure generation from archived outputs.

Supporting implementations for robustness and ablations are in `model_v12_robust.py`, `model_v12_ablations.py`, and `utils_v12.py`.

## Archived v1.1 files

`model.py`, `battery.py`, and `pipeline.py` are the archived v1.1/stage-1 implementation retained for the cross-version audit. They are not the final voluntary-arm specification. Do not use `model.py` as the v1.2 model.

The frozen v1.2 numerical authority is `../docs/V1_2_LOCKED_PARAMETERS.md`; the computational authority is `../docs/FINAL_COMPUTATIONAL_AUDIT.md`.
