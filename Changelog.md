# Changelog

## 1.0.0
Initial app release. Converted from the `cnv-backbone-purple-atlas` `cgp-qc-flags` **applet** into
a versioned, namespaced DNAnexus **app** (`org-emee_1`, `aws:eu-central-1`) for the `eggd_atlas_cnv`
somatic CNV workflow.

Conversion changes only: app metadata and an explicit `timeoutPolicy` (1 h). The flag logic and the
`qc_report.tsv` column order are unchanged (the plotter reads columns by name). `sigs_allocation` /
`cup_summary` remain optional and are left unlinked in v0.1 of the workflow (SIGS/CUPpa out of scope).

### Post-release fixes (PR #1 review)
- README: expanded Inputs section to bullet points listing source app for each input.
- `dxapp.json`: added `help` field to every inputSpec entry.
- `src/code.sh` / `resources/home/ubuntu/run_qc_flags.py`: extracted the Python QC logic from
  an inline heredoc into a bundled module (`resources/home/ubuntu/run_qc_flags.py`) to enable
  unit testing; added validation of required purity TSV columns (`REQUIRED_FIELDS` check); added
  `try/except ValueError` guard on purity/ploidy float conversion.
- `tests/`: added pytest suite (23 tests) covering PASS output, all six flag thresholds,
  `NO_TUMOR_BESTFIT` best-fit selection, `LOW_SNV_COUNT`, and error-path exits (empty TSV,
  missing columns, NA values, empty range TSV for NO_TUMOR).
- CI: added `pytest tests/ -v` step.
