# Changelog

## 1.0.0
Initial app release. Converted from the `cnv-backbone-purple-atlas` `cgp-qc-flags` **applet** into
a versioned, namespaced DNAnexus **app** (`org-emee_1`, `aws:eu-central-1`) for the `eggd_atlas_cnv`
somatic CNV workflow.

Conversion changes only: app metadata and an explicit `timeoutPolicy` (1 h). The flag logic and the
`qc_report.tsv` column order are unchanged (the plotter reads columns by name). `sigs_allocation` /
`cup_summary` remain optional and are left unlinked in v0.1 of the workflow (SIGS/CUPpa out of scope).
