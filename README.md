# eggd_cgp-qc-flags

Post-processes PURPLE output into a structured QC report for the
[`eggd_atlas_cnv`](https://github.com/eastgenomics/eggd_atlas_cnv) somatic CNV workflow.

## What it does
Reads `*.purple.purity.tsv` (+ `.purple.purity.range.tsv` for `NO_TUMOR`) and emits
`{sample_id}.qc_report.tsv` with flags: `PURITY_FLOOR`, `WGD_SUSPECT`, `FLAT_GENOME`, `WIDE_CI`,
`NO_TUMOR_BESTFIT`, `LOW_SNV_COUNT`.

## Inputs
`sample_id`, `purity_tsv`, `purity_range_tsv` (+ optional `sigs_allocation`, `cup_summary` — unlinked in v0.1).

## Outputs
`qc_report` (`{sample_id}.qc_report.tsv`, 16 columns) → purple_plotter.

## Notes
- Pure-Python post-processing; no `execDepends` beyond the system `python3`.
- Instance `mem1_ssd1_v2_x2`; timeout 1 h. Column order is fixed (plotter reads by name).
