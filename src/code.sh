#!/bin/bash
# eggd_cgp-qc-flags v1.0.0 (converted from cgp-qc-flags applet: metadata + timeoutPolicy; logic unchanged)
# Post-processes PURPLE output to emit structured QC flags.
# Sources of flags (from REFERENCE.md §11):
#   PURITY_FLOOR      purity < 0.20 (near noise floor of backbone panel)
#   WGD_SUSPECT       WGD=true AND purity < 0.35 (unreliable below ~35%)
#   FLAT_GENOME       diploidProportion > 0.85 (weak CN signal, purity anchor uncertain)
#   WIDE_CI           maxPurity - minPurity > 0.40 (wide confidence interval)
#   NO_TUMOR_BESTFIT  status=NO_TUMOR → surfaces best-fit purity from range TSV
#   LOW_SNV_COUNT     total allocated SNVs < 50 (SIGS/CUPpa SNV features unreliable)
set -euo pipefail

main() {
    case "${sample_id}" in
        *[!A-Za-z0-9._-]* | "" | .* | -* )
            echo "ERROR: unsafe sample_id '${sample_id}' (allowed: A-Za-z0-9._-, no leading '-'/'.')" >&2; exit 1 ;;
    esac
    echo "====================================================="
    echo " eggd_cgp-qc-flags: QC flag generation"
    echo " Sample  : ${sample_id}"
    echo "====================================================="

    echo "[1/4] Downloading inputs..."
    dx download "${purity_tsv}"       -o purity.tsv
    dx download "${purity_range_tsv}" -o purity_range.tsv

    # Explicit branches: absent optional inputs are allowed (empty stub), but a FAILED
    # download of a supplied input must abort under set -e (not fall through to touch).
    if [[ -n "${sigs_allocation:-}" ]]; then dx download "${sigs_allocation}" -o sigs.tsv; else touch sigs.tsv; fi
    if [[ -n "${cup_summary:-}" ]]; then dx download "${cup_summary}" -o cuppa.tsv; else touch cuppa.tsv; fi

    echo "[2/4] Computing QC flags..."
    python3 /home/ubuntu/run_qc_flags.py "${sample_id}" purity.tsv purity_range.tsv sigs.tsv cuppa.tsv

    echo "[3/4] Uploading QC report..."
    qc_report=$(dx upload "${sample_id}.qc_report.tsv" --brief)
    dx-jobutil-add-output qc_report "${qc_report}" --class=file

    echo "====================================================="
    echo " eggd_cgp-qc-flags DONE: ${sample_id}"
    echo "====================================================="
}
