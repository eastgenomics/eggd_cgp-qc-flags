#!/usr/bin/env python3
"""eggd_cgp-qc-flags: compute PURPLE QC flags and emit a structured report TSV.

Called by code.sh as:
    python3 /home/ubuntu/run_qc_flags.py <sample_id> <purity_tsv> <purity_range_tsv> <sigs_tsv> <cuppa_tsv>
The output file is written to <sample_id>.qc_report.tsv in the current directory.
"""
import sys
import csv

REQUIRED_FIELDS = [
    'purity', 'ploidy', 'status', 'wholeGenomeDuplication',
    'diploidProportion', 'minPurity', 'maxPurity',
]


def compute_qc_flags(sample_id, purity_file, range_file, sigs_file, cuppa_file, out_file):
    """Compute QC flags from PURPLE output files and write a TSV report.

    Parameters
    ----------
    sample_id   : str   — sample identifier
    purity_file : str   — path to *.purple.purity.tsv
    range_file  : str   — path to *.purple.purity.range.tsv
    sigs_file   : str   — path to *.sig.allocation.tsv (may be empty stub)
    cuppa_file  : str   — path to *.cuppa.pred_summ.tsv (may be empty stub)
    out_file    : str   — path to write the output *.qc_report.tsv

    Returns
    -------
    (flags, info) : (list[str], dict[str, str])
        flags — list of raised flag names (empty = PASS)
        info  — dict of string-formatted values written to the TSV

    Exits via sys.exit(1) on fatal errors: empty purity TSV, missing required
    columns, non-numeric purity/ploidy, or NO_TUMOR best-fit recovery failure.
    """
    flags = []
    info  = {}

    # ── 1. Read purity TSV ────────────────────────────────────────────────
    with open(purity_file) as f:
        reader = csv.DictReader(f, delimiter='\t')
        rows = list(reader)

    if not rows:
        print("ERROR: purity.tsv has no data row", file=sys.stderr)
        sys.exit(1)

    missing = [fld for fld in REQUIRED_FIELDS if fld not in rows[0]]
    if missing:
        print(f"ERROR: purity.tsv is missing required columns: {missing}", file=sys.stderr)
        sys.exit(1)

    d = rows[0]
    try:
        purity = float(d['purity'])
        ploidy = float(d['ploidy'])
    except ValueError as e:
        print(
            f"ERROR: purity/ploidy is not a valid number in purity.tsv "
            f"(purity={d['purity']!r}, ploidy={d['ploidy']!r}): {e}",
            file=sys.stderr,
        )
        sys.exit(1)

    status   = d['status']
    wgd      = d['wholeGenomeDuplication'].upper() == 'TRUE'
    dip_prop = float(d.get('diploidProportion', 0))
    min_pur  = float(d.get('minPurity', purity))
    max_pur  = float(d.get('maxPurity', purity))

    info['purity']   = f"{purity:.3f}"
    info['ploidy']   = f"{ploidy:.3f}"
    info['status']   = status
    info['wgd']      = 'TRUE' if wgd else 'FALSE'
    info['dip_prop'] = f"{dip_prop:.3f}"
    info['ci_width'] = f"{max_pur - min_pur:.3f}"

    if purity < 0.20:
        flags.append('PURITY_FLOOR')
    if wgd and purity < 0.35:
        flags.append('WGD_SUSPECT')
    if dip_prop > 0.85:
        flags.append('FLAT_GENOME')
    if (max_pur - min_pur) > 0.40:
        flags.append('WIDE_CI')

    # ── 2. NO_TUMOR best-fit from range TSV ──────────────────────────────
    bestfit_purity = None
    bestfit_ploidy = None
    if status == 'NO_TUMOR':
        flags.append('NO_TUMOR_BESTFIT')
        try:
            with open(range_file) as f:
                rdr = csv.DictReader(f, delimiter='\t')
                best_row = min(rdr, key=lambda r: float(r.get('score', '9999')))
                bestfit_purity = float(best_row['purity'])
                bestfit_ploidy = float(best_row['ploidy'])
                info['bestfit_purity'] = f"{bestfit_purity:.3f}"
                info['bestfit_ploidy'] = f"{bestfit_ploidy:.3f}"
                info['bestfit_score']  = best_row.get('score', 'N/A')
        except Exception as e:
            print(f"ERROR: NO_TUMOR best-fit recovery failed from {range_file}: {e}", file=sys.stderr)
            sys.exit(1)

    # ── 3. SIGS — check SNV count ─────────────────────────────────────────
    total_snvs = 0
    try:
        with open(sigs_file) as f:
            rdr = csv.DictReader(f, delimiter='\t')
            for row in rdr:
                try:
                    total_snvs += float(row.get('allocation', 0))
                except Exception:
                    pass
        info['total_snvs'] = f"{total_snvs:.0f}"
        if 0 < total_snvs < 50:
            flags.append('LOW_SNV_COUNT')
    except Exception:
        info['total_snvs'] = 'N/A'

    # ── 4. CUPpa — top prediction ─────────────────────────────────────────
    try:
        with open(cuppa_file) as f:
            rdr = csv.DictReader(f, delimiter='\t')
            rows_c = list(rdr)
        combined = [r for r in rows_c if r.get('clf_name') == 'dna_combined']
        if combined:
            r = combined[0]
            info['cuppa_top1']  = r.get('pred_class_1', 'N/A')
            info['cuppa_prob1'] = f"{float(r.get('pred_prob_1', 0)):.3f}"
            info['cuppa_top2']  = r.get('pred_class_2', 'N/A')
            info['cuppa_prob2'] = f"{float(r.get('pred_prob_2', 0)):.3f}"
    except Exception:
        info['cuppa_top1'] = 'N/A'

    # ── 5. Write report ───────────────────────────────────────────────────
    flag_str = ';'.join(flags) if flags else 'PASS'
    print(f"\n{'='*55}")
    print(f" QC Report: {sample_id}")
    print(f"{'='*55}")
    print(f"  Flags     : {flag_str}")
    print(f"  Purity    : {info['purity']} [{min_pur:.2f},{max_pur:.2f}]")
    print(f"  Ploidy    : {info['ploidy']}")
    print(f"  Status    : {info['status']}")
    print(f"  WGD       : {info['wgd']}")
    print(f"  DipProp   : {info['dip_prop']}")
    if bestfit_purity is not None:
        print(f"  BestFit   : purity={info['bestfit_purity']} ploidy={info['bestfit_ploidy']}")
    print(f"  SNVs      : {info.get('total_snvs', 'N/A')}")
    print(f"  CUPpa     : {info.get('cuppa_top1', 'N/A')} (p={info.get('cuppa_prob1', 'N/A')})")
    print(f"{'='*55}\n")

    header_fields = [
        'sample_id', 'flags', 'purity', 'ploidy', 'status', 'wgd', 'diploidProportion',
        'ci_width', 'bestfit_purity', 'bestfit_ploidy', 'bestfit_score',
        'total_snvs', 'cuppa_top1', 'cuppa_prob1', 'cuppa_top2', 'cuppa_prob2',
    ]
    with open(out_file, 'w') as f:
        f.write('\t'.join(header_fields) + '\n')
        out_row = [
            sample_id, flag_str,
            info['purity'], info['ploidy'], info['status'], info['wgd'],
            info['dip_prop'], info['ci_width'],
            info.get('bestfit_purity', ''), info.get('bestfit_ploidy', ''), info.get('bestfit_score', ''),
            info.get('total_snvs', ''), info.get('cuppa_top1', ''), info.get('cuppa_prob1', ''),
            info.get('cuppa_top2', ''), info.get('cuppa_prob2', ''),
        ]
        f.write('\t'.join(out_row) + '\n')

    print(f"Written: {out_file}")
    return flags, info


if __name__ == '__main__':
    if len(sys.argv) != 6:
        print(
            f"Usage: {sys.argv[0]} <sample_id> <purity_tsv> <purity_range_tsv> <sigs_tsv> <cuppa_tsv>",
            file=sys.stderr,
        )
        sys.exit(1)
    _sample_id   = sys.argv[1]
    _purity_file = sys.argv[2]
    _range_file  = sys.argv[3]
    _sigs_file   = sys.argv[4]
    _cuppa_file  = sys.argv[5]
    _out_file    = f"{_sample_id}.qc_report.tsv"
    compute_qc_flags(_sample_id, _purity_file, _range_file, _sigs_file, _cuppa_file, _out_file)
