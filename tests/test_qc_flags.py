"""Tests for run_qc_flags — QC flag computation logic for eggd_cgp-qc-flags.

conftest.py adds resources/home/ubuntu to sys.path so run_qc_flags can be
imported without a DNAnexus environment.
"""
import csv
import pytest

from run_qc_flags import compute_qc_flags, REQUIRED_FIELDS

# ── Fixtures ──────────────────────────────────────────────────────────────────

PURITY_HDR = '\t'.join([
    'purity', 'ploidy', 'status', 'wholeGenomeDuplication',
    'diploidProportion', 'minPurity', 'maxPurity',
])


def make_purity_tsv(tmp_path, *, purity=0.50, ploidy=2.0, status='NORMAL',
                    wgd='FALSE', dip_prop=0.40, min_pur=0.30, max_pur=0.70):
    p = tmp_path / 'purity.tsv'
    p.write_text(
        f"{PURITY_HDR}\n"
        f"{purity}\t{ploidy}\t{status}\t{wgd}\t{dip_prop}\t{min_pur}\t{max_pur}\n"
    )
    return str(p)


def make_range_tsv(tmp_path, rows=None):
    if rows is None:
        rows = [
            {'purity': '0.12', 'ploidy': '2.1', 'score': '0.05'},
            {'purity': '0.20', 'ploidy': '2.0', 'score': '0.10'},
        ]
    p = tmp_path / 'range.tsv'
    header = list(rows[0].keys())
    lines = ['\t'.join(header)] + ['\t'.join(r[h] for h in header) for r in rows]
    p.write_text('\n'.join(lines) + '\n')
    return str(p)


def make_empty_stub(tmp_path, name):
    p = tmp_path / name
    p.write_text('')
    return str(p)


def run(tmp_path, *, purity_kwargs=None, range_rows=None, sigs_text=None, cuppa_text=None):
    """Helper: build inputs and call compute_qc_flags; return (flags, info, out_tsv_path)."""
    pur  = make_purity_tsv(tmp_path, **(purity_kwargs or {}))
    rng  = make_range_tsv(tmp_path, range_rows)
    sigs = tmp_path / 'sigs.tsv'
    sigs.write_text(sigs_text or '')
    cup  = tmp_path / 'cuppa.tsv'
    cup.write_text(cuppa_text or '')
    out  = tmp_path / 'out.tsv'
    flags, info = compute_qc_flags('S1', pur, str(rng), str(sigs), str(cup), str(out))
    return flags, info, out


# ── PASS ──────────────────────────────────────────────────────────────────────

class TestPass:
    def test_clean_sample_emits_pass(self, tmp_path):
        flags, info, out = run(tmp_path)
        assert flags == []
        rows = list(csv.DictReader(out.open(), delimiter='\t'))
        assert rows[0]['flags'] == 'PASS'
        assert rows[0]['sample_id'] == 'S1'

    def test_output_has_16_columns(self, tmp_path):
        _, _, out = run(tmp_path)
        rows = list(csv.DictReader(out.open(), delimiter='\t'))
        assert len(rows[0]) == 16


# ── Individual flags ──────────────────────────────────────────────────────────

class TestFlags:
    def test_purity_floor_below_threshold(self, tmp_path):
        flags, _, _ = run(tmp_path, purity_kwargs=dict(purity=0.15, min_pur=0.10, max_pur=0.20))
        assert 'PURITY_FLOOR' in flags

    def test_no_purity_floor_at_exactly_0_20(self, tmp_path):
        flags, _, _ = run(tmp_path, purity_kwargs=dict(purity=0.20, min_pur=0.15, max_pur=0.25))
        assert 'PURITY_FLOOR' not in flags

    def test_wgd_suspect_low_purity(self, tmp_path):
        flags, _, _ = run(tmp_path, purity_kwargs=dict(purity=0.30, wgd='TRUE', min_pur=0.25, max_pur=0.35))
        assert 'WGD_SUSPECT' in flags

    def test_no_wgd_suspect_high_purity(self, tmp_path):
        flags, _, _ = run(tmp_path, purity_kwargs=dict(purity=0.80, wgd='TRUE', min_pur=0.75, max_pur=0.85))
        assert 'WGD_SUSPECT' not in flags

    def test_no_wgd_suspect_no_wgd(self, tmp_path):
        flags, _, _ = run(tmp_path, purity_kwargs=dict(purity=0.25, wgd='FALSE', min_pur=0.20, max_pur=0.30))
        assert 'WGD_SUSPECT' not in flags

    def test_flat_genome(self, tmp_path):
        flags, _, _ = run(tmp_path, purity_kwargs=dict(dip_prop=0.90, min_pur=0.40, max_pur=0.60))
        assert 'FLAT_GENOME' in flags

    def test_no_flat_genome_at_threshold(self, tmp_path):
        flags, _, _ = run(tmp_path, purity_kwargs=dict(dip_prop=0.85, min_pur=0.40, max_pur=0.60))
        assert 'FLAT_GENOME' not in flags

    def test_wide_ci(self, tmp_path):
        flags, _, _ = run(tmp_path, purity_kwargs=dict(min_pur=0.10, max_pur=0.55))
        assert 'WIDE_CI' in flags

    def test_no_wide_ci_at_threshold(self, tmp_path):
        flags, _, _ = run(tmp_path, purity_kwargs=dict(min_pur=0.30, max_pur=0.70))
        assert 'WIDE_CI' not in flags  # 0.40 is not > 0.40

    def test_no_tumor_bestfit(self, tmp_path):
        range_rows = [
            {'purity': '0.12', 'ploidy': '2.1', 'score': '0.05'},
            {'purity': '0.20', 'ploidy': '2.0', 'score': '0.10'},
        ]
        flags, info, out = run(
            tmp_path,
            purity_kwargs=dict(status='NO_TUMOR', min_pur=0.05, max_pur=0.20),
            range_rows=range_rows,
        )
        assert 'NO_TUMOR_BESTFIT' in flags
        # lowest-score row is the first one (score 0.05)
        assert info['bestfit_purity'] == '0.120'
        assert info['bestfit_ploidy'] == '2.100'
        rows = list(csv.DictReader(out.open(), delimiter='\t'))
        assert rows[0]['bestfit_purity'] == '0.120'

    def test_low_snv_count(self, tmp_path):
        sigs = 'sig\tallocation\nSBS1\t20\nSBS2\t15\n'
        flags, info, _ = run(tmp_path, sigs_text=sigs)
        assert 'LOW_SNV_COUNT' in flags
        assert info['total_snvs'] == '35'

    def test_no_low_snv_count_above_threshold(self, tmp_path):
        sigs = 'sig\tallocation\nSBS1\t40\nSBS2\t20\n'
        flags, _, _ = run(tmp_path, sigs_text=sigs)
        assert 'LOW_SNV_COUNT' not in flags

    def test_zero_snvs_not_flagged(self, tmp_path):
        """Zero SNVs = absent/empty sigs file, not flagged as LOW_SNV_COUNT."""
        flags, _, _ = run(tmp_path)
        assert 'LOW_SNV_COUNT' not in flags

    def test_multiple_flags_combined(self, tmp_path):
        """Sample below purity floor AND with WGD raises both flags."""
        flags, _, _ = run(
            tmp_path,
            purity_kwargs=dict(purity=0.15, wgd='TRUE', min_pur=0.10, max_pur=0.20),
        )
        assert 'PURITY_FLOOR' in flags
        assert 'WGD_SUSPECT' in flags


# ── Error handling ────────────────────────────────────────────────────────────

class TestErrors:
    def test_empty_purity_tsv_exits(self, tmp_path):
        empty = tmp_path / 'purity.tsv'
        empty.write_text('')
        with pytest.raises(SystemExit):
            compute_qc_flags(
                'S', str(empty), make_range_tsv(tmp_path),
                make_empty_stub(tmp_path, 's.tsv'), make_empty_stub(tmp_path, 'c.tsv'),
                str(tmp_path / 'o.tsv'),
            )

    def test_missing_required_column_exits(self, tmp_path):
        pur = tmp_path / 'purity.tsv'
        # Missing wholeGenomeDuplication, diploidProportion, minPurity, maxPurity
        pur.write_text('purity\tploidy\tstatus\n0.50\t2.0\tNORMAL\n')
        with pytest.raises(SystemExit):
            compute_qc_flags(
                'S', str(pur), make_range_tsv(tmp_path),
                make_empty_stub(tmp_path, 's.tsv'), make_empty_stub(tmp_path, 'c.tsv'),
                str(tmp_path / 'o.tsv'),
            )

    def test_na_purity_exits(self, tmp_path):
        pur = tmp_path / 'purity.tsv'
        pur.write_text(f"{PURITY_HDR}\nNA\t2.0\tNORMAL\tFALSE\t0.40\t0.30\t0.70\n")
        with pytest.raises(SystemExit):
            compute_qc_flags(
                'S', str(pur), make_range_tsv(tmp_path),
                make_empty_stub(tmp_path, 's.tsv'), make_empty_stub(tmp_path, 'c.tsv'),
                str(tmp_path / 'o.tsv'),
            )

    def test_na_ploidy_exits(self, tmp_path):
        pur = tmp_path / 'purity.tsv'
        pur.write_text(f"{PURITY_HDR}\n0.50\tNA\tNORMAL\tFALSE\t0.40\t0.30\t0.70\n")
        with pytest.raises(SystemExit):
            compute_qc_flags(
                'S', str(pur), make_range_tsv(tmp_path),
                make_empty_stub(tmp_path, 's.tsv'), make_empty_stub(tmp_path, 'c.tsv'),
                str(tmp_path / 'o.tsv'),
            )

    def test_no_tumor_empty_range_exits(self, tmp_path):
        pur = make_purity_tsv(tmp_path, status='NO_TUMOR', min_pur=0.05, max_pur=0.20)
        empty_range = tmp_path / 'range.tsv'
        empty_range.write_text('')
        with pytest.raises(SystemExit):
            compute_qc_flags(
                'S', pur, str(empty_range),
                make_empty_stub(tmp_path, 's.tsv'), make_empty_stub(tmp_path, 'c.tsv'),
                str(tmp_path / 'o.tsv'),
            )
