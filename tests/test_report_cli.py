"""
End-to-end tests for `tirmite pair --report`.

The point of these is that adding a report must not disturb anything else the
run writes, and that a failure while building one must not fail the run.
"""

import argparse
from pathlib import Path
import random

import pytest

from tirmite.cli.hmm_pair import main, validate_arguments

MODEL_LEN = 60
LEFT_MOTIF = 'CAGGGGTGCAAAAGTACAATTTGGCCCATTAGCTGACCAGTTGACCATGGACCTTAGCAT'
RIGHT_MOTIF = 'ATGCTAAGGTCCATGGTCAACTGGTCAGCTAATGGGCCAAATTGTACTTTTGCACCCCTG'

CONTIG_LEN = 40_000
ELEMENTS = [(2_000, 5_000), (12_000, 14_500), (25_000, 31_000)]


def revcomp(seq):
    return seq.translate(str.maketrans('ACGT', 'TGCA'))[::-1]


@pytest.fixture
def run_inputs(tmp_path):
    """A genome with three implanted elements, plus a matching nhmmer table."""
    rng = random.Random(11)
    seq = [rng.choice('ACGT') for _ in range(CONTIG_LEN)]

    rows = []
    for start, end in ELEMENTS:
        seq[start - 1 : start - 1 + MODEL_LEN] = list(LEFT_MOTIF)
        right_start = end - MODEL_LEN + 1
        seq[right_start - 1 : end] = list(revcomp(RIGHT_MOTIF))
        rows.append(('TIR_L', start, start + MODEL_LEN - 1, '+', 1, MODEL_LEN))
        rows.append(('TIR_R', right_start, end, '-', 1, MODEL_LEN))

    # An orphan left terminus, partially matching its model.
    orphan_start = 35_000
    seq[orphan_start - 1 : orphan_start - 1 + 40] = list(LEFT_MOTIF[20:])
    rows.append(('TIR_L', orphan_start, orphan_start + 39, '+', 21, MODEL_LEN))

    genome = tmp_path / 'genome.fa'
    genome.write_text('>chr1\n' + ''.join(seq) + '\n')

    # nhmmer --tblout: whitespace-delimited, columns read by index.
    lines = ['# target name accession query name accession hmmfrom hmm to ...']
    for model, start, end, strand, hmm_from, hmm_to in rows:
        ali_from, ali_to = (start, end) if strand == '+' else (end, start)
        lines.append(
            ' '.join(
                [
                    'chr1',
                    '-',
                    model,
                    '-',
                    str(hmm_from),
                    str(hmm_to),
                    str(ali_from),
                    str(ali_to),
                    str(ali_from),
                    str(ali_to),
                    str(CONTIG_LEN),
                    strand,
                    '1e-25',
                    '95.0',
                    '0.1',
                    '-',
                ]
            )
        )
    hits = tmp_path / 'hits.tbl'
    hits.write_text('\n'.join(lines) + '\n')

    lengths = tmp_path / 'lengths.txt'
    lengths.write_text(f'TIR_L\t{MODEL_LEN}\nTIR_R\t{MODEL_LEN}\n')

    pairing_map = tmp_path / 'pairs.tsv'
    pairing_map.write_text('TIR_L\tTIR_R\n')

    return {
        'genome': genome,
        'hits': hits,
        'lengths': lengths,
        'pairing_map': pairing_map,
    }


def make_args(run_inputs, outdir, **overrides):
    """Build a fully populated namespace for hmm_pair.main."""
    args = argparse.Namespace(
        loglevel='ERROR',
        logfile=None,
        genome=str(run_inputs['genome']),
        blastdb=None,
        nhmmer_file=str(run_inputs['hits']),
        left_nhmmer=None,
        right_nhmmer=None,
        blast_file=None,
        left_blast=None,
        right_blast=None,
        hmm_file=None,
        left_model=None,
        right_model=None,
        lengths_file=str(run_inputs['lengths']),
        query_len=None,
        maxeval=0.001,
        mincov=0.5,
        maxdist=None,
        max_offset=None,
        orientation='F,R',
        stable_reps=0,
        pairing_map=str(run_inputs['pairing_map']),
        outdir=str(outdir),
        prefix=None,
        nopairing=False,
        no_hits=False,
        no_elements=False,
        gff_out=True,
        gff_report='all',
        padlen=None,
        flanks=False,
        flanks_paired=False,
        flank_len=50,
        flank_max_offset=None,
        no_pad_flanks=False,
        extend_hits_to_model=False,
        insertion_site=False,
        tsd_length=0,
        tsd_length_map=None,
        tsd_in_model=False,
        tempdir=None,
        keep_temp=False,
        report=False,
        report_out=None,
        report_title=None,
        report_no_sequences=False,
        report_max_seq_mb=20.0,
        report_msa='anchor',
        report_msa_max_rows=500,
        report_max_hits=200000,
        report_max_rows=30,
    )
    for key, value in overrides.items():
        setattr(args, key, value)
    return args


def tree(root):
    """Map every file under `root` to its bytes, keyed by relative path."""
    return {
        str(path.relative_to(root)): path.read_bytes()
        for path in sorted(Path(root).rglob('*'))
        if path.is_file() and path.suffix != '.log'
    }


class TestReportIsWritten:
    def test_report_file_is_created(self, run_inputs, tmp_path):
        outdir = tmp_path / 'out'
        assert main(make_args(run_inputs, outdir, report=True)) == 0
        report = outdir / 'tirmite_pair_report.html'
        assert report.exists()
        html = report.read_text(encoding='utf-8')
        assert html.startswith('<!DOCTYPE html>')
        assert 'tirmite-report-data' in html

    def test_report_reflects_the_run(self, run_inputs, tmp_path):
        import json
        import re

        outdir = tmp_path / 'out'
        main(make_args(run_inputs, outdir, report=True))
        html = (outdir / 'tirmite_pair_report.html').read_text(encoding='utf-8')
        payload = json.loads(
            re.search(
                r'<script id="tirmite-report-data"[^>]*>(.*?)</script>', html, re.DOTALL
            ).group(1)
        )
        assert payload['hits']['n'] == 7
        assert len(payload['elements']) == len(ELEMENTS)
        assert [c['name'] for c in payload['contigs']] == ['chr1']
        # The contig length came from the genome, not from the hits.
        assert payload['contigs'][0]['length'] == CONTIG_LEN
        assert payload['contigs'][0]['length_source'] == 'source'
        assert payload['groups'][0]['group_id'] == 'TIR_L__TIR_R'
        # The orphan terminus is present and only partially matches its model.
        assert any(cov is not None and cov < 1 for cov in payload['hits']['model_cov'])

    def test_prefix_is_applied_to_the_report_name(self, run_inputs, tmp_path):
        outdir = tmp_path / 'out'
        main(make_args(run_inputs, outdir, report=True, prefix='run1'))
        assert (outdir / 'run1_tirmite_pair_report.html').exists()

    def test_report_out_sets_the_path(self, run_inputs, tmp_path):
        outdir = tmp_path / 'out'
        target = tmp_path / 'elsewhere' / 'custom.html'
        args = make_args(run_inputs, outdir, report_out=str(target))
        validate_arguments(args)
        assert args.report is True
        assert main(args) == 0
        assert target.exists()


class TestOtherOutputsAreUnchanged:
    def test_adding_report_changes_nothing_else(self, run_inputs, tmp_path):
        without = tmp_path / 'without'
        with_report = tmp_path / 'with'
        assert main(make_args(run_inputs, without)) == 0
        assert main(make_args(run_inputs, with_report, report=True)) == 0

        baseline = tree(without)
        after = tree(with_report)
        added = set(after) - set(baseline)
        assert added == {'tirmite_pair_report.html'}
        for name, content in baseline.items():
            assert after[name] == content, f'{name} changed when --report was added'


class TestDegradedModes:
    def test_no_elements_omits_sequences_and_warns(self, run_inputs, tmp_path, caplog):
        outdir = tmp_path / 'out'
        args = make_args(run_inputs, outdir, report=True, no_elements=True)
        validate_arguments(args)
        assert args.report_no_sequences is True
        assert 'element sequences cannot be embedded' in caplog.text.replace('\n', ' ')

        assert main(args) == 0
        html = (outdir / 'tirmite_pair_report.html').read_text(encoding='utf-8')
        assert '"embedded":false' in html.replace(' ', '')

    def test_nopairing_still_writes_a_hits_only_report(self, run_inputs, tmp_path):
        import json
        import re

        outdir = tmp_path / 'out'
        assert main(make_args(run_inputs, outdir, report=True, nopairing=True)) == 0
        html = (outdir / 'tirmite_pair_report.html').read_text(encoding='utf-8')
        payload = json.loads(
            re.search(
                r'<script id="tirmite-report-data"[^>]*>(.*?)</script>', html, re.DOTALL
            ).group(1)
        )
        assert payload['hits']['n'] == 7
        assert payload['elements'] == []
        assert payload['groups'] == []

    def test_msa_off_omits_the_panels(self, run_inputs, tmp_path):
        outdir = tmp_path / 'out'
        main(make_args(run_inputs, outdir, report=True, report_msa='off'))
        html = (outdir / 'tirmite_pair_report.html').read_text(encoding='utf-8')
        assert '"msa":[]' in html.replace(' ', '')

    def test_anchor_mode_builds_panels_without_mafft(self, run_inputs, tmp_path):
        import json
        import re

        outdir = tmp_path / 'out'
        main(make_args(run_inputs, outdir, report=True, report_msa='anchor'))
        html = (outdir / 'tirmite_pair_report.html').read_text(encoding='utf-8')
        payload = json.loads(
            re.search(
                r'<script id="tirmite-report-data"[^>]*>(.*?)</script>', html, re.DOTALL
            ).group(1)
        )
        assert [panel['aligner'] for panel in payload['msa']] == ['anchor', 'anchor']

    def test_no_report_flag_writes_nothing(self, run_inputs, tmp_path):
        outdir = tmp_path / 'out'
        main(make_args(run_inputs, outdir))
        assert not list(outdir.glob('*.html'))


class TestFailureIsolation:
    def test_a_failing_report_does_not_fail_the_run(
        self, run_inputs, tmp_path, monkeypatch
    ):
        import tirmite.report.render as render_module

        def boom(*args, **kwargs):
            raise RuntimeError('renderer exploded')

        monkeypatch.setattr(render_module, 'write_pair_report', boom)

        outdir = tmp_path / 'out'
        # The run's real outputs were written before the report was attempted,
        # so a broken visualisation must not turn into a failed analysis.
        assert main(make_args(run_inputs, outdir, report=True)) == 0
        assert (outdir / 'TIR_L_TIR_R').is_dir()
        assert not (outdir / 'tirmite_pair_report.html').exists()


class TestReportArgumentValidation:
    @pytest.mark.parametrize(
        'field,value',
        [
            ('report_max_seq_mb', 0),
            ('report_max_seq_mb', -1),
            ('report_max_hits', 0),
            ('report_max_rows', 0),
            ('report_msa_max_rows', 0),
        ],
    )
    def test_out_of_range_values_are_rejected(self, run_inputs, tmp_path, field, value):
        args = make_args(run_inputs, tmp_path, report=True, **{field: value})
        with pytest.raises(ValueError, match=field.replace('_', '-')):
            validate_arguments(args)

    def test_values_are_not_checked_without_report(self, run_inputs, tmp_path):
        # Nothing reads them, so an unused default should not block a run.
        args = make_args(run_inputs, tmp_path, report_max_hits=0)
        validate_arguments(args)

    def test_mafft_mode_without_mafft_is_rejected(
        self, run_inputs, tmp_path, monkeypatch
    ):
        monkeypatch.setattr('tirmite.runners.mafft.mafft_available', lambda: False)
        args = make_args(run_inputs, tmp_path, report=True, report_msa='mafft')
        with pytest.raises(ValueError, match='--report-msa auto'):
            validate_arguments(args)
