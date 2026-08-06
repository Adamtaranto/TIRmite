"""Tests for tirmite.report.render — HTML assembly.

These assert structural facts about the document rather than comparing it byte
for byte: a golden HTML file would fail on every wording or styling change
without telling anyone anything useful.
"""

from html.parser import HTMLParser
import json
import re

import pytest
from test_report_model import make_element, make_hit, make_index

from tirmite.core.pairing import PairingConfig
from tirmite.report.collect import PairReportAccumulator
from tirmite.report.render import inline_asset, render_report, write_pair_report

PAYLOAD_RE = re.compile(
    r'<script id="tirmite-report-data" type="application/json">(.*?)</script>',
    re.DOTALL,
)


class Wellformed(HTMLParser):
    """Parses a document and records the tags it saw, to catch broken markup."""

    def __init__(self):
        super().__init__(convert_charrefs=True)
        self.tags = []

    def handle_starttag(self, tag, attrs):
        self.tags.append(tag)


def build(hit_table_factory, **finalise_kwargs):
    """Build a small two-hit, one-element report."""
    hits = [
        make_hit(0, 'TIR_L', 'chr1', 100, 200),
        make_hit(1, 'TIR_R', 'chr1', 5000, 5100, strand='-'),
        make_hit(2, 'TIR_L', 'chr9', 10, 90),
    ]
    table = hit_table_factory(
        [
            {'model': 'TIR_L', 'target': 'chr1', 'hitStart': '100', 'hitEnd': '200'},
            {
                'model': 'TIR_R',
                'target': 'chr1',
                'hitStart': '5000',
                'hitEnd': '5100',
                'strand': '-',
            },
            {'model': 'TIR_L', 'target': 'chr9', 'hitStart': '10', 'hitEnd': '90'},
        ]
    )
    acc = PairReportAccumulator(
        tirmite_version='9.9.9',
        command='tirmite pair --report',
        title='Test report',
        hit_table=table,
        model_lengths={'TIR_L': 100, 'TIR_R': 100},
        contig_length=lambda name: 50_000,
    )
    index = make_index(hits)
    acc.add_group(
        left_feature='TIR_L',
        right_feature='TIR_R',
        config=PairingConfig(
            orientation='F,R', left_model='TIR_L', right_model='TIR_R'
        ),
        paired={'TIR_L': [{0, 1}]},
        hit_index=index,
        elements={'TIR_L': [make_element(hits[0], hits[1], seq='ACGTACGTAC')]},
    )
    acc.add_unpaired(index)
    return acc.finalise(**finalise_kwargs)


@pytest.fixture
def html(hit_table_factory):
    return render_report(build(hit_table_factory))


class TestSelfContained:
    def test_no_external_stylesheets_or_scripts(self, html):
        # A report is routinely opened from disk long after the run, so
        # anything fetched over the network would be missing when it matters.
        assert '<link rel="stylesheet"' not in html
        assert 'rel="stylesheet"' not in html
        assert '<script src=' not in html

    def test_nothing_is_fetched_over_the_network(self, html):
        # Any src/href pointing off-document is a resource the report would
        # try to load. The SVG and XLink namespace URIs are declarations, not
        # fetches, and are the only permitted absolute URLs.
        allowed = {'http://www.w3.org/2000/svg', 'http://www.w3.org/1999/xlink'}
        targets = re.findall(r'(?:src|href|xlink:href)="([^"]*)"', html)
        remote = [
            target
            for target in targets
            if target.startswith(('http://', 'https://', '//'))
            and target not in allowed
        ]
        assert remote == []
        assert '@import' not in html
        assert 'url(http' not in html

    def test_styles_and_scripts_are_inlined(self, html):
        assert '<style>' in html
        assert '.track svg' in html
        assert 'window.TIRmite' in html

    def test_document_is_wellformed(self, html):
        parser = Wellformed()
        parser.feed(html)
        for tag in ('html', 'head', 'body', 'style', 'script', 'table', 'dialog'):
            assert tag in parser.tags


class TestPayload:
    def test_exactly_one_payload_block(self, html):
        assert len(PAYLOAD_RE.findall(html)) == 1

    def test_payload_matches_the_report_data(self, hit_table_factory):
        data = build(hit_table_factory)
        payload = json.loads(PAYLOAD_RE.search(render_report(data)).group(1))
        assert payload == json.loads(data.to_json())

    def test_payload_carries_the_hits(self, html):
        payload = json.loads(PAYLOAD_RE.search(html).group(1))
        assert payload['hits']['n'] == 3
        assert payload['elements'][0]['left_uid'] == 0


class TestContent:
    def test_every_group_appears_in_the_legend(self, html):
        # A legend is always present: colour alone never carries identity.
        assert 'TIR_L &lt;-&gt; TIR_R' in html or 'TIR_L <-> TIR_R' in html

    def test_group_colours_are_exposed_as_custom_properties(self, html):
        assert '--group-0:' in html
        assert 'var(--group-0)' in html
        # Dark steps are selected, not derived from the light ones.
        assert 'prefers-color-scheme: dark' in html

    def test_contig_with_hits_is_present_and_one_without_is_not(self, html):
        payload = json.loads(PAYLOAD_RE.search(html).group(1))
        names = [c['name'] for c in payload['contigs']]
        assert 'chr1' in names
        assert 'chr9' in names
        assert 'chr_absent' not in html

    def test_stats_tables_are_rendered_and_sortable(self, html):
        assert html.count('data-sortable') >= 3
        assert 'Pairing groups' in html
        assert 'Terminus models' in html
        assert 'Sequences' in html

    def test_hit_and_element_tables_have_mount_points(self, html):
        # These are built in the browser from the payload: writing hundreds of
        # thousands of rows into the document would make the file unusable.
        assert 'id="all-hits-table"' in html
        assert 'id="elements-table"' in html
        assert 'All terminus hits' in html
        assert 'Predicted elements' in html

    def test_element_table_is_omitted_when_there_are_no_elements(
        self, hit_table_factory
    ):
        table = hit_table_factory([{'model': 'TIR', 'target': 'chr1'}])
        acc = PairReportAccumulator(hit_table=table, contig_length=lambda n: 1000)
        acc.add_unpaired(make_index([make_hit(0, 'TIR', 'chr1', 100, 200)]))
        html = render_report(acc.finalise())
        assert 'id="all-hits-table"' in html
        assert 'id="elements-table"' not in html

    def test_metadata_is_shown(self, html):
        assert 'Test report' in html
        assert '9.9.9' in html
        assert 'tirmite pair --report' in html


class TestSequences:
    def test_sequences_present_by_default(self, html):
        payload = json.loads(PAYLOAD_RE.search(html).group(1))
        assert payload['sequences']['embedded'] is True
        assert list(payload['sequences']['seq'].values()) == ['ACGTACGTAC']

    def test_sequences_absent_when_disabled(self, hit_table_factory):
        data = build(hit_table_factory, embed_sequences=False)
        payload = json.loads(PAYLOAD_RE.search(render_report(data)).group(1))
        assert payload['sequences']['embedded'] is False
        assert payload['sequences']['seq'] == {}
        assert 'ACGTACGTAC' not in render_report(data)


class TestElementPopup:
    """
    The popup's sequence tab is built in JavaScript, so these check the
    shipped script rather than the rendered document: the assertions are that
    the behaviour is present in what the report carries, not that it ran.
    """

    @pytest.fixture
    def script(self):
        return inline_asset('track-browser.js')

    def test_download_is_a_link_carrying_the_filename(self, script):
        # An anchor rather than a button, so it offers Save link as and names
        # the file it will write before you commit to the click.
        assert "createElement('a')" in script
        assert "download.className = 'button-link'" in script
        assert "element.element_id + '.fasta'" in script
        assert 'createObjectURL' in script

    def test_download_link_is_styled_like_the_buttons_beside_it(self):
        assert 'a.button-link {' in inline_asset('report.css')

    def test_panel_shows_the_fasta_record_not_bare_sequence(self, script):
        # What is on screen is what the copy and download actions produce.
        assert 'pre.textContent = record()' in script
        assert "return header() + '\\n' + T.wrapSequence(current(), 60)" in script

    def test_clicking_an_annotation_opens_the_popup(self, script):
        # The same openFeature an element name in the tables calls, so both
        # routes land on one popup rather than two near-identical ones.
        assert "addEventListener('click', function (event) {" in script
        assert 'openFeature(h)' in script

    def test_hit_testing_is_forgiving_and_shared_with_hover(self, script):
        # A hit is drawn at least MIN_GLYPH_PX wide -- two pixels at
        # whole-contig zoom -- which is not a target anyone can hit. Hover uses
        # the same test, so a summary appearing means a click there will work.
        assert 'var CLICK_SLOP_X' in script
        assert 'var CLICK_SLOP_Y' in script
        assert script.count('self.hitAt(event)') >= 2
        assert 'var h = this.hitAt(event);' in script

    def test_vertical_tolerance_cannot_reach_two_rows_at_once(self, script):
        # Grabbing the row above or below would be worse than grabbing
        # nothing, and a point within tolerance of both rows would make the
        # choice arbitrary. Twice the slop has to fit in the gap between rows.
        def constant(name):
            return int(re.search(r'var ' + name + r' = (\d+)', script).group(1))

        gap = constant('ROW_HEIGHT') - constant('GLYPH_HEIGHT')
        assert 2 * constant('CLICK_SLOP_Y') < gap

    def test_a_pan_does_not_open_a_popup(self, script):
        # A drag ends in a click, and with a forgiving radius that click could
        # land on a glyph the reader never aimed at.
        assert 'var DRAG_THRESHOLD_PX' in script
        assert 'if (travelled > DRAG_THRESHOLD_PX) return;' in script

    def test_blob_is_released_when_the_popup_goes_away(self, script):
        assert 'revokeObjectURL' in script
        # Closing, dismissing with Escape and opening the next feature all
        # have to run the teardown; only the first goes through closeModal.
        assert script.count('runModalCleanups()') >= 2
        assert "modal.addEventListener('close', runModalCleanups)" in script


class TestEscaping:
    def test_hostile_model_name_cannot_close_the_payload_element(
        self, hit_table_factory
    ):
        hits = [make_hit(0, 'a</script><img src=x>b', 'chr1', 100, 200)]
        table = hit_table_factory(
            [{'model': 'a</script><img src=x>b', 'target': 'chr1'}]
        )
        acc = PairReportAccumulator(hit_table=table, contig_length=lambda n: 1000)
        acc.add_unpaired(make_index(hits))
        html = render_report(acc.finalise())

        payload = PAYLOAD_RE.search(html)
        assert payload is not None
        # The name survives intact inside the JSON...
        assert json.loads(payload.group(1))['models'][0]['name'] == (
            'a</script><img src=x>b'
        )
        # ...but never appears as live markup.
        assert '<img src=x>' not in html

    def test_hostile_contig_name_is_escaped_in_the_body(self, hit_table_factory):
        name = 'chr<script>alert(1)</script>'
        hits = [make_hit(0, 'TIR', name, 100, 200)]
        table = hit_table_factory([{'model': 'TIR', 'target': name}])
        acc = PairReportAccumulator(hit_table=table, contig_length=lambda n: 1000)
        acc.add_unpaired(make_index(hits))
        html = render_report(acc.finalise())
        assert '<script>alert(1)</script>' not in html


class TestAssets:
    def test_missing_asset_names_the_packaging_problem(self):
        with pytest.raises(FileNotFoundError) as excinfo:
            inline_asset('does-not-exist.css')
        assert 'artifacts' in str(excinfo.value)

    def test_closing_tags_in_assets_are_neutralised(self, monkeypatch):
        # An asset containing '</script' would end its element early.
        import tirmite.report.render as render_module

        class FakeResource:
            def joinpath(self, *_parts):
                return self

            def read_text(self, encoding='utf-8'):
                return 'var a = "</script>";'

        monkeypatch.setattr(
            render_module.resources, 'files', lambda _pkg: FakeResource()
        )
        assert '</script>' not in inline_asset('anything.js')


class TestWriteReport:
    def test_writes_and_creates_parent_directories(self, hit_table_factory, tmp_path):
        target = tmp_path / 'nested' / 'report.html'
        written = write_pair_report(build(hit_table_factory), target)
        assert written == target
        assert target.read_text(encoding='utf-8').startswith('<!DOCTYPE html>')
