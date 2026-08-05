"""
Tests that the report's templates and assets are reachable as package data.

These are not .py files, so they ship only because pyproject.toml lists them
under the wheel's `artifacts`. Without this test a packaging mistake would
surface as a broken install rather than a failing build.
"""

from importlib import resources

import pytest

from tirmite.report import render

EXPECTED_ASSETS = [
    'report.css',
    'report-core.js',
    'track-browser.js',
    'stats-table.js',
    'hit-tables.js',
    'msa-panel.js',
]
EXPECTED_TEMPLATES = [
    '_base.html.j2',
    '_stats_table.html.j2',
    'pair_report.html.j2',
    'search_report.html.j2',
]


@pytest.mark.parametrize('name', EXPECTED_ASSETS)
def test_asset_is_readable(name):
    text = render.inline_asset(name)
    assert text.strip()


@pytest.mark.parametrize('name', EXPECTED_TEMPLATES)
def test_template_is_readable(name):
    text = (
        resources.files('tirmite.report')
        .joinpath('templates', name)
        .read_text(encoding='utf-8')
    )
    assert text.strip()


def test_renderer_loads_every_asset_it_declares():
    # Guards against an asset being added to the renderer's list but not to
    # the package.
    for name in render._STYLES + render._SCRIPTS:
        assert render.inline_asset(name).strip()


def test_templates_resolve_through_the_jinja_loader():
    env = render._environment()
    for name in EXPECTED_TEMPLATES:
        assert env.get_template(name) is not None


def test_declared_templates_match_the_shipped_directory():
    # Catches a template added to the package but never listed here, which
    # would otherwise go untested for shipping.
    shipped = {
        entry.name
        for entry in resources.files('tirmite.report').joinpath('templates').iterdir()
        if entry.name.endswith('.j2')
    }
    assert shipped == set(EXPECTED_TEMPLATES)


def test_declared_assets_match_the_shipped_directory():
    shipped = {
        entry.name
        for entry in resources.files('tirmite.report').joinpath('assets').iterdir()
        if entry.name.endswith(('.css', '.js'))
    }
    assert shipped == set(EXPECTED_ASSETS)
