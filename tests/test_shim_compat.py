"""Guards on the tirmite.tirmitetools backwards-compatibility shim.

tirmitetools.py used to hold every core algorithm. Its contents now live under
tirmite.core, and it re-exports them so that existing imports keep working.
These tests pin that contract: if a name is moved again without updating the
shim, they fail rather than the failure surfacing as a user's ImportError.
"""

import importlib

import pytest

import tirmite.tirmitetools as shim

# Which core module each name is expected to live in now. This doubles as
# documentation of the split.
EXPECTED_HOME = {
    'tirmite.core.parsers': [
        'convertAlign',
        'import_nhmmer',
        'import_BED',
        'import_blast',
        'detect_input_format',
    ],
    'tirmite.core.filters': ['filterHitsLen', 'filterHitsEval'],
    'tirmite.core.termini': [
        'flipTIRs',
        '_model_deficit',
        '_determine_terminus_type',
        'TerminusAssignment',
        '_pair_roles',
        'resolve_terminus',
    ],
    'tirmite.core.extraction': [
        'fetch_padded_hit',
        '_model_extension',
        'extractTIRs',
        'writeTIRs',
        'fetchElements',
        'writeElements',
        'writePairedTIRs',
    ],
    'tirmite.core.flanks': [
        'compute_flank_coordinates',
        'compute_inner_tsd_coordinates',
        'FlankResult',
        'extract_terminus_flank',
        'writeFlanks',
    ],
    'tirmite.core.tsd': [
        'hamming_distance',
        'load_tsd_length_map',
        'reconstruct_target_site',
        'compare_tsds',
        'format_interleaved_flanks',
        'writeTargetSites',
    ],
    'tirmite.core.output': ['fetchUnpaired', 'gffWrite'],
    'tirmite.core.pairing': [
        'table2dict',
        'parseHits',
        'isfirstUnpaired',
        'getPairs',
        'countUnpaired',
        'listunpaired',
        'iterateGetPairs',
        'VALID_ORIENTATION_CODES',
        'parse_orientation',
        'PairingConfig',
        'parseHitsGeneral',
        'inter_hit_distance',
        'candidate_separation',
        '_check_distance',
        '_find_candidates',
        'iterateGetPairsAsymmetric',
        'getPairsAsymmetric',
        'checkAsymmetricReciprocity',
        'iterateGetPairsCustom',
        'getPairsSymmetric',
        'checkSymmetricReciprocity',
        'countUnpairedAsymmetric',
        'listunpairedAsymmetric',
        'countUnpairedSymmetric',
        'listunpairedSymmetric',
    ],
}

ALL_NAMES = [name for names in EXPECTED_HOME.values() for name in names]

# (module, name) pairs, for per-name parametrisation.
HOMED_NAMES = [
    pytest.param(module, name, id=f'{module.rsplit(".", 1)[-1]}.{name}')
    for module, names in EXPECTED_HOME.items()
    for name in names
]


class TestShimSurface:
    """The shim must expose exactly the names it used to."""

    def test_all_is_complete(self):
        """__all__ lists every re-exported name and nothing else."""
        assert sorted(shim.__all__) == sorted(ALL_NAMES)

    def test_no_duplicate_entries_in_all(self):
        """A name listed twice would hide a bad merge."""
        assert len(shim.__all__) == len(set(shim.__all__))

    @pytest.mark.parametrize('name', ALL_NAMES)
    def test_attribute_access_style_works(self, name):
        """`import tirmite.tirmitetools as tirmite; tirmite.NAME` still works.

        This is how test_pairing_orientation.py and several CLI modules use it.
        """
        assert hasattr(shim, name)

    @pytest.mark.parametrize('name', ALL_NAMES)
    def test_from_import_style_works(self, name):
        """`from tirmite.tirmitetools import NAME` still works."""
        module = importlib.import_module('tirmite.tirmitetools')
        assert getattr(module, name) is not None


class TestNamesResolveToTheirNewHome:
    """Re-exports must be the same objects, not copies."""

    @pytest.mark.parametrize('module_name,name', HOMED_NAMES)
    def test_shim_name_is_the_core_object(self, module_name, name):
        """The shim binds the identical object defined in tirmite.core."""
        module = importlib.import_module(module_name)
        assert getattr(shim, name) is getattr(module, name)


class TestCoreModulesAreIndependentlyImportable:
    """Each core module must import on its own, without cycles."""

    @pytest.mark.parametrize('module_name', sorted(EXPECTED_HOME))
    def test_module_imports_alone(self, module_name):
        """Importing one core module does not require the others."""
        assert importlib.import_module(module_name) is not None

    def test_core_does_not_import_cli(self):
        """core must never depend on cli; the dependency runs one way only."""
        import pkgutil

        import tirmite.core

        offenders = []
        for info in pkgutil.iter_modules(tirmite.core.__path__):
            module = importlib.import_module(f'tirmite.core.{info.name}')
            source_file = module.__file__
            assert source_file is not None
            with open(source_file) as handle:
                text = handle.read()
            if 'from tirmite.cli' in text or 'import tirmite.cli' in text:
                offenders.append(info.name)

        assert offenders == [], f'core modules importing cli: {offenders}'
