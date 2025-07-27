from tirmite.utils.utils import cleanID, getTimestring


def test_cleanID():
    """Test the cleanID function properly cleans strings."""
    assert cleanID('test-string!') == 'teststring'
    assert cleanID('test string') == 'test_string'
    assert cleanID('test_string') == 'test_string'
