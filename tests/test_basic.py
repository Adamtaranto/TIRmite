import pytest
from tirmite.utils import cleanID, getTimestring

def test_cleanID():
    """Test the cleanID function properly cleans strings."""
    assert cleanID("test-string!") == "teststring"
    assert cleanID("test string") == "test_string"
    assert cleanID("test_string") == "test_string"

def test_getTimestring():
    """Test that getTimestring returns a string."""
    result = getTimestring()
    assert isinstance(result, str)
    assert len(result) > 0