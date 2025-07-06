from matchmaps._utils import _validate_environment
import pytest


@pytest.mark.parametrize("ccp4", (True, False))
def test_validate_environment_failure(ccp4):
    with pytest.raises(OSError):
        _validate_environment(ccp4=ccp4)
