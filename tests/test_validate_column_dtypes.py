import reciprocalspaceship as rs
import pytest

from matchmaps._utils import _validate_column_dtypes


def test_validate_column_dtypes(mtz_1rx2):
    _validate_column_dtypes(
        mtz_1rx2,
        ("FP", "SIGFP"),
        (rs.StructureFactorAmplitudeDtype, rs.StandardDeviationDtype),
    )
    _validate_column_dtypes(
        mtz_1rx2, ("FC", "PHIC"), (rs.StructureFactorAmplitudeDtype, rs.PhaseDtype)
    )


def test_validate_column_dtypes_errors(mtz_1rx2):
    with pytest.raises(ValueError):
        _validate_column_dtypes(
            mtz_1rx2, ("FP", "SIGFP"), (rs.StructureFactorAmplitudeDtype, rs.PhaseDtype)
        )

    with pytest.raises(KeyError):
        _validate_column_dtypes(
            mtz_1rx2,
            ("F", "SIGFP"),
            (rs.StructureFactorAmplitudeDtype, rs.StandardDeviationDtype),
        )
