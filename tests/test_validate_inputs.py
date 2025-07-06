from pathlib import Path

import pytest

from matchmaps._utils import _validate_inputs


@pytest.mark.parametrize(
    "inputs",
    (
        (
            "./",
            "./",
            None,
            "data/raw_data/1rx2_phases.mtz",
            "data/raw_data/1rx4_phases.mtz",
            "data/raw_data/1rx2.pdb",
        ),
        (
            "./",
            "./",
            ("data/raw_data/FOL.cif", "data/raw_data/NAP.cif"),
            "data/raw_data/1rx2_phases.mtz",
            "data/raw_data/1rx4_phases.mtz",
            "data/raw_data/1rx2.pdb",
        ),
        (
            "data/raw_data",
            "outputs",
            ("FOL.cif", "NAP.cif"),
            "1rx2_phases.mtz",
            "1rx4_phases.mtz",
            "1rx2.pdb",
        ),
    ),
)
def test_validate_inputs(inputs):
    (input_dir, output_dir, ligands, mtzoff, mtzon, pdboff) = _validate_inputs(*inputs)

    assert isinstance(input_dir, Path)
    assert isinstance(output_dir, Path)
    assert isinstance(mtzoff, Path)
    assert isinstance(mtzon, Path)
    assert isinstance(pdboff, Path)

    if ligands is not None:
        assert [isinstance(ligand, Path) for ligand in ligands]


def test_validate_inputs_missingfile():
    with pytest.raises(ValueError):
        _validate_inputs(
            "./",
            "./",
            None,
            "data/raw_data/1rx2_phases.mtz",
            "data/raw_data/1rx4_phases.mtz",
            "data/raw_data/1rx.pdb",
        )
