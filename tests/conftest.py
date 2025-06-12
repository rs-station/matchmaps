import pytest

import gemmi
import reciprocalspaceship as rs


@pytest.fixture
def phenix_outputs():
    """
    Return a 4-tuple of the mtz and pdb for off and on datasets output by phenix.refine
    Data was produced by running matchmaps on PDB IDs 1RX2 (off) and 1RX4 (on)
    """

    d = "data/mm_intermediates"

    mtzoff = rs.read_mtz(f"{d}/1rx2_phases_trunc_rbr_to_self_1_1.mtz")
    pdboff = gemmi.read_structure(f"{d}/1rx2_phases_trunc_rbr_to_self_1_1.pdb")

    mtzon = rs.read_mtz(
        f"{d}/1rx4_phases_scaled_truecell_rbr_to_1rx2_nospecialpositions_renumbered_1_1.mtz"
    )
    pdbon = gemmi.read_structure(
        f"{d}/1rx4_phases_scaled_truecell_rbr_to_1rx2_nospecialpositions_renumbered_1_1.pdb"
    )

    return mtzoff, pdboff, mtzon, pdbon


@pytest.fixture
def mtz_1rx2():
    return rs.read_mtz("data/raw_data/1rx2_phases.mtz")
