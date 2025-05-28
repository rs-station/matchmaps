import pytest

import gemmi
import reciprocalspaceship as rs


@pytest.fixture
def phenix_outputs():

    d = "data/mm_intermediates"

    mtzoff = rs.read_mtz(f"{d}/1rx2_phases_trunc_rbr_to_self_1_1.mtz")
    pdboff = gemmi.read_structure(f"{d}/1rx2_phases_trunc_rbr_to_self_1_1.pdb")

    mtzon = rs.read_mtz(f"{d}/1rx4_phases_scaled_truecell_rbr_to_1rx2_nospecialpositions_renumbered_1_1.mtz")
    pdbon = gemmi.read_structure(f"{d}/1rx4_phases_scaled_truecell_rbr_to_1rx2_nospecialpositions_renumbered_1_1.pdb")

    return mtzoff, pdboff, mtzon, pdbon