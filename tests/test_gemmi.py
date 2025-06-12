from pathlib import Path

import gemmi

from matchmaps._utils import make_floatgrid_from_mtz, _realspace_align_and_subtract


def test_gemmi_pipeline(phenix_outputs):
    mtzoff, pdboff, mtzon, pdbon = phenix_outputs

    spacing = 0.5
    F = "F-obs-filtered"
    SigF = "SIGF-obs-filtered"
    Phi = "PH2FOFCWT"
    spacegroup = "P1"
    dmin = max(
        min(mtzoff.compute_dHKL(inplace=True).dHKL),
        min(mtzon.compute_dHKL(inplace=True).dHKL),
    )
    alpha = 0

    fg_off = make_floatgrid_from_mtz(
        mtzoff,
        spacing=spacing,
        F=F,
        SigF=SigF,
        Phi=Phi,
        spacegroup=spacegroup,
        dmin=dmin,
        alpha=alpha,
    )
    fg_on = make_floatgrid_from_mtz(
        mtzon,
        spacing=spacing,
        F=F,
        SigF=SigF,
        Phi=Phi,
        spacegroup=spacegroup,
        dmin=dmin,
        alpha=alpha,
    )

    assert isinstance(fg_off, gemmi.FloatGrid)
    assert isinstance(fg_on, gemmi.FloatGrid)

    output_dir = Path("data/outputs")
    on_name = "on_test"
    off_name = "off_test"
    on_as_stationary = False
    rbr_gemmi = None
    radius = 5

    _realspace_align_and_subtract(
        fg_off=fg_off,
        fg_on=fg_on,
        pdboff=pdboff,
        pdbon=pdbon,
        output_dir=output_dir,
        on_name=on_name,
        off_name=off_name,
        on_as_stationary=on_as_stationary,
        selection=rbr_gemmi,
        radius=radius,
    )

    return
