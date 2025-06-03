"""Compute unbiased real space difference map."""

import sys
import time
from pathlib import Path

import gemmi
import reciprocalspaceship as rs

from matchmaps._parsers import matchmaps_parser
from matchmaps._phenix_utils import (
    rigid_body_refinement_wrapper,
    _renumber_waters,
    _custom_subprocess,
)
from matchmaps._utils import (
    _handle_special_positions,
    make_floatgrid_from_mtz,
    _realspace_align_and_subtract,
    _rbr_selection_parser,
    _clean_up_files,
    _validate_environment,
    _validate_inputs,
    _cif_or_mtz_to_mtz,
    _cif_or_pdb_to_pdb,
    _write_script,
    _validate_column_dtypes,
)


def compute_realspace_difference_map(
    pdboff: Path,
    mtzoff: Path,
    mtzon: Path,
    Foff: str,
    SigFoff: str,
    Fon: str,
    SigFon: str,
    ligands: list = None,
    dmin: float = None,
    spacing=0.5,
    on_as_stationary: bool = False,
    input_dir=Path("."),
    output_dir=Path("."),
    verbose=False,
    rbr_selections: list[str] = None,
    eff: str = None,
    keep_temp_files: str = None,
    radius: float = 5,
    alpha: float = 0,
    no_bss=False,
    phenix_version: str = None,
):
    """
    Compute a real-space difference map from mtzs.

    Parameters
    ----------
    pdboff : pathlib.Path
        Name of input .pdb file to use for phasing
    mtzoff : pathlib.Path
        Name of input .mtz file containing 'off' data
    mtzon : pathlib.Path
        Name of input .mtz file containing 'off' data
    Foff : string
        Column in mtzoff containing structure factor amplitudes
    SigFoff : string
        Column in mtzoff containing structure factor uncertainties
    Fon : string
        Column in mtzon containing structure factor amplitudes
    SigFon : string
        Column in mtzon containing structure factor uncertainties
    ligands : list of pathlib.Path
        Filename(s) of any .cif ligand restraint files necessary for phenix.refine
        by default None, meaning only the .pdb is required for refinement
    dmin : float, optional
        Minimum resolution (in Angstroms) reflections to be used in computing real-space maps from mtzs.
        If omitted, resolution cutoff is the maximum resolution from the lower-resolution input file.
    spacing : float, optional
        Approximate size of real-space voxels in Angstroms, by default 0.5 A
    on_as_stationary : bool, optional
        If True, align "off" data onto "on" data, by default False
    input_dir : pathlib.Path, optional
        Path to directory containing input files, by default "./" (current directory)
    output_dir : pathlib.Path, optional
        Path to directory to which output files should be written, by default "./" (current directory)
    verbose : bool, optional
        If True, print outputs of scaleit and phenix.refine, by default False
    rbr_selections : list of strings, optional
        Custom selections to provide to refinement.refine.sites.rigid_body=
        If omitted, then refinement.refine.sites.rigid_body=all, and the entire structure is refined as a single rigid body.
    eff : str, optional
        Name of a file containing a template .eff parameter file for phenix.refine.
        If omitted, the sensible built-in .eff template is used. If you need to change something,
        I recommend copying the template from the source code and editing that.
    keep_temp_files : str, optional
        If not None, the name of a subdirectory of the output_dir into which intermediate matchmaps files are moved upon program completion.
    radius : float, optional
        Maximum distance away from protein model to include voxels. Only applies to the "unmasked" difference map output.
    alpha : float, optional
        Alpha to use in error weighting of F-obs prior to Fourier Transform. Defaults to 0, e.g. no weighting.
    no_bss : bool, optional
        If True, skip bulk solvent scaling feature of phenix.refine
    phenix_version: str, optional
        Phenix version string to override the automatically detected version. I don't know why this would be necessary.
    """

    auto_phenix_version = _validate_environment(ccp4=True)

    if not phenix_version:
        phenix_version = auto_phenix_version

    output_dir_contents = list(output_dir.glob("*"))

    pdboff = _cif_or_pdb_to_pdb(pdboff, output_dir)

    mtzoff, off_name = _cif_or_mtz_to_mtz(mtzoff, output_dir)
    mtzon, on_name = _cif_or_mtz_to_mtz(mtzon, output_dir)

    # take in the list of rbr selections and parse them into phenix and gemmi selection formats
    # if rbr_groups = None, just returns (None, None)
    rbr_phenix, rbr_gemmi = _rbr_selection_parser(rbr_selections)

    ### scaleit
    mtzon_scaled = output_dir / (on_name + "_scaled.mtz")

    print(
        f"{time.strftime('%H:%M:%S')}: Running scaleit to scale 'on' data to 'off' data..."
    )

    _custom_subprocess(
        command="rs.scaleit",
        params=f"-r {mtzoff} {Foff} {SigFoff} -i {mtzon} {Fon} {SigFon} -o {mtzon_scaled} --ignore-isomorphism",
        verbose=verbose,
    )

    # now that scaleit has run, let's swap out the spacegroup from the scaled file
    mtzon_scaled_py = rs.read_mtz(str(mtzon_scaled))
    mtzon_original_py = rs.read_mtz(str(mtzon))
    mtzoff_original_py = rs.read_mtz(str(mtzoff))

    mtzoff_trunc = output_dir / (mtzoff.name.removesuffix(".mtz") + "_trunc.mtz")
    mtzon_scaled_truecell = output_dir / (
        mtzon_scaled.name.removesuffix(".mtz") + "_truecell.mtz"
    )

    mtzon_scaled_py.cell = mtzon_original_py.cell

    mtzoff_original_py.compute_dHKL(inplace=True)
    mtzon_scaled_py.compute_dHKL(inplace=True)

    # make resolutions match for mtzon_scaled_py and mtzon_original_py
    resolution = max(mtzoff_original_py["dHKL"].min(), mtzon_scaled_py["dHKL"].min())
    mtzoff_original_py = mtzoff_original_py.loc[mtzoff_original_py.dHKL >= resolution]
    mtzon_scaled_py = mtzon_scaled_py.loc[mtzon_scaled_py.dHKL >= resolution]

    mtzoff_original_py.write_mtz(str(mtzoff_trunc))
    mtzon_scaled_py.write_mtz(str(mtzon_scaled_truecell))

    # reset short nicknames to the latest files
    mtzon = mtzon_scaled_truecell
    mtzoff = mtzoff_trunc
    # done with cell swapping and resolution matching

    pdboff = _handle_special_positions(pdboff, output_dir)

    pdboff = _renumber_waters(pdboff, verbose)

    print(f"{time.strftime('%H:%M:%S')}: Running phenix.refine for the 'on' data...")

    nickname_on = rigid_body_refinement_wrapper(
        mtzon=mtzon,
        pdboff=pdboff,
        input_dir=input_dir,
        output_dir=output_dir,
        ligands=ligands,
        eff=eff,
        verbose=verbose,
        rbr_selections=rbr_phenix,
        no_bss=no_bss,
        phenix_style=phenix_version,
    )

    print(f"{time.strftime('%H:%M:%S')}: Running phenix.refine for the 'off' data...")

    nickname_off = rigid_body_refinement_wrapper(
        mtzon=mtzoff,
        pdboff=pdboff,
        input_dir=input_dir,
        output_dir=output_dir,
        ligands=ligands,
        eff=eff,
        verbose=verbose,
        rbr_selections=rbr_phenix,
        off_labels=f"{Foff},{SigFoff}",
        no_bss=no_bss,
        phenix_style=phenix_version,
    )

    # read back in the files created by phenix
    # these have knowable names
    mtzon = rs.read_mtz(f"{nickname_on}_1.mtz")
    mtzoff = rs.read_mtz(f"{nickname_off}_1.mtz")

    pdbon = gemmi.read_structure(f"{nickname_on}_1.pdb")
    pdboff = gemmi.read_structure(f"{nickname_off}_1.pdb")

    if dmin is None:
        dmin = max(
            min(mtzoff.compute_dHKL(inplace=True).dHKL),
            min(mtzon.compute_dHKL(inplace=True).dHKL),
        )

    print(f"{time.strftime('%H:%M:%S')}: Constructing FloatGrids from mtzs...")
    # hard-coding F, Phi because they're always phenix outputs
    # TO-DO: Figure out why phenix outputs are sometimes still split into (+) and (-) columns, even when I specify that anomalous=False
    # As a workaround, even anomalous files have a single 'F-obs-filtered' column, so I can always just use that.
    fg_off = make_floatgrid_from_mtz(
        mtzoff,
        spacing,
        F="F-obs-filtered",
        SigF="SIGF-obs-filtered",
        Phi="PH2FOFCWT",
        spacegroup="P1",
        dmin=dmin,
        alpha=alpha,
    )
    fg_on = make_floatgrid_from_mtz(
        mtzon,
        spacing,
        F="F-obs-filtered",
        SigF="SIGF-obs-filtered",
        Phi="PH2FOFCWT",
        spacegroup="P1",
        dmin=dmin,
        alpha=alpha,
    )

    if rbr_gemmi is None:
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

    else:  # run helper function separately for each selection
        for n, selection in enumerate(rbr_gemmi, start=1):
            on_name_rbr = on_name + "_rbrgroup" + str(n)
            off_name_rbr = off_name + "_rbrgroup" + str(n)

            _realspace_align_and_subtract(
                fg_off=fg_off.clone(),
                fg_on=fg_on.clone(),
                pdboff=pdboff,
                pdbon=pdbon,
                output_dir=output_dir,
                on_name=on_name_rbr,
                off_name=off_name_rbr,
                on_as_stationary=on_as_stationary,
                selection=selection,
                radius=radius,
            )

    print(f"{time.strftime('%H:%M:%S')}: Cleaning up files...")
    _clean_up_files(output_dir, output_dir_contents, keep_temp_files)

    print(f"{time.strftime('%H:%M:%S')}: Done!")

    return


def main():
    args = matchmaps_parser.parse_args()

    (input_dir, output_dir, ligands, mtzoff, mtzon, pdboff) = _validate_inputs(
        args.input_dir,
        args.output_dir,
        args.ligands,
        args.mtzoff[0],
        args.mtzon[0],
        args.pdboff,
    )

    _validate_column_dtypes(
        rs.read_mtz(str(mtzoff)),
        (args.mtzoff[1], args.mtzoff[2]),
        (rs.StructureFactorAmplitudeDtype, rs.StandardDeviationDtype),
    )
    _validate_column_dtypes(
        rs.read_mtz(str(mtzon)),
        (args.mtzon[1], args.mtzon[2]),
        (rs.StructureFactorAmplitudeDtype, rs.StandardDeviationDtype),
    )

    compute_realspace_difference_map(
        pdboff=pdboff,
        ligands=ligands,
        mtzoff=mtzoff,
        mtzon=mtzon,
        Foff=args.mtzoff[1],
        SigFoff=args.mtzoff[2],
        Fon=args.mtzon[1],
        SigFon=args.mtzon[2],
        input_dir=input_dir,
        output_dir=output_dir,
        verbose=args.verbose,
        rbr_selections=args.rbr_selections,
        eff=args.eff,
        dmin=args.dmin,
        spacing=args.spacing,
        radius=args.unmasked_radius,
        alpha=args.alpha,
        on_as_stationary=args.on_as_stationary,
        keep_temp_files=args.keep_temp_files,
        no_bss=args.no_bss,
        phenix_version=args.phenix_version,
    )

    if args.script:
        _write_script(
            utility="matchmaps",
            arguments=sys.argv[1:],
            script_name=args.script,
            phenix_version=args.phenix_version,
        )

    return


if __name__ == "__main__":
    main()
