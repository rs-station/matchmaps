"""Compute unbiased real space difference map from inputs in different spacegroups"""

import os
import sys
import time
from pathlib import Path

import gemmi
import reciprocalspaceship as rs

from matchmaps._parsers import matchmaps_mr_parser
from matchmaps._phenix_utils import (
    rigid_body_refinement_wrapper,
    phaser_wrapper,
    _remove_waters,
)
from matchmaps._utils import (
    _handle_special_positions,
    make_floatgrid_from_mtz,
    _realspace_align_and_subtract,
    _rbr_selection_parser,
    _restore_ligand_occupancy,
    _validate_environment,
    _validate_inputs,
    _clean_up_files,
    _cif_or_pdb_to_pdb,
    _cif_or_mtz_to_mtz,
    _write_script,
    _validate_column_dtypes,
)


def compute_mr_difference_map(
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
    on_as_stationary=False,
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
    Compute a real-space difference map from mtzs in different spacegroups.

    Parameters
    ----------
    pdboff : string
        Name of input .pdb file to use for phasing and as an MR search model
    mtzoff : string
        Name of input .mtz containing 'off data
    mtzon : string
        Name of input .mtz file containing 'off' data
    Foff : string
        Column in mtzoff containing structure factor amplitudes
    SigFoff : string
        Column in mtzoff containing structure factor uncertainties
    Fon : string
        Column in mtzon containing structure factor amplitudes
    SigFon : string
        Column in mtzon containing structure factor uncertainties
    ligands : list of strings
        Filename(s) of any .cif ligand restraint files necessary for phenix.refine
        by default None, meaning only the .pdb is required for refinement
    dmin : float, optional
        Minimum resolution (in Angstroms) reflections to be used in computing real-space maps from mtzs.
        If omitted, resolution cutoff is the maximum resolution from the lower-resolution input file.
    spacing : float, optional
        Approximate size of real-space voxels in Angstroms, by default 0.5 A
    on_as_stationary : bool, optional
        If True, align "off" data onto "on" data, by default False
        Note that this applies only to post-molecular-replacement refinement, not to molecular replacement itself.
    input_dir : str, optional
        Path to directory containing input files, by default "./" (current directory)
    output_dir : str, optional
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
    keep_temp_files : bool, optional
        If not None, the name of a subdirectory of the output_dir into which intermediate matchmaps files are moved upon program completion.
    radius : float, optional
        Maximum distance away from protein model to include voxels. Only applies to the "unmasked" difference map output.
    alpha : float, optional
        Alpha to use in error weighting of F-obs prior to Fourier Transform. Defaults to 0, e.g. no weighting.
     no_bss : bool, optional
        If True, skip bulk solvent scaling feature of phenix.refine
    """

    auto_phenix_version = _validate_environment(ccp4=False)

    if not phenix_version:
        phenix_version = auto_phenix_version

    output_dir_contents = list(output_dir.glob("*"))

    pdboff = _cif_or_pdb_to_pdb(pdboff, output_dir)

    mtzoff, off_name = _cif_or_mtz_to_mtz(mtzoff, output_dir)
    mtzon, on_name = _cif_or_mtz_to_mtz(mtzon, output_dir)

    # take in the list of rbr selections and parse them into phenix and gemmi selection formats
    # if rbr_groups = None, just returns (None, None)
    rbr_phenix, rbr_gemmi = _rbr_selection_parser(rbr_selections)

    # this is where scaling takes place in the usual pipeline, but that doesn't make sense with different-spacegroup inputs

    pdboff = _handle_special_positions(pdboff, output_dir)

    # write this function as a wrapper around phenix.pdbtools
    # modified pdboff already moved to output_dir by _handle_special_positions
    pdboff = _remove_waters(pdboff, output_dir, verbose)

    print(
        f"{time.strftime('%H:%M:%S')}: Running phenix.phaser to place 'off' model into 'on' data..."
    )

    phaser_nickname = phaser_wrapper(
        mtzfile=mtzon,
        pdb=pdboff,
        output_dir=output_dir,
        off_labels=f"{Fon},{SigFon}",
        phenix_style=phenix_version,
        eff=None,
        verbose=verbose,
    )

    # TO-DO: fix ligand occupancies in pdb_mr_to_on
    edited_mr_pdb = _restore_ligand_occupancy(
        pdb_to_be_restored=str(phaser_nickname) + ".1.pdb",
        original_pdb=pdboff,
        output_dir=output_dir,
    )

    print(f"{time.strftime('%H:%M:%S')}: Running phenix.refine for the 'on' data...")

    nickname_on = rigid_body_refinement_wrapper(
        mtzon=mtzon,
        pdboff=edited_mr_pdb,
        input_dir=input_dir,
        output_dir=output_dir,
        ligands=ligands,
        eff=eff,
        verbose=verbose,
        rbr_selections=rbr_phenix,
        off_labels=f"{Fon},{SigFon}",  # workaround for compatibility
        mr_on=True,
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

    # from here down I just copied over the stuff from the normal version

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
    args = matchmaps_mr_parser.parse_args()

    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    if not os.path.exists(args.input_dir):
        raise ValueError(f"Input directory '{args.input_dir}' does not exist")

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

    compute_mr_difference_map(
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
            utility="matchmaps.mr",
            arguments=sys.argv[1:],
            script_name=args.script,
            phenix_version=args.phenix_version,
        )

    return


if __name__ == "__main__":
    main()
