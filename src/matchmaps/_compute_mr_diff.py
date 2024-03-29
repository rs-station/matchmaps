"""Compute unbiased real space difference map from inputs in different spacegroups"""

import argparse
import os
import sys
import subprocess
import time
from functools import partial
from pathlib import Path

import gemmi
import numpy as np
import reciprocalspaceship as rs

from matchmaps._utils import (
    _handle_special_positions,
    make_floatgrid_from_mtz,
    rigid_body_refinement_wrapper,
    _realspace_align_and_subtract,
    _rbr_selection_parser,
    _remove_waters,
    _restore_ligand_occupancy,
    _validate_environment,
    _validate_inputs,
    phaser_wrapper,
    _clean_up_files,
    _cif_or_pdb_to_pdb,
    _cif_or_mtz_to_mtz,
    _write_script,
)


def compute_mr_difference_map(
    pdboff : Path,
    mtzoff : Path,
    mtzon : Path,
    Foff : str,
    SigFoff : str,
    Fon : str,
    SigFon : str,
    ligands: list = None,
    dmin : int = None,
    spacing = 0.5,
    on_as_stationary = False,
    input_dir = Path("."),
    output_dir = Path("."),
    verbose = False,
    rbr_selections : list[str] = None,
    eff : str = None,
    keep_temp_files: str = None,
    radius : float = 5,
    alpha : float = 0,
    no_bss = False,
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
    
    _validate_environment(ccp4=False)

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
    pdboff = _remove_waters(pdboff, output_dir)

    print(
        f"{time.strftime('%H:%M:%S')}: Running phenix.phaser to place 'off' model into 'on' data..."
    )

    phaser_nickname = phaser_wrapper(
        mtzfile=mtzon,
        pdb=pdboff,
        input_dir=input_dir,
        output_dir=output_dir,
        off_labels=f"{Fon},{SigFon}",
        eff=None,
        verbose=verbose,
    )

    # TO-DO: fix ligand occupancies in pdb_mr_to_on
    edited_mr_pdb = _restore_ligand_occupancy(
        pdb_to_be_restored= str(phaser_nickname) + ".1.pdb",
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
        mr_off=True,
        no_bss=no_bss,
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
        mtzoff, spacing, F="F-obs-filtered", SigF="SIGF-obs-filtered", Phi="PH2FOFCWT", spacegroup="P1", dmin=dmin, alpha=alpha,
    )
    fg_on = make_floatgrid_from_mtz(
        mtzon, spacing, F="F-obs-filtered", SigF="SIGF-obs-filtered", Phi="PH2FOFCWT", spacegroup="P1", dmin=dmin, alpha=alpha,
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


def parse_arguments():
    """Parse commandline arguments."""
    parser = argparse.ArgumentParser(
        description=(
            "Compute a real-space difference map between inputs in different space groups / crystal packings. "
            "You will need two MTZ files, which will be referred to throughout as 'on' and 'off', "
            "though they could also be light/dark, bound/apo, mutant/WT, hot/cold, etc. "
            "Each mtz will need to contain structure factor amplitudes and uncertainties; you will not need any phases. "
            "You will, however, need an input model (assumed to correspond with the 'off' state) which will be used to determine phases. "
            "The input file may be in .pdb or .cif format. "
            "Please note that phenix must be installed and active in your environment for this function to run. "
            ""
            "If your mtzoff and mtzon are in the same spacegroup and crystal packing, see the basic matchmaps utility. "
            "More information can be found online at https://rs-station.github.io/matchmaps/index.html"
        )
    )

    parser.add_argument(
        "--mtzoff",
        "-f",
        nargs=3,
        metavar=("mtzfileoff", "Foff", "SigFoff"),
        required=True,
        help=(
            "MTZ or sfCIF containing off/apo/ground/dark state data. "
            "Specified as [filename F SigF]"
        ),
    )

    parser.add_argument(
        "--mtzon",
        "-n",
        nargs=3,
        metavar=("mtzfileon", "Fon", "SigFon"),
        required=True,
        help=(
            "MTZ or SFCIF containing on/bound/excited/bright state data. "
            "Specified as [filename F SigF]"
            "This file may be in a different spacegroup / crystal packing than mtzoff"
        ),
    )

    parser.add_argument(
        "--pdboff",
        "-p",
        required=True,
        help=(
            "Reference PDB or mmCIF corresponding to the off/apo/ground/dark state. "
            "Used as a molecular replacement solution for mtzon and for rigid-body refinement of both input MTZs to generate phases."
            "Should match mtzoff well enough that molecular replacement is not necessary."
        ),
    )

    parser.add_argument(
        "--ligands",
        "-l",
        required=False,
        default=None,
        nargs="*",
        help=("Any .cif restraint files needed for refinement"),
    )

    parser.add_argument(
        "--input-dir",
        "-i",
        required=False,
        default="./",
        help="Path to input files. Optional, defaults to './' (current directory)",
    )

    parser.add_argument(
        "--output-dir",
        "-o",
        required=False,
        default="./",
        help="Path to which output files should be written. Optional, defaults to './' (current directory)",
    )

    parser.add_argument(
        "--on-as-stationary",
        required=False,
        action="store_true",
        default=False,
        help=(
            "Include this flag to align 'off' data onto 'on' data. By default, 'off' data is stationary and 'on' data is moved."
            "For matchmaps.mr, this only applies to the post-molecular-replacement alignment; "
            "all maps will be placed in the spacegroup of mtzoff."
        ),
    )
    
    parser.add_argument(
        "--no-bss",
        required=False,
        action="store_true",
        default=False,
        help=(
            "Include this flag to skip bulk solvent scaling in phenix.refine. By default, BSS is included."
        ),
    )

    parser.add_argument(
        "--spacing",
        "-s",
        required=False,
        type=float,
        default=0.5,
        help=(
            "Approximate voxel size in Angstroms for real-space maps. Defaults to 0.5 A. "
            "Value is approximate because there must be an integer number of voxels along each unit cell dimension"
        ),
    )

    parser.add_argument(
        "--dmin",
        required=False,
        type=float,
        default=None,
        help=(
            "Highest-resolution (in Angstroms) reflections to include in Fourier transform for FloatGrid creation. "
            "By default, cutoff is the resolution limit of the lower-resolution input MTZ. "
        ),
    )
    
    parser.add_argument(
        "--unmasked-radius",
        required=False,
        type=float,
        default=5,
        help=(
            "Maximum distance (in Anstroms) away from protein model to include voxels. Only applies to the 'unmasked' difference map output. "
            "Defaults to 5. "
            "Note that the regular difference map (e.g. the 'masked' version) is not affected by this parameter and maintains a solvent mask radius of 2 Angstroms."
        ),
    )

    parser.add_argument(
        "--alpha",
        required=False,
        type=float,
        default=0,
        help=(
            "Alpha to use for error weighting of F-obs prior to Fourier Transform. "
            "Weights are computed as: 1 / ((1+(alpha*(SigF^2)) / <SigF>^2). "
            "Default value is alpha=0, e.g., no weighting is performed. "
        )
    )

    parser.add_argument(
        "--verbose",
        "-v",
        required=False,
        action="store_true",
        default=False,
        help="Include this flag to print out phenix.phaser and phenix.refine outputs to the terminal. Useful for troubleshooting, but annoying; defaults to False.",
    )

    parser.add_argument(
        "--rbr-selections",
        "-r",
        required=False,
        default=None,
        nargs="*",
        help=(
            "Specification of multiple rigid-body groups for refinement. By default, everything is refined as one rigid-body group. "
            "For matchmaps.mr, everything will always be molecular replaced as a single rigid-body, but may then be refined as multiple rigid bodies."
        ),
    )

    parser.add_argument(
        "--eff",
        required=False,
        default=None,
        help=("Custom .eff template for running phenix.refine. "),
    )
    
    parser.add_argument(
        "--keep-temp-files",
        "-k",
        required=False,
        default=None,
        help=(
            "Do not delete intermediate matchmaps files, but rather place them in the supplied directory. "
            "This directory is created as a subdirectory of the supplied output-dir."
        )
    )

    parser.add_argument(
        "--script",
        required=False,
        default='run_matchmaps',
        help=(
            "Name for a file {script}.sh which can be run to repeat this command. "
            "By default, this file is called `run_matchmaps.sh`. "
            "Note that this file is written out in the current working directory, NOT the input or output directories"
        )
    )

    return parser


def main():
    parser = parse_arguments()
    args = parser.parse_args()

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
        no_bss = args.no_bss
    )
    
    if args.script:
        _write_script(
            utility = 'matchmaps.mr', 
            arguments = sys.argv[1:],
            script_name = args.script,
            )

    return


if __name__ == "__main__":
    main()
