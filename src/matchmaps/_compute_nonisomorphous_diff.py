"""Compute unbiased nonisomorphous difference map."""

import argparse
import os
import subprocess
import time
from functools import partial

import gemmi
import numpy as np
import reciprocalspaceship as rs

from matchmaps._utils import (
    _handle_special_positions,
    align_grids_from_model_transform,
    make_floatgrid_from_mtz,
    rigid_body_refinement_wrapper,
)


def compute_nonisomorphous_difference_map(
    pdboff,
    mtzoff,
    mtzon,
    Foff,
    SigFoff,
    Fon,
    SigFon,
    ligands=None,
    dmin=None,
    spacing=0.5,
    on_as_stationary=False,
    input_dir="./",
    output_dir="./",
    verbose=False,
    selection=None,
    eff=None,
):
    """
    _summary_.

    Parameters
    ----------
    pdboff : string
        Name of input .pdb file to use for phasing
    mtzoff : string
        Name of input .mtz file containing 'off' data
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
        _description_, by default False
    input_dir : str, optional
        Path to directory containing input files, by default "./" (current directory)
    output_dir : str, optional
        Path to directory to which output files should be written, by default "./" (current directory)
    verbose : bool, optional
        If True, print outputs of scaleit and phenix.refine, by default False
    selection : str, optional
        Custom selection to provide to refinement.refine.sites.rigid_body=
        If omitted, then refinement.refine.sites.rigid_body=all
        Custom selection is only necessary if special-position atoms need to be excluded from refinement
    eff : str, optional
        Name of a file containing a template .eff parameter file for phenix.refine.
        If omitted, the sensible built-in .eff template is used. If you need to change something,
        I recommend copying the template from the source code and editing that.
    """
    off_name = str(mtzoff.removesuffix(".mtz"))
    on_name = str(mtzon.removesuffix(".mtz"))

    # make sure directories have a trailing slash!
    if input_dir[-1] != "/":
        input_dir = input_dir + "/"

    if output_dir[-1] != "/":
        output_dir = output_dir + "/"

    mtzon_scaled = mtzon.removesuffix(".mtz") + "_scaled" + ".mtz"

    print(
        f"{time.strftime('%H:%M:%S')}: Running scaleit to scale 'on' data to 'off' data..."
    )

    subprocess.run(
        f"rs.scaleit -r {input_dir}/{mtzoff} {Foff} {SigFoff} -i {input_dir}/{mtzon} {Fon} {SigFon} -o {output_dir}/{mtzon_scaled}",
        shell=True,
        capture_output=(not verbose),
    )

    pdboff = _handle_special_positions(pdboff, input_dir, output_dir)

    mtzon = mtzon_scaled

    print(f"{time.strftime('%H:%M:%S')}: Running phenix.refine for the 'on' data...")

    nickname_on = rigid_body_refinement_wrapper(
        mtzon=mtzon,
        pdboff=pdboff,
        input_dir=input_dir,
        output_dir=output_dir,
        ligands=ligands,
        eff=eff,
        verbose=verbose,
        selection=selection,
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
        selection=selection,
        off_labels=f"{Foff},{SigFoff}",
    )

    # done with reciprocal space; transitioning to real space

    # read back in the files created by phenix
    # these have knowable names
    mtzon = rs.read_mtz(f"{output_dir}/{nickname_on}_1.mtz")
    mtzoff = rs.read_mtz(f"{output_dir}/{nickname_off}_1.mtz")

    pdbon = gemmi.read_structure(f"{output_dir}/{nickname_on}_1.pdb")
    pdboff = gemmi.read_structure(f"{output_dir}/{nickname_off}_1.pdb")

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
        mtzoff, spacing, F="F-obs-filtered", Phi="PH2FOFCWT", spacegroup="P1", dmin=dmin
    )
    fg_on = make_floatgrid_from_mtz(
        mtzon, spacing, F="F-obs-filtered", Phi="PH2FOFCWT", spacegroup="P1", dmin=dmin
    )

    # TO-DO: think more about scaling!

    print(f"{time.strftime('%H:%M:%S')}: Using models to rigid-body align maps...")
    if on_as_stationary:
        # rs.io.write_ccp4_map(fg_off.array, f'{output_dir}/off_before_transforming.map',
        #                      fg_off.unit_cell, fg_off.spacegroup)
        fg_off = align_grids_from_model_transform(fg_on, fg_off, pdbon, pdboff)
        fg_on = align_grids_from_model_transform(fg_on, fg_on, pdbon, pdbon)
        pdb = pdbon
        fg_ref = fg_on
    else:
        # rs.io.write_ccp4_map(fg_on.array, f'{output_dir}/on_before_transforming.map',
        #                      fg_on.unit_cell, fg_on.spacegroup)
        fg_on = align_grids_from_model_transform(fg_off, fg_on, pdboff, pdbon)
        fg_off = align_grids_from_model_transform(
            fg_off, fg_off, pdboff, pdboff
        )  # apply same masking sitch to both grids?
        pdb = pdboff
        fg_ref = fg_off

    print(f"{fg_off.array.mean()=}, {fg_on.array.mean()=}")

    # do this again, because transformation + carving can mess up scales:
    fg_on.normalize()
    fg_off.normalize()

    print(f"{fg_off.array.mean()=}, {fg_on.array.mean()=}")

    print(f"{time.strftime('%H:%M:%S')}: Writing files...")

    difference_array = fg_on.array - fg_off.array

    # all that's left is to mask out voxels that aren't near the model!
    # we can do this in gemmi
    fg_mask_only = fg_ref.clone()
    masker = gemmi.SolventMasker(gemmi.AtomicRadiiSet.Cctbx)
    masker.rprobe = 2  # this should do it, idk, no need to make this a user parameter

    masker.put_mask_on_float_grid(fg_mask_only, pdb[0])
    masked_difference_array = np.logical_not(fg_mask_only.array) * difference_array

    # and finally, write stuff out
    # use partial function to guarantee I'm always using the same and correct cell

    # coot refuses to render periodic boundaries for P1 maps with alpha=beta=gamma=90, sooooo
    if all(
        [
            angle == 90
            for angle in (
                fg_ref.unit_cell.alpha,
                fg_ref.unit_cell.beta,
                fg_ref.unit_cell.gamma,
            )
        ]
    ):
        fg_ref.unit_cell = gemmi.UnitCell(
            fg_ref.unit_cell.a,
            fg_ref.unit_cell.b,
            fg_ref.unit_cell.c,
            90.006,
            fg_ref.unit_cell.beta,
            fg_ref.unit_cell.gamma,
        )
        print("did silly angle thing")

    write_maps = partial(
        rs.io.write_ccp4_map, cell=fg_ref.unit_cell, spacegroup=fg_ref.spacegroup
    )

    write_maps(fg_on.array, f"{output_dir}/{on_name}.map")

    write_maps(fg_off.array, f"{output_dir}/{off_name}.map")

    write_maps(masked_difference_array, f"{output_dir}/{on_name}_minus_{off_name}.map")
    write_maps(
        difference_array, f"{output_dir}/{on_name}_minus_{off_name}_unmasked.map"
    )

    print(f"{time.strftime('%H:%M:%S')}: Done!")

    return


def parse_arguments():
    """Parse commandline arguments."""
    parser = argparse.ArgumentParser(
        description=(
            "Compute a difference map using non-isomorphous inputs. "
            "You will need two MTZ files, which will be referred to throughout as 'on' and 'off', "
            "though they could also be light/dark, bound/apo, mutant/WT, hot/cold, etc. "
            "Each mtz will need to contain structure factor amplitudes and uncertainties; you will not need any phases. "
            "You will, however, need an input model (assumed to correspond with the 'off' state) which will be used to determine phases. "
            "Please note that both ccp4 and phenix must be installed and active in your environment for this function to work. "
        )
    )

    parser.add_argument(
        "--mtzoff",
        "-f",
        nargs=3,
        metavar=("mtzfileoff", "Foff", "SigFoff"),
        required=True,
        help=(
            "MTZ containing off/apo/ground/dark state data. "
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
            "MTZ containing on/bound/excited/bright state data. "
            "Specified as [filename F SigF]"
        ),
    )

    parser.add_argument(
        "--pdboff",
        "-p",
        required=True,
        help=(
            "Reference pdb corresponding to the off/apo/ground/dark state. "
            "Used for rigid-body refinement of both input MTZs to generate phases."
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
        help="Path to input mtzs and pdb. Optional, defaults to './' (current directory)",
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
            "By default, cutoff is the resolution of the lower-resolution input MTZ. "
        ),
    )

    parser.add_argument(
        "--verbose",
        "-v",
        required=False,
        action="store_true",
        default=False,
        help="Include this flag to print out scaleit and phenix.refine outputs to the terminal. Useful for troubleshooting, but annoying; defaults to False.",
    )

    parser.add_argument(
        "--selection",
        required=False,
        default=None,
        help=(
            "Argument to phenix.refine's refinement.refine.sites.rigid_body. "
            "Use this flag if atoms on special positions need to be omitted from the atom selection. "
        ),
    )

    parser.add_argument(
        "--eff",
        required=False,
        default=None,
        help=("Custom .eff template for running phenix.refine. "),
    )

    return parser


def main():
    parser = parse_arguments()
    args = parser.parse_args()

    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    compute_nonisomorphous_difference_map(
        pdboff=args.pdboff,
        ligands=args.ligands,
        mtzoff=args.mtzoff[0],
        mtzon=args.mtzon[0],
        Foff=args.mtzoff[1],
        SigFoff=args.mtzoff[2],
        Fon=args.mtzon[1],
        SigFon=args.mtzon[2],
        input_dir=args.input_dir,
        output_dir=args.output_dir,
        verbose=args.verbose,
        selection=args.selection,
        eff=args.eff,
        dmin=args.dmin,
        spacing=args.spacing,
        on_as_stationary=args.on_as_stationary,
    )

    return


if __name__ == "__main__":
    main()
