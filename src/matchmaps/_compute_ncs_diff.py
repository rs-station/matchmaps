"""Compute internal real space difference map across non-crystallographic symmetry."""

import argparse
import os
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
    _renumber_waters,
    _ncs_align_and_subtract,
    _validate_environment,
    _validate_inputs,
    _clean_up_files,
    _cif_or_pdb_to_pdb,
)


def compute_ncs_difference_map(
    pdb,
    mtz,
    F,
    SigF=None,
    Phi=None,
    ligands=None,
    name=None,
    dmin=None,
    spacing=0.5,
    input_dir=Path("."),
    output_dir=Path("."),
    verbose=False,
    ncs_chains=None,
    refine_ncs_separately=False,
    eff=None,
    keep_temp_files=None
):
    _validate_environment(ccp4=False)
    
    output_dir_contents = list(output_dir.glob("*"))
    
    pdb = _cif_or_pdb_to_pdb(pdb, output_dir)
    
    rbr_phenix, rbr_gemmi = _rbr_selection_parser(ncs_chains)

    if Phi is None:  # do rigid-body refinement to get phases
        if not refine_ncs_separately:
            rbr_phenix = None

        pdb = _handle_special_positions(pdb, output_dir)

        pdb = _renumber_waters(pdb)

        print(f"{time.strftime('%H:%M:%S')}: Running phenix.refine...")

        nickname = rigid_body_refinement_wrapper(
            mtzon=mtz,
            pdboff=pdb,
            input_dir=input_dir,
            output_dir=output_dir,
            ligands=ligands,
            eff=eff,
            verbose=verbose,
            rbr_selections=rbr_phenix,
            off_labels=f"{F},{SigF}",
        )

        # use phenix names for columns when computing FloatGrid
        fname = "F-obs-filtered"
        phiname = "PH2FOFCWT"

        mtzfilename = f"{nickname}_1.mtz"
        pdbfilename = f"{nickname}_1.pdb"

    else:
        print(f"{time.strftime('%H:%M:%S')}: Using provided phases...")
        fname = F
        phiname = Phi

        mtzfilename = str(mtz)
        pdbfilename = str(pdb)

    # regardless of whether refinement was performed, read in these two files:
    mtz = rs.read_mtz(mtzfilename)
    pdb = gemmi.read_structure(pdbfilename)

    print(f"{time.strftime('%H:%M:%S')}: Constructing FloatGrid from mtz...")
    # hard-coding F, Phi because they're always phenix outputs
    # TO-DO: Figure out why phenix outputs are sometimes still split into (+) and (-) columns, even when I specify that anomalous=False
    # As a workaround, even anomalous files have a single 'F-obs-filtered' column, so I can always just use that.
    fg = make_floatgrid_from_mtz(mtz, spacing, F=fname, Phi=phiname, spacegroup="P1")

    _ncs_align_and_subtract(
        fg=fg,
        pdb=pdb,
        ncs_chains=ncs_chains,
        output_dir=output_dir,
        name=name,
    )
    
    print(f"{time.strftime('%H:%M:%S')}: Cleaning up files...")
    print(keep_temp_files)
    _clean_up_files(output_dir, output_dir_contents, keep_temp_files)

    print(f"{time.strftime('%H:%M:%S')}: Done!")

    return

def parse_arguments():
    """Parse commandline arguments."""
    parser = argparse.ArgumentParser(
        description=(
            "Compute an 'internal' real-space difference map between NCS-related molecules. "
            "You will need an MTZ file with structure factor amplitudes and optionally containing phases, and a PDB/CIF file."
            ""
            "Please note that phenix must be installed and active in your environment for this function to run. "
            ""
            "More information can be found online at https://rs-station.github.io/matchmaps/index.html"
        )
    )

    parser.add_argument(
        "--mtz",
        "-m",
        nargs="*",
        # metavar=("mtzfile", "F", "SigF"),
        required=True,
        help=(
            "MTZ file containing structure factor amplitudes. "
            "Specified as [filename F SigF] or [filename F]. "
            "SigF is not necessary if phases are also provided"
        ),
    )

    parser.add_argument(
        "--phases",
        required=False,
        default=None,
        help=(
            "Optional. Column in MTZ file containing phases. "
            "If phases are not provided, phases will be computed via rigid-body refinement of "
            "the provided model and structure factor amplitudes."
        ),
    )

    parser.add_argument(
        "--pdb",
        "-p",
        required=True,
        help=(
            "Reference pdb/cif. "
            "If phases are not provided, used for rigid-body refinement of input MTZ to generate phases."
        ),
    )

    parser.add_argument(
        "--ncs-chains",
        "-n",
        required=True,
        metavar=("fixed_chain", "moving_chain"),
        default=None,
        nargs=2,
        help=(
            "NCS chains to overlay and subtract, specified as [fixed_chain, moving_chain]."
            "E.g. to overlay chain C onto chain B, specify: --ncs-chains B C"
        ),
    )

    parser.add_argument(
        "--mapname",
        required=False,
        default="matchmaps_ncs",
        help=("Base filename for the output map files. "),
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
        help="Path to input mtz and pdb. Optional, defaults to './' (current directory)",
    )

    parser.add_argument(
        "--output-dir",
        "-o",
        required=False,
        default="./",
        help="Path to which output files should be written. Optional, defaults to './' (current directory)",
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
            "By default, cutoff is the resolution limit of the input MTZ. "
        ),
    )

    parser.add_argument(
        "--verbose",
        "-v",
        required=False,
        action="store_true",
        default=False,
        help="Include this flag to print out phenix.refine outputs to the terminal. Useful for troubleshooting, but annoying; defaults to False.",
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
    
    return parser


def main():
    parser = parse_arguments()
    args = parser.parse_args()

    (input_dir, output_dir, ligands, mtz, pdb) = _validate_inputs(
        args.input_dir,
        args.output_dir,
        args.ligands,
        args.mtz[0],
        args.pdb,
    )

    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    if not os.path.exists(args.input_dir):
        raise ValueError(f"Input directory '{args.input_dir}' does not exist")

    compute_ncs_difference_map(
        pdb=pdb,
        mtz=mtz,
        F=args.mtz[1],
        SigF=args.mtz[2] if (len(args.mtz) == 3) else None,
        Phi=args.phases,
        name=args.mapname,
        ligands=ligands,
        input_dir=input_dir,
        output_dir=output_dir,
        verbose=args.verbose,
        ncs_chains=args.ncs_chains,
        eff=args.eff,
        dmin=args.dmin,
        spacing=args.spacing,
        keep_temp_files=args.keep_temp_files,
    )

    return


if __name__ == "__main__":
    main()
