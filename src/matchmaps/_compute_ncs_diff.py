"""Compute internal real space difference map across non-crystallographic symmetry."""

import argparse
import os
import sys
import time
from pathlib import Path

import gemmi
import reciprocalspaceship as rs

from matchmaps._parsers import matchmaps_ncs_parser
from matchmaps._phenix_utils import rigid_body_refinement_wrapper, _renumber_waters
from matchmaps._utils import (
    _handle_special_positions,
    make_floatgrid_from_mtz,
    _rbr_selection_parser,
    _ncs_align_and_subtract,
    _validate_environment,
    _validate_inputs,
    _clean_up_files,
    _cif_or_pdb_to_pdb,
    _cif_or_mtz_to_mtz,
    _write_script,
    _validate_column_dtypes,
)


def compute_ncs_difference_map(
    pdb: Path,
    mtz: Path,
    F: str,
    SigF: str = None,
    Phi: str = None,
    ligands: list = None,
    name: str = None,
    dmin: float = None,
    spacing=0.5,
    input_dir=Path("."),
    output_dir=Path("."),
    verbose=False,
    ncs_chains: list[str] = None,
    refine_ncs_separately=False,
    eff: str = None,
    keep_temp_files: str = None,
    no_bss=False,
    phenix_version: str = None,
):
    """
    Compute an internal difference map across non-crystallographic symmetry

    Parameters
    ----------
    pdb : Path
        Name of input .pdb file to use for phasing
    mtz : Path
        Name of input .mtz file
    F : str
        Column in mtz containing structure factor amplitudes
    SigF : str, optional
        Column in mtz containing structure factor uncertainties (optional)
    Phi : str, optional
        Column in mtz containing structure factor phases (optional)
    ligands : list, optional
        Any necessary .cif restrain files for refinement, by default None
    name : str, optional
        Name prefix for the output maps, by default None
    dmin : int, optional
        Resolution cutoff for input mtz, by default None
    spacing : float, optional
        Approximate size of real-space voxels in Angstroms, by default 0.5 A
    input_dir : pathlib.Path, optional
        Path to directory containing input files, by default Path(".") (current directory),
    output_dir : pathlib.Path, optional
        Path to directory containing input files, by default Path(".") (current directory),
    verbose : bool, optional
        If True, print outputs of scaleit and phenix.refine, by default False
    ncs_chains : list[str], optional
        Chains in input pdb comparable by ncs
    refine_ncs_separately : bool, optional
        If True, refine each ncs selection as a separate rigid body, by default False
    eff : str, optional
        Name of a file containing a template .eff parameter file for phenix.refine.
        If omitted, the sensible built-in .eff template is used. If you need to change something,
        I recommend copying the template from the source code and editing that.
    keep_temp_files : str, optional
        If not None, the name of a subdirectory of the output_dir into which intermediate matchmaps files are moved upon program completion.
    no_bss : bool, optional
        If True, skip bulk solvent scaling feature of phenix.refine
    phenix_version: str, optional
        Phenix version string to override the automatically detected version. I don't know why this would be necessary.
    """
    auto_phenix_version = _validate_environment(ccp4=False)

    if phenix_version:
        pass
    else:
        phenix_version = auto_phenix_version

    output_dir_contents = list(output_dir.glob("*"))

    pdb = _cif_or_pdb_to_pdb(pdb, output_dir)

    mtz, _ = _cif_or_mtz_to_mtz(mtz, output_dir)

    rbr_phenix, rbr_gemmi = _rbr_selection_parser(ncs_chains)

    if Phi is None:  # do rigid-body refinement to get phases
        if not refine_ncs_separately:
            rbr_phenix = None

        pdb = _handle_special_positions(pdb, output_dir)

        pdb = _renumber_waters(pdb, verbose)

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
            no_bss=no_bss,
            phenix_style=phenix_version,
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

    _clean_up_files(output_dir, output_dir_contents, keep_temp_files)

    print(f"{time.strftime('%H:%M:%S')}: Done!")

    return


def main():
    args = matchmaps_ncs_parser.parse_args()

    (input_dir, output_dir, ligands, mtz, pdb) = _validate_inputs(
        args.input_dir,
        args.output_dir,
        args.ligands,
        args.mtz[0],
        args.pdb,
    )

    if args.phases:
        columns = (args.mtz[1], args.phases)
        dtypes = (rs.StructureFactorAmplitudeDtype, rs.PhaseDtype)
    else:
        columns = (args.mtz[1], args.mtz[2])
        dtypes = (rs.StructureFactorAmplitudeDtype, rs.StandardDeviationDtype)

    _validate_column_dtypes(rs.read_mtz(str(mtz)), columns, dtypes)

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
        no_bss=args.no_bss,
        phenix_version=args.phenix_version,
    )

    if args.script:
        _write_script(
            utility="matchmaps.ncs",
            arguments=sys.argv[1:],
            script_name=args.script,
            phenix_version=args.phenix_version,
        )

    return


if __name__ == "__main__":
    main()
