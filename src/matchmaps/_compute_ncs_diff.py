"""Compute internal real space difference map across non-crystallographic symmetry."""

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
    # align_grids_from_model_transform,
    make_floatgrid_from_mtz,
    rigid_body_refinement_wrapper,
    _realspace_align_and_subtract,
    _rbr_selection_parser,
    _renumber_waters,
)


def compute_ncs_difference_map(
    pdb,
    mtz,
    F,
    SigF,
    ligands=None,
    dmin=None,
    spacing=0.5,
    input_dir="./",
    output_dir="./",
    verbose=False,
    ncs_selections=None,
    refine_ncs_separately=False,
    eff=None,
):
    mtz_name = str(mtz.removesuffix(".mtz"))
    
    # make sure directories have a trailing slash!
    if input_dir[-1] != "/":
        input_dir = input_dir + "/"

    if output_dir[-1] != "/":
        output_dir = output_dir + "/"
        
    rbr_phenix, rbr_gemmi = _rbr_selection_parser(ncs_selections)
    
    if not refine_ncs_separately:
        rbr_phenix = None

    pdb = _handle_special_positions(pdb, input_dir, output_dir)

    pdb = _renumber_waters(pdb, output_dir)

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

    mtz = rs.read_mtz(f"{output_dir}/{nickname}_1.mtz")

    pdb = gemmi.read_structure(f"{output_dir}/{nickname}_1.pdb")

    if dmin is None:
        dmin = min(mtz.compute_dHKL(inplace=True).dHKL),
        
    print(f"{time.strftime('%H:%M:%S')}: Constructing FloatGrids from mtzs...")
    # hard-coding F, Phi because they're always phenix outputs
    # TO-DO: Figure out why phenix outputs are sometimes still split into (+) and (-) columns, even when I specify that anomalous=False
    # As a workaround, even anomalous files have a single 'F-obs-filtered' column, so I can always just use that.
    fg = make_floatgrid_from_mtz(
        mtz, spacing, F="F-obs-filtered", Phi="PH2FOFCWT", spacegroup="P1", dmin=dmin
    )
    
    _ncs_align_and_subtract(
       # this is the key utility! write this!
    )

def parse_arguments():
    
    parser = argparse.ArgumentParser()
    
    return parser

def main():
    parser = parse_arguments()
    args = parser.parse_args()

    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    compute_ncs_difference_map(
        
    )
    
    return

if __name__ == "__main__":
    main()