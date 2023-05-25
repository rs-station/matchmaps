"""Compute unbiased real space difference map from inputs in different spacegroups"""

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
    _remove_waters,
    _restore_ligand_occupancy,
    phaser_wrapper,
)


def compute_mr_difference_map(
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
    rbr_selections=None,
    eff=None,
):
    
    off_name = str(mtzoff.removesuffix(".mtz"))
    on_name = str(mtzon.removesuffix(".mtz"))
    
    # make sure directories have a trailing slash!
    if input_dir[-1] != "/":
        input_dir = input_dir + "/"

    if output_dir[-1] != "/":
        output_dir = output_dir + "/"

    # take in the list of rbr selections and parse them into phenix and gemmi selection formats
    # if rbr_groups = None, just returns (None, None)
    rbr_phenix, rbr_gemmi = _rbr_selection_parser(rbr_selections)

    # this is where scaling takes place in the usual pipeline, but that doesn't make sense with different-spacegroup inputs 
    # side note: I need to test the importance of scaling even in the normal case!! Might be more artifact than good, who knows
    
    pdboff = _handle_special_positions(pdboff, input_dir, output_dir)

    # write this function as a wrapper around phenix.pdbtools
    # modified pdboff already moved to output_dir by _handle_special_positions
    pdboff = _remove_waters(pdboff, output_dir)
    
    print(f"{time.strftime('%H:%M:%S')}: Running phenix.phaser to place 'off' model into 'on' data...")

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
    
    
    # mtzon = phaser_nickname + ".1.mtz" # jokes, original mtzon is fine!
    pdb_mr_to_on = phaser_nickname + ".1.pdb"
    
    # the refinement process *should* be identical. Waters are gone already
    # I just need to make sure that the phaser outputs go together
    print(f"{time.strftime('%H:%M:%S')}: Running phenix.refine for the 'on' data...")

    nickname_on = rigid_body_refinement_wrapper(
        mtzon=mtzon,
        pdboff=pdb_mr_to_on,
        input_dir=input_dir,
        output_dir=output_dir,
        ligands=ligands,
        eff=eff,
        verbose=verbose,
        rbr_selections=rbr_phenix,
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
    )
    
    # from here down I just copied over the stuff from the normal version
    # this should be proofread for compatibility but should all work
    
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
            )
    # print(f"{time.strftime('%H:%M:%S')}: Cleaning up files...")
        
    # _clean_up_files()

    print(f"{time.strftime('%H:%M:%S')}: Done!")

    return