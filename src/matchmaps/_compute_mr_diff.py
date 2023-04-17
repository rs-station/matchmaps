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
    _renumber_waters,
)


def compute_realspace_difference_map(
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

    