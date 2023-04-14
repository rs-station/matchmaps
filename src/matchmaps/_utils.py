"""
Utilities and helperfunctions used by matchmaps.

Functions make_floatgrid_from_mtz, rigid_body_refinement_wrapper, and align_grids_from_model_transform
are exported to python for use in prototyping and testing
"""

import glob
import shutil
import subprocess
import time
import re
from functools import partial

import gemmi
import numpy as np
import reciprocalspaceship as rs


def _rbr_selection_parser(rbr_selections):
    # end early and return nones if this feature isn't being used
    if rbr_selections is None:
        return None, None

    # validate input is a list of strings
    if not isinstance(rbr_selections, list):
        raise ValueError(
            f"rbr_selections must be a list of str, not {type(rbr_selections).__name__}"
        )
    if not all([isinstance(sel, str) for sel in rbr_selections]):
        raise ValueError(
            f"rbr_selections must be a list of str, not list of {[type(sel).__name__ for sel in rbr_selections]}"
        )

    # get rid of spaces, which have no semantic meaning to my parser
    rbr_selections = [g.replace(" ", "") for g in rbr_selections]
    phenix_selections = []
    gemmi_selections = []

    for group in rbr_selections:
        # check if there are multiple, separate selections
        # that make up this rigid-body selection
        if len(group.split(",")) > 1:
            subgroups = group.split(",")
            phegem = [_subparser(sub) for sub in subgroups]

            phenixtemp = [pg[0] for pg in phegem]

            phe = " or ".join([f"({p})" for p in phenixtemp])

            # gem = [pg[1] for pg in phegem]
            gem = phegem[0][1]

            phenix_selections.append(phe)
            gemmi_selections.append(gem)

        # just one selection in each
        else:
            phe, gem = _subparser(group)
            phenix_selections.append(phe)
            gemmi_selections.append(gem)

    return phenix_selections, gemmi_selections


def _subparser(selection):
    if selection.isalpha():  # just a chain name
        phe = "chain " + selection
        gem = selection

    else:
        chain = "".join(re.findall("[a-zA-Z]", selection))
        residues = selection.removeprefix(chain)

        if "-" in residues:
            start, end = [int(r) for r in residues.split("-")]

            phe = f"chain {chain} and resseq {start}:{end}"
            gem = f"{chain}/{start}-{end}"

        else:
            phe = f"chain {chain} and resid {residues}"
            gem = f"{chain}/{residues}"

    return phe, gem


def make_floatgrid_from_mtz(mtz, spacing, F, Phi, spacegroup="P1", dmin=None):
    """
    Make a gemmi.FloatGrid from an rs.DataSet.

    Parameters
    ----------
    mtz : rs.DataSet
        mtz data to be transformed into real space
    spacing : float
        Approximate voxel size desired (will be rounded as necessary to create integer grid dimensions)
    F : str, optional
        Column in mtz containing structure factor amplitudes to use for calculation
    Phi : str, optional
        Column in mtz containing phases to be used for calculation
    spacegroup : str, optional
        Spacegroup for the output FloatGrid. Defaults to P1.
    dmin: float, optional
        Highest resolution reflections to include in Fourier transform. Defaults to None, no cutoff.

    Returns
    -------
    float_grid : gemmi.FloatGrid
        Fourier transform of mtz, written out as a gemmi object containing a 3D voxel array
        and various other metadata and methods
    """
    # drop NAs in either of the specified columns
    # this has the secondary purpose of not silently modifying the input mtz
    new_mtz = mtz[~mtz[F].isnull()]
    new_mtz = new_mtz[~new_mtz[Phi].isnull()]

    new_mtz.hkl_to_asu(inplace=True)

    # apply user-provided resolution cutoff:
    if dmin is not None:
        new_mtz = new_mtz.compute_dHKL(inplace=True).loc[new_mtz.dHKL > dmin]

    # compute desired grid size based on given spacing
    gridsize = [
        int(dim // spacing) for dim in (new_mtz.cell.a, new_mtz.cell.b, new_mtz.cell.c)
    ]

    # perform FFT using the desired amplitudes and phases
    new_mtz["Fcomplex"] = new_mtz.to_structurefactor(F, Phi)
    reciprocal_grid = new_mtz.to_reciprocal_grid("Fcomplex", grid_size=gridsize)
    real_grid = np.real(np.fft.fftn(reciprocal_grid))

    # declare gemmi.FloatGrid object
    float_grid = gemmi.FloatGrid(*gridsize)
    float_grid.set_unit_cell(new_mtz.cell)

    if spacegroup is not None:
        float_grid.spacegroup = gemmi.find_spacegroup_by_name(spacegroup)
    else:
        float_grid.spacegroup = new_mtz.spacegroup

    # write real_grid into float_grid via buffer protocol
    temp = np.array(float_grid, copy=False)
    temp[:, :, :] = real_grid[:, :, :]

    # Enforce that mean=0, stdev=1 for voxel intensities
    float_grid.normalize()

    return float_grid


def rigid_body_refinement_wrapper(
    mtzon,
    pdboff,
    input_dir,
    output_dir,
    off_labels=None,
    ligands=None,
    eff=None,
    verbose=False,
    rbr_selections=None,
):
    # confirm that phenix is active in the command-line environment
    if shutil.which("phenix.refine") is None:
        raise OSError(
            "Cannot find executable, phenix.refine. Please set up your phenix environment."
        )

    if eff is None:
        eff_contents = """
refinement {
  crystal_symmetry {
    unit_cell = cell_parameters
    space_group = sg
  }
  input {
    pdb {
      file_name = pdb_input
    }
    xray_data {
      file_name = "mtz_input"
      labels = columns
      r_free_flags {
        generate=True
      }
      force_anomalous_flag_to_be_equal_to = False
    }
    monomers {
      ligands
    }
  }
  output {
    prefix = '''nickname'''
    serial = 1
    serial_format = "%d"
    job_title = '''nickname'''
    write_def_file = False
    write_eff_file = False
    write_geo_file = False
  }
  electron_density_maps {
    map_coefficients {
      map_type = "2mFo-DFc"
      mtz_label_amplitudes = "2FOFCWT"
      mtz_label_phases = "PH2FOFCWT"
    }
  }
  refine {
    strategy = *rigid_body
    sites {
      rigid_body_sites
    }
  }
  main {
    number_of_macro_cycles = 1
    nproc = 8
  }
}
    """
    else:
        with open(input_dir + eff) as file:
            eff_contents = file.read()

    if off_labels is None:
        nickname = f"{mtzon.removesuffix('.mtz')}_rbr_to_{pdboff.removesuffix('.pdb')}"
    else:
        nickname = f"{mtzon.removesuffix('.mtz')}_rbr_to_self"

    # check existing files because phenix doesn't like to overwrite things
    similar_files = glob.glob(f"{output_dir}/{nickname}_[0-9]_1.*")
    if len(similar_files) == 0:
        nickname += "_0"
    else:
        n = max([int(s.split("_")[-2]) for s in similar_files])
        nickname += f"_{n+1}"

    # read in mtz to access cell parameters and spacegroup
    mtz = rs.read_mtz((output_dir if (off_labels is None) else input_dir) + mtzon)
    cell_string = f"{mtz.cell.a} {mtz.cell.b} {mtz.cell.c} {mtz.cell.alpha} {mtz.cell.beta} {mtz.cell.gamma}"
    sg = mtz.spacegroup.short_name()

    # edit refinement template
    eff = f"{output_dir}/params_{nickname}.eff"

    params = {
        "sg": sg,
        "cell_parameters": cell_string,
        "pdb_input": output_dir + pdboff,
        "mtz_input": (output_dir if (off_labels is None) else input_dir) + mtzon,
        "nickname": output_dir + nickname,
    }

    if off_labels is None:
        params["columns"] = "FPH1,SIGFPH1"  # names from scaleit output
    else:
        params["columns"] = off_labels

    # if selection is not None:
    #     params["all"] = selection  # overwrite atom selection

    for key, value in params.items():
        eff_contents = eff_contents.replace(key, value)

    # either add ligands to .eff file or delete "ligands" placeholder
    if ligands is not None:
        ligand_string = "\n".join([f"file_name = '{input_dir}/{l}'" for l in ligands])
        eff_contents = eff_contents.replace("ligands", ligand_string)
    else:
        eff_contents = eff_contents.replace("ligands", "")

    if rbr_selections is not None:
        selection_string = "\n".join(
            [f"rigid_body = '{sel}'" for sel in rbr_selections]
        )
        eff_contents = eff_contents.replace("rigid_body_sites", selection_string)
    else:
        eff_contents = eff_contents.replace("rigid_body_sites", "rigid_body = all")

    # write out customized .eff file for use by phenix
    with open(eff, "w") as file:
        file.write(eff_contents)

    # run refinement!
    # print refinement output to terminal if user supplied the --verbose flag
    subprocess.run(
        f"phenix.refine {eff}",
        shell=True,
        capture_output=(not verbose),
    )

    return nickname


def _handle_special_positions(pdboff, input_dir, output_dir):
    """
    Check if any waters happen to sit on special positions, and if so, remove them.
    If any non-water atoms sit on special positions, throw a (hopefully helpful) error.

    Regardless of whether any special positions were found, copy this file over to output_dir

    Parameters
    ----------
    pdboff : str
        name of input pdb
    path : str
        Relative in/out path for files
    """
    pdb = gemmi.read_structure(input_dir + pdboff)

    # gemmi pdb heirarchy is nonsense, but this is what we've got
    for model in pdb:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    # check if atom is within 0.5 A of a special position
                    # returns a number > 0 if so
                    if pdb.cell.is_special_position(atom.pos, max_dist=0.5) > 0:
                        if residue.is_water():
                            print(
                                "   Input model contains water at special position. Removing water so as to not break rigid-body refinement."
                            )
                            print(
                                # "   If it is important that you keep this water and just omit it from your refinement selection, you many the --rbr-selections flag"
                            )
                            residue.remove_atom(
                                name=atom.name, altloc=atom.altloc, el=atom.element
                            )

                        else:
                            raise NotImplementedError(
                                """
Input model contains a non-water atom on a special position, which is not allowed in rigid-body refinement.
You may attempt to exclude this atom via a custom atom selection to the `--rbr-selections` argument, but this feature is not currently tested.

Alternatively, you can remove this atom from your structure altogether and try again.
"""
                            )

    pdboff_nospecialpositions = pdboff.removesuffix(".pdb") + "_nospecialpositions.pdb"

    pdb.write_pdb(output_dir + pdboff_nospecialpositions)

    return pdboff_nospecialpositions


def _renumber_waters(pdb, dir):
    """
    Call phenix.sort_hetatms to place waters onto the nearest protein chain. This ensures that rbr selections handle waters properly

    Parameters
    ----------
    pdb : str
        name of pdb file
    dir : str
        directory in which pdb file lives
    """

    pdb_renumbered = pdb.removesuffix(".pdb") + "_renumbered.pdb"

    subprocess.run(
        f"phenix.sort_hetatms file_name={dir}/{pdb} output_file={dir}/{pdb_renumbered}",
        shell=True,
        capture_output=True,
    )

    print(f"{time.strftime('%H:%M:%S')}: Moved waters to nearest protein chains...")

    return pdb_renumbered


def _realspace_align_and_subtract(
    output_dir,
    fg_off,
    fg_on,
    pdboff,
    pdbon,
    on_name,
    off_name,
    on_as_stationary,
    selection,
):
    """
    Take in two floatgrids and two pdbs. Apply the alignment transformation and subtract.

    If there were multiple selections used for rigid-body refinement, this method should be called separately for each selection.

    Parameters
    ----------
    output_dir : _type_
        _description_
    fg_off : _type_
        _description_
    fg_on : _type_
        _description_
    pdboff : _type_
        _description_
    pdbon : _type_
        _description_
    on_name : _type_
        _description_
    off_name : _type_
        _description_
    on_as_stationary : _type_
        _description_
    selection : str or list of str
        If not None, should be a valid gemmi selection string or a list of valid gemmi selection strings
    """

    if selection:
        print(
            f"{time.strftime('%H:%M:%S')}: Using models to rigid-body align maps for rigid-body selection {selection}..."
        )
    else:
        print(f"{time.strftime('%H:%M:%S')}: Using models to rigid-body align maps...")

    rs.io.write_ccp4_map(
        fg_on.array,
        output_dir + on_name + "_before.map",
        cell=fg_on.unit_cell,
        spacegroup=fg_on.spacegroup,
    )
    rs.io.write_ccp4_map(
        fg_off.array,
        output_dir + off_name + "_before.map",
        cell=fg_off.unit_cell,
        spacegroup=fg_off.spacegroup,
    )

    if on_as_stationary:
        fg_fixed = fg_on.clone()
        pdb_fixed = pdbon.clone()
    else:
        fg_fixed = fg_off.clone()
        pdb_fixed = pdboff.clone()

    fg_off = align_grids_from_model_transform(
        fg_fixed, fg_off, pdb_fixed, pdboff, selection
    )
    fg_on = align_grids_from_model_transform(
        fg_fixed, fg_on, pdb_fixed, pdbon, selection
    )

    # do this again, because transformation + carving can mess up scales:
    fg_on.normalize()
    fg_off.normalize()

    print(f"{time.strftime('%H:%M:%S')}: Writing files...")

    difference_array = fg_on.array - fg_off.array

    # all that's left is to mask out voxels that aren't near the model!
    # we can do this in gemmi
    fg_mask_only = fg_fixed.clone()
    masker = gemmi.SolventMasker(gemmi.AtomicRadiiSet.Cctbx)
    masker.rprobe = 2  # this should do it, idk, no need to make this a user parameter

    # if selection is None:
    if selection is None:
        pdb_for_mask = pdb_fixed[0]
    else:
        pdb_for_mask = _extract_pdb_selection(pdb_fixed[0], selection)

    masker.put_mask_on_float_grid(fg_mask_only, pdb_for_mask)
    masked_difference_array = np.logical_not(fg_mask_only.array) * difference_array

    # and finally, write stuff out

    # coot refuses to render periodic boundaries for P1 maps with alpha=beta=gamma=90, sooooo
    fg_fixed.unit_cell = _unit_cell_hack(fg_fixed.unit_cell)

    # use partial function to guarantee I'm always using the same and correct cell
    write_maps = partial(
        rs.io.write_ccp4_map, cell=fg_fixed.unit_cell, spacegroup=fg_fixed.spacegroup
    )

    write_maps(fg_on.array, f"{output_dir}/{on_name}.map")

    write_maps(fg_off.array, f"{output_dir}/{off_name}.map")

    write_maps(masked_difference_array, f"{output_dir}/{on_name}_minus_{off_name}.map")
    write_maps(
        difference_array, f"{output_dir}/{on_name}_minus_{off_name}_unmasked.map"
    )

    return


def _unit_cell_hack(cell):
    """
    Helper function to check if alpha, beta and gamma are all 90, and if so,
    set alpha to 90.006. Otherwise, do nothing. a, b, c, beta, and gamma are never changed.

    This is necessary because coot assumes that maps in P1 with alpha, beta, gamma all = 90 are
    EM maps, and do not have periodic boundaries. Setting alpha=90.006 is sufficient to convince
    coot to render periodic boundaries.

    Parameters
    ----------
    cell : gemmi.UnitCell

    Returns
    -------
    cell_with_hack : gemmi.UnitCell
    """
    if all(
        [
            angle == 90
            for angle in (
                cell.alpha,
                cell.beta,
                cell.gamma,
            )
        ]
    ):
        return gemmi.UnitCell(
            cell.a,
            cell.b,
            cell.c,
            90.006,
            cell.beta,
            cell.gamma,
        )
    else:
        return cell


def align_grids_from_model_transform(grid1, grid2, structure1, structure2, selection):
    """
    This function is basically just a wrapper around `gemmi.interpolate_grid_of_aligned_model2`, which is an amazing thing that exists!!.

    Parameters
    ----------
    grid1 : gemmi.FloatGrid
        _description_
    grid2 : gemmi.FloatGrid
        _description_
    structure1 : gemmi.Structure
        _description_
    structure2 : gemmi.Structure
        _description_
    selection : str
        a string sufficient to describe the rbr selection to gemmi
        If the rbr selection contained multiple chains / spans, can just be the first chain / span

    Returns
    -------
    grid2_out : gemmi.FloatGrid
        Aligned, interpolated, and trimmed grid
    """

    if selection is None:
        span1 = structure1[0][0].get_polymer()
        span2 = structure2[0][0].get_polymer()

        dest_model = structure1[0]

    else:
        sel = gemmi.Selection(selection)

        span1 = sel.copy_structure_selection(structure1)[0][0].get_polymer()
        span2 = sel.copy_structure_selection(structure2)[0][0].get_polymer()

        dest_model = sel.copy_structure_selection(structure1)[0]

    sup = gemmi.calculate_superposition(
        span1, span2, span1.check_polymer_type(), gemmi.SupSelect.CaP
    )
    transform = sup.transform.inverse()

    # clone a grid to hold the output
    grid2_out = (
        grid1.clone()
    )  # this makes sure that grid2_out has a voxel frame matching grid1

    gemmi.interpolate_grid_of_aligned_model2(
        dest=grid2_out,
        src=grid2,
        tr=transform,
        dest_model=structure1[0],  # dest_model,
        radius=3,
        order=2,
    )

    return grid2_out


def _extract_pdb_selection(pdb, selection):
    # return a gemmi.Model containing only the selected things
    # selection can be either a string or a list of strings

    sel = gemmi.Selection(selection)

    pdb_selection = sel.copy_model_selection(pdb)

    return pdb_selection


def _ncs_align_and_subtract(
    fg,
    pdb,
    ncs_chains,
    name,
    output_dir,
):
    # here, ncs_chains would be ["B", "C"] for example
    # first element is 'fixed', second is 'moving'

    fg.unit_cell = _unit_cell_hack(fg.unit_cell)

    sup = gemmi.calculate_superposition(
        pdb[0][ncs_chains[0]].get_polymer(),
        pdb[0][ncs_chains[1]].get_polymer(),
        gemmi.PolymerType.PeptideL,
        gemmi.SupSelect.CaP,
    )

    fg2 = fg.clone()
    fg2.fill(0)

    gemmi.interpolate_grid_of_aligned_model2(
        dest=fg2,
        src=fg,
        tr=sup.transform.inverse(),
        dest_model=pdb[0],
        radius=8,
        order=2,
    )

    model = gemmi.Model("dummy model")
    model.add_chain(pdb[0][ncs_chains[0]])

    fg_mask_only = fg.clone()
    fg_mask_only.fill(0)

    masker = gemmi.SolventMasker(gemmi.AtomicRadiiSet.Cctbx)
    masker.rprobe = 4  # generous mask to be sure to leave no islands
    masker.put_mask_on_float_grid(fg_mask_only, model)

    masked_arr1 = np.logical_not(fg_mask_only.array) * fg.array
    masked_arr2 = np.logical_not(fg_mask_only.array) * fg2.array

    mask = np.logical_not(fg_mask_only.array.astype(bool))

    masked_arr1[mask] = _quicknorm(masked_arr1[mask])
    masked_arr2[mask] = _quicknorm(masked_arr2[mask])

    fg.set_subarray(arr=masked_arr1, start=(0, 0, 0))
    fg2.set_subarray(arr=masked_arr2, start=(0, 0, 0))

    print(f"{time.strftime('%H:%M:%S')}: Writing files...")

    write_maps = partial(
        rs.io.write_ccp4_map, cell=fg.unit_cell, spacegroup=fg.spacegroup
    )

    write_maps(fg.array, f"{output_dir}/{name}_{ncs_chains[0]}.map")
    write_maps(fg2.array, f"{output_dir}/{name}_{ncs_chains[1]}.map")

    write_maps(
        fg2.array - fg.array,
        f"{output_dir}/{name}_{ncs_chains[1]}_minus_{ncs_chains[0]}.map",
    )

    print(f"{time.strftime('%H:%M:%S')}: Done!")

    return


def _quicknorm(array):
    return (array - array.mean()) / array.std()
