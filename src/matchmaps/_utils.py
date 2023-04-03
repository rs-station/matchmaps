"""
Utilities and helperfunctions used by matchmaps.

Functions make_floatgrid_from_mtz, rigid_body_refinement_wrapper, and align_grids_from_model_transform
are exported to python for use in prototyping and testing
"""

import glob
import shutil
import subprocess

import gemmi
import numpy as np
import reciprocalspaceship as rs


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
        Column in mtz containing structure factor amplitudes to use for calculation, by default "2FOFCWT"
    Phi : str, optional
        Column in mtz containing phases to be used for calculation, by default "PH2FOFCWT"
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
    selection=None,
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
      rigid_body = all
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

    if selection is not None:
        params["all"] = selection  # overwrite atom selection

    for key, value in params.items():
        eff_contents = eff_contents.replace(key, value)

    # either add ligands to .eff file or delete "ligands" placeholder
    if ligands is not None:
        ligand_string = "\n".join([f"file_name = '{input_dir}/{l}'" for l in ligands])
        eff_contents = eff_contents.replace("ligands", ligand_string)
    else:
        eff_contents = eff_contents.replace("ligands", "")

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

    # if off_labels is not None:
    #     print(f"                                ^^^ Use this file as --mtzoff for mapreg.register ")
    # else:
    #     print(f"                                ^^^ Use this file as --mtzon for mapreg.register ")

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
                                "Input model contains water at special position. Removing water so as to not break rigid-body refinement."
                            )
                            print(
                                "If it is important that you keep this water and just omit it from your refinement selection, use the --selection flag"
                            )
                            residue.remove_atom(
                                name=atom.name, altloc=atom.altloc, el=atom.element
                            )

                        else:
                            raise ValueError(
                                """
Input model contains a non-water atom on a special position.
Please use the --selection flag to supply an input suitable for rigid-body refinement by phenix.
"""
                            )

    pdboff_nospecialpositions = pdboff.removesuffix(".pdb") + "_nospecialpositions.pdb"

    pdb.write_pdb(output_dir + pdboff_nospecialpositions)

    return pdboff_nospecialpositions


def align_grids_from_model_transform(grid1, grid2, structure1, structure2):
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

    Returns
    -------
    grid2_out : gemmi.FloatGrid
        Aligned, interpolated, and trimmed grid
    """
    span1 = structure1[0]["A"].get_polymer()
    span2 = structure2[0]["A"].get_polymer()
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
        dest_model=structure1[0],
        radius=3,
        order=2,
    )

    return grid2_out
