"""
Utilities and helperfunctions used by matchmaps.

Functions make_floatgrid_from_mtz, rigid_body_refinement_wrapper, and align_grids_from_model_transform
are exported to python for use in prototyping and testing
"""

import os
import shutil
import subprocess
import time
import re
from functools import partial
from pathlib import Path
import gemmi
import numpy as np
import reciprocalspaceship as rs
from pandas.api.extensions import ExtensionDtype


def _validate_environment(ccp4):
    """
    Check if the environment contains phenix (and if necessary, ccp4) and throw a helpful error if not

    If the function runs successfully, it returns a string of form "X.XX" denoting the major and minor version of the phenix detected, e.g. "1.20" or "1.21".
    """

    if shutil.which("phenix.refine") is None:
        raise OSError(
            "phenix is not active in your current environment. Please activate phenix and try again."
            "\n"
            "For more information, see https://rs-station.github.io/matchmaps/quickstart.html#additional-dependencies"
        )
    else:
        phenix_version = _detect_phenix_version()

    if ccp4:
        if shutil.which("scaleit") is None:
            raise OSError(
                "ccp4 is not active in your current environment. Please activate ccp4 and try again."
                "\n"
                "For more information, see https://rs-station.github.io/matchmaps/quickstart.html#additional-dependencies"
            )

    print(
        f"Detected phenix {phenix_version} in your environment.",
        "\n",
        "If this is not the version you are using, please specify the version directly via the --phenix-version flag",
    )

    return phenix_version


def _detect_phenix_version():
    version_printout = subprocess.run(
        "phenix.version | grep Version", shell=True, capture_output=True
    )
    version_string = str(version_printout.stdout)
    # if version_string.find('21') > 0:
    #     raise NotImplementedError("It seems that you are using phenix 1.21, which is not yet supported by matchmaps"
    #                               "\n"
    #                               "Please use phenix 1.20 or earlier.")
    phenix_version = ".".join(version_string.split(": ")[1].split(".")[:-1])
    return phenix_version


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


def make_floatgrid_from_mtz(
    mtz, spacing, F, Phi, SigF=None, spacegroup="P1", dmin=None, alpha=0
):
    """
    Make a gemmi.FloatGrid from an rs.DataSet.

    Parameters
    ----------
    mtz : rs.DataSet
        mtz data to be transformed into real space
    spacing : float
        Approximate voxel size desired (will be rounded as necessary to create integer grid dimensions)
    F : str
        Column in mtz containing structure factor amplitudes to use for calculation
    Phi : str
        Column in mtz containing phases to be used for calculation
    SigF : str, optional
        Column in mtz containing uncertainties of structure factor amplitudes.
        If provided, output map is weighted using the alpha value provided
    spacegroup : str, optional
        Spacegroup for the output FloatGrid. Defaults to P1.
    dmin: float, optional
        Highest resolution reflections to include in Fourier transform. Defaults to None, no cutoff.
    alpha: float, optional
        Alpha-factor to use in error-weighting of the output map. If SigF is not provided, this parameter is ignored.

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

    # apply weighting
    # note: if alpha==0, then these numbers all just become 1, e.g. no weighting
    if SigF is not None:
        weights = 1 / (1 + (alpha * new_mtz[SigF] ** 2 / np.mean(new_mtz[SigF] ** 2)))
    else:
        weights = 1

    if "weighted_Fobs" in new_mtz.columns:
        raise NotImplementedError(
            "Error: mtz already contains a column named weighted_Fobs; email Dennis bugging him to support this"
        )
    new_mtz["weighted_Fobs"] = new_mtz[F] * weights

    # perform FFT using the desired amplitudes and phases
    new_mtz["Fcomplex"] = new_mtz.to_structurefactor("weighted_Fobs", Phi)
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


def _handle_special_positions(pdboff, output_dir):
    """
    Check if any waters happen to sit on special positions, and if so, remove them.
    If any non-water atoms sit on special positions, throw a (hopefully helpful) error.

    Regardless of whether any special positions were found, copy this file over to output_dir

    Additionally, regardless of whether the input is a .pdb or .cif file, write the output to a .pdb file for downstream use.

    Parameters
    ----------
    pdboff : pathlib.Path
        name of input pdb (or mmcif)
    output_dir : str
    """
    pdb = gemmi.read_structure(str(pdboff))

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

    # if pdboff.suffix in ('.cif', '.CIF'):
    #     pdboff_nospecialpositions = output_dir / (pdboff.name.lower().removesuffix(".cif") + "_nospecialpositions.pdb")
    #     pdb.make_mmcif_document().write_file(str(pdboff_nospecialpositions))
    # else:
    pdboff_nospecialpositions = output_dir / (
        pdboff.name.removesuffix(".pdb") + "_nospecialpositions.pdb"
    )
    pdb.write_pdb(str(pdboff_nospecialpositions))

    return pdboff_nospecialpositions


def _restore_ligand_occupancy(
    pdb_to_be_restored,
    original_pdb,
    output_dir,
):
    # grab occupancies of all HETATMs in original_pdb
    with open(original_pdb, "r") as o:
        original = o.readlines()
    original_hetatm = []
    for l in original:
        if ("HETATM" in l) and (not "REMARK" in l):
            original_hetatm.append(l)
    original_occs = [h[56:60] for h in original_hetatm]

    with open(pdb_to_be_restored, "r") as p:
        pdb = p.readlines()
    n = 0
    for i in range(len(pdb)):
        if ("HETATM" in pdb[i]) and (not "REMARK" in pdb[i]):
            pdb[i] = pdb[i][:56] + original_occs[n] + pdb[i][60:]
            n += 1

    edited_pdb = output_dir / (
        original_pdb.name.removesuffix(".pdb") + "_restorehetatms.pdb"
    )

    with open(edited_pdb, "w") as output:
        output.write("".join(pdb))

    return edited_pdb


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
    radius,
):
    """
    Take in two floatgrids and two pdbs. Apply the alignment transformation and subtract.

    If there were multiple selections used for rigid-body refinement, this method should be called separately for each selection.

    Parameters
    ----------
    output_dir : pathlib.Path
        Where to write output files
    fg_off : gemmi.FloatGrid
        FloatGrid containing electron density for "off" data
    fg_on : gemmi.FloatGrid
        FloatGrid containing electron density for "on" data
    pdboff : gemmi.Structure
        Reference structure, rigid-body-refined to the "off" data
    pdbon : gemmi.Structure
        Reference structure, rigid-body-refined to the "on" data
    on_name : str
        Original base name of the input "on" data
    off_name : str
        Original base name of the input "off" data
    on_as_stationary : bool
        If True, align "off" data onto "on" data. If False, align "on" data onto "off" data.
    selection : str or list of str
        If not None, should be a valid gemmi selection string or a list of valid gemmi selection strings
    radius : float
        Maximum distance away from protein model to include voxels. Only applies to the "unmasked" output.
    """

    if selection:
        print(
            f"{time.strftime('%H:%M:%S')}: Using models to rigid-body align maps for rigid-body selection {selection}..."
        )
    else:
        print(f"{time.strftime('%H:%M:%S')}: Using models to rigid-body align maps...")

    # write out "raw" maps, prior to any alignment
    rs.io.write_ccp4_map(
        fg_on.array,
        str(output_dir / (on_name + "_before.map")),
        cell=fg_on.unit_cell,
        spacegroup=fg_on.spacegroup,
    )
    rs.io.write_ccp4_map(
        fg_off.array,
        str(output_dir / (off_name + "_before.map")),
        cell=fg_off.unit_cell,
        spacegroup=fg_off.spacegroup,
    )

    # pick which grid/model is "fixed"
    if on_as_stationary:
        fg_fixed = fg_on.clone()
        pdb_fixed = pdbon.clone()
    else:
        fg_fixed = fg_off.clone()
        pdb_fixed = pdboff.clone()

    # Use models to align grids
    # align_grids_from_model_transform is a light wrapper around gemmi.interpolate_grid_of_aligned_model2 (previously called gemmi.interpolate_grid_of_aligned_model2)
    fg_off = align_grids_from_model_transform(
        fg_fixed, fg_off, pdb_fixed, pdboff, selection, radius=radius
    )
    fg_on = align_grids_from_model_transform(
        fg_fixed, fg_on, pdb_fixed, pdbon, selection, radius=radius
    )

    # These grids contain mostly zeros (as governed by the radius supplied above)
    # Store the locations of these zeros as a boolean array so we can ignore them
    loose_mask = np.invert(fg_on.array == 0)

    # Normalize the non-zero values in each grid
    fg_on.array[loose_mask] = _quicknorm(fg_on.array[loose_mask])
    fg_off.array[loose_mask] = _quicknorm(fg_off.array[loose_mask])

    # Subtract
    difference_array = fg_on.array - fg_off.array

    # Normalize again, this time on the differences
    # This ensures that the map looks the same in units of e/A^3 and sigma
    difference_array[loose_mask] = _quicknorm(difference_array[loose_mask])

    # Use gemmi to compute a second, more stringent mask with a radius of 2
    fg_mask_only = fg_fixed.clone()
    masker = gemmi.SolventMasker(gemmi.AtomicRadiiSet.Cctbx)
    masker.rprobe = (
        2  # the user's --unmasked-radius parameter does NOT change this value
    )

    if selection is None:
        pdb_for_mask = pdb_fixed[0]
    else:
        pdb_for_mask = _extract_pdb_selection(pdb_fixed[0], selection)

    masker.put_mask_on_float_grid(fg_mask_only, pdb_for_mask)

    tight_mask = np.logical_not(fg_mask_only.array).astype(bool)

    # apply tighter mask to a the difference array
    masked_difference_array = tight_mask * difference_array

    # NOT doing this final normalization for now:
    # masked_difference_array[tight_mask] = _quicknorm(masked_difference_array[tight_mask])

    # coot refuses to render periodic boundaries for P1 maps with alpha=beta=gamma=90, sooooo
    fg_fixed.unit_cell = _unit_cell_hack(fg_fixed.unit_cell)

    # use partial function to guarantee I'm always using the same and correct cell
    write_maps = partial(
        rs.io.write_ccp4_map, cell=fg_fixed.unit_cell, spacegroup=fg_fixed.spacegroup
    )

    # write files!
    print(f"{time.strftime('%H:%M:%S')}: Writing files...")

    write_maps(fg_on.array, str(output_dir / (on_name + ".map")))

    write_maps(fg_off.array, str(output_dir / (off_name + ".map")))

    write_maps(
        masked_difference_array,
        str(output_dir / (on_name + "_minus_" + off_name + ".map")),
    )
    write_maps(
        difference_array,
        str(output_dir / (on_name + "_minus_" + off_name + "_unmasked.map")),
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


def align_grids_from_model_transform(
    grid1, grid2, structure1, structure2, selection, radius
):
    """
    This function is basically just a wrapper around `interpolate_grid_around_model` (previously called `gemmi.interpolate_grid_of_aligned_model2`),
    which is an amazing thing that exists!!.

    Parameters
    ----------
    grid1 : gemmi.FloatGrid
        FloatGrid containing reference (static) electron density
    grid2 : gemmi.FloatGrid
        FloatGrid containing electron density to be moved
    structure1 : gemmi.Structure
        Structure containing reference (static) protein model.
    structure2 : gemmi.Structure
        Structure which, when aligned to structure1, reveals the necessary transformation to align grid2 onto grid1.
    selection : str
        a string sufficient to describe the rbr selection to gemmi
        If the rbr selection contained multiple chains / spans, can just be the first chain / span
    radius : float
        Voxels up to this far away (in Angstroms) from the protein model will be kept in the "unmasked" output map

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

    grid2_out.fill(0)

    gemmi.interpolate_grid_around_model(
        dest=grid2_out,
        src=grid2,
        tr=transform,
        dest_model=structure1[0],  # dest_model,
        radius=radius,
        order=1,
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

    gemmi.interpolate_grid_around_model(
        dest=fg2,
        src=fg,
        tr=sup.transform.inverse(),
        dest_model=pdb[0],
        radius=8,
        order=1,
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

    write_maps(fg.array, str(output_dir / (name + "_" + ncs_chains[0] + ".map")))
    write_maps(
        fg2.array,
        str(
            output_dir
            / (name + "_" + ncs_chains[1] + "_onto_" + ncs_chains[0] + ".map")
        ),
    )

    write_maps(
        fg2.array - fg.array,
        str(
            output_dir
            / (name + "_" + ncs_chains[1] + "_minus_" + ncs_chains[0] + ".map")
        ),
    )

    print(f"{time.strftime('%H:%M:%S')}: Done!")

    return


def _quicknorm(array):
    return (array - array.mean()) / array.std()


def _validate_inputs(
    input_dir: str,
    output_dir: str,
    ligands,
    *files,
):
    # use pathlib to validate input files and directories
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    input_dir = Path(input_dir)

    if not input_dir.exists():
        raise ValueError(f"Input directory '{input_dir}' does not exist")

    if ligands is not None:
        ligands = [Path(ligand) for ligand in ligands]

        ligands = [
            ligand if ligand.is_absolute() else input_dir / ligand for ligand in ligands
        ]

        for ligand in ligands:
            if not ligand.exists():
                raise ValueError(f"Input ligand '{ligand}' does not exist")

    files = [Path(f) for f in files]

    files = [f if f.is_absolute() else input_dir / f for f in files]

    for f in files:
        if not f.exists():
            raise ValueError(f"Input file '{f}' does not exist")

    return input_dir, output_dir, ligands, *files


def _validate_column_dtypes(
    mtz: rs.DataSet,
    columns: tuple[str, ...],
    dtypes: tuple[ExtensionDtype, ...],  # ????
):
    """
    Validate dtypes of input columns
    """

    for column, dtype in zip(columns, dtypes):
        try:
            if not dtype.is_dtype(mtz[column].dtype):
                raise ValueError(
                    f"MTZ column '{column}' has dtype '{mtz[column].dtype}'; expected dtype '{dtype.name}'"
                )
        except KeyError:
            raise KeyError(f"Column '{column}' does not exist in the input MTZ")


def _cif_or_mtz_to_mtz(input_file, output_dir):
    """
    Return a File() object and also the original file name as a string

    Parameters
    ----------
    input_file : pathlib.Path
        _description_
    output_dir : pathlib.Path
        _description_

    Returns
    -------
    (output_file, filename)

    """

    if input_file.suffix.lower() == ".mtz":
        output_file = output_dir / (input_file.name)

        try:
            shutil.copy(input_file, output_file)
        except shutil.SameFileError:
            print(
                f"Note: Input file '{input_file.name}' is located in the output directory '{output_dir}'\n"
                "      This is ok, but I recommend directing outputs elsewhere using the --output-dir option"
            )
            pass

    elif input_file.suffix.lower() == ".cif":
        output_file = output_dir / (
            input_file.name.lower().removesuffix(".cif") + ".mtz"
        )

        reflections = rs.read_cif(str(input_file))

        reflections.write_mtz(str(output_file))

    else:
        raise ValueError(
            f"Invalid file type {input_file.suffix} for starting model, must be '.mtz' or '.cif'"
        )

    return (output_file, input_file.name.removesuffix(input_file.suffix))


def _cif_or_pdb_to_pdb(input_file, output_dir):
    """
    _summary_

    Parameters
    ----------
    input_file : pathlib.Path
    output_dir : pathlib.Path

    Returns
    -------
    pathlib.Path
        path to output file

    Raises
    ------
    ValueError
        _description_
    """

    if input_file.suffix.lower() == ".pdb":
        output_file = output_dir / (input_file.name)

        try:
            shutil.copy(input_file, output_file)
        except shutil.SameFileError:
            print(
                f"Note: Input file '{input_file.name}' is located in the output directory '{output_dir}'\n"
                "      This is ok, but I recommend directing outputs elsewhere using the --output-dir option"
            )
            pass

    elif input_file.suffix.lower() == ".cif":
        output_file = output_dir / (
            input_file.name.lower().removesuffix(".cif") + ".pdb"
        )

        structure = gemmi.read_structure(str(input_file))
        structure.write_pdb(str(output_file))

    else:
        raise ValueError(
            f"Invalid file type {input_file.suffix} for starting model, must be '.pdb' or '.cif'"
        )

    return output_file


def _clean_up_files(output_dir, old_files, keep_temp_files):
    candidate_files = []
    for suffix in ("eff", "pdb", "mtz", "log", "cif"):
        candidate_files.extend(list(output_dir.glob("*" + suffix)))

    files_to_delete = set(candidate_files) - set(old_files)

    if keep_temp_files is not None:
        new_dir = output_dir / keep_temp_files
        new_dir.mkdir(parents=True, exist_ok=True)
        for f in files_to_delete:
            f.rename(new_dir / f.name)
    else:
        for f in files_to_delete:
            os.remove(f)

    return


def _write_script(utility, arguments, script_name, phenix_version):
    from matchmaps import __version__ as version

    if phenix_version is None:
        phenix_version = _detect_phenix_version()

    contents = f"""#!/bin/bash

# This file was produced on {time.strftime('%c')}
# Using matchmaps version {version} and phenix version {phenix_version}
#
# The command below was originally run in the following directory:
# {os.getcwd()}

{utility} {' '.join(arguments)}

"""

    with open(script_name + ".sh", "w") as file:
        file.write(contents)

    return
