import shutil
import subprocess
import time
from pathlib import Path

import reciprocalspaceship as rs


def _auto_eff_refinement_template(phenix_style: str):
    if phenix_style == "1.20":
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
        map_coefficients {
          map_type = "mFo-DFc"
          mtz_label_amplitudes = "FOFCWT"
          mtz_label_phases = "PHFOFCWT"
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
        bulk_solvent_and_scale=bss
        nqh_flips=False
      }
    }
        """
    elif phenix_style == "1.21":
        eff_contents = """
data_manager {
  model {
    file = pdb_input
  }
  miller_array {
    file = mtz_input
    labels {
      name = columns
    }
  }
  fmodel {
    xray_data {
        r_free_flags {
            generate = True
        }
        force_anomalous_flag_to_be_equal_to = False
    }
  }
}
refinement {
  crystal_symmetry {
    unit_cell = cell_parameters
    space_group = sg
  }
  output {
    write_def_file = False
    write_eff_file = False
    write_geo_file = False
  }
  electron_density_maps {
    apply_default_maps = False
    map_coefficients {
      map_type = "2mFo-DFc"
      mtz_label_amplitudes = "2FOFCWT"
      mtz_label_phases = "PH2FOFCWT"
    }
    map_coefficients {
      map_type = "mFo-DFc"
      mtz_label_amplitudes = "FOFCWT"
      mtz_label_phases = "PHFOFCWT"
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
    bulk_solvent_and_scale=bss
    nqh_flips=False
  }
  rigid_body {
    bulk_solvent_and_scale=bss
    high_resolution=None
  }
}
output {
  prefix = '''nickname'''
  serial = 1
  serial_format = "%d"
}
    """
    else:
        raise NotImplementedError("Unsupported phenix version")

    return eff_contents


def rigid_body_refinement_wrapper(
    mtzon,
    pdboff,
    input_dir,
    output_dir,
    phenix_style,
    off_labels=None,
    ligands=None,
    eff=None,
    verbose=False,
    rbr_selections=None,
    mr_on=False,
    no_bss=False,
):
    if eff is None:
        eff_contents = _auto_eff_refinement_template(phenix_style=phenix_style)
    else:
        with open(input_dir + eff) as file:
            eff_contents = file.read()

    if (off_labels is None) or mr_on:
        nickname = f"{mtzon.name.removesuffix('.mtz')}_rbr_to_{pdboff.name.removesuffix('.pdb')}"
    else:
        nickname = f"{mtzon.name.removesuffix('.mtz')}_rbr_to_self"
    ####
    # update this logic in the future if matchmaps.mr changes
    # mtz_location = input_dir if (mr_on or mr_off) else output_dir
    ####
    similar_files = list(output_dir.glob(f"{nickname}_[0-9]_1.*"))
    if len(similar_files) == 0:
        nickname += "_0"
    else:
        nums = []
        for s in similar_files:
            try:
                nums.append(int(str(s).split("_")[-2]))
            except ValueError:
                pass
        nickname += f"_{max(nums) + 1}"
    # read in mtz to access cell parameters and spacegroup
    cell_string, sg = _parse_mtz(mtzfile=str(mtzon))

    # name for modified refinement file
    eff = output_dir / f"params_{nickname}.eff"
    params = {
        "sg": sg,
        "cell_parameters": cell_string,
        "bss": str(not no_bss),
        "pdb_input": str(pdboff),
        "mtz_input": str(mtzon),
        "nickname": str(output_dir / nickname),
    }
    if off_labels is None:
        params["columns"] = "FPH1,SIGFPH1"  # names from scaleit output
    else:
        params["columns"] = off_labels  # user-provided column nanes
    # if selection is not None:
    #     params["all"] = selection  # overwrite atom selection
    for key, value in params.items():
        eff_contents = eff_contents.replace(key, value)
    # either add ligands to .eff file or delete "ligands" placeholder
    if ligands is not None:
        ligand_string = "\n".join([f"file_name = '{l}'" for l in ligands])
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

    _custom_subprocess(
        command="phenix.refine",
        params=str(eff),
        verbose=verbose,
    )

    return output_dir / nickname


def _auto_eff_phaser_template(phenix_style):
    if (phenix_style == "1.20") or (phenix_style == "1.21"):
        eff_contents = """
phaser {
  mode = ANO CCA EP_AUTO *MR_AUTO MR_FRF MR_FTF MR_PAK MR_RNP NMAXYZ SCEDS
  hklin = mtz_input
  labin = labels
  model = pdb_input
  model_identity = 100
  component_copies = 1
  search_copies = 1
  chain_type = *protein dna rna
  crystal_symmetry {
    unit_cell = cell_parameters
    space_group = sg
  }
  keywords {
    general {
      root = '''nickname'''
      title = '''matchmaps_MR'''
      mute = None
      xyzout = True
      xyzout_ensemble = True
      hklout = True
      jobs = 6
    }
  }
}
        """
    elif phenix_style == "1.21":
        eff_contents = """"""

    else:
        raise NotImplementedError(f"Phenix version {phenix_style} not supported")

    return eff_contents


def phaser_wrapper(
    mtzfile,
    pdb,
    output_dir,
    off_labels,
    phenix_style,
    eff=None,
    verbose=False,
):
    """
    Handle simple phaser run from the command line

    Args:
        phenix_style:
    """

    # this should never be needed; environment is already checked
    if shutil.which("phenix.phaser") is None:
        raise OSError(
            "Cannot find executable, phenix.phaser. Please set up your phenix environment."
        )

    if eff is None:
        eff_contents = _auto_eff_phaser_template(phenix_style=phenix_style)

    else:
        raise NotImplementedError("Custom phaser specifications are not yet supported")

    nickname = f"{mtzfile.name.removesuffix('.mtz')}_phased_with_{pdb.name.removesuffix('.pdb')}"

    similar_files = list(output_dir.glob(f"{nickname}_*"))
    if len(similar_files) == 0:
        nickname += "_0"
    else:
        nums = []
        for s in similar_files:
            try:
                nums.append(int(str(s).split("_")[-1].split(".")[0]))
            except ValueError:
                pass
        nickname += f"_{max(nums) + 1}"

    cell_string, sg = _parse_mtz(mtzfile)

    eff = output_dir / f"params_{nickname}.eff"

    params = {
        "sg": sg,
        "cell_parameters": cell_string,
        "pdb_input": str(pdb),
        "mtz_input": str(mtzfile),
        "nickname": str(output_dir / nickname),
        "labels": off_labels,  # should be prepackaged as a string
    }

    for key, value in params.items():
        eff_contents = eff_contents.replace(key, value)

    with open(eff, "w") as file:
        file.write(eff_contents)

    _custom_subprocess(command="phenix.phaser", params=str(eff), verbose=verbose)

    return output_dir / nickname


def _parse_mtz(mtzfile):
    mtz = rs.read_mtz(str(mtzfile))
    cell_string = f"{mtz.cell.a} {mtz.cell.b} {mtz.cell.c} {mtz.cell.alpha} {mtz.cell.beta} {mtz.cell.gamma}"
    sg = mtz.spacegroup.short_name()
    return cell_string, sg


def _renumber_waters(pdb, verbose):
    """
    Call phenix.sort_hetatms to place waters onto the nearest protein chain.
    This ensures that rbr selections handle waters properly

    Parameters
    ----------
    pdb : str
        name of pdb file
    dir : str
        directory in which pdb file lives
    verbose:
    """

    pdb_renumbered = Path(str(pdb).removesuffix(".pdb") + "_renumbered.pdb")

    _custom_subprocess(
        command="phenix.sort_hetatms",
        params=f"file_name={pdb} output_file={pdb_renumbered}",
        verbose=verbose,
    )

    print(f"{time.strftime('%H:%M:%S')}: Moved waters to nearest protein chains...")

    return pdb_renumbered


def _remove_waters(pdb, output_dir, verbose):
    pdb_dry = pdb.name.removesuffix(".pdb") + "_dry"

    _custom_subprocess(
        command="phenix.pdbtools",
        params=f"{pdb} remove='water' \
            output.prefix='{output_dir}/' \
            output.suffix='{pdb_dry}'",
        verbose=verbose,
    )

    return output_dir / (pdb_dry + ".pdb")


def _custom_subprocess(command, params, verbose, shell=True):
    subproc = subprocess.run(
        command + " " + params, shell=shell, capture_output=(not verbose)
    )

    if subproc.returncode != 0:
        error = f"matchmaps encountered an error while running {command}" + (
            "\n              Try again in --verbose mode for more a more detailed error message"
            if (not verbose)
            else ""
        )

        raise RuntimeError(error)

    return
