"""
Dictionaries for arguments parsed by one or more matchmaps utilities
"""

matchmaps_description = (
            "Compute a real-space difference map. "
            "You will need two MTZ files, which will be referred to throughout as 'on' and 'off', "
            "though they could also be light/dark, bound/apo, mutant/WT, hot/cold, etc. "
            "Each mtz will need to contain structure factor amplitudes and uncertainties; you will not need any phases. "
            "You will, however, need an input model (assumed to correspond with the 'off' state) which will be used to determine phases. "
            "The input file may be in .pdb or .cif format. "
            "Please note that both ccp4 and phenix must be installed and active in your environment for this function to run. "
            ""
            "More information can be found online at https://rs-station.github.io/matchmaps/index.html"
        )

matchmaps_mr_description = ()

matchmaps_ncs_description = ()

common_args = (
    # args used by all three utilities
    (("--ligands", "-l"),  {
        "help":"Any .cif restraint files needed for refinement",
        "required": False,
        "default": None,
        "nargs": "*",
    }),

    (("--input-dir", "-i"), {
        "help": "Path to input files. Optional, defaults to './' (current directory)",
        "required": False,
        "default": "./",
    }),

    (("--output-dir", "-o"), {
        "help": "Path to which output files should be written. Optional, defaults to './' (current directory)",
        "required": False,
        "default": "./",
    }),

)