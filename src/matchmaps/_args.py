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

matchmaps_mr_description = (
    "Compute a real-space difference map between inputs in different space groups / crystal packings. "
    "You will need two MTZ files, which will be referred to throughout as 'on' and 'off', "
    "though they could also be light/dark, bound/apo, mutant/WT, hot/cold, etc. "
    "Each mtz will need to contain structure factor amplitudes and uncertainties; you will not need any phases. "
    "You will, however, need an input model (assumed to correspond with the 'off' state) which will be used to determine phases. "
    "The input file may be in .pdb or .cif format. "
    "Please note that phenix must be installed and active in your environment for this function to run. "
    ""
    "If your mtzoff and mtzon are in the same spacegroup and crystal packing, see the basic matchmaps utility. "
    "More information can be found online at https://rs-station.github.io/matchmaps/index.html"
)

matchmaps_ncs_description = (
    "Compute an 'internal' real-space difference map between NCS-related molecules. "
    "You will need an MTZ file with structure factor amplitudes and optionally containing phases, and a PDB/CIF file. "
    ""
    "Please note that phenix must be installed and active in your environment for this function to run. "
    ""
    "More information can be found online at https://rs-station.github.io/matchmaps/index.html"
)

matchmaps_diagnose_description = (
    "Plot the isomorphism of two datasets by resolution bin; "
    "this plot is analogous to Figure 1a from the matchmaps paper. "
    ""
    "By default, the plot will open in a separate window, but you can save it with the `--filename` flag"
)

base_and_mr_args = (
    (
        ("--mtzoff", "-f"),
        {
            "nargs": 3,
            "metavar": ("mtzfileoff", "Foff", "SigFoff"),
            "required": True,
            "help": (
                "MTZ or sfCIF containing off/apo/ground/dark state data. "
                "Specified as [filename F SigF]"
            ),
        },
    ),
    (
        ("--mtzon", "-n"),
        {
            "nargs": 3,
            "metavar": ("mtzfileon", "Fon", "SigFon"),
            "required": True,
            "help": (
                "MTZ or sfCIF containing on/bound/excited/bright state data. "
                "Specified as [filename F SigF]"
            ),
        },
    ),
    (
        ("--pdboff", "-p"),
        {
            "required": True,
            "help": (
                "Reference PDB or mmCIF corresponding to the off/apo/ground/dark state. "
                "Used for rigid-body refinement of both input MTZs to generate phases."
            ),
        },
    ),
    (
        ("--on-as-stationary",),
        {
            "required": False,
            "action": "store_true",
            "default": False,
            "help": (
                "Include this flag to align 'off' data onto 'on' data. By default, 'off' data is stationary and 'on' data is moved."
            ),
        },
    ),
    (
        ("--alpha",),
        {
            "required": False,
            "type": float,
            "default": 0,
            "help": (
                "Alpha to use for error weighting of F-obs prior to Fourier Transform. "
                "Weights are computed as: 1 / ((1+(alpha*(SigF^2)) / <SigF>^2). "
                "Default value is alpha=0, e.g., no weighting is performed. "
            ),
        },
    ),
    (
        ("--unmasked-radius",),
        {
            "required": False,
            "type": float,
            "default": 5,
            "help": (
                "Maximum distance (in Anstroms) away from protein model to include voxels. Only applies to the 'unmasked' difference map output. "
                "Defaults to 5. "
                "Note that the regular difference map (e.g. the 'masked' version) is not affected by this parameter and maintains a solvent mask radius of 2 Angstroms."
            ),
        },
    ),
    (
        ("--rbr-selections", "-r"),
        {
            "required": False,
            "default": None,
            "nargs": "*",
            "help": (
                "Specification of multiple rigid-body groups for refinement. By default, everything is refined as one rigid-body group. "
            ),
        },
    ),
)

ncs_args = (
    (
        ("--mtz", "-m"),
        {
            "nargs": "*",
            "required": True,
            "help": (
                "MTZ or sfCIF file containing structure factor amplitudes. "
                "Specified as [filename F SigF] or [filename F]. "
                "SigF is not necessary if phases are also provided"
            ),
        },
    ),
    (
        ("--phases",),
        {
            "required": False,
            "default": None,
            "help": (
                "Optional. Column in MTZ/sfCIF file containing phases. "
                "If phases are not provided, phases will be computed via rigid-body refinement of "
                "the provided model and structure factor amplitudes."
            ),
        },
    ),
    (
        ("--pdb", "-p"),
        {
            "required": True,
            "help": (
                "Reference PDB or mmCIF. "
                "If phases are not provided, used for rigid-body refinement of input MTZ/sfCIF to generate phases."
            ),
        },
    ),
    (
        ("--ncs-chains", "-n"),
        {
            "required": True,
            "metavar": ("fixed_chain", "moving_chain"),
            "default": None,
            "nargs": 2,
            "help": (
                "NCS chains to overlay and subtract, specified as [fixed_chain, moving_chain]. "
                "E.g. to overlay chain C onto chain B, specify: `--ncs-chains B C`"
            ),
        },
    ),
    (
        ("--mapname",),
        {
            "required": False,
            "default": "matchmaps_ncs",
            "help": "Base filename for the output map files. ",
        },
    ),
)

common_matchmaps_args = (
    # args used by all three utilities
    (
        ("--ligands", "-l"),
        {
            "help": "Any .cif restraint files needed for refinement",
            "required": False,
            "default": None,
            "nargs": "*",
        },
    ),
    (
        ("--spacing", "-s"),
        {
            "help": (
                "Approximate voxel size in Angstroms for real-space maps. Defaults to 0.5 A. "
                "Value is approximate because there must be an integer number of voxels along each unit cell dimension"
            ),
            "required": False,
            "default": 0.5,
            "type": float,
        },
    ),
    (
        ("--no-bss",),
        {
            "help": (
                "Include this flag to skip bulk solvent scaling in phenix.refine. By default, BSS is included."
            ),
            "required": False,
            "action": "store_true",
            "default": False,
        },
    ),
    (
        ("--verbose", "-v"),
        {
            "help": "Include this flag to print out scaleit and phenix.refine outputs to the terminal. Useful for troubleshooting, but annoying; defaults to False.",
            "required": False,
            "action": "store_true",
            "default": False,
        },
    ),
    (
        ("--eff",),
        {
            "help": "Custom .eff template for running phenix.refine. ",
            "required": False,
            "default": None,
        },
    ),
    (
        ("--keep-temp-files", "-k"),
        {
            "required": False,
            "default": None,
            "help": (
                "Do not delete intermediate matchmaps files, but rather place them in the supplied directory. "
                "This directory is created as a subdirectory of the supplied `output-dir`."
            ),
        },
    ),
    (
        ("--script",),
        {
            "required": False,
            "default": "run_matchmaps",
            "help": (
                "Name for a file {script}.sh which can be run to repeat this command. "
                "By default, this file is called `run_matchmaps.sh`. "
                "Note that this file is written out in the current working directory, NOT the input or output directories"
            ),
        },
    ),
    (
        ("--phenix-version",),
        {
            "required": False,
            "help": (
                "Specify phenix version as a string, e.g. '1.20'. "
                "If omitted, matchmaps will attempt to automatically detect the version in use "
                "by analyzing the output of phenix.version"
            ),
        },
    ),
)

global_args = (
    (
        ("--input-dir", "-i"),
        {
            "help": "Path to input files. Optional, defaults to './' (current directory)",
            "required": False,
            "default": "./",
        },
    ),
    (
        ("--output-dir", "-o"),
        {
            "help": "Path to which output files should be written. Optional, defaults to './' (current directory)",
            "required": False,
            "default": "./",
        },
    ),
    (
        ("--dmin",),
        {
            "help": (
                "Highest-resolution (in Angstroms) reflections to include in Fourier transform for FloatGrid creation. "
                "By default, cutoff is the resolution limit of the lower-resolution input MTZ. "
            ),
            "required": False,
            "type": float,
            "default": None,
        },
    ),
)

diagnose_args = (
    # args used by only matchmaps.diagnose
    # these differ from the standard mtzoff/mtzon because the third argument is phase, not SigF!
    (
        ("--mtzoff", "-f"),
        {
            "nargs": 3,
            "metavar": ("mtzfileoff", "Foff", "Phioff"),
            "required": True,
            "help": (
                "MTZ or sfCIF containing off/apo/ground/dark state data. "
                "Specified as [filename F Phi] (no need for SigF here!)"
            ),
        },
    ),
    (
        ("--mtzon", "-n"),
        {
            "nargs": 3,
            "metavar": ("mtzfileon", "Fon", "Phioff"),
            "required": True,
            "help": (
                "MTZ or sfCIF containing on/bound/excited/bright state data. "
                "Specified as [filename F Phi] (no need for SigF here!)"
            ),
        },
    ),
    (
        ("--filename",),
        {
            "required": False,
            "help": (
                "Filename for the plot to save; if it has an extension, it should be compatible with plt.savefig. "
                "If omitted, plot will be opened interactively but no file will be saved. "
                "File will be saved in the `--output-dir` directory (current directory by default)"
            ),
        },
    ),
    (
        ("--title", "-t"),
        {
            "required": False,
            "default": "Produced by matchmaps.diagnose",
            "help": "Title for the plot",
        },
    ),
    (
        ("--bins", "-b"),
        {
            "required": False,
            "default": 15,
            "help": "Numbers of resolution bins to compute correlation",
            "type": int,
        },
    ),
)
