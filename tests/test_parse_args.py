from matchmaps._parsers import (
    matchmaps_parser,
    matchmaps_mr_parser,
    matchmaps_ncs_parser,
    matchmaps_diagnose_parser,
)

import pytest


@pytest.mark.parametrize(
    "cli_extras",
    (
        ("--ligands", "data/raw_data/FOL.cif", "data/raw_data/NAP.cif"),
        ("--script", "script"),
        ("--dmin", "3"),
        ("--on-as-stationary",),
    ),
)
def test_mm_mr_full_paths(cli_extras):
    args = (
        "--mtzoff",
        "data/raw_data/1rx2_phases.mtz",
        "FP",
        "SIGFP",
        "--mtzon",
        "data/raw_data/1rx4_phases.mtz",
        "FP",
        "SIGFP",
        "--pdboff",
        "data/raw_data/1rx4.pdb",
    ) + cli_extras

    matchmaps_parser.parse_args(args)

    matchmaps_mr_parser.parse_args(args)

    return


@pytest.mark.parametrize(
    "cli_extras",
    (
        ("--ligands", "FOL.cif", "NAP.cif"),
        ("--script", "script"),
        ("--dmin", "3"),
        ("--on-as-stationary",),
    ),
)
def test_mm_mr_provide_dir(cli_extras):
    args = (
        "--input-dir",
        "data/raw_data",
        "--mtzoff",
        "1rx2_phases.mtz",
        "FP",
        "SIGFP",
        "--mtzon",
        "1rx4_phases.mtz",
        "FP",
        "SIGFP",
        "--pdboff",
        "1rx4.pdb",
    ) + cli_extras

    matchmaps_parser.parse_args(args)

    matchmaps_mr_parser.parse_args(args)

    return


@pytest.mark.parametrize(
    "cli_args",
    (
        (
            "--mtzon",
            "data/raw_data/1rx2_phases.mtz",
            "FP",
            "SIGFP",
        ),
        ("--ligands", "data/raw_data/FOL.cif", "data/raw_data/NAP.cif"),
        ("--script", "script"),
        ("--dmin", "3"),
    ),
)
def test_mm_missingargs(cli_args):
    with pytest.raises(SystemExit):
        matchmaps_parser.parse_args(cli_args)


@pytest.mark.parametrize(
    "cli_extras",
    (
        ("--ligands", "data/raw_data/FOL.cif", "data/raw_data/NAP.cif"),
        ("--script", "script"),
        ("--dmin", "3"),
    ),
)
def test_ncs_full_paths(cli_extras):
    args = (
        "--mtz",
        "data/raw_data/1rx2_phases.mtz",
        "FP",
        "SIGFP",
        "--pdb",
        "data/raw_data/1rx4.pdb",
        "--ncs-chains",
        "A",
        "B",
    ) + cli_extras

    matchmaps_ncs_parser.parse_args(args)

    return


@pytest.mark.parametrize("cli_extras", (("--filename", "name.png"), ("--dmin", "3")))
def test_diagnose(cli_extras):
    args = (
        "--mtzoff",
        "data/raw_data/1rx2_phases.mtz",
        "FP",
        "PHIC",
        "--mtzon",
        "data/raw_data/1rx4_phases.mtz",
        "FP",
        "PHIC",
    ) + cli_extras

    matchmaps_diagnose_parser.parse_args(args)
