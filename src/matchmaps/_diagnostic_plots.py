"""Make diagnostic plots similar to those in Figure 1 of the MatchMaps paper"""

from pathlib import Path

import reciprocalspaceship as rs
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams

# to do list:
# [  ] write function for making plot
# [  ] write main function which directs command-line arguments into the plot
# [  ] write argument parser (in the _parsers and _args files), reusing code where possible
# [  ] write a useful docstring
# [  ] make sure matchmaps functions aren't broken!


def compute_diagnostic_plots(
    mtzoff: Path,
    mtzon: Path,
    Foff: str,
    Fon: str,
    Phioff: str,
    Phion: str,
    bins: int = 15,
    dmin: float = None,
    input_dir=Path("."),
    output_dir=Path("."),
    filename=None,
    title = 'Produced by matchmaps.diagnose',
):
    rcParams.update({"figure.autolayout": True})
    rcParams.update({"font.size": 16})

    offmtz = rs.read_mtz(str(input_dir / mtzoff))
    onmtz = rs.read_mtz(str(input_dir / mtzon))

    offmtz.canonicalize_phases(inplace=True)
    onmtz.canonicalize_phases(inplace=True)

    mtz = offmtz.merge(onmtz, on=["H", "K", "L"], suffixes=("_off", "_on"))

    mtz = mtz.label_centrics(inplace=True).query("not CENTRIC").dropna()

    if dmin:
        mtz.compute_dHKL(inplace=True)
        mtz = mtz.loc[mtz.dHKL >= dmin]

    for newname, oldname in zip(
        ("cosphi_off", "cosphi_on"), (f"{Phioff}_off", f"{Phion}_on")
    ):

        mtz[newname] = np.cos(mtz[oldname] * np.pi / 180)

    plt.figure(figsize=(8, 7))
    plt.grid(linestyle="--")
    plt.title(title)
    plt.hlines(0, 0, bins-1, color='gray', linestyle='-')
    plt.xlabel(r"Resolution bins ($\AA$)")
    plt.ylabel(r"$CC_{1/2}$")

    edges = cc_datasets(
        mtz,
        "cosphi_off",
        "cosphi_on",
        marker="o",
        color="#1f77b4",
        linestyle="-",
        label="cos(Phase) correlation",
        bins=bins,
    )
    cc_datasets(
        mtz,
        f"{Foff}_off",
        f"{Fon}_on",
        bins=edges,
        color="#be4f30",
        marker="o",
        linestyle="--",
        label="Amplitude correlation",
        return_edges=False
    )

    plt.legend(loc="upper right")
    if filename:
        plt.savefig(output_dir / filename, dpi=600)
    plt.show()
    return mtz


def cc_datasets(
    mtz,
    column1,
    column2,
    bins=15,
    return_edges=True,
    method="spearman",
    label=None,
    **kwargs,
):
    mtz, labels, edges = mtz.assign_resolution_bins(
        bins, return_labels=True, return_edges=True
    )
    nbins = len(edges) - 1
    grouper = mtz.groupby(["bin"])[[column1, column2]]

    result = (
        grouper.corr(method=method)
        .unstack()[(column1, column2)]
        .to_frame()
        .reset_index()
    )
    plt.plot(result.bin, result[column1], label=label, **kwargs)
    plt.xticks(range(nbins), labels, rotation=45, ha="right", rotation_mode="anchor")

    if return_edges:
        return edges
    else:
        return


def main():

    print("entered diagnostic plots")

    compute_diagnostic_plots()

    return


if __name__ == "__main__":
    main()
