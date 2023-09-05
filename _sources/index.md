# `matchmaps` Documentation

Welcome to the docs for `matchmaps`, a python package for aligning and subtracting non-isomorphous crystallographic datasets. Do you have a pair of datasets for which an isomorphous difference maps doesn't quite work due to changes in your unit cell dimensions? If so, this is the package for you.

The [quickstart guide](quickstart.md) contains installation instructions and sample usage of the basic `matchmaps` utility, along with other useful tips.

If you'd like to learn more about `matchmaps` and see the package in action, please check out our pre-print!
> [MatchMaps: Non-isomorphous difference maps for X-ray crystallography](https://www.biorxiv.org/content/10.1101/2023.09.01.555333v1.full.pdf+html)

This software is part of the [Reciprocal Space Station](https://rs-station.github.io/) family of open-source crystallography software and was conceived in [Doeke Hekstra's Lab](https://hekstralab.fas.harvard.edu/) at Harvard by [Dennis Brookner](https://dennisbrookner.github.io/)

## The `matchmaps` functions in brief

`matchmaps` consists of three command-line utilities. A walkthrough of the core utility can be found [here](quickstart.md), and all options for all three utilites are documented [here](cli.md). All three utilities produce a real-space difference map and require a single starting model for phasing, but differ in the types of reflection data they require. Briefly:

 - **`matchmaps`** takes in two mtzs containing similar data and which are nearly isomorphous and computes an unbiased real-space difference map between them
 - **`matchmaps.mr`** takes in two mtzs containing similar data but which are in different spacegroups (or the same spacegroup but different crystal packing) and computes an unbiased real-space difference map between them.
 - **`matchmaps.ncs`** takes in a single mtz and computes an internal difference map across a defined non-crystallographic symmetry present in the data.


```{eval-rst}
.. toctree::
   :maxdepth: 1
   :hidden:

   Quickstart guide <quickstart>
   Command-line options <cli>
   About the algorithm <about>

```
