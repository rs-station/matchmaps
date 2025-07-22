# `matchmaps` Documentation

Welcome to the docs for `matchmaps`, a python package for aligning and subtracting non-isomorphous crystallographic datasets. Do you have a pair of datasets for which an isomorphous difference maps doesn't quite work due to changes in your unit cell dimensions? If so, this is the package for you.

The [quickstart guide](quickstart.md) contains installation instructions and sample usage of the basic `matchmaps` utility, along with other useful tips.

If you'd like to learn more about `matchmaps`, you can find our exploration of the package [here](about.md). For even more information, and see the package in action, please check out our paper in the Journal of Applied Crystallography!

> [MatchMaps: non-isomorphous difference maps for X-ray crystallography](https://journals.iucr.org/j/issues/2024/03/00/ei5112/index.html)

This software is part of the [Reciprocal Space Station](https://rs-station.github.io/) family of open-source crystallography software and was conceived in [Doeke Hekstra's Lab](https://hekstralab.fas.harvard.edu/) at Harvard by [Dennis Brookner](https://dennisbrookner.github.io/)

## The `matchmaps` functions in brief

`matchmaps` consists of three main command-line utilities. A walkthrough of the core utility can be found [here](quickstart.md), and all options for all three utilites are documented [here](cli.md). All three utilities produce a real-space difference map and require a single starting model for phasing, but differ in the types of reflection data they require. Briefly:

 - **`matchmaps`** takes in two mtzs containing similar data and which are nearly isomorphous and computes an unbiased real-space difference map between them
 - **`matchmaps.mr`** takes in two mtzs containing similar data but which are in different spacegroups (or the same spacegroup but different crystal packing) and computes an unbiased real-space difference map between them.
 - **`matchmaps.ncs`** takes in a single mtz and computes an internal difference map across a defined non-crystallographic symmetry present in the data.

Additionally, `matchmaps` provides a utility called `matchmaps.diagnose`, which will produce a plot to help evaluate the isomorphism (or lack thereof!) between a pair of datasets. This is described in more detail [here](diagnose.md).

Last and also least, you can run `matchmaps.version` to check which version of `matchmaps` is installed.

```{eval-rst}
.. toctree::
   :maxdepth: 1
   :hidden:

   Quickstart guide <quickstart>
   Should I use matchmaps? <diagnose>
   Troubleshooting and advanced usage <troubleshooting>
   Visualizing results <visualization>
   About the algorithm <about>
   Full command-line API <cli>


```
