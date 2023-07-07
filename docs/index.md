# matchmaps's Documentation

Welcome to the docs for `matchmaps`, a python package for aligning and subtracting non-isomorphous crystallographic datasets. Do you have a pair of datasets for which an isomorphous difference maps doesn't quite work due to changes in your unit cell dimensions? If so, this is the package for you.

The [quickstart guide](quickstart.md) contains installation instructions and sample usage of the basic `matchmaps` utility, along with other useful tips.

## The `matchmaps` functions in brief

`matchmaps` consists of three command-line utilities. All options for all three utilites are documented [here](cli.md). All three produce a real-space difference map and require a single starting model for phasing, but differ in the types of reflection data they require. Briefly:

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
