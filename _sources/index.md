# matchmaps's Documentation

Welcome to the docs for `matchmaps`, a python package for aligning and subtracting non-isomorphous crystallographic datasets. Do you have a pair of datasets for which an isomorphous difference maps doesn't quite work due to changes in your unit cell dimensions? If so, this is the package for you.

## Installation

The python dependencies for `matchmaps` can be installed via `pip`.
```bash
pip install matchmaps
```
I recommend that you use a package manager such as [`conda`](https://docs.conda.io/en/latest/) and always install into a fresh environment, e.g.:

```bash
conda create -n my-matchmaps-env
conda activate my-matchmaps-env
pip install matchmaps
```

### Additional dependencies

Though `matchmaps` is a python package, it relies on two pieces of external software that are not (yet!) `pip`-installable. If they do become `pip`-installable in the future, I will excitedly update this package and save you the trouble. For the time being, you will need to install:

 - [ccp4](https://www.ccp4.ac.uk/download/#os=mac)
 - [phenix](https://phenix-online.org/documentation/install-setup-run.html)

When actually using `mapreg` in the command-line, you'll need to have both ccp4 and phenix active. Doing that will look something like:
```bash
source /path/to/phenix/phenix_env.sh
/path/to/ccp4/start
```

At this point, you should be good to go! Please [file an issue on github](https://github.com/dennisbrookner/matchmaps/issues) is this is not working.

## Using `matchmaps`

`matchmaps` consists of three command-line utilities. Usage for all three is documented [here](cli.md). All three produce a real-space difference map and require a single starting model for phasing, but differ in the types of reflection data they require. Briefly:

 - **`matchmaps`** takes in two mtzs containing similar data and which are nearly isomorphous and computes an unbiased real-space difference map between them
 - **`matchmaps.mr`** takes in two mtzs containing similar data but which are in different spacegroups (or the same spacegroup but different crystal packing) and computes an unbiased real-space difference map between them.
 - **`matchmaps.mr`** takes in a single mtz and computes an internal difference map across a defined non-crystallographic symmetry present in the data.

```{eval-rst}
.. toctree::
   :maxdepth: 1
   :hidden:

   Command-line usage <cli>
   About the algorithm <about>
```
