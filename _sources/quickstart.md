# Getting started with `matchmaps`

On this page, we'll explore how to use the basic `matchmaps` utility and examine its outputs.  Full documation of all options for all three command-line utilities can be found [here](cli.md) or by typing the command plus `--help` into the command line.

## Installation

`matchmaps` and its python dependencies can be installed via `pip`:
```bash
pip install matchmaps
```
I recommend that you use a package manager such as [`conda`](https://docs.conda.io/en/latest/) and always install into a fresh environment, e.g.:

```bash
conda create -n my-matchmaps-env python=3.9
conda activate my-matchmaps-env
pip install matchmaps
```

If you would like to use the latest version of `matchmaps` that hasn't yet made it to PyPI, you can alternatively install directly from GitHub:
```bash
pip install git+https://github.com/rs-station/matchmaps.git
```

### Additional dependencies

Though `matchmaps` is a python package, it relies on two pieces of external software that are not (yet!) `pip`-installable. If they do become `pip`-installable in the future, I will excitedly update this package and save you the trouble. For the time being, you will need to install:

 - [ccp4](https://www.ccp4.ac.uk/download/#os=mac)
 - [phenix](https://phenix-online.org/documentation/install-setup-run.html)

When actually using `matchmaps` in the command-line, you'll need to have both ccp4 and phenix active. Doing that will look something like:
```bash
source /path/to/phenix/phenix_env.sh
/path/to/ccp4/start
```

At this point, you should be good to go! Please [file an issue on github](https://github.com/dennisbrookner/matchmaps/issues) is this is not working.

## Input files

To run `matchmaps`, you will need:
 - one `.pdb` (or `.cif`) file containing a refined structural model corresponding to your "off" data.
 - two `.mtz` (or `.cif`) files corresponding to your "on" and "off" data respectively.

You will also need to know the names of the columns in these `mtz`/`cif`s containing your observed structure factor amplitudes and uncertainties. Depending on what software you used to produce these files, this may be something like `FP`/`SIGFP`, `Fobs`/`SIGFobs`, or similar. If you don't know these off-hand and your input is an `.mtz` file, you can figure it out using [`reciprocalspaceship`](https://rs-station.github.io/reciprocalspaceship/)'s `rs.mtzdump` utility, which is installed along with `matchmaps`. You can do this right in the command-line as:
```bash
rs.mtzdump mymtz.mtz
```
which will print a summary of the contents of the `.mtz` file.

Finally, if your structure contains any ligands or solvent for which a restraint file (`.cif`) is needed, you will need those as well.

#### A note on "on" and "off" data

Throughout this documentation, we will assume to be working with a pair of datasets that differ by some perturbation. These datasets could be apo/bound, light/dark, warm/cold, WT/mutant, etc. We will generalize these perturbations as representing either "off" or "on" data. Importantly, these datasets are not created equal! Your "off" data should be the data for which you have refined structural coordinates. For your "on" data, you do not need to provide a corresponding structure. This is the data which you hope to visualize in a model-bias-free manner.

## Basic usage

If the above files are in your current directory, and you would like to write output files into your current directory, then you only need three arguments: `--mtzoff`, `--mtzon`, and `--pdboff`. For example:

```bash
matchmaps --mtzoff apo_data.mtz Fobs SIGFobs --mtzon bound_data.mtz Fobs SIGFobs --pdboff apo.pdb
```

Any necessary ligand restraints can be added via the `--ligands` flag, e.g.:

```bash
matchmaps --mtzoff apo_data.mtz Fobs SIGFobs \
    --mtzon bound_data.mtz Fobs SIGFobs \
    --pdboff apo.pdb \
    --ligands weird_solvent_1.cif weird_solvent_2.cif
```
 
If you'd like read or write files from somewhere other than your current directory, you can! There are three ways to specify input files:
 - Provide relative paths directly for all input files
 - Provide only file names, and add the `--input-dir` option to specify where those files live. If you do this, the same `--input-dir` will be preprended to all filenames, so your files should all live in the same place.
 - Provide absolute paths to input files. For any input supplied as an absolute path, the `--input-dir` will be ignored.

To direct output files to a specific directory, use the `--output-dir` flag. By default, all of the temporary files created by `matchmaps` will be deleted when the program finishes. Only the `.map` files described below are kept. If you would like to keep all files, you may additionally supply a `--keep-temp-files` directory, which will be created inside of `--output-dir`.

### Examples 

Supply an input directory:
```bash
matchmaps --mtzoff apo_data.mtz Fobs SIGFobs \
    --mtzon bound_data.mtz Fobs SIGFobs \
    --pdboff apo.pdb \
    --input-dir analysis/matchmaps \
    --output-dir ../data/myproject
```

Supply relative inputs; keep all files:
```bash
matchmaps --mtzoff input_files/apo_data.mtz Fobs SIGFobs \
    --mtzon input_files/bound_data.mtz Fobs SIGFobs \
    --pdboff different_dir_with_input_files/apo.pdb \
    --output-dir ../data/myproject \
    --keep-temp-files mmfiles
```

Supply a mix of relative and absolute inputs:
```bash
matchmaps --mtzoff apo_data.mtz Fobs SIGFobs \
    --mtzon bound_data.mtz Fobs SIGFobs \
    --pdboff /complicated/absolute/path/to/apo.pdb \
    --input-dir input_files \
    --output-dir ../data/myproject \
```

### Running `matchmaps` from a script
After your `matchmaps` run completes successfully, it will write out a file (called `run_matchmaps.sh` unless you specify a different name to the `--script` flag) which can be run to reproduce exactly the same command. Note that, to ensure the compatibility of input/output paths, this script is written to your ***current working directory***.

If you'd then like to run `matchmaps` again with slightly different parameters, you can use this script as a starting point. No need to remember exactly which parameters you used the first time!

## Other useful options

 - `--on-as-stationary`: The `matchmaps` algorithm always involves an alignment in real-space of the "on" and "off" maps. By default, the "off" map is stationary, and the "on" map is moved. This is typically desired, such that everything lines up with your "off" structural model. However, say that your structures are "apo" and "bound", and you would like to line up your maps with a "bound" structure (which you never have to supply to `matchmaps`!). In this case, you could use the `--on-as-stationary` flag.
 - `--dmin`: The input `mtz` files are truncated to equal resolution by default. If you would like, the `mtz`s may be truncated even more stringently.
 - `--unmasked-radius`: How far away from the protein model do you expect to see difference signal? Use this flag to change the behavior of the `_unmasked` difference map output to show more or less signal far from the protein. Defaults to 5 A. See [below](#important-map-outputs) for more details.
 - `--no-bss`: If included, skip the bulk solvent scaling step of phenix.refine. Like `--unmasked-radius`, this option may be useful in situtations where you expect signal far away from your protein model. For example, bulk solvent scaling may "flatten" or otherwise alter signal for an unmodeled bound ligand.
 - `--spacing`: This flag defines the approximate size of the voxels in your real-space maps. The default (0.5 A) is fine for most purposes. For making figures in PyMOL, you might want finer spacing (0.25 A or so); this comes at a cost of much larger file size. If your computer/coot is being really slow, you could consider increasing the spacing.
 - `--verbose`: Use this option to print out to the terminal all of the log output from CCP4 and phenix. This is disabled by default because it's very annoying, but it can be useful for debugging purposes.
 - `--rbr-selections`: When doing rigid-body refinement, refine as multiple explicitly defined rigid bodies rather than a single rigid body containing everything. This flag is admittedly a little finnicky; please [file an issue](https://github.com/dennisbrookner/matchmaps/issues) if you have any trouble.

Note that most of the command-line options have short and long versions, e.g. `-i` vs. `--input-dir`. For clarity, the long names have been used exclusively on this page. The [full documentation](cli.md) lits all short and long options.

## Output files

Below is a quick tour of the output files that `matchmaps` will produce and what you might want to do with them.

### Important `.map` outputs

Let's assume that your input files are called `off.mtz` and `on.mtz`. The following files created by `matchmaps` may be of interest:

 - `on_minus_off.map`: This is your difference map! It should contain positive and negative signal in the vicinity (>= 2 Angstroms) of your protein model.
 - `on_minus_off_unmasked.map`: The same as the previous difference map, but before the 2 A solvent mask was applied. This file can be useful if you're expecting signal (perhaps a bound ligand) far away from your protein model. By default, this map contains signal up to 5 A away from the protein model; this radius can be changed with the `--unmasked-radius` flag for `matchmaps` and `matchmaps.mr` utilites. 
 - `on.map` / `off.map`: The real-space maps which are subtracted to produce the above difference maps. It is a good idea to open these files and inspect them. They should be generally aligned in space. Any interesting signal you expect to see in a difference map may also be apparent by inspecting these maps. Remember that both of these maps were computed using the "off" model, so structural features of the "off" data are likely to be more prominent.
 - `on_before.map` / `off_before.map`: The real-space maps, prior to alignment. These maps may be useful a) if you're curious how much alignment was necessary, or b) to troubleshoot where in the pipeline something went wrong.

Additionally, `matchmaps` produces ~15 other files which by default are deleted after the program finishes. If you would like to keep these files, you can use the `--keep-temp-files` flag described above.

Note that if you re-run `matchmaps` into the same output directory, the `.map` output files ***will*** be overwritten. I recommend directing each `matchmaps` run in to a unique, informatively-named output directory.

### Working with `.map` files in Coot

If you have experience with cryo-EM, you're probably familiar with `.map`/`.CCP4`/`.mrc` real-space map files. But if you're a crystallographer, you might not have worked with these at all! Crystallographers typically just open `.mtz` files directly in a software (like coot or PyMOL) that does the Fourier transform and computes the map "on the fly". But here you are, working with `matchmaps`, which by design produces real-space outputs.

`.map` files should be opened in Coot via "File > Open Map...":

![Open Map...](images/openmap.png)

If the map you are opening is a difference map, it is essential that you check the option "Is Difference Map" in the bottom left corner of the "Select Map..." dialog. If you don't do this here, there is no way to do it later; you'll have to close and reopen the map.

![Is difference map](images/isdifferencemap.png)

Happy difference mapping! From here, working with a `.map` file should be no different than typical work in Coot with an `.mtz` file.
