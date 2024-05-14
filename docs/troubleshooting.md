# Troubleshooting and advanced usage

If you've read through the [quickstart guide](quickstart.md) and you still have questions or issues, read on for more information.

<!-- ## Some useful `matchmaps` options
Below is an overview of some parameters which aren't strictly necessary, but may come in handy. Note that the full list of parameters for each `matchmaps` utility can be found [here](cli.md). -->

## `matchmaps` with bound ligands.
The `matchmaps` algorithm always involves an alignment in real-space of the "on" and "off" maps. By default, the "off" map is stationary, and the "on" map is moved. This is typically desired, such that everything lines up with your "off" structural model. However, say that your structures are "apo" and "bound", and you would like to line up your maps with a "bound" structure (which you never have to supply to `matchmaps`!). In this case, you could use the `--on-as-stationary` flag.

 - `--unmasked-radius`: How far away from the protein model do you expect to see difference signal? Use this flag to change the behavior of the `_unmasked` difference map output to show more or less signal far from the protein. Defaults to 5 A. See [below](#important-map-outputs) for more details.
 - `--no-bss`: If included, skip the bulk solvent scaling step of phenix.refine. Like `--unmasked-radius`, this option may be useful in situtations where you expect signal far away from your protein model. For example, bulk solvent scaling may "flatten" or otherwise alter signal for an unmodeled bound ligand.

## My difference map is noisy
 - `--dmin`: The input `mtz` files are truncated to equal resolution by default. If you would like, the `mtz`s may be truncated even more stringently. This could be useful if you expect your highest-resolution reflections to be more noise than signal.
 - `--alpha`: Rather than just truncating at a particular resolution, you might want to error-weight your reflections. You can 

## I have multiple protein chains in my model

 - `--rbr-selections`: When doing rigid-body refinement, refine as multiple explicitly defined rigid bodies rather than a single rigid body containing everything. This flag is admittedly a little finnicky; please [file an issue](https://github.com/rs-station/matchmaps/issues)  if you have any trouble.


## Miscellaneous useful options

 - `--spacing`: This flag defines the approximate size of the voxels in your real-space maps. The default (0.5 A) is fine for most purposes. For making figures in PyMOL, you might want finer spacing (0.25 A or so); this comes at a cost of much larger file size. If your computer/coot is being really slow, you could consider increasing the spacing.
 - `--verbose`: Use this option to print out to the terminal all of the log output from CCP4 and phenix. This is disabled by default because it's very annoying, but it can be useful for debugging purposes.