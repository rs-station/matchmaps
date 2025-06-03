# About the `matchmaps` algorithm

If you want to learn more about the idea behind `matchmaps`, along with some examples, please check out paper in the Journal of Applied Crystallography!

> [MatchMaps: non-isomorphous difference maps for X-ray crystallography](https://journals.iucr.org/j/issues/2024/03/00/ei5112/index.html)

If you're looking for a user guide, you can find that [here](quickstart.md). But if you're looking for more details about how `matchmaps` works, read on!

## Abstract

Conformational change mediates the biological functions of macromolecules. Crystallographic measurements can map these changes with extraordinary sensitivity as a function of mutations, ligands, and time. A popular method for detecting structural differences between crystallographic datasets is the isomorphous difference map. Isomorphous difference maps combine the phases of a chosen reference state with the observed changes in structure factor amplitudes to yield a map of changes in electron density. Such maps are much more sensitive to conformational change than structure refinement is, and are unbiased in the sense that observed differences do not depend on refinement of the perturbed state. However, even modest changes in unit cell properties can render isomorphous difference maps useless. This is unnecessary. Here we describe a generalized procedure for calculating observed difference maps that retains the high sensitivity to conformational change and avoids structure refinement of the perturbed state. We have implemented this procedure in an open-source python package, MatchMaps, that can be run in any software environment supporting PHENIX and CCP4. Through examples, we show that MatchMaps "rescues" observed difference electron density maps for poorly-isomorphous crystals, corrects artifacts in nominally isomorphous difference maps, and extends to detecting differences across copies within the asymmetric unit, or across altogether different crystal forms.

## Algorithm overview

  1. Place both sets of structure factor amplitudes on a common scale using CCP4’s `SCALEIT` utility and truncate the data to the same resolution range.
  2. Generate phases for each dataset via the `phenix.refine` program. For both datasets, the OFF starting model is used, and only rigid-body refinement is permitted to prevent the introduction of model bias.
  3. Fourier-transform each set of complex structure factors into a real-space electron density map using the python packages `reciprocalspaceship` and `gemmi`.
  4. Compute the translation and rotation necessary to overlay the two rigid-body refined models. Apply this translation-rotation to the ON real-space map such that it overlays with the OFF map. These computations are carried out using `gemmi`.
  5. Subtract real-space maps voxel-wise.

## Algorithmic details

### Scaling

Scaling includes fitting both an overall scale factor and an anisotropic B-factor. `matchmaps` performs scaling via the [CCP4 `SCALEIT` utility](https://www.ccp4.ac.uk/html/scaleit.html), assisted by the [`rs-booster` `rs.scaleit` utility](https://rs-station.github.io/rs-booster/misc.html#rs-scaleit).

### Refinement

Refinement is performed via `phenix.refine`. `matchmaps` makes use of a custom `.eff` parameter template which can be found in full in the [source code](https://github.com/rs-station/matchmaps/blob/d59fa78c2f549904d0042e637262ef6c5171d355/src/matchmaps/_utils.py#L216).

By default, refinement includes bulk-solvent scaling, as this produces the best refinement results. However, in some cases, you may expect your ON data to include interesting signal far away from the OFF model, in regions outside the solvent mask. A common example of this would be if your ON data includes a bound ligand. In such situations, we recommend that bulk-solvent scaling be deactivated. You can find instructions for turning off bulk-solvent scaling (and for changing the solvent mask; see [below](#solvent-masking)) [here](quickstart.md#other-useful-options)

The user also has the option to perform rigid-body refinement on multiple different selections. For example, if your protein model contains multiple chains, those chains may move slightly relative to each other; you may not be interested in visualizing this shift in your difference map. In this case, you can specify the model selections that should be refined separately. Instuctions for doing so can be found [here](quickstart.md#other-useful-options). Note that any non-macromolecule atoms in your model will be renumbered to belong to the nearest macromolecule chain using the `phenix.sort_hetatms` utility.

More generally, this refinement can be fully customized by providing a `.eff` file. If you're interested in doing this, I recommend using the [template found in the source code](https://github.com/rs-station/matchmaps/blob/d59fa78c2f549904d0042e637262ef6c5171d355/src/matchmaps/_utils.py#L216) as a starting point. Don't hesitate to [file an issue on github](https://github.com/dennisbrookner/matchmaps/issues) if you have any issues.

### The Fourier transform

Conversion of structure factors into a real-space electron density grid is handled by `matchmaps` in python, with the help of the `reciprocalspaceship` and `gemmi` packages. This approach allows for maximum flexibility and minimal "black-box code."

`matchmaps` uses the "F-obs-filtered" column for structure factor amplitudes, and the "PH2FOFCWT" column for structure factor phases. Optionally, structure factor amplitudes can be error-weighted (using "SIGF-obs-filtered" uncertainties) via the formula described in Equation 7 here: [*reciprocalspaceship*: a Python library for crystallographic data analysis](https://scripts.iucr.org/cgi-bin/paper?S160057672100755X). If you are interested in using different MTZ columns, please let me know by [filing an issue on github](https://github.com/dennisbrookner/matchmaps/issues) and this feature could be added.

By default, both input datasets are truncated to matching resolution. Note that this is a "refine, then truncate" approach - both refinements make use of the full resolution of the dataset, and the data is only truncated afterwards. Optionally, you may provide an explicit resolution cut to be applied to both datasets; this (or error-weighting) may be useful if you believe your high-resolution reflections are noisy.

By default, the real-space voxels in the output maps are approximately 0.5 Angstrom cubes. If you're planning on visualizing your maps in Coot, I recommend keeping this default; if you're planning on visualizing your maps in PyMOL, I recommend 0.25-Angstrom spacing. Find more details [here](quickstart.md#other-useful-options).

### Real-space alignment

Conveniently, real-space alignment of the two rigid-body-refined models is very easy, because they're exactly the same! `matchmaps` use C-alphas, but any atom selection would do. Then, the transformation defining this alignment is applied to the ON real-space grid, such that it aligns with the OFF real-space grid. (You can reverse this an align the OFF grid to the ON grid if you would like, as described [here](quickstart.md#other-useful-options)). Grid alignment is handled conveniently via the `gemmi` method `interpolate_grid_of_aligned_model2`.

Note that if you have specified multiple rigid-body selections (as described [above](#refinement)), then this real-space alignment will be performed separately for each selection.

### Solvent masking

The real-space alignment performed by `matchmaps` will align the molecule chain(s) of interest, but will likely **mis-align** any symmetry-related molecules. For this reason, it is essential to mask the symmetry-related chains out of your final maps. The main difference map output by `matchmaps` is solvent masked with a 2 Angstrom radius - pretty tight.

However, as discussed [above](#refinement), you may expect your ON data to include signal far away from your OFF model. To account for this possibility, `matchmaps` produces a second difference map with a more generous solvent-masking radius. This more generous radius defaults to 5 Angstroms and [can be changed](quickstart.md#other-useful-options). The file containing this generously-masked difference map will have `_unmasked` at the end of its name; `matchmaps` output files are described [here](quickstart.md#important-map-outputs).

Note that real-space normalization of the two maps is critical for producing a "balanced" difference map. This normalization is performed not on across the unit cell or the asymmetric unit, but rather on the remaining non-zero voxels after applying the more generous solvent mask.

## Variants of `matchmaps`

In addition to the core `matchmaps` utility, two additional command-line functions are available: `matchmaps.mr` and `matchmaps.ncs`. The algorithm is modified slightly for each utility.

### `matchmaps.mr`

`matchmaps.mr` can support two input reflection files that are in different crystal packings or spacegroups. Accordingly, step 1 above (scaling) is replaced with a round of molecular replacement (using `phenix.phaser`) wherein the OFF starting model is used as a molecular replacement solution for the ON reflections. The algorithm proceeds essentially identically from this point.

One small but important change in `matchmaps.mr` is that all ordered water molecules in the OFF starting model are discarded. This is important, as the ordered water molecules in one spacegroup/crystal packing are often not appropriate in a different spacegroup/crystal packing.

### `matchmaps.ncs`

`matchmaps.ncs` creates internal difference maps across a non-crystallographic symmetry operation. This function requires a number of modifications to the core algorithm:
 - Step 1 is omitted; there is only one input reflection file, so no scaling is necessary.
 - Step 2 is optional. If phases are not included, `phenix.refine` is run to generate phases. However, if phases are provided, they are used and refinement is skipped.
 - In step 4, there is only one rigid-body refined model. Alignment is performed on the two regions of the model related by NCS (e.g., two different protein chains).

## Questions or thoughts?
If you have any questions, thoughts, or suggestions regarding this algorithm and software, we would love to hear from you! The easiest way to get in touch is by [opening an issue on GitHub](https://github.com/rs-station/matchmaps/issues)
