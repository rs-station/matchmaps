# About the `matchmaps` algorithm

If you want to learn more about the idea behind `matchmaps`, along with some examples, please check out our new pre-print!  
  
> [MatchMaps: Non-isomorphous difference maps for X-ray crystallography](https://www.biorxiv.org/content/10.1101/2023.09.01.555333v1) 

## Excerpts from the pre-print

### Abstract
Conformational change mediates the biological functions of proteins. Crystallographic measurements can map these changes with extraordinary sensitivity as a function of mutations, ligands, and time. The isomorphous difference map remains the gold standard for detecting structural differences between datasets. Isomorphous difference maps combine the phases of a chosen reference state with the observed changes in structure factor amplitudes to yield a map of changes in electron density. Such maps are much more sensitive to conformational change than structure refinement is, and are unbiased in the sense that observed differences do not depend on refinement of the perturbed state. However, even minute changes in unit cell dimensions can render isomorphous difference maps useless. This is unnecessary. Here we describe a generalized procedure for calculating observed difference maps that retains the high sensitivity to conformational change and avoids structure refinement of the perturbed state. We have implemented this procedure in an open-source python package, MatchMaps, that can be run in any software environment supporting PHENIX and CCP4. Through examples, we show that MatchMaps “rescues” observed difference electron density maps for near-isomorphous crystals, corrects artifacts in nominally isomorphous difference maps, and extends to detecting differences across copies within the asymmetric unit, or across altogether different crystal forms.

### Algorithm description

  1. Place both sets of structure factor amplitudes on a common scale using CCP4’s `SCALEIT` utility and truncate the data to the same resolution range.
  2. Generate phases for each dataset via the `phenix.refine` program. For each dataset, the OFF starting model is used, and only rigid-body refinement is permitted to prevent the introduction of model bias.
  3. Fourier-transform each set of complex structure factors into a real-space electron density map using the python packages `reciprocalspaceship` and `gemmi`.
  4. Compute the translation and rotation necessary to overlay the two rigid-body refined models. Apply this translation-rotation to the ON real-space map such that it overlays with the OFF map. These computations are carried out using `gemmi`.
  5. Place both real-space maps on a common scale.
  6. Subtract real-space maps voxel-wise.
  7. Apply a solvent mask to the final difference map.

## `matchmaps` variants

In addition to the core `matchmaps` utility, two additional command-line functions are available: `matchmaps.mr` and `matchmaps.ncs`. The algorithm is modified slightly for each utility.

### `matchmaps.mr`

`matchmaps.mr` can support two input reflection files that are in different crystal packings or spacegroups. Accordingly, step 1 above is replaced with a round of molecular replacement (using `phenix.phaser`) whereing the OFF starting model is used as a molecular replacement solution for the ON reflections. The algorithm proceeds essentially identically from this point.

One small but important change in `matchmaps.mr` is that all ordered water molecules in the OFF starting model are discarded. This is important, as the ordered water molecules in one spacegroup/crystal packing are often not appropriate in a different spacegroup/crystal packing.

### `matchmaps.ncs`

`matchmaps.ncs` creates internal difference maps across a non-crystallographic symmetry operation. This function requires a number of modifications to the core algorithm:
 - Step 1 is omitted; there is only one input reflection file, so no scaling is necessary.
 - Step 2 is optional. If phases are not included, `phenix.refine` is run to generate phases. However, if phases are provided, they are used and refinement is skipped.
 - In step 4, there is only one rigid-body refined model. Alignment is performed on the two regions of the model related by NCS (e.g., two different protein chains).

## Questions or thoughts?
If you have any questions, thoughts, or suggestions regarding this algorithm and software, we would love to hear from you! The easiest way to get in touch is by [opening an issue on GitHub](https://github.com/rs-station/matchmaps/issues)