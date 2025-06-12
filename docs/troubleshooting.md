# Troubleshooting and advanced usage

If you've read through the [quickstart guide](quickstart.md) and you still have questions or issues, read on for more information.

## Solvent masking and bulk-solvent scaling
The goal of `matchmaps` is to visualize difference density around the protein molecule. However, the question remains: what does "around" mean exactly? To address this ambiguity, `matchmaps` writes out two difference maps. The main difference map is solvent masked pretty strictly and only contains signal within 2 A of the protein model. The second difference map, which will include the suffix `_unmasked`, is solvent masked more laxly. The default masking radius for the `_unmasked` map is 5 A, but you can change this value using the `--unmasked-radius` flag. This might be useful if you expect to see signal farther away from your protein model, e.g. a bound ligand.

Another important parameter for visualizing difference density far from your protein model is bulk-solvent scaling (BSS). BSS is a typical feature of crystallographic refinement, and including it will typically result in better refinement statistics and nicer maps. By default, `phenix.refine` will perform BSS prior to rigid-body refinement. However, if you expect to see signal far away from your protein model, you may find that BSS will "flatten" or otherwise alter this signal. You can disable BSS using the `--no-bss` flag. I strongly recommend this when analyzing an apo/bound pair via `matchmaps`.

Also, for using `matchmaps` with bound ligands, see the `--on-as-stationary` flag [below](#miscellaneous-useful-options).

### Symmetry-related molecules

A side effect of the `matchmaps` real-space alignment approach is that while the ON and OFF models will end up aligned, the *symmetry mates* of the ON and OFF models will necessarily end up *misaligned*. This is a big reason that solvent masking is so important: otherwise, difference density for the misaligned symmetry mates will dominate your maps. Unfortunately, even after strict solvent masking, your final difference map is likely to contain this artifactual signal. The good news is that this artifactual signal is pretty easy to identify and disregard. The bad news is that if you're using a feature like Coot's "find difference peaks," you're likely to see these a lot.

In the medium-to-long term, I aim for the next generation of `matchmaps` masking approaches to resolve this issue more elegantly. If you have ideas about this, or if you just want to bug me about it, I would defintely welcome an [issue](https://github.com/rs-station/matchmaps/issues) and/or pull request on GitHub! In the shorter term, if this signal is becoming a large issue, you could try applying an even stricter solvent radius (say, 1 A) to the `--unmasked-radius` flag and see if that helps.

## Resolution cuts and error weighting
By default, `matchmaps` will simply truncate the two input reflection files to equal resolution. However, if you expect the highest-resolution reflections to still be noisy after this, or if your difference maps look very noisy, you might consider cutting resolution even further. You can do this by providing a resolution cut to the `--dmin` flag.

Alternatively, rather than just truncating at a particular resolution, you can apply an error weighting to the reflections. Error weighting is performed immediately prior to the Fourier transform. Weights are computed via the following formula:

```{eval-rst}
.. math::

   \frac{1}{1 + \alpha \frac{{(\sigma F)}^2}{{<\sigma F>}^2}}
```

where $\sigma F$ is the error estimate for the reflection, $<\sigma F>$ is the average error estimate across all reflections, and $\alpha$ is a weighting parameter. By default, $\alpha = 0$, e.g. no weighting. You can use the `--alpha` flag to supply your own value and thus include error-weighting.

## Multiple protein chains

You may have multiple protein chains in your model. This presents you with some interesting opportunities for difference maps, but also some potential headaches.

### Refining chains individually

One possibility is that your `matchmaps` difference map contains essentially no signal for one protein chain, but strong signal throughout the other chain indicating a global motion. This signal is "real" -- it tells you that the relative packing of these two chains together is different in your two datasets -- but it probably doesn't make for a very useful map. Instead, you might consider using the `--rbr-selections` flag to rigid-body refine each chain separately. If your chains are A and B, then you would use the flag as `--rbr-selections A B`. `matchmaps` will then produce a difference map specific to each chain. This flag is admittedly a little finnicky; please [file an issue](https://github.com/rs-station/matchmaps/issues)  if you have any trouble.

### Comparing chains to each other

Assuming your data is some sort of homo-multimer / non-crystallographic symmetry, another option is to use `matchmaps.ncs` to compare the protein chains within a dataset against each other. If you have multiple datasets, of course, you could still run `matchmaps.ncs` on each dataset and see how the difference map changes. See the full documentation for `matchmaps.ncs` [here](cli.md#matchmaps-ncs)

## Miscellaneous useful options

 - `--on-as-stationary`: The `matchmaps` algorithm always involves an alignment in real-space of the "on" and "off" maps. By default, the "off" map is stationary, and the "on" map is moved. This is typically desired, such that everything lines up with your "off" structural model. However, say that your structures are "apo" and "bound", and you would like to line up your maps with a "bound" structure (which you never have to supply to `matchmaps`!). In this case, you could use the `--on-as-stationary` flag.
 - `--spacing`: This flag defines the approximate size of the voxels in your real-space maps. The default (0.5 A) is fine for most purposes. For making figures in PyMOL, you might want finer spacing (0.25 A or so); this comes at a cost of much larger file size. If your computer/coot is being really slow, you could consider increasing the spacing.
 - `--verbose`: Use this option to print out to the terminal all of the log output from CCP4 and phenix. This is disabled by default because it's very annoying, but it can be useful for debugging purposes.
 - `--eff`: The `matchmaps` source code contains a hard-coded `.eff` template which is modified, written to a file, and passed to `phenix.refine`. For most cases, this `.eff` template should do the trick. However, if there's something specific that you would like phenix to do, you can pass your own custom `.eff` template to `matchmaps` via the `--eff` flag. There is a lot of potential for error here, because the code has very specific expectations for what the `.eff` file contains. If you're interested in trying this out, I would recommend that you use the [`.eff` template in the source code](https://github.com/rs-station/matchmaps/blob/7531ff1b13da91b01ede273fa9b1f5a99d72a5ca/src/matchmaps/_utils.py#L227) as a starting point. Don't hesitate to [file an issue on GitHub](https://github.com/rs-station/matchmaps/issues) if anything isn't working. Depending on what you're trying to do, we may decide to try and implement your desired functionality directly, so you don't need to provide a custom `.eff` in the long run.
