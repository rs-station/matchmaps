# Should I use `matchmaps`? Diagnosing non-isomorphism

You have two crystallographic datasets that you'd like to compare with an Fo-Fo difference map. Should you use `matchmaps`, or will an isomorphous difference map do the trick? To help answer this question, `matchmaps` provides the utility `matchmaps.diagnose`, which creates a plot similar to the following (Figure 1a from the [MatchMaps paper](https://journals.iucr.org/j/issues/2024/03/00/ei5112/index.html))

![Figure 1a: resolution dependence of inter-dataset correlation](images/figure.jpg)

## Interpreting this plot

The plot above shows three curves for each dataset. Two of them report on a 

## Using `matchmaps.diagnose`