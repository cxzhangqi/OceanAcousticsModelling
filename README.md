# OceanAcousticsModelling
This package is an implementation of a number of textbooks I study in ocean acoustics modelling.

Notes:
* I still haven't figured out how to produce Julia packages yet.

## Bottom Loss
I've started off simple, replicating the bottom loss curves in Jensen et al [[1]](#JensenEtAl).

The bottom loss is dependent on the complex-valued sound speed which includes the volume attenuation in its imaginary component.
![](img/BottomLoss_Parameters.png)

It was then simple enough to replicate the bottom loss for various bottom sediments.
![](img/BottomLoss_Types.png)

## Ray Tracing
The Eikonal equation is solved using Julia's [DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl) using the time variable in the solvers for the arc-length.

This first example shows a Gaussian bathymetry and sinusoidal altimetry. The sound speed profile is also range-dependent as a parabola with maximums at the ocean boundaries, which shifts with the changing boundary depths.

![](img/../RayTrace_FirstExample.png)

## References
> <a name="JensenEtAl">[1]</a> Jensen, F. B., Kuperman, W. A., Porter, M. B., & Schmidt, H. (2011). Computational ocean acoustics. Springer Science & Business Media.