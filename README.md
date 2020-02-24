# Noise Extras

## Exponential Simplex
#### OpenSimplex2S_Exp.java
OpenSimplex2S modified to have an exponential-esque distribution as inspired by http://jcgt.org/published/0004/02/01/. Preview image generated using expBase=64 and expMin=0.1. Left-hand side is a single octave, and right hand side is fBm / summed octaves with lacunarity=(1/persistence)=sqrt(2).

![Simplex-Cellular Noise](images/OpenSimplex2S_Exp.png?raw=true)

## Simplex-Cellular Noise
#### SimplexCellularNoise.java
Cellular noise implemented using the lattice of Simplex noise. May or may not look nicer than traditional cellular noise. Cellular noise on a square grid is usually already less obvious than gradient/value noise on a square grid.

![Simplex-Cellular Noise](images/SimplexCellularNoise.png?raw=true)

## Simplex-Gabor Noise
#### SimlexGaborNoise.java
Simplex-style noise modified to use a Gabor-like sinusoid kernel with variable direction and magnitude.

![Simplex-Gabor Noise](images/SimplexGaborNoise.png?raw=true)

## Simplex-Gabor Noise (Isotropic version)
#### SimlexGaborNoiseIso.java
Simplex-style noise modified to use a Gabor-like sinusoid kernel with variable magnitude and randomized direction.

![Simplex-Gabor Noise Isotropic](images/SimplexGaborNoiseIso.png?raw=true)

## Simplex-Style Noise with Neighborhood-Figure Kernel
#### NeighborFigureKernelSimplexStyleNoise.glsl
Experimental simplex-style noise with a kernel that conforms to the lattice neighborhood figure instead of a sphere, to reduce the bubbly appearance, but without introducing discontinuities or accounting for more lattice points. May have a lower degree of differentiability.

![Space Tiling Kernel Simplex Style Noise](images/NeighborFigureKernelSimplexStyleNoise.png?raw=true)

https://www.shadertoy.com/view/Wl3Gzl
