# Noise Extras

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
