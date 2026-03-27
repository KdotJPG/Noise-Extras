# Noise Extras

## 2D-tileable 3D Simplex-type noise

#### OpenSimplex2S_ImprovedXZPlanes_TileableXZ.java

3D OpenSimplex2S noise that achieves periodicity in a 2D rectangle by slightly distorting the grid.

![Image composed of four seamless quadrants of repeating smooth noise](https://github.com/user-attachments/assets/8ff13ed7-94fe-4573-b7e7-84a65b290058)

## Simplex-Gabor Noise
#### SimlexGaborNoise.java
Simplex-style noise modified to use a Gabor-like sinusoid kernel with variable direction and magnitude.

![Smooth noise image with slightly varied diagonal waves packed tightly together](https://github.com/user-attachments/assets/4f69f61d-99cf-4d95-a02b-c5e7f6c77f5e)

## Simplex-Gabor Noise (Isotropic version)
#### SimlexGaborNoiseIso.java
Simplex-style noise modified to use a Gabor-like sinusoid kernel with variable magnitude and randomized direction.

![Smooth noise image with highly variational packed wave angles](https://github.com/user-attachments/assets/2f9102c6-4507-4d05-ad41-5f687d44d00b)

https://www.shadertoy.com/view/Wl3Gzl

## Simplex-Style Noise with Neighborhood-Figure Kernel
#### NeighborFigureKernelSimplexStyleNoise.glsl
Experimental simplex-type noise with a kernel that conforms to the lattice neighborhood figure instead of a sphere, to reduce the bubbly appearance without introducing discontinuities or accounting for more lattice points. May have a lower degree of differentiability.

![Noise image more relaxed than 3D simplex noise, but not as much so as OpenSimplex2S](https://github.com/user-attachments/assets/0a9d82d3-646a-4021-987a-a871261d3497)

## Exponential Simplex
#### OpenSimplex2S_Exp.java
OpenSimplex2S modified to have an exponential-esque distribution as inspired by http://jcgt.org/published/0004/02/01/. Preview image generated using expBase=64 and expMin=0.1. Left-hand side is a single octave, and right hand side is fBm / summed octaves with lacunarity=(1/persistence)=sqrt(2).

![Smooth simplex-reminiscent noise with wide variation in extrema magnitude](https://github.com/user-attachments/assets/23033874-51eb-4c69-97a8-2bd7d2eae8c4)

## Simplex-Cellular Noise
#### SimplexCellularNoise.java
Cellular noise using the lattice of Simplex noise. May or may not look nicer than traditional cellular noise. Cellular noise on a square grid is usually already less obvious than gradient/value noise on a square grid.

![Randomized honeycomb-like image reminiscent of familiar cellular noise](https://github.com/user-attachments/assets/dcfa94e9-8e4b-4231-9047-41eeeda8cf29)

## Multi-seedable Noise

#### Multi-seedable Noise/*

Simplex, domain-rotated Perlin, and Cellular noise implementations capable of [evaluating multiple noise seeds at once](https://noiseposti.ng/posts/2024-01-20-Evaluating-Multiple-Noise-Seeds-At-Once.html) for improved performance.

## Vector-valued Gradient Noise (Orthogonal Gradient Matrix)

#### Vector-Valued Gradient Noise/**

OpemSimplex2S implementations that output vectors with more spherical probability density by sampling orthogonal gradient vectors at each vertex. Useful for domain warping and other directional needs.

## Vector-valued Gradient Noise (Gradient Outer Product)

#### OpenSimplex2S_3D_WithVectorOutputNoise.java

OpemSimplex2S implementation that outputs vectors with more spherical probability density by sampling a second vector at each vertex to determine its output contribution direction. Also useful for domain warping and other directional needs.

## Cached Column Noise

#### CachedColumnNoise/**

3D OpenSimplex2, OpenSimplex2S, and domain-rotated Perlin implementations that improve the performance of vertical column sampling (y) for fixed horizontal (x,z) coordinates.

