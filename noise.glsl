/////////////// K.jpg's Simplex-Style Re-oriented 4Point BCC Noise ///////////////
///////////////////// Experimantal non-spherical kernel test ////////////////////

// Inspired by Stefan Gustavson's noise
vec4 permute(vec4 t) {
    return t * (t * 34.0 + 133.0);
}

// Gradient set is a normalized expanded rhombic dodecahedron
vec3 grad(float hash) {
    
    // Random vertex of a cube, +/- 1 each
    vec3 cube = mod(floor(hash / vec3(1.0, 2.0, 4.0)), 2.0) * 2.0 - 1.0;
    
    // Random edge of the three edges connected to that vertex
    // Also a cuboctahedral vertex
    // And corresponds to the face of its dual, the rhombic dodecahedron
    vec3 cuboct = cube;
    cuboct[int(hash / 16.0)] = 0.0;
    
    // In a funky way, pick one of the four points on the rhombic face
    float type = mod(floor(hash / 8.0), 2.0);
    vec3 rhomb = (1.0 - type) * cube + type * (cuboct + cross(cube, cuboct));
    
    // Expand it so that the new edges are the same length
    // as the existing ones
    vec3 grad = cuboct * 1.22474487139 + rhomb;
    
    // To make all gradients the same length, we only need to shorten the
    // second type of vector. We also put in the whole noise scale constant.
    // The compiler should reduce it into the existing floats. I think.
    grad *= (1.0 - 0.042942436724648037 * type) * 2.0; //TODO we can find a better normalization constant probably
    
    return grad;
}

// BCC lattice split up into 2 cube lattices
float bccNoiseBase(vec3 X) {
    
    // First half-lattice, closest edge
    vec3 v1 = round(X);
    vec3 d1 = X - v1;
    vec3 score1 = abs(d1);
    vec3 dir1 = step(max(score1.yzx, score1.zxy), score1);
    vec3 v2 = v1 + dir1 * sign(d1);
    vec3 d2 = X - v2;
    
    // Second half-lattice, closest edge
    vec3 X2 = X + 144.5;
    vec3 v3 = round(X2);
    vec3 d3 = X2 - v3;
    vec3 score2 = abs(d3);
    vec3 dir2 = step(max(score2.yzx, score2.zxy), score2);
    vec3 v4 = v3 + dir2 * sign(d3);
    vec3 d4 = X2 - v4;
    
    // Gradient hashes for the four points, two from each half-lattice
    vec4 hashes = permute(mod(vec4(v1.x, v2.x, v3.x, v4.x), 289.0));
    hashes = permute(mod(hashes + vec4(v1.y, v2.y, v3.y, v4.y), 289.0));
    hashes = mod(permute(mod(hashes + vec4(v1.z, v2.z, v3.z, v4.z), 289.0)), 48.0);
    
    // One 1D bump function between each pair of opposing faces on the lattice neighborhood figure.
    // All multiplied together, then squared at the end.
    // Should be possible in any number of dimensions.
    // May not achieve as high a degree of continuity.
    vec3 d1c = d1.xxy + d1.yzz; vec3 d2c = d2.xxy + d2.yzz;
    vec3 d3c = d3.xxy + d3.yzz; vec3 d4c = d4.xxy + d4.yzz;
    vec3 d1d = d1.xxy - d1.yzz; vec3 d2d = d2.xxy - d2.yzz;
    vec3 d3d = d3.xxy - d3.yzz; vec3 d4d = d4.xxy - d4.yzz;
    d1c = 1.0 - d1c * d1c; d2c = 1.0 - d2c * d2c;
    d3c = 1.0 - d3c * d3c; d4c = 1.0 - d4c * d4c;
    d1d = 1.0 - d1d * d1d; d2d = 1.0 - d2d * d2d;
    d3d = 1.0 - d3d * d3d; d4d = 1.0 - d4d * d4d;
    vec4 b = vec4(d1c.x * d1c.y * d1c.z, d2c.x * d2c.y * d2c.z, d3c.x * d3c.y * d3c.z, d4c.x * d4c.y * d4c.z);
    b *= vec4(d1d.x * d1d.y * d1d.z, d2d.x * d2d.y * d2d.z, d3d.x * d3d.y * d3d.z, d4d.x * d4d.y * d4d.z);
    vec4 bb = b * b;
    
    // Gradient extrapolations
    vec3 g1 = grad(hashes.x); vec3 g2 = grad(hashes.y);
    vec3 g3 = grad(hashes.z); vec3 g4 = grad(hashes.w);
    vec4 extrapolations = vec4(dot(d1, g1), dot(d2, g2), dot(d3, g3), dot(d4, g4));
    
    return dot(bb, extrapolations);
}

// Use this if you don't want Z to look different from X and Y
float bccNoiseClassic(vec3 X) {
    
    // Rotate around the main diagonal. Not a skew transform.
    return bccNoiseBase(dot(X, vec3(2.0/3.0)) - X);
}

// Use this if you want to show X and Y in a plane, and use Z for time, etc.
float bccNoisePlaneFirst(vec3 X) {
    
    // Rotate so Z points down the main diagonal. Not a skew transform.
    mat3 orthonormalMap = mat3(
        0.788675134594813, -0.211324865405187, -0.577350269189626,
        -0.211324865405187, 0.788675134594813, -0.577350269189626,
        0.577350269189626, 0.577350269189626, 0.577350269189626);
    
    return bccNoiseBase(orthonormalMap * X);
}


//////////////////////////////// End noise code ////////////////////////////////
