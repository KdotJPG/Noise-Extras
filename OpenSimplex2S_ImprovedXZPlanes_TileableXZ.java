public class OpenSimplex2S_ImprovedXZPlanes_TileableXZ {
    
    // For larger repeat periods combined with higher frequencies, it is possible that pre-multiplying the primes
    // can create overflow in the 64-bit integer coordinate position values, breaking tileability. For more modest
    // use cases, you can enable this for a ~3% performance improvement (according to rudimentary tests on my end).
    // Note that, after moving the pre-prime multiplication to after the trailing zero removal step, I found that
    // cases which break when PREPRIME=true but tile when it's false, became more elusive but were still possible.
    // For example, frequency=0.125 and repeatX=repeatZ=6619035. In general, I have not seen this occur below
    // frequency*repeat = 65536 or so (e.g. frequency=0.125, repeat=524288).
    private static final boolean PREPRIME = false;
    
    private static final int PRIME_X = 1091;
    private static final int PRIME_Y = 30869;
    private static final int PRIME_Z = 1879;
    private static final int PREPRIME_X = PREPRIME ? PRIME_X : 1;
    private static final int PREPRIME_Z = PREPRIME ? PRIME_Z : 1;
    private static final int POSTPRIME_X = !PREPRIME ? PRIME_X : 1;
    private static final int POSTPRIME_Z = !PREPRIME ? PRIME_Z : 1;

    private double m00, m02, m20, m22, m1;
    private static final double M1 = 0.577350269189626;
    
    float i00, i02, i20, i22;
    private static final float I1 = 0.577350269189626f;
    
    long ii00, ii02, ii20, ii22;
    long iModX, iModZ;
    
    float rSquared;
    float approxNormalizer;
    
    private final float[] di;
    private final long[] dbp;
    
    /*
     * Constructor
     */
    
    public OpenSimplex2S_ImprovedXZPlanes_TileableXZ(double xzFrequency, double yFrequency, int repeatX, int repeatZ) {
        
        // Generate initial transform matrix in X and Z, by plugging (repeatX, 0) and (0, repeatZ), multiplied
        // by the frequency, into the 2D / triangular part of the base rotation matrix. Round the results to
        // integers, and divide each vector by repeatX and repeatZ respectively so that when they're used in the
        // adjusted matrix, (repeatX, 0) and (0, repeatZ) produce these integer results which are close to what
        // they would produce ordinarily. One could instead find the true closest points, but this works fine.
        // The resulting matrix includes the frequency.
        double SX0 = repeatX * -0.211324865405187;
        double S0Z = repeatZ * -0.211324865405187;
        long im00 = fastRoundL(xzFrequency * (repeatX + SX0)), im02 = fastRoundL(xzFrequency * S0Z);
        long im20 = fastRoundL(xzFrequency * SX0), im22 = fastRoundL(xzFrequency * (repeatZ + S0Z));
        m00 = im00 / (double)repeatX; m02 = im02 / (double)repeatZ;
        m20 = im20 / (double)repeatX; m22 = im22 / (double)repeatZ;
        m1 = M1 * yFrequency;
        
        // Integer inverse matrix. When we invert the coordinates of the lattice vertices, which are in the
        // adjusted rotated coordinate space, we want them to stay as integers to keep modulo+hashing accurate.
        // We can achieve this by inverting the general form of the adjusted rotation matrix with the common
        // denominator (matrixScale) multiplied out, as well as multiplying repeatX and repeatZ out of the it
        // from that matrix. Note that, since ii01 = -(ii00 + ii02), ii01 is not stored explicitly. Same for i21.
        // A -1 is shifted to matrixScale, so that the results tend towards positive instead of negative.
        // I must have made an error somewhere else for this to be necessary, but it works now.
        // https://www.wolframalpha.com/input/?i=inverse+%5B%5Ba%2Cq%2Cb%5D%2C%5B-a-c%2Cq%2C-b-d%5D%2C%5Bc%2Cq%2Cd%5D%5D
        ii00 = (im02 + 2L * im22) * repeatX;
        ii02 = -(im22 + 2L * im02) * repeatX;
        ii20 = -(im00 + 2L * im20) * repeatZ;
        ii22 = (im20 + 2L * im00) * repeatZ;
        long matrixScale = -3 * (im02 * im20 - im00 * im22);
        
        // True inverse matrix (with frequency factored out). The above, finally dividing by matrixScale, gives
        // us the inverse of m00..m22. Multiplying by xzFrequency counters the original multiplication by it.
        double inverseMatrixScale = xzFrequency / matrixScale;
        i00 = (float)(ii00 * inverseMatrixScale);
        i02 = (float)(ii02 * inverseMatrixScale);
        i20 = (float)(ii20 * inverseMatrixScale);
        i22 = (float)(ii22 * inverseMatrixScale);
        
        // Modulos for the X/Z coordinates returned by that matrix.
        iModX = matrixScale * repeatX;
        iModZ = matrixScale * repeatZ;
        
        // Certain repeat periods at certain frequencies create trailing zeroes, which decrease the quality of
        // the hash, and lead to unsightly patterns in the noise. Remove them where we can. An alternative to this
        // could be to always transform back to lattice space after the modulo, but this avoids that runtime cost.
        while (((ii00 | ii02 | iModX) & 1) == 0 && (ii00 | ii02 | iModX) != 0) {
            ii00 >>= 1; ii02 >>= 1; iModX >>= 1;
        }
        while (((ii20 | ii22 | iModZ) & 1) == 0 && (ii20 | ii22 | iModZ) != 0) {
            ii20 >>= 1; ii22 >>= 1; iModZ >>= 1;
        }
        
        // Apply hash primes to inverse matrix (if option set). We can modulo with them already in place.
        ii00 *= PREPRIME_X; ii02 *= PREPRIME_X;
        ii20 *= PREPRIME_Z; ii22 *= PREPRIME_Z;
        iModX *= PREPRIME_X; iModZ *= PREPRIME_Z;
        
        // Avoid discontinuities as the transform slightly changes the distance between lattice vertices.
        float rSquaredA = I1 * I1 * 2.25f;
        float rSquaredB = i00 * i00 + i20 * i20 + I1 * I1 * 0.25f;
        float rSquaredC = i02 * i02 + i22 * i22 + I1 * I1 * 0.25f;
        float rSquaredD = (i00 + i02) * (i00 + i02) + (i20 + i22) * (i20 + i22) + I1 * I1 * 0.25f;
        rSquared = Math.min(Math.min(rSquaredA, rSquaredB), Math.min(rSquaredC, rSquaredD));
        
        // Don't let that push us too far away from a range of -1 to 1 (it will not be exact either way).
        double approxNormalizerDouble = 0.75 / rSquared;
        approxNormalizerDouble *= approxNormalizerDouble;
        approxNormalizerDouble *= approxNormalizerDouble;
        approxNormalizer = (float)approxNormalizerDouble;
        
        // Lookup tables, customized for the repeat period.
        di = new float[8 * 16 * 4];
        dbp = new long[8 * 16 * 4];
        for (int i = 0; i < 8; i++) {
            int i1, j1, k1, i2, j2, k2;
            i1 = (i >> 0) & 1; j1 = (i >> 1) & 1; k1 = (i >> 2) & 1;
            i2 = i1 ^ 1; j2 = j1 ^ 1; k2 = k1 ^ 1;
            
            int index = i * 16;
            setLatticePointLookup(index++, i1, j1, k1, 0, 1);
            setLatticePointLookup(index++, i1 + i2, j1 + j2, k1 + k2, 1, 2);
            setLatticePointLookup(index++, i1 ^ 1, j1, k1, 0, 5);
            setLatticePointLookup(index++, i1, j1 ^ 1, k1 ^ 1, 0, 4);
            setLatticePointLookup(index++, i1 + (i2 ^ 1), j1 + j2, k1 + k2, 1, 6);
            setLatticePointLookup(index++, i1 + i2, j1 + (j2 ^ 1), k1 + (k2 ^ 1), 1, 6);
            setLatticePointLookup(index++, i1, j1 ^ 1, k1, 0, 9);
            setLatticePointLookup(index++, i1 ^ 1, j1, k1 ^ 1, 0, 8);
            setLatticePointLookup(index++, i1 + i2, j1 + (j2 ^ 1), k1 + k2, 1, 10);
            setLatticePointLookup(index++, i1 + (i2 ^ 1), j1 + j2, k1 + (k2 ^ 1), 1, 10);
            setLatticePointLookup(index++, i1, j1, k1 ^ 1, 0, 13);
            setLatticePointLookup(index++, i1 ^ 1, j1 ^ 1, k1, 0, 12);
            setLatticePointLookup(index++, i1 + i2, j1 + j2, k1 + (k2 ^ 1), 1, 14);
            setLatticePointLookup(index++, i1 + (i2 ^ 1), j1 + (j2 ^ 1), k1 + k2, 1, 14);
        }
    }
    
    private void setLatticePointLookup(int index, int xrv, int yrv, int zrv, int lattice, int nextOnSuccess) {
        index *= 4;
        
        float dxr = lattice * 0.5f - xrv;
        float dyr = lattice * 0.5f - yrv;
        float dzr = lattice * 0.5f - zrv;
        float rxyi = dxr - dyr;
        float rzyi = dzr - dyr;
        di[index + 0] = i00 * rxyi + i02 * rzyi;
        di[index + 1] = I1 * (dxr + dyr + dzr);
        di[index + 2] = i20 * rxyi + i22 * rzyi;
        
        int xyrv = xrv - yrv;
        int zyrv = zrv - yrv;
        dbp[index + 0] = (ii00 * xyrv + ii02 * zyrv);
        dbp[index + 1] = ((xrv + yrv + zrv) * 2 + lattice) * PRIME_Y;
        dbp[index + 2] = (ii20 * xyrv + ii22 * zyrv);
        dbp[index + 3] = nextOnSuccess * 4;
    }
    
    /*
     * Evaluators
     */
    
    public double noise2(long seed, double x, double z) {
        return noise3(seed, x, 0.8660254037844386, z);
    }
    
    public double noise3(long seed, double x, double y, double z) {
        
        // Approximate rotation, adjusted for tiling
        double yy = y * m1;
        double rx = m00 * x + m02 * z;
        double rz = m20 * x + m22 * z;
        double ry = yy - rx - rz;
        rx += yy; rz += yy;
        
        // Relative coordinates like normal
        int rxb = fastFloor(rx);
        int ryb = fastFloor(ry);
        int rzb = fastFloor(rz);
        float rxi = (float)(rx - rxb);
        float ryi = (float)(ry - ryb);
        float rzi = (float)(rz - rzb);
        
        // Undo entire approximate rotation on relative cell coordinates
        float rxyi = rxi - ryi;
        float rzyi = rzi - ryi;
        float xi = i00 * rxyi + i02 * rzyi;
        float yi = I1 * (rxi + ryi + rzi);
        float zi = i20 * rxyi + i22 * rzyi;
        
        // Undo entire approximate rotation on base vertex coordinates,
        // pre-multiplied by a constant so they stay as integers,
        // and pre-multiplied by primes to speed up hashing,
        // since the modulo can still be performed after prime multiplication.
        int rxyb = rxb - ryb;
        int rzyb = rzb - ryb;
        long xbp = (ii00 * rxyb + ii02 * rzyb);
        long ybp = (rxb + ryb + rzb) * (2 * PRIME_Y);
        long zbp = (ii20 * rxyb + ii22 * rzyb);
        
        // Vertex loop.
        int index = ((int)(rxi + 63.5f) & 64) | ((int)(ryi + 127.5f) & 128) | ((int)(rzi + 255.5f) & 256);
        int i = 0;
        float value = 0;
        do {
            float dx = xi + di[index | i | 0];
            float dy = yi + di[index | i | 1];
            float dz = zi + di[index | i | 2];
            float a = dx * dx + dy * dy + dz * dz;
            
            if (a < rSquared) {
                // Convert to actual spherical bump function
                a -= rSquared;
                a *= a; a *= a;
                
                // Actual primed modulo-able offsets for this vertex
                long xvp = xbp + dbp[index | i | 0];
                long yvp = ybp + dbp[index | i | 1];
                long zvp = zbp + dbp[index | i | 2];
                
                // Local modulo. Could be substituted for a true modulo if needed.
                // This is perfectly fine if you only ever evaluate the noise inside (or slightly outside) the repeat boundary.
                if (xvp < 0) xvp += iModX;
                else if (xvp >= iModX) xvp -= iModX;
                if (zvp < 0) zvp += iModZ;
                else if (zvp >= iModZ) zvp -= iModZ;
                
                // Hash
                long hash = seed ^ (xvp * POSTPRIME_X) ^ yvp ^ (zvp * POSTPRIME_Z);
                hash *= 2325943009213694033L;
                hash ^= hash >> 30;
                
                // Pseudo-modulo on index to effectively pick from 0-47
                // (repeated like 0012 3345 6678...)
                // This has nothing to do with noise tiling,
                // just clean gradient picking from a non-power-of-two sized set.
                long giL = hash & 0x3FF_FFFF_FFFF_FFFFL;
                giL *= 0x555_5555_5555_5555L;
                int gi = (int)(giL >> 56) & 0xFC;
                
                // Gradient, multiply, and add.
                value += a * (dx * GRAD3[gi + 0] + dy * GRAD3[gi + 1] + dz * GRAD3[gi + 2]);
                
                // Next on success
                i = (int)dbp[index | i | 3];
            } else i += 4;
        } while (i < 56);
        
        return value * approxNormalizer;
    }
    
    /*
     * Utility
     */
    
    private static int fastFloor(double x) {
        int xi = (int)x;
        return x < xi ? xi - 1 : xi;
    }
    
    private static long fastRoundL(double x) {
        return x > 0 ? (long)(x + 0.5) : (long)(x - 0.5);
    }
    
    /*
     * Gradients
     */
    
    public static final double N3 = 0.2781926117527186;
    static float[] GRAD3;
    static {
        double[] grad3 = {
            -2.22474487139,      -2.22474487139,      -1.0,                 0.0,
            -2.22474487139,      -2.22474487139,       1.0,                 0.0,
            -3.0862664687972017, -1.1721513422464978,  0.0,                 0.0,
            -1.1721513422464978, -3.0862664687972017,  0.0,                 0.0,
            -2.22474487139,      -1.0,                -2.22474487139,       0.0,
            -2.22474487139,       1.0,                -2.22474487139,       0.0,
            -1.1721513422464978,  0.0,                -3.0862664687972017,  0.0,
            -3.0862664687972017,  0.0,                -1.1721513422464978,  0.0,
            -2.22474487139,      -1.0,                 2.22474487139,       0.0,
            -2.22474487139,       1.0,                 2.22474487139,       0.0,
            -3.0862664687972017,  0.0,                 1.1721513422464978,  0.0,
            -1.1721513422464978,  0.0,                 3.0862664687972017,  0.0,
            -2.22474487139,       2.22474487139,      -1.0,                 0.0,
            -2.22474487139,       2.22474487139,       1.0,                 0.0,
            -1.1721513422464978,  3.0862664687972017,  0.0,                 0.0,
            -3.0862664687972017,  1.1721513422464978,  0.0,                 0.0,
            -1.0,                -2.22474487139,      -2.22474487139,       0.0,
             1.0,                -2.22474487139,      -2.22474487139,       0.0,
             0.0,                -3.0862664687972017, -1.1721513422464978,  0.0,
             0.0,                -1.1721513422464978, -3.0862664687972017,  0.0,
            -1.0,                -2.22474487139,       2.22474487139,       0.0,
             1.0,                -2.22474487139,       2.22474487139,       0.0,
             0.0,                -1.1721513422464978,  3.0862664687972017,  0.0,
             0.0,                -3.0862664687972017,  1.1721513422464978,  0.0,
            -1.0,                 2.22474487139,      -2.22474487139,       0.0,
             1.0,                 2.22474487139,      -2.22474487139,       0.0,
             0.0,                 1.1721513422464978, -3.0862664687972017,  0.0,
             0.0,                 3.0862664687972017, -1.1721513422464978,  0.0,
            -1.0,                 2.22474487139,       2.22474487139,       0.0,
             1.0,                 2.22474487139,       2.22474487139,       0.0,
             0.0,                 3.0862664687972017,  1.1721513422464978,  0.0,
             0.0,                 1.1721513422464978,  3.0862664687972017,  0.0,
             2.22474487139,      -2.22474487139,      -1.0,                 0.0,
             2.22474487139,      -2.22474487139,       1.0,                 0.0,
             1.1721513422464978, -3.0862664687972017,  0.0,                 0.0,
             3.0862664687972017, -1.1721513422464978,  0.0,                 0.0,
             2.22474487139,      -1.0,                -2.22474487139,       0.0,
             2.22474487139,       1.0,                -2.22474487139,       0.0,
             3.0862664687972017,  0.0,                -1.1721513422464978,  0.0,
             1.1721513422464978,  0.0,                -3.0862664687972017,  0.0,
             2.22474487139,      -1.0,                 2.22474487139,       0.0,
             2.22474487139,       1.0,                 2.22474487139,       0.0,
             1.1721513422464978,  0.0,                 3.0862664687972017,  0.0,
             3.0862664687972017,  0.0,                 1.1721513422464978,  0.0,
             2.22474487139,       2.22474487139,      -1.0,                 0.0,
             2.22474487139,       2.22474487139,       1.0,                 0.0,
             3.0862664687972017,  1.1721513422464978,  0.0,                 0.0,
             1.1721513422464978,  3.0862664687972017,  0.0,                 0.0
        };
        
        // Copy into GRAD, rotated, and with the first of every 3 repeated like 0012 3345 6678.
        GRAD3 = new float[grad3.length * 4 / 3];
        for (int i = 0, j = 0; i < grad3.length; i += 4*3, j += 4*4) {
            for (int k = 0; k < 16; k += 4) {
                int k2 = k == 0 ? 0 : k - 4;
                double gxr = grad3[i + k2 + 0];
                double gyr = grad3[i + k2 + 1];
                double gzr = grad3[i + k2 + 2];
                double s2 = (gxr + gzr) * -0.211324865405187f;
                double yy = gyr * 0.577350269189626f;
                GRAD3[j + k + 0] = (float)((gxr + s2 - yy) / N3);
                GRAD3[j + k + 1] = (float)((gxr + gzr + gyr) * 0.577350269189626 / N3);
                GRAD3[j + k + 2] = (float)((gzr + s2 - yy) / N3);
            }
        }
    }
}
