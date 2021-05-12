/**
 * Domain-rotated Perlin noise, speed-optimized for column generation.
 * Lock in an X/Z position, and call getForY(y) for each successive Y position.
 * It remembers enough of the current state of the noise to achieve a ~1.8x speedup.
 *
 * Incorporates "Improve XZ Planes" domain rotation to eliminate square bias from X/Z planes.
 * Y points up the main diagonal of the noise grid.
 *
 * Not recommended for 2D noise where X and Y are both horizontal, as X and Z have different.
 * characteristics compared to Y. I may later release noise better suited to this purpose.
 *
 * @author K.jpg
 */

public class RotoperlinCachedColumnNoise {
    private static final double ROOT3 = 1.7320508075688772;
    private static final double ROOT3OVER3 = 0.577350269189626;

    private static final int PERMUTATION_TABLE_SIZE_EXPONENT = 8;

    private static final int PSIZE = 1 << PERMUTATION_TABLE_SIZE_EXPONENT;
    private static final int PMASK = PSIZE - 1;
    
    private short[] perm;

    public RotoperlinCachedColumnNoise(long seed) {
        perm = new short[PSIZE];
        short[] source = new short[PSIZE];
        for (short i = 0; i < PSIZE; i++)
            source[i] = i;
        for (int i = PSIZE - 1; i >= 0; i--) {
            seed = seed * 6364136223846793005L + 1442695040888963407L;
            int r = (int)((seed + 31) % (i + 1));
            if (r < 0)
                r += (i + 1);
            perm[i] = source[r];
            source[r] = source[i];
        }
    }
    
    public ColumnGenerator columnGenerator() {
        return new ColumnGenerator();
    }
    
    public class ColumnGenerator {
    
        private double x, z;
        private double localMinY, localMaxY, y0;
        private double xri0, yri0, zri0;
        private double g000b, g001b, g010b, g011b, g100b, g101b, g110b, g111b;
        private double g000v, g001v, g010v, g011v, g100v, g101v, g110v, g111v;
        
        public ColumnGenerator() {
            localMinY = Double.POSITIVE_INFINITY;
        }
        
        public void setXZ(double x, double z) {
            this.x = x;
            this.z = z;
            localMinY = Double.POSITIVE_INFINITY;
        }
        
        public double getForY(double y) {
            if (y < localMinY || y > localMaxY) {
                
                // Domain rotation
                double xz = x + z;
                double s2 = xz * -0.211324865405187;
                double yy = y * ROOT3OVER3;
                double xr = x + (s2 + yy);
                double zr = z + (s2 + yy);
                double yr = xz * -ROOT3OVER3 + yy;
                
                // Cube base and bounds
                int xrb = fastFloor(xr), yrb = fastFloor(yr), zrb = fastFloor(zr);
                double xri = xr - xrb, yri = yr - yrb, zri = zr - zrb;
                
                // Set the current coordinates since things will be relative to them in the dependent step.
                xri0 = xri;
                yri0 = yri;
                zri0 = zri;
                
                // Set bounds of Y that will stay inside the cube.
                localMinY = y - ROOT3 * (xri < zri ? (xri < yri ? xri : yri) : (zri < yri ? zri : yri));
                localMaxY = y + ROOT3 * (1 - (xri >= zri ? (xri >= yri ? xri : yri) : (zri >= yri ? zri : yri)));
                y0 = y;
                
                // Gradient vector lookup indices
                int gi000 =  perm[perm[perm[xrb & PMASK] ^ (yrb & PMASK)] ^ (zrb & PMASK)] & 0xFC;
                int gi001 =  perm[perm[perm[xrb & PMASK] ^ (yrb & PMASK)] ^ ((zrb + 1) & PMASK)] & 0xFC;
                int gi010 =  perm[perm[perm[xrb & PMASK] ^ ((yrb + 1) & PMASK)] ^ (zrb & PMASK)] & 0xFC;
                int gi011 =  perm[perm[perm[xrb & PMASK] ^ ((yrb + 1) & PMASK)] ^ ((zrb + 1) & PMASK)] & 0xFC;
                int gi100 =  perm[perm[perm[(xrb + 1) & PMASK] ^ (yrb & PMASK)] ^ (zrb & PMASK)] & 0xFC;
                int gi101 =  perm[perm[perm[(xrb + 1) & PMASK] ^ (yrb & PMASK)] ^ ((zrb + 1) & PMASK)] & 0xFC;
                int gi110 =  perm[perm[perm[(xrb + 1) & PMASK] ^ ((yrb + 1) & PMASK)] ^ (zrb & PMASK)] & 0xFC;
                int gi111 =  perm[perm[perm[(xrb + 1) & PMASK] ^ ((yrb + 1) & PMASK)] ^ ((zrb + 1) & PMASK)] & 0xFC;
                
                // Gradients vectors: base value at this Y coordinate, as well as slope in this Y coordinate (rotated)
                g000b = GRADIENTS_3D[gi000 | 0] * xri + GRADIENTS_3D[gi000 | 1] * yri + GRADIENTS_3D[gi000 | 2] * zri; g000v = GRADIENTS_3D[gi000 | 3];
                g100b = GRADIENTS_3D[gi100 | 0] * (xri - 1) + GRADIENTS_3D[gi100 | 1] * yri + GRADIENTS_3D[gi100 | 2] * zri; g100v = GRADIENTS_3D[gi100 | 3];
                g001b = GRADIENTS_3D[gi001 | 0] * xri + GRADIENTS_3D[gi001 | 1] * yri + GRADIENTS_3D[gi001 | 2] * (zri - 1); g001v = GRADIENTS_3D[gi001 | 3];
                g101b = GRADIENTS_3D[gi101 | 0] * (xri - 1) + GRADIENTS_3D[gi101 | 1] * yri + GRADIENTS_3D[gi101 | 2] * (zri - 1); g101v = GRADIENTS_3D[gi101 | 3];
                g010b = GRADIENTS_3D[gi010 | 0] * xri + GRADIENTS_3D[gi010 | 1] * (yri - 1) + GRADIENTS_3D[gi010 | 2] * zri; g010v = GRADIENTS_3D[gi010 | 3];
                g110b = GRADIENTS_3D[gi110 | 0] * (xri - 1) + GRADIENTS_3D[gi110 | 1] * (yri - 1) + GRADIENTS_3D[gi110 | 2] * zri; g110v = GRADIENTS_3D[gi110 | 3];
                g011b = GRADIENTS_3D[gi011 | 0] * xri + GRADIENTS_3D[gi011 | 1] * (yri - 1) + GRADIENTS_3D[gi011 | 2] * (zri - 1); g011v = GRADIENTS_3D[gi011 | 3];
                g111b = GRADIENTS_3D[gi111 | 0] * (xri - 1) + GRADIENTS_3D[gi111 | 1] * (yri - 1) + GRADIENTS_3D[gi111 | 2] * (zri - 1); g111v = GRADIENTS_3D[gi111 | 3];
                
            }
            
            double t = y - y0;
            double g000 = g000b + g000v * t;
            double g100 = g100b + g100v * t;
            double g001 = g001b + g001v * t;
            double g101 = g101b + g101v * t;
            double g010 = g010b + g010v * t;
            double g110 = g110b + g110v * t;
            double g011 = g011b + g011v * t;
            double g111 = g111b + g111v * t;
            double fadeX = fadeCurve(xri0 + t * ROOT3OVER3);
            double fadeY = fadeCurve(yri0 + t * ROOT3OVER3);
            double fadeZ = fadeCurve(zri0 + t * ROOT3OVER3);
            double g00Z = fadeZ * (g001 - g000) + g000;
            double g01Z = fadeZ * (g011 - g010) + g010;
            double g10Z = fadeZ * (g101 - g100) + g100;
            double g11Z = fadeZ * (g111 - g110) + g110;
            double g0YZ = fadeY * (g01Z - g00Z) + g00Z;
            double g1YZ = fadeY * (g11Z - g10Z) + g10Z;
            double gXYZ = fadeX * (g1YZ - g0YZ) + g0YZ;
            
            return gXYZ;
        }
    }

    private static double fadeCurve(double t) {
        return t * t * t * (10 + t * (-15 + t * 6));
    }

    private static int fastFloor(double x) {
        int xi = (int)x;
        return x < xi ? xi - 1 : xi;
    }

    private static final double[] GRADIENTS_3D = {
        0, 1, 1, 0,  0,-1, 1, 0,  0, 1,-1, 0,  0,-1,-1, 0,
        1, 0, 1, 0, -1, 0, 1, 0,  1, 0,-1, 0, -1, 0,-1, 0,
        1, 1, 0, 0, -1, 1, 0, 0,  1,-1, 0, 0, -1,-1, 0, 0,
        0, 1, 1, 0,  0,-1, 1, 0,  0, 1,-1, 0,  0,-1,-1, 0,
        1, 0, 1, 0, -1, 0, 1, 0,  1, 0,-1, 0, -1, 0,-1, 0,
        1, 1, 0, 0, -1, 1, 0, 0,  1,-1, 0, 0, -1,-1, 0, 0,
        0, 1, 1, 0,  0,-1, 1, 0,  0, 1,-1, 0,  0,-1,-1, 0,
        1, 0, 1, 0, -1, 0, 1, 0,  1, 0,-1, 0, -1, 0,-1, 0,
        1, 1, 0, 0, -1, 1, 0, 0,  1,-1, 0, 0, -1,-1, 0, 0,
        0, 1, 1, 0,  0,-1, 1, 0,  0, 1,-1, 0,  0,-1,-1, 0,
        1, 0, 1, 0, -1, 0, 1, 0,  1, 0,-1, 0, -1, 0,-1, 0,
        1, 1, 0, 0, -1, 1, 0, 0,  1,-1, 0, 0, -1,-1, 0, 0,
        0, 1, 1, 0,  0,-1, 1, 0,  0, 1,-1, 0,  0,-1,-1, 0,
        1, 0, 1, 0, -1, 0, 1, 0,  1, 0,-1, 0, -1, 0,-1, 0,
        1, 1, 0, 0, -1, 1, 0, 0,  1,-1, 0, 0, -1,-1, 0, 0,
        1, 1, 0, 0,  0,-1, 1, 0, -1, 1, 0, 0,  0,-1,-1, 0
    };

    public static final double NORMALIZING_MULTIPLIER_3D = 1.0 / 1.0363538112118038;
    static {
        for (int i = 0; i < GRADIENTS_3D.length; i += 4) {
            GRADIENTS_3D[i] *= NORMALIZING_MULTIPLIER_3D;
        }
        for (int i = 0; i < GRADIENTS_3D.length; i += 4) {
            double inRotatedVerticalDirection = 0;
            for (int j = 0; j < 3; j++) {
                inRotatedVerticalDirection += GRADIENTS_3D[i + j];
                GRADIENTS_3D[i | j] *= NORMALIZING_MULTIPLIER_3D;
            }
            GRADIENTS_3D[i | 3] = inRotatedVerticalDirection * (ROOT3OVER3 * NORMALIZING_MULTIPLIER_3D);
        }
    }
    
}
