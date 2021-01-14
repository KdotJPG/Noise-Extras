public class DomainRotatedClassicNoise {

    private static final int PSIZE = 2048;
    private static final int PMASK = 2047;

    private short[] perm;

    public DomainRotatedClassicNoise(long seed) {
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

    public DomainRotatedClassicNoise(java.util.Random rand) {
        perm = new short[PSIZE];
        short[] source = new short[PSIZE];
        for (short i = 0; i < PSIZE; i++)
            source[i] = i;
        for (int i = PSIZE - 1; i >= 0; i--) {
            int r = rand.nextInt(i + 1);
            if (r < 0)
                r += (i + 1);
            perm[i] = source[r];
            source[r] = source[i];
        }
    }
    
    public double noise2(double x, double y) {
        // No way to hide grid bias in Perlin without using the 3D noise.
        // The value chosen for Y is far away from the origin, and between plane slices.
        return noise3(x, (PSIZE + 1) * 0.8660254037844386, y);
    }

    public double noise3(double x, double y, double z) {

        // Re-orient the cubic lattice without skewing, to make Y look down <1,1,1>.
        // This hides the vast majority of the square alignment characteristic of Perlin, in X/Z planes.
        // It also reduces the variation in characteristic of the different horizontal slices.
        double xz = x + z;
        double s2 = xz * -0.211324865405187;
        double yy = y * 0.577350269189626;
        double xr = x + s2 + yy; double zr = z + s2 + yy;
        double yr = xz * -0.577350269189626 + yy;

        // The rest is a modified Perlin.
        int xrb = fastFloor(xr), yrb = fastFloor(yr), zrb = fastFloor(zr);
        double xri = xr - xrb, yri = yr - yrb, zri = zr - zrb;
        double g000 = grad3(4 * perm[perm[perm[xrb & PMASK] ^ (yrb & PMASK)] ^ (zrb & PMASK)], xri, yri, zri);
        double g001 = grad3(4 * perm[perm[perm[xrb & PMASK] ^ (yrb & PMASK)] ^ ((zrb + 1) & PMASK)], xri, yri, zri - 1);
        double g010 = grad3(4 * perm[perm[perm[xrb & PMASK] ^ ((yrb + 1) & PMASK)] ^ (zrb & PMASK)], xri, yri - 1, zri);
        double g011 = grad3(4 * perm[perm[perm[xrb & PMASK] ^ ((yrb + 1) & PMASK)] ^ ((zrb + 1) & PMASK)], xri, yri - 1, zri - 1);
        double g100 = grad3(4 * perm[perm[perm[(xrb + 1) & PMASK] ^ (yrb & PMASK)] ^ (zrb & PMASK)], xri - 1, yri, zri);
        double g101 = grad3(4 * perm[perm[perm[(xrb + 1) & PMASK] ^ (yrb & PMASK)] ^ ((zrb + 1) & PMASK)], xri - 1, yri, zri - 1);
        double g110 = grad3(4 * perm[perm[perm[(xrb + 1) & PMASK] ^ ((yrb + 1) & PMASK)] ^ (zrb & PMASK)], xri - 1, yri - 1, zri);
        double g111 = grad3(4 * perm[perm[perm[(xrb + 1) & PMASK] ^ ((yrb + 1) & PMASK)] ^ ((zrb + 1) & PMASK)], xri - 1, yri - 1, zri - 1);
        double fadeX = fadeCurve(xri);
        double fadeY = fadeCurve(yri);
        double fadeZ = fadeCurve(zri);
        double g00Z = (1 - fadeZ) * g000 + fadeZ * g001;
        double g01Z = (1 - fadeZ) * g010 + fadeZ * g011;
        double g10Z = (1 - fadeZ) * g100 + fadeZ * g101;
        double g11Z = (1 - fadeZ) * g110 + fadeZ * g111;
        double g0YZ = (1 - fadeY) * g00Z + fadeY * g01Z;
        double g1YZ = (1 - fadeY) * g10Z + fadeY * g11Z;
        double gXYZ = (1 - fadeX) * g0YZ + fadeX * g1YZ;

        return gXYZ;
    }

    public double noise3WithDerivatives(Vector3 derivatives, double x, double y, double z) {

        // Re-orient the cubic lattice without skewing, to make Y look down <1,1,1>.
        // This hides the vast majority of the square alignment characteristic of Perlin, in X/Z planes.
        // It also reduces the variation in characteristic of the different horizontal slices.
        double xz = x + z;
        double s2 = xz * -0.211324865405187;
        double yy = y * 0.577350269189626;
        double xr = x + s2 + yy; double zr = z + s2 + yy;
        double yr = xz * -0.577350269189626 + yy;

        // The rest is a modified Perlin.
        int xrb = fastFloor(xr), yrb = fastFloor(yr), zrb = fastFloor(zr);
        double xri = xr - xrb, yri = yr - yrb, zri = zr - zrb;
        
        int gi000 = 4 * perm[perm[perm[xrb & PMASK] ^ (yrb & PMASK)] ^ (zrb & PMASK)];
        int gi001 = 4 * perm[perm[perm[xrb & PMASK] ^ (yrb & PMASK)] ^ ((zrb + 1) & PMASK)];
        int gi010 = 4 * perm[perm[perm[xrb & PMASK] ^ ((yrb + 1) & PMASK)] ^ (zrb & PMASK)];
        int gi011 = 4 * perm[perm[perm[xrb & PMASK] ^ ((yrb + 1) & PMASK)] ^ ((zrb + 1) & PMASK)];
        int gi100 = 4 * perm[perm[perm[(xrb + 1) & PMASK] ^ (yrb & PMASK)] ^ (zrb & PMASK)];
        int gi101 = 4 * perm[perm[perm[(xrb + 1) & PMASK] ^ (yrb & PMASK)] ^ ((zrb + 1) & PMASK)];
        int gi110 = 4 * perm[perm[perm[(xrb + 1) & PMASK] ^ ((yrb + 1) & PMASK)] ^ (zrb & PMASK)];
        int gi111 = 4 * perm[perm[perm[(xrb + 1) & PMASK] ^ ((yrb + 1) & PMASK)] ^ ((zrb + 1) & PMASK)];
        double g000 = grad3(gi000, xri, yri, zri);
        double g001 = grad3(gi001, xri, yri, zri - 1);
        double g010 = grad3(gi010, xri, yri - 1, zri);
        double g011 = grad3(gi011, xri, yri - 1, zri - 1);
        double g100 = grad3(gi100, xri - 1, yri, zri);
        double g101 = grad3(gi101, xri - 1, yri, zri - 1);
        double g110 = grad3(gi110, xri - 1, yri - 1, zri);
        double g111 = grad3(gi111, xri - 1, yri - 1, zri - 1);
        double fadeX = fadeCurve(xri);
        double fadeY = fadeCurve(yri);
        double fadeZ = fadeCurve(zri);
        double g00Z = (1 - fadeZ) * g000 + fadeZ * g001;
        double g01Z = (1 - fadeZ) * g010 + fadeZ * g011;
        double g10Z = (1 - fadeZ) * g100 + fadeZ * g101;
        double g11Z = (1 - fadeZ) * g110 + fadeZ * g111;
        double g0YZ = (1 - fadeY) * g00Z + fadeY * g01Z;
        double g1YZ = (1 - fadeY) * g10Z + fadeY * g11Z;
        double gXYZ = (1 - fadeX) * g0YZ + fadeX * g1YZ;
        
        double dFadeX = dFadeCurve(xri);
        double dFadeY = dFadeCurve(yri);
        double dFadeZ = dFadeCurve(zri);
        double dXg00Z = (1 - fadeZ) * GRADIENTS_3D[gi000 + 0] + fadeZ * GRADIENTS_3D[gi001 + 0];
        double dYg00Z = (1 - fadeZ) * GRADIENTS_3D[gi000 + 1] + fadeZ * GRADIENTS_3D[gi001 + 1];
        double dZg00Z = (1 - fadeZ) * GRADIENTS_3D[gi000 + 2] + fadeZ * GRADIENTS_3D[gi001 + 2] + (-dFadeZ * g000 + dFadeZ * g001);
        double dXg01Z = (1 - fadeZ) * GRADIENTS_3D[gi010 + 0] + fadeZ * GRADIENTS_3D[gi011 + 0];
        double dYg01Z = (1 - fadeZ) * GRADIENTS_3D[gi010 + 1] + fadeZ * GRADIENTS_3D[gi011 + 1];
        double dZg01Z = (1 - fadeZ) * GRADIENTS_3D[gi010 + 2] + fadeZ * GRADIENTS_3D[gi011 + 2] + (-dFadeZ * g010 + dFadeZ * g011);
        double dXg10Z = (1 - fadeZ) * GRADIENTS_3D[gi100 + 0] + fadeZ * GRADIENTS_3D[gi101 + 0];
        double dYg10Z = (1 - fadeZ) * GRADIENTS_3D[gi100 + 1] + fadeZ * GRADIENTS_3D[gi101 + 1];
        double dZg10Z = (1 - fadeZ) * GRADIENTS_3D[gi100 + 2] + fadeZ * GRADIENTS_3D[gi101 + 2] + (-dFadeZ * g100 + dFadeZ * g101);
        double dXg11Z = (1 - fadeZ) * GRADIENTS_3D[gi110 + 0] + fadeZ * GRADIENTS_3D[gi111 + 0];
        double dYg11Z = (1 - fadeZ) * GRADIENTS_3D[gi110 + 1] + fadeZ * GRADIENTS_3D[gi111 + 1];
        double dZg11Z = (1 - fadeZ) * GRADIENTS_3D[gi110 + 2] + fadeZ * GRADIENTS_3D[gi111 + 2] + (-dFadeZ * g110 + dFadeZ * g111);
        double dXg0YZ = (1 - fadeY) * dXg00Z + fadeY * dXg01Z;
        double dYg0YZ = (1 - fadeY) * dYg00Z + fadeY * dYg01Z + (-dFadeY * g00Z + dFadeY * g01Z);
        double dZg0YZ = (1 - fadeY) * dZg00Z + fadeY * dZg01Z;
        double dXg1YZ = (1 - fadeY) * dXg10Z + fadeY * dXg11Z;
        double dYg1YZ = (1 - fadeY) * dYg10Z + fadeY * dYg11Z + (-dFadeY * g10Z + dFadeY * g11Z);
        double dZg1YZ = (1 - fadeY) * dZg10Z + fadeY * dZg11Z;
        double dXgXYZ = (1 - fadeX) * dXg0YZ + fadeX * dXg1YZ + (-dFadeX * g0YZ + dFadeX * g1YZ);
        double dYgXYZ = (1 - fadeX) * dYg0YZ + fadeX * dYg1YZ;
        double dZgXYZ = (1 - fadeX) * dZg0YZ + fadeX * dZg1YZ;
        double xzd = dXgXYZ + dZgXYZ;
        double s2d = xzd * -0.211324865405187;
        double yyd = dYgXYZ * 0.577350269189626;
        derivatives.x = dXgXYZ + s2d - yyd;
        derivatives.z = dZgXYZ + s2d - yyd;
        derivatives.y = xzd * 0.577350269189626 + yyd;
        
        return gXYZ;
    }

	/*
	 * Utility
	 */

    private static int fastFloor(double x) {
        int xi = (int)x;
        return x < xi ? xi - 1 : xi;
    }

    private static double grad3(int gi, double dx, double dy, double dz) {
        return GRADIENTS_3D[gi + 0] * dx + GRADIENTS_3D[gi + 1] * dy + GRADIENTS_3D[gi + 2] * dz;
    }

    private static double fadeCurve(double t) {
        return t * t * t * (10 + t * (-15 + t * 6));
    }

    private static double dFadeCurve(double t) {
        return t * t * (30 + t * (-60 + t * 30));
    }

	/*
	 * Gradients
	 */

    public static final double NORMALIZING_MULTIPLIER = 0.964921414852142333984375;
    private static final double[] GRADIENTS_3D;
    static {
        double[] grad3 = {  
             1,  1,  0,  0,
            -1,  1,  0,  0,
             1, -1,  0,  0,
            -1, -1,  0,  0,
             1,  0,  1,  0,
            -1,  0,  1,  0,
             1,  0, -1,  0,
            -1,  0, -1,  0,
             0,  1,  1,  0,
             0, -1,  1,  0,
             0,  1, -1,  0,
             0, -1, -1,  0
        };
        GRADIENTS_3D = new double[PSIZE * 4];
        for (int i = 0; i < grad3.length; i++) {
            grad3[i] *= NORMALIZING_MULTIPLIER;
        }
        for (int i = 0; i < PSIZE; i++) {
            int iSrc = i % (grad3.length / 4);
            for (int j = 0; j < 3; j++) {
                GRADIENTS_3D[i * 4 + j] = grad3[iSrc * 4 + j];
            }
            double xz = GRADIENTS_3D[i * 4 + 0] + GRADIENTS_3D[i * 4 + 2];
            double s2 = xz * -0.211324865405187;
            double yy = GRADIENTS_3D[i * 4 + 1] * 0.577350269189626;
        }
    }

	/*
	 * Inner Classes
	 */
    
    public static class Vector3 {
        public double x, y, z;
    }
}
