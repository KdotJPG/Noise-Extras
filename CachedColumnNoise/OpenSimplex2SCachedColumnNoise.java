/**
 * OpenSimplex2S noise, speed-optimized for column generation.
 * Lock in an X/Z position, and call getForY(y) for each successive Y position.
 * It remembers enough of the current state of the noise to achieve a ~1.8x speedup.
 *
 * Uses "Improve XZ Planes" domain rotation to reduce subtle diagonal patterns.
 *
 * Not recommended for 2D noise where X and Y are both horizontal, as there is some bias.
 * I may later release noise better suited to this purpose.
 *
 * @author K.jpg
 */

public class OpenSimplex2SCachedColumnNoise {

    private static final double ROOT3 = 1.7320508075688772;
    private static final double ROOT3OVER3 = 0.577350269189626;

    private static final int PERMUTATION_TABLE_SIZE_EXPONENT = 8;

    private static final int PSIZE = 1 << PERMUTATION_TABLE_SIZE_EXPONENT;
    private static final int PMASK = PSIZE - 1;
    
    private short[] perm;

    public OpenSimplex2SCachedColumnNoise(long seed) {
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
    
        HalfColumnGenerator cA, cB;
        
        public ColumnGenerator() {
            cA = new HalfColumnGenerator();
            cB = new HalfColumnGenerator();
        }
        
        public void setXZ(double x, double z) {
            cA.setXZ(x, z);
            cB.setXZ(x, z);
        }
        
        public double getForY(double y) {
            return cA.getForY(y) + cB.getForY(y + ((PSIZE + 1) * (ROOT3/2)));
        }
    
    }
    
    private class HalfColumnGenerator {
    
        private double xs, zs, xzr;
        private double localMinY, localMaxY, y0;
        private int xrb, yrb, zrb;
        private double xi0, zi0;
        private double g000b, g100b, g010b, g001b;
        private double g000v, g100v, g010v, g001v;
        private double xzFalloff000, xzFalloff100, xzFalloff010, xzFalloff001;
        private double yb000, yb100, yb010, yb001;
        
        private double yBound000, yBound100, yBound010, yBound001;
        private double yMin000, yMin100, yMin010, yMin001;
        private double yMax000, yMax100, yMax010, yMax001;
        
        protected HalfColumnGenerator() {
            localMinY = Double.POSITIVE_INFINITY;
        }
        
        public void setXZ(double x, double z) {
            
            // Domain rotation, start
            double xz = x + z;
            double s2 = xz * -0.211324865405187;
            this.xs = x + s2;
            this.zs = z + s2;
            this.xzr = xz * -ROOT3OVER3;
            
            localMinY = Double.POSITIVE_INFINITY;
        }
        
        public double getForY(double y) {
            
            if (y < localMinY || y > localMaxY) {
                
                // Domain rotation, finish
                double yy = y * ROOT3OVER3;
                double xr = xs + yy;
                double zr = zs + yy;
                double yr = xzr + yy;
                
                // Cube base and bounds
                xrb = fastFloor(xr); yrb = fastFloor(yr); zrb = fastFloor(zr);
                double xri = xr - xrb, yri = yr - yrb, zri = zr - zrb;
                
                // Set bounds of Y that will stay inside the cube.
                localMinY = y - ROOT3 * (xri < zri ? (xri < yri ? xri : yri) : (zri < yri ? zri : yri));
                localMaxY = y + ROOT3 * (1 - (xri >= zri ? (xri >= yri ? xri : yri) : (zri >= yri ? zri : yri)));
                y0 = (yrb + xrb + zrb) * ROOT3OVER3;
                
                // Remove the (world-space) Y coordinate so we can consider everything relative to the base cube vertex.
                double flatten = (xri + yri + zri) * -0.3333333333333333;
                double xri0 = xri + flatten;
                double yri0 = yri + flatten;
                double zri0 = zri + flatten;
                
                // Planar dividers between contributing vertices, intersections with the current column.
                yBound000 = 1.5 * ROOT3OVER3;
                yBound100 = (1.5 * ROOT3OVER3) + (2 * ROOT3) * xri0;
                yBound010 = (1.5 * ROOT3OVER3) + (2 * ROOT3) * yri0;
                yBound001 = (1.5 * ROOT3OVER3) + (2 * ROOT3) * zri0;
                
                double s2i = (xri0 + zri0) * -0.211324865405187 - yri0 * ROOT3OVER3;
                xi0 = xri0 + s2i;
                zi0 = zri0 + s2i;
                xzFalloff000 = 0.75 - xi0 * xi0 - zi0 * zi0;
                
                // Trigger the compares to happen regardless
                yMin000 = yMin100 = yMin010 = yMin001 = Double.POSITIVE_INFINITY;
            }
            double t = y - y0;
            
            if (t < yMin000 || t > yMax000) {
                if (t > yBound000) { // 111
                    yb000 = ROOT3;
                    yMin000 = yBound000;
                    yMax000 = Double.POSITIVE_INFINITY;
                    int gi = gradIndex(xrb + 1, yrb + 1, zrb + 1);
                    g000b = RGRADIENTS_3D[gi | 0] * xi0 + RGRADIENTS_3D[gi | 1] * zi0;
                    g000v = RGRADIENTS_3D[gi | 2];
                } else { // 000
                    yb000 = 0;
                    yMin000 = Double.NEGATIVE_INFINITY;
                    yMax000 = yBound000;
                    int gi = gradIndex(xrb, yrb, zrb);
                    g000b = RGRADIENTS_3D[gi | 0] * xi0 + RGRADIENTS_3D[gi | 1] * zi0;
                    g000v = RGRADIENTS_3D[gi | 2];
                }
            }
            
            if (t < yMin100 || t > yMax100) {
                if (t > yBound100) { // 011
                    yb100 = 2 * ROOT3OVER3;
                    yMin100 = yBound100;
                    yMax100 = Double.POSITIVE_INFINITY;
                    double xi = xi0 + 0.788675134594813, zi = zi0 - 0.211324865405187;
                    int gi = gradIndex(xrb, yrb + 1, zrb + 1);
                    g100b = RGRADIENTS_3D[gi | 0] * xi + RGRADIENTS_3D[gi | 1] * zi;
                    g100v = RGRADIENTS_3D[gi | 2];
                    xzFalloff100 = 0.75 - xi * xi - zi * zi;
                } else { // 100
                    yb100 = ROOT3OVER3;
                    yMin100 = Double.NEGATIVE_INFINITY;
                    yMax100 = yBound100;
                    double xi = xi0 - 0.788675134594813, zi = zi0 + 0.211324865405187;
                    int gi = gradIndex(xrb + 1, yrb, zrb);
                    g100b = RGRADIENTS_3D[gi | 0] * xi + RGRADIENTS_3D[gi | 1] * zi;
                    g100v = RGRADIENTS_3D[gi | 2];
                    xzFalloff100 = 0.75 - xi * xi - zi * zi;
                }
            }
            
            if (t < yMin010 || t > yMax010) {
                if (t > yBound010) { // 101
                    yb010 = 2 * ROOT3OVER3;
                    yMin010 = yBound010;
                    yMax010 = Double.POSITIVE_INFINITY;
                    double xi = xi0 - 0.57735026918962, zi = zi0 - 0.57735026918962;
                    int gi = gradIndex(xrb + 1, yrb, zrb + 1);
                    g010b = RGRADIENTS_3D[gi | 0] * xi + RGRADIENTS_3D[gi | 1] * zi;
                    g010v = RGRADIENTS_3D[gi | 2];
                    xzFalloff010 = 0.75 - xi * xi - zi * zi;
                } else { // 010
                    yb010 = ROOT3OVER3;
                    yMin010 = Double.NEGATIVE_INFINITY;
                    yMax010 = yBound010;
                    double xi = xi0 + 0.577350269189626, zi = zi0 + 0.577350269189626;
                    int gi = gradIndex(xrb, yrb + 1, zrb);
                    g010b = RGRADIENTS_3D[gi | 0] * xi + RGRADIENTS_3D[gi | 1] * zi;
                    g010v = RGRADIENTS_3D[gi | 2];
                    xzFalloff010 = 0.75 - xi * xi - zi * zi;
                }
            }
            
            if (t < yMin001 || t > yMax001) {
                if (t > yBound001) { // 110
                    yb001 = 2 * ROOT3OVER3;
                    yMin001 = yBound001;
                    yMax001 = Double.POSITIVE_INFINITY;
                    double xi = xi0 - 0.211324865405187, zi = zi0 + 0.788675134594813;
                    int gi = gradIndex(xrb + 1, yrb + 1, zrb);
                    g001b = RGRADIENTS_3D[gi | 0] * xi + RGRADIENTS_3D[gi | 1] * zi;
                    g001v = RGRADIENTS_3D[gi | 2];
                    xzFalloff001 = 0.75 - xi * xi - zi * zi;
                } else { // 001
                    yb001 = ROOT3OVER3;
                    yMin001 = Double.NEGATIVE_INFINITY;
                    yMax001 = yBound001;
                    double xi = xi0 + 0.211324865405187, zi = zi0 - 0.788675134594813;
                    int gi = gradIndex(xrb, yrb, zrb + 1);
                    g001b = RGRADIENTS_3D[gi | 0] * xi + RGRADIENTS_3D[gi | 1] * zi;
                    g001v = RGRADIENTS_3D[gi | 2];
                    xzFalloff001 = 0.75 - xi * xi - zi * zi;
                }
            }
            
            double value = 0;
            
            double dy000 = t - yb000;
            double dy000Sq = dy000 * dy000;
            if (xzFalloff000 > dy000Sq) {
                double falloff = xzFalloff000 - dy000Sq;
                falloff *= falloff;
                value = falloff * falloff * (g000b + g000v * dy000);
            }
            
            double dy100 = t - yb100;
            double dy100Sq = dy100 * dy100;
            if (xzFalloff100 > dy100Sq) {
                double falloff = xzFalloff100 - dy100Sq;
                falloff *= falloff;
                value += falloff * falloff * (g100b + g100v * dy100);
            }
            
            double dy010 = t - yb010;
            double dy010Sq = dy010 * dy010;
            if (xzFalloff010 > dy010Sq) {
                double falloff = xzFalloff010 - dy010Sq;
                falloff *= falloff;
                value += falloff * falloff * (g010b + g010v * dy010);
            }
            
            double dy001 = t - yb001;
            double dy001Sq = dy001 * dy001;
            if (xzFalloff001 > dy001Sq) {
                double falloff = xzFalloff001 - dy001Sq;
                falloff *= falloff;
                value += falloff * falloff * (g001b + g001v * dy001);
            }
            
            return value;
        }
        
    }

    private static int fastFloor(double x) {
        int xi = (int)x;
        return x < xi ? xi - 1 : xi;
    }
    
    private int gradIndex(int xrv, int yrv, int zrv) {
        return perm[perm[perm[xrv & PMASK] ^ (yrv & PMASK)] ^ (zrv & PMASK)] & 0xFC;
    }

    private static final double[] RGRADIENTS_3D = {
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

    private static final double NORMALIZING_MULTIPLIER_3D = 9.046026385208288;
    static {
        for (int i = 0; i < RGRADIENTS_3D.length; i += 4) {
            double gx = RGRADIENTS_3D[i | 0];
            double gy = RGRADIENTS_3D[i | 1];
            double gz = RGRADIENTS_3D[i | 2];
            
            // Rotate and apply noise normalizing constant
            // Store as XZY so X and Z are next to each other when we read them at the same time.
			double s2 = (gx + gz) * -0.211324865405187 + gy * 0.577350269189626;
            RGRADIENTS_3D[i | 0] = (gx + s2) * NORMALIZING_MULTIPLIER_3D;
            RGRADIENTS_3D[i | 1] = (gz + s2) * NORMALIZING_MULTIPLIER_3D;
            RGRADIENTS_3D[i | 2] = (gy - gx - gz) * (ROOT3OVER3 * NORMALIZING_MULTIPLIER_3D);
        }
    }
    
}
