/**
 * OpenSimplex2S noise with
 * - orthonormal-frame-based vector-valued gradient noise endpoints, and
 * - high granularity gradient sampling
 */

public class OpenSimplex2S_WithVectorValuedSamplers_Fancy {

    private static final long PRIME_X = 0x5205402B9270C86FL;
    private static final long PRIME_Y = 0x598CD327003817B5L;
    private static final long PRIME_Z = 0x5BCC226E9FA0BACBL;
    private static final long HASH_MULTIPLIER = 0x53A3F72DEEC546F5L;
    private static final long SEED_FLIP_3D = -0x52D547B2E96ED629L;

    private static final double ROOT2OVER2 = 0.7071067811865476;
    private static final double SKEW_2D = 0.366025403784439;
    private static final double UNSKEW_2D = -0.21132486540518713;

    private static final double ROOT3OVER3 = 0.577350269189626;
    private static final double FALLBACK_ROTATE3 = 2.0 / 3.0;
    private static final double ROTATE3_ORTHOGONALIZER = UNSKEW_2D;

    private static final int N_SINES_360_EXP = 16;
    private static final int N_ROTOR_SLIDES_3D_90_EXP = N_SINES_360_EXP - 2;
    private static final int N_VECTOR_SCALES_3D_180_EXP = N_SINES_360_EXP - 1;
    private static final int N_SINES_360 = 1 << N_SINES_360_EXP;
    private static final int COSINE_OFFSET = N_SINES_360 / 4;
    private static final int N_ROTOR_SLIDES_3D_90 = 1 << N_ROTOR_SLIDES_3D_90_EXP;
    private static final int N_VECTOR_SCALES_3D_180 = 1 << N_VECTOR_SCALES_3D_180_EXP;

    private static final float NORMALIZER_2D = (float)(1.0 / 0.05518041029835581);
    private static final float NORMALIZER_3D = (float)(1.0 / 0.08677623918785461);

    private static final float RSQUARED_2D = 2.0f / 3.0f;
    private static final float RSQUARED_3D = 3.0f / 4.0f;

    /*
     * Noise Evaluators
     */

    /**
     * 2D OpenSimplex2S/SuperSimplex noise, standard lattice orientation.
     * Value range: [-1, 1]
     */
    public static float noise2(long seed, double x, double y) {

        // Get points for A2* lattice
        double s = SKEW_2D * (x + y);
        double xs = x + s, ys = y + s;

        return noise2_UnskewedBase(seed, xs, ys);
    }

    /**
     * 2D OpenSimplex2S/SuperSimplex noise, with Y pointing down the main diagonal.
     * Consider this if one coordinate is vertical or equatorial.
     * Not advised if both X and Y are horizontal.
     * Value range: [-1, 1]
     */
    public static float noise2_ImproveX(long seed, double x, double y) {

        // Skew transform and rotation baked into one.
        double xx = x * ROOT2OVER2;
        double yy = y * (ROOT2OVER2 * (1 + 2 * SKEW_2D));

        return noise2_UnskewedBase(seed, yy + xx, yy - xx);
    }

    /**
     * 2D OpenSimplex2S/SuperSimplex noise unskewed base.
     * Value range: [-1, 1]
     */
    private static float noise2_UnskewedBase(long seed, double xs, double ys) {

        // Get base points and offsets.
        int xsb = fastFloor(xs), ysb = fastFloor(ys);
        float xi = (float)(xs - xsb), yi = (float)(ys - ysb);

        // Prime pre-multiplication for hash.
        long xsbp = xsb * PRIME_X, ysbp = ysb * PRIME_Y;

        // Unskew.
        float t = (xi + yi) * (float)UNSKEW_2D;
        float dx0 = xi + t, dy0 = yi + t;

        // First vertex.
        float a0 = RSQUARED_2D - dx0 * dx0 - dy0 * dy0;
        float value = (a0 * a0) * (a0 * a0) * grad(seed, xsbp, ysbp, dx0, dy0);

        // Second vertex.
        float a1 = (float)(2 * (1 + 2 * UNSKEW_2D) * (1 / UNSKEW_2D + 2)) * t + ((float)(-2 * (1 + 2 * UNSKEW_2D) * (1 + 2 * UNSKEW_2D)) + a0);
        float dx1 = dx0 - (float)(1 + 2 * UNSKEW_2D);
        float dy1 = dy0 - (float)(1 + 2 * UNSKEW_2D);
        value += (a1 * a1) * (a1 * a1) * grad(seed, xsbp + PRIME_X, ysbp + PRIME_Y, dx1, dy1);

        // Third and fourth vertices.
        // Nested conditionals were faster than compact bit logic/arithmetic.
        float xmyi = xi - yi;
        if (t < UNSKEW_2D) {
            if (xi + xmyi > 1) {
                float dx2 = dx0 - (float)(3 * UNSKEW_2D + 2);
                float dy2 = dy0 - (float)(3 * UNSKEW_2D + 1);
                float a2 = RSQUARED_2D - dx2 * dx2 - dy2 * dy2;
                if (a2 > 0) {
                    value += (a2 * a2) * (a2 * a2) * grad(seed, xsbp + (PRIME_X << 1), ysbp + PRIME_Y, dx2, dy2);
                }
            }
            else
            {
                float dx2 = dx0 - (float)UNSKEW_2D;
                float dy2 = dy0 - (float)(UNSKEW_2D + 1);
                float a2 = RSQUARED_2D - dx2 * dx2 - dy2 * dy2;
                if (a2 > 0) {
                    value += (a2 * a2) * (a2 * a2) * grad(seed, xsbp, ysbp + PRIME_Y, dx2, dy2);
                }
            }

            if (yi - xmyi > 1) {
                float dx3 = dx0 - (float)(3 * UNSKEW_2D + 1);
                float dy3 = dy0 - (float)(3 * UNSKEW_2D + 2);
                float a3 = RSQUARED_2D - dx3 * dx3 - dy3 * dy3;
                if (a3 > 0) {
                    value += (a3 * a3) * (a3 * a3) * grad(seed, xsbp + PRIME_X, ysbp + (PRIME_Y << 1), dx3, dy3);
                }
            }
            else
            {
                float dx3 = dx0 - (float)(UNSKEW_2D + 1);
                float dy3 = dy0 - (float)UNSKEW_2D;
                float a3 = RSQUARED_2D - dx3 * dx3 - dy3 * dy3;
                if (a3 > 0) {
                    value += (a3 * a3) * (a3 * a3) * grad(seed, xsbp + PRIME_X, ysbp, dx3, dy3);
                }
            }
        }
        else
        {
            if (xi + xmyi < 0) {
                float dx2 = dx0 + (float)(1 + UNSKEW_2D);
                float dy2 = dy0 + (float)UNSKEW_2D;
                float a2 = RSQUARED_2D - dx2 * dx2 - dy2 * dy2;
                if (a2 > 0) {
                    value += (a2 * a2) * (a2 * a2) * grad(seed, xsbp - PRIME_X, ysbp, dx2, dy2);
                }
            }
            else
            {
                float dx2 = dx0 - (float)(UNSKEW_2D + 1);
                float dy2 = dy0 - (float)UNSKEW_2D;
                float a2 = RSQUARED_2D - dx2 * dx2 - dy2 * dy2;
                if (a2 > 0) {
                    value += (a2 * a2) * (a2 * a2) * grad(seed, xsbp + PRIME_X, ysbp, dx2, dy2);
                }
            }

            if (yi < xmyi) {
                float dx3 = dx0 + (float)UNSKEW_2D;
                float dy3 = dy0 + (float)(UNSKEW_2D + 1);
                float a3 = RSQUARED_2D - dx3 * dx3 - dy3 * dy3;
                if (a3 > 0) {
                    value += (a3 * a3) * (a3 * a3) * grad(seed, xsbp, ysbp - PRIME_Y, dx3, dy3);
                }
            }
            else
            {
                float dx3 = dx0 - (float)UNSKEW_2D;
                float dy3 = dy0 - (float)(UNSKEW_2D + 1);
                float a3 = RSQUARED_2D - dx3 * dx3 - dy3 * dy3;
                if (a3 > 0) {
                    value += (a3 * a3) * (a3 * a3) * grad(seed, xsbp, ysbp + PRIME_Y, dx3, dy3);
                }
            }
        }

        return value * NORMALIZER_2D;
    }

    /**
     * Vector-valued 2D OpenSimplex2S/SuperSimplex noise, standard lattice orientation.
     * Vector magnitude range: [0, 1]
     */
    public static void vectorValuedNoise2(long seed, double x, double y, Vec2 destination) {

        // Get points for A2* lattice
        double s = SKEW_2D * (x + y);
        double xs = x + s, ys = y + s;

        vectorValuedNoise2_UnskewedBase(seed, xs, ys, destination);
    }

    /**
     * Vector-valued 2D OpenSimplex2S/SuperSimplex noise, with Y pointing down the main diagonal.
     * Consider this if one coordinate is vertical or equatorial.
     * Not advised if both X and Y are horizontal.
     * Magnitude range: [0, 1]
     */
    public static void vectorValuedNoise2_ImproveX(long seed, double x, double y, OpenSimplex2S_WithFrameVectorValuedEvaluators.Vec2 destination) {

        // Skew transform and rotation baked into one.
        double xx = x * ROOT2OVER2;
        double yy = y * (ROOT2OVER2 * (1 + 2 * SKEW_2D));

        vectorValuedNoise2_UnskewedBase(seed, yy + xx, yy - xx, destination);
    }

    /**
     * Vector-valued 2D OpenSimplex2S/SuperSimplex vector-output noise.
     * Vector magnitude range: [0, 1]
     */
    private static void vectorValuedNoise2_UnskewedBase(long seed, double xs, double ys, Vec2 destination) {
        destination.clear();

        // Get base points and offsets.
        int xsb = fastFloor(xs), ysb = fastFloor(ys);
        float xi = (float)(xs - xsb), yi = (float)(ys - ysb);

        // Prime pre-multiplication for hash.
        long xsbp = xsb * PRIME_X, ysbp = ysb * PRIME_Y;

        // Unskew.
        float t = (xi + yi) * (float)UNSKEW_2D;
        float dx0 = xi + t, dy0 = yi + t;

        // First vertex.
        float a0 = RSQUARED_2D - dx0 * dx0 - dy0 * dy0;
        gradFrame(seed, xsbp, ysbp, dx0, dy0, (a0 * a0) * (a0 * a0), destination);

        // Second vertex.
        float a1 = (float)(2 * (1 + 2 * UNSKEW_2D) * (1 / UNSKEW_2D + 2)) * t + ((float)(-2 * (1 + 2 * UNSKEW_2D) * (1 + 2 * UNSKEW_2D)) + a0);
        float dx1 = dx0 - (float)(1 + 2 * UNSKEW_2D);
        float dy1 = dy0 - (float)(1 + 2 * UNSKEW_2D);
        gradFrame(seed, xsbp + PRIME_X, ysbp + PRIME_Y, dx1, dy1, (a1 * a1) * (a1 * a1), destination);

        // Third and fourth vertices.
        // Nested conditionals were faster than compact bit logic/arithmetic.
        float xmyi = xi - yi;
        if (t < UNSKEW_2D) {
            if (xi + xmyi > 1) {
                float dx2 = dx0 - (float)(3 * UNSKEW_2D + 2);
                float dy2 = dy0 - (float)(3 * UNSKEW_2D + 1);
                float a2 = RSQUARED_2D - dx2 * dx2 - dy2 * dy2;
                if (a2 > 0) {
                    gradFrame(seed, xsbp + (PRIME_X << 1), ysbp + PRIME_Y, dx2, dy2, (a2 * a2) * (a2 * a2), destination);
                }
            }
            else
            {
                float dx2 = dx0 - (float)UNSKEW_2D;
                float dy2 = dy0 - (float)(UNSKEW_2D + 1);
                float a2 = RSQUARED_2D - dx2 * dx2 - dy2 * dy2;
                if (a2 > 0) {
                    gradFrame(seed, xsbp, ysbp + PRIME_Y, dx2, dy2, (a2 * a2) * (a2 * a2), destination);
                }
            }

            if (yi - xmyi > 1) {
                float dx3 = dx0 - (float)(3 * UNSKEW_2D + 1);
                float dy3 = dy0 - (float)(3 * UNSKEW_2D + 2);
                float a3 = RSQUARED_2D - dx3 * dx3 - dy3 * dy3;
                if (a3 > 0) {
                    gradFrame(seed, xsbp + PRIME_X, ysbp + (PRIME_Y << 1), dx3, dy3, (a3 * a3) * (a3 * a3), destination);
                }
            }
            else
            {
                float dx3 = dx0 - (float)(UNSKEW_2D + 1);
                float dy3 = dy0 - (float)UNSKEW_2D;
                float a3 = RSQUARED_2D - dx3 * dx3 - dy3 * dy3;
                if (a3 > 0) {
                    gradFrame(seed, xsbp + PRIME_X, ysbp, dx3, dy3, (a3 * a3) * (a3 * a3), destination);
                }
            }
        }
        else
        {
            if (xi + xmyi < 0) {
                float dx2 = dx0 + (float)(1 + UNSKEW_2D);
                float dy2 = dy0 + (float)UNSKEW_2D;
                float a2 = RSQUARED_2D - dx2 * dx2 - dy2 * dy2;
                if (a2 > 0) {
                    gradFrame(seed, xsbp - PRIME_X, ysbp, dx2, dy2, (a2 * a2) * (a2 * a2), destination);
                }
            }
            else
            {
                float dx2 = dx0 - (float)(UNSKEW_2D + 1);
                float dy2 = dy0 - (float)UNSKEW_2D;
                float a2 = RSQUARED_2D - dx2 * dx2 - dy2 * dy2;
                if (a2 > 0) {
                    gradFrame(seed, xsbp + PRIME_X, ysbp, dx2, dy2, (a2 * a2) * (a2 * a2), destination);
                }
            }

            if (yi < xmyi) {
                float dx3 = dx0 + (float)UNSKEW_2D;
                float dy3 = dy0 + (float)(UNSKEW_2D + 1);
                float a3 = RSQUARED_2D - dx3 * dx3 - dy3 * dy3;
                if (a3 > 0) {
                    gradFrame(seed, xsbp, ysbp - PRIME_Y, dx3, dy3, (a3 * a3) * (a3 * a3), destination);
                }
            }
            else
            {
                float dx3 = dx0 - (float)UNSKEW_2D;
                float dy3 = dy0 - (float)(UNSKEW_2D + 1);
                float a3 = RSQUARED_2D - dx3 * dx3 - dy3 * dy3;
                if (a3 > 0) {
                    gradFrame(seed, xsbp, ysbp + PRIME_Y, dx3, dy3, (a3 * a3) * (a3 * a3), destination);
                }
            }
        }

        destination.mul(NORMALIZER_2D);
    }

    /**
     * 3D OpenSimplex2S/SuperSimplex noise, with better visual isotropy in (X, Y).
     * Recommended for 3D terrain and time-varied animations.
     * The Z coordinate should always be the "different" coordinate in whatever your use case is.
     * If Y is vertical in world coordinates, call noise3_ImproveXZ(x, z, Y) or use noise3_XZBeforeY.
     * If Z is vertical in world coordinates, call noise3_ImproveXZ(x, y, Z).
     * For a time varied animation, call noise3_ImproveXY(x, y, T).
     * Value range: [-1, 1]
     */
    public static float noise3_ImproveXY(long seed, double x, double y, double z) {

        // Re-orient the cubic lattices without skewing, so Z points up the main lattice diagonal,
        // and the planes formed by XY are moved far out of alignment with the cube faces.
        // Orthonormal rotation. Not a skew transform.
        double xy = x + y;
        double s2 = xy * ROTATE3_ORTHOGONALIZER;
        double zz = z * ROOT3OVER3;
        double xr = x + s2 + zz;
        double yr = y + s2 + zz;
        double zr = xy * -ROOT3OVER3 + zz;

        return noise3_UnrotatedBase(seed, xr, yr, zr);
    }

    /**
     * 3D OpenSimplex2S/SuperSimplex noise, with better visual isotropy in (X, Z).
     * Recommended for 3D terrain and time-varied animations.
     * The Y coordinate should always be the "different" coordinate in whatever your use case is.
     * If Y is vertical in world coordinates, call noise3_ImproveXZ(x, Y, z).
     * If Z is vertical in world coordinates, call noise3_ImproveXZ(x, Z, y) or use noise3_ImproveXY.
     * For a time varied animation, call noise3_ImproveXZ(x, T, y) or use noise3_ImproveXY.
     * Value range: [-1, 1]
     */
    public static float noise3_ImproveXZ(long seed, double x, double y, double z) {

        // Re-orient the cubic lattices without skewing, so Y points up the main lattice diagonal,
        // and the planes formed by XZ are moved far out of alignment with the cube faces.
        // Orthonormal rotation. Not a skew transform.
        double xz = x + z;
        double s2 = xz * ROTATE3_ORTHOGONALIZER;
        double yy = y * ROOT3OVER3;
        double xr = x + s2 + yy;
        double zr = z + s2 + yy;
        double yr = xz * -ROOT3OVER3 + yy;

        return noise3_UnrotatedBase(seed, xr, yr, zr);
    }

    /**
     * 3D OpenSimplex2S/SuperSimplex noise, fallback lattice orientation.
     * Use noise3_ImproveXY or noise3_ImproveXZ instead, wherever appropriate.
     * Value range: [-1, 1]
     */
    public static float noise3_Fallback(long seed, double x, double y, double z) {

        // Re-orient the cubic lattices via rotation, to produce a familiar look.
        // Orthonormal rotation. Not a skew transform.
        double r = FALLBACK_ROTATE3 * (x + y + z);
        double xr = r - x, yr = r - y, zr = r - z;

        return noise3_UnrotatedBase(seed, xr, yr, zr);
    }

    /**
     * Generates unrotated base 3D noise using the BCC lattice.
     * Value range: [-1, 1]
     */
    private static float noise3_UnrotatedBase(long seed, double xr, double yr, double zr) {

        // Get base points and offsets.
        int xrb = fastFloor(xr), yrb = fastFloor(yr), zrb = fastFloor(zr);
        float xi = (float)(xr - xrb), yi = (float)(yr - yrb), zi = (float)(zr - zrb);

        // Prime pre-multiplication for hash. Also flip seed for second lattice copy.
        long xrbp = xrb * PRIME_X, yrbp = yrb * PRIME_Y, zrbp = zrb * PRIME_Z;
        long seed2 = seed ^ -SEED_FLIP_3D;

        // -1 if positive, 0 if negative.
        int xNMask = (int)(-0.5f - xi), yNMask = (int)(-0.5f - yi), zNMask = (int)(-0.5f - zi);

        // First vertex.
        float x0 = xi + xNMask;
        float y0 = yi + yNMask;
        float z0 = zi + zNMask;
        float a0 = RSQUARED_3D - x0 * x0 - y0 * y0 - z0 * z0;
        float value = (a0 * a0) * (a0 * a0) * grad(seed,
                xrbp + (xNMask & PRIME_X), yrbp + (yNMask & PRIME_Y), zrbp + (zNMask & PRIME_Z), x0, y0, z0);

        // Second vertex.
        float x1 = xi - 0.5f;
        float y1 = yi - 0.5f;
        float z1 = zi - 0.5f;
        float a1 = RSQUARED_3D - x1 * x1 - y1 * y1 - z1 * z1;
        value += (a1 * a1) * (a1 * a1) * grad(seed2,
                xrbp + PRIME_X, yrbp + PRIME_Y, zrbp + PRIME_Z, x1, y1, z1);

        // Shortcuts for building the remaining falloffs.
        // Derived by subtracting the polynomials with the offsets plugged in.
        float xAFlipMask0 = ((xNMask | 1) << 1) * x1;
        float yAFlipMask0 = ((yNMask | 1) << 1) * y1;
        float zAFlipMask0 = ((zNMask | 1) << 1) * z1;
        float xAFlipMask1 = (-2 - (xNMask << 2)) * x1 - 1.0f;
        float yAFlipMask1 = (-2 - (yNMask << 2)) * y1 - 1.0f;
        float zAFlipMask1 = (-2 - (zNMask << 2)) * z1 - 1.0f;

        boolean skip5 = false;
        float a2 = xAFlipMask0 + a0;
        if (a2 > 0) {
            float x2 = x0 - (xNMask | 1);
            float y2 = y0;
            float z2 = z0;
            value += (a2 * a2) * (a2 * a2) * grad(seed,
                    xrbp + (~xNMask & PRIME_X), yrbp + (yNMask & PRIME_Y), zrbp + (zNMask & PRIME_Z), x2, y2, z2);
        }
        else
        {
            float a3 = yAFlipMask0 + zAFlipMask0 + a0;
            if (a3 > 0) {
                float x3 = x0;
                float y3 = y0 - (yNMask | 1);
                float z3 = z0 - (zNMask | 1);
                value += (a3 * a3) * (a3 * a3) * grad(seed,
                        xrbp + (xNMask & PRIME_X), yrbp + (~yNMask & PRIME_Y), zrbp + (~zNMask & PRIME_Z), x3, y3, z3);
            }

            float a4 = xAFlipMask1 + a1;
            if (a4 > 0) {
                float x4 = (xNMask | 1) + x1;
                float y4 = y1;
                float z4 = z1;
                value += (a4 * a4) * (a4 * a4) * grad(seed2,
                        xrbp + (xNMask & (PRIME_X * 2)), yrbp + PRIME_Y, zrbp + PRIME_Z, x4, y4, z4);
                skip5 = true;
            }
        }

        boolean skip9 = false;
        float a6 = yAFlipMask0 + a0;
        if (a6 > 0) {
            float x6 = x0;
            float y6 = y0 - (yNMask | 1);
            float z6 = z0;
            value += (a6 * a6) * (a6 * a6) * grad(seed,
                    xrbp + (xNMask & PRIME_X), yrbp + (~yNMask & PRIME_Y), zrbp + (zNMask & PRIME_Z), x6, y6, z6);
        }
        else
        {
            float a7 = xAFlipMask0 + zAFlipMask0 + a0;
            if (a7 > 0) {
                float x7 = x0 - (xNMask | 1);
                float y7 = y0;
                float z7 = z0 - (zNMask | 1);
                value += (a7 * a7) * (a7 * a7) * grad(seed,
                        xrbp + (~xNMask & PRIME_X), yrbp + (yNMask & PRIME_Y), zrbp + (~zNMask & PRIME_Z), x7, y7, z7);
            }

            float a8 = yAFlipMask1 + a1;
            if (a8 > 0) {
                float x8 = x1;
                float y8 = (yNMask | 1) + y1;
                float z8 = z1;
                value += (a8 * a8) * (a8 * a8) * grad(seed2,
                        xrbp + PRIME_X, yrbp + (yNMask & (PRIME_Y << 1)), zrbp + PRIME_Z, x8, y8, z8);
                skip9 = true;
            }
        }

        boolean skipD = false;
        float aA = zAFlipMask0 + a0;
        if (aA > 0) {
            float xA = x0;
            float yA = y0;
            float zA = z0 - (zNMask | 1);
            value += (aA * aA) * (aA * aA) * grad(seed,
                    xrbp + (xNMask & PRIME_X), yrbp + (yNMask & PRIME_Y), zrbp + (~zNMask & PRIME_Z), xA, yA, zA);
        }
        else
        {
            float aB = xAFlipMask0 + yAFlipMask0 + a0;
            if (aB > 0) {
                float xB = x0 - (xNMask | 1);
                float yB = y0 - (yNMask | 1);
                float zB = z0;
                value += (aB * aB) * (aB * aB) * grad(seed,
                        xrbp + (~xNMask & PRIME_X), yrbp + (~yNMask & PRIME_Y), zrbp + (zNMask & PRIME_Z), xB, yB, zB);
            }

            float aC = zAFlipMask1 + a1;
            if (aC > 0) {
                float xC = x1;
                float yC = y1;
                float zC = (zNMask | 1) + z1;
                value += (aC * aC) * (aC * aC) * grad(seed2,
                        xrbp + PRIME_X, yrbp + PRIME_Y, zrbp + (zNMask & (PRIME_Z << 1)), xC, yC, zC);
                skipD = true;
            }
        }

        if (!skip5) {
            float a5 = yAFlipMask1 + zAFlipMask1 + a1;
            if (a5 > 0) {
                float x5 = x1;
                float y5 = (yNMask | 1) + y1;
                float z5 = (zNMask | 1) + z1;
                value += (a5 * a5) * (a5 * a5) * grad(seed2,
                        xrbp + PRIME_X, yrbp + (yNMask & (PRIME_Y << 1)), zrbp + (zNMask & (PRIME_Z << 1)), x5, y5, z5);
            }
        }

        if (!skip9) {
            float a9 = xAFlipMask1 + zAFlipMask1 + a1;
            if (a9 > 0) {
                float x9 = (xNMask | 1) + x1;
                float y9 = y1;
                float z9 = (zNMask | 1) + z1;
                value += (a9 * a9) * (a9 * a9) * grad(seed2,
                        xrbp + (xNMask & (PRIME_X * 2)), yrbp + PRIME_Y, zrbp + (zNMask & (PRIME_Z << 1)), x9, y9, z9);
            }
        }

        if (!skipD) {
            float aD = xAFlipMask1 + yAFlipMask1 + a1;
            if (aD > 0) {
                float xD = (xNMask | 1) + x1;
                float yD = (yNMask | 1) + y1;
                float zD = z1;
                value += (aD * aD) * (aD * aD) * grad(seed2,
                        xrbp + (xNMask & (PRIME_X << 1)), yrbp + (yNMask & (PRIME_Y << 1)), zrbp + PRIME_Z, xD, yD, zD);
            }
        }

        return value * NORMALIZER_3D;
    }

    /**
     * Vector-valued 3D OpenSimplex2S/SuperSimplex noise, with better visual isotropy in (X, Y).
     * Recommended for 3D terrain and time-varied animations.
     * The Z coordinate should always be the "different" coordinate in whatever your use case is.
     * If Y is vertical in world coordinates, call noise3_ImproveXZ(x, z, Y) or use noise3_XZBeforeY.
     * If Z is vertical in world coordinates, call noise3_ImproveXZ(x, y, Z).
     * For a time varied animation, call noise3_ImproveXY(x, y, T).
     * Vector magnitude range: [0, 1]
     */
    public static void vectorValuedNoise3_ImproveXY(long seed, double x, double y, double z, Vec3 destination) {

        // Re-orient the cubic lattices without skewing, so Z points up the main lattice diagonal,
        // and the planes formed by XY are moved far out of alignment with the cube faces.
        // Orthonormal rotation. Not a skew transform.
        double xy = x + y;
        double s2 = xy * ROTATE3_ORTHOGONALIZER;
        double zz = z * ROOT3OVER3;
        double xr = x + s2 + zz;
        double yr = y + s2 + zz;
        double zr = xy * -ROOT3OVER3 + zz;

        vectorValuedNoise3_UnrotatedBase(seed, xr, yr, zr, destination);
    }

    /**
     * Vector-valued 3D OpenSimplex2S/SuperSimplex noise, with better visual isotropy in (X, Z).
     * Recommended for 3D terrain and time-varied animations.
     * The Y coordinate should always be the "different" coordinate in whatever your use case is.
     * If Y is vertical in world coordinates, call noise3_ImproveXZ(x, Y, z).
     * If Z is vertical in world coordinates, call noise3_ImproveXZ(x, Z, y) or use noise3_ImproveXY.
     * For a time varied animation, call noise3_ImproveXZ(x, T, y) or use noise3_ImproveXY.
     * Vector magnitude range: [0, 1]
     */
    public static void vectorValuedNoise3_ImproveXZ(long seed, double x, double y, double z, Vec3 destination) {

        // Re-orient the cubic lattices without skewing, so Y points up the main lattice diagonal,
        // and the planes formed by XZ are moved far out of alignment with the cube faces.
        // Orthonormal rotation. Not a skew transform.
        double xz = x + z;
        double s2 = xz * ROTATE3_ORTHOGONALIZER;
        double yy = y * ROOT3OVER3;
        double xr = x + s2 + yy;
        double zr = z + s2 + yy;
        double yr = xz * -ROOT3OVER3 + yy;

        vectorValuedNoise3_UnrotatedBase(seed, xr, yr, zr, destination);
    }

    /**
     * Vector-valued 3D OpenSimplex2S/SuperSimplex noise, fallback lattice orientation.
     * Use noise3_ImproveXY or noise3_ImproveXZ instead, wherever appropriate.
     * Vector magnitude range: [0, 1]
     */
    public static void vectorValuedNoise3_Fallback(long seed, double x, double y, double z, Vec3 destination) {

        // Re-orient the cubic lattices via rotation, to produce a familiar look.
        // Orthonormal rotation. Not a skew transform.
        double r = FALLBACK_ROTATE3 * (x + y + z);
        double xr = r - x, yr = r - y, zr = r - z;

        vectorValuedNoise3_UnrotatedBase(seed, xr, yr, zr, destination);
    }

    /**
     * Generates unrotated base vector-valued 3D noise using the BCC lattice.
     * Vector magnitude range: [0, 1]
     */
    private static void vectorValuedNoise3_UnrotatedBase(long seed, double xr, double yr, double zr, Vec3 destination) {
        destination.clear();

        // Get base points and offsets.
        int xrb = fastFloor(xr), yrb = fastFloor(yr), zrb = fastFloor(zr);
        float xi = (float)(xr - xrb), yi = (float)(yr - yrb), zi = (float)(zr - zrb);

        // Prime pre-multiplication for hash. Also flip seed for second lattice copy.
        long xrbp = xrb * PRIME_X, yrbp = yrb * PRIME_Y, zrbp = zrb * PRIME_Z;
        long seed2 = seed ^ -SEED_FLIP_3D;

        // -1 if positive, 0 if negative.
        int xNMask = (int)(-0.5f - xi), yNMask = (int)(-0.5f - yi), zNMask = (int)(-0.5f - zi);

        // First vertex.
        float x0 = xi + xNMask;
        float y0 = yi + yNMask;
        float z0 = zi + zNMask;
        float a0 = RSQUARED_3D - x0 * x0 - y0 * y0 - z0 * z0;
        gradFrame(seed, xrbp + (xNMask & PRIME_X), yrbp + (yNMask & PRIME_Y), zrbp + (zNMask & PRIME_Z), x0, y0, z0, (a0 * a0) * (a0 * a0), destination);

        // Second vertex.
        float x1 = xi - 0.5f;
        float y1 = yi - 0.5f;
        float z1 = zi - 0.5f;
        float a1 = RSQUARED_3D - x1 * x1 - y1 * y1 - z1 * z1;
        gradFrame(seed2, xrbp + PRIME_X, yrbp + PRIME_Y, zrbp + PRIME_Z, x1, y1, z1, (a1 * a1) * (a1 * a1), destination);

        // Shortcuts for building the remaining falloffs.
        // Derived by subtracting the polynomials with the offsets plugged in.
        float xAFlipMask0 = ((xNMask | 1) << 1) * x1;
        float yAFlipMask0 = ((yNMask | 1) << 1) * y1;
        float zAFlipMask0 = ((zNMask | 1) << 1) * z1;
        float xAFlipMask1 = (-2 - (xNMask << 2)) * x1 - 1.0f;
        float yAFlipMask1 = (-2 - (yNMask << 2)) * y1 - 1.0f;
        float zAFlipMask1 = (-2 - (zNMask << 2)) * z1 - 1.0f;

        boolean skip5 = false;
        float a2 = xAFlipMask0 + a0;
        if (a2 > 0) {
            float x2 = x0 - (xNMask | 1);
            float y2 = y0;
            float z2 = z0;
            gradFrame(seed, xrbp + (~xNMask & PRIME_X), yrbp + (yNMask & PRIME_Y), zrbp + (zNMask & PRIME_Z), x2, y2, z2, (a2 * a2) * (a2 * a2), destination);
        }
        else
        {
            float a3 = yAFlipMask0 + zAFlipMask0 + a0;
            if (a3 > 0) {
                float x3 = x0;
                float y3 = y0 - (yNMask | 1);
                float z3 = z0 - (zNMask | 1);
                gradFrame(seed, xrbp + (xNMask & PRIME_X), yrbp + (~yNMask & PRIME_Y), zrbp + (~zNMask & PRIME_Z), x3, y3, z3, (a3 * a3) * (a3 * a3), destination);
            }

            float a4 = xAFlipMask1 + a1;
            if (a4 > 0) {
                float x4 = (xNMask | 1) + x1;
                float y4 = y1;
                float z4 = z1;
                gradFrame(seed2, xrbp + (xNMask & (PRIME_X * 2)), yrbp + PRIME_Y, zrbp + PRIME_Z, x4, y4, z4, (a4 * a4) * (a4 * a4), destination);
                skip5 = true;
            }
        }

        boolean skip9 = false;
        float a6 = yAFlipMask0 + a0;
        if (a6 > 0) {
            float x6 = x0;
            float y6 = y0 - (yNMask | 1);
            float z6 = z0;
            gradFrame(seed, xrbp + (xNMask & PRIME_X), yrbp + (~yNMask & PRIME_Y), zrbp + (zNMask & PRIME_Z), x6, y6, z6, (a6 * a6) * (a6 * a6), destination);
        }
        else
        {
            float a7 = xAFlipMask0 + zAFlipMask0 + a0;
            if (a7 > 0) {
                float x7 = x0 - (xNMask | 1);
                float y7 = y0;
                float z7 = z0 - (zNMask | 1);
                gradFrame(seed, xrbp + (~xNMask & PRIME_X), yrbp + (yNMask & PRIME_Y), zrbp + (~zNMask & PRIME_Z), x7, y7, z7, (a7 * a7) * (a7 * a7), destination);
            }

            float a8 = yAFlipMask1 + a1;
            if (a8 > 0) {
                float x8 = x1;
                float y8 = (yNMask | 1) + y1;
                float z8 = z1;
                gradFrame(seed2, xrbp + PRIME_X, yrbp + (yNMask & (PRIME_Y << 1)), zrbp + PRIME_Z, x8, y8, z8, (a8 * a8) * (a8 * a8), destination);
                skip9 = true;
            }
        }

        boolean skipD = false;
        float aA = zAFlipMask0 + a0;
        if (aA > 0) {
            float xA = x0;
            float yA = y0;
            float zA = z0 - (zNMask | 1);
            gradFrame(seed, xrbp + (xNMask & PRIME_X), yrbp + (yNMask & PRIME_Y), zrbp + (~zNMask & PRIME_Z), xA, yA, zA, (aA * aA) * (aA * aA), destination);
        }
        else
        {
            float aB = xAFlipMask0 + yAFlipMask0 + a0;
            if (aB > 0) {
                float xB = x0 - (xNMask | 1);
                float yB = y0 - (yNMask | 1);
                float zB = z0;
                gradFrame(seed, xrbp + (~xNMask & PRIME_X), yrbp + (~yNMask & PRIME_Y), zrbp + (zNMask & PRIME_Z), xB, yB, zB, (aB * aB) * (aB * aB), destination);
            }

            float aC = zAFlipMask1 + a1;
            if (aC > 0) {
                float xC = x1;
                float yC = y1;
                float zC = (zNMask | 1) + z1;
                gradFrame(seed2, xrbp + PRIME_X, yrbp + PRIME_Y, zrbp + (zNMask & (PRIME_Z << 1)), xC, yC, zC, (aC * aC) * (aC * aC), destination);
                skipD = true;
            }
        }

        if (!skip5) {
            float a5 = yAFlipMask1 + zAFlipMask1 + a1;
            if (a5 > 0) {
                float x5 = x1;
                float y5 = (yNMask | 1) + y1;
                float z5 = (zNMask | 1) + z1;
                gradFrame(seed2, xrbp + PRIME_X, yrbp + (yNMask & (PRIME_Y << 1)), zrbp + (zNMask & (PRIME_Z << 1)), x5, y5, z5, (a5 * a5) * (a5 * a5), destination);
            }
        }

        if (!skip9) {
            float a9 = xAFlipMask1 + zAFlipMask1 + a1;
            if (a9 > 0) {
                float x9 = (xNMask | 1) + x1;
                float y9 = y1;
                float z9 = (zNMask | 1) + z1;
                gradFrame(seed2, xrbp + (xNMask & (PRIME_X * 2)), yrbp + PRIME_Y, zrbp + (zNMask & (PRIME_Z << 1)), x9, y9, z9, (a9 * a9) * (a9 * a9), destination);
            }
        }

        if (!skipD) {
            float aD = xAFlipMask1 + yAFlipMask1 + a1;
            if (aD > 0) {
                float xD = (xNMask | 1) + x1;
                float yD = (yNMask | 1) + y1;
                float zD = z1;
                gradFrame(seed2, xrbp + (xNMask & (PRIME_X << 1)), yrbp + (yNMask & (PRIME_Y << 1)), zrbp + PRIME_Z, xD, yD, zD, (aD * aD) * (aD * aD), destination);
            }
        }

        destination.mul(NORMALIZER_3D);
    }

    /*
     * Utility
     */

    private static float grad(long seed, long xsvp, long ysvp, float dx, float dy) {
        long hash = seed ^ xsvp ^ ysvp;
        hash *= HASH_MULTIPLIER;
        hash ^= hash >> (64 - N_SINES_360_EXP);

        float gx = SINES_360[ (int)hash                  & (N_SINES_360 - 1)];
        float gy = SINES_360[((int)hash + COSINE_OFFSET) & (N_SINES_360 - 1)];

        return gx * dx + gy * dy;
    }

    private static void gradFrame(long seed, long xsvp, long ysvp, float dx, float dy, float falloff, Vec2 destination) {
        long hash = seed ^ xsvp ^ ysvp;
        hash *= HASH_MULTIPLIER;
        hash ^= hash >> (64 - (N_SINES_360_EXP + 1));

        float a = SINES_360[ (int)hash                  & (N_SINES_360 - 1)];
        float b = SINES_360[((int)hash + COSINE_OFFSET) & (N_SINES_360 - 1)];
        float rampX = a * dx + b * dy;
        float rampY = b * dx - a * dy;
        rampY = ((hash & N_SINES_360) != 0) ? -rampY : rampY;

        destination.add(falloff * rampX, falloff * rampY);
    }

    private static float grad(long seed, long xrvp, long yrvp, long zrvp, float dx, float dy, float dz) {
        long hash = seed ^ xrvp ^ yrvp ^ zrvp;
        hash *= HASH_MULTIPLIER;
        hash ^= hash >> (64 - (N_SINES_360_EXP + N_VECTOR_SCALES_3D_180_EXP));

        // Unit-spherical point sampling.
        int iScale = (int)(hash >> (N_SINES_360_EXP - 1)) & ((N_VECTOR_SCALES_3D_180 - 1) << 1);
        float scale = VECTOR_SCALES_3D_180[iScale];
        float gx = SINES_360[ (int)hash                  & (N_SINES_360 - 1)] * scale;
        float gy = VECTOR_SCALES_3D_180[iScale | 1];
        float gz = SINES_360[((int)hash + COSINE_OFFSET) & (N_SINES_360 - 1)] * scale;

        return gx * dx + gy * dy + gz * dz;
    }

    private static void gradFrame(long seed, long xrvp, long yrvp, long zrvp, float dx, float dy, float dz, float falloff, Vec3 destination) {
        long hash = seed ^ xrvp ^ yrvp ^ zrvp;
        hash *= HASH_MULTIPLIER;
        hash ^= hash >> (64 - (N_SINES_360_EXP + N_SINES_360_EXP + N_ROTOR_SLIDES_3D_90_EXP + 1));

        // Rotor/quaternion sampling.
        // Rotors info/video: https://marctenbosch.com/quaternions/
        // Sampling based on https://stackoverflow.com/a/56794499 but with scaled sin/cos instead of rejection.
        int iSlide = (int)hash & ((N_ROTOR_SLIDES_3D_90 - 1) << 1);
        float slide0 = ROTOR_SLIDES_3D_90[iSlide    ];
        float slide1 = ROTOR_SLIDES_3D_90[iSlide | 1];
        int hAB = (int)(hash >> (N_ROTOR_SLIDES_3D_90_EXP - 1));
        int hCD = (int)(hash >> (N_ROTOR_SLIDES_3D_90_EXP - 1 + N_SINES_360_EXP));
        float a = SINES_360[ hAB                  & (N_SINES_360 - 1)] * slide0;
        float b = SINES_360[ hCD                  & (N_SINES_360 - 1)] * slide1;
        float c = SINES_360[(hAB + COSINE_OFFSET) & (N_SINES_360 - 1)] * slide0;
        float d = SINES_360[(hCD + COSINE_OFFSET) & (N_SINES_360 - 1)] * slide1;

        // Apply rotation to displacement vector.
        // https://blog.molecular-matters.com/2013/05/24/a-faster-quaternion-vector-multiplication/
        float tx = 2 * (b * dz - c * dy);
        float ty = 2 * (c * dx - a * dz);
        float tz = 2 * (a * dy - b * dx);
        float rampX = dx + d * tx + (b * tz - c * ty);
        float rampY = dy + d * ty + (c * tx - a * tz);
        float rampZ = dz + d * tz + (a * ty - b * tx);

        // To get all possible orthogonal frames (within precision),
        // and not just those that follow the right-hand rule,
        // it suffices to optionally flip just one of the components.
        rampY = ((hash & 1) != 0) ? -rampY : rampY;

        destination.add(falloff * rampX, falloff * rampY, falloff * rampZ);
    }

    private static int fastFloor(double x) {
        int xi = (int)x;
        return x < xi ? xi - 1 : xi;
    }

    /*
     * Gradient and frame construction tables
     */

    private static final float[] SINES_360;
    private static final float[] ROTOR_SLIDES_3D_90;
    private static final float[] VECTOR_SCALES_3D_180;
    static {
        SINES_360 = new float[N_SINES_360];
        for (int i = 0; i < N_SINES_360; i++) {
            double angle = (i + 0.5) * (Math.PI * 2.0 / N_SINES_360);
            SINES_360[i] = (float)Math.sin(angle);
        }

        ROTOR_SLIDES_3D_90 = new float[N_ROTOR_SLIDES_3D_90 * 2];
        for (int i = 0; i < N_ROTOR_SLIDES_3D_90; i++) {
            double t = (i + 0.5) * (1.0 / N_ROTOR_SLIDES_3D_90);
            ROTOR_SLIDES_3D_90[i * 2    ] = (float)Math.sqrt(1 - t);
            ROTOR_SLIDES_3D_90[i * 2 + 1] = (float)Math.sqrt(t);
        }

        VECTOR_SCALES_3D_180 = new float[N_VECTOR_SCALES_3D_180 * 2];
        for (int i = 0; i < N_VECTOR_SCALES_3D_180; i++) {
            double t = (i + (-N_VECTOR_SCALES_3D_180 / 2.0 + 0.5)) * (2.0 / N_VECTOR_SCALES_3D_180);
            VECTOR_SCALES_3D_180[i * 2    ] = (float)Math.sqrt(1 - t * t);
            VECTOR_SCALES_3D_180[i * 2 + 1] = (float)t;
        }
    }

    public static class Vec2 {
        public float x, y;
        public void set(float x, float y) {
            this.x = x;
            this.y = y;
        }
        public void add(float x, float y) {
            this.x += x;
            this.y += y;
        }
        public void mul(float t) {
            this.x *= t;
            this.y *= t;
        }
        public void clear() {
            this.x = this.y = 0;
        }
    }

    public static class Vec3 {
        public float x, y, z;
        public void set(float x, float y, float z) {
            this.x = x;
            this.y = y;
            this.z = z;
        }
        public void add(float x, float y, float z) {
            this.x += x;
            this.y += y;
            this.z += z;
        }
        public void mul(float t) {
            this.x *= t;
            this.y *= t;
            this.z *= t;
        }
        public void clear() {
            this.x = this.y = this.z = 0;
        }
    }
}
