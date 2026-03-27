/**
 * OpenSimplex2S noise with
 * - orthonormal-frame-based vector-valued gradient noise endpoints, and
 * - uniform table indexing
 */

public class OpenSimplex2S_WithVectorValuedSamplers {

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

    private static final int N_GRADIENTS_2D_BASE_EXPONENT = 3;
    private static final int N_GRADIENTS_2D_MULTIPLIER = 3;
    private static final int N_GRADIENTS_2D_MULTIPLIER_COVERAGE_EXPONENT = 2;
    private static final int N_GRADIENTS_2D = (1 << N_GRADIENTS_2D_BASE_EXPONENT) * N_GRADIENTS_2D_MULTIPLIER;
    private static final int GRADIENTS_2D_BITS_CONSIDERED = 31 - N_GRADIENTS_2D_MULTIPLIER_COVERAGE_EXPONENT;
    private static final int GRADIENTS_2D_HASH_MASK = (1 << GRADIENTS_2D_BITS_CONSIDERED) - 1;
    private static final int GRADIENT_VECTOR_2D_LENGTH_EXPONENT = 1;
    private static final int GRADIENT_VECTOR_2D_LENGTH = 1 << GRADIENT_VECTOR_2D_LENGTH_EXPONENT;
    private static final int GRADIENTS_2D_BIT_SHIFT = GRADIENTS_2D_BITS_CONSIDERED - N_GRADIENTS_2D_BASE_EXPONENT - GRADIENT_VECTOR_2D_LENGTH_EXPONENT;
    private static final int GRADIENTS_2D_SELECTION_MASK = ((1 << (N_GRADIENTS_2D_BASE_EXPONENT + N_GRADIENTS_2D_MULTIPLIER_COVERAGE_EXPONENT)) - 1) << GRADIENT_VECTOR_2D_LENGTH_EXPONENT;

    private static final int N_GRADIENTS_3D_BASE_EXPONENT = 4;
    private static final int N_GRADIENTS_3D_MULTIPLIER = 3;
    private static final int N_GRADIENTS_3D_MULTIPLIER_COVERAGE_EXPONENT = 2;
    private static final int N_GRADIENTS_3D = (1 << N_GRADIENTS_3D_BASE_EXPONENT) * N_GRADIENTS_3D_MULTIPLIER;
    private static final int GRADIENTS_3D_BITS_CONSIDERED = 31 - N_GRADIENTS_3D_MULTIPLIER_COVERAGE_EXPONENT;
    private static final int GRADIENTS_3D_HASH_MASK = (1 << GRADIENTS_3D_BITS_CONSIDERED) - 1;
    private static final int GRADIENT_VECTOR_3D_LENGTH = 3;
    private static final int GRADIENT_VECTOR_3D_LENGTH_PADDED_EXPONENT = 2;
    private static final int GRADIENT_VECTOR_3D_LENGTH_PADDED = 1 << GRADIENT_VECTOR_3D_LENGTH_PADDED_EXPONENT;
    private static final int GRADIENTS_3D_BIT_SHIFT = GRADIENTS_3D_BITS_CONSIDERED - N_GRADIENTS_3D_BASE_EXPONENT - GRADIENT_VECTOR_3D_LENGTH_PADDED_EXPONENT;
    private static final int GRADIENTS_3D_SELECTION_MASK = ((1 << (N_GRADIENTS_3D_BASE_EXPONENT + N_GRADIENTS_3D_MULTIPLIER_COVERAGE_EXPONENT)) - 1) << GRADIENT_VECTOR_3D_LENGTH_PADDED_EXPONENT;

    private static final int N_FRAMES_2D_BASE_EXPONENT = 4;
    private static final int N_FRAMES_2D_MULTIPLIER = 3;
    private static final int N_FRAMES_2D_MULTIPLIER_COVERAGE_EXPONENT = 2;
    private static final int N_FRAMES_2D = (1 << N_FRAMES_2D_BASE_EXPONENT) * N_FRAMES_2D_MULTIPLIER;
    private static final int FRAMES_2D_BITS_CONSIDERED = 31 - N_FRAMES_2D_MULTIPLIER_COVERAGE_EXPONENT;
    private static final int FRAMES_2D_HASH_MASK = (1 << FRAMES_2D_BITS_CONSIDERED) - 1;
    private static final int FRAME_MATRIX_2D_LENGTH_EXPONENT = 2;
    private static final int FRAME_MATRIX_2D_LENGTH = 1 << FRAME_MATRIX_2D_LENGTH_EXPONENT;
    private static final int FRAMES_2D_BIT_SHIFT = FRAMES_2D_BITS_CONSIDERED - N_FRAMES_2D_BASE_EXPONENT - FRAME_MATRIX_2D_LENGTH_EXPONENT;
    private static final int FRAMES_2D_SELECTION_MASK = ((1 << (N_FRAMES_2D_BASE_EXPONENT + N_FRAMES_2D_MULTIPLIER_COVERAGE_EXPONENT)) - 1) << FRAME_MATRIX_2D_LENGTH_EXPONENT;

    private static final int N_FRAMES_3D_BASE_EXPONENT = 7;
    private static final int N_FRAMES_3D_MULTIPLIER = 9;
    private static final int N_FRAMES_3D_MULTIPLIER_COVERAGE_EXPONENT = 4;
    private static final int N_FRAMES_3D = (1 << N_FRAMES_3D_BASE_EXPONENT) * N_FRAMES_3D_MULTIPLIER;
    private static final int FRAMES_3D_BITS_CONSIDERED = 31 - N_FRAMES_3D_MULTIPLIER_COVERAGE_EXPONENT;
    private static final int FRAMES_3D_HASH_MASK = (1 << FRAMES_3D_BITS_CONSIDERED) - 1;
    private static final int FRAME_MATRIX_3D_LENGTH = 9;
    private static final int FRAMES_3D_BIT_SHIFT = FRAMES_3D_BITS_CONSIDERED - N_FRAMES_3D_BASE_EXPONENT;

    private static final double GRADIENT_NORMALIZATION_DIVISOR_2D = 0.05481866495625118;
    private static final double GRADIENT_NORMALIZATION_DIVISOR_3D = 0.2781926117527186;
    
    private static final double FRAME_NORMALIZATION_DIVISOR_2D = 0.05518041029835581;
    private static final double FRAME_NORMALIZATION_DIVISOR_3D = 0.08677623918785461;

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
                float dx2 = dx0 + (float)UNSKEW_2D;
                float dy2 = dy0 + (float)(UNSKEW_2D + 1);
                float a2 = RSQUARED_2D - dx2 * dx2 - dy2 * dy2;
                if (a2 > 0) {
                    value += (a2 * a2) * (a2 * a2) * grad(seed, xsbp, ysbp - PRIME_Y, dx2, dy2);
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
        }

        return value;
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
    public static void vectorValuedNoise2_ImproveX(long seed, double x, double y, Vec2 destination) {

        // Skew transform and rotation baked into one.
        double xx = x * ROOT2OVER2;
        double yy = y * (ROOT2OVER2 * (1 + 2 * SKEW_2D));

        vectorValuedNoise2_UnskewedBase(seed, yy + xx, yy - xx, destination);
    }

    /**
     * Vector-valued 2D OpenSimplex2S/SuperSimplex noise unskewed base.
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

        // Evaluate both lattices to form a BCC lattice.
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
        double s2 = xz * -0.211324865405187;
        double yy = y * ROOT3OVER3;
        double xr = x + s2 + yy;
        double zr = z + s2 + yy;
        double yr = xz * -ROOT3OVER3 + yy;

        // Evaluate both lattices to form a BCC lattice.
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

        // Evaluate both lattices to form a BCC lattice.
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
        long seed2 = seed ^ SEED_FLIP_3D;

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

        return value;
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
    }

    /*
     * Utility
     */

    private static float grad(long seed, long xsvp, long ysvp, float dx, float dy) {
        long hash = seed ^ xsvp ^ ysvp;
        hash *= HASH_MULTIPLIER;
        hash ^= hash >> (64 - GRADIENTS_2D_BITS_CONSIDERED);

        int index = (int)hash & GRADIENTS_2D_HASH_MASK;
        index = (index * N_GRADIENTS_2D_MULTIPLIER) >> GRADIENTS_2D_BIT_SHIFT;
        index &= GRADIENTS_2D_SELECTION_MASK;

        return GRADIENTS_2D[index | 0] * dx + GRADIENTS_2D[index | 1] * dy;
    }

    private static void gradFrame(long seed, long xrvp, long yrvp, float dx, float dy, float falloff, Vec2 destination) {
        long hash = seed ^ xrvp ^ yrvp;
        hash *= HASH_MULTIPLIER;
        hash ^= hash >> (64 - FRAMES_2D_BITS_CONSIDERED);

        int index = (int)hash & FRAMES_2D_HASH_MASK;
        index = (index * N_FRAMES_2D_MULTIPLIER) >> FRAMES_2D_BIT_SHIFT;
        index &= FRAMES_2D_SELECTION_MASK;

        float rampX = FRAMES_2D[index | 0] * dx + FRAMES_2D[index | 1] * dy;
        float rampY = FRAMES_2D[index | 2] * dx + FRAMES_2D[index | 3] * dy;

        destination.add(falloff * rampX, falloff * rampY);
    }

    private static float grad(long seed, long xrvp, long yrvp, long zrvp, float dx, float dy, float dz) {
        long hash = (seed ^ xrvp) ^ (yrvp ^ zrvp);
        hash *= HASH_MULTIPLIER;
        hash ^= hash >> (64 - GRADIENTS_3D_BITS_CONSIDERED);

        int index = (int)hash & GRADIENTS_3D_HASH_MASK;
        index = (index * N_GRADIENTS_3D_MULTIPLIER) >> GRADIENTS_3D_BIT_SHIFT;
        index &= GRADIENTS_3D_SELECTION_MASK;

        return GRADIENTS_3D[index | 0] * dx + GRADIENTS_3D[index | 1] * dy + GRADIENTS_3D[index | 2] * dz;
    }

    private static void gradFrame(long seed, long xrvp, long yrvp, long zrvp, float dx, float dy, float dz, float falloff, Vec3 destination) {
        long hash = seed ^ xrvp ^ yrvp ^ zrvp;
        hash *= HASH_MULTIPLIER;
        hash ^= hash >> (64 - FRAMES_3D_BITS_CONSIDERED);

        int index = (int)hash & FRAMES_3D_HASH_MASK;
        index = (index * N_FRAMES_3D_MULTIPLIER) >> FRAMES_3D_BIT_SHIFT;
        index *= FRAME_MATRIX_3D_LENGTH;

        float rampX = FRAMES_3D[index + 0] * dx + FRAMES_3D[index + 1] * dy + FRAMES_3D[index + 2] * dz;
        float rampY = FRAMES_3D[index + 3] * dx + FRAMES_3D[index + 4] * dy + FRAMES_3D[index + 5] * dz;
        float rampZ = FRAMES_3D[index + 6] * dx + FRAMES_3D[index + 7] * dy + FRAMES_3D[index + 8] * dz;

        destination.add(falloff * rampX, falloff * rampY, falloff * rampZ);
    }

    private static int fastFloor(double x) {
        int xi = (int)x;
        return x < xi ? xi - 1 : xi;
    }

    /*
     * Lookup Tables & Gradients
     */

    private static final float[] GRADIENTS_2D;
    private static final float[] GRADIENTS_3D;
    private static final float[] FRAMES_2D;
    private static final float[] FRAMES_3D;
    static {

        GRADIENTS_2D = new float[N_GRADIENTS_2D * GRADIENT_VECTOR_2D_LENGTH];
        double[] grad2 = {
                 0.38268343236509,   0.923879532511287,
                 0.923879532511287,  0.38268343236509,
                 0.923879532511287, -0.38268343236509,
                 0.38268343236509,  -0.923879532511287,
                -0.38268343236509,  -0.923879532511287,
                -0.923879532511287, -0.38268343236509,
                -0.923879532511287,  0.38268343236509,
                -0.38268343236509,   0.923879532511287,
                 0.130526192220052,  0.99144486137381,
                 0.608761429008721,  0.793353340291235,
                 0.793353340291235,  0.608761429008721,
                 0.99144486137381,   0.130526192220051,
                 0.99144486137381,  -0.130526192220051,
                 0.793353340291235, -0.60876142900872,
                 0.608761429008721, -0.793353340291235,
                 0.130526192220052, -0.99144486137381,
                -0.130526192220052, -0.99144486137381,
                -0.608761429008721, -0.793353340291235,
                -0.793353340291235, -0.608761429008721,
                -0.99144486137381,  -0.130526192220052,
                -0.99144486137381,   0.130526192220051,
                -0.793353340291235,  0.608761429008721,
                -0.608761429008721,  0.793353340291235,
                -0.130526192220052,  0.99144486137381,
        };
        for (int i = 0; i < GRADIENTS_2D.length; i++) {
            GRADIENTS_2D[i] = (float)(grad2[i] / GRADIENT_NORMALIZATION_DIVISOR_2D);
        }

        GRADIENTS_3D = new float[N_GRADIENTS_3D * GRADIENT_VECTOR_3D_LENGTH_PADDED];
        double[] grad3 = {
                 2.22474487139,       2.22474487139,      -1.0,
                 2.22474487139,       2.22474487139,       1.0,
                 3.0862664687972017,  1.1721513422464978,  0.0,
                 1.1721513422464978,  3.0862664687972017,  0.0,
                -2.22474487139,       2.22474487139,      -1.0,
                -2.22474487139,       2.22474487139,       1.0,
                -1.1721513422464978,  3.0862664687972017,  0.0,
                -3.0862664687972017,  1.1721513422464978,  0.0,
                -1.0,                -2.22474487139,      -2.22474487139,
                 1.0,                -2.22474487139,      -2.22474487139,
                 0.0,                -3.0862664687972017, -1.1721513422464978,
                 0.0,                -1.1721513422464978, -3.0862664687972017,
                -1.0,                -2.22474487139,       2.22474487139,
                 1.0,                -2.22474487139,       2.22474487139,
                 0.0,                -1.1721513422464978,  3.0862664687972017,
                 0.0,                -3.0862664687972017,  1.1721513422464978,
                -2.22474487139,      -2.22474487139,      -1.0,
                -2.22474487139,      -2.22474487139,       1.0,
                -3.0862664687972017, -1.1721513422464978,  0.0,
                -1.1721513422464978, -3.0862664687972017,  0.0,
                -2.22474487139,      -1.0,                -2.22474487139,
                -2.22474487139,       1.0,                -2.22474487139,
                -1.1721513422464978,  0.0,                -3.0862664687972017,
                -3.0862664687972017,  0.0,                -1.1721513422464978,
                -2.22474487139,      -1.0,                 2.22474487139,
                -2.22474487139,       1.0,                 2.22474487139,
                -3.0862664687972017,  0.0,                 1.1721513422464978,
                -1.1721513422464978,  0.0,                 3.0862664687972017,
                -1.0,                 2.22474487139,      -2.22474487139,
                 1.0,                 2.22474487139,      -2.22474487139,
                 0.0,                 1.1721513422464978, -3.0862664687972017,
                 0.0,                 3.0862664687972017, -1.1721513422464978,
                -1.0,                 2.22474487139,       2.22474487139,
                 1.0,                 2.22474487139,       2.22474487139,
                 0.0,                 3.0862664687972017,  1.1721513422464978,
                 0.0,                 1.1721513422464978,  3.0862664687972017,
                 2.22474487139,      -2.22474487139,      -1.0,
                 2.22474487139,      -2.22474487139,       1.0,
                 1.1721513422464978, -3.0862664687972017,  0.0,
                 3.0862664687972017, -1.1721513422464978,  0.0,
                 2.22474487139,      -1.0,                -2.22474487139,
                 2.22474487139,       1.0,                -2.22474487139,
                 3.0862664687972017,  0.0,                -1.1721513422464978,
                 1.1721513422464978,  0.0,                -3.0862664687972017,
                 2.22474487139,      -1.0,                 2.22474487139,
                 2.22474487139,       1.0,                 2.22474487139,
                 1.1721513422464978,  0.0,                 3.0862664687972017,
                 3.0862664687972017,  0.0,                 1.1721513422464978,
        };
        for (int i = 0; i < grad3.length; i++) {
            grad3[i] /= GRADIENT_NORMALIZATION_DIVISOR_3D;
        }
        for (int indexBase = 0; indexBase < N_GRADIENTS_3D; indexBase++) {
            int indexSource = indexBase * GRADIENT_VECTOR_3D_LENGTH;
            int indexDestination = indexBase * GRADIENT_VECTOR_3D_LENGTH_PADDED;
            for (int iAxis = 0; iAxis < GRADIENT_VECTOR_3D_LENGTH; iAxis++) {
                GRADIENTS_3D[indexDestination + iAxis] = (float)grad3[indexSource + iAxis];
            }
        }
        
        FRAMES_2D = new float[N_FRAMES_2D * FRAME_MATRIX_2D_LENGTH];
        double[] frame2 = {
                 0.79335334029123539,  0.60876142900872054,  0.60876142900872088, -0.79335334029123517,
                 0.60876142900872088, -0.79335334029123505,  0.79335334029123539,  0.60876142900872066,
                 0.79335334029123539, -0.60876142900872054,  0.60876142900872088,  0.79335334029123517,
                 0.60876142900872088,  0.79335334029123505,  0.79335334029123539, -0.60876142900872066,
                 0.38268343236508978,  0.92387953251128674,  0.92387953251128685, -0.38268343236508989,
                 0.92387953251128674, -0.38268343236508978,  0.38268343236508973,  0.92387953251128696,
                 0.99144486137381038, -0.13052619222005157,  0.13052619222005185,  0.99144486137381038,
                 0.13052619222005191,  0.99144486137381027,  0.99144486137381049, -0.13052619222005168,
                -0.13052619222005146,  0.99144486137381049,  0.99144486137381049,  0.13052619222005157,
                 0.99144486137381038,  0.13052619222005163, -0.13052619222005157,  0.99144486137381049,
                 0.92387953251128674,  0.38268343236508984, -0.38268343236508973,  0.92387953251128674,
                -0.38268343236508962,  0.92387953251128674,  0.92387953251128685,  0.38268343236508978,
                -0.60876142900872054,  0.79335334029123550,  0.79335334029123539,  0.60876142900872077,
                 0.79335334029123528,  0.60876142900872077, -0.60876142900872066,  0.79335334029123550,
                 0.60876142900872077,  0.79335334029123517, -0.79335334029123517,  0.60876142900872088,
                -0.79335334029123505,  0.60876142900872088,  0.60876142900872088,  0.79335334029123517,
                -0.92387953251128674,  0.38268343236509006,  0.38268343236509006,  0.92387953251128674,
                 0.38268343236509006,  0.92387953251128674, -0.92387953251128674,  0.38268343236509006,
                 0.13052619222005191,  0.99144486137381027, -0.99144486137381038,  0.13052619222005191,
                -0.99144486137381038,  0.13052619222005202,  0.13052619222005191,  0.99144486137381038,
                -0.99144486137381072, -0.13052619222005124, -0.13052619222005152,  0.99144486137381060,
                -0.13052619222005163,  0.99144486137381060, -0.99144486137381072, -0.13052619222005124,
                -0.38268343236508962,  0.92387953251128674, -0.92387953251128696, -0.38268343236508967,
                -0.92387953251128696, -0.38268343236508956, -0.38268343236508956,  0.92387953251128696,
                -0.79335334029123550, -0.60876142900872021, -0.60876142900872032,  0.79335334029123550,
                -0.60876142900872032,  0.79335334029123539, -0.79335334029123550, -0.60876142900872032,
                -0.79335334029123517,  0.60876142900872077, -0.60876142900872088, -0.79335334029123494,
                -0.60876142900872088, -0.79335334029123483, -0.79335334029123517,  0.60876142900872088,
                -0.38268343236509056, -0.92387953251128663, -0.92387953251128674,  0.38268343236509017,
                -0.92387953251128674,  0.38268343236509017, -0.38268343236509056, -0.92387953251128663,
                -0.99144486137381038,  0.13052619222005191, -0.13052619222005185, -0.99144486137381038,
                -0.13052619222005191, -0.99144486137381038, -0.99144486137381049,  0.13052619222005191,
                 0.13052619222005118, -0.99144486137381049, -0.99144486137381049, -0.13052619222005124,
                -0.99144486137381049, -0.13052619222005130,  0.13052619222005130, -0.99144486137381060,
                -0.92387953251128696, -0.38268343236508928,  0.38268343236508928, -0.92387953251128718,
                 0.38268343236508917, -0.92387953251128707, -0.92387953251128707, -0.38268343236508923,
                 0.60876142900871999, -0.79335334029123572, -0.79335334029123594, -0.60876142900871988,
                -0.79335334029123583, -0.60876142900871988,  0.60876142900872010, -0.79335334029123572,
                -0.60876142900872088, -0.79335334029123505,  0.79335334029123505, -0.60876142900872099,
                 0.79335334029123494, -0.60876142900872099, -0.60876142900872099, -0.79335334029123505,
                 0.92387953251128652, -0.38268343236509056, -0.38268343236509023, -0.92387953251128674,
                -0.38268343236509023, -0.92387953251128674,  0.92387953251128652, -0.38268343236509056,
                -0.13052619222005257, -0.99144486137381027,  0.99144486137381038, -0.13052619222005257,
                 0.99144486137381038, -0.13052619222005263, -0.13052619222005257, -0.99144486137381038,
                 0.99144486137381027,  0.13052619222005118,  0.13052619222005107, -0.99144486137381049,
                 0.13052619222005118, -0.99144486137381049,  0.99144486137381027,  0.13052619222005130,
                 0.38268343236508917, -0.92387953251128696,  0.92387953251128718,  0.38268343236508923,
                 0.92387953251128718,  0.38268343236508912,  0.38268343236508917, -0.92387953251128718
        };
        for (int i = 0; i < FRAMES_2D.length; i++) {
            FRAMES_2D[i] = (float)(frame2[i] / FRAME_NORMALIZATION_DIVISOR_2D);
        }
        
        FRAMES_3D = new float[N_FRAMES_3D * FRAME_MATRIX_3D_LENGTH];
        populateFrames3D_Part1();
        populateFrames3D_Part2();
        populateFrames3D_Part3();
        populateFrames3D_Part4();
    }

    private static void populateFrames3D(int startingIndex, double[] unnormalizedValues) {
        for (int i = 0; i < unnormalizedValues.length; i++) {
            FRAMES_3D[startingIndex + i] = (float)(unnormalizedValues[i] / FRAME_NORMALIZATION_DIVISOR_3D);
        }
    }

    private static void populateFrames3D_Part1() {
        populateFrames3D(0, new double[] {
             0.37627949404651073, -0.21387143172584849,  0.90148363992456226,  0.72510924905369167,  0.67365669613657153, -0.14283988479115478, -0.57674111982491216,  0.70742184475876546,  0.40856335403608102,
             0.72510924905369167,  0.67365669613657131, -0.14283988479115467,  0.37627949404651057, -0.21387143172584838,  0.90148363992456226, -0.57674111982491194,  0.70742184475876568,  0.40856335403608085,
             0.92650620200843870,  0.36957471874357806,  0.07071622802089356, -0.02651441186298362, -0.12334556714490336,  0.99200950450550773, -0.37534416687016486,  0.92097795757081369,  0.10448137664308749,
            -0.02651441186298353, -0.12334556714490325,  0.99200950450550751,  0.92650620200843870,  0.36957471874357806,  0.07071622802089356, -0.37534416687016486,  0.92097795757081369,  0.10448137664308760,
             0.37627949404651051,  0.90148363992456249, -0.21387143172584849,  0.72510924905369190, -0.14283988479115489,  0.67365669613657153, -0.57674111982491194,  0.40856335403608113,  0.70742184475876546,
             0.72510924905369190, -0.14283988479115489,  0.67365669613657153,  0.37627949404651040,  0.90148363992456249, -0.21387143172584849, -0.57674111982491194,  0.40856335403608113,  0.70742184475876546,
             0.92650620200843870,  0.07071622802089345,  0.36957471874357817, -0.02651441186298364,  0.99200950450550751, -0.12334556714490336, -0.37534416687016475,  0.10448137664308760,  0.92097795757081391,
            -0.02651441186298367,  0.99200950450550751, -0.12334556714490347,  0.92650620200843892,  0.07071622802089339,  0.36957471874357828, -0.37534416687016497,  0.10448137664308771,  0.92097795757081369,
            -0.37534416687016453,  0.10448137664308754,  0.92097795757081369,  0.92650620200843892,  0.07071622802089333,  0.36957471874357828, -0.02651441186298392,  0.99200950450550773, -0.12334556714490325,
             0.92650620200843870,  0.07071622802089345,  0.36957471874357822, -0.37534416687016475,  0.10448137664308749,  0.92097795757081369, -0.02651441186298376,  0.99200950450550751, -0.12334556714490325,
             0.72510924905369167, -0.14283988479115500,  0.67365669613657175, -0.57674111982491216,  0.40856335403608113,  0.70742184475876546,  0.37627949404651057,  0.90148363992456226, -0.21387143172584849,
            -0.57674111982491194,  0.40856335403608113,  0.70742184475876546,  0.72510924905369167, -0.14283988479115500,  0.67365669613657175,  0.37627949404651045,  0.90148363992456226, -0.21387143172584849,
            -0.37534416687016486,  0.92097795757081413,  0.10448137664308732,  0.92650620200843892,  0.36957471874357795,  0.07071622802089356, -0.02651441186298370, -0.12334556714490325,  0.99200950450550751,
             0.92650620200843870,  0.36957471874357795,  0.07071622802089356, -0.37534416687016486,  0.92097795757081413,  0.10448137664308732, -0.02651441186298359, -0.12334556714490325,  0.99200950450550751,
             0.72510924905369145,  0.67365669613657153, -0.14283988479115467, -0.57674111982491194,  0.70742184475876568,  0.40856335403608079,  0.37627949404651079, -0.21387143172584872,  0.90148363992456226,
            -0.57674111982491194,  0.70742184475876568,  0.40856335403608085,  0.72510924905369167,  0.67365669613657153, -0.14283988479115467,  0.37627949404651062, -0.21387143172584849,  0.90148363992456226,
            -0.02651441186298362, -0.12334556714490336,  0.99200950450550751, -0.37534416687016497,  0.92097795757081391,  0.10448137664308749,  0.92650620200843870,  0.36957471874357795,  0.07071622802089361,
            -0.37534416687016497,  0.92097795757081391,  0.10448137664308749, -0.02651441186298356, -0.12334556714490336,  0.99200950450550751,  0.92650620200843870,  0.36957471874357795,  0.07071622802089361,
            -0.57674111982491194,  0.70742184475876568,  0.40856335403608091,  0.37627949404651051, -0.21387143172584849,  0.90148363992456249,  0.72510924905369167,  0.67365669613657153, -0.14283988479115478,
             0.37627949404651051, -0.21387143172584849,  0.90148363992456249, -0.57674111982491216,  0.70742184475876568,  0.40856335403608079,  0.72510924905369190,  0.67365669613657131, -0.14283988479115467,
            -0.02651441186298384,  0.99200950450550751, -0.12334556714490325, -0.37534416687016486,  0.10448137664308754,  0.92097795757081391,  0.92650620200843892,  0.07071622802089356,  0.36957471874357806,
            -0.37534416687016475,  0.10448137664308771,  0.92097795757081369, -0.02651441186298370,  0.99200950450550751, -0.12334556714490325,  0.92650620200843870,  0.07071622802089345,  0.36957471874357822,
            -0.57674111982491194,  0.40856335403608102,  0.70742184475876546,  0.37627949404651051,  0.90148363992456249, -0.21387143172584860,  0.72510924905369167, -0.14283988479115467,  0.67365669613657153,
             0.37627949404651040,  0.90148363992456226, -0.21387143172584849, -0.57674111982491194,  0.40856335403608102,  0.70742184475876546,  0.72510924905369167, -0.14283988479115467,  0.67365669613657153,
            -0.78867513459481253,  0.61476523277246187,  0.00743240513391599,  0.57735026918962573,  0.74472292782631933, -0.33474531727338774,  0.21132486540518666,  0.25971420705077969,  0.94227932796886937,
             0.57735026918962562,  0.74472292782631933, -0.33474531727338763, -0.78867513459481253,  0.61476523277246187,  0.00743240513391608,  0.21132486540518683,  0.25971420705077974,  0.94227932796886926,
             0.21132486540518713,  0.81975403804892744, -0.53230171741593757, -0.78867513459481264,  0.46470301232724576,  0.40254520541901562,  0.57735026918962551,  0.33474531727338763,  0.74472292782631966,
            -0.78867513459481264,  0.46470301232724587,  0.40254520541901562,  0.21132486540518713,  0.81975403804892744, -0.53230171741593768,  0.57735026918962540,  0.33474531727338763,  0.74472292782631966,
             0.78867513459481264,  0.61476523277246142,  0.00743240513391624, -0.57735026918962573,  0.74472292782631966, -0.33474531727338785, -0.21132486540518700,  0.25971420705077991,  0.94227932796886937,
            -0.57735026918962562,  0.74472292782631966, -0.33474531727338785,  0.78867513459481275,  0.61476523277246142,  0.00743240513391630, -0.21132486540518711,  0.25971420705077986,  0.94227932796886926,
            -0.21132486540518700,  0.81975403804892766, -0.53230171741593757,  0.78867513459481275,  0.46470301232724542,  0.40254520541901573, -0.57735026918962573,  0.33474531727338785,  0.74472292782631944,
             0.78867513459481264,  0.46470301232724548,  0.40254520541901567, -0.21132486540518705,  0.81975403804892766, -0.53230171741593757, -0.57735026918962562,  0.33474531727338785,  0.74472292782631966,
            -0.57735026918962573,  0.33474531727338808,  0.74472292782631933, -0.21132486540518711,  0.81975403804892766, -0.53230171741593768,  0.78867513459481264,  0.46470301232724531,  0.40254520541901595,
            -0.21132486540518705,  0.81975403804892766, -0.53230171741593757, -0.57735026918962595,  0.33474531727338791,  0.74472292782631944,  0.78867513459481275,  0.46470301232724542,  0.40254520541901573,
            -0.57735026918962573,  0.74472292782631955, -0.33474531727338763, -0.21132486540518713,  0.25971420705077980,  0.94227932796886937,  0.78867513459481275,  0.61476523277246153,  0.00743240513391605,
            -0.21132486540518722,  0.25971420705077980,  0.94227932796886926, -0.57735026918962562,  0.74472292782631955, -0.33474531727338763,  0.78867513459481275,  0.61476523277246153,  0.00743240513391616,
             0.57735026918962573,  0.33474531727338774,  0.74472292782631944,  0.21132486540518708,  0.81975403804892744, -0.53230171741593757, -0.78867513459481275,  0.46470301232724576,  0.40254520541901584,
             0.21132486540518716,  0.81975403804892744, -0.53230171741593746,  0.57735026918962562,  0.33474531727338774,  0.74472292782631944, -0.78867513459481275,  0.46470301232724576,  0.40254520541901573,
             0.57735026918962573,  0.74472292782631921, -0.33474531727338752,  0.21132486540518705,  0.25971420705077969,  0.94227932796886948, -0.78867513459481264,  0.61476523277246198,  0.00743240513391588,
             0.21132486540518700,  0.25971420705077969,  0.94227932796886937,  0.57735026918962595,  0.74472292782631944, -0.33474531727338763, -0.78867513459481275,  0.61476523277246187,  0.00743240513391610,
            -0.78867513459481264,  0.46470301232724587,  0.40254520541901562,  0.57735026918962573,  0.33474531727338763,  0.74472292782631966,  0.21132486540518694,  0.81975403804892744, -0.53230171741593757,
             0.57735026918962562,  0.33474531727338763,  0.74472292782631966, -0.78867513459481275,  0.46470301232724587,  0.40254520541901556,  0.21132486540518705,  0.81975403804892744, -0.53230171741593746,
             0.21132486540518694,  0.25971420705077969,  0.94227932796886937, -0.78867513459481275,  0.61476523277246187,  0.00743240513391613,  0.57735026918962573,  0.74472292782631944, -0.33474531727338763,
            -0.78867513459481264,  0.61476523277246176,  0.00743240513391619,  0.21132486540518700,  0.25971420705077969,  0.94227932796886937,  0.57735026918962562,  0.74472292782631944, -0.33474531727338785,
             0.78867513459481253,  0.46470301232724542,  0.40254520541901584, -0.57735026918962573,  0.33474531727338797,  0.74472292782631955, -0.21132486540518672,  0.81975403804892766, -0.53230171741593757,
            -0.57735026918962562,  0.33474531727338797,  0.74472292782631944,  0.78867513459481253,  0.46470301232724542,  0.40254520541901578, -0.21132486540518688,  0.81975403804892755, -0.53230171741593746,
            -0.21132486540518719,  0.25971420705077980,  0.94227932796886937,  0.78867513459481264,  0.61476523277246153,  0.00743240513391624, -0.57735026918962551,  0.74472292782631966, -0.33474531727338785,
             0.78867513459481264,  0.61476523277246142,  0.00743240513391627, -0.21132486540518719,  0.25971420705077980,  0.94227932796886948, -0.57735026918962540,  0.74472292782631966, -0.33474531727338785,
             0.78867513459481253,  0.00743240513391608,  0.61476523277246153, -0.57735026918962562, -0.33474531727338763,  0.74472292782631966, -0.21132486540518694,  0.94227932796886926,  0.25971420705077980,
            -0.57735026918962551, -0.33474531727338763,  0.74472292782631966,  0.78867513459481253,  0.00743240513391613,  0.61476523277246153, -0.21132486540518711,  0.94227932796886926,  0.25971420705077980,
            -0.21132486540518688, -0.53230171741593746,  0.81975403804892755,  0.78867513459481264,  0.40254520541901562,  0.46470301232724548, -0.57735026918962573,  0.74472292782631955,  0.33474531727338785,
             0.78867513459481264,  0.40254520541901562,  0.46470301232724553, -0.21132486540518697, -0.53230171741593757,  0.81975403804892755, -0.57735026918962573,  0.74472292782631955,  0.33474531727338785,
            -0.78867513459481264,  0.00743240513391619,  0.61476523277246176,  0.57735026918962584, -0.33474531727338763,  0.74472292782631944,  0.21132486540518683,  0.94227932796886926,  0.25971420705077974,
             0.57735026918962584, -0.33474531727338763,  0.74472292782631944, -0.78867513459481264,  0.00743240513391624,  0.61476523277246176,  0.21132486540518683,  0.94227932796886915,  0.25971420705077974,
             0.21132486540518722, -0.53230171741593746,  0.81975403804892744, -0.78867513459481275,  0.40254520541901562,  0.46470301232724570,  0.57735026918962551,  0.74472292782631955,  0.33474531727338774,
            -0.78867513459481275,  0.40254520541901556,  0.46470301232724576,  0.21132486540518722, -0.53230171741593746,  0.81975403804892744,  0.57735026918962551,  0.74472292782631955,  0.33474531727338769,
             0.57735026918962551,  0.74472292782631933,  0.33474531727338774,  0.21132486540518733, -0.53230171741593746,  0.81975403804892755, -0.78867513459481287,  0.40254520541901584,  0.46470301232724565,
             0.21132486540518716, -0.53230171741593746,  0.81975403804892744,  0.57735026918962551,  0.74472292782631944,  0.33474531727338780, -0.78867513459481264,  0.40254520541901573,  0.46470301232724565,
             0.57735026918962584, -0.33474531727338763,  0.74472292782631944,  0.21132486540518686,  0.94227932796886937,  0.25971420705077969, -0.78867513459481264,  0.00743240513391610,  0.61476523277246187,
             0.21132486540518691,  0.94227932796886926,  0.25971420705077969,  0.57735026918962584, -0.33474531727338763,  0.74472292782631944, -0.78867513459481264,  0.00743240513391616,  0.61476523277246187,
            -0.57735026918962584,  0.74472292782631944,  0.33474531727338791, -0.21132486540518691, -0.53230171741593757,  0.81975403804892766,  0.78867513459481264,  0.40254520541901573,  0.46470301232724542,
            -0.21132486540518697, -0.53230171741593746,  0.81975403804892766, -0.57735026918962584,  0.74472292782631944,  0.33474531727338785,  0.78867513459481264,  0.40254520541901573,  0.46470301232724548,
            -0.57735026918962551, -0.33474531727338752,  0.74472292782631955, -0.21132486540518738,  0.94227932796886926,  0.25971420705077974,  0.78867513459481287,  0.00743240513391605,  0.61476523277246165,
            -0.21132486540518722,  0.94227932796886926,  0.25971420705077986, -0.57735026918962551, -0.33474531727338763,  0.74472292782631944,  0.78867513459481264,  0.00743240513391610,  0.61476523277246165,
             0.78867513459481264,  0.40254520541901567,  0.46470301232724553, -0.57735026918962584,  0.74472292782631944,  0.33474531727338785, -0.21132486540518688, -0.53230171741593746,  0.81975403804892755,
            -0.57735026918962584,  0.74472292782631944,  0.33474531727338785,  0.78867513459481264,  0.40254520541901562,  0.46470301232724553, -0.21132486540518688, -0.53230171741593735,  0.81975403804892755,
            -0.21132486540518727,  0.94227932796886926,  0.25971420705077980,  0.78867513459481275,  0.00743240513391627,  0.61476523277246153, -0.57735026918962551, -0.33474531727338774,  0.74472292782631955,
             0.78867513459481275,  0.00743240513391630,  0.61476523277246153, -0.21132486540518727,  0.94227932796886926,  0.25971420705077980, -0.57735026918962551, -0.33474531727338774,  0.74472292782631966,
            -0.78867513459481253,  0.40254520541901578,  0.46470301232724570,  0.57735026918962562,  0.74472292782631944,  0.33474531727338769,  0.21132486540518688, -0.53230171741593746,  0.81975403804892744,
             0.57735026918962551,  0.74472292782631944,  0.33474531727338769, -0.78867513459481253,  0.40254520541901573,  0.46470301232724576,  0.21132486540518705, -0.53230171741593746,  0.81975403804892744,
             0.21132486540518683,  0.94227932796886926,  0.25971420705077974, -0.78867513459481264,  0.00743240513391621,  0.61476523277246176,  0.57735026918962573, -0.33474531727338774,  0.74472292782631944,
            -0.78867513459481264,  0.00743240513391624,  0.61476523277246176,  0.21132486540518691,  0.94227932796886937,  0.25971420705077974,  0.57735026918962573, -0.33474531727338774,  0.74472292782631944,
            -0.37627949404651051,  0.90148363992456215, -0.21387143172584822, -0.72510924905369167, -0.14283988479115428,  0.67365669613657153,  0.57674111982491172,  0.40856335403608096,  0.70742184475876546,
            -0.72510924905369156, -0.14283988479115417,  0.67365669613657142, -0.37627949404651051,  0.90148363992456215, -0.21387143172584822,  0.57674111982491161,  0.40856335403608074,  0.70742184475876546,
            -0.92650620200843881,  0.07071622802089400,  0.36957471874357806,  0.02651441186298364,  0.99200950450550718, -0.12334556714490297,  0.37534416687016459,  0.10448137664308754,  0.92097795757081369,
             0.02651441186298359,  0.99200950450550718, -0.12334556714490297, -0.92650620200843881,  0.07071622802089395,  0.36957471874357806,  0.37534416687016481,  0.10448137664308765,  0.92097795757081369,
            -0.37627949404651051, -0.21387143172584799,  0.90148363992456215, -0.72510924905369167,  0.67365669613657164, -0.14283988479115450,  0.57674111982491172,  0.70742184475876513,  0.40856335403608113,
            -0.72510924905369167,  0.67365669613657153, -0.14283988479115450, -0.37627949404651051, -0.21387143172584799,  0.90148363992456215,  0.57674111982491172,  0.70742184475876513,  0.40856335403608113,
            -0.92650620200843870,  0.36957471874357845,  0.07071622802089356,  0.02651441186298351, -0.12334556714490308,  0.99200950450550740,  0.37534416687016470,  0.92097795757081335,  0.10448137664308776,
             0.02651441186298342, -0.12334556714490308,  0.99200950450550740, -0.92650620200843870,  0.36957471874357850,  0.07071622802089350,  0.37534416687016470,  0.92097795757081335,  0.10448137664308782,
             0.37534416687016448,  0.92097795757081335,  0.10448137664308771, -0.92650620200843881,  0.36957471874357839,  0.07071622802089356,  0.02651441186298384, -0.12334556714490297,  0.99200950450550740,
            -0.92650620200843881,  0.36957471874357839,  0.07071622802089356,  0.37534416687016459,  0.92097795757081358,  0.10448137664308765,  0.02651441186298367, -0.12334556714490319,  0.99200950450550751,
            -0.72510924905369145,  0.67365669613657175, -0.14283988479115461,  0.57674111982491183,  0.70742184475876524,  0.40856335403608113, -0.37627949404651073, -0.21387143172584810,  0.90148363992456215,
             0.57674111982491161,  0.70742184475876524,  0.40856335403608107, -0.72510924905369145,  0.67365669613657175, -0.14283988479115439, -0.37627949404651062, -0.21387143172584810,  0.90148363992456215,
             0.37534416687016448,  0.10448137664308743,  0.92097795757081380, -0.92650620200843881,  0.07071622802089395,  0.36957471874357806,  0.02651441186298384,  0.99200950450550729, -0.12334556714490297,
            -0.92650620200843858,  0.07071622802089395,  0.36957471874357811,  0.37534416687016448,  0.10448137664308743,  0.92097795757081358,  0.02651441186298367,  0.99200950450550729, -0.12334556714490297,
            -0.72510924905369145, -0.14283988479115417,  0.67365669613657153,  0.57674111982491183,  0.40856335403608079,  0.70742184475876568, -0.37627949404651073,  0.90148363992456215, -0.21387143172584822,
             0.57674111982491183,  0.40856335403608079,  0.70742184475876568, -0.72510924905369156, -0.14283988479115439,  0.67365669613657153, -0.37627949404651062,  0.90148363992456237, -0.21387143172584833,
             0.02651441186298356,  0.99200950450550718, -0.12334556714490297,  0.37534416687016470,  0.10448137664308754,  0.92097795757081369, -0.92650620200843870,  0.07071622802089406,  0.36957471874357806,
             0.37534416687016470,  0.10448137664308760,  0.92097795757081369,  0.02651441186298356,  0.99200950450550718, -0.12334556714490297, -0.92650620200843870,  0.07071622802089406,  0.36957471874357806,
             0.57674111982491172,  0.40856335403608074,  0.70742184475876568, -0.37627949404651040,  0.90148363992456226, -0.21387143172584822, -0.72510924905369167, -0.14283988479115417,  0.67365669613657142,
            -0.37627949404651034,  0.90148363992456226, -0.21387143172584822,  0.57674111982491172,  0.40856335403608068,  0.70742184475876568, -0.72510924905369167, -0.14283988479115417,  0.67365669613657131,
             0.02651441186298362, -0.12334556714490297,  0.99200950450550740,  0.37534416687016470,  0.92097795757081347,  0.10448137664308771, -0.92650620200843870,  0.36957471874357822,  0.07071622802089372,
             0.37534416687016459,  0.92097795757081335,  0.10448137664308776,  0.02651441186298359, -0.12334556714490297,  0.99200950450550740, -0.92650620200843858,  0.36957471874357845,  0.07071622802089372,
             0.57674111982491183,  0.70742184475876524,  0.40856335403608113, -0.37627949404651057, -0.21387143172584799,  0.90148363992456215, -0.72510924905369156,  0.67365669613657164, -0.14283988479115450,
            -0.37627949404651051, -0.21387143172584799,  0.90148363992456215,  0.57674111982491183,  0.70742184475876524,  0.40856335403608113, -0.72510924905369178,  0.67365669613657153, -0.14283988479115450,
            -0.21387143172584849,  0.37627949404651095,  0.90148363992456249,  0.67365669613657153,  0.72510924905369167, -0.14283988479115489,  0.70742184475876557, -0.57674111982491216,  0.40856335403608096,
             0.67365669613657131,  0.72510924905369145, -0.14283988479115467, -0.21387143172584849,  0.37627949404651084,  0.90148363992456249,  0.70742184475876568, -0.57674111982491194,  0.40856335403608079,
             0.36957471874357795,  0.92650620200843892,  0.07071622802089361, -0.12334556714490325, -0.02651441186298345,  0.99200950450550773,  0.92097795757081369, -0.37534416687016509,  0.10448137664308732,
            -0.12334556714490325, -0.02651441186298342,  0.99200950450550773,  0.36957471874357795,  0.92650620200843892,  0.07071622802089361,  0.92097795757081369, -0.37534416687016520,  0.10448137664308749,
             0.90148363992456226,  0.37627949404651029, -0.21387143172584849, -0.14283988479115489,  0.72510924905369190,  0.67365669613657175,  0.40856335403608118, -0.57674111982491194,  0.70742184475876546,
            -0.14283988479115489,  0.72510924905369190,  0.67365669613657175,  0.90148363992456226,  0.37627949404651018, -0.21387143172584849,  0.40856335403608118, -0.57674111982491194,  0.70742184475876546,
             0.07071622802089339,  0.92650620200843892,  0.36957471874357828,  0.99200950450550751, -0.02651441186298381, -0.12334556714490336,  0.10448137664308771, -0.37534416687016486,  0.92097795757081369,
             0.99200950450550751, -0.02651441186298381, -0.12334556714490347,  0.07071622802089339,  0.92650620200843892,  0.36957471874357828,  0.10448137664308771, -0.37534416687016486,  0.92097795757081369,
             0.10448137664308765, -0.37534416687016453,  0.92097795757081369,  0.07071622802089333,  0.92650620200843914,  0.36957471874357828,  0.99200950450550751, -0.02651441186298420, -0.12334556714490336,
             0.07071622802089333,  0.92650620200843892,  0.36957471874357833,  0.10448137664308754, -0.37534416687016464,  0.92097795757081391,  0.99200950450550773, -0.02651441186298392, -0.12334556714490347,
            -0.14283988479115500,  0.72510924905369190,  0.67365669613657186,  0.40856335403608118, -0.57674111982491194,  0.70742184475876546,  0.90148363992456226,  0.37627949404651040, -0.21387143172584860,
             0.40856335403608113, -0.57674111982491194,  0.70742184475876546, -0.14283988479115489,  0.72510924905369190,  0.67365669613657186,  0.90148363992456226,  0.37627949404651029, -0.21387143172584860,
             0.92097795757081413, -0.37534416687016509,  0.10448137664308721,  0.36957471874357789,  0.92650620200843892,  0.07071622802089361, -0.12334556714490325, -0.02651441186298353,  0.99200950450550773,
             0.36957471874357795,  0.92650620200843870,  0.07071622802089361,  0.92097795757081391, -0.37534416687016509,  0.10448137664308721, -0.12334556714490313, -0.02651441186298342,  0.99200950450550773,
             0.67365669613657142,  0.72510924905369145, -0.14283988479115467,  0.70742184475876568, -0.57674111982491216,  0.40856335403608079, -0.21387143172584849,  0.37627949404651106,  0.90148363992456249,
             0.70742184475876568, -0.57674111982491216,  0.40856335403608074,  0.67365669613657153,  0.72510924905369145, -0.14283988479115489, -0.21387143172584872,  0.37627949404651079,  0.90148363992456249,
            -0.12334556714490325, -0.02651441186298339,  0.99200950450550751,  0.92097795757081391, -0.37534416687016509,  0.10448137664308738,  0.36957471874357789,  0.92650620200843870,  0.07071622802089367,
             0.92097795757081391, -0.37534416687016509,  0.10448137664308738, -0.12334556714490325, -0.02651441186298334,  0.99200950450550751,  0.36957471874357789,  0.92650620200843870,  0.07071622802089367,
             0.70742184475876568, -0.57674111982491216,  0.40856335403608079, -0.21387143172584849,  0.37627949404651068,  0.90148363992456249,  0.67365669613657131,  0.72510924905369167, -0.14283988479115467,
            -0.21387143172584849,  0.37627949404651068,  0.90148363992456249,  0.70742184475876568, -0.57674111982491216,  0.40856335403608079,  0.67365669613657131,  0.72510924905369167, -0.14283988479115467,
             0.99200950450550751, -0.02651441186298412, -0.12334556714490336,  0.10448137664308760, -0.37534416687016475,  0.92097795757081391,  0.07071622802089350,  0.92650620200843914,  0.36957471874357811,
             0.10448137664308771, -0.37534416687016453,  0.92097795757081369,  0.99200950450550751, -0.02651441186298401, -0.12334556714490347,  0.07071622802089345,  0.92650620200843870,  0.36957471874357828,
             0.40856335403608113, -0.57674111982491194,  0.70742184475876546,  0.90148363992456226,  0.37627949404651029, -0.21387143172584860, -0.14283988479115467,  0.72510924905369190,  0.67365669613657175,
             0.90148363992456226,  0.37627949404651029, -0.21387143172584860,  0.40856335403608113, -0.57674111982491216,  0.70742184475876546, -0.14283988479115467,  0.72510924905369212,  0.67365669613657153,
             0.61476523277246187, -0.78867513459481253,  0.00743240513391605,  0.74472292782631944,  0.57735026918962573, -0.33474531727338774,  0.25971420705077969,  0.21132486540518672,  0.94227932796886926,
             0.74472292782631944,  0.57735026918962551, -0.33474531727338763,  0.61476523277246176, -0.78867513459481242,  0.00743240513391608,  0.25971420705077974,  0.21132486540518683,  0.94227932796886926,
             0.81975403804892766,  0.21132486540518711, -0.53230171741593757,  0.46470301232724565, -0.78867513459481253,  0.40254520541901573,  0.33474531727338769,  0.57735026918962540,  0.74472292782631955,
             0.46470301232724576, -0.78867513459481242,  0.40254520541901573,  0.81975403804892766,  0.21132486540518711, -0.53230171741593768,  0.33474531727338763,  0.57735026918962529,  0.74472292782631966,
             0.61476523277246153,  0.78867513459481242,  0.00743240513391616,  0.74472292782631966, -0.57735026918962551, -0.33474531727338774,  0.25971420705077980, -0.21132486540518694,  0.94227932796886926,
             0.74472292782631966, -0.57735026918962551, -0.33474531727338774,  0.61476523277246153,  0.78867513459481242,  0.00743240513391621,  0.25971420705077980, -0.21132486540518700,  0.94227932796886926,
             0.81975403804892755, -0.21132486540518694, -0.53230171741593757,  0.46470301232724553,  0.78867513459481264,  0.40254520541901562,  0.33474531727338785, -0.57735026918962562,  0.74472292782631955,
             0.46470301232724553,  0.78867513459481264,  0.40254520541901562,  0.81975403804892766, -0.21132486540518702, -0.53230171741593757,  0.33474531727338774, -0.57735026918962551,  0.74472292782631966,
             0.33474531727338791, -0.57735026918962551,  0.74472292782631933,  0.81975403804892766, -0.21132486540518708, -0.53230171741593757,  0.46470301232724542,  0.78867513459481242,  0.40254520541901595,
             0.81975403804892755, -0.21132486540518700, -0.53230171741593757,  0.33474531727338785, -0.57735026918962573,  0.74472292782631955,  0.46470301232724553,  0.78867513459481264,  0.40254520541901573,
             0.74472292782631944, -0.57735026918962562, -0.33474531727338763,  0.25971420705077980, -0.21132486540518711,  0.94227932796886937,  0.61476523277246165,  0.78867513459481264,  0.00743240513391599,
             0.25971420705077974, -0.21132486540518711,  0.94227932796886926,  0.74472292782631955, -0.57735026918962551, -0.33474531727338763,  0.61476523277246165,  0.78867513459481264,  0.00743240513391616,
             0.33474531727338780,  0.57735026918962562,  0.74472292782631944,  0.81975403804892744,  0.21132486540518705, -0.53230171741593757,  0.46470301232724565, -0.78867513459481264,  0.40254520541901584,
             0.81975403804892755,  0.21132486540518705, -0.53230171741593746,  0.33474531727338774,  0.57735026918962551,  0.74472292782631944,  0.46470301232724565, -0.78867513459481264,  0.40254520541901573,
             0.74472292782631944,  0.57735026918962551, -0.33474531727338752,  0.25971420705077969,  0.21132486540518702,  0.94227932796886937,  0.61476523277246187, -0.78867513459481242,  0.00743240513391594,
             0.25971420705077974,  0.21132486540518694,  0.94227932796886937,  0.74472292782631944,  0.57735026918962573, -0.33474531727338774,  0.61476523277246176, -0.78867513459481264,  0.00743240513391610,
             0.46470301232724570, -0.78867513459481242,  0.40254520541901573,  0.33474531727338769,  0.57735026918962551,  0.74472292782631955,  0.81975403804892744,  0.21132486540518688, -0.53230171741593746,
             0.33474531727338769,  0.57735026918962551,  0.74472292782631955,  0.46470301232724576, -0.78867513459481242,  0.40254520541901562,  0.81975403804892744,  0.21132486540518694, -0.53230171741593746,
             0.25971420705077974,  0.21132486540518688,  0.94227932796886937,  0.61476523277246176, -0.78867513459481264,  0.00743240513391624,  0.74472292782631944,  0.57735026918962562, -0.33474531727338774,
             0.61476523277246176, -0.78867513459481264,  0.00743240513391624,  0.25971420705077969,  0.21132486540518697,  0.94227932796886937,  0.74472292782631955,  0.57735026918962551, -0.33474531727338785,
             0.46470301232724542,  0.78867513459481253,  0.40254520541901584,  0.33474531727338785, -0.57735026918962573,  0.74472292782631955,  0.81975403804892766, -0.21132486540518677, -0.53230171741593746,
             0.33474531727338785, -0.57735026918962551,  0.74472292782631944,  0.46470301232724553,  0.78867513459481242,  0.40254520541901578,  0.81975403804892755, -0.21132486540518688, -0.53230171741593746,
             0.25971420705077963, -0.21132486540518716,  0.94227932796886937,  0.61476523277246165,  0.78867513459481253,  0.00743240513391613,  0.74472292782631966, -0.57735026918962540, -0.33474531727338774,
             0.61476523277246153,  0.78867513459481242,  0.00743240513391616,  0.25971420705077969, -0.21132486540518716,  0.94227932796886948,  0.74472292782631966, -0.57735026918962529, -0.33474531727338785,
             0.00743240513391599,  0.78867513459481264,  0.61476523277246165, -0.33474531727338774, -0.57735026918962573,  0.74472292782631966,  0.94227932796886937, -0.21132486540518700,  0.25971420705077963,
            -0.33474531727338763, -0.57735026918962551,  0.74472292782631966,  0.00743240513391608,  0.78867513459481264,  0.61476523277246176,  0.94227932796886926, -0.21132486540518711,  0.25971420705077958,
            -0.53230171741593757, -0.21132486540518702,  0.81975403804892766,  0.40254520541901573,  0.78867513459481275,  0.46470301232724553,  0.74472292782631966, -0.57735026918962562,  0.33474531727338774,
             0.40254520541901562,  0.78867513459481264,  0.46470301232724559, -0.53230171741593757, -0.21132486540518705,  0.81975403804892766,  0.74472292782631966, -0.57735026918962562,  0.33474531727338769,
             0.00743240513391610, -0.78867513459481275,  0.61476523277246176, -0.33474531727338774,  0.57735026918962573,  0.74472292782631955,  0.94227932796886937,  0.21132486540518694,  0.25971420705077963,
            -0.33474531727338774,  0.57735026918962584,  0.74472292782631955,  0.00743240513391621, -0.78867513459481287,  0.61476523277246176,  0.94227932796886926,  0.21132486540518700,  0.25971420705077963,
            -0.53230171741593757,  0.21132486540518708,  0.81975403804892766,  0.40254520541901573, -0.78867513459481264,  0.46470301232724565,  0.74472292782631966,  0.57735026918962551,  0.33474531727338769,
             0.40254520541901562, -0.78867513459481264,  0.46470301232724565, -0.53230171741593757,  0.21132486540518716,  0.81975403804892766,  0.74472292782631966,  0.57735026918962551,  0.33474531727338769,
             0.74472292782631944,  0.57735026918962562,  0.33474531727338774, -0.53230171741593757,  0.21132486540518716,  0.81975403804892766,  0.40254520541901595, -0.78867513459481275,  0.46470301232724553,
            -0.53230171741593757,  0.21132486540518722,  0.81975403804892766,  0.74472292782631955,  0.57735026918962584,  0.33474531727338763,  0.40254520541901573, -0.78867513459481287,  0.46470301232724559,
            -0.33474531727338774,  0.57735026918962595,  0.74472292782631966,  0.94227932796886948,  0.21132486540518697,  0.25971420705077952,  0.00743240513391599, -0.78867513459481287,  0.61476523277246176,
             0.94227932796886937,  0.21132486540518697,  0.25971420705077952, -0.33474531727338774,  0.57735026918962595,  0.74472292782631966,  0.00743240513391610, -0.78867513459481287,  0.61476523277246176,
             0.74472292782631955, -0.57735026918962595,  0.33474531727338769, -0.53230171741593768, -0.21132486540518702,  0.81975403804892777,  0.40254520541901584,  0.78867513459481287,  0.46470301232724553,
            -0.53230171741593757, -0.21132486540518702,  0.81975403804892777,  0.74472292782631955, -0.57735026918962595,  0.33474531727338769,  0.40254520541901573,  0.78867513459481287,  0.46470301232724553,
            -0.33474531727338763, -0.57735026918962562,  0.74472292782631955,  0.94227932796886937, -0.21132486540518722,  0.25971420705077963,  0.00743240513391594,  0.78867513459481275,  0.61476523277246176,
             0.94227932796886937, -0.21132486540518727,  0.25971420705077958, -0.33474531727338774, -0.57735026918962584,  0.74472292782631966,  0.00743240513391610,  0.78867513459481287,  0.61476523277246165,
             0.40254520541901573,  0.78867513459481275,  0.46470301232724553,  0.74472292782631955, -0.57735026918962573,  0.33474531727338774, -0.53230171741593757, -0.21132486540518700,  0.81975403804892766,
             0.74472292782631955, -0.57735026918962584,  0.33474531727338774,  0.40254520541901562,  0.78867513459481287,  0.46470301232724553, -0.53230171741593746, -0.21132486540518705,  0.81975403804892766,
             0.94227932796886937, -0.21132486540518713,  0.25971420705077963,  0.00743240513391616,  0.78867513459481264,  0.61476523277246165, -0.33474531727338785, -0.57735026918962551,  0.74472292782631966,
             0.00743240513391621,  0.78867513459481264,  0.61476523277246165,  0.94227932796886937, -0.21132486540518722,  0.25971420705077963, -0.33474531727338785, -0.57735026918962551,  0.74472292782631966,
             0.40254520541901584, -0.78867513459481264,  0.46470301232724559,  0.74472292782631955,  0.57735026918962573,  0.33474531727338763, -0.53230171741593757,  0.21132486540518694,  0.81975403804892766,
             0.74472292782631944,  0.57735026918962551,  0.33474531727338763,  0.40254520541901578, -0.78867513459481264,  0.46470301232724559, -0.53230171741593746,  0.21132486540518705,  0.81975403804892766,
             0.94227932796886937,  0.21132486540518697,  0.25971420705077963,  0.00743240513391616, -0.78867513459481275,  0.61476523277246176, -0.33474531727338785,  0.57735026918962562,  0.74472292782631955,
             0.00743240513391624, -0.78867513459481264,  0.61476523277246176,  0.94227932796886937,  0.21132486540518700,  0.25971420705077958, -0.33474531727338785,  0.57735026918962562,  0.74472292782631966,
             0.90148363992456226, -0.37627949404651062, -0.21387143172584822, -0.14283988479115450, -0.72510924905369178,  0.67365669613657164,  0.40856335403608096,  0.57674111982491183,  0.70742184475876535,
            -0.14283988479115428, -0.72510924905369156,  0.67365669613657153,  0.90148363992456226, -0.37627949404651062, -0.21387143172584822,  0.40856335403608079,  0.57674111982491161,  0.70742184475876546,
             0.07071622802089383, -0.92650620200843881,  0.36957471874357828,  0.99200950450550740,  0.02651441186298367, -0.12334556714490319,  0.10448137664308754,  0.37534416687016470,  0.92097795757081369,
             0.99200950450550729,  0.02651441186298362, -0.12334556714490319,  0.07071622802089383, -0.92650620200843881,  0.36957471874357828,  0.10448137664308760,  0.37534416687016470,  0.92097795757081358,
            -0.21387143172584822, -0.37627949404651062,  0.90148363992456237,  0.67365669613657164, -0.72510924905369178, -0.14283988479115450,  0.70742184475876524,  0.57674111982491183,  0.40856335403608091,
             0.67365669613657164, -0.72510924905369178, -0.14283988479115450, -0.21387143172584822, -0.37627949404651057,  0.90148363992456237,  0.70742184475876546,  0.57674111982491183,  0.40856335403608091,
             0.36957471874357839, -0.92650620200843881,  0.07071622802089372, -0.12334556714490330,  0.02651441186298351,  0.99200950450550762,  0.92097795757081358,  0.37534416687016481,  0.10448137664308754,
            -0.12334556714490319,  0.02651441186298345,  0.99200950450550751,  0.36957471874357839, -0.92650620200843870,  0.07071622802089361,  0.92097795757081358,  0.37534416687016481,  0.10448137664308765,
             0.92097795757081358,  0.37534416687016459,  0.10448137664308754,  0.36957471874357839, -0.92650620200843881,  0.07071622802089372, -0.12334556714490308,  0.02651441186298378,  0.99200950450550751,
             0.36957471874357833, -0.92650620200843870,  0.07071622802089383,  0.92097795757081358,  0.37534416687016470,  0.10448137664308749, -0.12334556714490319,  0.02651441186298351,  0.99200950450550740,
             0.67365669613657175, -0.72510924905369156, -0.14283988479115461,  0.70742184475876524,  0.57674111982491183,  0.40856335403608102, -0.21387143172584822, -0.37627949404651073,  0.90148363992456237,
             0.70742184475876524,  0.57674111982491172,  0.40856335403608102,  0.67365669613657175, -0.72510924905369156, -0.14283988479115450, -0.21387143172584810, -0.37627949404651062,  0.90148363992456226,
             0.10448137664308743,  0.37534416687016459,  0.92097795757081380,  0.07071622802089389, -0.92650620200843881,  0.36957471874357817,  0.99200950450550740,  0.02651441186298378, -0.12334556714490319,
             0.07071622802089395, -0.92650620200843870,  0.36957471874357817,  0.10448137664308743,  0.37534416687016459,  0.92097795757081369,  0.99200950450550729,  0.02651441186298373, -0.12334556714490308,
            -0.14283988479115439, -0.72510924905369156,  0.67365669613657164,  0.40856335403608079,  0.57674111982491183,  0.70742184475876546,  0.90148363992456226, -0.37627949404651073, -0.21387143172584833,
             0.40856335403608085,  0.57674111982491172,  0.70742184475876535, -0.14283988479115439, -0.72510924905369167,  0.67365669613657175,  0.90148363992456237, -0.37627949404651040, -0.21387143172584822,
             0.99200950450550740,  0.02651441186298370, -0.12334556714490319,  0.10448137664308754,  0.37534416687016481,  0.92097795757081369,  0.07071622802089389, -0.92650620200843881,  0.36957471874357828,
             0.10448137664308754,  0.37534416687016481,  0.92097795757081369,  0.99200950450550740,  0.02651441186298364, -0.12334556714490319,  0.07071622802089378, -0.92650620200843881,  0.36957471874357828,
             0.40856335403608079,  0.57674111982491183,  0.70742184475876546,  0.90148363992456249, -0.37627949404651040, -0.21387143172584844, -0.14283988479115439, -0.72510924905369178,  0.67365669613657164,
             0.90148363992456237, -0.37627949404651040, -0.21387143172584833,  0.40856335403608079,  0.57674111982491172,  0.70742184475876557, -0.14283988479115439, -0.72510924905369178,  0.67365669613657153,
            -0.12334556714490308,  0.02651441186298370,  0.99200950450550740,  0.92097795757081369,  0.37534416687016481,  0.10448137664308754,  0.36957471874357822, -0.92650620200843881,  0.07071622802089383,
             0.92097795757081347,  0.37534416687016459,  0.10448137664308765, -0.12334556714490308,  0.02651441186298370,  0.99200950450550740,  0.36957471874357839, -0.92650620200843858,  0.07071622802089372,
             0.70742184475876535,  0.57674111982491183,  0.40856335403608091, -0.21387143172584822, -0.37627949404651062,  0.90148363992456237,  0.67365669613657164, -0.72510924905369167, -0.14283988479115450,
            -0.21387143172584810, -0.37627949404651051,  0.90148363992456237,  0.70742184475876535,  0.57674111982491183,  0.40856335403608091,  0.67365669613657153, -0.72510924905369167, -0.14283988479115439,
            -0.21387143172584838,  0.90148363992456226,  0.37627949404651062,  0.67365669613657153, -0.14283988479115456,  0.72510924905369167,  0.70742184475876524,  0.40856335403608068, -0.57674111982491216,
             0.67365669613657142, -0.14283988479115434,  0.72510924905369145, -0.21387143172584849,  0.90148363992456204,  0.37627949404651051,  0.70742184475876535,  0.40856335403608041, -0.57674111982491194,
             0.36957471874357811,  0.07071622802089383,  0.92650620200843892, -0.12334556714490347,  0.99200950450550729, -0.02651441186298367,  0.92097795757081369,  0.10448137664308726, -0.37534416687016503,
            -0.12334556714490336,  0.99200950450550729, -0.02651441186298362,  0.36957471874357806,  0.07071622802089383,  0.92650620200843881,  0.92097795757081347,  0.10448137664308726, -0.37534416687016503,
             0.90148363992456226, -0.21387143172584849,  0.37627949404651040, -0.14283988479115478,  0.67365669613657175,  0.72510924905369190,  0.40856335403608079,  0.70742184475876502, -0.57674111982491205,
            -0.14283988479115478,  0.67365669613657164,  0.72510924905369190,  0.90148363992456249, -0.21387143172584838,  0.37627949404651040,  0.40856335403608068,  0.70742184475876502, -0.57674111982491194,
             0.07071622802089356,  0.36957471874357839,  0.92650620200843892,  0.99200950450550751, -0.12334556714490347, -0.02651441186298384,  0.10448137664308738,  0.92097795757081324, -0.37534416687016492,
             0.99200950450550751, -0.12334556714490336, -0.02651441186298373,  0.07071622802089345,  0.36957471874357845,  0.92650620200843892,  0.10448137664308743,  0.92097795757081324, -0.37534416687016503,
             0.10448137664308732,  0.92097795757081324, -0.37534416687016470,  0.07071622802089345,  0.36957471874357839,  0.92650620200843892,  0.99200950450550751, -0.12334556714490347, -0.02651441186298406,
             0.07071622802089345,  0.36957471874357845,  0.92650620200843892,  0.10448137664308721,  0.92097795757081324, -0.37534416687016481,  0.99200950450550751, -0.12334556714490347, -0.02651441186298389,
            -0.14283988479115478,  0.67365669613657186,  0.72510924905369167,  0.40856335403608079,  0.70742184475876502, -0.57674111982491216,  0.90148363992456226, -0.21387143172584838,  0.37627949404651051,
             0.40856335403608079,  0.70742184475876502, -0.57674111982491205, -0.14283988479115478,  0.67365669613657175,  0.72510924905369178,  0.90148363992456226, -0.21387143172584838,  0.37627949404651040,
             0.92097795757081369,  0.10448137664308699, -0.37534416687016492,  0.36957471874357806,  0.07071622802089389,  0.92650620200843892, -0.12334556714490336,  0.99200950450550729, -0.02651441186298384,
             0.36957471874357806,  0.07071622802089389,  0.92650620200843881,  0.92097795757081369,  0.10448137664308710, -0.37534416687016503, -0.12334556714490336,  0.99200950450550729, -0.02651441186298373,
             0.67365669613657153, -0.14283988479115445,  0.72510924905369145,  0.70742184475876546,  0.40856335403608046, -0.57674111982491216, -0.21387143172584860,  0.90148363992456226,  0.37627949404651073,
             0.70742184475876546,  0.40856335403608041, -0.57674111982491216,  0.67365669613657164, -0.14283988479115445,  0.72510924905369156, -0.21387143172584860,  0.90148363992456226,  0.37627949404651062,
            -0.12334556714490347,  0.99200950450550729, -0.02651441186298373,  0.92097795757081369,  0.10448137664308715, -0.37534416687016514,  0.36957471874357806,  0.07071622802089383,  0.92650620200843881,
             0.92097795757081369,  0.10448137664308721, -0.37534416687016514, -0.12334556714490358,  0.99200950450550729, -0.02651441186298367,  0.36957471874357817,  0.07071622802089389,  0.92650620200843870,
             0.70742184475876524,  0.40856335403608046, -0.57674111982491216, -0.21387143172584860,  0.90148363992456226,  0.37627949404651051,  0.67365669613657153, -0.14283988479115445,  0.72510924905369167,
            -0.21387143172584860,  0.90148363992456226,  0.37627949404651040,  0.70742184475876546,  0.40856335403608041, -0.57674111982491216,  0.67365669613657142, -0.14283988479115434,  0.72510924905369178,
             0.99200950450550729, -0.12334556714490336, -0.02651441186298392,  0.10448137664308726,  0.92097795757081347, -0.37534416687016492,  0.07071622802089367,  0.36957471874357817,  0.92650620200843892,
             0.10448137664308743,  0.92097795757081324, -0.37534416687016470,  0.99200950450550729, -0.12334556714490325, -0.02651441186298384,  0.07071622802089350,  0.36957471874357845,  0.92650620200843870,
             0.40856335403608074,  0.70742184475876502, -0.57674111982491216,  0.90148363992456226, -0.21387143172584838,  0.37627949404651040, -0.14283988479115478,  0.67365669613657153,  0.72510924905369178,
             0.90148363992456226, -0.21387143172584838,  0.37627949404651029,  0.40856335403608079,  0.70742184475876502, -0.57674111982491205, -0.14283988479115467,  0.67365669613657153,  0.72510924905369178,
             0.61476523277246187,  0.00743240513391605, -0.78867513459481253,  0.74472292782631944, -0.33474531727338758,  0.57735026918962573,  0.25971420705077974,  0.94227932796886926,  0.21132486540518666,
             0.74472292782631933, -0.33474531727338758,  0.57735026918962562,  0.61476523277246176,  0.00743240513391610, -0.78867513459481253,  0.25971420705077974,  0.94227932796886926,  0.21132486540518677,
             0.81975403804892744, -0.53230171741593746,  0.21132486540518716,  0.46470301232724570,  0.40254520541901562, -0.78867513459481264,  0.33474531727338774,  0.74472292782631966,  0.57735026918962540,
             0.46470301232724576,  0.40254520541901562, -0.78867513459481264,  0.81975403804892755, -0.53230171741593746,  0.21132486540518722,  0.33474531727338763,  0.74472292782631966,  0.57735026918962529,
             0.61476523277246153,  0.00743240513391624,  0.78867513459481242,  0.74472292782631955, -0.33474531727338769, -0.57735026918962551,  0.25971420705077986,  0.94227932796886926, -0.21132486540518694,
             0.74472292782631966, -0.33474531727338780, -0.57735026918962551,  0.61476523277246153,  0.00743240513391630,  0.78867513459481253,  0.25971420705077986,  0.94227932796886926, -0.21132486540518711,
             0.81975403804892755, -0.53230171741593746, -0.21132486540518688,  0.46470301232724548,  0.40254520541901567,  0.78867513459481264,  0.33474531727338785,  0.74472292782631966, -0.57735026918962573,
             0.46470301232724553,  0.40254520541901562,  0.78867513459481264,  0.81975403804892755, -0.53230171741593746, -0.21132486540518700,  0.33474531727338785,  0.74472292782631966, -0.57735026918962562,
             0.33474531727338797,  0.74472292782631921, -0.57735026918962562,  0.81975403804892766, -0.53230171741593746, -0.21132486540518705,  0.46470301232724542,  0.40254520541901601,  0.78867513459481264,
             0.81975403804892755, -0.53230171741593746, -0.21132486540518694,  0.33474531727338791,  0.74472292782631944, -0.57735026918962573,  0.46470301232724548,  0.40254520541901573,  0.78867513459481264,
             0.74472292782631944, -0.33474531727338758, -0.57735026918962551,  0.25971420705077980,  0.94227932796886948, -0.21132486540518716,  0.61476523277246165,  0.00743240513391602,  0.78867513459481264,
             0.25971420705077986,  0.94227932796886926, -0.21132486540518722,  0.74472292782631944, -0.33474531727338758, -0.57735026918962551,  0.61476523277246165,  0.00743240513391613,  0.78867513459481264,
             0.33474531727338780,  0.74472292782631944,  0.57735026918962551,  0.81975403804892744, -0.53230171741593746,  0.21132486540518711,  0.46470301232724565,  0.40254520541901589, -0.78867513459481264,
             0.81975403804892744, -0.53230171741593746,  0.21132486540518716,  0.33474531727338780,  0.74472292782631944,  0.57735026918962551,  0.46470301232724565,  0.40254520541901578, -0.78867513459481264,
             0.74472292782631933, -0.33474531727338735,  0.57735026918962562,  0.25971420705077969,  0.94227932796886926,  0.21132486540518700,  0.61476523277246187,  0.00743240513391591, -0.78867513459481264,
             0.25971420705077974,  0.94227932796886926,  0.21132486540518688,  0.74472292782631944, -0.33474531727338758,  0.57735026918962573,  0.61476523277246176,  0.00743240513391619, -0.78867513459481264,
             0.46470301232724576,  0.40254520541901567, -0.78867513459481242,  0.33474531727338774,  0.74472292782631966,  0.57735026918962551,  0.81975403804892744, -0.53230171741593746,  0.21132486540518688,
             0.33474531727338763,  0.74472292782631966,  0.57735026918962551,  0.46470301232724576,  0.40254520541901562, -0.78867513459481253,  0.81975403804892744, -0.53230171741593746,  0.21132486540518705,
             0.25971420705077974,  0.94227932796886926,  0.21132486540518683,  0.61476523277246176,  0.00743240513391624, -0.78867513459481264,  0.74472292782631944, -0.33474531727338769,  0.57735026918962573,
             0.61476523277246176,  0.00743240513391627, -0.78867513459481264,  0.25971420705077974,  0.94227932796886948,  0.21132486540518694,  0.74472292782631944, -0.33474531727338769,  0.57735026918962562,
             0.46470301232724548,  0.40254520541901584,  0.78867513459481253,  0.33474531727338785,  0.74472292782631944, -0.57735026918962573,  0.81975403804892755, -0.53230171741593746, -0.21132486540518672,
             0.33474531727338797,  0.74472292782631944, -0.57735026918962562,  0.46470301232724548,  0.40254520541901584,  0.78867513459481253,  0.81975403804892755, -0.53230171741593746, -0.21132486540518683,
             0.25971420705077980,  0.94227932796886948, -0.21132486540518722,  0.61476523277246153,  0.00743240513391627,  0.78867513459481264,  0.74472292782631955, -0.33474531727338769, -0.57735026918962540,
             0.61476523277246153,  0.00743240513391627,  0.78867513459481264,  0.25971420705077974,  0.94227932796886948, -0.21132486540518727,  0.74472292782631966, -0.33474531727338780, -0.57735026918962529,
             0.00743240513391610,  0.61476523277246142,  0.78867513459481264, -0.33474531727338780,  0.74472292782631966, -0.57735026918962573,  0.94227932796886948,  0.25971420705077991, -0.21132486540518700,
            -0.33474531727338769,  0.74472292782631955, -0.57735026918962540,  0.00743240513391613,  0.61476523277246142,  0.78867513459481253,  0.94227932796886926,  0.25971420705078002, -0.21132486540518716,
            -0.53230171741593746,  0.81975403804892766, -0.21132486540518688,  0.40254520541901584,  0.46470301232724537,  0.78867513459481264,  0.74472292782631966,  0.33474531727338797, -0.57735026918962573,
             0.40254520541901573,  0.46470301232724542,  0.78867513459481264, -0.53230171741593768,  0.81975403804892766, -0.21132486540518691,  0.74472292782631988,  0.33474531727338797, -0.57735026918962573,
             0.00743240513391608,  0.61476523277246198, -0.78867513459481253, -0.33474531727338769,  0.74472292782631921,  0.57735026918962584,  0.94227932796886948,  0.25971420705077980,  0.21132486540518672,
            -0.33474531727338769,  0.74472292782631921,  0.57735026918962584,  0.00743240513391613,  0.61476523277246187, -0.78867513459481275,  0.94227932796886948,  0.25971420705077986,  0.21132486540518688,
            -0.53230171741593746,  0.81975403804892744,  0.21132486540518727,  0.40254520541901567,  0.46470301232724587, -0.78867513459481287,  0.74472292782631966,  0.33474531727338763,  0.57735026918962540,
             0.40254520541901562,  0.46470301232724587, -0.78867513459481275, -0.53230171741593746,  0.81975403804892744,  0.21132486540518730,  0.74472292782631966,  0.33474531727338763,  0.57735026918962540,
             0.74472292782631944,  0.33474531727338780,  0.57735026918962551, -0.53230171741593746,  0.81975403804892744,  0.21132486540518733,  0.40254520541901589,  0.46470301232724576, -0.78867513459481287,
            -0.53230171741593746,  0.81975403804892732,  0.21132486540518727,  0.74472292782631966,  0.33474531727338774,  0.57735026918962551,  0.40254520541901584,  0.46470301232724587, -0.78867513459481275,
            -0.33474531727338769,  0.74472292782631933,  0.57735026918962595,  0.94227932796886948,  0.25971420705077974,  0.21132486540518675,  0.00743240513391596,  0.61476523277246198, -0.78867513459481275,
             0.94227932796886948,  0.25971420705077986,  0.21132486540518680, -0.33474531727338758,  0.74472292782631921,  0.57735026918962595,  0.00743240513391613,  0.61476523277246187, -0.78867513459481287,
             0.74472292782631966,  0.33474531727338797, -0.57735026918962595, -0.53230171741593768,  0.81975403804892755, -0.21132486540518680,  0.40254520541901595,  0.46470301232724537,  0.78867513459481275,
            -0.53230171741593746,  0.81975403804892744, -0.21132486540518686,  0.74472292782631944,  0.33474531727338802, -0.57735026918962595,  0.40254520541901578,  0.46470301232724542,  0.78867513459481287,
            -0.33474531727338758,  0.74472292782631944, -0.57735026918962551,  0.94227932796886948,  0.25971420705077991, -0.21132486540518738,  0.00743240513391602,  0.61476523277246153,  0.78867513459481287,
             0.94227932796886948,  0.25971420705077997, -0.21132486540518733, -0.33474531727338769,  0.74472292782631955, -0.57735026918962551,  0.00743240513391608,  0.61476523277246142,  0.78867513459481275,
             0.40254520541901584,  0.46470301232724537,  0.78867513459481253,  0.74472292782631966,  0.33474531727338802, -0.57735026918962584, -0.53230171741593746,  0.81975403804892744, -0.21132486540518677,
             0.74472292782631966,  0.33474531727338802, -0.57735026918962584,  0.40254520541901578,  0.46470301232724542,  0.78867513459481275, -0.53230171741593746,  0.81975403804892744, -0.21132486540518694,
             0.94227932796886948,  0.25971420705077986, -0.21132486540518733,  0.00743240513391624,  0.61476523277246142,  0.78867513459481287, -0.33474531727338780,  0.74472292782631966, -0.57735026918962540,
             0.00743240513391633,  0.61476523277246142,  0.78867513459481275,  0.94227932796886926,  0.25971420705077980, -0.21132486540518736, -0.33474531727338780,  0.74472292782631966, -0.57735026918962540,
             0.40254520541901584,  0.46470301232724587, -0.78867513459481264,  0.74472292782631966,  0.33474531727338763,  0.57735026918962573, -0.53230171741593768,  0.81975403804892744,  0.21132486540518694,
             0.74472292782631966,  0.33474531727338774,  0.57735026918962540,  0.40254520541901578,  0.46470301232724587, -0.78867513459481253, -0.53230171741593746,  0.81975403804892721,  0.21132486540518711,
             0.94227932796886948,  0.25971420705077969,  0.21132486540518683,  0.00743240513391610,  0.61476523277246187, -0.78867513459481264, -0.33474531727338769,  0.74472292782631933,  0.57735026918962573,
             0.00743240513391616,  0.61476523277246187, -0.78867513459481264,  0.94227932796886971,  0.25971420705077969,  0.21132486540518686, -0.33474531727338791,  0.74472292782631933,  0.57735026918962573,
             0.90148363992456226, -0.21387143172584827, -0.37627949404651062, -0.14283988479115434,  0.67365669613657153, -0.72510924905369145,  0.40856335403608074,  0.70742184475876524,  0.57674111982491172,
            -0.14283988479115423,  0.67365669613657153, -0.72510924905369145,  0.90148363992456226, -0.21387143172584827, -0.37627949404651062,  0.40856335403608057,  0.70742184475876535,  0.57674111982491172,
             0.07071622802089400,  0.36957471874357817, -0.92650620200843870,  0.99200950450550729, -0.12334556714490302,  0.02651441186298359,  0.10448137664308738,  0.92097795757081347,  0.37534416687016475,
             0.99200950450550707, -0.12334556714490302,  0.02651441186298353,  0.07071622802089389,  0.36957471874357817, -0.92650620200843870,  0.10448137664308754,  0.92097795757081347,  0.37534416687016486,
            -0.21387143172584805,  0.90148363992456226, -0.37627949404651045,  0.67365669613657175, -0.14283988479115456, -0.72510924905369167,  0.70742184475876502,  0.40856335403608091,  0.57674111982491172,
             0.67365669613657164, -0.14283988479115456, -0.72510924905369167, -0.21387143172584805,  0.90148363992456226, -0.37627949404651045,  0.70742184475876502,  0.40856335403608085,  0.57674111982491172,
             0.36957471874357850,  0.07071622802089367, -0.92650620200843870, -0.12334556714490313,  0.99200950450550729,  0.02651441186298362,  0.92097795757081324,  0.10448137664308760,  0.37534416687016464,
            -0.12334556714490313,  0.99200950450550729,  0.02651441186298359,  0.36957471874357850,  0.07071622802089356, -0.92650620200843870,  0.92097795757081324,  0.10448137664308771,  0.37534416687016475,
             0.92097795757081324,  0.10448137664308760,  0.37534416687016453,  0.36957471874357850,  0.07071622802089361, -0.92650620200843870, -0.12334556714490302,  0.99200950450550729,  0.02651441186298392,
             0.36957471874357839,  0.07071622802089378, -0.92650620200843847,  0.92097795757081347,  0.10448137664308743,  0.37534416687016464, -0.12334556714490325,  0.99200950450550729,  0.02651441186298359,
             0.67365669613657186, -0.14283988479115456, -0.72510924905369167,  0.70742184475876502,  0.40856335403608096,  0.57674111982491194, -0.21387143172584827,  0.90148363992456226, -0.37627949404651068,
             0.70742184475876502,  0.40856335403608091,  0.57674111982491172,  0.67365669613657175, -0.14283988479115445, -0.72510924905369145, -0.21387143172584816,  0.90148363992456226, -0.37627949404651045,
             0.10448137664308721,  0.92097795757081369,  0.37534416687016475,  0.07071622802089400,  0.36957471874357811, -0.92650620200843870,  0.99200950450550729, -0.12334556714490313,  0.02651441186298381,
             0.07071622802089400,  0.36957471874357817, -0.92650620200843847,  0.10448137664308732,  0.92097795757081347,  0.37534416687016464,  0.99200950450550729, -0.12334556714490313,  0.02651441186298359,
            -0.14283988479115423,  0.67365669613657153, -0.72510924905369145,  0.40856335403608057,  0.70742184475876546,  0.57674111982491194,  0.90148363992456204, -0.21387143172584827, -0.37627949404651079,
             0.40856335403608068,  0.70742184475876524,  0.57674111982491172, -0.14283988479115445,  0.67365669613657164, -0.72510924905369145,  0.90148363992456226, -0.21387143172584827, -0.37627949404651045,
             0.99200950450550707, -0.12334556714490313,  0.02651441186298359,  0.10448137664308738,  0.92097795757081369,  0.37534416687016475,  0.07071622802089406,  0.36957471874357817, -0.92650620200843870,
             0.10448137664308743,  0.92097795757081369,  0.37534416687016475,  0.99200950450550707, -0.12334556714490325,  0.02651441186298359,  0.07071622802089406,  0.36957471874357822, -0.92650620200843847,
             0.40856335403608057,  0.70742184475876546,  0.57674111982491194,  0.90148363992456226, -0.21387143172584827, -0.37627949404651051, -0.14283988479115423,  0.67365669613657153, -0.72510924905369145,
             0.90148363992456226, -0.21387143172584827, -0.37627949404651045,  0.40856335403608057,  0.70742184475876546,  0.57674111982491194, -0.14283988479115423,  0.67365669613657131, -0.72510924905369167,
            -0.12334556714490313,  0.99200950450550729,  0.02651441186298378,  0.92097795757081347,  0.10448137664308754,  0.37534416687016464,  0.36957471874357833,  0.07071622802089378, -0.92650620200843870,
             0.92097795757081324,  0.10448137664308760,  0.37534416687016464, -0.12334556714490325,  0.99200950450550729,  0.02651441186298378,  0.36957471874357850,  0.07071622802089372, -0.92650620200843847,
             0.70742184475876502,  0.40856335403608091,  0.57674111982491194, -0.21387143172584816,  0.90148363992456204, -0.37627949404651045,  0.67365669613657175, -0.14283988479115445, -0.72510924905369167,
            -0.21387143172584805,  0.90148363992456204, -0.37627949404651040,  0.70742184475876524,  0.40856335403608091,  0.57674111982491194,  0.67365669613657153, -0.14283988479115445, -0.72510924905369167,
        });
    }

    private static void populateFrames3D_Part2() {
        populateFrames3D(2592, new double[] {
            -0.00743240513391594,  0.61476523277246187, -0.78867513459481287,  0.33474531727338758,  0.74472292782631944,  0.57735026918962584, -0.94227932796886948,  0.25971420705077963,  0.21132486540518694,
             0.33474531727338758,  0.74472292782631944,  0.57735026918962584, -0.00743240513391605,  0.61476523277246187, -0.78867513459481287, -0.94227932796886926,  0.25971420705077974,  0.21132486540518700,
             0.53230171741593746,  0.81975403804892766,  0.21132486540518711, -0.40254520541901562,  0.46470301232724570, -0.78867513459481287, -0.74472292782631966,  0.33474531727338763,  0.57735026918962573,
            -0.40254520541901556,  0.46470301232724581, -0.78867513459481275,  0.53230171741593746,  0.81975403804892766,  0.21132486540518711, -0.74472292782631988,  0.33474531727338752,  0.57735026918962551,
            -0.00743240513391624,  0.61476523277246153,  0.78867513459481298,  0.33474531727338769,  0.74472292782631955, -0.57735026918962573, -0.94227932796886926,  0.25971420705077986, -0.21132486540518716,
             0.33474531727338780,  0.74472292782631966, -0.57735026918962573, -0.00743240513391630,  0.61476523277246153,  0.78867513459481298, -0.94227932796886926,  0.25971420705077986, -0.21132486540518722,
             0.53230171741593746,  0.81975403804892766, -0.21132486540518711, -0.40254520541901584,  0.46470301232724548,  0.78867513459481287, -0.74472292782631966,  0.33474531727338780, -0.57735026918962573,
            -0.40254520541901573,  0.46470301232724548,  0.78867513459481287,  0.53230171741593746,  0.81975403804892766, -0.21132486540518716, -0.74472292782631966,  0.33474531727338780, -0.57735026918962573,
            -0.74472292782631944,  0.33474531727338797, -0.57735026918962584,  0.53230171741593768,  0.81975403804892766, -0.21132486540518716, -0.40254520541901601,  0.46470301232724531,  0.78867513459481287,
             0.53230171741593746,  0.81975403804892766, -0.21132486540518713, -0.74472292782631944,  0.33474531727338785, -0.57735026918962584, -0.40254520541901595,  0.46470301232724548,  0.78867513459481298,
             0.33474531727338758,  0.74472292782631966, -0.57735026918962595, -0.94227932796886948,  0.25971420705077974, -0.21132486540518702, -0.00743240513391608,  0.61476523277246165,  0.78867513459481287,
            -0.94227932796886926,  0.25971420705077986, -0.21132486540518716,  0.33474531727338769,  0.74472292782631955, -0.57735026918962584, -0.00743240513391619,  0.61476523277246153,  0.78867513459481287,
            -0.74472292782631944,  0.33474531727338769,  0.57735026918962595,  0.53230171741593746,  0.81975403804892755,  0.21132486540518697, -0.40254520541901584,  0.46470301232724570, -0.78867513459481287,
             0.53230171741593746,  0.81975403804892744,  0.21132486540518711, -0.74472292782631966,  0.33474531727338774,  0.57735026918962584, -0.40254520541901573,  0.46470301232724570, -0.78867513459481287,
             0.33474531727338747,  0.74472292782631933,  0.57735026918962584, -0.94227932796886948,  0.25971420705077963,  0.21132486540518711, -0.00743240513391591,  0.61476523277246198, -0.78867513459481287,
            -0.94227932796886948,  0.25971420705077969,  0.21132486540518708,  0.33474531727338758,  0.74472292782631944,  0.57735026918962584, -0.00743240513391596,  0.61476523277246176, -0.78867513459481298,
            -0.40254520541901567,  0.46470301232724570, -0.78867513459481298, -0.74472292782631966,  0.33474531727338774,  0.57735026918962573,  0.53230171741593746,  0.81975403804892744,  0.21132486540518711,
            -0.74472292782631966,  0.33474531727338763,  0.57735026918962573, -0.40254520541901562,  0.46470301232724576, -0.78867513459481298,  0.53230171741593746,  0.81975403804892744,  0.21132486540518716,
            -0.94227932796886948,  0.25971420705077963,  0.21132486540518705, -0.00743240513391610,  0.61476523277246176, -0.78867513459481287,  0.33474531727338769,  0.74472292782631944,  0.57735026918962573,
            -0.00743240513391619,  0.61476523277246187, -0.78867513459481287, -0.94227932796886948,  0.25971420705077963,  0.21132486540518711,  0.33474531727338780,  0.74472292782631944,  0.57735026918962573,
            -0.40254520541901595,  0.46470301232724548,  0.78867513459481287, -0.74472292782631944,  0.33474531727338780, -0.57735026918962584,  0.53230171741593746,  0.81975403804892766, -0.21132486540518700,
            -0.74472292782631944,  0.33474531727338785, -0.57735026918962584, -0.40254520541901584,  0.46470301232724542,  0.78867513459481287,  0.53230171741593746,  0.81975403804892755, -0.21132486540518705,
            -0.94227932796886948,  0.25971420705077969, -0.21132486540518716, -0.00743240513391627,  0.61476523277246153,  0.78867513459481287,  0.33474531727338780,  0.74472292782631966, -0.57735026918962573,
            -0.00743240513391635,  0.61476523277246153,  0.78867513459481275, -0.94227932796886926,  0.25971420705077969, -0.21132486540518716,  0.33474531727338791,  0.74472292782631977, -0.57735026918962551,
            -0.90148363992456237, -0.21387143172584849,  0.37627949404651073,  0.14283988479115439,  0.67365669613657153,  0.72510924905369178, -0.40856335403608085,  0.70742184475876546, -0.57674111982491216,
             0.14283988479115417,  0.67365669613657153,  0.72510924905369145, -0.90148363992456226, -0.21387143172584849,  0.37627949404651073, -0.40856335403608057,  0.70742184475876546, -0.57674111982491194,
            -0.07071622802089395,  0.36957471874357817,  0.92650620200843892, -0.99200950450550729, -0.12334556714490325, -0.02651441186298364, -0.10448137664308754,  0.92097795757081369, -0.37534416687016503,
            -0.99200950450550729, -0.12334556714490313, -0.02651441186298356, -0.07071622802089395,  0.36957471874357817,  0.92650620200843881, -0.10448137664308754,  0.92097795757081369, -0.37534416687016503,
             0.21387143172584833,  0.90148363992456226,  0.37627949404651051, -0.67365669613657186, -0.14283988479115467,  0.72510924905369190, -0.70742184475876524,  0.40856335403608096, -0.57674111982491216,
            -0.67365669613657175, -0.14283988479115489,  0.72510924905369190,  0.21387143172584810,  0.90148363992456249,  0.37627949404651045, -0.70742184475876524,  0.40856335403608096, -0.57674111982491216,
            -0.36957471874357845,  0.07071622802089356,  0.92650620200843892,  0.12334556714490319,  0.99200950450550751, -0.02651441186298373, -0.92097795757081347,  0.10448137664308760, -0.37534416687016492,
             0.12334556714490319,  0.99200950450550751, -0.02651441186298373, -0.36957471874357850,  0.07071622802089350,  0.92650620200843892, -0.92097795757081347,  0.10448137664308760, -0.37534416687016492,
            -0.92097795757081347,  0.10448137664308743, -0.37534416687016481, -0.36957471874357850,  0.07071622802089350,  0.92650620200843892,  0.12334556714490319,  0.99200950450550751, -0.02651441186298401,
            -0.36957471874357850,  0.07071622802089356,  0.92650620200843903, -0.92097795757081358,  0.10448137664308749, -0.37534416687016503,  0.12334556714490319,  0.99200950450550751, -0.02651441186298378,
            -0.67365669613657198, -0.14283988479115478,  0.72510924905369178, -0.70742184475876513,  0.40856335403608102, -0.57674111982491227,  0.21387143172584822,  0.90148363992456226,  0.37627949404651062,
            -0.70742184475876513,  0.40856335403608091, -0.57674111982491205, -0.67365669613657198, -0.14283988479115478,  0.72510924905369178,  0.21387143172584822,  0.90148363992456249,  0.37627949404651051,
            -0.10448137664308726,  0.92097795757081391, -0.37534416687016503, -0.07071622802089406,  0.36957471874357806,  0.92650620200843903, -0.99200950450550740, -0.12334556714490325, -0.02651441186298384,
            -0.07071622802089406,  0.36957471874357817,  0.92650620200843881, -0.10448137664308726,  0.92097795757081391, -0.37534416687016503, -0.99200950450550740, -0.12334556714490336, -0.02651441186298367,
             0.14283988479115428,  0.67365669613657164,  0.72510924905369156, -0.40856335403608068,  0.70742184475876557, -0.57674111982491216, -0.90148363992456237, -0.21387143172584849,  0.37627949404651084,
            -0.40856335403608068,  0.70742184475876546, -0.57674111982491227,  0.14283988479115439,  0.67365669613657153,  0.72510924905369178, -0.90148363992456237, -0.21387143172584849,  0.37627949404651062,
            -0.99200950450550751, -0.12334556714490325, -0.02651441186298373, -0.10448137664308732,  0.92097795757081369, -0.37534416687016514, -0.07071622802089400,  0.36957471874357811,  0.92650620200843892,
            -0.10448137664308743,  0.92097795757081391, -0.37534416687016514, -0.99200950450550729, -0.12334556714490336, -0.02651441186298364, -0.07071622802089400,  0.36957471874357811,  0.92650620200843892,
            -0.40856335403608074,  0.70742184475876546, -0.57674111982491216, -0.90148363992456237, -0.21387143172584849,  0.37627949404651051,  0.14283988479115428,  0.67365669613657153,  0.72510924905369167,
            -0.90148363992456237, -0.21387143172584849,  0.37627949404651051, -0.40856335403608068,  0.70742184475876557, -0.57674111982491216,  0.14283988479115428,  0.67365669613657153,  0.72510924905369167,
             0.12334556714490319,  0.99200950450550751, -0.02651441186298395, -0.92097795757081358,  0.10448137664308749, -0.37534416687016503, -0.36957471874357833,  0.07071622802089361,  0.92650620200843892,
            -0.92097795757081335,  0.10448137664308760, -0.37534416687016470,  0.12334556714490308,  0.99200950450550751, -0.02651441186298389, -0.36957471874357861,  0.07071622802089356,  0.92650620200843870,
            -0.70742184475876524,  0.40856335403608091, -0.57674111982491216,  0.21387143172584810,  0.90148363992456226,  0.37627949404651045, -0.67365669613657164, -0.14283988479115456,  0.72510924905369178,
             0.21387143172584810,  0.90148363992456226,  0.37627949404651040, -0.70742184475876524,  0.40856335403608091, -0.57674111982491205, -0.67365669613657164, -0.14283988479115456,  0.72510924905369178,
             0.21387143172584827,  0.90148363992456226, -0.37627949404651062, -0.67365669613657153, -0.14283988479115450, -0.72510924905369178, -0.70742184475876524,  0.40856335403608096,  0.57674111982491216,
            -0.67365669613657153, -0.14283988479115417, -0.72510924905369156,  0.21387143172584827,  0.90148363992456226, -0.37627949404651068, -0.70742184475876535,  0.40856335403608068,  0.57674111982491205,
            -0.36957471874357817,  0.07071622802089389, -0.92650620200843892,  0.12334556714490325,  0.99200950450550729,  0.02651441186298364, -0.92097795757081347,  0.10448137664308754,  0.37534416687016503,
             0.12334556714490325,  0.99200950450550729,  0.02651441186298356, -0.36957471874357822,  0.07071622802089383, -0.92650620200843881, -0.92097795757081347,  0.10448137664308760,  0.37534416687016503,
            -0.90148363992456226, -0.21387143172584822, -0.37627949404651051,  0.14283988479115456,  0.67365669613657175, -0.72510924905369190, -0.40856335403608085,  0.70742184475876524,  0.57674111982491216,
             0.14283988479115467,  0.67365669613657153, -0.72510924905369190, -0.90148363992456226, -0.21387143172584810, -0.37627949404651051, -0.40856335403608085,  0.70742184475876524,  0.57674111982491205,
            -0.07071622802089361,  0.36957471874357839, -0.92650620200843881, -0.99200950450550751, -0.12334556714490319,  0.02651441186298376, -0.10448137664308754,  0.92097795757081358,  0.37534416687016492,
            -0.99200950450550729, -0.12334556714490308,  0.02651441186298370, -0.07071622802089367,  0.36957471874357839, -0.92650620200843892, -0.10448137664308765,  0.92097795757081347,  0.37534416687016503,
            -0.10448137664308738,  0.92097795757081347,  0.37534416687016470, -0.07071622802089367,  0.36957471874357839, -0.92650620200843892, -0.99200950450550751, -0.12334556714490308,  0.02651441186298395,
            -0.07071622802089367,  0.36957471874357839, -0.92650620200843903, -0.10448137664308749,  0.92097795757081358,  0.37534416687016503, -0.99200950450550751, -0.12334556714490319,  0.02651441186298378,
             0.14283988479115467,  0.67365669613657175, -0.72510924905369167, -0.40856335403608091,  0.70742184475876524,  0.57674111982491227, -0.90148363992456226, -0.21387143172584822, -0.37627949404651062,
            -0.40856335403608085,  0.70742184475876524,  0.57674111982491205,  0.14283988479115456,  0.67365669613657175, -0.72510924905369178, -0.90148363992456226, -0.21387143172584810, -0.37627949404651051,
            -0.92097795757081369,  0.10448137664308743,  0.37534416687016492, -0.36957471874357817,  0.07071622802089389, -0.92650620200843903,  0.12334556714490325,  0.99200950450550740,  0.02651441186298378,
            -0.36957471874357822,  0.07071622802089395, -0.92650620200843881, -0.92097795757081369,  0.10448137664308743,  0.37534416687016503,  0.12334556714490325,  0.99200950450550729,  0.02651441186298367,
            -0.67365669613657175, -0.14283988479115428, -0.72510924905369145, -0.70742184475876546,  0.40856335403608079,  0.57674111982491216,  0.21387143172584849,  0.90148363992456226, -0.37627949404651073,
            -0.70742184475876546,  0.40856335403608079,  0.57674111982491227, -0.67365669613657153, -0.14283988479115439, -0.72510924905369178,  0.21387143172584838,  0.90148363992456237, -0.37627949404651062,
             0.12334556714490325,  0.99200950450550740,  0.02651441186298373, -0.92097795757081369,  0.10448137664308749,  0.37534416687016514, -0.36957471874357822,  0.07071622802089395, -0.92650620200843892,
            -0.92097795757081369,  0.10448137664308760,  0.37534416687016514,  0.12334556714490325,  0.99200950450550729,  0.02651441186298367, -0.36957471874357822,  0.07071622802089389, -0.92650620200843881,
            -0.70742184475876546,  0.40856335403608079,  0.57674111982491205,  0.21387143172584849,  0.90148363992456237, -0.37627949404651057, -0.67365669613657153, -0.14283988479115439, -0.72510924905369167,
             0.21387143172584827,  0.90148363992456226, -0.37627949404651051, -0.70742184475876546,  0.40856335403608079,  0.57674111982491216, -0.67365669613657142, -0.14283988479115428, -0.72510924905369178,
            -0.99200950450550729, -0.12334556714490308,  0.02651441186298384, -0.10448137664308749,  0.92097795757081369,  0.37534416687016503, -0.07071622802089378,  0.36957471874357822, -0.92650620200843892,
            -0.10448137664308749,  0.92097795757081335,  0.37534416687016481, -0.99200950450550729, -0.12334556714490308,  0.02651441186298387, -0.07071622802089372,  0.36957471874357850, -0.92650620200843881,
            -0.40856335403608091,  0.70742184475876524,  0.57674111982491216, -0.90148363992456226, -0.21387143172584810, -0.37627949404651045,  0.14283988479115445,  0.67365669613657164, -0.72510924905369178,
            -0.90148363992456226, -0.21387143172584810, -0.37627949404651040, -0.40856335403608085,  0.70742184475876535,  0.57674111982491205,  0.14283988479115445,  0.67365669613657153, -0.72510924905369178,
            -0.61476523277246187,  0.00743240513391594,  0.78867513459481253, -0.74472292782631944, -0.33474531727338758, -0.57735026918962551, -0.25971420705077963,  0.94227932796886948, -0.21132486540518683,
            -0.74472292782631944, -0.33474531727338758, -0.57735026918962540, -0.61476523277246176,  0.00743240513391599,  0.78867513459481242, -0.25971420705077974,  0.94227932796886926, -0.21132486540518700,
            -0.81975403804892766, -0.53230171741593746, -0.21132486540518697, -0.46470301232724559,  0.40254520541901562,  0.78867513459481264, -0.33474531727338774,  0.74472292782631966, -0.57735026918962551,
            -0.46470301232724570,  0.40254520541901567,  0.78867513459481253, -0.81975403804892744, -0.53230171741593768, -0.21132486540518713, -0.33474531727338774,  0.74472292782631988, -0.57735026918962540,
            -0.61476523277246165,  0.00743240513391624, -0.78867513459481242, -0.74472292782631966, -0.33474531727338769,  0.57735026918962551, -0.25971420705077974,  0.94227932796886926,  0.21132486540518694,
            -0.74472292782631955, -0.33474531727338780,  0.57735026918962551, -0.61476523277246165,  0.00743240513391630, -0.78867513459481264, -0.25971420705077980,  0.94227932796886926,  0.21132486540518705,
            -0.81975403804892766, -0.53230171741593746,  0.21132486540518702, -0.46470301232724553,  0.40254520541901584, -0.78867513459481264, -0.33474531727338780,  0.74472292782631966,  0.57735026918962551,
            -0.46470301232724553,  0.40254520541901573, -0.78867513459481264, -0.81975403804892766, -0.53230171741593746,  0.21132486540518711, -0.33474531727338774,  0.74472292782631966,  0.57735026918962540,
            -0.33474531727338791,  0.74472292782631944,  0.57735026918962562, -0.81975403804892755, -0.53230171741593746,  0.21132486540518711, -0.46470301232724548,  0.40254520541901595, -0.78867513459481264,
            -0.81975403804892766, -0.53230171741593746,  0.21132486540518711, -0.33474531727338785,  0.74472292782631944,  0.57735026918962573, -0.46470301232724548,  0.40254520541901589, -0.78867513459481264,
            -0.74472292782631955, -0.33474531727338780,  0.57735026918962584, -0.25971420705077963,  0.94227932796886948,  0.21132486540518702, -0.61476523277246176,  0.00743240513391613, -0.78867513459481264,
            -0.25971420705077974,  0.94227932796886926,  0.21132486540518700, -0.74472292782631955, -0.33474531727338769,  0.57735026918962562, -0.61476523277246165,  0.00743240513391624, -0.78867513459481253,
            -0.33474531727338774,  0.74472292782631966, -0.57735026918962584, -0.81975403804892766, -0.53230171741593746, -0.21132486540518697, -0.46470301232724553,  0.40254520541901578,  0.78867513459481264,
            -0.81975403804892755, -0.53230171741593746, -0.21132486540518694, -0.33474531727338774,  0.74472292782631966, -0.57735026918962562, -0.46470301232724565,  0.40254520541901567,  0.78867513459481253,
            -0.74472292782631944, -0.33474531727338747, -0.57735026918962562, -0.25971420705077974,  0.94227932796886948, -0.21132486540518705, -0.61476523277246176,  0.00743240513391596,  0.78867513459481264,
            -0.25971420705077963,  0.94227932796886948, -0.21132486540518705, -0.74472292782631944, -0.33474531727338758, -0.57735026918962573, -0.61476523277246187,  0.00743240513391602,  0.78867513459481264,
            -0.46470301232724565,  0.40254520541901567,  0.78867513459481242, -0.33474531727338769,  0.74472292782631966, -0.57735026918962551, -0.81975403804892755, -0.53230171741593746, -0.21132486540518688,
            -0.33474531727338774,  0.74472292782631966, -0.57735026918962551, -0.46470301232724570,  0.40254520541901562,  0.78867513459481264, -0.81975403804892744, -0.53230171741593746, -0.21132486540518700,
            -0.25971420705077969,  0.94227932796886948, -0.21132486540518697, -0.61476523277246176,  0.00743240513391608,  0.78867513459481264, -0.74472292782631944, -0.33474531727338769, -0.57735026918962551,
            -0.61476523277246176,  0.00743240513391619,  0.78867513459481264, -0.25971420705077969,  0.94227932796886948, -0.21132486540518705, -0.74472292782631955, -0.33474531727338780, -0.57735026918962540,
            -0.46470301232724548,  0.40254520541901595, -0.78867513459481253, -0.33474531727338780,  0.74472292782631944,  0.57735026918962551, -0.81975403804892766, -0.53230171741593746,  0.21132486540518688,
            -0.33474531727338785,  0.74472292782631944,  0.57735026918962540, -0.46470301232724553,  0.40254520541901595, -0.78867513459481242, -0.81975403804892755, -0.53230171741593746,  0.21132486540518705,
            -0.25971420705077969,  0.94227932796886948,  0.21132486540518702, -0.61476523277246165,  0.00743240513391630, -0.78867513459481264, -0.74472292782631955, -0.33474531727338780,  0.57735026918962551,
            -0.61476523277246165,  0.00743240513391624, -0.78867513459481253, -0.25971420705077980,  0.94227932796886971,  0.21132486540518719, -0.74472292782631955, -0.33474531727338791,  0.57735026918962540,
            -0.00743240513391588, -0.78867513459481287,  0.61476523277246165,  0.33474531727338763,  0.57735026918962584,  0.74472292782631966, -0.94227932796886937,  0.21132486540518702,  0.25971420705077936,
             0.33474531727338752,  0.57735026918962562,  0.74472292782631944, -0.00743240513391596, -0.78867513459481287,  0.61476523277246176, -0.94227932796886926,  0.21132486540518713,  0.25971420705077952,
             0.53230171741593746,  0.21132486540518713,  0.81975403804892766, -0.40254520541901551, -0.78867513459481287,  0.46470301232724548, -0.74472292782631955,  0.57735026918962573,  0.33474531727338752,
            -0.40254520541901551, -0.78867513459481264,  0.46470301232724553,  0.53230171741593746,  0.21132486540518716,  0.81975403804892766, -0.74472292782631977,  0.57735026918962551,  0.33474531727338752,
            -0.00743240513391630,  0.78867513459481287,  0.61476523277246165,  0.33474531727338785, -0.57735026918962584,  0.74472292782631955, -0.94227932796886926, -0.21132486540518697,  0.25971420705077941,
             0.33474531727338797, -0.57735026918962584,  0.74472292782631966, -0.00743240513391638,  0.78867513459481287,  0.61476523277246153, -0.94227932796886926, -0.21132486540518702,  0.25971420705077941,
             0.53230171741593746, -0.21132486540518716,  0.81975403804892766, -0.40254520541901573,  0.78867513459481287,  0.46470301232724537, -0.74472292782631944, -0.57735026918962562,  0.33474531727338763,
            -0.40254520541901573,  0.78867513459481287,  0.46470301232724542,  0.53230171741593757, -0.21132486540518722,  0.81975403804892766, -0.74472292782631944, -0.57735026918962562,  0.33474531727338758,
            -0.74472292782631921, -0.57735026918962573,  0.33474531727338758,  0.53230171741593757, -0.21132486540518727,  0.81975403804892766, -0.40254520541901612,  0.78867513459481287,  0.46470301232724531,
             0.53230171741593746, -0.21132486540518725,  0.81975403804892766, -0.74472292782631921, -0.57735026918962584,  0.33474531727338752, -0.40254520541901589,  0.78867513459481287,  0.46470301232724531,
             0.33474531727338785, -0.57735026918962595,  0.74472292782631966, -0.94227932796886937, -0.21132486540518700,  0.25971420705077941, -0.00743240513391619,  0.78867513459481287,  0.61476523277246165,
            -0.94227932796886926, -0.21132486540518700,  0.25971420705077947,  0.33474531727338785, -0.57735026918962595,  0.74472292782631955, -0.00743240513391635,  0.78867513459481287,  0.61476523277246153,
            -0.74472292782631966,  0.57735026918962595,  0.33474531727338747,  0.53230171741593757,  0.21132486540518700,  0.81975403804892766, -0.40254520541901562, -0.78867513459481287,  0.46470301232724548,
             0.53230171741593746,  0.21132486540518700,  0.81975403804892766, -0.74472292782631966,  0.57735026918962595,  0.33474531727338752, -0.40254520541901545, -0.78867513459481287,  0.46470301232724553,
             0.33474531727338741,  0.57735026918962573,  0.74472292782631944, -0.94227932796886937,  0.21132486540518727,  0.25971420705077936, -0.00743240513391569, -0.78867513459481287,  0.61476523277246176,
            -0.94227932796886926,  0.21132486540518725,  0.25971420705077936,  0.33474531727338741,  0.57735026918962584,  0.74472292782631955, -0.00743240513391591, -0.78867513459481287,  0.61476523277246176,
            -0.40254520541901551, -0.78867513459481287,  0.46470301232724542, -0.74472292782631966,  0.57735026918962584,  0.33474531727338752,  0.53230171741593746,  0.21132486540518697,  0.81975403804892766,
            -0.74472292782631977,  0.57735026918962584,  0.33474531727338741, -0.40254520541901539, -0.78867513459481287,  0.46470301232724553,  0.53230171741593746,  0.21132486540518702,  0.81975403804892766,
            -0.94227932796886926,  0.21132486540518716,  0.25971420705077941, -0.00743240513391605, -0.78867513459481287,  0.61476523277246165,  0.33474531727338763,  0.57735026918962562,  0.74472292782631944,
            -0.00743240513391608, -0.78867513459481287,  0.61476523277246165, -0.94227932796886937,  0.21132486540518722,  0.25971420705077941,  0.33474531727338763,  0.57735026918962562,  0.74472292782631944,
            -0.40254520541901595,  0.78867513459481287,  0.46470301232724542, -0.74472292782631944, -0.57735026918962584,  0.33474531727338747,  0.53230171741593757, -0.21132486540518702,  0.81975403804892766,
            -0.74472292782631933, -0.57735026918962562,  0.33474531727338763, -0.40254520541901584,  0.78867513459481287,  0.46470301232724537,  0.53230171741593746, -0.21132486540518713,  0.81975403804892755,
            -0.94227932796886926, -0.21132486540518713,  0.25971420705077941, -0.00743240513391630,  0.78867513459481287,  0.61476523277246165,  0.33474531727338774, -0.57735026918962573,  0.74472292782631955,
            -0.00743240513391630,  0.78867513459481264,  0.61476523277246153, -0.94227932796886926, -0.21132486540518716,  0.25971420705077947,  0.33474531727338797, -0.57735026918962551,  0.74472292782631955,
            -0.90148363992456204,  0.37627949404651051, -0.21387143172584844,  0.14283988479115445,  0.72510924905369190,  0.67365669613657142, -0.40856335403608085, -0.57674111982491216,  0.70742184475876524,
             0.14283988479115434,  0.72510924905369167,  0.67365669613657142, -0.90148363992456204,  0.37627949404651051, -0.21387143172584844, -0.40856335403608052, -0.57674111982491194,  0.70742184475876524,
            -0.07071622802089378,  0.92650620200843892,  0.36957471874357800, -0.99200950450550707, -0.02651441186298384, -0.12334556714490330, -0.10448137664308743, -0.37534416687016492,  0.92097795757081347,
            -0.99200950450550707, -0.02651441186298376, -0.12334556714490341, -0.07071622802089378,  0.92650620200843892,  0.36957471874357806, -0.10448137664308749, -0.37534416687016492,  0.92097795757081335,
             0.21387143172584805,  0.37627949404651062,  0.90148363992456215, -0.67365669613657153,  0.72510924905369167, -0.14283988479115484, -0.70742184475876480, -0.57674111982491216,  0.40856335403608074,
            -0.67365669613657153,  0.72510924905369178, -0.14283988479115484,  0.21387143172584827,  0.37627949404651051,  0.90148363992456237, -0.70742184475876502, -0.57674111982491194,  0.40856335403608063,
            -0.36957471874357833,  0.92650620200843881,  0.07071622802089350,  0.12334556714490325, -0.02651441186298362,  0.99200950450550729, -0.92097795757081324, -0.37534416687016503,  0.10448137664308738,
             0.12334556714490325, -0.02651441186298359,  0.99200950450550729, -0.36957471874357828,  0.92650620200843881,  0.07071622802089339, -0.92097795757081324, -0.37534416687016503,  0.10448137664308743,
            -0.92097795757081324, -0.37534416687016492,  0.10448137664308732, -0.36957471874357822,  0.92650620200843903,  0.07071622802089339,  0.12334556714490325, -0.02651441186298387,  0.99200950450550740,
            -0.36957471874357828,  0.92650620200843881,  0.07071622802089339, -0.92097795757081347, -0.37534416687016503,  0.10448137664308726,  0.12334556714490336, -0.02651441186298370,  0.99200950450550740,
            -0.67365669613657164,  0.72510924905369156, -0.14283988479115484, -0.70742184475876502, -0.57674111982491216,  0.40856335403608079,  0.21387143172584838,  0.37627949404651073,  0.90148363992456215,
            -0.70742184475876502, -0.57674111982491205,  0.40856335403608068, -0.67365669613657153,  0.72510924905369156, -0.14283988479115484,  0.21387143172584816,  0.37627949404651051,  0.90148363992456215,
            -0.10448137664308721, -0.37534416687016481,  0.92097795757081358, -0.07071622802089383,  0.92650620200843892,  0.36957471874357795, -0.99200950450550729, -0.02651441186298398, -0.12334556714490341,
            -0.07071622802089389,  0.92650620200843881,  0.36957471874357806, -0.10448137664308726, -0.37534416687016481,  0.92097795757081358, -0.99200950450550707, -0.02651441186298376, -0.12334556714490341,
             0.14283988479115445,  0.72510924905369167,  0.67365669613657142, -0.40856335403608063, -0.57674111982491227,  0.70742184475876535, -0.90148363992456204,  0.37627949404651062, -0.21387143172584866,
            -0.40856335403608057, -0.57674111982491205,  0.70742184475876535,  0.14283988479115456,  0.72510924905369178,  0.67365669613657153, -0.90148363992456226,  0.37627949404651045, -0.21387143172584866,
            -0.99200950450550685, -0.02651441186298389, -0.12334556714490341, -0.10448137664308738, -0.37534416687016492,  0.92097795757081358, -0.07071622802089400,  0.92650620200843892,  0.36957471874357800,
            -0.10448137664308738, -0.37534416687016503,  0.92097795757081358, -0.99200950450550707, -0.02651441186298373, -0.12334556714490363, -0.07071622802089378,  0.92650620200843870,  0.36957471874357811,
            -0.40856335403608052, -0.57674111982491205,  0.70742184475876524, -0.90148363992456204,  0.37627949404651040, -0.21387143172584855,  0.14283988479115434,  0.72510924905369178,  0.67365669613657131,
            -0.90148363992456204,  0.37627949404651034, -0.21387143172584855, -0.40856335403608057, -0.57674111982491205,  0.70742184475876535,  0.14283988479115434,  0.72510924905369178,  0.67365669613657131,
             0.12334556714490313, -0.02651441186298376,  0.99200950450550718, -0.92097795757081324, -0.37534416687016514,  0.10448137664308732, -0.36957471874357800,  0.92650620200843892,  0.07071622802089356,
            -0.92097795757081324, -0.37534416687016492,  0.10448137664308732,  0.12334556714490313, -0.02651441186298373,  0.99200950450550718, -0.36957471874357833,  0.92650620200843870,  0.07071622802089356,
            -0.70742184475876502, -0.57674111982491216,  0.40856335403608074,  0.21387143172584816,  0.37627949404651062,  0.90148363992456204, -0.67365669613657142,  0.72510924905369167, -0.14283988479115473,
             0.21387143172584816,  0.37627949404651051,  0.90148363992456215, -0.70742184475876502, -0.57674111982491216,  0.40856335403608068, -0.67365669613657131,  0.72510924905369167, -0.14283988479115461,
             0.21387143172584827, -0.37627949404651062,  0.90148363992456204, -0.67365669613657153, -0.72510924905369178, -0.14283988479115445, -0.70742184475876524,  0.57674111982491216,  0.40856335403608102,
            -0.67365669613657153, -0.72510924905369145, -0.14283988479115423,  0.21387143172584838, -0.37627949404651062,  0.90148363992456226, -0.70742184475876546,  0.57674111982491194,  0.40856335403608063,
            -0.36957471874357817, -0.92650620200843892,  0.07071622802089378,  0.12334556714490325,  0.02651441186298364,  0.99200950450550729, -0.92097795757081369,  0.37534416687016503,  0.10448137664308743,
             0.12334556714490325,  0.02651441186298356,  0.99200950450550729, -0.36957471874357828, -0.92650620200843881,  0.07071622802089378, -0.92097795757081347,  0.37534416687016503,  0.10448137664308754,
            -0.90148363992456226, -0.37627949404651051, -0.21387143172584816,  0.14283988479115467, -0.72510924905369167,  0.67365669613657153, -0.40856335403608091,  0.57674111982491194,  0.70742184475876524,
             0.14283988479115467, -0.72510924905369190,  0.67365669613657153, -0.90148363992456226, -0.37627949404651045, -0.21387143172584816, -0.40856335403608085,  0.57674111982491216,  0.70742184475876524,
            -0.07071622802089356, -0.92650620200843881,  0.36957471874357828, -0.99200950450550751,  0.02651441186298373, -0.12334556714490313, -0.10448137664308760,  0.37534416687016492,  0.92097795757081347,
            -0.99200950450550729,  0.02651441186298367, -0.12334556714490313, -0.07071622802089361, -0.92650620200843892,  0.36957471874357839, -0.10448137664308771,  0.37534416687016503,  0.92097795757081347,
            -0.10448137664308738,  0.37534416687016470,  0.92097795757081347, -0.07071622802089367, -0.92650620200843892,  0.36957471874357828, -0.99200950450550751,  0.02651441186298395, -0.12334556714490302,
            -0.07071622802089356, -0.92650620200843881,  0.36957471874357839, -0.10448137664308749,  0.37534416687016492,  0.92097795757081347, -0.99200950450550751,  0.02651441186298373, -0.12334556714490325,
             0.14283988479115467, -0.72510924905369167,  0.67365669613657164, -0.40856335403608091,  0.57674111982491227,  0.70742184475876524, -0.90148363992456226, -0.37627949404651062, -0.21387143172584827,
            -0.40856335403608085,  0.57674111982491205,  0.70742184475876524,  0.14283988479115467, -0.72510924905369167,  0.67365669613657175, -0.90148363992456226, -0.37627949404651051, -0.21387143172584816,
            -0.92097795757081369,  0.37534416687016492,  0.10448137664308743, -0.36957471874357817, -0.92650620200843903,  0.07071622802089378,  0.12334556714490325,  0.02651441186298384,  0.99200950450550729,
            -0.36957471874357822, -0.92650620200843881,  0.07071622802089389, -0.92097795757081369,  0.37534416687016492,  0.10448137664308738,  0.12334556714490325,  0.02651441186298373,  0.99200950450550729,
            -0.67365669613657175, -0.72510924905369145, -0.14283988479115434, -0.70742184475876546,  0.57674111982491216,  0.40856335403608079,  0.21387143172584849, -0.37627949404651073,  0.90148363992456204,
            -0.70742184475876546,  0.57674111982491205,  0.40856335403608068, -0.67365669613657153, -0.72510924905369167, -0.14283988479115445,  0.21387143172584849, -0.37627949404651051,  0.90148363992456226,
             0.12334556714490325,  0.02651441186298370,  0.99200950450550729, -0.92097795757081369,  0.37534416687016492,  0.10448137664308749, -0.36957471874357817, -0.92650620200843870,  0.07071622802089389,
            -0.92097795757081369,  0.37534416687016514,  0.10448137664308754,  0.12334556714490325,  0.02651441186298364,  0.99200950450550729, -0.36957471874357822, -0.92650620200843892,  0.07071622802089383,
            -0.70742184475876546,  0.57674111982491205,  0.40856335403608079,  0.21387143172584849, -0.37627949404651051,  0.90148363992456226, -0.67365669613657153, -0.72510924905369167, -0.14283988479115434,
             0.21387143172584827, -0.37627949404651051,  0.90148363992456226, -0.70742184475876546,  0.57674111982491216,  0.40856335403608068, -0.67365669613657131, -0.72510924905369178, -0.14283988479115434,
            -0.99200950450550729,  0.02651441186298378, -0.12334556714490302, -0.10448137664308749,  0.37534416687016503,  0.92097795757081347, -0.07071622802089378, -0.92650620200843892,  0.36957471874357806,
            -0.10448137664308749,  0.37534416687016470,  0.92097795757081324, -0.99200950450550751,  0.02651441186298378, -0.12334556714490313, -0.07071622802089367, -0.92650620200843870,  0.36957471874357845,
            -0.40856335403608091,  0.57674111982491216,  0.70742184475876524, -0.90148363992456226, -0.37627949404651045, -0.21387143172584827,  0.14283988479115456, -0.72510924905369178,  0.67365669613657164,
            -0.90148363992456226, -0.37627949404651040, -0.21387143172584816, -0.40856335403608079,  0.57674111982491205,  0.70742184475876524,  0.14283988479115445, -0.72510924905369178,  0.67365669613657153,
            -0.61476523277246187,  0.78867513459481287,  0.00743240513391585, -0.74472292782631955, -0.57735026918962584, -0.33474531727338752, -0.25971420705077958, -0.21132486540518702,  0.94227932796886926,
            -0.74472292782631944, -0.57735026918962562, -0.33474531727338741, -0.61476523277246187,  0.78867513459481287,  0.00743240513391594, -0.25971420705077963, -0.21132486540518713,  0.94227932796886915,
            -0.81975403804892766, -0.21132486540518708, -0.53230171741593746, -0.46470301232724565,  0.78867513459481287,  0.40254520541901551, -0.33474531727338763, -0.57735026918962573,  0.74472292782631966,
            -0.46470301232724576,  0.78867513459481264,  0.40254520541901539, -0.81975403804892766, -0.21132486540518708, -0.53230171741593746, -0.33474531727338758, -0.57735026918962562,  0.74472292782631966,
            -0.61476523277246165, -0.78867513459481275,  0.00743240513391633, -0.74472292782631966,  0.57735026918962584, -0.33474531727338774, -0.25971420705077974,  0.21132486540518691,  0.94227932796886904,
            -0.74472292782631966,  0.57735026918962573, -0.33474531727338785, -0.61476523277246153, -0.78867513459481275,  0.00743240513391635, -0.25971420705077980,  0.21132486540518702,  0.94227932796886904,
            -0.81975403804892766,  0.21132486540518719, -0.53230171741593746, -0.46470301232724548, -0.78867513459481287,  0.40254520541901573, -0.33474531727338774,  0.57735026918962562,  0.74472292782631944,
            -0.46470301232724548, -0.78867513459481275,  0.40254520541901562, -0.81975403804892766,  0.21132486540518719, -0.53230171741593757, -0.33474531727338780,  0.57735026918962551,  0.74472292782631955,
            -0.33474531727338785,  0.57735026918962562,  0.74472292782631921, -0.81975403804892766,  0.21132486540518730, -0.53230171741593746, -0.46470301232724537, -0.78867513459481287,  0.40254520541901606,
            -0.81975403804892766,  0.21132486540518725, -0.53230171741593757, -0.33474531727338780,  0.57735026918962584,  0.74472292782631933, -0.46470301232724542, -0.78867513459481287,  0.40254520541901584,
            -0.74472292782631966,  0.57735026918962595, -0.33474531727338763, -0.25971420705077963,  0.21132486540518700,  0.94227932796886937, -0.61476523277246165, -0.78867513459481287,  0.00743240513391608,
            -0.25971420705077974,  0.21132486540518700,  0.94227932796886915, -0.74472292782631955,  0.57735026918962595, -0.33474531727338763, -0.61476523277246165, -0.78867513459481287,  0.00743240513391624,
            -0.33474531727338763, -0.57735026918962595,  0.74472292782631944, -0.81975403804892766, -0.21132486540518700, -0.53230171741593757, -0.46470301232724570,  0.78867513459481287,  0.40254520541901573,
            -0.81975403804892755, -0.21132486540518700, -0.53230171741593735, -0.33474531727338774, -0.57735026918962595,  0.74472292782631944, -0.46470301232724565,  0.78867513459481287,  0.40254520541901556,
            -0.74472292782631944, -0.57735026918962562, -0.33474531727338741, -0.25971420705077958, -0.21132486540518730,  0.94227932796886926, -0.61476523277246187,  0.78867513459481287,  0.00743240513391574,
            -0.25971420705077958, -0.21132486540518725,  0.94227932796886937, -0.74472292782631944, -0.57735026918962584, -0.33474531727338752, -0.61476523277246187,  0.78867513459481287,  0.00743240513391596,
            -0.46470301232724565,  0.78867513459481275,  0.40254520541901551, -0.33474531727338769, -0.57735026918962584,  0.74472292782631955, -0.81975403804892755, -0.21132486540518691, -0.53230171741593724,
            -0.33474531727338769, -0.57735026918962573,  0.74472292782631966, -0.46470301232724570,  0.78867513459481275,  0.40254520541901545, -0.81975403804892744, -0.21132486540518702, -0.53230171741593724,
            -0.25971420705077958, -0.21132486540518719,  0.94227932796886926, -0.61476523277246187,  0.78867513459481287,  0.00743240513391610, -0.74472292782631955, -0.57735026918962562, -0.33474531727338763,
            -0.61476523277246187,  0.78867513459481275,  0.00743240513391616, -0.25971420705077963, -0.21132486540518719,  0.94227932796886937, -0.74472292782631944, -0.57735026918962551, -0.33474531727338774,
            -0.46470301232724542, -0.78867513459481287,  0.40254520541901595, -0.33474531727338774,  0.57735026918962584,  0.74472292782631933, -0.81975403804892766,  0.21132486540518702, -0.53230171741593746,
            -0.33474531727338780,  0.57735026918962562,  0.74472292782631921, -0.46470301232724542, -0.78867513459481287,  0.40254520541901584, -0.81975403804892766,  0.21132486540518713, -0.53230171741593735,
            -0.25971420705077963,  0.21132486540518708,  0.94227932796886926, -0.61476523277246165, -0.78867513459481287,  0.00743240513391633, -0.74472292782631966,  0.57735026918962573, -0.33474531727338785,
            -0.61476523277246153, -0.78867513459481264,  0.00743240513391638, -0.25971420705077969,  0.21132486540518708,  0.94227932796886926, -0.74472292782631966,  0.57735026918962562, -0.33474531727338785,
             0.61476523277246187, -0.00743240513391594, -0.78867513459481253,  0.74472292782631944,  0.33474531727338758,  0.57735026918962551,  0.25971420705077963, -0.94227932796886948,  0.21132486540518683,
             0.74472292782631944,  0.33474531727338758,  0.57735026918962540,  0.61476523277246176, -0.00743240513391599, -0.78867513459481242,  0.25971420705077974, -0.94227932796886926,  0.21132486540518700,
             0.81975403804892766,  0.53230171741593746,  0.21132486540518697,  0.46470301232724559, -0.40254520541901562, -0.78867513459481264,  0.33474531727338774, -0.74472292782631966,  0.57735026918962551,
             0.46470301232724570, -0.40254520541901567, -0.78867513459481253,  0.81975403804892744,  0.53230171741593768,  0.21132486540518713,  0.33474531727338774, -0.74472292782631988,  0.57735026918962540,
             0.61476523277246165, -0.00743240513391624,  0.78867513459481242,  0.74472292782631966,  0.33474531727338769, -0.57735026918962551,  0.25971420705077974, -0.94227932796886926, -0.21132486540518694,
             0.74472292782631955,  0.33474531727338780, -0.57735026918962551,  0.61476523277246165, -0.00743240513391630,  0.78867513459481264,  0.25971420705077980, -0.94227932796886926, -0.21132486540518705,
             0.81975403804892766,  0.53230171741593746, -0.21132486540518702,  0.46470301232724553, -0.40254520541901584,  0.78867513459481264,  0.33474531727338780, -0.74472292782631966, -0.57735026918962551,
             0.46470301232724553, -0.40254520541901573,  0.78867513459481264,  0.81975403804892766,  0.53230171741593746, -0.21132486540518711,  0.33474531727338774, -0.74472292782631966, -0.57735026918962540,
             0.33474531727338791, -0.74472292782631944, -0.57735026918962562,  0.81975403804892755,  0.53230171741593746, -0.21132486540518711,  0.46470301232724548, -0.40254520541901595,  0.78867513459481264,
             0.81975403804892766,  0.53230171741593746, -0.21132486540518711,  0.33474531727338785, -0.74472292782631944, -0.57735026918962573,  0.46470301232724548, -0.40254520541901589,  0.78867513459481264,
             0.74472292782631955,  0.33474531727338780, -0.57735026918962584,  0.25971420705077963, -0.94227932796886948, -0.21132486540518702,  0.61476523277246176, -0.00743240513391613,  0.78867513459481264,
             0.25971420705077974, -0.94227932796886926, -0.21132486540518700,  0.74472292782631955,  0.33474531727338769, -0.57735026918962562,  0.61476523277246165, -0.00743240513391624,  0.78867513459481253,
             0.33474531727338774, -0.74472292782631966,  0.57735026918962584,  0.81975403804892766,  0.53230171741593746,  0.21132486540518697,  0.46470301232724553, -0.40254520541901578, -0.78867513459481264,
             0.81975403804892755,  0.53230171741593746,  0.21132486540518694,  0.33474531727338774, -0.74472292782631966,  0.57735026918962562,  0.46470301232724565, -0.40254520541901567, -0.78867513459481253,
             0.74472292782631944,  0.33474531727338747,  0.57735026918962562,  0.25971420705077974, -0.94227932796886948,  0.21132486540518705,  0.61476523277246176, -0.00743240513391596, -0.78867513459481264,
             0.25971420705077963, -0.94227932796886948,  0.21132486540518705,  0.74472292782631944,  0.33474531727338758,  0.57735026918962573,  0.61476523277246187, -0.00743240513391602, -0.78867513459481264,
             0.46470301232724565, -0.40254520541901567, -0.78867513459481242,  0.33474531727338769, -0.74472292782631966,  0.57735026918962551,  0.81975403804892755,  0.53230171741593746,  0.21132486540518688,
             0.33474531727338774, -0.74472292782631966,  0.57735026918962551,  0.46470301232724570, -0.40254520541901562, -0.78867513459481264,  0.81975403804892744,  0.53230171741593746,  0.21132486540518700,
             0.25971420705077969, -0.94227932796886948,  0.21132486540518697,  0.61476523277246176, -0.00743240513391608, -0.78867513459481264,  0.74472292782631944,  0.33474531727338769,  0.57735026918962551,
             0.61476523277246176, -0.00743240513391619, -0.78867513459481264,  0.25971420705077969, -0.94227932796886948,  0.21132486540518705,  0.74472292782631955,  0.33474531727338780,  0.57735026918962540,
             0.46470301232724548, -0.40254520541901595,  0.78867513459481253,  0.33474531727338780, -0.74472292782631944, -0.57735026918962551,  0.81975403804892766,  0.53230171741593746, -0.21132486540518688,
             0.33474531727338785, -0.74472292782631944, -0.57735026918962540,  0.46470301232724553, -0.40254520541901595,  0.78867513459481242,  0.81975403804892755,  0.53230171741593746, -0.21132486540518705,
             0.25971420705077969, -0.94227932796886948, -0.21132486540518702,  0.61476523277246165, -0.00743240513391630,  0.78867513459481264,  0.74472292782631955,  0.33474531727338780, -0.57735026918962551,
             0.61476523277246165, -0.00743240513391624,  0.78867513459481253,  0.25971420705077980, -0.94227932796886971, -0.21132486540518719,  0.74472292782631955,  0.33474531727338791, -0.57735026918962540,
            -0.21387143172584827, -0.90148363992456226,  0.37627949404651062,  0.67365669613657153,  0.14283988479115450,  0.72510924905369178,  0.70742184475876524, -0.40856335403608096, -0.57674111982491216,
             0.67365669613657153,  0.14283988479115417,  0.72510924905369156, -0.21387143172584827, -0.90148363992456226,  0.37627949404651068,  0.70742184475876535, -0.40856335403608068, -0.57674111982491205,
             0.36957471874357817, -0.07071622802089389,  0.92650620200843892, -0.12334556714490325, -0.99200950450550729, -0.02651441186298364,  0.92097795757081347, -0.10448137664308754, -0.37534416687016503,
            -0.12334556714490325, -0.99200950450550729, -0.02651441186298356,  0.36957471874357822, -0.07071622802089383,  0.92650620200843881,  0.92097795757081347, -0.10448137664308760, -0.37534416687016503,
             0.90148363992456226,  0.21387143172584822,  0.37627949404651051, -0.14283988479115456, -0.67365669613657175,  0.72510924905369190,  0.40856335403608085, -0.70742184475876524, -0.57674111982491216,
            -0.14283988479115467, -0.67365669613657153,  0.72510924905369190,  0.90148363992456226,  0.21387143172584810,  0.37627949404651051,  0.40856335403608085, -0.70742184475876524, -0.57674111982491205,
             0.07071622802089361, -0.36957471874357839,  0.92650620200843881,  0.99200950450550751,  0.12334556714490319, -0.02651441186298376,  0.10448137664308754, -0.92097795757081358, -0.37534416687016492,
             0.99200950450550729,  0.12334556714490308, -0.02651441186298370,  0.07071622802089367, -0.36957471874357839,  0.92650620200843892,  0.10448137664308765, -0.92097795757081347, -0.37534416687016503,
             0.10448137664308738, -0.92097795757081347, -0.37534416687016470,  0.07071622802089367, -0.36957471874357839,  0.92650620200843892,  0.99200950450550751,  0.12334556714490308, -0.02651441186298395,
             0.07071622802089367, -0.36957471874357839,  0.92650620200843903,  0.10448137664308749, -0.92097795757081358, -0.37534416687016503,  0.99200950450550751,  0.12334556714490319, -0.02651441186298378,
            -0.14283988479115467, -0.67365669613657175,  0.72510924905369167,  0.40856335403608091, -0.70742184475876524, -0.57674111982491227,  0.90148363992456226,  0.21387143172584822,  0.37627949404651062,
             0.40856335403608085, -0.70742184475876524, -0.57674111982491205, -0.14283988479115456, -0.67365669613657175,  0.72510924905369178,  0.90148363992456226,  0.21387143172584810,  0.37627949404651051,
             0.92097795757081369, -0.10448137664308743, -0.37534416687016492,  0.36957471874357817, -0.07071622802089389,  0.92650620200843903, -0.12334556714490325, -0.99200950450550740, -0.02651441186298378,
             0.36957471874357822, -0.07071622802089395,  0.92650620200843881,  0.92097795757081369, -0.10448137664308743, -0.37534416687016503, -0.12334556714490325, -0.99200950450550729, -0.02651441186298367,
             0.67365669613657175,  0.14283988479115428,  0.72510924905369145,  0.70742184475876546, -0.40856335403608079, -0.57674111982491216, -0.21387143172584849, -0.90148363992456226,  0.37627949404651073,
             0.70742184475876546, -0.40856335403608079, -0.57674111982491227,  0.67365669613657153,  0.14283988479115439,  0.72510924905369178, -0.21387143172584838, -0.90148363992456237,  0.37627949404651062,
            -0.12334556714490325, -0.99200950450550740, -0.02651441186298373,  0.92097795757081369, -0.10448137664308749, -0.37534416687016514,  0.36957471874357822, -0.07071622802089395,  0.92650620200843892,
             0.92097795757081369, -0.10448137664308760, -0.37534416687016514, -0.12334556714490325, -0.99200950450550729, -0.02651441186298367,  0.36957471874357822, -0.07071622802089389,  0.92650620200843881,
             0.70742184475876546, -0.40856335403608079, -0.57674111982491205, -0.21387143172584849, -0.90148363992456237,  0.37627949404651057,  0.67365669613657153,  0.14283988479115439,  0.72510924905369167,
            -0.21387143172584827, -0.90148363992456226,  0.37627949404651051,  0.70742184475876546, -0.40856335403608079, -0.57674111982491216,  0.67365669613657142,  0.14283988479115428,  0.72510924905369178,
             0.99200950450550729,  0.12334556714490308, -0.02651441186298384,  0.10448137664308749, -0.92097795757081369, -0.37534416687016503,  0.07071622802089378, -0.36957471874357822,  0.92650620200843892,
             0.10448137664308749, -0.92097795757081335, -0.37534416687016481,  0.99200950450550729,  0.12334556714490308, -0.02651441186298387,  0.07071622802089372, -0.36957471874357850,  0.92650620200843881,
             0.40856335403608091, -0.70742184475876524, -0.57674111982491216,  0.90148363992456226,  0.21387143172584810,  0.37627949404651045, -0.14283988479115445, -0.67365669613657164,  0.72510924905369178,
             0.90148363992456226,  0.21387143172584810,  0.37627949404651040,  0.40856335403608085, -0.70742184475876535, -0.57674111982491205, -0.14283988479115445, -0.67365669613657153,  0.72510924905369178,
             0.90148363992456237,  0.21387143172584849, -0.37627949404651073, -0.14283988479115439, -0.67365669613657153, -0.72510924905369178,  0.40856335403608085, -0.70742184475876546,  0.57674111982491216,
            -0.14283988479115417, -0.67365669613657153, -0.72510924905369145,  0.90148363992456226,  0.21387143172584849, -0.37627949404651073,  0.40856335403608057, -0.70742184475876546,  0.57674111982491194,
             0.07071622802089395, -0.36957471874357817, -0.92650620200843892,  0.99200950450550729,  0.12334556714490325,  0.02651441186298364,  0.10448137664308754, -0.92097795757081369,  0.37534416687016503,
             0.99200950450550729,  0.12334556714490313,  0.02651441186298356,  0.07071622802089395, -0.36957471874357817, -0.92650620200843881,  0.10448137664308754, -0.92097795757081369,  0.37534416687016503,
            -0.21387143172584833, -0.90148363992456226, -0.37627949404651051,  0.67365669613657186,  0.14283988479115467, -0.72510924905369190,  0.70742184475876524, -0.40856335403608096,  0.57674111982491216,
             0.67365669613657175,  0.14283988479115489, -0.72510924905369190, -0.21387143172584810, -0.90148363992456249, -0.37627949404651045,  0.70742184475876524, -0.40856335403608096,  0.57674111982491216,
             0.36957471874357845, -0.07071622802089356, -0.92650620200843892, -0.12334556714490319, -0.99200950450550751,  0.02651441186298373,  0.92097795757081347, -0.10448137664308760,  0.37534416687016492,
            -0.12334556714490319, -0.99200950450550751,  0.02651441186298373,  0.36957471874357850, -0.07071622802089350, -0.92650620200843892,  0.92097795757081347, -0.10448137664308760,  0.37534416687016492,
             0.92097795757081347, -0.10448137664308743,  0.37534416687016481,  0.36957471874357850, -0.07071622802089350, -0.92650620200843892, -0.12334556714490319, -0.99200950450550751,  0.02651441186298401,
             0.36957471874357850, -0.07071622802089356, -0.92650620200843903,  0.92097795757081358, -0.10448137664308749,  0.37534416687016503, -0.12334556714490319, -0.99200950450550751,  0.02651441186298378,
             0.67365669613657198,  0.14283988479115478, -0.72510924905369178,  0.70742184475876513, -0.40856335403608102,  0.57674111982491227, -0.21387143172584822, -0.90148363992456226, -0.37627949404651062,
             0.70742184475876513, -0.40856335403608091,  0.57674111982491205,  0.67365669613657198,  0.14283988479115478, -0.72510924905369178, -0.21387143172584822, -0.90148363992456249, -0.37627949404651051,
             0.10448137664308726, -0.92097795757081391,  0.37534416687016503,  0.07071622802089406, -0.36957471874357806, -0.92650620200843903,  0.99200950450550740,  0.12334556714490325,  0.02651441186298384,
             0.07071622802089406, -0.36957471874357817, -0.92650620200843881,  0.10448137664308726, -0.92097795757081391,  0.37534416687016503,  0.99200950450550740,  0.12334556714490336,  0.02651441186298367,
            -0.14283988479115428, -0.67365669613657164, -0.72510924905369156,  0.40856335403608068, -0.70742184475876557,  0.57674111982491216,  0.90148363992456237,  0.21387143172584849, -0.37627949404651084,
             0.40856335403608068, -0.70742184475876546,  0.57674111982491227, -0.14283988479115439, -0.67365669613657153, -0.72510924905369178,  0.90148363992456237,  0.21387143172584849, -0.37627949404651062,
             0.99200950450550751,  0.12334556714490325,  0.02651441186298373,  0.10448137664308732, -0.92097795757081369,  0.37534416687016514,  0.07071622802089400, -0.36957471874357811, -0.92650620200843892,
             0.10448137664308743, -0.92097795757081391,  0.37534416687016514,  0.99200950450550729,  0.12334556714490336,  0.02651441186298364,  0.07071622802089400, -0.36957471874357811, -0.92650620200843892,
             0.40856335403608074, -0.70742184475876546,  0.57674111982491216,  0.90148363992456237,  0.21387143172584849, -0.37627949404651051, -0.14283988479115428, -0.67365669613657153, -0.72510924905369167,
             0.90148363992456237,  0.21387143172584849, -0.37627949404651051,  0.40856335403608068, -0.70742184475876557,  0.57674111982491216, -0.14283988479115428, -0.67365669613657153, -0.72510924905369167,
            -0.12334556714490319, -0.99200950450550751,  0.02651441186298395,  0.92097795757081358, -0.10448137664308749,  0.37534416687016503,  0.36957471874357833, -0.07071622802089361, -0.92650620200843892,
             0.92097795757081335, -0.10448137664308760,  0.37534416687016470, -0.12334556714490308, -0.99200950450550751,  0.02651441186298389,  0.36957471874357861, -0.07071622802089356, -0.92650620200843870,
             0.70742184475876524, -0.40856335403608091,  0.57674111982491216, -0.21387143172584810, -0.90148363992456226, -0.37627949404651045,  0.67365669613657164,  0.14283988479115456, -0.72510924905369178,
            -0.21387143172584810, -0.90148363992456226, -0.37627949404651040,  0.70742184475876524, -0.40856335403608091,  0.57674111982491205,  0.67365669613657164,  0.14283988479115456, -0.72510924905369178,
             0.00743240513391594, -0.61476523277246187,  0.78867513459481287, -0.33474531727338758, -0.74472292782631944, -0.57735026918962584,  0.94227932796886948, -0.25971420705077963, -0.21132486540518694,
            -0.33474531727338758, -0.74472292782631944, -0.57735026918962584,  0.00743240513391605, -0.61476523277246187,  0.78867513459481287,  0.94227932796886926, -0.25971420705077974, -0.21132486540518700,
            -0.53230171741593746, -0.81975403804892766, -0.21132486540518711,  0.40254520541901562, -0.46470301232724570,  0.78867513459481287,  0.74472292782631966, -0.33474531727338763, -0.57735026918962573,
             0.40254520541901556, -0.46470301232724581,  0.78867513459481275, -0.53230171741593746, -0.81975403804892766, -0.21132486540518711,  0.74472292782631988, -0.33474531727338752, -0.57735026918962551,
             0.00743240513391624, -0.61476523277246153, -0.78867513459481298, -0.33474531727338769, -0.74472292782631955,  0.57735026918962573,  0.94227932796886926, -0.25971420705077986,  0.21132486540518716,
            -0.33474531727338780, -0.74472292782631966,  0.57735026918962573,  0.00743240513391630, -0.61476523277246153, -0.78867513459481298,  0.94227932796886926, -0.25971420705077986,  0.21132486540518722,
            -0.53230171741593746, -0.81975403804892766,  0.21132486540518711,  0.40254520541901584, -0.46470301232724548, -0.78867513459481287,  0.74472292782631966, -0.33474531727338780,  0.57735026918962573,
             0.40254520541901573, -0.46470301232724548, -0.78867513459481287, -0.53230171741593746, -0.81975403804892766,  0.21132486540518716,  0.74472292782631966, -0.33474531727338780,  0.57735026918962573,
             0.74472292782631944, -0.33474531727338797,  0.57735026918962584, -0.53230171741593768, -0.81975403804892766,  0.21132486540518716,  0.40254520541901601, -0.46470301232724531, -0.78867513459481287,
            -0.53230171741593746, -0.81975403804892766,  0.21132486540518713,  0.74472292782631944, -0.33474531727338785,  0.57735026918962584,  0.40254520541901595, -0.46470301232724548, -0.78867513459481298,
            -0.33474531727338758, -0.74472292782631966,  0.57735026918962595,  0.94227932796886948, -0.25971420705077974,  0.21132486540518702,  0.00743240513391608, -0.61476523277246165, -0.78867513459481287,
             0.94227932796886926, -0.25971420705077986,  0.21132486540518716, -0.33474531727338769, -0.74472292782631955,  0.57735026918962584,  0.00743240513391619, -0.61476523277246153, -0.78867513459481287,
             0.74472292782631944, -0.33474531727338769, -0.57735026918962595, -0.53230171741593746, -0.81975403804892755, -0.21132486540518697,  0.40254520541901584, -0.46470301232724570,  0.78867513459481287,
            -0.53230171741593746, -0.81975403804892744, -0.21132486540518711,  0.74472292782631966, -0.33474531727338774, -0.57735026918962584,  0.40254520541901573, -0.46470301232724570,  0.78867513459481287,
            -0.33474531727338747, -0.74472292782631933, -0.57735026918962584,  0.94227932796886948, -0.25971420705077963, -0.21132486540518711,  0.00743240513391591, -0.61476523277246198,  0.78867513459481287,
             0.94227932796886948, -0.25971420705077969, -0.21132486540518708, -0.33474531727338758, -0.74472292782631944, -0.57735026918962584,  0.00743240513391596, -0.61476523277246176,  0.78867513459481298,
             0.40254520541901567, -0.46470301232724570,  0.78867513459481298,  0.74472292782631966, -0.33474531727338774, -0.57735026918962573, -0.53230171741593746, -0.81975403804892744, -0.21132486540518711,
             0.74472292782631966, -0.33474531727338763, -0.57735026918962573,  0.40254520541901562, -0.46470301232724576,  0.78867513459481298, -0.53230171741593746, -0.81975403804892744, -0.21132486540518716,
             0.94227932796886948, -0.25971420705077963, -0.21132486540518705,  0.00743240513391610, -0.61476523277246176,  0.78867513459481287, -0.33474531727338769, -0.74472292782631944, -0.57735026918962573,
             0.00743240513391619, -0.61476523277246187,  0.78867513459481287,  0.94227932796886948, -0.25971420705077963, -0.21132486540518711, -0.33474531727338780, -0.74472292782631944, -0.57735026918962573,
             0.40254520541901595, -0.46470301232724548, -0.78867513459481287,  0.74472292782631944, -0.33474531727338780,  0.57735026918962584, -0.53230171741593746, -0.81975403804892766,  0.21132486540518700,
             0.74472292782631944, -0.33474531727338785,  0.57735026918962584,  0.40254520541901584, -0.46470301232724542, -0.78867513459481287, -0.53230171741593746, -0.81975403804892755,  0.21132486540518705,
             0.94227932796886948, -0.25971420705077969,  0.21132486540518716,  0.00743240513391627, -0.61476523277246153, -0.78867513459481287, -0.33474531727338780, -0.74472292782631966,  0.57735026918962573,
             0.00743240513391635, -0.61476523277246153, -0.78867513459481275,  0.94227932796886926, -0.25971420705077969,  0.21132486540518716, -0.33474531727338791, -0.74472292782631977,  0.57735026918962551,
        });
    }

    private static void populateFrames3D_Part3() {
        populateFrames3D(5184, new double[] {
            -0.78867513459481287, -0.00743240513391580,  0.61476523277246176,  0.57735026918962562,  0.33474531727338763,  0.74472292782631944,  0.21132486540518713, -0.94227932796886948,  0.25971420705077936,
             0.57735026918962551,  0.33474531727338763,  0.74472292782631944, -0.78867513459481287, -0.00743240513391588,  0.61476523277246176,  0.21132486540518741, -0.94227932796886937,  0.25971420705077941,
             0.21132486540518691,  0.53230171741593768,  0.81975403804892766, -0.78867513459481287, -0.40254520541901551,  0.46470301232724553,  0.57735026918962584, -0.74472292782631988,  0.33474531727338741,
            -0.78867513459481287, -0.40254520541901539,  0.46470301232724553,  0.21132486540518700,  0.53230171741593768,  0.81975403804892766,  0.57735026918962584, -0.74472292782631988,  0.33474531727338741,
             0.78867513459481275, -0.00743240513391630,  0.61476523277246153, -0.57735026918962595,  0.33474531727338808,  0.74472292782631966, -0.21132486540518669, -0.94227932796886937,  0.25971420705077941,
            -0.57735026918962595,  0.33474531727338808,  0.74472292782631966,  0.78867513459481264, -0.00743240513391630,  0.61476523277246153, -0.21132486540518680, -0.94227932796886926,  0.25971420705077947,
            -0.21132486540518738,  0.53230171741593768,  0.81975403804892766,  0.78867513459481320, -0.40254520541901584,  0.46470301232724531, -0.57735026918962562, -0.74472292782631944,  0.33474531727338763,
             0.78867513459481298, -0.40254520541901578,  0.46470301232724531, -0.21132486540518744,  0.53230171741593779,  0.81975403804892766, -0.57735026918962551, -0.74472292782631966,  0.33474531727338752,
            -0.57735026918962551, -0.74472292782631933,  0.33474531727338763, -0.21132486540518747,  0.53230171741593768,  0.81975403804892766,  0.78867513459481287, -0.40254520541901612,  0.46470301232724526,
            -0.21132486540518736,  0.53230171741593768,  0.81975403804892766, -0.57735026918962562, -0.74472292782631944,  0.33474531727338763,  0.78867513459481287, -0.40254520541901589,  0.46470301232724531,
            -0.57735026918962595,  0.33474531727338785,  0.74472292782631966, -0.21132486540518688, -0.94227932796886948,  0.25971420705077941,  0.78867513459481287, -0.00743240513391608,  0.61476523277246153,
            -0.21132486540518680, -0.94227932796886948,  0.25971420705077947, -0.57735026918962606,  0.33474531727338797,  0.74472292782631966,  0.78867513459481287, -0.00743240513391624,  0.61476523277246153,
             0.57735026918962595, -0.74472292782631966,  0.33474531727338747,  0.21132486540518688,  0.53230171741593768,  0.81975403804892766, -0.78867513459481287, -0.40254520541901573,  0.46470301232724548,
             0.21132486540518680,  0.53230171741593768,  0.81975403804892766,  0.57735026918962606, -0.74472292782631977,  0.33474531727338741, -0.78867513459481287, -0.40254520541901556,  0.46470301232724553,
             0.57735026918962551,  0.33474531727338752,  0.74472292782631944,  0.21132486540518747, -0.94227932796886948,  0.25971420705077936, -0.78867513459481287, -0.00743240513391569,  0.61476523277246176,
             0.21132486540518736, -0.94227932796886948,  0.25971420705077936,  0.57735026918962562,  0.33474531727338763,  0.74472292782631944, -0.78867513459481287, -0.00743240513391591,  0.61476523277246176,
            -0.78867513459481275, -0.40254520541901551,  0.46470301232724553,  0.57735026918962595, -0.74472292782631988,  0.33474531727338741,  0.21132486540518669,  0.53230171741593757,  0.81975403804892766,
             0.57735026918962595, -0.74472292782631988,  0.33474531727338741, -0.78867513459481264, -0.40254520541901551,  0.46470301232724553,  0.21132486540518680,  0.53230171741593746,  0.81975403804892766,
             0.21132486540518738, -0.94227932796886948,  0.25971420705077941, -0.78867513459481320, -0.00743240513391596,  0.61476523277246176,  0.57735026918962562,  0.33474531727338763,  0.74472292782631944,
            -0.78867513459481298, -0.00743240513391602,  0.61476523277246176,  0.21132486540518744, -0.94227932796886960,  0.25971420705077936,  0.57735026918962551,  0.33474531727338785,  0.74472292782631955,
             0.78867513459481287, -0.40254520541901601,  0.46470301232724526, -0.57735026918962562, -0.74472292782631944,  0.33474531727338758, -0.21132486540518713,  0.53230171741593768,  0.81975403804892766,
            -0.57735026918962551, -0.74472292782631944,  0.33474531727338763,  0.78867513459481287, -0.40254520541901595,  0.46470301232724531, -0.21132486540518741,  0.53230171741593757,  0.81975403804892766,
            -0.21132486540518691, -0.94227932796886948,  0.25971420705077941,  0.78867513459481287, -0.00743240513391633,  0.61476523277246153, -0.57735026918962584,  0.33474531727338808,  0.74472292782631966,
             0.78867513459481287, -0.00743240513391638,  0.61476523277246153, -0.21132486540518700, -0.94227932796886948,  0.25971420705077947, -0.57735026918962584,  0.33474531727338808,  0.74472292782631966,
             0.37627949404651062, -0.90148363992456226, -0.21387143172584844,  0.72510924905369156,  0.14283988479115456,  0.67365669613657153, -0.57674111982491205, -0.40856335403608068,  0.70742184475876502,
             0.72510924905369145,  0.14283988479115445,  0.67365669613657142,  0.37627949404651062, -0.90148363992456204, -0.21387143172584844, -0.57674111982491183, -0.40856335403608063,  0.70742184475876524,
             0.92650620200843870, -0.07071622802089378,  0.36957471874357817, -0.02651441186298359, -0.99200950450550707, -0.12334556714490341, -0.37534416687016492, -0.10448137664308738,  0.92097795757081335,
            -0.02651441186298356, -0.99200950450550707, -0.12334556714490352,  0.92650620200843881, -0.07071622802089367,  0.36957471874357817, -0.37534416687016514, -0.10448137664308743,  0.92097795757081335,
             0.37627949404651040,  0.21387143172584827,  0.90148363992456215,  0.72510924905369167, -0.67365669613657153, -0.14283988479115484, -0.57674111982491194, -0.70742184475876502,  0.40856335403608074,
             0.72510924905369167, -0.67365669613657153, -0.14283988479115461,  0.37627949404651040,  0.21387143172584838,  0.90148363992456226, -0.57674111982491183, -0.70742184475876513,  0.40856335403608052,
             0.92650620200843870, -0.36957471874357822,  0.07071622802089356, -0.02651441186298370,  0.12334556714490336,  0.99200950450550718, -0.37534416687016481, -0.92097795757081347,  0.10448137664308732,
            -0.02651441186298367,  0.12334556714490336,  0.99200950450550729,  0.92650620200843870, -0.36957471874357828,  0.07071622802089350, -0.37534416687016492, -0.92097795757081324,  0.10448137664308738,
            -0.37534416687016459, -0.92097795757081347,  0.10448137664308721,  0.92650620200843881, -0.36957471874357828,  0.07071622802089350, -0.02651441186298398,  0.12334556714490336,  0.99200950450550740,
             0.92650620200843870, -0.36957471874357828,  0.07071622802089339, -0.37534416687016470, -0.92097795757081347,  0.10448137664308726, -0.02651441186298387,  0.12334556714490336,  0.99200950450550740,
             0.72510924905369167, -0.67365669613657153, -0.14283988479115484, -0.57674111982491205, -0.70742184475876513,  0.40856335403608068,  0.37627949404651051,  0.21387143172584838,  0.90148363992456215,
            -0.57674111982491194, -0.70742184475876502,  0.40856335403608068,  0.72510924905369167, -0.67365669613657164, -0.14283988479115484,  0.37627949404651045,  0.21387143172584838,  0.90148363992456215,
            -0.37534416687016492, -0.10448137664308726,  0.92097795757081358,  0.92650620200843881, -0.07071622802089372,  0.36957471874357806, -0.02651441186298376, -0.99200950450550729, -0.12334556714490341,
             0.92650620200843870, -0.07071622802089383,  0.36957471874357806, -0.37534416687016492, -0.10448137664308721,  0.92097795757081358, -0.02651441186298370, -0.99200950450550729, -0.12334556714490341,
             0.72510924905369134,  0.14283988479115456,  0.67365669613657153, -0.57674111982491205, -0.40856335403608057,  0.70742184475876524,  0.37627949404651073, -0.90148363992456226, -0.21387143172584866,
            -0.57674111982491194, -0.40856335403608057,  0.70742184475876535,  0.72510924905369145,  0.14283988479115456,  0.67365669613657153,  0.37627949404651062, -0.90148363992456226, -0.21387143172584866,
            -0.02651441186298362, -0.99200950450550707, -0.12334556714490341, -0.37534416687016492, -0.10448137664308726,  0.92097795757081358,  0.92650620200843870, -0.07071622802089383,  0.36957471874357800,
            -0.37534416687016492, -0.10448137664308726,  0.92097795757081335, -0.02651441186298362, -0.99200950450550729, -0.12334556714490352,  0.92650620200843858, -0.07071622802089372,  0.36957471874357822,
            -0.57674111982491194, -0.40856335403608063,  0.70742184475876524,  0.37627949404651045, -0.90148363992456226, -0.21387143172584844,  0.72510924905369156,  0.14283988479115456,  0.67365669613657142,
             0.37627949404651040, -0.90148363992456226, -0.21387143172584855, -0.57674111982491194, -0.40856335403608057,  0.70742184475876524,  0.72510924905369167,  0.14283988479115445,  0.67365669613657131,
            -0.02651441186298389,  0.12334556714490336,  0.99200950450550718, -0.37534416687016481, -0.92097795757081347,  0.10448137664308721,  0.92650620200843881, -0.36957471874357817,  0.07071622802089367,
            -0.37534416687016470, -0.92097795757081324,  0.10448137664308732, -0.02651441186298384,  0.12334556714490325,  0.99200950450550718,  0.92650620200843858, -0.36957471874357822,  0.07071622802089356,
            -0.57674111982491194, -0.70742184475876502,  0.40856335403608057,  0.37627949404651034,  0.21387143172584827,  0.90148363992456215,  0.72510924905369167, -0.67365669613657153, -0.14283988479115461,
             0.37627949404651029,  0.21387143172584827,  0.90148363992456226, -0.57674111982491205, -0.70742184475876524,  0.40856335403608057,  0.72510924905369190, -0.67365669613657142, -0.14283988479115461,
            -0.37627949404651062,  0.21387143172584827,  0.90148363992456226, -0.72510924905369156, -0.67365669613657175, -0.14283988479115445,  0.57674111982491205, -0.70742184475876524,  0.40856335403608079,
            -0.72510924905369156, -0.67365669613657164, -0.14283988479115423, -0.37627949404651062,  0.21387143172584827,  0.90148363992456226,  0.57674111982491205, -0.70742184475876524,  0.40856335403608063,
            -0.92650620200843881, -0.36957471874357828,  0.07071622802089389,  0.02651441186298359,  0.12334556714490325,  0.99200950450550729,  0.37534416687016503, -0.92097795757081347,  0.10448137664308749,
             0.02651441186298353,  0.12334556714490336,  0.99200950450550707, -0.92650620200843881, -0.36957471874357828,  0.07071622802089378,  0.37534416687016514, -0.92097795757081347,  0.10448137664308760,
            -0.37627949404651051, -0.90148363992456226, -0.21387143172584816, -0.72510924905369178,  0.14283988479115456,  0.67365669613657175,  0.57674111982491194, -0.40856335403608079,  0.70742184475876513,
            -0.72510924905369178,  0.14283988479115456,  0.67365669613657164, -0.37627949404651040, -0.90148363992456249, -0.21387143172584827,  0.57674111982491194, -0.40856335403608074,  0.70742184475876524,
            -0.92650620200843870, -0.07071622802089367,  0.36957471874357833,  0.02651441186298367, -0.99200950450550729, -0.12334556714490325,  0.37534416687016481, -0.10448137664308754,  0.92097795757081347,
             0.02651441186298364, -0.99200950450550751, -0.12334556714490325, -0.92650620200843870, -0.07071622802089367,  0.36957471874357839,  0.37534416687016492, -0.10448137664308754,  0.92097795757081347,
             0.37534416687016470, -0.10448137664308738,  0.92097795757081347, -0.92650620200843881, -0.07071622802089367,  0.36957471874357839,  0.02651441186298395, -0.99200950450550751, -0.12334556714490313,
            -0.92650620200843881, -0.07071622802089367,  0.36957471874357839,  0.37534416687016481, -0.10448137664308726,  0.92097795757081347,  0.02651441186298389, -0.99200950450550751, -0.12334556714490325,
            -0.72510924905369167,  0.14283988479115467,  0.67365669613657175,  0.57674111982491216, -0.40856335403608079,  0.70742184475876524, -0.37627949404651062, -0.90148363992456249, -0.21387143172584827,
             0.57674111982491205, -0.40856335403608079,  0.70742184475876513, -0.72510924905369167,  0.14283988479115467,  0.67365669613657175, -0.37627949404651051, -0.90148363992456249, -0.21387143172584827,
             0.37534416687016492, -0.92097795757081369,  0.10448137664308738, -0.92650620200843892, -0.36957471874357828,  0.07071622802089383,  0.02651441186298378,  0.12334556714490336,  0.99200950450550729,
            -0.92650620200843881, -0.36957471874357828,  0.07071622802089395,  0.37534416687016492, -0.92097795757081369,  0.10448137664308732,  0.02651441186298367,  0.12334556714490336,  0.99200950450550729,
            -0.72510924905369145, -0.67365669613657175, -0.14283988479115434,  0.57674111982491205, -0.70742184475876546,  0.40856335403608068, -0.37627949404651073,  0.21387143172584849,  0.90148363992456226,
             0.57674111982491205, -0.70742184475876546,  0.40856335403608068, -0.72510924905369156, -0.67365669613657175, -0.14283988479115445, -0.37627949404651073,  0.21387143172584849,  0.90148363992456226,
             0.02651441186298370,  0.12334556714490325,  0.99200950450550729,  0.37534416687016503, -0.92097795757081369,  0.10448137664308738, -0.92650620200843870, -0.36957471874357828,  0.07071622802089395,
             0.37534416687016503, -0.92097795757081369,  0.10448137664308743,  0.02651441186298359,  0.12334556714490336,  0.99200950450550729, -0.92650620200843870, -0.36957471874357833,  0.07071622802089389,
             0.57674111982491194, -0.70742184475876546,  0.40856335403608074, -0.37627949404651051,  0.21387143172584827,  0.90148363992456226, -0.72510924905369156, -0.67365669613657153, -0.14283988479115445,
            -0.37627949404651045,  0.21387143172584838,  0.90148363992456226,  0.57674111982491194, -0.70742184475876546,  0.40856335403608068, -0.72510924905369167, -0.67365669613657153, -0.14283988479115434,
             0.02651441186298381, -0.99200950450550729, -0.12334556714490325,  0.37534416687016481, -0.10448137664308738,  0.92097795757081347, -0.92650620200843881, -0.07071622802089383,  0.36957471874357828,
             0.37534416687016481, -0.10448137664308743,  0.92097795757081324,  0.02651441186298381, -0.99200950450550729, -0.12334556714490313, -0.92650620200843881, -0.07071622802089378,  0.36957471874357845,
             0.57674111982491205, -0.40856335403608079,  0.70742184475876524, -0.37627949404651040, -0.90148363992456226, -0.21387143172584816, -0.72510924905369178,  0.14283988479115445,  0.67365669613657153,
            -0.37627949404651034, -0.90148363992456249, -0.21387143172584805,  0.57674111982491205, -0.40856335403608079,  0.70742184475876524, -0.72510924905369190,  0.14283988479115445,  0.67365669613657153,
             0.78867513459481287, -0.61476523277246187,  0.00743240513391580, -0.57735026918962562, -0.74472292782631966, -0.33474531727338763, -0.21132486540518713, -0.25971420705077947,  0.94227932796886948,
            -0.57735026918962551, -0.74472292782631966, -0.33474531727338774,  0.78867513459481287, -0.61476523277246187,  0.00743240513391585, -0.21132486540518741, -0.25971420705077947,  0.94227932796886948,
            -0.21132486540518691, -0.81975403804892766, -0.53230171741593768,  0.78867513459481287, -0.46470301232724565,  0.40254520541901551, -0.57735026918962584, -0.33474531727338752,  0.74472292782631988,
             0.78867513459481275, -0.46470301232724570,  0.40254520541901539, -0.21132486540518697, -0.81975403804892766, -0.53230171741593768, -0.57735026918962573, -0.33474531727338752,  0.74472292782631988,
            -0.78867513459481275, -0.61476523277246165,  0.00743240513391630,  0.57735026918962595, -0.74472292782631988, -0.33474531727338808,  0.21132486540518669, -0.25971420705077952,  0.94227932796886926,
             0.57735026918962595, -0.74472292782631988, -0.33474531727338808, -0.78867513459481275, -0.61476523277246165,  0.00743240513391641,  0.21132486540518686, -0.25971420705077952,  0.94227932796886926,
             0.21132486540518741, -0.81975403804892777, -0.53230171741593768, -0.78867513459481320, -0.46470301232724542,  0.40254520541901584,  0.57735026918962562, -0.33474531727338769,  0.74472292782631944,
            -0.78867513459481298, -0.46470301232724542,  0.40254520541901578,  0.21132486540518744, -0.81975403804892788, -0.53230171741593779,  0.57735026918962551, -0.33474531727338763,  0.74472292782631966,
             0.57735026918962562, -0.33474531727338774,  0.74472292782631933,  0.21132486540518747, -0.81975403804892777, -0.53230171741593768, -0.78867513459481287, -0.46470301232724537,  0.40254520541901612,
             0.21132486540518736, -0.81975403804892788, -0.53230171741593768,  0.57735026918962562, -0.33474531727338774,  0.74472292782631933, -0.78867513459481287, -0.46470301232724537,  0.40254520541901595,
             0.57735026918962606, -0.74472292782631966, -0.33474531727338797,  0.21132486540518680, -0.25971420705077952,  0.94227932796886948, -0.78867513459481287, -0.61476523277246165,  0.00743240513391608,
             0.21132486540518680, -0.25971420705077952,  0.94227932796886948,  0.57735026918962606, -0.74472292782631966, -0.33474531727338797, -0.78867513459481287, -0.61476523277246165,  0.00743240513391624,
            -0.57735026918962606, -0.33474531727338758,  0.74472292782631977, -0.21132486540518680, -0.81975403804892777, -0.53230171741593768,  0.78867513459481287, -0.46470301232724559,  0.40254520541901573,
            -0.21132486540518680, -0.81975403804892777, -0.53230171741593768, -0.57735026918962606, -0.33474531727338758,  0.74472292782631977,  0.78867513459481287, -0.46470301232724559,  0.40254520541901556,
            -0.57735026918962562, -0.74472292782631955, -0.33474531727338752, -0.21132486540518747, -0.25971420705077952,  0.94227932796886948,  0.78867513459481287, -0.61476523277246187,  0.00743240513391569,
            -0.21132486540518736, -0.25971420705077947,  0.94227932796886948, -0.57735026918962562, -0.74472292782631955, -0.33474531727338752,  0.78867513459481287, -0.61476523277246198,  0.00743240513391585,
             0.78867513459481275, -0.46470301232724565,  0.40254520541901551, -0.57735026918962595, -0.33474531727338747,  0.74472292782631988, -0.21132486540518669, -0.81975403804892777, -0.53230171741593746,
            -0.57735026918962595, -0.33474531727338747,  0.74472292782631988,  0.78867513459481275, -0.46470301232724570,  0.40254520541901539, -0.21132486540518686, -0.81975403804892777, -0.53230171741593746,
            -0.21132486540518741, -0.25971420705077952,  0.94227932796886948,  0.78867513459481320, -0.61476523277246187,  0.00743240513391596, -0.57735026918962562, -0.74472292782631966, -0.33474531727338763,
             0.78867513459481298, -0.61476523277246187,  0.00743240513391602, -0.21132486540518744, -0.25971420705077947,  0.94227932796886960, -0.57735026918962551, -0.74472292782631966, -0.33474531727338785,
            -0.78867513459481287, -0.46470301232724537,  0.40254520541901601,  0.57735026918962562, -0.33474531727338769,  0.74472292782631944,  0.21132486540518713, -0.81975403804892788, -0.53230171741593768,
             0.57735026918962551, -0.33474531727338769,  0.74472292782631955, -0.78867513459481287, -0.46470301232724542,  0.40254520541901595,  0.21132486540518741, -0.81975403804892788, -0.53230171741593768,
             0.21132486540518691, -0.25971420705077958,  0.94227932796886948, -0.78867513459481287, -0.61476523277246165,  0.00743240513391633,  0.57735026918962584, -0.74472292782631977, -0.33474531727338808,
            -0.78867513459481275, -0.61476523277246153,  0.00743240513391638,  0.21132486540518697, -0.25971420705077958,  0.94227932796886948,  0.57735026918962573, -0.74472292782631977, -0.33474531727338808,
             0.61476523277246187, -0.78867513459481287, -0.00743240513391585,  0.74472292782631955,  0.57735026918962584,  0.33474531727338752,  0.25971420705077958,  0.21132486540518702, -0.94227932796886926,
             0.74472292782631944,  0.57735026918962562,  0.33474531727338741,  0.61476523277246187, -0.78867513459481287, -0.00743240513391594,  0.25971420705077963,  0.21132486540518713, -0.94227932796886915,
             0.81975403804892766,  0.21132486540518708,  0.53230171741593746,  0.46470301232724565, -0.78867513459481287, -0.40254520541901551,  0.33474531727338763,  0.57735026918962573, -0.74472292782631966,
             0.46470301232724576, -0.78867513459481264, -0.40254520541901539,  0.81975403804892766,  0.21132486540518708,  0.53230171741593746,  0.33474531727338758,  0.57735026918962562, -0.74472292782631966,
             0.61476523277246165,  0.78867513459481275, -0.00743240513391633,  0.74472292782631966, -0.57735026918962584,  0.33474531727338774,  0.25971420705077974, -0.21132486540518691, -0.94227932796886904,
             0.74472292782631966, -0.57735026918962573,  0.33474531727338785,  0.61476523277246153,  0.78867513459481275, -0.00743240513391635,  0.25971420705077980, -0.21132486540518702, -0.94227932796886904,
             0.81975403804892766, -0.21132486540518719,  0.53230171741593746,  0.46470301232724548,  0.78867513459481287, -0.40254520541901573,  0.33474531727338774, -0.57735026918962562, -0.74472292782631944,
             0.46470301232724548,  0.78867513459481275, -0.40254520541901562,  0.81975403804892766, -0.21132486540518719,  0.53230171741593757,  0.33474531727338780, -0.57735026918962551, -0.74472292782631955,
             0.33474531727338785, -0.57735026918962562, -0.74472292782631921,  0.81975403804892766, -0.21132486540518730,  0.53230171741593746,  0.46470301232724537,  0.78867513459481287, -0.40254520541901606,
             0.81975403804892766, -0.21132486540518725,  0.53230171741593757,  0.33474531727338780, -0.57735026918962584, -0.74472292782631933,  0.46470301232724542,  0.78867513459481287, -0.40254520541901584,
             0.74472292782631966, -0.57735026918962595,  0.33474531727338763,  0.25971420705077963, -0.21132486540518700, -0.94227932796886937,  0.61476523277246165,  0.78867513459481287, -0.00743240513391608,
             0.25971420705077974, -0.21132486540518700, -0.94227932796886915,  0.74472292782631955, -0.57735026918962595,  0.33474531727338763,  0.61476523277246165,  0.78867513459481287, -0.00743240513391624,
             0.33474531727338763,  0.57735026918962595, -0.74472292782631944,  0.81975403804892766,  0.21132486540518700,  0.53230171741593757,  0.46470301232724570, -0.78867513459481287, -0.40254520541901573,
             0.81975403804892755,  0.21132486540518700,  0.53230171741593735,  0.33474531727338774,  0.57735026918962595, -0.74472292782631944,  0.46470301232724565, -0.78867513459481287, -0.40254520541901556,
             0.74472292782631944,  0.57735026918962562,  0.33474531727338741,  0.25971420705077958,  0.21132486540518730, -0.94227932796886926,  0.61476523277246187, -0.78867513459481287, -0.00743240513391574,
             0.25971420705077958,  0.21132486540518725, -0.94227932796886937,  0.74472292782631944,  0.57735026918962584,  0.33474531727338752,  0.61476523277246187, -0.78867513459481287, -0.00743240513391596,
             0.46470301232724565, -0.78867513459481275, -0.40254520541901551,  0.33474531727338769,  0.57735026918962584, -0.74472292782631955,  0.81975403804892755,  0.21132486540518691,  0.53230171741593724,
             0.33474531727338769,  0.57735026918962573, -0.74472292782631966,  0.46470301232724570, -0.78867513459481275, -0.40254520541901545,  0.81975403804892744,  0.21132486540518702,  0.53230171741593724,
             0.25971420705077958,  0.21132486540518719, -0.94227932796886926,  0.61476523277246187, -0.78867513459481287, -0.00743240513391610,  0.74472292782631955,  0.57735026918962562,  0.33474531727338763,
             0.61476523277246187, -0.78867513459481275, -0.00743240513391616,  0.25971420705077963,  0.21132486540518719, -0.94227932796886937,  0.74472292782631944,  0.57735026918962551,  0.33474531727338774,
             0.46470301232724542,  0.78867513459481287, -0.40254520541901595,  0.33474531727338774, -0.57735026918962584, -0.74472292782631933,  0.81975403804892766, -0.21132486540518702,  0.53230171741593746,
             0.33474531727338780, -0.57735026918962562, -0.74472292782631921,  0.46470301232724542,  0.78867513459481287, -0.40254520541901584,  0.81975403804892766, -0.21132486540518713,  0.53230171741593735,
             0.25971420705077963, -0.21132486540518708, -0.94227932796886926,  0.61476523277246165,  0.78867513459481287, -0.00743240513391633,  0.74472292782631966, -0.57735026918962573,  0.33474531727338785,
             0.61476523277246153,  0.78867513459481264, -0.00743240513391638,  0.25971420705077969, -0.21132486540518708, -0.94227932796886926,  0.74472292782631966, -0.57735026918962562,  0.33474531727338785,
            -0.21387143172584827,  0.37627949404651062, -0.90148363992456204,  0.67365669613657153,  0.72510924905369178,  0.14283988479115445,  0.70742184475876524, -0.57674111982491216, -0.40856335403608102,
             0.67365669613657153,  0.72510924905369145,  0.14283988479115423, -0.21387143172584838,  0.37627949404651062, -0.90148363992456226,  0.70742184475876546, -0.57674111982491194, -0.40856335403608063,
             0.36957471874357817,  0.92650620200843892, -0.07071622802089378, -0.12334556714490325, -0.02651441186298364, -0.99200950450550729,  0.92097795757081369, -0.37534416687016503, -0.10448137664308743,
            -0.12334556714490325, -0.02651441186298356, -0.99200950450550729,  0.36957471874357828,  0.92650620200843881, -0.07071622802089378,  0.92097795757081347, -0.37534416687016503, -0.10448137664308754,
             0.90148363992456226,  0.37627949404651051,  0.21387143172584816, -0.14283988479115467,  0.72510924905369167, -0.67365669613657153,  0.40856335403608091, -0.57674111982491194, -0.70742184475876524,
            -0.14283988479115467,  0.72510924905369190, -0.67365669613657153,  0.90148363992456226,  0.37627949404651045,  0.21387143172584816,  0.40856335403608085, -0.57674111982491216, -0.70742184475876524,
             0.07071622802089356,  0.92650620200843881, -0.36957471874357828,  0.99200950450550751, -0.02651441186298373,  0.12334556714490313,  0.10448137664308760, -0.37534416687016492, -0.92097795757081347,
             0.99200950450550729, -0.02651441186298367,  0.12334556714490313,  0.07071622802089361,  0.92650620200843892, -0.36957471874357839,  0.10448137664308771, -0.37534416687016503, -0.92097795757081347,
             0.10448137664308738, -0.37534416687016470, -0.92097795757081347,  0.07071622802089367,  0.92650620200843892, -0.36957471874357828,  0.99200950450550751, -0.02651441186298395,  0.12334556714490302,
             0.07071622802089356,  0.92650620200843881, -0.36957471874357839,  0.10448137664308749, -0.37534416687016492, -0.92097795757081347,  0.99200950450550751, -0.02651441186298373,  0.12334556714490325,
            -0.14283988479115467,  0.72510924905369167, -0.67365669613657164,  0.40856335403608091, -0.57674111982491227, -0.70742184475876524,  0.90148363992456226,  0.37627949404651062,  0.21387143172584827,
             0.40856335403608085, -0.57674111982491205, -0.70742184475876524, -0.14283988479115467,  0.72510924905369167, -0.67365669613657175,  0.90148363992456226,  0.37627949404651051,  0.21387143172584816,
             0.92097795757081369, -0.37534416687016492, -0.10448137664308743,  0.36957471874357817,  0.92650620200843903, -0.07071622802089378, -0.12334556714490325, -0.02651441186298384, -0.99200950450550729,
             0.36957471874357822,  0.92650620200843881, -0.07071622802089389,  0.92097795757081369, -0.37534416687016492, -0.10448137664308738, -0.12334556714490325, -0.02651441186298373, -0.99200950450550729,
             0.67365669613657175,  0.72510924905369145,  0.14283988479115434,  0.70742184475876546, -0.57674111982491216, -0.40856335403608079, -0.21387143172584849,  0.37627949404651073, -0.90148363992456204,
             0.70742184475876546, -0.57674111982491205, -0.40856335403608068,  0.67365669613657153,  0.72510924905369167,  0.14283988479115445, -0.21387143172584849,  0.37627949404651051, -0.90148363992456226,
            -0.12334556714490325, -0.02651441186298370, -0.99200950450550729,  0.92097795757081369, -0.37534416687016492, -0.10448137664308749,  0.36957471874357817,  0.92650620200843870, -0.07071622802089389,
             0.92097795757081369, -0.37534416687016514, -0.10448137664308754, -0.12334556714490325, -0.02651441186298364, -0.99200950450550729,  0.36957471874357822,  0.92650620200843892, -0.07071622802089383,
             0.70742184475876546, -0.57674111982491205, -0.40856335403608079, -0.21387143172584849,  0.37627949404651051, -0.90148363992456226,  0.67365669613657153,  0.72510924905369167,  0.14283988479115434,
            -0.21387143172584827,  0.37627949404651051, -0.90148363992456226,  0.70742184475876546, -0.57674111982491216, -0.40856335403608068,  0.67365669613657131,  0.72510924905369178,  0.14283988479115434,
             0.99200950450550729, -0.02651441186298378,  0.12334556714490302,  0.10448137664308749, -0.37534416687016503, -0.92097795757081347,  0.07071622802089378,  0.92650620200843892, -0.36957471874357806,
             0.10448137664308749, -0.37534416687016470, -0.92097795757081324,  0.99200950450550751, -0.02651441186298378,  0.12334556714490313,  0.07071622802089367,  0.92650620200843870, -0.36957471874357845,
             0.40856335403608091, -0.57674111982491216, -0.70742184475876524,  0.90148363992456226,  0.37627949404651045,  0.21387143172584827, -0.14283988479115456,  0.72510924905369178, -0.67365669613657164,
             0.90148363992456226,  0.37627949404651040,  0.21387143172584816,  0.40856335403608079, -0.57674111982491205, -0.70742184475876524, -0.14283988479115445,  0.72510924905369178, -0.67365669613657153,
             0.90148363992456204, -0.37627949404651051,  0.21387143172584844, -0.14283988479115445, -0.72510924905369190, -0.67365669613657142,  0.40856335403608085,  0.57674111982491216, -0.70742184475876524,
            -0.14283988479115434, -0.72510924905369167, -0.67365669613657142,  0.90148363992456204, -0.37627949404651051,  0.21387143172584844,  0.40856335403608052,  0.57674111982491194, -0.70742184475876524,
             0.07071622802089378, -0.92650620200843892, -0.36957471874357800,  0.99200950450550707,  0.02651441186298384,  0.12334556714490330,  0.10448137664308743,  0.37534416687016492, -0.92097795757081347,
             0.99200950450550707,  0.02651441186298376,  0.12334556714490341,  0.07071622802089378, -0.92650620200843892, -0.36957471874357806,  0.10448137664308749,  0.37534416687016492, -0.92097795757081335,
            -0.21387143172584805, -0.37627949404651062, -0.90148363992456215,  0.67365669613657153, -0.72510924905369167,  0.14283988479115484,  0.70742184475876480,  0.57674111982491216, -0.40856335403608074,
             0.67365669613657153, -0.72510924905369178,  0.14283988479115484, -0.21387143172584827, -0.37627949404651051, -0.90148363992456237,  0.70742184475876502,  0.57674111982491194, -0.40856335403608063,
             0.36957471874357833, -0.92650620200843881, -0.07071622802089350, -0.12334556714490325,  0.02651441186298362, -0.99200950450550729,  0.92097795757081324,  0.37534416687016503, -0.10448137664308738,
            -0.12334556714490325,  0.02651441186298359, -0.99200950450550729,  0.36957471874357828, -0.92650620200843881, -0.07071622802089339,  0.92097795757081324,  0.37534416687016503, -0.10448137664308743,
             0.92097795757081324,  0.37534416687016492, -0.10448137664308732,  0.36957471874357822, -0.92650620200843903, -0.07071622802089339, -0.12334556714490325,  0.02651441186298387, -0.99200950450550740,
             0.36957471874357828, -0.92650620200843881, -0.07071622802089339,  0.92097795757081347,  0.37534416687016503, -0.10448137664308726, -0.12334556714490336,  0.02651441186298370, -0.99200950450550740,
             0.67365669613657164, -0.72510924905369156,  0.14283988479115484,  0.70742184475876502,  0.57674111982491216, -0.40856335403608079, -0.21387143172584838, -0.37627949404651073, -0.90148363992456215,
             0.70742184475876502,  0.57674111982491205, -0.40856335403608068,  0.67365669613657153, -0.72510924905369156,  0.14283988479115484, -0.21387143172584816, -0.37627949404651051, -0.90148363992456215,
             0.10448137664308721,  0.37534416687016481, -0.92097795757081358,  0.07071622802089383, -0.92650620200843892, -0.36957471874357795,  0.99200950450550729,  0.02651441186298398,  0.12334556714490341,
             0.07071622802089389, -0.92650620200843881, -0.36957471874357806,  0.10448137664308726,  0.37534416687016481, -0.92097795757081358,  0.99200950450550707,  0.02651441186298376,  0.12334556714490341,
            -0.14283988479115445, -0.72510924905369167, -0.67365669613657142,  0.40856335403608063,  0.57674111982491227, -0.70742184475876535,  0.90148363992456204, -0.37627949404651062,  0.21387143172584866,
             0.40856335403608057,  0.57674111982491205, -0.70742184475876535, -0.14283988479115456, -0.72510924905369178, -0.67365669613657153,  0.90148363992456226, -0.37627949404651045,  0.21387143172584866,
             0.99200950450550685,  0.02651441186298389,  0.12334556714490341,  0.10448137664308738,  0.37534416687016492, -0.92097795757081358,  0.07071622802089400, -0.92650620200843892, -0.36957471874357800,
             0.10448137664308738,  0.37534416687016503, -0.92097795757081358,  0.99200950450550707,  0.02651441186298373,  0.12334556714490363,  0.07071622802089378, -0.92650620200843870, -0.36957471874357811,
             0.40856335403608052,  0.57674111982491205, -0.70742184475876524,  0.90148363992456204, -0.37627949404651040,  0.21387143172584855, -0.14283988479115434, -0.72510924905369178, -0.67365669613657131,
             0.90148363992456204, -0.37627949404651034,  0.21387143172584855,  0.40856335403608057,  0.57674111982491205, -0.70742184475876535, -0.14283988479115434, -0.72510924905369178, -0.67365669613657131,
            -0.12334556714490313,  0.02651441186298376, -0.99200950450550718,  0.92097795757081324,  0.37534416687016514, -0.10448137664308732,  0.36957471874357800, -0.92650620200843892, -0.07071622802089356,
             0.92097795757081324,  0.37534416687016492, -0.10448137664308732, -0.12334556714490313,  0.02651441186298373, -0.99200950450550718,  0.36957471874357833, -0.92650620200843870, -0.07071622802089356,
             0.70742184475876502,  0.57674111982491216, -0.40856335403608074, -0.21387143172584816, -0.37627949404651062, -0.90148363992456204,  0.67365669613657142, -0.72510924905369167,  0.14283988479115473,
            -0.21387143172584816, -0.37627949404651051, -0.90148363992456215,  0.70742184475876502,  0.57674111982491216, -0.40856335403608068,  0.67365669613657131, -0.72510924905369167,  0.14283988479115461,
             0.00743240513391588,  0.78867513459481287, -0.61476523277246165, -0.33474531727338763, -0.57735026918962584, -0.74472292782631966,  0.94227932796886937, -0.21132486540518702, -0.25971420705077936,
            -0.33474531727338752, -0.57735026918962562, -0.74472292782631944,  0.00743240513391596,  0.78867513459481287, -0.61476523277246176,  0.94227932796886926, -0.21132486540518713, -0.25971420705077952,
            -0.53230171741593746, -0.21132486540518713, -0.81975403804892766,  0.40254520541901551,  0.78867513459481287, -0.46470301232724548,  0.74472292782631955, -0.57735026918962573, -0.33474531727338752,
             0.40254520541901551,  0.78867513459481264, -0.46470301232724553, -0.53230171741593746, -0.21132486540518716, -0.81975403804892766,  0.74472292782631977, -0.57735026918962551, -0.33474531727338752,
             0.00743240513391630, -0.78867513459481287, -0.61476523277246165, -0.33474531727338785,  0.57735026918962584, -0.74472292782631955,  0.94227932796886926,  0.21132486540518697, -0.25971420705077941,
            -0.33474531727338797,  0.57735026918962584, -0.74472292782631966,  0.00743240513391638, -0.78867513459481287, -0.61476523277246153,  0.94227932796886926,  0.21132486540518702, -0.25971420705077941,
            -0.53230171741593746,  0.21132486540518716, -0.81975403804892766,  0.40254520541901573, -0.78867513459481287, -0.46470301232724537,  0.74472292782631944,  0.57735026918962562, -0.33474531727338763,
             0.40254520541901573, -0.78867513459481287, -0.46470301232724542, -0.53230171741593757,  0.21132486540518722, -0.81975403804892766,  0.74472292782631944,  0.57735026918962562, -0.33474531727338758,
             0.74472292782631921,  0.57735026918962573, -0.33474531727338758, -0.53230171741593757,  0.21132486540518727, -0.81975403804892766,  0.40254520541901612, -0.78867513459481287, -0.46470301232724531,
            -0.53230171741593746,  0.21132486540518725, -0.81975403804892766,  0.74472292782631921,  0.57735026918962584, -0.33474531727338752,  0.40254520541901589, -0.78867513459481287, -0.46470301232724531,
            -0.33474531727338785,  0.57735026918962595, -0.74472292782631966,  0.94227932796886937,  0.21132486540518700, -0.25971420705077941,  0.00743240513391619, -0.78867513459481287, -0.61476523277246165,
             0.94227932796886926,  0.21132486540518700, -0.25971420705077947, -0.33474531727338785,  0.57735026918962595, -0.74472292782631955,  0.00743240513391635, -0.78867513459481287, -0.61476523277246153,
             0.74472292782631966, -0.57735026918962595, -0.33474531727338747, -0.53230171741593757, -0.21132486540518700, -0.81975403804892766,  0.40254520541901562,  0.78867513459481287, -0.46470301232724548,
            -0.53230171741593746, -0.21132486540518700, -0.81975403804892766,  0.74472292782631966, -0.57735026918962595, -0.33474531727338752,  0.40254520541901545,  0.78867513459481287, -0.46470301232724553,
            -0.33474531727338741, -0.57735026918962573, -0.74472292782631944,  0.94227932796886937, -0.21132486540518727, -0.25971420705077936,  0.00743240513391569,  0.78867513459481287, -0.61476523277246176,
             0.94227932796886926, -0.21132486540518725, -0.25971420705077936, -0.33474531727338741, -0.57735026918962584, -0.74472292782631955,  0.00743240513391591,  0.78867513459481287, -0.61476523277246176,
             0.40254520541901551,  0.78867513459481287, -0.46470301232724542,  0.74472292782631966, -0.57735026918962584, -0.33474531727338752, -0.53230171741593746, -0.21132486540518697, -0.81975403804892766,
             0.74472292782631977, -0.57735026918962584, -0.33474531727338741,  0.40254520541901539,  0.78867513459481287, -0.46470301232724553, -0.53230171741593746, -0.21132486540518702, -0.81975403804892766,
             0.94227932796886926, -0.21132486540518716, -0.25971420705077941,  0.00743240513391605,  0.78867513459481287, -0.61476523277246165, -0.33474531727338763, -0.57735026918962562, -0.74472292782631944,
             0.00743240513391608,  0.78867513459481287, -0.61476523277246165,  0.94227932796886937, -0.21132486540518722, -0.25971420705077941, -0.33474531727338763, -0.57735026918962562, -0.74472292782631944,
             0.40254520541901595, -0.78867513459481287, -0.46470301232724542,  0.74472292782631944,  0.57735026918962584, -0.33474531727338747, -0.53230171741593757,  0.21132486540518702, -0.81975403804892766,
             0.74472292782631933,  0.57735026918962562, -0.33474531727338763,  0.40254520541901584, -0.78867513459481287, -0.46470301232724537, -0.53230171741593746,  0.21132486540518713, -0.81975403804892755,
             0.94227932796886926,  0.21132486540518713, -0.25971420705077941,  0.00743240513391630, -0.78867513459481287, -0.61476523277246165, -0.33474531727338774,  0.57735026918962573, -0.74472292782631955,
             0.00743240513391630, -0.78867513459481264, -0.61476523277246153,  0.94227932796886926,  0.21132486540518716, -0.25971420705077947, -0.33474531727338797,  0.57735026918962551, -0.74472292782631955,
            -0.78867513459481287,  0.61476523277246187, -0.00743240513391580,  0.57735026918962562,  0.74472292782631966,  0.33474531727338763,  0.21132486540518713,  0.25971420705077947, -0.94227932796886948,
             0.57735026918962551,  0.74472292782631966,  0.33474531727338774, -0.78867513459481287,  0.61476523277246187, -0.00743240513391585,  0.21132486540518741,  0.25971420705077947, -0.94227932796886948,
             0.21132486540518691,  0.81975403804892766,  0.53230171741593768, -0.78867513459481287,  0.46470301232724565, -0.40254520541901551,  0.57735026918962584,  0.33474531727338752, -0.74472292782631988,
            -0.78867513459481275,  0.46470301232724570, -0.40254520541901539,  0.21132486540518697,  0.81975403804892766,  0.53230171741593768,  0.57735026918962573,  0.33474531727338752, -0.74472292782631988,
             0.78867513459481275,  0.61476523277246165, -0.00743240513391630, -0.57735026918962595,  0.74472292782631988,  0.33474531727338808, -0.21132486540518669,  0.25971420705077952, -0.94227932796886926,
            -0.57735026918962595,  0.74472292782631988,  0.33474531727338808,  0.78867513459481275,  0.61476523277246165, -0.00743240513391641, -0.21132486540518686,  0.25971420705077952, -0.94227932796886926,
            -0.21132486540518741,  0.81975403804892777,  0.53230171741593768,  0.78867513459481320,  0.46470301232724542, -0.40254520541901584, -0.57735026918962562,  0.33474531727338769, -0.74472292782631944,
             0.78867513459481298,  0.46470301232724542, -0.40254520541901578, -0.21132486540518744,  0.81975403804892788,  0.53230171741593779, -0.57735026918962551,  0.33474531727338763, -0.74472292782631966,
            -0.57735026918962562,  0.33474531727338774, -0.74472292782631933, -0.21132486540518747,  0.81975403804892777,  0.53230171741593768,  0.78867513459481287,  0.46470301232724537, -0.40254520541901612,
            -0.21132486540518736,  0.81975403804892788,  0.53230171741593768, -0.57735026918962562,  0.33474531727338774, -0.74472292782631933,  0.78867513459481287,  0.46470301232724537, -0.40254520541901595,
            -0.57735026918962606,  0.74472292782631966,  0.33474531727338797, -0.21132486540518680,  0.25971420705077952, -0.94227932796886948,  0.78867513459481287,  0.61476523277246165, -0.00743240513391608,
            -0.21132486540518680,  0.25971420705077952, -0.94227932796886948, -0.57735026918962606,  0.74472292782631966,  0.33474531727338797,  0.78867513459481287,  0.61476523277246165, -0.00743240513391624,
             0.57735026918962606,  0.33474531727338758, -0.74472292782631977,  0.21132486540518680,  0.81975403804892777,  0.53230171741593768, -0.78867513459481287,  0.46470301232724559, -0.40254520541901573,
             0.21132486540518680,  0.81975403804892777,  0.53230171741593768,  0.57735026918962606,  0.33474531727338758, -0.74472292782631977, -0.78867513459481287,  0.46470301232724559, -0.40254520541901556,
             0.57735026918962562,  0.74472292782631955,  0.33474531727338752,  0.21132486540518747,  0.25971420705077952, -0.94227932796886948, -0.78867513459481287,  0.61476523277246187, -0.00743240513391569,
             0.21132486540518736,  0.25971420705077947, -0.94227932796886948,  0.57735026918962562,  0.74472292782631955,  0.33474531727338752, -0.78867513459481287,  0.61476523277246198, -0.00743240513391585,
            -0.78867513459481275,  0.46470301232724565, -0.40254520541901551,  0.57735026918962595,  0.33474531727338747, -0.74472292782631988,  0.21132486540518669,  0.81975403804892777,  0.53230171741593746,
             0.57735026918962595,  0.33474531727338747, -0.74472292782631988, -0.78867513459481275,  0.46470301232724570, -0.40254520541901539,  0.21132486540518686,  0.81975403804892777,  0.53230171741593746,
             0.21132486540518741,  0.25971420705077952, -0.94227932796886948, -0.78867513459481320,  0.61476523277246187, -0.00743240513391596,  0.57735026918962562,  0.74472292782631966,  0.33474531727338763,
            -0.78867513459481298,  0.61476523277246187, -0.00743240513391602,  0.21132486540518744,  0.25971420705077947, -0.94227932796886960,  0.57735026918962551,  0.74472292782631966,  0.33474531727338785,
             0.78867513459481287,  0.46470301232724537, -0.40254520541901601, -0.57735026918962562,  0.33474531727338769, -0.74472292782631944, -0.21132486540518713,  0.81975403804892788,  0.53230171741593768,
            -0.57735026918962551,  0.33474531727338769, -0.74472292782631955,  0.78867513459481287,  0.46470301232724542, -0.40254520541901595, -0.21132486540518741,  0.81975403804892788,  0.53230171741593768,
            -0.21132486540518691,  0.25971420705077958, -0.94227932796886948,  0.78867513459481287,  0.61476523277246165, -0.00743240513391633, -0.57735026918962584,  0.74472292782631977,  0.33474531727338808,
             0.78867513459481275,  0.61476523277246153, -0.00743240513391638, -0.21132486540518697,  0.25971420705077958, -0.94227932796886948, -0.57735026918962573,  0.74472292782631977,  0.33474531727338808,
             0.37627949404651062, -0.21387143172584827, -0.90148363992456226,  0.72510924905369156,  0.67365669613657175,  0.14283988479115445, -0.57674111982491205,  0.70742184475876524, -0.40856335403608079,
             0.72510924905369156,  0.67365669613657164,  0.14283988479115423,  0.37627949404651062, -0.21387143172584827, -0.90148363992456226, -0.57674111982491205,  0.70742184475876524, -0.40856335403608063,
             0.92650620200843881,  0.36957471874357828, -0.07071622802089389, -0.02651441186298359, -0.12334556714490325, -0.99200950450550729, -0.37534416687016503,  0.92097795757081347, -0.10448137664308749,
            -0.02651441186298353, -0.12334556714490336, -0.99200950450550707,  0.92650620200843881,  0.36957471874357828, -0.07071622802089378, -0.37534416687016514,  0.92097795757081347, -0.10448137664308760,
             0.37627949404651051,  0.90148363992456226,  0.21387143172584816,  0.72510924905369178, -0.14283988479115456, -0.67365669613657175, -0.57674111982491194,  0.40856335403608079, -0.70742184475876513,
             0.72510924905369178, -0.14283988479115456, -0.67365669613657164,  0.37627949404651040,  0.90148363992456249,  0.21387143172584827, -0.57674111982491194,  0.40856335403608074, -0.70742184475876524,
             0.92650620200843870,  0.07071622802089367, -0.36957471874357833, -0.02651441186298367,  0.99200950450550729,  0.12334556714490325, -0.37534416687016481,  0.10448137664308754, -0.92097795757081347,
            -0.02651441186298364,  0.99200950450550751,  0.12334556714490325,  0.92650620200843870,  0.07071622802089367, -0.36957471874357839, -0.37534416687016492,  0.10448137664308754, -0.92097795757081347,
            -0.37534416687016470,  0.10448137664308738, -0.92097795757081347,  0.92650620200843881,  0.07071622802089367, -0.36957471874357839, -0.02651441186298395,  0.99200950450550751,  0.12334556714490313,
             0.92650620200843881,  0.07071622802089367, -0.36957471874357839, -0.37534416687016481,  0.10448137664308726, -0.92097795757081347, -0.02651441186298389,  0.99200950450550751,  0.12334556714490325,
             0.72510924905369167, -0.14283988479115467, -0.67365669613657175, -0.57674111982491216,  0.40856335403608079, -0.70742184475876524,  0.37627949404651062,  0.90148363992456249,  0.21387143172584827,
            -0.57674111982491205,  0.40856335403608079, -0.70742184475876513,  0.72510924905369167, -0.14283988479115467, -0.67365669613657175,  0.37627949404651051,  0.90148363992456249,  0.21387143172584827,
            -0.37534416687016492,  0.92097795757081369, -0.10448137664308738,  0.92650620200843892,  0.36957471874357828, -0.07071622802089383, -0.02651441186298378, -0.12334556714490336, -0.99200950450550729,
             0.92650620200843881,  0.36957471874357828, -0.07071622802089395, -0.37534416687016492,  0.92097795757081369, -0.10448137664308732, -0.02651441186298367, -0.12334556714490336, -0.99200950450550729,
             0.72510924905369145,  0.67365669613657175,  0.14283988479115434, -0.57674111982491205,  0.70742184475876546, -0.40856335403608068,  0.37627949404651073, -0.21387143172584849, -0.90148363992456226,
            -0.57674111982491205,  0.70742184475876546, -0.40856335403608068,  0.72510924905369156,  0.67365669613657175,  0.14283988479115445,  0.37627949404651073, -0.21387143172584849, -0.90148363992456226,
            -0.02651441186298370, -0.12334556714490325, -0.99200950450550729, -0.37534416687016503,  0.92097795757081369, -0.10448137664308738,  0.92650620200843870,  0.36957471874357828, -0.07071622802089395,
            -0.37534416687016503,  0.92097795757081369, -0.10448137664308743, -0.02651441186298359, -0.12334556714490336, -0.99200950450550729,  0.92650620200843870,  0.36957471874357833, -0.07071622802089389,
            -0.57674111982491194,  0.70742184475876546, -0.40856335403608074,  0.37627949404651051, -0.21387143172584827, -0.90148363992456226,  0.72510924905369156,  0.67365669613657153,  0.14283988479115445,
             0.37627949404651045, -0.21387143172584838, -0.90148363992456226, -0.57674111982491194,  0.70742184475876546, -0.40856335403608068,  0.72510924905369167,  0.67365669613657153,  0.14283988479115434,
            -0.02651441186298381,  0.99200950450550729,  0.12334556714490325, -0.37534416687016481,  0.10448137664308738, -0.92097795757081347,  0.92650620200843881,  0.07071622802089383, -0.36957471874357828,
            -0.37534416687016481,  0.10448137664308743, -0.92097795757081324, -0.02651441186298381,  0.99200950450550729,  0.12334556714490313,  0.92650620200843881,  0.07071622802089378, -0.36957471874357845,
            -0.57674111982491205,  0.40856335403608079, -0.70742184475876524,  0.37627949404651040,  0.90148363992456226,  0.21387143172584816,  0.72510924905369178, -0.14283988479115445, -0.67365669613657153,
             0.37627949404651034,  0.90148363992456249,  0.21387143172584805, -0.57674111982491205,  0.40856335403608079, -0.70742184475876524,  0.72510924905369190, -0.14283988479115445, -0.67365669613657153,
            -0.37627949404651062,  0.90148363992456226,  0.21387143172584844, -0.72510924905369156, -0.14283988479115456, -0.67365669613657153,  0.57674111982491205,  0.40856335403608068, -0.70742184475876502,
            -0.72510924905369145, -0.14283988479115445, -0.67365669613657142, -0.37627949404651062,  0.90148363992456204,  0.21387143172584844,  0.57674111982491183,  0.40856335403608063, -0.70742184475876524,
            -0.92650620200843870,  0.07071622802089378, -0.36957471874357817,  0.02651441186298359,  0.99200950450550707,  0.12334556714490341,  0.37534416687016492,  0.10448137664308738, -0.92097795757081335,
             0.02651441186298356,  0.99200950450550707,  0.12334556714490352, -0.92650620200843881,  0.07071622802089367, -0.36957471874357817,  0.37534416687016514,  0.10448137664308743, -0.92097795757081335,
            -0.37627949404651040, -0.21387143172584827, -0.90148363992456215, -0.72510924905369167,  0.67365669613657153,  0.14283988479115484,  0.57674111982491194,  0.70742184475876502, -0.40856335403608074,
            -0.72510924905369167,  0.67365669613657153,  0.14283988479115461, -0.37627949404651040, -0.21387143172584838, -0.90148363992456226,  0.57674111982491183,  0.70742184475876513, -0.40856335403608052,
            -0.92650620200843870,  0.36957471874357822, -0.07071622802089356,  0.02651441186298370, -0.12334556714490336, -0.99200950450550718,  0.37534416687016481,  0.92097795757081347, -0.10448137664308732,
             0.02651441186298367, -0.12334556714490336, -0.99200950450550729, -0.92650620200843870,  0.36957471874357828, -0.07071622802089350,  0.37534416687016492,  0.92097795757081324, -0.10448137664308738,
             0.37534416687016459,  0.92097795757081347, -0.10448137664308721, -0.92650620200843881,  0.36957471874357828, -0.07071622802089350,  0.02651441186298398, -0.12334556714490336, -0.99200950450550740,
            -0.92650620200843870,  0.36957471874357828, -0.07071622802089339,  0.37534416687016470,  0.92097795757081347, -0.10448137664308726,  0.02651441186298387, -0.12334556714490336, -0.99200950450550740,
            -0.72510924905369167,  0.67365669613657153,  0.14283988479115484,  0.57674111982491205,  0.70742184475876513, -0.40856335403608068, -0.37627949404651051, -0.21387143172584838, -0.90148363992456215,
             0.57674111982491194,  0.70742184475876502, -0.40856335403608068, -0.72510924905369167,  0.67365669613657164,  0.14283988479115484, -0.37627949404651045, -0.21387143172584838, -0.90148363992456215,
             0.37534416687016492,  0.10448137664308726, -0.92097795757081358, -0.92650620200843881,  0.07071622802089372, -0.36957471874357806,  0.02651441186298376,  0.99200950450550729,  0.12334556714490341,
            -0.92650620200843870,  0.07071622802089383, -0.36957471874357806,  0.37534416687016492,  0.10448137664308721, -0.92097795757081358,  0.02651441186298370,  0.99200950450550729,  0.12334556714490341,
            -0.72510924905369134, -0.14283988479115456, -0.67365669613657153,  0.57674111982491205,  0.40856335403608057, -0.70742184475876524, -0.37627949404651073,  0.90148363992456226,  0.21387143172584866,
             0.57674111982491194,  0.40856335403608057, -0.70742184475876535, -0.72510924905369145, -0.14283988479115456, -0.67365669613657153, -0.37627949404651062,  0.90148363992456226,  0.21387143172584866,
             0.02651441186298362,  0.99200950450550707,  0.12334556714490341,  0.37534416687016492,  0.10448137664308726, -0.92097795757081358, -0.92650620200843870,  0.07071622802089383, -0.36957471874357800,
             0.37534416687016492,  0.10448137664308726, -0.92097795757081335,  0.02651441186298362,  0.99200950450550729,  0.12334556714490352, -0.92650620200843858,  0.07071622802089372, -0.36957471874357822,
             0.57674111982491194,  0.40856335403608063, -0.70742184475876524, -0.37627949404651045,  0.90148363992456226,  0.21387143172584844, -0.72510924905369156, -0.14283988479115456, -0.67365669613657142,
            -0.37627949404651040,  0.90148363992456226,  0.21387143172584855,  0.57674111982491194,  0.40856335403608057, -0.70742184475876524, -0.72510924905369167, -0.14283988479115445, -0.67365669613657131,
             0.02651441186298389, -0.12334556714490336, -0.99200950450550718,  0.37534416687016481,  0.92097795757081347, -0.10448137664308721, -0.92650620200843881,  0.36957471874357817, -0.07071622802089367,
             0.37534416687016470,  0.92097795757081324, -0.10448137664308732,  0.02651441186298384, -0.12334556714490325, -0.99200950450550718, -0.92650620200843858,  0.36957471874357822, -0.07071622802089356,
             0.57674111982491194,  0.70742184475876502, -0.40856335403608057, -0.37627949404651034, -0.21387143172584827, -0.90148363992456215, -0.72510924905369167,  0.67365669613657153,  0.14283988479115461,
            -0.37627949404651029, -0.21387143172584827, -0.90148363992456226,  0.57674111982491205,  0.70742184475876524, -0.40856335403608057, -0.72510924905369190,  0.67365669613657142,  0.14283988479115461,
             0.78867513459481287,  0.00743240513391580, -0.61476523277246176, -0.57735026918962562, -0.33474531727338763, -0.74472292782631944, -0.21132486540518713,  0.94227932796886948, -0.25971420705077936,
            -0.57735026918962551, -0.33474531727338763, -0.74472292782631944,  0.78867513459481287,  0.00743240513391588, -0.61476523277246176, -0.21132486540518741,  0.94227932796886937, -0.25971420705077941,
            -0.21132486540518691, -0.53230171741593768, -0.81975403804892766,  0.78867513459481287,  0.40254520541901551, -0.46470301232724553, -0.57735026918962584,  0.74472292782631988, -0.33474531727338741,
             0.78867513459481287,  0.40254520541901539, -0.46470301232724553, -0.21132486540518700, -0.53230171741593768, -0.81975403804892766, -0.57735026918962584,  0.74472292782631988, -0.33474531727338741,
            -0.78867513459481275,  0.00743240513391630, -0.61476523277246153,  0.57735026918962595, -0.33474531727338808, -0.74472292782631966,  0.21132486540518669,  0.94227932796886937, -0.25971420705077941,
             0.57735026918962595, -0.33474531727338808, -0.74472292782631966, -0.78867513459481264,  0.00743240513391630, -0.61476523277246153,  0.21132486540518680,  0.94227932796886926, -0.25971420705077947,
             0.21132486540518738, -0.53230171741593768, -0.81975403804892766, -0.78867513459481320,  0.40254520541901584, -0.46470301232724531,  0.57735026918962562,  0.74472292782631944, -0.33474531727338763,
            -0.78867513459481298,  0.40254520541901578, -0.46470301232724531,  0.21132486540518744, -0.53230171741593779, -0.81975403804892766,  0.57735026918962551,  0.74472292782631966, -0.33474531727338752,
             0.57735026918962551,  0.74472292782631933, -0.33474531727338763,  0.21132486540518747, -0.53230171741593768, -0.81975403804892766, -0.78867513459481287,  0.40254520541901612, -0.46470301232724526,
             0.21132486540518736, -0.53230171741593768, -0.81975403804892766,  0.57735026918962562,  0.74472292782631944, -0.33474531727338763, -0.78867513459481287,  0.40254520541901589, -0.46470301232724531,
             0.57735026918962595, -0.33474531727338785, -0.74472292782631966,  0.21132486540518688,  0.94227932796886948, -0.25971420705077941, -0.78867513459481287,  0.00743240513391608, -0.61476523277246153,
             0.21132486540518680,  0.94227932796886948, -0.25971420705077947,  0.57735026918962606, -0.33474531727338797, -0.74472292782631966, -0.78867513459481287,  0.00743240513391624, -0.61476523277246153,
            -0.57735026918962595,  0.74472292782631966, -0.33474531727338747, -0.21132486540518688, -0.53230171741593768, -0.81975403804892766,  0.78867513459481287,  0.40254520541901573, -0.46470301232724548,
            -0.21132486540518680, -0.53230171741593768, -0.81975403804892766, -0.57735026918962606,  0.74472292782631977, -0.33474531727338741,  0.78867513459481287,  0.40254520541901556, -0.46470301232724553,
            -0.57735026918962551, -0.33474531727338752, -0.74472292782631944, -0.21132486540518747,  0.94227932796886948, -0.25971420705077936,  0.78867513459481287,  0.00743240513391569, -0.61476523277246176,
            -0.21132486540518736,  0.94227932796886948, -0.25971420705077936, -0.57735026918962562, -0.33474531727338763, -0.74472292782631944,  0.78867513459481287,  0.00743240513391591, -0.61476523277246176,
             0.78867513459481275,  0.40254520541901551, -0.46470301232724553, -0.57735026918962595,  0.74472292782631988, -0.33474531727338741, -0.21132486540518669, -0.53230171741593757, -0.81975403804892766,
            -0.57735026918962595,  0.74472292782631988, -0.33474531727338741,  0.78867513459481264,  0.40254520541901551, -0.46470301232724553, -0.21132486540518680, -0.53230171741593746, -0.81975403804892766,
            -0.21132486540518738,  0.94227932796886948, -0.25971420705077941,  0.78867513459481320,  0.00743240513391596, -0.61476523277246176, -0.57735026918962562, -0.33474531727338763, -0.74472292782631944,
             0.78867513459481298,  0.00743240513391602, -0.61476523277246176, -0.21132486540518744,  0.94227932796886960, -0.25971420705077936, -0.57735026918962551, -0.33474531727338785, -0.74472292782631955,
            -0.78867513459481287,  0.40254520541901601, -0.46470301232724526,  0.57735026918962562,  0.74472292782631944, -0.33474531727338758,  0.21132486540518713, -0.53230171741593768, -0.81975403804892766,
             0.57735026918962551,  0.74472292782631944, -0.33474531727338763, -0.78867513459481287,  0.40254520541901595, -0.46470301232724531,  0.21132486540518741, -0.53230171741593757, -0.81975403804892766,
             0.21132486540518691,  0.94227932796886948, -0.25971420705077941, -0.78867513459481287,  0.00743240513391633, -0.61476523277246153,  0.57735026918962584, -0.33474531727338808, -0.74472292782631966,
            -0.78867513459481287,  0.00743240513391638, -0.61476523277246153,  0.21132486540518700,  0.94227932796886948, -0.25971420705077947,  0.57735026918962584, -0.33474531727338808, -0.74472292782631966,
        });
    }

    private static void populateFrames3D_Part4() {
        populateFrames3D(7776, new double[] {
             0.21387143172584838, -0.90148363992456226, -0.37627949404651062, -0.67365669613657153,  0.14283988479115456, -0.72510924905369167, -0.70742184475876524, -0.40856335403608068,  0.57674111982491216,
            -0.67365669613657142,  0.14283988479115434, -0.72510924905369145,  0.21387143172584849, -0.90148363992456204, -0.37627949404651051, -0.70742184475876535, -0.40856335403608041,  0.57674111982491194,
            -0.36957471874357811, -0.07071622802089383, -0.92650620200843892,  0.12334556714490347, -0.99200950450550729,  0.02651441186298367, -0.92097795757081369, -0.10448137664308726,  0.37534416687016503,
             0.12334556714490336, -0.99200950450550729,  0.02651441186298362, -0.36957471874357806, -0.07071622802089383, -0.92650620200843881, -0.92097795757081347, -0.10448137664308726,  0.37534416687016503,
            -0.90148363992456226,  0.21387143172584849, -0.37627949404651040,  0.14283988479115478, -0.67365669613657175, -0.72510924905369190, -0.40856335403608079, -0.70742184475876502,  0.57674111982491205,
             0.14283988479115478, -0.67365669613657164, -0.72510924905369190, -0.90148363992456249,  0.21387143172584838, -0.37627949404651040, -0.40856335403608068, -0.70742184475876502,  0.57674111982491194,
            -0.07071622802089356, -0.36957471874357839, -0.92650620200843892, -0.99200950450550751,  0.12334556714490347,  0.02651441186298384, -0.10448137664308738, -0.92097795757081324,  0.37534416687016492,
            -0.99200950450550751,  0.12334556714490336,  0.02651441186298373, -0.07071622802089345, -0.36957471874357845, -0.92650620200843892, -0.10448137664308743, -0.92097795757081324,  0.37534416687016503,
            -0.10448137664308732, -0.92097795757081324,  0.37534416687016470, -0.07071622802089345, -0.36957471874357839, -0.92650620200843892, -0.99200950450550751,  0.12334556714490347,  0.02651441186298406,
            -0.07071622802089345, -0.36957471874357845, -0.92650620200843892, -0.10448137664308721, -0.92097795757081324,  0.37534416687016481, -0.99200950450550751,  0.12334556714490347,  0.02651441186298389,
             0.14283988479115478, -0.67365669613657186, -0.72510924905369167, -0.40856335403608079, -0.70742184475876502,  0.57674111982491216, -0.90148363992456226,  0.21387143172584838, -0.37627949404651051,
            -0.40856335403608079, -0.70742184475876502,  0.57674111982491205,  0.14283988479115478, -0.67365669613657175, -0.72510924905369178, -0.90148363992456226,  0.21387143172584838, -0.37627949404651040,
            -0.92097795757081369, -0.10448137664308699,  0.37534416687016492, -0.36957471874357806, -0.07071622802089389, -0.92650620200843892,  0.12334556714490336, -0.99200950450550729,  0.02651441186298384,
            -0.36957471874357806, -0.07071622802089389, -0.92650620200843881, -0.92097795757081369, -0.10448137664308710,  0.37534416687016503,  0.12334556714490336, -0.99200950450550729,  0.02651441186298373,
            -0.67365669613657153,  0.14283988479115445, -0.72510924905369145, -0.70742184475876546, -0.40856335403608046,  0.57674111982491216,  0.21387143172584860, -0.90148363992456226, -0.37627949404651073,
            -0.70742184475876546, -0.40856335403608041,  0.57674111982491216, -0.67365669613657164,  0.14283988479115445, -0.72510924905369156,  0.21387143172584860, -0.90148363992456226, -0.37627949404651062,
             0.12334556714490347, -0.99200950450550729,  0.02651441186298373, -0.92097795757081369, -0.10448137664308715,  0.37534416687016514, -0.36957471874357806, -0.07071622802089383, -0.92650620200843881,
            -0.92097795757081369, -0.10448137664308721,  0.37534416687016514,  0.12334556714490358, -0.99200950450550729,  0.02651441186298367, -0.36957471874357817, -0.07071622802089389, -0.92650620200843870,
            -0.70742184475876524, -0.40856335403608046,  0.57674111982491216,  0.21387143172584860, -0.90148363992456226, -0.37627949404651051, -0.67365669613657153,  0.14283988479115445, -0.72510924905369167,
             0.21387143172584860, -0.90148363992456226, -0.37627949404651040, -0.70742184475876546, -0.40856335403608041,  0.57674111982491216, -0.67365669613657142,  0.14283988479115434, -0.72510924905369178,
            -0.99200950450550729,  0.12334556714490336,  0.02651441186298392, -0.10448137664308726, -0.92097795757081347,  0.37534416687016492, -0.07071622802089367, -0.36957471874357817, -0.92650620200843892,
            -0.10448137664308743, -0.92097795757081324,  0.37534416687016470, -0.99200950450550729,  0.12334556714490325,  0.02651441186298384, -0.07071622802089350, -0.36957471874357845, -0.92650620200843870,
            -0.40856335403608074, -0.70742184475876502,  0.57674111982491216, -0.90148363992456226,  0.21387143172584838, -0.37627949404651040,  0.14283988479115478, -0.67365669613657153, -0.72510924905369178,
            -0.90148363992456226,  0.21387143172584838, -0.37627949404651029, -0.40856335403608079, -0.70742184475876502,  0.57674111982491205,  0.14283988479115467, -0.67365669613657153, -0.72510924905369178,
            -0.61476523277246187, -0.00743240513391605,  0.78867513459481253, -0.74472292782631944,  0.33474531727338758, -0.57735026918962573, -0.25971420705077974, -0.94227932796886926, -0.21132486540518666,
            -0.74472292782631933,  0.33474531727338758, -0.57735026918962562, -0.61476523277246176, -0.00743240513391610,  0.78867513459481253, -0.25971420705077974, -0.94227932796886926, -0.21132486540518677,
            -0.81975403804892744,  0.53230171741593746, -0.21132486540518716, -0.46470301232724570, -0.40254520541901562,  0.78867513459481264, -0.33474531727338774, -0.74472292782631966, -0.57735026918962540,
            -0.46470301232724576, -0.40254520541901562,  0.78867513459481264, -0.81975403804892755,  0.53230171741593746, -0.21132486540518722, -0.33474531727338763, -0.74472292782631966, -0.57735026918962529,
            -0.61476523277246153, -0.00743240513391624, -0.78867513459481242, -0.74472292782631955,  0.33474531727338769,  0.57735026918962551, -0.25971420705077986, -0.94227932796886926,  0.21132486540518694,
            -0.74472292782631966,  0.33474531727338780,  0.57735026918962551, -0.61476523277246153, -0.00743240513391630, -0.78867513459481253, -0.25971420705077986, -0.94227932796886926,  0.21132486540518711,
            -0.81975403804892755,  0.53230171741593746,  0.21132486540518688, -0.46470301232724548, -0.40254520541901567, -0.78867513459481264, -0.33474531727338785, -0.74472292782631966,  0.57735026918962573,
            -0.46470301232724553, -0.40254520541901562, -0.78867513459481264, -0.81975403804892755,  0.53230171741593746,  0.21132486540518700, -0.33474531727338785, -0.74472292782631966,  0.57735026918962562,
            -0.33474531727338797, -0.74472292782631921,  0.57735026918962562, -0.81975403804892766,  0.53230171741593746,  0.21132486540518705, -0.46470301232724542, -0.40254520541901601, -0.78867513459481264,
            -0.81975403804892755,  0.53230171741593746,  0.21132486540518694, -0.33474531727338791, -0.74472292782631944,  0.57735026918962573, -0.46470301232724548, -0.40254520541901573, -0.78867513459481264,
            -0.74472292782631944,  0.33474531727338758,  0.57735026918962551, -0.25971420705077980, -0.94227932796886948,  0.21132486540518716, -0.61476523277246165, -0.00743240513391602, -0.78867513459481264,
            -0.25971420705077986, -0.94227932796886926,  0.21132486540518722, -0.74472292782631944,  0.33474531727338758,  0.57735026918962551, -0.61476523277246165, -0.00743240513391613, -0.78867513459481264,
            -0.33474531727338780, -0.74472292782631944, -0.57735026918962551, -0.81975403804892744,  0.53230171741593746, -0.21132486540518711, -0.46470301232724565, -0.40254520541901589,  0.78867513459481264,
            -0.81975403804892744,  0.53230171741593746, -0.21132486540518716, -0.33474531727338780, -0.74472292782631944, -0.57735026918962551, -0.46470301232724565, -0.40254520541901578,  0.78867513459481264,
            -0.74472292782631933,  0.33474531727338735, -0.57735026918962562, -0.25971420705077969, -0.94227932796886926, -0.21132486540518700, -0.61476523277246187, -0.00743240513391591,  0.78867513459481264,
            -0.25971420705077974, -0.94227932796886926, -0.21132486540518688, -0.74472292782631944,  0.33474531727338758, -0.57735026918962573, -0.61476523277246176, -0.00743240513391619,  0.78867513459481264,
            -0.46470301232724576, -0.40254520541901567,  0.78867513459481242, -0.33474531727338774, -0.74472292782631966, -0.57735026918962551, -0.81975403804892744,  0.53230171741593746, -0.21132486540518688,
            -0.33474531727338763, -0.74472292782631966, -0.57735026918962551, -0.46470301232724576, -0.40254520541901562,  0.78867513459481253, -0.81975403804892744,  0.53230171741593746, -0.21132486540518705,
            -0.25971420705077974, -0.94227932796886926, -0.21132486540518683, -0.61476523277246176, -0.00743240513391624,  0.78867513459481264, -0.74472292782631944,  0.33474531727338769, -0.57735026918962573,
            -0.61476523277246176, -0.00743240513391627,  0.78867513459481264, -0.25971420705077974, -0.94227932796886948, -0.21132486540518694, -0.74472292782631944,  0.33474531727338769, -0.57735026918962562,
            -0.46470301232724548, -0.40254520541901584, -0.78867513459481253, -0.33474531727338785, -0.74472292782631944,  0.57735026918962573, -0.81975403804892755,  0.53230171741593746,  0.21132486540518672,
            -0.33474531727338797, -0.74472292782631944,  0.57735026918962562, -0.46470301232724548, -0.40254520541901584, -0.78867513459481253, -0.81975403804892755,  0.53230171741593746,  0.21132486540518683,
            -0.25971420705077980, -0.94227932796886948,  0.21132486540518722, -0.61476523277246153, -0.00743240513391627, -0.78867513459481264, -0.74472292782631955,  0.33474531727338769,  0.57735026918962540,
            -0.61476523277246153, -0.00743240513391627, -0.78867513459481264, -0.25971420705077974, -0.94227932796886948,  0.21132486540518727, -0.74472292782631966,  0.33474531727338780,  0.57735026918962529,
            -0.00743240513391610, -0.61476523277246142, -0.78867513459481264,  0.33474531727338780, -0.74472292782631966,  0.57735026918962573, -0.94227932796886948, -0.25971420705077991,  0.21132486540518700,
             0.33474531727338769, -0.74472292782631955,  0.57735026918962540, -0.00743240513391613, -0.61476523277246142, -0.78867513459481253, -0.94227932796886926, -0.25971420705078002,  0.21132486540518716,
             0.53230171741593746, -0.81975403804892766,  0.21132486540518688, -0.40254520541901584, -0.46470301232724537, -0.78867513459481264, -0.74472292782631966, -0.33474531727338797,  0.57735026918962573,
            -0.40254520541901573, -0.46470301232724542, -0.78867513459481264,  0.53230171741593768, -0.81975403804892766,  0.21132486540518691, -0.74472292782631988, -0.33474531727338797,  0.57735026918962573,
            -0.00743240513391608, -0.61476523277246198,  0.78867513459481253,  0.33474531727338769, -0.74472292782631921, -0.57735026918962584, -0.94227932796886948, -0.25971420705077980, -0.21132486540518672,
             0.33474531727338769, -0.74472292782631921, -0.57735026918962584, -0.00743240513391613, -0.61476523277246187,  0.78867513459481275, -0.94227932796886948, -0.25971420705077986, -0.21132486540518688,
             0.53230171741593746, -0.81975403804892744, -0.21132486540518727, -0.40254520541901567, -0.46470301232724587,  0.78867513459481287, -0.74472292782631966, -0.33474531727338763, -0.57735026918962540,
            -0.40254520541901562, -0.46470301232724587,  0.78867513459481275,  0.53230171741593746, -0.81975403804892744, -0.21132486540518730, -0.74472292782631966, -0.33474531727338763, -0.57735026918962540,
            -0.74472292782631944, -0.33474531727338780, -0.57735026918962551,  0.53230171741593746, -0.81975403804892744, -0.21132486540518733, -0.40254520541901589, -0.46470301232724576,  0.78867513459481287,
             0.53230171741593746, -0.81975403804892732, -0.21132486540518727, -0.74472292782631966, -0.33474531727338774, -0.57735026918962551, -0.40254520541901584, -0.46470301232724587,  0.78867513459481275,
             0.33474531727338769, -0.74472292782631933, -0.57735026918962595, -0.94227932796886948, -0.25971420705077974, -0.21132486540518675, -0.00743240513391596, -0.61476523277246198,  0.78867513459481275,
            -0.94227932796886948, -0.25971420705077986, -0.21132486540518680,  0.33474531727338758, -0.74472292782631921, -0.57735026918962595, -0.00743240513391613, -0.61476523277246187,  0.78867513459481287,
            -0.74472292782631966, -0.33474531727338797,  0.57735026918962595,  0.53230171741593768, -0.81975403804892755,  0.21132486540518680, -0.40254520541901595, -0.46470301232724537, -0.78867513459481275,
             0.53230171741593746, -0.81975403804892744,  0.21132486540518686, -0.74472292782631944, -0.33474531727338802,  0.57735026918962595, -0.40254520541901578, -0.46470301232724542, -0.78867513459481287,
             0.33474531727338758, -0.74472292782631944,  0.57735026918962551, -0.94227932796886948, -0.25971420705077991,  0.21132486540518738, -0.00743240513391602, -0.61476523277246153, -0.78867513459481287,
            -0.94227932796886948, -0.25971420705077997,  0.21132486540518733,  0.33474531727338769, -0.74472292782631955,  0.57735026918962551, -0.00743240513391608, -0.61476523277246142, -0.78867513459481275,
            -0.40254520541901584, -0.46470301232724537, -0.78867513459481253, -0.74472292782631966, -0.33474531727338802,  0.57735026918962584,  0.53230171741593746, -0.81975403804892744,  0.21132486540518677,
            -0.74472292782631966, -0.33474531727338802,  0.57735026918962584, -0.40254520541901578, -0.46470301232724542, -0.78867513459481275,  0.53230171741593746, -0.81975403804892744,  0.21132486540518694,
            -0.94227932796886948, -0.25971420705077986,  0.21132486540518733, -0.00743240513391624, -0.61476523277246142, -0.78867513459481287,  0.33474531727338780, -0.74472292782631966,  0.57735026918962540,
            -0.00743240513391633, -0.61476523277246142, -0.78867513459481275, -0.94227932796886926, -0.25971420705077980,  0.21132486540518736,  0.33474531727338780, -0.74472292782631966,  0.57735026918962540,
            -0.40254520541901584, -0.46470301232724587,  0.78867513459481264, -0.74472292782631966, -0.33474531727338763, -0.57735026918962573,  0.53230171741593768, -0.81975403804892744, -0.21132486540518694,
            -0.74472292782631966, -0.33474531727338774, -0.57735026918962540, -0.40254520541901578, -0.46470301232724587,  0.78867513459481253,  0.53230171741593746, -0.81975403804892721, -0.21132486540518711,
            -0.94227932796886948, -0.25971420705077969, -0.21132486540518683, -0.00743240513391610, -0.61476523277246187,  0.78867513459481264,  0.33474531727338769, -0.74472292782631933, -0.57735026918962573,
            -0.00743240513391616, -0.61476523277246187,  0.78867513459481264, -0.94227932796886971, -0.25971420705077969, -0.21132486540518686,  0.33474531727338791, -0.74472292782631933, -0.57735026918962573,
            -0.90148363992456226,  0.21387143172584827,  0.37627949404651062,  0.14283988479115434, -0.67365669613657153,  0.72510924905369145, -0.40856335403608074, -0.70742184475876524, -0.57674111982491172,
             0.14283988479115423, -0.67365669613657153,  0.72510924905369145, -0.90148363992456226,  0.21387143172584827,  0.37627949404651062, -0.40856335403608057, -0.70742184475876535, -0.57674111982491172,
            -0.07071622802089400, -0.36957471874357817,  0.92650620200843870, -0.99200950450550729,  0.12334556714490302, -0.02651441186298359, -0.10448137664308738, -0.92097795757081347, -0.37534416687016475,
            -0.99200950450550707,  0.12334556714490302, -0.02651441186298353, -0.07071622802089389, -0.36957471874357817,  0.92650620200843870, -0.10448137664308754, -0.92097795757081347, -0.37534416687016486,
             0.21387143172584805, -0.90148363992456226,  0.37627949404651045, -0.67365669613657175,  0.14283988479115456,  0.72510924905369167, -0.70742184475876502, -0.40856335403608091, -0.57674111982491172,
            -0.67365669613657164,  0.14283988479115456,  0.72510924905369167,  0.21387143172584805, -0.90148363992456226,  0.37627949404651045, -0.70742184475876502, -0.40856335403608085, -0.57674111982491172,
            -0.36957471874357850, -0.07071622802089367,  0.92650620200843870,  0.12334556714490313, -0.99200950450550729, -0.02651441186298362, -0.92097795757081324, -0.10448137664308760, -0.37534416687016464,
             0.12334556714490313, -0.99200950450550729, -0.02651441186298359, -0.36957471874357850, -0.07071622802089356,  0.92650620200843870, -0.92097795757081324, -0.10448137664308771, -0.37534416687016475,
            -0.92097795757081324, -0.10448137664308760, -0.37534416687016453, -0.36957471874357850, -0.07071622802089361,  0.92650620200843870,  0.12334556714490302, -0.99200950450550729, -0.02651441186298392,
            -0.36957471874357839, -0.07071622802089378,  0.92650620200843847, -0.92097795757081347, -0.10448137664308743, -0.37534416687016464,  0.12334556714490325, -0.99200950450550729, -0.02651441186298359,
            -0.67365669613657186,  0.14283988479115456,  0.72510924905369167, -0.70742184475876502, -0.40856335403608096, -0.57674111982491194,  0.21387143172584827, -0.90148363992456226,  0.37627949404651068,
            -0.70742184475876502, -0.40856335403608091, -0.57674111982491172, -0.67365669613657175,  0.14283988479115445,  0.72510924905369145,  0.21387143172584816, -0.90148363992456226,  0.37627949404651045,
            -0.10448137664308721, -0.92097795757081369, -0.37534416687016475, -0.07071622802089400, -0.36957471874357811,  0.92650620200843870, -0.99200950450550729,  0.12334556714490313, -0.02651441186298381,
            -0.07071622802089400, -0.36957471874357817,  0.92650620200843847, -0.10448137664308732, -0.92097795757081347, -0.37534416687016464, -0.99200950450550729,  0.12334556714490313, -0.02651441186298359,
             0.14283988479115423, -0.67365669613657153,  0.72510924905369145, -0.40856335403608057, -0.70742184475876546, -0.57674111982491194, -0.90148363992456204,  0.21387143172584827,  0.37627949404651079,
            -0.40856335403608068, -0.70742184475876524, -0.57674111982491172,  0.14283988479115445, -0.67365669613657164,  0.72510924905369145, -0.90148363992456226,  0.21387143172584827,  0.37627949404651045,
            -0.99200950450550707,  0.12334556714490313, -0.02651441186298359, -0.10448137664308738, -0.92097795757081369, -0.37534416687016475, -0.07071622802089406, -0.36957471874357817,  0.92650620200843870,
            -0.10448137664308743, -0.92097795757081369, -0.37534416687016475, -0.99200950450550707,  0.12334556714490325, -0.02651441186298359, -0.07071622802089406, -0.36957471874357822,  0.92650620200843847,
            -0.40856335403608057, -0.70742184475876546, -0.57674111982491194, -0.90148363992456226,  0.21387143172584827,  0.37627949404651051,  0.14283988479115423, -0.67365669613657153,  0.72510924905369145,
            -0.90148363992456226,  0.21387143172584827,  0.37627949404651045, -0.40856335403608057, -0.70742184475876546, -0.57674111982491194,  0.14283988479115423, -0.67365669613657131,  0.72510924905369167,
             0.12334556714490313, -0.99200950450550729, -0.02651441186298378, -0.92097795757081347, -0.10448137664308754, -0.37534416687016464, -0.36957471874357833, -0.07071622802089378,  0.92650620200843870,
            -0.92097795757081324, -0.10448137664308760, -0.37534416687016464,  0.12334556714490325, -0.99200950450550729, -0.02651441186298378, -0.36957471874357850, -0.07071622802089372,  0.92650620200843847,
            -0.70742184475876502, -0.40856335403608091, -0.57674111982491194,  0.21387143172584816, -0.90148363992456204,  0.37627949404651045, -0.67365669613657175,  0.14283988479115445,  0.72510924905369167,
             0.21387143172584805, -0.90148363992456204,  0.37627949404651040, -0.70742184475876524, -0.40856335403608091, -0.57674111982491194, -0.67365669613657153,  0.14283988479115445,  0.72510924905369167,
             0.21387143172584849, -0.37627949404651095, -0.90148363992456249, -0.67365669613657153, -0.72510924905369167,  0.14283988479115489, -0.70742184475876557,  0.57674111982491216, -0.40856335403608096,
            -0.67365669613657131, -0.72510924905369145,  0.14283988479115467,  0.21387143172584849, -0.37627949404651084, -0.90148363992456249, -0.70742184475876568,  0.57674111982491194, -0.40856335403608079,
            -0.36957471874357795, -0.92650620200843892, -0.07071622802089361,  0.12334556714490325,  0.02651441186298345, -0.99200950450550773, -0.92097795757081369,  0.37534416687016509, -0.10448137664308732,
             0.12334556714490325,  0.02651441186298342, -0.99200950450550773, -0.36957471874357795, -0.92650620200843892, -0.07071622802089361, -0.92097795757081369,  0.37534416687016520, -0.10448137664308749,
            -0.90148363992456226, -0.37627949404651029,  0.21387143172584849,  0.14283988479115489, -0.72510924905369190, -0.67365669613657175, -0.40856335403608118,  0.57674111982491194, -0.70742184475876546,
             0.14283988479115489, -0.72510924905369190, -0.67365669613657175, -0.90148363992456226, -0.37627949404651018,  0.21387143172584849, -0.40856335403608118,  0.57674111982491194, -0.70742184475876546,
            -0.07071622802089339, -0.92650620200843892, -0.36957471874357828, -0.99200950450550751,  0.02651441186298381,  0.12334556714490336, -0.10448137664308771,  0.37534416687016486, -0.92097795757081369,
            -0.99200950450550751,  0.02651441186298381,  0.12334556714490347, -0.07071622802089339, -0.92650620200843892, -0.36957471874357828, -0.10448137664308771,  0.37534416687016486, -0.92097795757081369,
            -0.10448137664308765,  0.37534416687016453, -0.92097795757081369, -0.07071622802089333, -0.92650620200843914, -0.36957471874357828, -0.99200950450550751,  0.02651441186298420,  0.12334556714490336,
            -0.07071622802089333, -0.92650620200843892, -0.36957471874357833, -0.10448137664308754,  0.37534416687016464, -0.92097795757081391, -0.99200950450550773,  0.02651441186298392,  0.12334556714490347,
             0.14283988479115500, -0.72510924905369190, -0.67365669613657186, -0.40856335403608118,  0.57674111982491194, -0.70742184475876546, -0.90148363992456226, -0.37627949404651040,  0.21387143172584860,
            -0.40856335403608113,  0.57674111982491194, -0.70742184475876546,  0.14283988479115489, -0.72510924905369190, -0.67365669613657186, -0.90148363992456226, -0.37627949404651029,  0.21387143172584860,
            -0.92097795757081413,  0.37534416687016509, -0.10448137664308721, -0.36957471874357789, -0.92650620200843892, -0.07071622802089361,  0.12334556714490325,  0.02651441186298353, -0.99200950450550773,
            -0.36957471874357795, -0.92650620200843870, -0.07071622802089361, -0.92097795757081391,  0.37534416687016509, -0.10448137664308721,  0.12334556714490313,  0.02651441186298342, -0.99200950450550773,
            -0.67365669613657142, -0.72510924905369145,  0.14283988479115467, -0.70742184475876568,  0.57674111982491216, -0.40856335403608079,  0.21387143172584849, -0.37627949404651106, -0.90148363992456249,
            -0.70742184475876568,  0.57674111982491216, -0.40856335403608074, -0.67365669613657153, -0.72510924905369145,  0.14283988479115489,  0.21387143172584872, -0.37627949404651079, -0.90148363992456249,
             0.12334556714490325,  0.02651441186298339, -0.99200950450550751, -0.92097795757081391,  0.37534416687016509, -0.10448137664308738, -0.36957471874357789, -0.92650620200843870, -0.07071622802089367,
            -0.92097795757081391,  0.37534416687016509, -0.10448137664308738,  0.12334556714490325,  0.02651441186298334, -0.99200950450550751, -0.36957471874357789, -0.92650620200843870, -0.07071622802089367,
            -0.70742184475876568,  0.57674111982491216, -0.40856335403608079,  0.21387143172584849, -0.37627949404651068, -0.90148363992456249, -0.67365669613657131, -0.72510924905369167,  0.14283988479115467,
             0.21387143172584849, -0.37627949404651068, -0.90148363992456249, -0.70742184475876568,  0.57674111982491216, -0.40856335403608079, -0.67365669613657131, -0.72510924905369167,  0.14283988479115467,
            -0.99200950450550751,  0.02651441186298412,  0.12334556714490336, -0.10448137664308760,  0.37534416687016475, -0.92097795757081391, -0.07071622802089350, -0.92650620200843914, -0.36957471874357811,
            -0.10448137664308771,  0.37534416687016453, -0.92097795757081369, -0.99200950450550751,  0.02651441186298401,  0.12334556714490347, -0.07071622802089345, -0.92650620200843870, -0.36957471874357828,
            -0.40856335403608113,  0.57674111982491194, -0.70742184475876546, -0.90148363992456226, -0.37627949404651029,  0.21387143172584860,  0.14283988479115467, -0.72510924905369190, -0.67365669613657175,
            -0.90148363992456226, -0.37627949404651029,  0.21387143172584860, -0.40856335403608113,  0.57674111982491216, -0.70742184475876546,  0.14283988479115467, -0.72510924905369212, -0.67365669613657153,
            -0.61476523277246187,  0.78867513459481253, -0.00743240513391605, -0.74472292782631944, -0.57735026918962573,  0.33474531727338774, -0.25971420705077969, -0.21132486540518672, -0.94227932796886926,
            -0.74472292782631944, -0.57735026918962551,  0.33474531727338763, -0.61476523277246176,  0.78867513459481242, -0.00743240513391608, -0.25971420705077974, -0.21132486540518683, -0.94227932796886926,
            -0.81975403804892766, -0.21132486540518711,  0.53230171741593757, -0.46470301232724565,  0.78867513459481253, -0.40254520541901573, -0.33474531727338769, -0.57735026918962540, -0.74472292782631955,
            -0.46470301232724576,  0.78867513459481242, -0.40254520541901573, -0.81975403804892766, -0.21132486540518711,  0.53230171741593768, -0.33474531727338763, -0.57735026918962529, -0.74472292782631966,
            -0.61476523277246153, -0.78867513459481242, -0.00743240513391616, -0.74472292782631966,  0.57735026918962551,  0.33474531727338774, -0.25971420705077980,  0.21132486540518694, -0.94227932796886926,
            -0.74472292782631966,  0.57735026918962551,  0.33474531727338774, -0.61476523277246153, -0.78867513459481242, -0.00743240513391621, -0.25971420705077980,  0.21132486540518700, -0.94227932796886926,
            -0.81975403804892755,  0.21132486540518694,  0.53230171741593757, -0.46470301232724553, -0.78867513459481264, -0.40254520541901562, -0.33474531727338785,  0.57735026918962562, -0.74472292782631955,
            -0.46470301232724553, -0.78867513459481264, -0.40254520541901562, -0.81975403804892766,  0.21132486540518702,  0.53230171741593757, -0.33474531727338774,  0.57735026918962551, -0.74472292782631966,
            -0.33474531727338791,  0.57735026918962551, -0.74472292782631933, -0.81975403804892766,  0.21132486540518708,  0.53230171741593757, -0.46470301232724542, -0.78867513459481242, -0.40254520541901595,
            -0.81975403804892755,  0.21132486540518700,  0.53230171741593757, -0.33474531727338785,  0.57735026918962573, -0.74472292782631955, -0.46470301232724553, -0.78867513459481264, -0.40254520541901573,
            -0.74472292782631944,  0.57735026918962562,  0.33474531727338763, -0.25971420705077980,  0.21132486540518711, -0.94227932796886937, -0.61476523277246165, -0.78867513459481264, -0.00743240513391599,
            -0.25971420705077974,  0.21132486540518711, -0.94227932796886926, -0.74472292782631955,  0.57735026918962551,  0.33474531727338763, -0.61476523277246165, -0.78867513459481264, -0.00743240513391616,
            -0.33474531727338780, -0.57735026918962562, -0.74472292782631944, -0.81975403804892744, -0.21132486540518705,  0.53230171741593757, -0.46470301232724565,  0.78867513459481264, -0.40254520541901584,
            -0.81975403804892755, -0.21132486540518705,  0.53230171741593746, -0.33474531727338774, -0.57735026918962551, -0.74472292782631944, -0.46470301232724565,  0.78867513459481264, -0.40254520541901573,
            -0.74472292782631944, -0.57735026918962551,  0.33474531727338752, -0.25971420705077969, -0.21132486540518702, -0.94227932796886937, -0.61476523277246187,  0.78867513459481242, -0.00743240513391594,
            -0.25971420705077974, -0.21132486540518694, -0.94227932796886937, -0.74472292782631944, -0.57735026918962573,  0.33474531727338774, -0.61476523277246176,  0.78867513459481264, -0.00743240513391610,
            -0.46470301232724570,  0.78867513459481242, -0.40254520541901573, -0.33474531727338769, -0.57735026918962551, -0.74472292782631955, -0.81975403804892744, -0.21132486540518688,  0.53230171741593746,
            -0.33474531727338769, -0.57735026918962551, -0.74472292782631955, -0.46470301232724576,  0.78867513459481242, -0.40254520541901562, -0.81975403804892744, -0.21132486540518694,  0.53230171741593746,
            -0.25971420705077974, -0.21132486540518688, -0.94227932796886937, -0.61476523277246176,  0.78867513459481264, -0.00743240513391624, -0.74472292782631944, -0.57735026918962562,  0.33474531727338774,
            -0.61476523277246176,  0.78867513459481264, -0.00743240513391624, -0.25971420705077969, -0.21132486540518697, -0.94227932796886937, -0.74472292782631955, -0.57735026918962551,  0.33474531727338785,
            -0.46470301232724542, -0.78867513459481253, -0.40254520541901584, -0.33474531727338785,  0.57735026918962573, -0.74472292782631955, -0.81975403804892766,  0.21132486540518677,  0.53230171741593746,
            -0.33474531727338785,  0.57735026918962551, -0.74472292782631944, -0.46470301232724553, -0.78867513459481242, -0.40254520541901578, -0.81975403804892755,  0.21132486540518688,  0.53230171741593746,
            -0.25971420705077963,  0.21132486540518716, -0.94227932796886937, -0.61476523277246165, -0.78867513459481253, -0.00743240513391613, -0.74472292782631966,  0.57735026918962540,  0.33474531727338774,
            -0.61476523277246153, -0.78867513459481242, -0.00743240513391616, -0.25971420705077969,  0.21132486540518716, -0.94227932796886948, -0.74472292782631966,  0.57735026918962529,  0.33474531727338785,
            -0.00743240513391599, -0.78867513459481264, -0.61476523277246165,  0.33474531727338774,  0.57735026918962573, -0.74472292782631966, -0.94227932796886937,  0.21132486540518700, -0.25971420705077963,
             0.33474531727338763,  0.57735026918962551, -0.74472292782631966, -0.00743240513391608, -0.78867513459481264, -0.61476523277246176, -0.94227932796886926,  0.21132486540518711, -0.25971420705077958,
             0.53230171741593757,  0.21132486540518702, -0.81975403804892766, -0.40254520541901573, -0.78867513459481275, -0.46470301232724553, -0.74472292782631966,  0.57735026918962562, -0.33474531727338774,
            -0.40254520541901562, -0.78867513459481264, -0.46470301232724559,  0.53230171741593757,  0.21132486540518705, -0.81975403804892766, -0.74472292782631966,  0.57735026918962562, -0.33474531727338769,
            -0.00743240513391610,  0.78867513459481275, -0.61476523277246176,  0.33474531727338774, -0.57735026918962573, -0.74472292782631955, -0.94227932796886937, -0.21132486540518694, -0.25971420705077963,
             0.33474531727338774, -0.57735026918962584, -0.74472292782631955, -0.00743240513391621,  0.78867513459481287, -0.61476523277246176, -0.94227932796886926, -0.21132486540518700, -0.25971420705077963,
             0.53230171741593757, -0.21132486540518708, -0.81975403804892766, -0.40254520541901573,  0.78867513459481264, -0.46470301232724565, -0.74472292782631966, -0.57735026918962551, -0.33474531727338769,
            -0.40254520541901562,  0.78867513459481264, -0.46470301232724565,  0.53230171741593757, -0.21132486540518716, -0.81975403804892766, -0.74472292782631966, -0.57735026918962551, -0.33474531727338769,
            -0.74472292782631944, -0.57735026918962562, -0.33474531727338774,  0.53230171741593757, -0.21132486540518716, -0.81975403804892766, -0.40254520541901595,  0.78867513459481275, -0.46470301232724553,
             0.53230171741593757, -0.21132486540518722, -0.81975403804892766, -0.74472292782631955, -0.57735026918962584, -0.33474531727338763, -0.40254520541901573,  0.78867513459481287, -0.46470301232724559,
             0.33474531727338774, -0.57735026918962595, -0.74472292782631966, -0.94227932796886948, -0.21132486540518697, -0.25971420705077952, -0.00743240513391599,  0.78867513459481287, -0.61476523277246176,
            -0.94227932796886937, -0.21132486540518697, -0.25971420705077952,  0.33474531727338774, -0.57735026918962595, -0.74472292782631966, -0.00743240513391610,  0.78867513459481287, -0.61476523277246176,
            -0.74472292782631955,  0.57735026918962595, -0.33474531727338769,  0.53230171741593768,  0.21132486540518702, -0.81975403804892777, -0.40254520541901584, -0.78867513459481287, -0.46470301232724553,
             0.53230171741593757,  0.21132486540518702, -0.81975403804892777, -0.74472292782631955,  0.57735026918962595, -0.33474531727338769, -0.40254520541901573, -0.78867513459481287, -0.46470301232724553,
             0.33474531727338763,  0.57735026918962562, -0.74472292782631955, -0.94227932796886937,  0.21132486540518722, -0.25971420705077963, -0.00743240513391594, -0.78867513459481275, -0.61476523277246176,
            -0.94227932796886937,  0.21132486540518727, -0.25971420705077958,  0.33474531727338774,  0.57735026918962584, -0.74472292782631966, -0.00743240513391610, -0.78867513459481287, -0.61476523277246165,
            -0.40254520541901573, -0.78867513459481275, -0.46470301232724553, -0.74472292782631955,  0.57735026918962573, -0.33474531727338774,  0.53230171741593757,  0.21132486540518700, -0.81975403804892766,
            -0.74472292782631955,  0.57735026918962584, -0.33474531727338774, -0.40254520541901562, -0.78867513459481287, -0.46470301232724553,  0.53230171741593746,  0.21132486540518705, -0.81975403804892766,
            -0.94227932796886937,  0.21132486540518713, -0.25971420705077963, -0.00743240513391616, -0.78867513459481264, -0.61476523277246165,  0.33474531727338785,  0.57735026918962551, -0.74472292782631966,
            -0.00743240513391621, -0.78867513459481264, -0.61476523277246165, -0.94227932796886937,  0.21132486540518722, -0.25971420705077963,  0.33474531727338785,  0.57735026918962551, -0.74472292782631966,
            -0.40254520541901584,  0.78867513459481264, -0.46470301232724559, -0.74472292782631955, -0.57735026918962573, -0.33474531727338763,  0.53230171741593757, -0.21132486540518694, -0.81975403804892766,
            -0.74472292782631944, -0.57735026918962551, -0.33474531727338763, -0.40254520541901578,  0.78867513459481264, -0.46470301232724559,  0.53230171741593746, -0.21132486540518705, -0.81975403804892766,
            -0.94227932796886937, -0.21132486540518697, -0.25971420705077963, -0.00743240513391616,  0.78867513459481275, -0.61476523277246176,  0.33474531727338785, -0.57735026918962562, -0.74472292782631955,
            -0.00743240513391624,  0.78867513459481264, -0.61476523277246176, -0.94227932796886937, -0.21132486540518700, -0.25971420705077958,  0.33474531727338785, -0.57735026918962562, -0.74472292782631966,
            -0.90148363992456226,  0.37627949404651062,  0.21387143172584822,  0.14283988479115450,  0.72510924905369178, -0.67365669613657164, -0.40856335403608096, -0.57674111982491183, -0.70742184475876535,
             0.14283988479115428,  0.72510924905369156, -0.67365669613657153, -0.90148363992456226,  0.37627949404651062,  0.21387143172584822, -0.40856335403608079, -0.57674111982491161, -0.70742184475876546,
            -0.07071622802089383,  0.92650620200843881, -0.36957471874357828, -0.99200950450550740, -0.02651441186298367,  0.12334556714490319, -0.10448137664308754, -0.37534416687016470, -0.92097795757081369,
            -0.99200950450550729, -0.02651441186298362,  0.12334556714490319, -0.07071622802089383,  0.92650620200843881, -0.36957471874357828, -0.10448137664308760, -0.37534416687016470, -0.92097795757081358,
             0.21387143172584822,  0.37627949404651062, -0.90148363992456237, -0.67365669613657164,  0.72510924905369178,  0.14283988479115450, -0.70742184475876524, -0.57674111982491183, -0.40856335403608091,
            -0.67365669613657164,  0.72510924905369178,  0.14283988479115450,  0.21387143172584822,  0.37627949404651057, -0.90148363992456237, -0.70742184475876546, -0.57674111982491183, -0.40856335403608091,
            -0.36957471874357839,  0.92650620200843881, -0.07071622802089372,  0.12334556714490330, -0.02651441186298351, -0.99200950450550762, -0.92097795757081358, -0.37534416687016481, -0.10448137664308754,
             0.12334556714490319, -0.02651441186298345, -0.99200950450550751, -0.36957471874357839,  0.92650620200843870, -0.07071622802089361, -0.92097795757081358, -0.37534416687016481, -0.10448137664308765,
            -0.92097795757081358, -0.37534416687016459, -0.10448137664308754, -0.36957471874357839,  0.92650620200843881, -0.07071622802089372,  0.12334556714490308, -0.02651441186298378, -0.99200950450550751,
            -0.36957471874357833,  0.92650620200843870, -0.07071622802089383, -0.92097795757081358, -0.37534416687016470, -0.10448137664308749,  0.12334556714490319, -0.02651441186298351, -0.99200950450550740,
            -0.67365669613657175,  0.72510924905369156,  0.14283988479115461, -0.70742184475876524, -0.57674111982491183, -0.40856335403608102,  0.21387143172584822,  0.37627949404651073, -0.90148363992456237,
            -0.70742184475876524, -0.57674111982491172, -0.40856335403608102, -0.67365669613657175,  0.72510924905369156,  0.14283988479115450,  0.21387143172584810,  0.37627949404651062, -0.90148363992456226,
            -0.10448137664308743, -0.37534416687016459, -0.92097795757081380, -0.07071622802089389,  0.92650620200843881, -0.36957471874357817, -0.99200950450550740, -0.02651441186298378,  0.12334556714490319,
            -0.07071622802089395,  0.92650620200843870, -0.36957471874357817, -0.10448137664308743, -0.37534416687016459, -0.92097795757081369, -0.99200950450550729, -0.02651441186298373,  0.12334556714490308,
             0.14283988479115439,  0.72510924905369156, -0.67365669613657164, -0.40856335403608079, -0.57674111982491183, -0.70742184475876546, -0.90148363992456226,  0.37627949404651073,  0.21387143172584833,
            -0.40856335403608085, -0.57674111982491172, -0.70742184475876535,  0.14283988479115439,  0.72510924905369167, -0.67365669613657175, -0.90148363992456237,  0.37627949404651040,  0.21387143172584822,
            -0.99200950450550740, -0.02651441186298370,  0.12334556714490319, -0.10448137664308754, -0.37534416687016481, -0.92097795757081369, -0.07071622802089389,  0.92650620200843881, -0.36957471874357828,
            -0.10448137664308754, -0.37534416687016481, -0.92097795757081369, -0.99200950450550740, -0.02651441186298364,  0.12334556714490319, -0.07071622802089378,  0.92650620200843881, -0.36957471874357828,
            -0.40856335403608079, -0.57674111982491183, -0.70742184475876546, -0.90148363992456249,  0.37627949404651040,  0.21387143172584844,  0.14283988479115439,  0.72510924905369178, -0.67365669613657164,
            -0.90148363992456237,  0.37627949404651040,  0.21387143172584833, -0.40856335403608079, -0.57674111982491172, -0.70742184475876557,  0.14283988479115439,  0.72510924905369178, -0.67365669613657153,
             0.12334556714490308, -0.02651441186298370, -0.99200950450550740, -0.92097795757081369, -0.37534416687016481, -0.10448137664308754, -0.36957471874357822,  0.92650620200843881, -0.07071622802089383,
            -0.92097795757081347, -0.37534416687016459, -0.10448137664308765,  0.12334556714490308, -0.02651441186298370, -0.99200950450550740, -0.36957471874357839,  0.92650620200843858, -0.07071622802089372,
            -0.70742184475876535, -0.57674111982491183, -0.40856335403608091,  0.21387143172584822,  0.37627949404651062, -0.90148363992456237, -0.67365669613657164,  0.72510924905369167,  0.14283988479115450,
             0.21387143172584810,  0.37627949404651051, -0.90148363992456237, -0.70742184475876535, -0.57674111982491183, -0.40856335403608091, -0.67365669613657153,  0.72510924905369167,  0.14283988479115439,
            -0.37627949404651073,  0.21387143172584849, -0.90148363992456226, -0.72510924905369167, -0.67365669613657153,  0.14283988479115478,  0.57674111982491216, -0.70742184475876546, -0.40856335403608102,
            -0.72510924905369167, -0.67365669613657131,  0.14283988479115467, -0.37627949404651057,  0.21387143172584838, -0.90148363992456226,  0.57674111982491194, -0.70742184475876568, -0.40856335403608085,
            -0.92650620200843870, -0.36957471874357806, -0.07071622802089356,  0.02651441186298362,  0.12334556714490336, -0.99200950450550773,  0.37534416687016486, -0.92097795757081369, -0.10448137664308749,
             0.02651441186298353,  0.12334556714490325, -0.99200950450550751, -0.92650620200843870, -0.36957471874357806, -0.07071622802089356,  0.37534416687016486, -0.92097795757081369, -0.10448137664308760,
            -0.37627949404651051, -0.90148363992456249,  0.21387143172584849, -0.72510924905369190,  0.14283988479115489, -0.67365669613657153,  0.57674111982491194, -0.40856335403608113, -0.70742184475876546,
            -0.72510924905369190,  0.14283988479115489, -0.67365669613657153, -0.37627949404651040, -0.90148363992456249,  0.21387143172584849,  0.57674111982491194, -0.40856335403608113, -0.70742184475876546,
            -0.92650620200843870, -0.07071622802089345, -0.36957471874357817,  0.02651441186298364, -0.99200950450550751,  0.12334556714490336,  0.37534416687016475, -0.10448137664308760, -0.92097795757081391,
             0.02651441186298367, -0.99200950450550751,  0.12334556714490347, -0.92650620200843892, -0.07071622802089339, -0.36957471874357828,  0.37534416687016497, -0.10448137664308771, -0.92097795757081369,
             0.37534416687016453, -0.10448137664308754, -0.92097795757081369, -0.92650620200843892, -0.07071622802089333, -0.36957471874357828,  0.02651441186298392, -0.99200950450550773,  0.12334556714490325,
            -0.92650620200843870, -0.07071622802089345, -0.36957471874357822,  0.37534416687016475, -0.10448137664308749, -0.92097795757081369,  0.02651441186298376, -0.99200950450550751,  0.12334556714490325,
            -0.72510924905369167,  0.14283988479115500, -0.67365669613657175,  0.57674111982491216, -0.40856335403608113, -0.70742184475876546, -0.37627949404651057, -0.90148363992456226,  0.21387143172584849,
             0.57674111982491194, -0.40856335403608113, -0.70742184475876546, -0.72510924905369167,  0.14283988479115500, -0.67365669613657175, -0.37627949404651045, -0.90148363992456226,  0.21387143172584849,
             0.37534416687016486, -0.92097795757081413, -0.10448137664308732, -0.92650620200843892, -0.36957471874357795, -0.07071622802089356,  0.02651441186298370,  0.12334556714490325, -0.99200950450550751,
            -0.92650620200843870, -0.36957471874357795, -0.07071622802089356,  0.37534416687016486, -0.92097795757081413, -0.10448137664308732,  0.02651441186298359,  0.12334556714490325, -0.99200950450550751,
            -0.72510924905369145, -0.67365669613657153,  0.14283988479115467,  0.57674111982491194, -0.70742184475876568, -0.40856335403608079, -0.37627949404651079,  0.21387143172584872, -0.90148363992456226,
             0.57674111982491194, -0.70742184475876568, -0.40856335403608085, -0.72510924905369167, -0.67365669613657153,  0.14283988479115467, -0.37627949404651062,  0.21387143172584849, -0.90148363992456226,
             0.02651441186298362,  0.12334556714490336, -0.99200950450550751,  0.37534416687016497, -0.92097795757081391, -0.10448137664308749, -0.92650620200843870, -0.36957471874357795, -0.07071622802089361,
             0.37534416687016497, -0.92097795757081391, -0.10448137664308749,  0.02651441186298356,  0.12334556714490336, -0.99200950450550751, -0.92650620200843870, -0.36957471874357795, -0.07071622802089361,
             0.57674111982491194, -0.70742184475876568, -0.40856335403608091, -0.37627949404651051,  0.21387143172584849, -0.90148363992456249, -0.72510924905369167, -0.67365669613657153,  0.14283988479115478,
            -0.37627949404651051,  0.21387143172584849, -0.90148363992456249,  0.57674111982491216, -0.70742184475876568, -0.40856335403608079, -0.72510924905369190, -0.67365669613657131,  0.14283988479115467,
             0.02651441186298384, -0.99200950450550751,  0.12334556714490325,  0.37534416687016486, -0.10448137664308754, -0.92097795757081391, -0.92650620200843892, -0.07071622802089356, -0.36957471874357806,
             0.37534416687016475, -0.10448137664308771, -0.92097795757081369,  0.02651441186298370, -0.99200950450550751,  0.12334556714490325, -0.92650620200843870, -0.07071622802089345, -0.36957471874357822,
             0.57674111982491194, -0.40856335403608102, -0.70742184475876546, -0.37627949404651051, -0.90148363992456249,  0.21387143172584860, -0.72510924905369167,  0.14283988479115467, -0.67365669613657153,
            -0.37627949404651040, -0.90148363992456226,  0.21387143172584849,  0.57674111982491194, -0.40856335403608102, -0.70742184475876546, -0.72510924905369167,  0.14283988479115467, -0.67365669613657153,
             0.78867513459481253, -0.61476523277246187, -0.00743240513391599, -0.57735026918962573, -0.74472292782631933,  0.33474531727338774, -0.21132486540518666, -0.25971420705077969, -0.94227932796886937,
            -0.57735026918962562, -0.74472292782631933,  0.33474531727338763,  0.78867513459481253, -0.61476523277246187, -0.00743240513391608, -0.21132486540518683, -0.25971420705077974, -0.94227932796886926,
            -0.21132486540518713, -0.81975403804892744,  0.53230171741593757,  0.78867513459481264, -0.46470301232724576, -0.40254520541901562, -0.57735026918962551, -0.33474531727338763, -0.74472292782631966,
             0.78867513459481264, -0.46470301232724587, -0.40254520541901562, -0.21132486540518713, -0.81975403804892744,  0.53230171741593768, -0.57735026918962540, -0.33474531727338763, -0.74472292782631966,
            -0.78867513459481264, -0.61476523277246142, -0.00743240513391624,  0.57735026918962573, -0.74472292782631966,  0.33474531727338785,  0.21132486540518700, -0.25971420705077991, -0.94227932796886937,
             0.57735026918962562, -0.74472292782631966,  0.33474531727338785, -0.78867513459481275, -0.61476523277246142, -0.00743240513391630,  0.21132486540518711, -0.25971420705077986, -0.94227932796886926,
             0.21132486540518700, -0.81975403804892766,  0.53230171741593757, -0.78867513459481275, -0.46470301232724542, -0.40254520541901573,  0.57735026918962573, -0.33474531727338785, -0.74472292782631944,
            -0.78867513459481264, -0.46470301232724548, -0.40254520541901567,  0.21132486540518705, -0.81975403804892766,  0.53230171741593757,  0.57735026918962562, -0.33474531727338785, -0.74472292782631966,
             0.57735026918962573, -0.33474531727338808, -0.74472292782631933,  0.21132486540518711, -0.81975403804892766,  0.53230171741593768, -0.78867513459481264, -0.46470301232724531, -0.40254520541901595,
             0.21132486540518705, -0.81975403804892766,  0.53230171741593757,  0.57735026918962595, -0.33474531727338791, -0.74472292782631944, -0.78867513459481275, -0.46470301232724542, -0.40254520541901573,
             0.57735026918962573, -0.74472292782631955,  0.33474531727338763,  0.21132486540518713, -0.25971420705077980, -0.94227932796886937, -0.78867513459481275, -0.61476523277246153, -0.00743240513391605,
             0.21132486540518722, -0.25971420705077980, -0.94227932796886926,  0.57735026918962562, -0.74472292782631955,  0.33474531727338763, -0.78867513459481275, -0.61476523277246153, -0.00743240513391616,
            -0.57735026918962573, -0.33474531727338774, -0.74472292782631944, -0.21132486540518708, -0.81975403804892744,  0.53230171741593757,  0.78867513459481275, -0.46470301232724576, -0.40254520541901584,
            -0.21132486540518716, -0.81975403804892744,  0.53230171741593746, -0.57735026918962562, -0.33474531727338774, -0.74472292782631944,  0.78867513459481275, -0.46470301232724576, -0.40254520541901573,
            -0.57735026918962573, -0.74472292782631921,  0.33474531727338752, -0.21132486540518705, -0.25971420705077969, -0.94227932796886948,  0.78867513459481264, -0.61476523277246198, -0.00743240513391588,
            -0.21132486540518700, -0.25971420705077969, -0.94227932796886937, -0.57735026918962595, -0.74472292782631944,  0.33474531727338763,  0.78867513459481275, -0.61476523277246187, -0.00743240513391610,
             0.78867513459481264, -0.46470301232724587, -0.40254520541901562, -0.57735026918962573, -0.33474531727338763, -0.74472292782631966, -0.21132486540518694, -0.81975403804892744,  0.53230171741593757,
            -0.57735026918962562, -0.33474531727338763, -0.74472292782631966,  0.78867513459481275, -0.46470301232724587, -0.40254520541901556, -0.21132486540518705, -0.81975403804892744,  0.53230171741593746,
            -0.21132486540518694, -0.25971420705077969, -0.94227932796886937,  0.78867513459481275, -0.61476523277246187, -0.00743240513391613, -0.57735026918962573, -0.74472292782631944,  0.33474531727338763,
             0.78867513459481264, -0.61476523277246176, -0.00743240513391619, -0.21132486540518700, -0.25971420705077969, -0.94227932796886937, -0.57735026918962562, -0.74472292782631944,  0.33474531727338785,
            -0.78867513459481253, -0.46470301232724542, -0.40254520541901584,  0.57735026918962573, -0.33474531727338797, -0.74472292782631955,  0.21132486540518672, -0.81975403804892766,  0.53230171741593757,
             0.57735026918962562, -0.33474531727338797, -0.74472292782631944, -0.78867513459481253, -0.46470301232724542, -0.40254520541901578,  0.21132486540518688, -0.81975403804892755,  0.53230171741593746,
             0.21132486540518719, -0.25971420705077980, -0.94227932796886937, -0.78867513459481264, -0.61476523277246153, -0.00743240513391624,  0.57735026918962551, -0.74472292782631966,  0.33474531727338785,
            -0.78867513459481264, -0.61476523277246142, -0.00743240513391627,  0.21132486540518719, -0.25971420705077980, -0.94227932796886948,  0.57735026918962540, -0.74472292782631966,  0.33474531727338785,
            -0.78867513459481253, -0.00743240513391608, -0.61476523277246153,  0.57735026918962562,  0.33474531727338763, -0.74472292782631966,  0.21132486540518694, -0.94227932796886926, -0.25971420705077980,
             0.57735026918962551,  0.33474531727338763, -0.74472292782631966, -0.78867513459481253, -0.00743240513391613, -0.61476523277246153,  0.21132486540518711, -0.94227932796886926, -0.25971420705077980,
             0.21132486540518688,  0.53230171741593746, -0.81975403804892755, -0.78867513459481264, -0.40254520541901562, -0.46470301232724548,  0.57735026918962573, -0.74472292782631955, -0.33474531727338785,
            -0.78867513459481264, -0.40254520541901562, -0.46470301232724553,  0.21132486540518697,  0.53230171741593757, -0.81975403804892755,  0.57735026918962573, -0.74472292782631955, -0.33474531727338785,
             0.78867513459481264, -0.00743240513391619, -0.61476523277246176, -0.57735026918962584,  0.33474531727338763, -0.74472292782631944, -0.21132486540518683, -0.94227932796886926, -0.25971420705077974,
            -0.57735026918962584,  0.33474531727338763, -0.74472292782631944,  0.78867513459481264, -0.00743240513391624, -0.61476523277246176, -0.21132486540518683, -0.94227932796886915, -0.25971420705077974,
            -0.21132486540518722,  0.53230171741593746, -0.81975403804892744,  0.78867513459481275, -0.40254520541901562, -0.46470301232724570, -0.57735026918962551, -0.74472292782631955, -0.33474531727338774,
             0.78867513459481275, -0.40254520541901556, -0.46470301232724576, -0.21132486540518722,  0.53230171741593746, -0.81975403804892744, -0.57735026918962551, -0.74472292782631955, -0.33474531727338769,
            -0.57735026918962551, -0.74472292782631933, -0.33474531727338774, -0.21132486540518733,  0.53230171741593746, -0.81975403804892755,  0.78867513459481287, -0.40254520541901584, -0.46470301232724565,
            -0.21132486540518716,  0.53230171741593746, -0.81975403804892744, -0.57735026918962551, -0.74472292782631944, -0.33474531727338780,  0.78867513459481264, -0.40254520541901573, -0.46470301232724565,
            -0.57735026918962584,  0.33474531727338763, -0.74472292782631944, -0.21132486540518686, -0.94227932796886937, -0.25971420705077969,  0.78867513459481264, -0.00743240513391610, -0.61476523277246187,
            -0.21132486540518691, -0.94227932796886926, -0.25971420705077969, -0.57735026918962584,  0.33474531727338763, -0.74472292782631944,  0.78867513459481264, -0.00743240513391616, -0.61476523277246187,
             0.57735026918962584, -0.74472292782631944, -0.33474531727338791,  0.21132486540518691,  0.53230171741593757, -0.81975403804892766, -0.78867513459481264, -0.40254520541901573, -0.46470301232724542,
             0.21132486540518697,  0.53230171741593746, -0.81975403804892766,  0.57735026918962584, -0.74472292782631944, -0.33474531727338785, -0.78867513459481264, -0.40254520541901573, -0.46470301232724548,
             0.57735026918962551,  0.33474531727338752, -0.74472292782631955,  0.21132486540518738, -0.94227932796886926, -0.25971420705077974, -0.78867513459481287, -0.00743240513391605, -0.61476523277246165,
             0.21132486540518722, -0.94227932796886926, -0.25971420705077986,  0.57735026918962551,  0.33474531727338763, -0.74472292782631944, -0.78867513459481264, -0.00743240513391610, -0.61476523277246165,
            -0.78867513459481264, -0.40254520541901567, -0.46470301232724553,  0.57735026918962584, -0.74472292782631944, -0.33474531727338785,  0.21132486540518688,  0.53230171741593746, -0.81975403804892755,
             0.57735026918962584, -0.74472292782631944, -0.33474531727338785, -0.78867513459481264, -0.40254520541901562, -0.46470301232724553,  0.21132486540518688,  0.53230171741593735, -0.81975403804892755,
             0.21132486540518727, -0.94227932796886926, -0.25971420705077980, -0.78867513459481275, -0.00743240513391627, -0.61476523277246153,  0.57735026918962551,  0.33474531727338774, -0.74472292782631955,
            -0.78867513459481275, -0.00743240513391630, -0.61476523277246153,  0.21132486540518727, -0.94227932796886926, -0.25971420705077980,  0.57735026918962551,  0.33474531727338774, -0.74472292782631966,
             0.78867513459481253, -0.40254520541901578, -0.46470301232724570, -0.57735026918962562, -0.74472292782631944, -0.33474531727338769, -0.21132486540518688,  0.53230171741593746, -0.81975403804892744,
            -0.57735026918962551, -0.74472292782631944, -0.33474531727338769,  0.78867513459481253, -0.40254520541901573, -0.46470301232724576, -0.21132486540518705,  0.53230171741593746, -0.81975403804892744,
            -0.21132486540518683, -0.94227932796886926, -0.25971420705077974,  0.78867513459481264, -0.00743240513391621, -0.61476523277246176, -0.57735026918962573,  0.33474531727338774, -0.74472292782631944,
             0.78867513459481264, -0.00743240513391624, -0.61476523277246176, -0.21132486540518691, -0.94227932796886937, -0.25971420705077974, -0.57735026918962573,  0.33474531727338774, -0.74472292782631944,
             0.37627949404651051, -0.90148363992456215,  0.21387143172584822,  0.72510924905369167,  0.14283988479115428, -0.67365669613657153, -0.57674111982491172, -0.40856335403608096, -0.70742184475876546,
             0.72510924905369156,  0.14283988479115417, -0.67365669613657142,  0.37627949404651051, -0.90148363992456215,  0.21387143172584822, -0.57674111982491161, -0.40856335403608074, -0.70742184475876546,
             0.92650620200843881, -0.07071622802089400, -0.36957471874357806, -0.02651441186298364, -0.99200950450550718,  0.12334556714490297, -0.37534416687016459, -0.10448137664308754, -0.92097795757081369,
            -0.02651441186298359, -0.99200950450550718,  0.12334556714490297,  0.92650620200843881, -0.07071622802089395, -0.36957471874357806, -0.37534416687016481, -0.10448137664308765, -0.92097795757081369,
             0.37627949404651051,  0.21387143172584799, -0.90148363992456215,  0.72510924905369167, -0.67365669613657164,  0.14283988479115450, -0.57674111982491172, -0.70742184475876513, -0.40856335403608113,
             0.72510924905369167, -0.67365669613657153,  0.14283988479115450,  0.37627949404651051,  0.21387143172584799, -0.90148363992456215, -0.57674111982491172, -0.70742184475876513, -0.40856335403608113,
             0.92650620200843870, -0.36957471874357845, -0.07071622802089356, -0.02651441186298351,  0.12334556714490308, -0.99200950450550740, -0.37534416687016470, -0.92097795757081335, -0.10448137664308776,
            -0.02651441186298342,  0.12334556714490308, -0.99200950450550740,  0.92650620200843870, -0.36957471874357850, -0.07071622802089350, -0.37534416687016470, -0.92097795757081335, -0.10448137664308782,
            -0.37534416687016448, -0.92097795757081335, -0.10448137664308771,  0.92650620200843881, -0.36957471874357839, -0.07071622802089356, -0.02651441186298384,  0.12334556714490297, -0.99200950450550740,
             0.92650620200843881, -0.36957471874357839, -0.07071622802089356, -0.37534416687016459, -0.92097795757081358, -0.10448137664308765, -0.02651441186298367,  0.12334556714490319, -0.99200950450550751,
             0.72510924905369145, -0.67365669613657175,  0.14283988479115461, -0.57674111982491183, -0.70742184475876524, -0.40856335403608113,  0.37627949404651073,  0.21387143172584810, -0.90148363992456215,
            -0.57674111982491161, -0.70742184475876524, -0.40856335403608107,  0.72510924905369145, -0.67365669613657175,  0.14283988479115439,  0.37627949404651062,  0.21387143172584810, -0.90148363992456215,
            -0.37534416687016448, -0.10448137664308743, -0.92097795757081380,  0.92650620200843881, -0.07071622802089395, -0.36957471874357806, -0.02651441186298384, -0.99200950450550729,  0.12334556714490297,
             0.92650620200843858, -0.07071622802089395, -0.36957471874357811, -0.37534416687016448, -0.10448137664308743, -0.92097795757081358, -0.02651441186298367, -0.99200950450550729,  0.12334556714490297,
             0.72510924905369145,  0.14283988479115417, -0.67365669613657153, -0.57674111982491183, -0.40856335403608079, -0.70742184475876568,  0.37627949404651073, -0.90148363992456215,  0.21387143172584822,
            -0.57674111982491183, -0.40856335403608079, -0.70742184475876568,  0.72510924905369156,  0.14283988479115439, -0.67365669613657153,  0.37627949404651062, -0.90148363992456237,  0.21387143172584833,
            -0.02651441186298356, -0.99200950450550718,  0.12334556714490297, -0.37534416687016470, -0.10448137664308754, -0.92097795757081369,  0.92650620200843870, -0.07071622802089406, -0.36957471874357806,
            -0.37534416687016470, -0.10448137664308760, -0.92097795757081369, -0.02651441186298356, -0.99200950450550718,  0.12334556714490297,  0.92650620200843870, -0.07071622802089406, -0.36957471874357806,
            -0.57674111982491172, -0.40856335403608074, -0.70742184475876568,  0.37627949404651040, -0.90148363992456226,  0.21387143172584822,  0.72510924905369167,  0.14283988479115417, -0.67365669613657142,
             0.37627949404651034, -0.90148363992456226,  0.21387143172584822, -0.57674111982491172, -0.40856335403608068, -0.70742184475876568,  0.72510924905369167,  0.14283988479115417, -0.67365669613657131,
            -0.02651441186298362,  0.12334556714490297, -0.99200950450550740, -0.37534416687016470, -0.92097795757081347, -0.10448137664308771,  0.92650620200843870, -0.36957471874357822, -0.07071622802089372,
            -0.37534416687016459, -0.92097795757081335, -0.10448137664308776, -0.02651441186298359,  0.12334556714490297, -0.99200950450550740,  0.92650620200843858, -0.36957471874357845, -0.07071622802089372,
            -0.57674111982491183, -0.70742184475876524, -0.40856335403608113,  0.37627949404651057,  0.21387143172584799, -0.90148363992456215,  0.72510924905369156, -0.67365669613657164,  0.14283988479115450,
             0.37627949404651051,  0.21387143172584799, -0.90148363992456215, -0.57674111982491183, -0.70742184475876524, -0.40856335403608113,  0.72510924905369178, -0.67365669613657153,  0.14283988479115450
        });
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
