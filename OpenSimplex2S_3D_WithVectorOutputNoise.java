/**
 * OpenSimplex2S 3D with 2D and 3D unit-length-bound vector-output evaluator endpoints.
 */

public class OpenSimplex2S_3D_WithVectorOutputNoise {

    private static final long PRIME_X = 0x5205402B9270C86FL;
    private static final long PRIME_Y = 0x598CD327003817B5L;
    private static final long PRIME_Z = 0x5BCC226E9FA0BACBL;
    private static final long HASH_MULTIPLIER = 0x53A3F72DEEC546F5L;
    private static final long SEED_FLIP_3D = -0x52D547B2E96ED629L;

    private static final double ROOT3OVER3 = 0.577350269189626;
    private static final double FALLBACK_ROTATE3 = 2.0 / 3.0;
    private static final double ROTATE3_ORTHOGONALIZER = -0.21132486540518713;

    private static final int N_GRADS_3D_EXPONENT = 8;
    private static final int N_GRADS_3D = 1 << N_GRADS_3D_EXPONENT;
    private static final int GRAD_HASH_MASK_3D = (N_GRADS_3D - 1) << 2;

    private static final int N_OUTPUT_VEC2 = 256;
    private static final int VEC2_HASH_SHIFT_3D = N_GRADS_3D_EXPONENT + 1;
    private static final int VEC2_HASH_MASK_3D = (N_OUTPUT_VEC2 - 1) << 1;

    private static final int N_OUTPUT_VEC3 = 256;
    private static final int VEC3_HASH_SHIFT_3D = N_GRADS_3D_EXPONENT;
    private static final int VEC3_HASH_MASK_3D = (N_OUTPUT_VEC3 - 1) << 2;

    private static final double NORMALIZER_3D = 0.2781926117527186;

    private static final float RSQUARED_3D = 3.0f / 4.0f;

    /*
     * Ordinary 3D noise
     */

    /**
     * 3D OpenSimplex2S/SuperSimplex noise, with better visual isotropy in (X, Z).
     * Recommended for 3D terrain and time-varied animations.
     * The Y coordinate should always be the "different" coordinate in whatever your use case is.
     * If Y is vertical in world coordinates, call noise3_ImproveXZ(x, Y, z).
     * If Z is vertical in world coordinates, call noise3_ImproveXZ(x, Z, y) or use noise3_ImproveXY.
     * For a time varied animation, call noise3_ImproveXZ(x, T, y) or use noise3_ImproveXY.
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
     * Generate overlapping cubic lattices for 3D OpenSimplex2S noise.
     * Lookup table implementation inspired by DigitalShadow.
     * It was actually faster to narrow down the points in the loop itself,
     * than to build up the index with enough info to isolate 8 points.
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

    /*
     * Unit-Vec3-Output 3D noise
     */

    /**
     * 3D OpenSimplex2S/SuperSimplex noise, with better visual isotropy in (X, Z).
     * Unit-bounded-vec2-output version.
     */
    public static Vector2f vec2Noise3_ImproveXZ(long seed, double x, double y, double z) {

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
        return vec2Noise3_UnrotatedBase(seed, xr, yr, zr);
    }

    /**
     * Generate overlapping cubic lattices for 3D OpenSimplex2 noise.
     * Unit-bounded-vec2-output version.
     */
    private static Vector2f vec2Noise3_UnrotatedBase(long seed, double xr, double yr, double zr) {

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
        int h0 = v2Hash(seed, xrbp + (xNMask & PRIME_X), yrbp + (yNMask & PRIME_Y), zrbp + (zNMask & PRIME_Z));
        float v0 = (a0 * a0) * (a0 * a0) * grad(h0, x0, y0, z0);
        int vi0 = (h0 >> VEC2_HASH_SHIFT_3D) & VEC2_HASH_MASK_3D;
        float vx = v0 * OUTPUT_VEC2[vi0 | 0];
        float vy = v0 * OUTPUT_VEC2[vi0 | 1];

        // Second vertex.
        float x1 = xi - 0.5f;
        float y1 = yi - 0.5f;
        float z1 = zi - 0.5f;
        float a1 = RSQUARED_3D - x1 * x1 - y1 * y1 - z1 * z1;
        int h1 = v2Hash(seed2, xrbp + PRIME_X, yrbp + PRIME_Y, zrbp + PRIME_Z);
        float v1 = (a1 * a1) * (a1 * a1) * grad(h1, x1, y1, z1);
        int vi1 = (h1 >> VEC2_HASH_SHIFT_3D) & VEC2_HASH_MASK_3D;
        vx += v1 * OUTPUT_VEC2[vi1 | 0];
        vy += v1 * OUTPUT_VEC2[vi1 | 1];

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
            int h2 = v2Hash(seed, xrbp + (~xNMask & PRIME_X), yrbp + (yNMask & PRIME_Y), zrbp + (zNMask & PRIME_Z));
            float v2 = (a2 * a2) * (a2 * a2) * grad(h2, x2, y2, z2);
            int vi2 = (h2 >> VEC2_HASH_SHIFT_3D) & VEC2_HASH_MASK_3D;
            vx += v2 * OUTPUT_VEC2[vi2 | 0];
            vy += v2 * OUTPUT_VEC2[vi2 | 1];
        }
        else
        {
            float a3 = yAFlipMask0 + zAFlipMask0 + a0;
            if (a3 > 0) {
                float x3 = x0;
                float y3 = y0 - (yNMask | 1);
                float z3 = z0 - (zNMask | 1);
                int h3 = v2Hash(seed, xrbp + (xNMask & PRIME_X), yrbp + (~yNMask & PRIME_Y), zrbp + (~zNMask & PRIME_Z));
                float v3 = (a3 * a3) * (a3 * a3) * grad(h3, x3, y3, z3);
                int vi3 = (h3 >> VEC2_HASH_SHIFT_3D) & VEC2_HASH_MASK_3D;
                vx += v3 * OUTPUT_VEC2[vi3 | 0];
                vy += v3 * OUTPUT_VEC2[vi3 | 1];
            }

            float a4 = xAFlipMask1 + a1;
            if (a4 > 0) {
                float x4 = (xNMask | 1) + x1;
                float y4 = y1;
                float z4 = z1;
                int h4 = v2Hash(seed2, xrbp + (xNMask & (PRIME_X * 2)), yrbp + PRIME_Y, zrbp + PRIME_Z);
                float v4 = (a4 * a4) * (a4 * a4) * grad(h4, x4, y4, z4);
                int vi4 = (h4 >> VEC2_HASH_SHIFT_3D) & VEC2_HASH_MASK_3D;
                vx += v4 * OUTPUT_VEC2[vi4 | 0];
                vy += v4 * OUTPUT_VEC2[vi4 | 1];
                skip5 = true;
            }
        }

        boolean skip9 = false;
        float a6 = yAFlipMask0 + a0;
        if (a6 > 0) {
            float x6 = x0;
            float y6 = y0 - (yNMask | 1);
            float z6 = z0;
            int h6 = v2Hash(seed, xrbp + (xNMask & PRIME_X), yrbp + (~yNMask & PRIME_Y), zrbp + (zNMask & PRIME_Z));
            float v6 = (a6 * a6) * (a6 * a6) * grad(h6, x6, y6, z6);
            int vi6 = (h6 >> VEC2_HASH_SHIFT_3D) & VEC2_HASH_MASK_3D;
            vx += v6 * OUTPUT_VEC2[vi6 | 0];
            vy += v6 * OUTPUT_VEC2[vi6 | 1];
        }
        else
        {
            float a7 = xAFlipMask0 + zAFlipMask0 + a0;
            if (a7 > 0) {
                float x7 = x0 - (xNMask | 1);
                float y7 = y0;
                float z7 = z0 - (zNMask | 1);
                int h7 = v2Hash(seed, xrbp + (~xNMask & PRIME_X), yrbp + (yNMask & PRIME_Y), zrbp + (~zNMask & PRIME_Z));
                float v7 = (a7 * a7) * (a7 * a7) * grad(h7, x7, y7, z7);
                int vi7 = (h7 >> VEC2_HASH_SHIFT_3D) & VEC2_HASH_MASK_3D;
                vx += v7 * OUTPUT_VEC2[vi7 | 0];
                vy += v7 * OUTPUT_VEC2[vi7 | 1];
            }

            float a8 = yAFlipMask1 + a1;
            if (a8 > 0) {
                float x8 = x1;
                float y8 = (yNMask | 1) + y1;
                float z8 = z1;
                int h8 = v2Hash(seed2, xrbp + PRIME_X, yrbp + (yNMask & (PRIME_Y << 1)), zrbp + PRIME_Z);
                float v8 = (a8 * a8) * (a8 * a8) * grad(h8, x8, y8, z8);
                int vi8 = (h8 >> VEC2_HASH_SHIFT_3D) & VEC2_HASH_MASK_3D;
                vx += v8 * OUTPUT_VEC2[vi8 | 0];
                vy += v8 * OUTPUT_VEC2[vi8 | 1];
                skip9 = true;
            }
        }

        boolean skipD = false;
        float aA = zAFlipMask0 + a0;
        if (aA > 0) {
            float xA = x0;
            float yA = y0;
            float zA = z0 - (zNMask | 1);
            int hA = v2Hash(seed, xrbp + (xNMask & PRIME_X), yrbp + (yNMask & PRIME_Y), zrbp + (~zNMask & PRIME_Z));
            float vA = (aA * aA) * (aA * aA) * grad(hA, xA, yA, zA);
            int viA = (hA >> VEC2_HASH_SHIFT_3D) & VEC2_HASH_MASK_3D;
            vx += vA * OUTPUT_VEC2[viA | 0];
            vy += vA * OUTPUT_VEC2[viA | 1];
        }
        else
        {
            float aB = xAFlipMask0 + yAFlipMask0 + a0;
            if (aB > 0) {
                float xB = x0 - (xNMask | 1);
                float yB = y0 - (yNMask | 1);
                float zB = z0;
                int hB = v2Hash(seed, xrbp + (~xNMask & PRIME_X), yrbp + (~yNMask & PRIME_Y), zrbp + (zNMask & PRIME_Z));
                float vB = (aB * aB) * (aB * aB) * grad(hB, xB, yB, zB);
                int viB = (hB >> VEC2_HASH_SHIFT_3D) & VEC2_HASH_MASK_3D;
                vx += vB * OUTPUT_VEC2[viB | 0];
                vy += vB * OUTPUT_VEC2[viB | 1];
            }

            float aC = zAFlipMask1 + a1;
            if (aC > 0) {
                float xC = x1;
                float yC = y1;
                float zC = (zNMask | 1) + z1;
                int hC = v2Hash(seed2, xrbp + PRIME_X, yrbp + PRIME_Y, zrbp + (zNMask & (PRIME_Z << 1)));
                float vC = (aC * aC) * (aC * aC) * grad(hC, xC, yC, zC);
                int viC = (hC >> VEC2_HASH_SHIFT_3D) & VEC2_HASH_MASK_3D;
                vx += vC * OUTPUT_VEC2[viC | 0];
                vy += vC * OUTPUT_VEC2[viC | 1];
                skipD = true;
            }
        }

        if (!skip5) {
            float a5 = yAFlipMask1 + zAFlipMask1 + a1;
            if (a5 > 0) {
                float x5 = x1;
                float y5 = (yNMask | 1) + y1;
                float z5 = (zNMask | 1) + z1;
                int h5 = v2Hash(seed2, xrbp + PRIME_X, yrbp + (yNMask & (PRIME_Y << 1)), zrbp + (zNMask & (PRIME_Z << 1)));
                float v5 = (a5 * a5) * (a5 * a5) * grad(h5, x5, y5, z5);
                int vi5 = (h5 >> VEC2_HASH_SHIFT_3D) & VEC2_HASH_MASK_3D;
                vx += v5 * OUTPUT_VEC2[vi5 | 0];
                vy += v5 * OUTPUT_VEC2[vi5 | 1];
            }
        }

        if (!skip9) {
            float a9 = xAFlipMask1 + zAFlipMask1 + a1;
            if (a9 > 0) {
                float x9 = (xNMask | 1) + x1;
                float y9 = y1;
                float z9 = (zNMask | 1) + z1;
                int h9 = v2Hash(seed2, xrbp + (xNMask & (PRIME_X * 2)), yrbp + PRIME_Y, zrbp + (zNMask & (PRIME_Z << 1)));
                float v9 = (a9 * a9) * (a9 * a9) * grad(h9, x9, y9, z9);
                int vi9 = (h9 >> VEC2_HASH_SHIFT_3D) & VEC2_HASH_MASK_3D;
                vx += v9 * OUTPUT_VEC2[vi9 | 0];
                vy += v9 * OUTPUT_VEC2[vi9 | 1];
            }
        }

        if (!skipD) {
            float aD = xAFlipMask1 + yAFlipMask1 + a1;
            if (aD > 0) {
                float xD = (xNMask | 1) + x1;
                float yD = (yNMask | 1) + y1;
                float zD = z1;
                int hD = v2Hash(seed2, xrbp + (xNMask & (PRIME_X << 1)), yrbp + (yNMask & (PRIME_Y << 1)), zrbp + PRIME_Z);
                float vD = (aD * aD) * (aD * aD) * grad(hD, xD, yD, zD);
                int viD = (hD >> VEC2_HASH_SHIFT_3D) & VEC2_HASH_MASK_3D;
                vx += vD * OUTPUT_VEC2[viD | 0];
                vy += vD * OUTPUT_VEC2[viD | 1];
            }
        }

        return new Vector2f(vx, vy);
    }

    /*
     * Unit-Vec3-Output 3D noise
     */

    /**
     * 3D OpenSimplex2S/SuperSimplex noise, with better visual isotropy in (X, Z).
     * Unit-bounded-vec3-output version.
     */
    public static Vector3f vec3Noise3_ImproveXZ(long seed, double x, double y, double z) {

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
        return vec3Noise3_UnrotatedBase(seed, xr, yr, zr);
    }

    /**
     * Generate overlapping cubic lattices for 3D OpenSimplex2 noise.
     * Unit-bounded-vec3-output version.
     */
    private static Vector3f vec3Noise3_UnrotatedBase(long seed, double xr, double yr, double zr) {

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
        int h0 = v3Hash(seed, xrbp + (xNMask & PRIME_X), yrbp + (yNMask & PRIME_Y), zrbp + (zNMask & PRIME_Z));
        float v0 = (a0 * a0) * (a0 * a0) * grad(h0, x0, y0, z0);
        int vi0 = (h0 >> VEC3_HASH_SHIFT_3D) & VEC3_HASH_MASK_3D;
        float vx = v0 * OUTPUT_VEC3[vi0 | 0];
        float vy = v0 * OUTPUT_VEC3[vi0 | 1];
        float vz = v0 * OUTPUT_VEC3[vi0 | 2];

        // Second vertex.
        float x1 = xi - 0.5f;
        float y1 = yi - 0.5f;
        float z1 = zi - 0.5f;
        float a1 = RSQUARED_3D - x1 * x1 - y1 * y1 - z1 * z1;
        int h1 = v3Hash(seed2, xrbp + PRIME_X, yrbp + PRIME_Y, zrbp + PRIME_Z);
        float v1 = (a1 * a1) * (a1 * a1) * grad(h1, x1, y1, z1);
        int vi1 = (h1 >> VEC3_HASH_SHIFT_3D) & VEC3_HASH_MASK_3D;
        vx += v1 * OUTPUT_VEC3[vi1 | 0];
        vy += v1 * OUTPUT_VEC3[vi1 | 1];
        vz += v1 * OUTPUT_VEC3[vi1 | 2];

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
            int h2 = v3Hash(seed, xrbp + (~xNMask & PRIME_X), yrbp + (yNMask & PRIME_Y), zrbp + (zNMask & PRIME_Z));
            float v2 = (a2 * a2) * (a2 * a2) * grad(h2, x2, y2, z2);
            int vi2 = (h2 >> VEC3_HASH_SHIFT_3D) & VEC3_HASH_MASK_3D;
            vx += v2 * OUTPUT_VEC3[vi2 | 0];
            vy += v2 * OUTPUT_VEC3[vi2 | 1];
            vz += v2 * OUTPUT_VEC3[vi2 | 2];
        }
        else
        {
            float a3 = yAFlipMask0 + zAFlipMask0 + a0;
            if (a3 > 0) {
                float x3 = x0;
                float y3 = y0 - (yNMask | 1);
                float z3 = z0 - (zNMask | 1);
                int h3 = v3Hash(seed, xrbp + (xNMask & PRIME_X), yrbp + (~yNMask & PRIME_Y), zrbp + (~zNMask & PRIME_Z));
                float v3 = (a3 * a3) * (a3 * a3) * grad(h3, x3, y3, z3);
                int vi3 = (h3 >> VEC3_HASH_SHIFT_3D) & VEC3_HASH_MASK_3D;
                vx += v3 * OUTPUT_VEC3[vi3 | 0];
                vy += v3 * OUTPUT_VEC3[vi3 | 1];
                vz += v3 * OUTPUT_VEC3[vi3 | 2];
            }

            float a4 = xAFlipMask1 + a1;
            if (a4 > 0) {
                float x4 = (xNMask | 1) + x1;
                float y4 = y1;
                float z4 = z1;
                int h4 = v3Hash(seed2, xrbp + (xNMask & (PRIME_X * 2)), yrbp + PRIME_Y, zrbp + PRIME_Z);
                float v4 = (a4 * a4) * (a4 * a4) * grad(h4, x4, y4, z4);
                int vi4 = (h4 >> VEC3_HASH_SHIFT_3D) & VEC3_HASH_MASK_3D;
                vx += v4 * OUTPUT_VEC3[vi4 | 0];
                vy += v4 * OUTPUT_VEC3[vi4 | 1];
                vz += v4 * OUTPUT_VEC3[vi4 | 2];
                skip5 = true;
            }
        }

        boolean skip9 = false;
        float a6 = yAFlipMask0 + a0;
        if (a6 > 0) {
            float x6 = x0;
            float y6 = y0 - (yNMask | 1);
            float z6 = z0;
            int h6 = v3Hash(seed, xrbp + (xNMask & PRIME_X), yrbp + (~yNMask & PRIME_Y), zrbp + (zNMask & PRIME_Z));
            float v6 = (a6 * a6) * (a6 * a6) * grad(h6, x6, y6, z6);
            int vi6 = (h6 >> VEC3_HASH_SHIFT_3D) & VEC3_HASH_MASK_3D;
            vx += v6 * OUTPUT_VEC3[vi6 | 0];
            vy += v6 * OUTPUT_VEC3[vi6 | 1];
            vz += v6 * OUTPUT_VEC3[vi6 | 2];
        }
        else
        {
            float a7 = xAFlipMask0 + zAFlipMask0 + a0;
            if (a7 > 0) {
                float x7 = x0 - (xNMask | 1);
                float y7 = y0;
                float z7 = z0 - (zNMask | 1);
                int h7 = v3Hash(seed, xrbp + (~xNMask & PRIME_X), yrbp + (yNMask & PRIME_Y), zrbp + (~zNMask & PRIME_Z));
                float v7 = (a7 * a7) * (a7 * a7) * grad(h7, x7, y7, z7);
                int vi7 = (h7 >> VEC3_HASH_SHIFT_3D) & VEC3_HASH_MASK_3D;
                vx += v7 * OUTPUT_VEC3[vi7 | 0];
                vy += v7 * OUTPUT_VEC3[vi7 | 1];
                vz += v7 * OUTPUT_VEC3[vi7 | 2];
            }

            float a8 = yAFlipMask1 + a1;
            if (a8 > 0) {
                float x8 = x1;
                float y8 = (yNMask | 1) + y1;
                float z8 = z1;
                int h8 = v3Hash(seed2, xrbp + PRIME_X, yrbp + (yNMask & (PRIME_Y << 1)), zrbp + PRIME_Z);
                float v8 = (a8 * a8) * (a8 * a8) * grad(h8, x8, y8, z8);
                int vi8 = (h8 >> VEC3_HASH_SHIFT_3D) & VEC3_HASH_MASK_3D;
                vx += v8 * OUTPUT_VEC3[vi8 | 0];
                vy += v8 * OUTPUT_VEC3[vi8 | 1];
                vz += v8 * OUTPUT_VEC3[vi8 | 2];
                skip9 = true;
            }
        }

        boolean skipD = false;
        float aA = zAFlipMask0 + a0;
        if (aA > 0) {
            float xA = x0;
            float yA = y0;
            float zA = z0 - (zNMask | 1);
            int hA = v3Hash(seed, xrbp + (xNMask & PRIME_X), yrbp + (yNMask & PRIME_Y), zrbp + (~zNMask & PRIME_Z));
            float vA = (aA * aA) * (aA * aA) * grad(hA, xA, yA, zA);
            int viA = (hA >> VEC3_HASH_SHIFT_3D) & VEC3_HASH_MASK_3D;
            vx += vA * OUTPUT_VEC3[viA | 0];
            vy += vA * OUTPUT_VEC3[viA | 1];
            vz += vA * OUTPUT_VEC3[viA | 2];
        }
        else
        {
            float aB = xAFlipMask0 + yAFlipMask0 + a0;
            if (aB > 0) {
                float xB = x0 - (xNMask | 1);
                float yB = y0 - (yNMask | 1);
                float zB = z0;
                int hB = v3Hash(seed, xrbp + (~xNMask & PRIME_X), yrbp + (~yNMask & PRIME_Y), zrbp + (zNMask & PRIME_Z));
                float vB = (aB * aB) * (aB * aB) * grad(hB, xB, yB, zB);
                int viB = (hB >> VEC3_HASH_SHIFT_3D) & VEC3_HASH_MASK_3D;
                vx += vB * OUTPUT_VEC3[viB | 0];
                vy += vB * OUTPUT_VEC3[viB | 1];
                vz += vB * OUTPUT_VEC3[viB | 2];
            }

            float aC = zAFlipMask1 + a1;
            if (aC > 0) {
                float xC = x1;
                float yC = y1;
                float zC = (zNMask | 1) + z1;
                int hC = v3Hash(seed2, xrbp + PRIME_X, yrbp + PRIME_Y, zrbp + (zNMask & (PRIME_Z << 1)));
                float vC = (aC * aC) * (aC * aC) * grad(hC, xC, yC, zC);
                int viC = (hC >> VEC3_HASH_SHIFT_3D) & VEC3_HASH_MASK_3D;
                vx += vC * OUTPUT_VEC3[viC | 0];
                vy += vC * OUTPUT_VEC3[viC | 1];
                vz += vC * OUTPUT_VEC3[viC | 2];
                skipD = true;
            }
        }

        if (!skip5) {
            float a5 = yAFlipMask1 + zAFlipMask1 + a1;
            if (a5 > 0) {
                float x5 = x1;
                float y5 = (yNMask | 1) + y1;
                float z5 = (zNMask | 1) + z1;
                int h5 = v3Hash(seed2, xrbp + PRIME_X, yrbp + (yNMask & (PRIME_Y << 1)), zrbp + (zNMask & (PRIME_Z << 1)));
                float v5 = (a5 * a5) * (a5 * a5) * grad(h5, x5, y5, z5);
                int vi5 = (h5 >> VEC3_HASH_SHIFT_3D) & VEC3_HASH_MASK_3D;
                vx += v5 * OUTPUT_VEC3[vi5 | 0];
                vy += v5 * OUTPUT_VEC3[vi5 | 1];
                vz += v5 * OUTPUT_VEC3[vi5 | 2];
            }
        }

        if (!skip9) {
            float a9 = xAFlipMask1 + zAFlipMask1 + a1;
            if (a9 > 0) {
                float x9 = (xNMask | 1) + x1;
                float y9 = y1;
                float z9 = (zNMask | 1) + z1;
                int h9 = v3Hash(seed2, xrbp + (xNMask & (PRIME_X * 2)), yrbp + PRIME_Y, zrbp + (zNMask & (PRIME_Z << 1)));
                float v9 = (a9 * a9) * (a9 * a9) * grad(h9, x9, y9, z9);
                int vi9 = (h9 >> VEC3_HASH_SHIFT_3D) & VEC3_HASH_MASK_3D;
                vx += v9 * OUTPUT_VEC3[vi9 | 0];
                vy += v9 * OUTPUT_VEC3[vi9 | 1];
                vz += v9 * OUTPUT_VEC3[vi9 | 2];
            }
        }

        if (!skipD) {
            float aD = xAFlipMask1 + yAFlipMask1 + a1;
            if (aD > 0) {
                float xD = (xNMask | 1) + x1;
                float yD = (yNMask | 1) + y1;
                float zD = z1;
                int hD = v3Hash(seed2, xrbp + (xNMask & (PRIME_X << 1)), yrbp + (yNMask & (PRIME_Y << 1)), zrbp + PRIME_Z);
                float vD = (aD * aD) * (aD * aD) * grad(hD, xD, yD, zD);
                int viD = (hD >> VEC3_HASH_SHIFT_3D) & VEC3_HASH_MASK_3D;
                vx += vD * OUTPUT_VEC3[viD | 0];
                vy += vD * OUTPUT_VEC3[viD | 1];
                vz += vD * OUTPUT_VEC3[viD | 2];
            }
        }

        return new Vector3f(vx, vy, vz);
    }

    /*
     * Utility
     */

    private static float grad(long seed, long xrvp, long yrvp, long zrvp, float dx, float dy, float dz) {
        long hash = (seed ^ xrvp) ^ (yrvp ^ zrvp);
        hash *= HASH_MULTIPLIER;
        hash ^= hash >> (64 - N_GRADS_3D_EXPONENT - 2);
        return grad((int)hash, dx, dy, dz);
    }

    private static float grad(int hash, float dx, float dy, float dz) {
        int gi = (int)hash & GRAD_HASH_MASK_3D;
        return GRADIENTS_3D[gi | 0] * dx + GRADIENTS_3D[gi | 1] * dy + GRADIENTS_3D[gi | 2] * dz;
    }

    private static int v2Hash(long seed, long xrvp, long yrvp, long zrvp) {
        long hash = (seed ^ xrvp) ^ (yrvp ^ zrvp);
        hash *= HASH_MULTIPLIER;
        hash ^= hash >> (64 - N_GRADS_3D_EXPONENT + N_OUTPUT_VEC2 - 2);
        return (int)hash;
    }

    private static int v3Hash(long seed, long xrvp, long yrvp, long zrvp) {
        long hash = (seed ^ xrvp) ^ (yrvp ^ zrvp);
        hash *= HASH_MULTIPLIER;
        hash ^= hash >> (64 - N_GRADS_3D_EXPONENT + N_OUTPUT_VEC3 - 2);
        return (int)hash;
    }

    private static int fastFloor(double x) {
        int xi = (int)x;
        return x < xi ? xi - 1 : xi;
    }

    public record Vector2f(float x, float y) { }
    public record Vector3f(float x, float y, float z) { }

    /*
     * Lookup Tables & Gradients
     */

    private static float[] GRADIENTS_3D;
    static {

        GRADIENTS_3D = new float[N_GRADS_3D * 4];
        float[] grad3 = {
                2.22474487139f,       2.22474487139f,      -1.0f,                 0.0f,
                2.22474487139f,       2.22474487139f,       1.0f,                 0.0f,
                3.0862664687972017f,  1.1721513422464978f,  0.0f,                 0.0f,
                1.1721513422464978f,  3.0862664687972017f,  0.0f,                 0.0f,
                -2.22474487139f,       2.22474487139f,      -1.0f,                 0.0f,
                -2.22474487139f,       2.22474487139f,       1.0f,                 0.0f,
                -1.1721513422464978f,  3.0862664687972017f,  0.0f,                 0.0f,
                -3.0862664687972017f,  1.1721513422464978f,  0.0f,                 0.0f,
                -1.0f,                -2.22474487139f,      -2.22474487139f,       0.0f,
                1.0f,                -2.22474487139f,      -2.22474487139f,       0.0f,
                0.0f,                -3.0862664687972017f, -1.1721513422464978f,  0.0f,
                0.0f,                -1.1721513422464978f, -3.0862664687972017f,  0.0f,
                -1.0f,                -2.22474487139f,       2.22474487139f,       0.0f,
                1.0f,                -2.22474487139f,       2.22474487139f,       0.0f,
                0.0f,                -1.1721513422464978f,  3.0862664687972017f,  0.0f,
                0.0f,                -3.0862664687972017f,  1.1721513422464978f,  0.0f,
                //--------------------------------------------------------------------//
                -2.22474487139f,      -2.22474487139f,      -1.0f,                 0.0f,
                -2.22474487139f,      -2.22474487139f,       1.0f,                 0.0f,
                -3.0862664687972017f, -1.1721513422464978f,  0.0f,                 0.0f,
                -1.1721513422464978f, -3.0862664687972017f,  0.0f,                 0.0f,
                -2.22474487139f,      -1.0f,                -2.22474487139f,       0.0f,
                -2.22474487139f,       1.0f,                -2.22474487139f,       0.0f,
                -1.1721513422464978f,  0.0f,                -3.0862664687972017f,  0.0f,
                -3.0862664687972017f,  0.0f,                -1.1721513422464978f,  0.0f,
                -2.22474487139f,      -1.0f,                 2.22474487139f,       0.0f,
                -2.22474487139f,       1.0f,                 2.22474487139f,       0.0f,
                -3.0862664687972017f,  0.0f,                 1.1721513422464978f,  0.0f,
                -1.1721513422464978f,  0.0f,                 3.0862664687972017f,  0.0f,
                -1.0f,                 2.22474487139f,      -2.22474487139f,       0.0f,
                1.0f,                 2.22474487139f,      -2.22474487139f,       0.0f,
                0.0f,                 1.1721513422464978f, -3.0862664687972017f,  0.0f,
                0.0f,                 3.0862664687972017f, -1.1721513422464978f,  0.0f,
                -1.0f,                 2.22474487139f,       2.22474487139f,       0.0f,
                1.0f,                 2.22474487139f,       2.22474487139f,       0.0f,
                0.0f,                 3.0862664687972017f,  1.1721513422464978f,  0.0f,
                0.0f,                 1.1721513422464978f,  3.0862664687972017f,  0.0f,
                2.22474487139f,      -2.22474487139f,      -1.0f,                 0.0f,
                2.22474487139f,      -2.22474487139f,       1.0f,                 0.0f,
                1.1721513422464978f, -3.0862664687972017f,  0.0f,                 0.0f,
                3.0862664687972017f, -1.1721513422464978f,  0.0f,                 0.0f,
                2.22474487139f,      -1.0f,                -2.22474487139f,       0.0f,
                2.22474487139f,       1.0f,                -2.22474487139f,       0.0f,
                3.0862664687972017f,  0.0f,                -1.1721513422464978f,  0.0f,
                1.1721513422464978f,  0.0f,                -3.0862664687972017f,  0.0f,
                2.22474487139f,      -1.0f,                 2.22474487139f,       0.0f,
                2.22474487139f,       1.0f,                 2.22474487139f,       0.0f,
                1.1721513422464978f,  0.0f,                 3.0862664687972017f,  0.0f,
                3.0862664687972017f,  0.0f,                 1.1721513422464978f,  0.0f,
        };
        for (int i = 0; i < grad3.length; i++) {
            grad3[i] = (float)(grad3[i] / NORMALIZER_3D);
        }
        for (int i = 0, j = 0; i < GRADIENTS_3D.length; i++, j++) {
            if (j == grad3.length) j = 0;
            GRADIENTS_3D[i] = grad3[j];
        }
    }

    // 256 evenly-spaced angles, rotated a half-step
    private static final float[] OUTPUT_VEC2 = {
            0.012271538285719925f, 0.9999247018391445f, 0.03680722294135883f, 0.9993223845883495f, 0.06132073630220858f, 0.9981181129001492f, 0.0857973123444399f, 0.996312612182778f,
            0.11022220729388306f, 0.9939069700023561f, 0.13458070850712617f, 0.99090263542778f, 0.15885814333386145f, 0.9873014181578584f, 0.18303988795514095f, 0.9831054874312163f,
            0.20711137619221856f, 0.9783173707196277f, 0.2310581082806711f, 0.9729399522055602f, 0.25486565960451457f, 0.9669764710448521f, 0.27851968938505306f, 0.9604305194155658f,
            0.3020059493192281f, 0.9533060403541939f, 0.325310292162263f, 0.9456073253805213f, 0.34841868024943456f, 0.937339011912575f, 0.37131719395183754f, 0.9285060804732156f,
            0.3939920400610481f, 0.9191138516900578f, 0.4164295600976372f, 0.9091679830905224f, 0.43861623853852766f, 0.8986744656939538f, 0.46053871095824f, 0.8876396204028539f,
            0.4821837720791227f, 0.8760700941954066f, 0.5035383837257176f, 0.8639728561215868f, 0.5245896826784688f, 0.8513551931052652f, 0.5453249884220464f, 0.8382247055548381f,
            0.5657318107836131f, 0.8245893027850253f, 0.5857978574564389f, 0.8104571982525948f, 0.6055110414043255f, 0.7958369046088836f, 0.6248594881423863f, 0.7807372285720945f,
            0.6438315428897914f, 0.765167265622459f, 0.6624157775901718f, 0.7491363945234594f, 0.6806009977954529f, 0.7326542716724129f, 0.6983762494089728f, 0.7157308252838187f,
            0.7157308252838186f, 0.6983762494089729f, 0.7326542716724128f, 0.6806009977954531f, 0.7491363945234593f, 0.6624157775901718f, 0.7651672656224588f, 0.6438315428897915f,
            0.7807372285720944f, 0.6248594881423865f, 0.7958369046088835f, 0.6055110414043255f, 0.8104571982525948f, 0.585797857456439f, 0.8245893027850252f, 0.5657318107836132f,
            0.838224705554838f, 0.5453249884220465f, 0.8513551931052651f, 0.5245896826784691f, 0.8639728561215867f, 0.5035383837257176f, 0.8760700941954066f, 0.48218377207912283f,
            0.8876396204028538f, 0.46053871095824017f, 0.8986744656939538f, 0.4386162385385277f, 0.9091679830905223f, 0.4164295600976373f, 0.9191138516900578f, 0.3939920400610481f,
            0.9285060804732155f, 0.3713171939518376f, 0.9373390119125748f, 0.34841868024943473f, 0.9456073253805213f, 0.325310292162263f, 0.9533060403541938f, 0.3020059493192282f,
            0.9604305194155658f, 0.27851968938505306f, 0.9669764710448521f, 0.2548656596045146f, 0.9729399522055601f, 0.23105810828067128f, 0.9783173707196277f, 0.20711137619221856f,
            0.9831054874312163f, 0.18303988795514106f, 0.9873014181578583f, 0.1588581433338616f, 0.99090263542778f, 0.13458070850712622f, 0.9939069700023561f, 0.11022220729388318f,
            0.996312612182778f, 0.0857973123444401f, 0.9981181129001492f, 0.06132073630220865f, 0.9993223845883495f, 0.03680722294135899f, 0.9999247018391445f, 0.012271538285719944f,
            0.9999247018391445f, -0.012271538285719823f, 0.9993223845883495f, -0.036807222941358644f, 0.9981181129001492f, -0.06132073630220853f, 0.996312612182778f, -0.08579731234443976f,
            0.9939069700023561f, -0.11022220729388306f, 0.99090263542778f, -0.1345807085071261f, 0.9873014181578584f, -0.15885814333386128f, 0.9831054874312163f, -0.18303988795514092f,
            0.9783173707196277f, -0.20711137619221845f, 0.9729399522055602f, -0.23105810828067094f, 0.9669764710448521f, -0.2548656596045145f, 0.9604305194155659f, -0.27851968938505295f,
            0.9533060403541939f, -0.30200594931922786f, 0.9456073253805214f, -0.32531029216226287f, 0.937339011912575f, -0.3484186802494344f, 0.9285060804732156f, -0.3713171939518375f,
            0.9191138516900578f, -0.393992040061048f, 0.9091679830905225f, -0.416429560097637f, 0.8986744656939539f, -0.4386162385385274f, 0.8876396204028539f, -0.46053871095824006f,
            0.8760700941954066f, -0.4821837720791227f, 0.8639728561215868f, -0.5035383837257175f, 0.8513551931052652f, -0.5245896826784687f, 0.8382247055548382f, -0.5453249884220462f,
            0.8245893027850255f, -0.5657318107836129f, 0.8104571982525948f, -0.5857978574564389f, 0.7958369046088836f, -0.6055110414043254f, 0.7807372285720946f, -0.6248594881423862f,
            0.7651672656224591f, -0.6438315428897913f, 0.7491363945234596f, -0.6624157775901716f, 0.7326542716724128f, -0.680600997795453f, 0.7157308252838187f, -0.6983762494089728f,
            0.6983762494089729f, -0.7157308252838186f, 0.6806009977954532f, -0.7326542716724127f, 0.662415777590172f, -0.7491363945234591f, 0.6438315428897914f, -0.765167265622459f,
            0.6248594881423863f, -0.7807372285720945f, 0.6055110414043257f, -0.7958369046088835f, 0.585797857456439f, -0.8104571982525947f, 0.5657318107836135f, -0.8245893027850251f,
            0.5453249884220464f, -0.8382247055548381f, 0.524589682678469f, -0.8513551931052652f, 0.5035383837257177f, -0.8639728561215867f, 0.4821837720791229f, -0.8760700941954065f,
            0.4605387109582402f, -0.8876396204028538f, 0.43861623853852794f, -0.8986744656939537f, 0.41642956009763715f, -0.9091679830905224f, 0.39399204006104815f, -0.9191138516900578f,
            0.3713171939518377f, -0.9285060804732155f, 0.3484186802494348f, -0.9373390119125748f, 0.32531029216226326f, -0.9456073253805212f, 0.30200594931922803f, -0.9533060403541939f,
            0.27851968938505317f, -0.9604305194155658f, 0.2548656596045147f, -0.9669764710448521f, 0.23105810828067133f, -0.9729399522055601f, 0.20711137619221884f, -0.9783173707196275f,
            0.18303988795514134f, -0.9831054874312163f, 0.15885814333386147f, -0.9873014181578584f, 0.13458070850712628f, -0.99090263542778f, 0.11022220729388324f, -0.9939069700023561f,
            0.08579731234444016f, -0.996312612182778f, 0.06132073630220893f, -0.9981181129001492f, 0.03680722294135883f, -0.9993223845883495f, 0.012271538285720007f, -0.9999247018391445f,
            -0.012271538285719762f, -0.9999247018391445f, -0.03680722294135858f, -0.9993223845883495f, -0.061320736302208245f, -0.9981181129001492f, -0.08579731234443992f, -0.996312612182778f,
            -0.110222207293883f, -0.9939069700023561f, -0.13458070850712606f, -0.99090263542778f, -0.15885814333386122f, -0.9873014181578584f, -0.18303988795514065f, -0.9831054874312164f,
            -0.2071113761922186f, -0.9783173707196277f, -0.23105810828067108f, -0.9729399522055602f, -0.25486565960451446f, -0.9669764710448522f, -0.2785196893850529f, -0.9604305194155659f,
            -0.3020059493192278f, -0.953306040354194f, -0.3253102921622626f, -0.9456073253805215f, -0.34841868024943456f, -0.937339011912575f, -0.37131719395183743f, -0.9285060804732156f,
            -0.39399204006104793f, -0.9191138516900578f, -0.41642956009763693f, -0.9091679830905225f, -0.4386162385385273f, -0.898674465693954f, -0.46053871095824006f, -0.8876396204028539f,
            -0.48218377207912266f, -0.8760700941954066f, -0.5035383837257175f, -0.8639728561215868f, -0.5245896826784687f, -0.8513551931052653f, -0.5453249884220461f, -0.8382247055548382f,
            -0.5657318107836129f, -0.8245893027850255f, -0.5857978574564389f, -0.8104571982525948f, -0.6055110414043254f, -0.7958369046088836f, -0.6248594881423862f, -0.7807372285720946f,
            -0.6438315428897913f, -0.7651672656224591f, -0.6624157775901715f, -0.7491363945234596f, -0.680600997795453f, -0.7326542716724128f, -0.6983762494089728f, -0.7157308252838187f,
            -0.7157308252838185f, -0.698376249408973f, -0.7326542716724126f, -0.6806009977954532f, -0.749136394523459f, -0.662415777590172f, -0.765167265622459f, -0.6438315428897915f,
            -0.7807372285720944f, -0.6248594881423865f, -0.7958369046088835f, -0.6055110414043257f, -0.8104571982525952f, -0.5857978574564383f, -0.8245893027850256f, -0.5657318107836128f,
            -0.8382247055548383f, -0.545324988422046f, -0.8513551931052653f, -0.5245896826784686f, -0.8639728561215869f, -0.5035383837257174f, -0.8760700941954067f, -0.48218377207912255f,
            -0.887639620402854f, -0.4605387109582399f, -0.8986744656939538f, -0.4386162385385276f, -0.9091679830905224f, -0.4164295600976372f, -0.9191138516900577f, -0.3939920400610482f,
            -0.9285060804732155f, -0.37131719395183777f, -0.9373390119125752f, -0.348418680249434f, -0.9456073253805215f, -0.3253102921622624f, -0.953306040354194f, -0.30200594931922764f,
            -0.9604305194155659f, -0.2785196893850528f, -0.9669764710448522f, -0.25486565960451435f, -0.9729399522055602f, -0.23105810828067094f, -0.9783173707196277f, -0.20711137619221848f,
            -0.9831054874312163f, -0.18303988795514095f, -0.9873014181578583f, -0.15885814333386153f, -0.99090263542778f, -0.13458070850712636f, -0.9939069700023561f, -0.11022220729388242f,
            -0.9963126121827781f, -0.08579731234443934f, -0.9981181129001493f, -0.061320736302208106f, -0.9993223845883495f, -0.03680722294135844f, -0.9999247018391445f, -0.012271538285719624f,
            -0.9999247018391445f, 0.012271538285720144f, -0.9993223845883495f, 0.036807222941358964f, -0.9981181129001492f, 0.06132073630220863f, -0.996312612182778f, 0.08579731234443985f,
            -0.9939069700023561f, 0.11022220729388293f, -0.99090263542778f, 0.13458070850712597f, -0.9873014181578583f, 0.15885814333386203f, -0.9831054874312163f, 0.18303988795514148f,
            -0.9783173707196275f, 0.20711137619221898f, -0.9729399522055601f, 0.23105810828067147f, -0.9669764710448521f, 0.25486565960451485f, -0.9604305194155658f, 0.2785196893850533f,
            -0.9533060403541939f, 0.30200594931922814f, -0.9456073253805213f, 0.3253102921622629f, -0.937339011912575f, 0.3484186802494345f, -0.9285060804732156f, 0.3713171939518374f,
            -0.9191138516900574f, 0.3939920400610487f, -0.9091679830905222f, 0.4164295600976377f, -0.8986744656939536f, 0.43861623853852805f, -0.8876396204028537f, 0.4605387109582404f,
            -0.8760700941954065f, 0.482183772079123f, -0.8639728561215866f, 0.5035383837257178f, -0.8513551931052651f, 0.5245896826784691f, -0.838224705554838f, 0.5453249884220465f,
            -0.8245893027850253f, 0.5657318107836131f, -0.8104571982525949f, 0.5857978574564388f, -0.7958369046088837f, 0.6055110414043253f, -0.780737228572094f, 0.6248594881423869f,
            -0.7651672656224586f, 0.6438315428897918f, -0.749136394523459f, 0.6624157775901721f, -0.7326542716724126f, 0.6806009977954534f, -0.7157308252838185f, 0.698376249408973f,
            -0.6983762494089727f, 0.7157308252838188f, -0.680600997795453f, 0.7326542716724129f, -0.6624157775901718f, 0.7491363945234594f, -0.6438315428897915f, 0.7651672656224588f,
            -0.6248594881423865f, 0.7807372285720944f, -0.6055110414043257f, 0.7958369046088833f, -0.5857978574564384f, 0.8104571982525951f, -0.5657318107836128f, 0.8245893027850255f,
            -0.5453249884220461f, 0.8382247055548383f, -0.5245896826784686f, 0.8513551931052653f, -0.5035383837257174f, 0.8639728561215868f, -0.4821837720791226f, 0.8760700941954067f,
            -0.46053871095823995f, 0.8876396204028539f, -0.43861623853852766f, 0.8986744656939538f, -0.41642956009763726f, 0.9091679830905224f, -0.39399204006104827f, 0.9191138516900577f,
            -0.3713171939518378f, 0.9285060804732155f, -0.34841868024943407f, 0.9373390119125751f, -0.32531029216226254f, 0.9456073253805215f, -0.3020059493192277f, 0.953306040354194f,
            -0.27851968938505284f, 0.9604305194155659f, -0.2548656596045144f, 0.9669764710448522f, -0.231058108280671f, 0.9729399522055602f, -0.20711137619221853f, 0.9783173707196277f,
            -0.183039887955141f, 0.9831054874312163f, -0.15885814333386158f, 0.9873014181578583f, -0.13458070850712642f, 0.99090263542778f, -0.11022220729388249f, 0.9939069700023561f,
            -0.0857973123444394f, 0.996312612182778f, -0.06132073630220817f, 0.9981181129001492f, -0.036807222941358506f, 0.9993223845883495f, -0.012271538285719684f, 0.9999247018391445f,
    };

    // 256 vectors generated using fixed-radius Poisson rejection sampling, but with randomly-rotated tetrahedra on a sphere.
    // The algorithm was re-run with adjusted spacing until a result with exactly 256 vectors was returned.
    private static final float[] OUTPUT_VEC3 = {
            -0.7855810203217722f, -0.5120242863985095f, 0.34740983096092787f, 0f, 0.44195007914860773f, -0.5220252645635574f, -0.7294996577777727f, 0f, 0.6523541316318193f, 0.15368963425027363f, 0.7421681637384562f, 0f, -0.30872319045865493f, 0.8803599167117933f, -0.36007833692161156f, 0f,
            0.6898583244214584f, 0.545276495380974f, -0.47620272553972254f, 0f, -0.31541683712510465f, 0.49487059528212773f, 0.8097007550836961f, 0f, 0.39922526013793186f, -0.8737379241824829f, 0.27785109593645013f, 0f, -0.7736667474342858f, -0.1664091664806189f, -0.6113491254804239f, 0f,
            0.7639719174602074f, 0.47770928933593376f, 0.4337519385712658f, 0f, -0.8625240947058035f, 0.39877897332167134f, 0.31149240197547f, 0f, 0.028345693290678365f, -0.9683321718294965f, 0.24805105658279553f, 0f, 0.0702064839549177f, 0.09184390917189147f, -0.9932953971295313f, 0f,
            0.17898593643336067f, 0.2770098767964981f, -0.9440495552227447f, 0f, 0.6988643568353057f, -0.6519001791695181f, 0.2943038687206214f, 0f, -0.9012969720852568f, -0.42954498183127393f, 0.05616828903855023f, 0f, 0.023446678816590575f, 0.8044352842042939f, 0.593577397463573f, 0f,
            -0.09078368654499426f, 0.9688670514753966f, 0.2303366206721613f, 0f, -0.8265312435094146f, -0.30975471280930633f, -0.4699980014800116f, 0f, 0.12608836044662863f, -0.5313821467246918f, 0.8376960901795706f, 0f, 0.7912265696077803f, -0.12773019194139865f, -0.5980347093717203f, 0f,
            0.9658556025219398f, 0.2508753898995745f, -0.06468766358214734f, 0f, -0.19584262894626087f, -0.3372687578399138f, 0.9208123857076091f, 0f, -0.2038408673117624f, -0.7160505583854263f, -0.6676230213597398f, 0f, -0.5661721062639166f, 0.8024439263257656f, -0.18850170076572195f, 0f,
            0.47478812322636377f, -0.04041326833958074f, -0.8791717726276818f, 0f, 0.20628690920529813f, -0.8235499246852277f, 0.5284044214817669f, 0f, -0.9860697317427367f, 0.057649242062784954f, -0.15602259141032665f, 0f, 0.30499469931107476f, 0.8063139509620236f, 0.5067899425562415f, 0f,
            0.3091079818184893f, 0.717919025397776f, -0.6237343413249009f, 0f, -0.9918703518605176f, 0.06017270755592631f, 0.1121269386253013f, 0f, 0.4435980305916518f, 0.11671666630171577f, 0.8885932742614148f, 0f, 0.23916433945037646f, -0.894808399255418f, -0.37698587156181523f, 0f,
            0.7643205245554541f, -0.35453878820241386f, 0.5386245291510954f, 0f, -0.557928635916593f, -0.7620082322224426f, -0.32872342668188015f, 0f, 0.35318344853967965f, 0.5156096460607112f, -0.7806459790242742f, 0f, -0.5595753371785408f, 0.6009373743641454f, 0.5707448765550588f, 0f,
            -0.37360693198037276f, -0.5527496306873113f, 0.7449065083292359f, 0f, 0.5390642001480362f, 0.7514311590063892f, 0.3804747053262182f, 0f, -0.7496142794659885f, 0.37151761026693697f, -0.5477710263259928f, 0f, 0.584157011298325f, -0.5701991385860148f, -0.5776101873294613f, 0f,
            0.07181287991979796f, 0.9335803476879151f, -0.351099194941453f, 0f, -0.95649664821889f, -0.20564815413689758f, 0.206937185266727f, 0f, 0.5471219409615322f, -0.08600273116743627f, 0.8326230311192621f, 0f, 0.3375618273375599f, -0.6419294623835812f, -0.688461021444536f, 0f,
            -0.3202365165496523f, -0.6268257800731106f, -0.7103083942231575f, 0f, 0.573072425743534f, 0.7075490759450512f, -0.41347587593682633f, 0f, 0.5332815615333012f, -0.5075548385167832f, 0.6767561318723356f, 0f, -0.7861174707271829f, 0.4268315426448425f, 0.4470281382876482f, 0f,
            0.9908212464937525f, 0.09780232581794611f, -0.09331646457708814f, 0f, -0.2147102203879096f, -0.9198781820278473f, 0.3282128143295414f, 0f, -0.43459554780442916f, 0.14675795161651578f, -0.8885881011277905f, 0f, -0.34151547830141377f, 0.6753179045933854f, 0.6536917513753374f, 0f,
            -0.7339345902809137f, 0.631095810794388f, -0.25111370887898843f, 0f, 0.49174971818087737f, 0.3592955942314991f, 0.7931512407068791f, 0f, 0.6327216030521022f, -0.09800524484662071f, -0.7681525532169583f, 0f, -0.390536730952066f, -0.8923861601792663f, 0.22611502138906758f, 0f,
            0.1121912831058285f, 0.5066428688698842f, -0.8548251981653098f, 0f, 0.046733537582263125f, -0.9815045200791781f, -0.18564712098276553f, 0f, -0.8875261072386896f, 0.2208715392951393f, 0.40436761999340887f, 0f, 0.728601286550598f, 0.2539901119141545f, 0.6361046991546666f, 0f,
            0.5950138811230042f, 0.7809811791796447f, 0.18980747887823315f, 0f, -0.49524927285502407f, -0.2515844962554219f, 0.8315247434565465f, 0f, 0.5538740555879639f, -0.7745887945982852f, -0.30534526331567124f, 0f, -0.6536386638559438f, 0.24519211167406232f, -0.7159869590191084f, 0f,
            -0.41830816885294014f, 0.12591490972994984f, 0.8995352752274913f, 0f, -0.08183753692570106f, 0.8458999380022647f, -0.5270255330035721f, 0f, -0.4663702234432518f, -0.740583703280193f, -0.4837712198149931f, 0f, 0.966515929221893f, -0.23123114445202145f, 0.11126147759107384f, 0f,
            -0.4265547642405082f, -0.8264639685736173f, 0.3674348129305f, 0f, 0.12589025143907626f, -0.10043331534459055f, -0.9869472092070073f, 0f, 0.8886873323636982f, 0.13892244603318246f, 0.43696153066883486f, 0f, -0.5880228195622664f, 0.7879748378850253f, 0.18255086560767242f, 0f,
            0.23617714375311633f, -0.7632615689695902f, -0.6013752024341371f, 0f, -0.9491498510349248f, 0.2702564216845885f, -0.16148073203519378f, 0f, 0.6039769693397199f, 0.7738542303104692f, -0.1906867870037102f, 0f, 0.10899573794208851f, -0.2808490830254677f, 0.9535427214730411f, 0f,
            0.7924276754540488f, 0.5408601639552588f, -0.2820082662277532f, 0f, 0.23132144400111276f, -0.9723239921804654f, -0.03280920260357504f, 0f, -0.7647094586838532f, 0.18175285382180814f, -0.6182114071473258f, 0f, -0.25903966077130836f, 0.24971097440339848f, 0.9330288759786541f, 0f,
            0.4397601676837157f, -0.6997968706307749f, 0.562934574150642f, 0f, 0.2834257067216901f, -0.11174766674169158f, -0.9524611948773023f, 0f, 0.2701171594997633f, 0.9066853669862031f, 0.3239727850250293f, 0f, -0.9933030339051692f, -0.09514082961373649f, 0.06555383570163097f, 0f,
            0.9349322364581333f, -0.2645457623526797f, -0.23646829143171166f, 0f, 0.008172544375040336f, 0.5321457550704652f, 0.8466133148485928f, 0f, -0.5565375134322379f, -0.8209553246012699f, 0.12766499579558016f, 0f, -0.3865672674009357f, 0.5533553318834843f, -0.7378100192124614f, 0f,
            -0.14771665249184782f, 0.5637174045046185f, 0.8126515110643546f, 0f, -0.8193273523953906f, 0.01999101980440576f, -0.5729773544774168f, 0f, 0.1897289121296018f, -0.9417767529293446f, 0.2775959069294394f, 0f, 0.7773150927576365f, 0.3580683286203203f, -0.5172700635163773f, 0f,
            0.9544041768521164f, 0.15931966561958918f, 0.2524478388778491f, 0f, -0.03845851220995278f, -0.6195723266631504f, -0.7839968589680747f, 0f, -0.43068993391588156f, 0.8697216922601873f, -0.24101941588928022f, 0f, -0.48525573072628214f, -0.4094690312166259f, 0.772568435979506f, 0f,
            0.062257144698645545f, 0.7156351982942329f, -0.6956941216485487f, 0f, -0.27099815783342807f, -0.8608682546638503f, -0.4306574585014225f, 0f, -0.6811964779845073f, 0.2591093217137989f, 0.6847143329758167f, 0f, 0.8899374911192899f, -0.11387626534418144f, 0.44163724717415465f, 0f,
            -0.09083418400031985f, 0.42071554527569943f, 0.9026336914719985f, 0f, 0.9170592331424601f, -0.38684959714203954f, -0.09669411615607605f, 0f, -0.14594129544681889f, 0.6923292134443365f, -0.7066692284901407f, 0f, -0.6802837536953215f, -0.7261951615779964f, -0.09927034682578173f, 0f,
            0.7374614901681634f, 0.6528283855425644f, 0.17310589114423774f, 0f, -0.010589427919937458f, -0.6904336463720835f, 0.7233182176426067f, 0f, 0.14900921006274737f, -0.4447395653449917f, -0.8831777705157801f, 0f, -0.8758812723109735f, 0.48234482617451074f, -0.013246338271064364f, 0f,
            -0.14426168944564297f, -0.32212999945076104f, -0.935639261901798f, 0f, 0.7311833787802975f, -0.5322101589550988f, 0.42675896393680635f, 0f, -0.8437581416720136f, -0.11198175022904902f, 0.5249116935045639f, 0f, 0.256836452337359f, 0.966321908634909f, -0.016031395539572165f, 0f,
            -0.5876935057953628f, -0.34346629798712786f, -0.7325621102629956f, 0f, 0.2531633987218102f, 0.9480619259819005f, -0.1925795369480124f, 0f, 0.826014920275886f, -0.5609238792814593f, 0.055351180054814854f, 0f, -0.49148481320233345f, -0.04367174871331346f, 0.8697904671561931f, 0f,
            -0.2992582955840803f, 0.88667955334857f, 0.35248239984103f, 0f, 0.9918158346799334f, -0.08056651463281023f, 0.09904739672534676f, 0f, -0.24562792887498533f, -0.07462644755182765f, -0.9664873583665643f, 0f, -0.44692961022086775f, -0.731486591163932f, 0.5149575618001876f, 0f,
            0.16481613315703178f, 0.3856526101410991f, 0.9078037819609043f, 0f, 0.8471118688582433f, 0.02334587018665747f, -0.5309015464139489f, 0f, -0.7016312141434594f, 0.5372623419886053f, -0.4680414674155422f, 0f, -0.31029678787181564f, -0.9462608223163619f, 0.09113923186858677f, 0f,
            -0.7201841829396932f, -0.658399461099545f, 0.21873475322251684f, 0f, -0.12488518359422617f, 0.8256613276636618f, 0.5501700309170936f, 0f, -0.04757203321794458f, 0.23585192146494155f, -0.9706239090383081f, 0f, 0.892641399751864f, -0.4031137880290584f, 0.20171912489869762f, 0f,
            0.974923160633712f, 0.19705317905806916f, 0.10341603107391945f, 0f, -0.31364375410952955f, 0.3278507216103314f, -0.8911461719873267f, 0f, -0.5120791324869708f, 0.4618598712687528f, 0.7241963969690837f, 0f, -0.14920027403721162f, -0.9867637719371533f, 0.06353374394432351f, 0f,
            -0.7755427508331851f, -0.26363766379442694f, 0.5736101671510999f, 0f, 0.36062604936339443f, -0.8084659632807343f, -0.46511465117442263f, 0f, 0.7152662341439137f, 0.40280728609827815f, 0.5710871251917031f, 0f, -0.300349532674123f, 0.6692963409768831f, -0.6795826411683804f, 0f,
            -0.8811356449101773f, -0.42403254136346724f, 0.20927584460124762f, 0f, -0.059035622266311136f, 0.4756159418990514f, -0.8776697961734286f, 0f, 0.706186351602365f, -0.7062286148323401f, -0.050418056314529974f, 0f, 0.2339849155741233f, 0.6546452142967559f, 0.718812007886711f, 0f,
            0.905669019138563f, 0.14492719158467698f, -0.39844665503925564f, 0f, -0.03739078737640134f, -0.9057588653322484f, 0.4221407429891251f, 0f, -0.6937015806648842f, 0.062219433330326486f, -0.7175701074453245f, 0f, -0.17457665109727752f, 0.6986122404172451f, 0.6938760194954549f, 0f,
            -0.9713121712088488f, -0.19366594171334378f, -0.13800785876839003f, 0f, 0.5137987576728388f, -0.27715753015451955f, -0.8119079628212142f, 0f, 0.331804206073051f, -0.5089581885664506f, 0.7942716985537063f, 0f, 0.1257092074629589f, 0.9797816604343138f, 0.15564412303589784f, 0f,
            0.3292249339658975f, -0.060348617990981096f, 0.9423210637366275f, 0f, -0.6831406673539048f, 0.727067116992706f, -0.06849989777493613f, 0f, -0.41280264895025337f, -0.8713220765494749f, -0.26531455281828215f, 0f, 0.7667183823382606f, 0.20460357754775013f, -0.6085066131434091f, 0f,
            -0.27831871592211327f, -0.953631748696043f, -0.11456517903094776f, 0f, 0.30026987719073683f, 0.3676449293831255f, -0.8801563535819853f, 0f, 0.7523949387768516f, 0.051056691786989586f, 0.6567305918919463f, 0f, -0.7743461000454751f, 0.534930127525928f, 0.33799094072098673f, 0f,
            0.46682087323242466f, 0.8782241497890452f, -0.1039260075333665f, 0f, -0.05244580833327629f, -0.23687514549758343f, 0.9701234986504389f, 0f, 0.5093345173266459f, -0.7081545460533599f, -0.4889739137864027f, 0f, -0.9237095822257945f, 0.06680554176189814f, -0.37722357733066963f, 0f,
            0.8635864148735739f, -0.44175482333032456f, -0.24304563380196953f, 0f, -0.45768667986861056f, 0.31688564944551223f, -0.8307264220225236f, 0f, 0.18156114860202402f, 0.780060303900726f, 0.5987833260843364f, 0f, -0.5874608836069876f, -0.6551911300159137f, 0.4749887297401567f, 0f,
            -0.6015712437089378f, 0.07741823611655728f, 0.7950587748456607f, 0f, 0.6255294979116325f, 0.7798931309536391f, -0.021898665113450344f, 0f, -0.5504358200423933f, -0.009361595640110087f, -0.8348249927627521f, 0f, 0.5264775658396986f, -0.8479497714300862f, 0.06166488303054183f, 0f,
            0.8339384681196594f, 0.4130400291907359f, -0.36598437900592384f, 0f, 0.235653272091136f, -0.6188640421353756f, 0.7493162434543944f, 0f, -0.46291313673801915f, -0.512943821956297f, -0.7229108266945755f, 0f, -0.6066786034727764f, 0.7187678349009367f, 0.3395789622461052f, 0f,
            -0.18492443560944202f, 0.9471244707213824f, -0.26221782943817995f, 0f, 0.8104995806690649f, -0.3280418298260722f, -0.4852617722627987f, 0f, -0.7853098392670248f, -0.5713097964426467f, -0.2385237364269253f, 0f, 0.15973469420740194f, -0.04777284445266355f, 0.9860033381279039f, 0f,
            -0.5846642868962936f, -0.4761197915443017f, 0.6568695576199732f, 0f, 0.9385378171725712f, 0.022778786445822787f, 0.34442400123249367f, 0f, -0.33191181632901495f, 0.934966122636414f, -0.12519143622230866f, 0f, -0.021961713947262804f, -0.4816251175379351f, -0.8761021226301583f, 0f,
            0.07881018032278211f, -0.7431175534688199f, -0.6645037676371792f, 0f, 0.47965074791085205f, 0.8071591161945931f, -0.34413561450760016f, 0f, 0.4067426201429073f, -0.2845610111537385f, 0.8680930087786932f, 0f, -0.9652035483765415f, 0.22051944842796525f, 0.14054637336608597f, 0f,
            0.45388387180939127f, 0.8493354231594679f, 0.269478700230359f, 0f, -0.9335360656174874f, -0.0072998300541518785f, 0.35840916097751047f, 0f, 0.5051561413049904f, -0.7797638172643385f, 0.36984545716717937f, 0f, -0.02550394749689448f, -0.062271775840977574f, -0.997733318375049f, 0f,
            -0.4092342171937601f, -0.0689431817225483f, -0.9098209676479114f, 0f, -0.7070425172350995f, -0.1332827569313765f, 0.6944973617852163f, 0f, 0.41165510243749687f, 0.9043465009352507f, 0.11268310824317235f, 0f, 0.7046216319913629f, -0.7021205622813258f, 0.10264049761952293f, 0f,
            -0.04732692985811512f, -0.43660908160697626f, 0.8984056275249603f, 0f, -0.2253853814632334f, 0.9702313712833114f, 0.08861442320738756f, 0f, -0.6520310976904612f, -0.4385384800961924f, -0.6184977357432229f, 0f, 0.9247434090118098f, -0.09508380958014283f, -0.36852231498912497f, 0f,
            0.7148113862648666f, -0.6408204529844312f, -0.27998898032409814f, 0f, 0.08817855466403214f, 0.19141711219558266f, 0.9775397852036847f, 0f, -0.8975827324986974f, -0.4018205910507603f, -0.1813434612214888f, 0f, 0.09459279156979862f, 0.851223931839609f, -0.5162073436580978f, 0f,
            0.38834439758730316f, -0.9200927943061323f, 0.051165210138218044f, 0f, 0.6681588581089581f, 0.6630612634046952f, 0.33751074249970575f, 0f, -0.8265686889408728f, 0.04461183159577979f, 0.5610650469813965f, 0f, -0.22993456675538845f, 0.21241969930565738f, -0.9497409996193203f, 0f,
            0.015661445499380405f, -0.29950139650432955f, -0.9539673121322483f, 0f, -0.7435326225126199f, 0.655670227637051f, 0.13136130271082064f, 0f, 0.8715521103783915f, 0.43441076044747673f, 0.22734161542141176f, 0f, -0.1436809333651521f, -0.7905795915801981f, 0.5952643940000159f, 0f,
            -0.817954133682781f, -0.5576167017996086f, 0.14147313902426017f, 0f, -0.23282141392641648f, 0.7984287565582466f, -0.5552528342279365f, 0f, 0.6957020891383141f, -0.5422550792034123f, -0.4711242216726776f, 0f, 0.3550734584708834f, 0.3014430244447743f, 0.884903916876354f, 0f,
            -0.3575134676849902f, -0.20022364322700392f, 0.9121922018504427f, 0f, 0.03779196186883926f, -0.8435244043543135f, -0.5357595980258353f, 0f, 0.9191295790607819f, 0.38584937185206164f, 0.07950521452659119f, 0f, -0.5994080732446307f, 0.6578986757292559f, -0.4559378183511987f, 0f,
            -0.6551222390919279f, -0.27776599810299163f, 0.7026100640789509f, 0f, -0.04915221319891831f, 0.9904665365247827f, -0.12868605193357033f, 0f, -0.21958429393379097f, -0.45926867116707687f, -0.8607293567330127f, 0f, 0.9238587462246373f, -0.253431867254714f, 0.28680534458763207f, 0f,
            -0.015744452317840597f, 0.9962813606878287f, 0.08470869239469932f, 0f, 0.7017725764385561f, -0.26734227648519915f, -0.6603357919744335f, 0f, -0.8931428665102982f, -0.3219913321345033f, -0.31403407781961473f, 0f, 0.20711474238958252f, -0.4069477520681264f, 0.8896611773993487f, 0f,
            -0.859295202830418f, 0.17308020503557464f, -0.4813055131800116f, 0f, -0.14530676286026736f, 0.09207305092013009f, 0.9850931417695127f, 0f, 0.6883365676916832f, 0.6610687603034884f, -0.2986316522561567f, 0f, 0.3162653979990022f, -0.9262220162591931f, -0.20515597633334434f, 0f,
            0.9231930920545681f, 0.2894435914641443f, -0.25285751351910263f, 0f, -0.43748359371402923f, 0.7173883541375309f, 0.5421826763913641f, 0f, 0.05014618509275305f, -0.8410951154047953f, 0.53855767282886f, 0f, -0.535855683433292f, -0.16573683019687993f, -0.8278828357011215f, 0f,
            0.5581665568888212f, -0.8217993717430312f, -0.11443726392061947f, 0f, -0.06281050728721442f, 0.22763938298436667f, 0.9717176295039727f, 0f, -0.9166882864336425f, -0.16641023805806285f, -0.36330458046219427f, 0f, 0.4213322368320357f, 0.7605702268167275f, -0.49397578512115875f, 0f,
            -0.6892017032413691f, 0.4268724214048745f, 0.5854749764875824f, 0f, -0.4126567283550981f, -0.24312428688850635f, -0.8778411050231256f, 0f, 0.3496796586736489f, -0.8250596378448894f, 0.44384764312671576f, 0f, 0.7521787729228183f, 0.6413115033285213f, -0.15148151459117262f, 0f,
            -0.4640280001194681f, 0.835481757740562f, 0.2943607439654067f, 0f, 0.5824682057697544f, 0.20194193672788777f, -0.7873691913314572f, 0f, -0.680399664847503f, -0.6866390026981048f, -0.25609212414511556f, 0f, 0.5619594591972169f, -0.3507846917703449f, 0.749100571511166f, 0f,
            -0.9562752416031773f, -0.28348095771831183f, 0.07194587484973432f, 0f, 0.11484365097654696f, 0.8968595424816005f, 0.4271462242488764f, 0f, 0.25997007902783403f, 0.054211107230928125f, -0.9640937267003964f, 0f, 0.5814615115987963f, -0.6675896919942167f, 0.46500162760178576f, 0f,
            0.8125786605926297f, -0.5030334858429766f, 0.29440318012914163f, 0f, 0.2784776203859986f, 0.8164474350964072f, -0.5058298139380969f, 0f, -0.5333305076141349f, 0.2702532573534887f, 0.8015745420969146f, 0f, -0.5577257733644934f, -0.5836672066069193f, -0.5901479082879593f, 0f,
            -0.802702296285481f, 0.4311004041290742f, -0.4120940003176034f, 0f, 0.5168288385619022f, 0.6827770537768403f, 0.5164333901545998f, 0f, 0.5794161049992462f, -0.38238247138332826f, -0.7197642828359022f, 0f, -0.29354264727566753f, -0.7314949865225862f, 0.6154248929989058f, 0f,
    };

}
