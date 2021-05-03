/**
 * K.jpg's OpenSimplex 2, smooth variant ("SuperSimplex")
 * 
 * WIP to replace implementation in OpenSimplex2 repo.
 * 2D and 3D functions replaced with statelessly-seedable,
 * highly optimized versions. Based on work I did for FastNoiseLite.
 * 4D noise will come before I move this to the main repo.
 */

using System.Runtime.CompilerServices;

namespace Noise
{
    public class OpenSimplex2S_Stateless
    {
        private const int PRIME_X = 7691;
        private const int PRIME_Y = 10163;
        private const int PRIME_Z = 4231;

        /*
         * Noise Evaluators
         */

        /**
         * 2D OpenSimplex2S/SuperSimplex noise, standard lattice orientation.
         */
        public static float Noise2(long seed, double x, double y)
        {

            // Get points for A2* lattice
            double s = 0.366025403784439 * (x + y);
            double xs = x + s, ys = y + s;

            return noise2_Base(seed, xs, ys);
        }

        /**
         * 2D OpenSimplex2S/SuperSimplex noise, with Y pointing down the main diagonal.
         * Might be better for a 2D sandbox style game, where Y is vertical.
         * Probably slightly less optimal for heightmaps or continent maps,
         * unless your map is centered around an equator. It's a slight
         * difference, but the option is here to make it easy.
         */
        public static float Noise2_ImproveXAxis(long seed, double x, double y)
        {

            // Skew transform and rotation baked into one.
            double xx = x * 0.7071067811865476;
            double yy = y * 1.224744871380249;

            return noise2_Base(seed, yy + xx, yy - xx);
        }

        /**
         * 2D  OpenSimplex2S/SuperSimplex noise base.
         */
        private static float noise2_Base(long seed, double xs, double ys)
        {
            // 2D OpenSimplex2S case is a modified 2D simplex noise.

            const double G2 = -0.21132486540518713;

            int xsb = FastFloor(xs);
            int ysb = FastFloor(ys);
            float xi = (float)(xs - xsb);
            float yi = (float)(ys - ysb);

            int xsbp = xsb * PRIME_X;
            int ysbp = ysb * PRIME_Y;

            float t = (xi + yi) * (float)G2;
            float x0 = xi + t;
            float y0 = yi + t;

            float a0 = (2.0f / 3.0f) - x0 * x0 - y0 * y0;
            float value = (a0 * a0) * (a0 * a0) * Grad(seed, xsbp, ysbp, x0, y0);

            float a1 = (float)(2 * (1 + 2 * G2) * (1 / G2 + 2)) * t + ((float)(-2 * (1 + 2 * G2) * (1 + 2 * G2)) + a0);
            float x1 = x0 - (float)(1 + 2 * G2);
            float y1 = y0 - (float)(1 + 2 * G2);
            value += (a1 * a1) * (a1 * a1) * Grad(seed, xsbp + PRIME_X, ysbp + PRIME_Y, x1, y1);

            // Nested conditionals were faster than compact bit logic/arithmetic.
            float xmyi = xi - yi;
            if (t < G2)
            {
                if (xi + xmyi > 1)
                {
                    float x2 = x0 - (float)(3 * G2 + 2);
                    float y2 = y0 - (float)(3 * G2 + 1);
                    float a2 = (2.0f / 3.0f) - x2 * x2 - y2 * y2;
                    if (a2 > 0)
                    {
                        value += (a2 * a2) * (a2 * a2) * Grad(seed, xsbp + (PRIME_X << 1), ysbp + PRIME_Y, x2, y2);
                    }
                }
                else
                {
                    float x2 = x0 - (float)G2;
                    float y2 = y0 - (float)(G2 + 1);
                    float a2 = (2.0f / 3.0f) - x2 * x2 - y2 * y2;
                    if (a2 > 0)
                    {
                        value += (a2 * a2) * (a2 * a2) * Grad(seed, xsbp, ysbp + PRIME_Y, x2, y2);
                    }
                }

                if (yi - xmyi > 1)
                {
                    float x3 = x0 - (float)(3 * G2 + 1);
                    float y3 = y0 - (float)(3 * G2 + 2);
                    float a3 = (2.0f / 3.0f) - x3 * x3 - y3 * y3;
                    if (a3 > 0)
                    {
                        value += (a3 * a3) * (a3 * a3) * Grad(seed, xsbp + PRIME_X, ysbp + (PRIME_Y << 1), x3, y3);
                    }
                }
                else
                {
                    float x3 = x0 - (float)(G2 + 1);
                    float y3 = y0 - (float)G2;
                    float a3 = (2.0f / 3.0f) - x3 * x3 - y3 * y3;
                    if (a3 > 0)
                    {
                        value += (a3 * a3) * (a3 * a3) * Grad(seed, xsbp + PRIME_X, ysbp, x3, y3);
                    }
                }
            }
            else
            {
                if (xi + xmyi < 0)
                {
                    float x2 = x0 + (float)(1 + G2);
                    float y2 = y0 + (float)G2;
                    float a2 = (2.0f / 3.0f) - x2 * x2 - y2 * y2;
                    if (a2 > 0)
                    {
                        value += (a2 * a2) * (a2 * a2) * Grad(seed, xsbp - PRIME_X, ysbp, x2, y2);
                    }
                }
                else
                {
                    float x2 = x0 - (float)(G2 + 1);
                    float y2 = y0 - (float)G2;
                    float a2 = (2.0f / 3.0f) - x2 * x2 - y2 * y2;
                    if (a2 > 0)
                    {
                        value += (a2 * a2) * (a2 * a2) * Grad(seed, xsbp + PRIME_X, ysbp, x2, y2);
                    }
                }

                if (yi < xmyi)
                {
                    float x2 = x0 + (float)G2;
                    float y2 = y0 + (float)(G2 + 1);
                    float a2 = (2.0f / 3.0f) - x2 * x2 - y2 * y2;
                    if (a2 > 0)
                    {
                        value += (a2 * a2) * (a2 * a2) * Grad(seed, xsbp, ysbp - PRIME_Y, x2, y2);
                    }
                }
                else
                {
                    float x2 = x0 - (float)G2;
                    float y2 = y0 - (float)(G2 + 1);
                    float a2 = (2.0f / 3.0f) - x2 * x2 - y2 * y2;
                    if (a2 > 0)
                    {
                        value += (a2 * a2) * (a2 * a2) * Grad(seed, xsbp, ysbp + PRIME_Y, x2, y2);
                    }
                }
            }

            return value;
        }

        /**
         * 3D Re-oriented 8-point BCC noise, classic orientation
         * Proper substitute for what 3D SuperSimplex "should" be,
         * in light of Forbidden Formulae.
         * Use noise3_XYBeforeZ or noise3_XZBeforeY instead, wherever appropriate.
         */
        public static float Noise3_Classic(long seed, double x, double y, double z)
        {

            // Re-orient the cubic lattices via rotation, to produce the expected look on cardinal planar slices.
            // If texturing objects that don't tend to have cardinal plane faces, you could even remove this.
            // Orthonormal rotation. Not a skew transform.
            double r = (2.0 / 3.0) * (x + y + z);
            double xr = r - x, yr = r - y, zr = r - z;

            // Evaluate both lattices to form a BCC lattice.
            return noise3_BCC(seed, xr, yr, zr);
        }

        /**
         * 3D Re-oriented 8-point BCC noise, with better visual isotropy in (X, Y).
         * Recommended for 3D terrain and time-varied animations.
         * The Z coordinate should always be the "different" coordinate in your use case.
         * If Y is vertical in world coordinates, call noise3_XYBeforeZ(x, z, Y) or use noise3_XZBeforeY.
         * If Z is vertical in world coordinates, call noise3_XYBeforeZ(x, y, Z).
         * For a time varied animation, call noise3_XYBeforeZ(x, y, T).
         */
        public static float Noise3_ImproveXYPlanes(long seed, double x, double y, double z)
        {

            // Re-orient the cubic lattices without skewing, to make X and Y triangular like 2D.
            // Orthonormal rotation. Not a skew transform.
            double xy = x + y;
            double s2 = xy * -0.211324865405187;
            double zz = z * 0.577350269189626;
            double xr = x + s2 - zz, yr = y + s2 - zz;
            double zr = xy * 0.577350269189626 + zz;

            // Evaluate both lattices to form a BCC lattice.
            return noise3_BCC(seed, xr, yr, zr);
        }

        /**
         * 3D Re-oriented 8-point BCC noise, with better visual isotropy in (X, Z).
         * Recommended for 3D terrain and time-varied animations.
         * The Y coordinate should always be the "different" coordinate in your use case.
         * If Y is vertical in world coordinates, call noise3_XZBeforeY(x, Y, z).
         * If Z is vertical in world coordinates, call noise3_XZBeforeY(x, Z, y) or use noise3_XYBeforeZ.
         * For a time varied animation, call noise3_XZBeforeY(x, T, y) or use noise3_XYBeforeZ.
         */
        public static float Noise3_ImproveXZPlanes(long seed, double x, double y, double z)
        {

            // Re-orient the cubic lattices without skewing, to make X and Z triangular like 2D.
            // Orthonormal rotation. Not a skew transform.
            double xz = x + z;
            double s2 = xz * -0.211324865405187;
            double yy = y * 0.577350269189626;
            double xr = x + s2 - yy; double zr = z + s2 - yy;
            double yr = xz * 0.577350269189626 + yy;

            // Evaluate both lattices to form a BCC lattice.
            return noise3_BCC(seed, xr, yr, zr);
        }

        /**
         * Generate overlapping cubic lattices for 3D Re-oriented BCC noise.
         * Lookup table implementation inspired by DigitalShadow.
         * It was actually faster to narrow down the points in the loop itself,
         * than to build up the index with enough info to isolate 8 points.
         */
        private static float noise3_BCC(long seed, double xr, double yr, double zr)
        {
            int xrbp = FastFloor(xr);
            int yrbp = FastFloor(yr);
            int zrbp = FastFloor(zr);
            float xi = (float)(xr - xrbp);
            float yi = (float)(yr - yrbp);
            float zi = (float)(zr - zrbp);

            xrbp *= PRIME_X;
            yrbp *= PRIME_Y;
            zrbp *= PRIME_Z;
            long seed2 = seed + 1293373;

            int xNMask = (int)(-0.5f - xi);
            int yNMask = (int)(-0.5f - yi);
            int zNMask = (int)(-0.5f - zi);

            float x0 = xi + xNMask;
            float y0 = yi + yNMask;
            float z0 = zi + zNMask;
            float a0 = 0.75f - x0 * x0 - y0 * y0 - z0 * z0;
            float value = (a0 * a0) * (a0 * a0) * Grad(seed,
                xrbp + (xNMask & PRIME_X), yrbp + (yNMask & PRIME_Y), zrbp + (zNMask & PRIME_Z), x0, y0, z0);

            float x1 = xi - 0.5f;
            float y1 = yi - 0.5f;
            float z1 = zi - 0.5f;
            float a1 = 0.75f - x1 * x1 - y1 * y1 - z1 * z1;
            value += (a1 * a1) * (a1 * a1) * Grad(seed2,
                xrbp + PRIME_X, yrbp + PRIME_Y, zrbp + PRIME_Z, x1, y1, z1);

            float xAFlipMask0 = ((xNMask | 1) << 1) * x1;
            float yAFlipMask0 = ((yNMask | 1) << 1) * y1;
            float zAFlipMask0 = ((zNMask | 1) << 1) * z1;
            float xAFlipMask1 = (-2 - (xNMask << 2)) * x1 - 1.0f;
            float yAFlipMask1 = (-2 - (yNMask << 2)) * y1 - 1.0f;
            float zAFlipMask1 = (-2 - (zNMask << 2)) * z1 - 1.0f;

            bool skip5 = false;
            float a2 = xAFlipMask0 + a0;
            if (a2 > 0)
            {
                float x2 = x0 - (xNMask | 1);
                float y2 = y0;
                float z2 = z0;
                value += (a2 * a2) * (a2 * a2) * Grad(seed,
                    xrbp + (~xNMask & PRIME_X), yrbp + (yNMask & PRIME_Y), zrbp + (zNMask & PRIME_Z), x2, y2, z2);
            }
            else
            {
                float a3 = yAFlipMask0 + zAFlipMask0 + a0;
                if (a3 > 0)
                {
                    float x3 = x0;
                    float y3 = y0 - (yNMask | 1);
                    float z3 = z0 - (zNMask | 1);
                    value += (a3 * a3) * (a3 * a3) * Grad(seed,
                        xrbp + (xNMask & PRIME_X), yrbp + (~yNMask & PRIME_Y), zrbp + (~zNMask & PRIME_Z), x3, y3, z3);
                }

                float a4 = xAFlipMask1 + a1;
                if (a4 > 0)
                {
                    float x4 = (xNMask | 1) + x1;
                    float y4 = y1;
                    float z4 = z1;
                    value += (a4 * a4) * (a4 * a4) * Grad(seed2,
                        xrbp + (xNMask & (PRIME_X * 2)), yrbp + PRIME_Y, zrbp + PRIME_Z, x4, y4, z4);
                    skip5 = true;
                }
            }

            bool skip9 = false;
            float a6 = yAFlipMask0 + a0;
            if (a6 > 0)
            {
                float x6 = x0;
                float y6 = y0 - (yNMask | 1);
                float z6 = z0;
                value += (a6 * a6) * (a6 * a6) * Grad(seed,
                    xrbp + (xNMask & PRIME_X), yrbp + (~yNMask & PRIME_Y), zrbp + (zNMask & PRIME_Z), x6, y6, z6);
            }
            else
            {
                float a7 = xAFlipMask0 + zAFlipMask0 + a0;
                if (a7 > 0)
                {
                    float x7 = x0 - (xNMask | 1);
                    float y7 = y0;
                    float z7 = z0 - (zNMask | 1);
                    value += (a7 * a7) * (a7 * a7) * Grad(seed,
                        xrbp + (~xNMask & PRIME_X), yrbp + (yNMask & PRIME_Y), zrbp + (~zNMask & PRIME_Z), x7, y7, z7);
                }

                float a8 = yAFlipMask1 + a1;
                if (a8 > 0)
                {
                    float x8 = x1;
                    float y8 = (yNMask | 1) + y1;
                    float z8 = z1;
                    value += (a8 * a8) * (a8 * a8) * Grad(seed2,
                        xrbp + PRIME_X, yrbp + (yNMask & (PRIME_Y << 1)), zrbp + PRIME_Z, x8, y8, z8);
                    skip9 = true;
                }
            }

            bool skipD = false;
            float aA = zAFlipMask0 + a0;
            if (aA > 0)
            {
                float xA = x0;
                float yA = y0;
                float zA = z0 - (zNMask | 1);
                value += (aA * aA) * (aA * aA) * Grad(seed,
                    xrbp + (xNMask & PRIME_X), yrbp + (yNMask & PRIME_Y), zrbp + (~zNMask & PRIME_Z), xA, yA, zA);
            }
            else
            {
                float aB = xAFlipMask0 + yAFlipMask0 + a0;
                if (aB > 0)
                {
                    float xB = x0 - (xNMask | 1);
                    float yB = y0 - (yNMask | 1);
                    float zB = z0;
                    value += (aB * aB) * (aB * aB) * Grad(seed,
                        xrbp + (~xNMask & PRIME_X), yrbp + (~yNMask & PRIME_Y), zrbp + (zNMask & PRIME_Z), xB, yB, zB);
                }

                float aC = zAFlipMask1 + a1;
                if (aC > 0)
                {
                    float xC = x1;
                    float yC = y1;
                    float zC = (zNMask | 1) + z1;
                    value += (aC * aC) * (aC * aC) * Grad(seed2,
                        xrbp + PRIME_X, yrbp + PRIME_Y, zrbp + (zNMask & (PRIME_Z << 1)), xC, yC, zC);
                    skipD = true;
                }
            }

            if (!skip5)
            {
                float a5 = yAFlipMask1 + zAFlipMask1 + a1;
                if (a5 > 0)
                {
                    float x5 = x1;
                    float y5 = (yNMask | 1) + y1;
                    float z5 = (zNMask | 1) + z1;
                    value += (a5 * a5) * (a5 * a5) * Grad(seed2,
                        xrbp + PRIME_X, yrbp + (yNMask & (PRIME_Y << 1)), zrbp + (zNMask & (PRIME_Z << 1)), x5, y5, z5);
                }
            }

            if (!skip9)
            {
                float a9 = xAFlipMask1 + zAFlipMask1 + a1;
                if (a9 > 0)
                {
                    float x9 = (xNMask | 1) + x1;
                    float y9 = y1;
                    float z9 = (zNMask | 1) + z1;
                    value += (a9 * a9) * (a9 * a9) * Grad(seed2,
                        xrbp + (xNMask & (PRIME_X * 2)), yrbp + PRIME_Y, zrbp + (zNMask & (PRIME_Z << 1)), x9, y9, z9);
                }
            }

            if (!skipD)
            {
                float aD = xAFlipMask1 + yAFlipMask1 + a1;
                if (aD > 0)
                {
                    float xD = (xNMask | 1) + x1;
                    float yD = (yNMask | 1) + y1;
                    float zD = z1;
                    value += (aD * aD) * (aD * aD) * Grad(seed2,
                        xrbp + (xNMask & (PRIME_X << 1)), yrbp + (yNMask & (PRIME_Y << 1)), zrbp + PRIME_Z, xD, yD, zD);
                }
            }

            return value;
        }

        /*
         * Utility
         */

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private static float Grad(long seed, int xrvp, int yrvp, float dx, float dy)
        {
            long hash = seed ^ xrvp ^ yrvp;
            hash *= 0x5555555555555555L;
            hash ^= hash >> 56;
            int gi = (int)hash & 0xFE;
            return GRADIENTS_2D[gi | 0] * dx + GRADIENTS_2D[gi | 1] * dy;
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private static float Grad(long seed, int xrvp, int yrvp, int zrvp, float dx, float dy, float dz)
        {
            long hash = (seed ^ xrvp) ^ (yrvp ^ zrvp);
            hash *= 0x5555555555555555L;
            hash ^= hash >> 54;
            int gi = (int)hash & 0x3FC;
            return GRADIENTS_3D[gi | 0] * dx + GRADIENTS_3D[gi | 1] * dy + GRADIENTS_3D[gi | 2] * dz;
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private static int FastFloor(double x)
        {
            int xi = (int)x;
            return x < xi ? xi - 1 : xi;
        }

        /*
         * Lookup Tables & Gradients
         */

        private const double N2 = 0.05481866495625118;
        private const double N3 = 0.2781926117527186;

        private const int N_GRADS_2D = 128;
        private const int N_GRADS_3D = 256;

        private static float[] GRADIENTS_2D;
        private static float[] GRADIENTS_3D;

        static OpenSimplex2S()
        {
            
            int[] grad2FillIndices = new int[] { 1, 4, 7, 10, 13, 16, 19, 22 };
            GRADIENTS_2D = new float[N_GRADS_2D * 2];
            float[] grad2 = {
                 0.130526192220052f,  0.99144486137381f,
                 0.38268343236509f,   0.923879532511287f,
                 0.608761429008721f,  0.793353340291235f,
                 0.793353340291235f,  0.608761429008721f,
                 0.923879532511287f,  0.38268343236509f,
                 0.99144486137381f,   0.130526192220051f,
                 0.99144486137381f,  -0.130526192220051f,
                 0.923879532511287f, -0.38268343236509f,
                 0.793353340291235f, -0.60876142900872f,
                 0.608761429008721f, -0.793353340291235f,
                 0.38268343236509f,  -0.923879532511287f,
                 0.130526192220052f, -0.99144486137381f,
                -0.130526192220052f, -0.99144486137381f,
                -0.38268343236509f,  -0.923879532511287f,
                -0.608761429008721f, -0.793353340291235f,
                -0.793353340291235f, -0.608761429008721f,
                -0.923879532511287f, -0.38268343236509f,
                -0.99144486137381f,  -0.130526192220052f,
                -0.99144486137381f,   0.130526192220051f,
                -0.923879532511287f,  0.38268343236509f,
                -0.793353340291235f,  0.608761429008721f,
                -0.608761429008721f,  0.793353340291235f,
                -0.38268343236509f,   0.923879532511287f,
                -0.130526192220052f,  0.99144486137381f
            };
            for (int i = 0; i < grad2.Length; i++)
            {
                grad2[i] = (float)(grad2[i] / N2);
            }
            {
                int i = 0;
                while (i < GRADIENTS_2D.Length - grad2.Length + 1)
                {
                    for (int j = 0; j < grad2.Length; j++)
                    {
                        GRADIENTS_2D[i + j] = grad2[j];
                    }
                    i += grad2.Length;
                }
                for (int f = 0; f < grad2FillIndices.Length; f++)
                {
                    int gi = grad2FillIndices[f] * 2;
                    GRADIENTS_2D[i | 0] = grad2[gi | 0];
                    GRADIENTS_2D[i | 1] = grad2[gi | 1];
                    i += 2;
                }
            }

            int[] grad3FillIndices = new int[] { 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 44, 45, 46, 47 };
            GRADIENTS_3D = new float[N_GRADS_3D * 4];
            float[] grad3 = {
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
                 2.22474487139f,       2.22474487139f,      -1.0f,                 0.0f,
                 2.22474487139f,       2.22474487139f,       1.0f,                 0.0f,
                 3.0862664687972017f,  1.1721513422464978f,  0.0f,                 0.0f,
                 1.1721513422464978f,  3.0862664687972017f,  0.0f,                 0.0f
            };
            for (int i = 0; i < grad3.Length; i++)
            {
                grad3[i] = (float)(grad3[i] / N3);
            }
            {
                int i = 0;
                while (i < GRADIENTS_3D.Length - grad3.Length + 1)
                {
                    for (int j = 0; j < grad3.Length; j++)
                    {
                        GRADIENTS_3D[i + j] = grad3[j];
                    }
                    i += grad3.Length;
                }
                for (int f = 0; f < grad3FillIndices.Length; f++)
                {
                    int gi = grad3FillIndices[f] * 4;
                    GRADIENTS_3D[i | 0] = grad3[gi | 0];
                    GRADIENTS_3D[i | 1] = grad3[gi | 1];
                    GRADIENTS_3D[i | 2] = grad3[gi | 2];
                    i += 4;
                }
            }
        }

    }
}
