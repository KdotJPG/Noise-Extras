using System.Runtime.CompilerServices;

public class RotatedPerlinNoise {
    private const long PRIME_X = 0x5205402B9270C86FL;
    private const long PRIME_Y = 0x598CD327003817B5L;
    private const long PRIME_Z = 0x5BCC226E9FA0BACBL;
    private const long HashMultiplier = 0x53A3F72DEEC546F5L;

    private const double Root3Over3 = 0.577350269189626;
    private const double RotationOrthogonalizer3D = -0.21132486540518713;

    private const double NormalizationDivisor3D = 2.742445288166158;

    public static double Noise3_ImproveXZ(long seed, double x, double y, double z) {
            
        // Rotation to conceal grid and tune noise for XZ horizontal and Y vertical (or time).
        double xPlusZ = x + z;
        double yRescaled = y * Root3Over3; // No skew here.
        double xzOffset = xPlusZ * RotationOrthogonalizer3D + yRescaled;
        x += xzOffset;
        z += xzOffset;
        y = yRescaled - xPlusZ * Root3Over3;

        int xBase = FastFloor(x), yBase = FastFloor(y), zBase = FastFloor(z);
        double xDelta = x - xBase, yDelta = y - yBase, zDelta = z - zBase;
        long xBasePrimed = xBase * PRIME_X, yBasePrimed = yBase * PRIME_Y, zBasePrimed = zBase * PRIME_Z;
        double inFadeX = Fade(xDelta), inFadeY = Fade(yDelta), inFadeZ = Fade(zDelta);
        double outFadeX = 1 - inFadeX, outFadeY = 1 - inFadeY, outFadeZ = 1 - inFadeZ;

        double grad000 = Gradient(seed, xBasePrimed,          yBasePrimed,            zBasePrimed,           xDelta,     yDelta,     zDelta    );
        double grad001 = Gradient(seed, xBasePrimed,          yBasePrimed,            zBasePrimed + PRIME_Z, xDelta,     yDelta,     zDelta - 1);
        double grad00Z = outFadeZ * grad000 + inFadeZ * grad001;
        double grad010 = Gradient(seed, xBasePrimed,           yBasePrimed + PRIME_Y, zBasePrimed,           xDelta,     yDelta - 1, zDelta    );
        double grad011 = Gradient(seed, xBasePrimed,           yBasePrimed + PRIME_Y, zBasePrimed + PRIME_Z, xDelta,     yDelta - 1, zDelta - 1);
        double grad01Z = outFadeZ * grad010 + inFadeZ * grad011;
        double grad0YZ = outFadeY * grad00Z + inFadeY * grad01Z;
        double grad100 = Gradient(seed, xBasePrimed + PRIME_X, yBasePrimed,           zBasePrimed,           xDelta - 1, yDelta,     zDelta    );
        double grad101 = Gradient(seed, xBasePrimed + PRIME_X, yBasePrimed,           zBasePrimed + PRIME_Z, xDelta - 1, yDelta,     zDelta - 1);
        double grad10Z = outFadeZ * grad100 + inFadeZ * grad101;
        double grad110 = Gradient(seed, xBasePrimed + PRIME_X, yBasePrimed + PRIME_Y, zBasePrimed,           xDelta - 1, yDelta - 1, zDelta    );
        double grad111 = Gradient(seed, xBasePrimed + PRIME_X, yBasePrimed + PRIME_Y, zBasePrimed + PRIME_Z, xDelta - 1, yDelta - 1, zDelta - 1);
        double grad11Z = outFadeZ * grad110 + inFadeZ * grad111;
        double grad1YZ = outFadeY * grad10Z + inFadeY * grad11Z;
        double gradXYZ = outFadeX * grad0YZ + inFadeX * grad1YZ;

        return gradXYZ;
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static double Gradient(long seed, long xPrimed, long yPrimed, long zPrimed, double xDelta, double yDelta, double zDelta) {
        long hash = seed ^ xPrimed ^ yPrimed ^ zPrimed;
        hash *= HashMultiplier;
        hash ^= hash >> (64 - Gradients3DBitsConsidered);

        int index = (int)hash & Gradients3DHashMask;
        index = (index * NGradients3DMultiplier) >> Gradients3DBitShift;
        index &= Gradients3DSelectionMask;

        var gradient = Gradients3D[index];
        return gradient.X * xDelta + gradient.Y * yDelta + gradient.Z * zDelta;
    }

    public static void Noise3_ImproveXZ(Span<long> seeds, Span<double> destination, double x, double y, double z) {

        // Rotation to conceal grid and tune noise for XZ horizontal and Y vertical (or time).
        double xPlusZ = x + z;
        double yRescaled = y * Root3Over3; // No skew here.
        double xzOffset = xPlusZ * RotationOrthogonalizer3D + yRescaled;
        x += xzOffset;
        z += xzOffset;
        y = yRescaled - xPlusZ * Root3Over3;

        int xBase = FastFloor(x), yBase = FastFloor(y), zBase = FastFloor(z);
        double xDelta = x - xBase, yDelta = y - yBase, zDelta = z - zBase;
        long xBasePrimed = xBase * PRIME_X, yBasePrimed = yBase * PRIME_Y, zBasePrimed = zBase * PRIME_Z;
        double inFadeX = Fade(xDelta), inFadeY = Fade(yDelta), inFadeZ = Fade(zDelta);
        double outFadeX = 1 - inFadeX, outFadeY = 1 - inFadeY, outFadeZ = 1 - inFadeZ;

        for (int i = 0; i < seeds.Length; i++) destination[i] = 0;
        Gradient(seeds, destination, xBasePrimed,           yBasePrimed,           zBasePrimed,           xDelta,     yDelta,     zDelta,     outFadeX * outFadeY * outFadeZ);
        Gradient(seeds, destination, xBasePrimed,           yBasePrimed,           zBasePrimed + PRIME_Z, xDelta,     yDelta,     zDelta - 1, outFadeX * outFadeY * inFadeZ );
        Gradient(seeds, destination, xBasePrimed,           yBasePrimed + PRIME_Y, zBasePrimed,           xDelta,     yDelta - 1, zDelta,     outFadeX * inFadeY  * outFadeZ);
        Gradient(seeds, destination, xBasePrimed,           yBasePrimed + PRIME_Y, zBasePrimed + PRIME_Z, xDelta,     yDelta - 1, zDelta - 1, outFadeX * inFadeY  * inFadeZ );
        Gradient(seeds, destination, xBasePrimed + PRIME_X, yBasePrimed,           zBasePrimed,           xDelta - 1, yDelta,     zDelta,     inFadeX  * outFadeY * outFadeZ);
        Gradient(seeds, destination, xBasePrimed + PRIME_X, yBasePrimed,           zBasePrimed + PRIME_Z, xDelta - 1, yDelta,     zDelta - 1, inFadeX  * outFadeY * inFadeZ );
        Gradient(seeds, destination, xBasePrimed + PRIME_X, yBasePrimed + PRIME_Y, zBasePrimed,           xDelta - 1, yDelta - 1, zDelta,     inFadeX  * inFadeY  * outFadeZ);
        Gradient(seeds, destination, xBasePrimed + PRIME_X, yBasePrimed + PRIME_Y, zBasePrimed + PRIME_Z, xDelta - 1, yDelta - 1, zDelta - 1, inFadeX  * inFadeY  * inFadeZ );
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static void Gradient(
            Span<long> seeds, Span<double> destination, long xPrimed, long yPrimed, long zPrimed, double xDelta, double yDelta, double zDelta, double falloff) {
        long vertexHashComponent = xPrimed ^ yPrimed ^ zPrimed;
        double xDeltaWithFalloff = xDelta * falloff, yDeltaWithFalloff = yDelta * falloff, zDeltaWithFalloff = zDelta * falloff;

        for (int i = 0; i < seeds.Length; i++) {
            long hash = seeds[i] ^ vertexHashComponent;
            hash *= HashMultiplier;
            hash ^= hash >> (64 - Gradients3DBitsConsidered);

            int index = (int)hash & Gradients3DHashMask;
            index = (index * NGradients3DMultiplier) >> Gradients3DBitShift;
            index &= Gradients3DSelectionMask;

            var gradient = Gradients3D[index];
            destination[i] += gradient.X * xDeltaWithFalloff + gradient.Y * yDeltaWithFalloff + gradient.Z * zDeltaWithFalloff;
        }
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static int FastFloor(double value) => value < (int)value ? (int)value - 1 : (int)value;

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static double Fade(double value) => value * value * value * (10 + value * (-15 + value * 6));

    private const int NGradients3DBaseExponent = 4;
    private const int NGradients3DMultiplier = 3;
    private const int NGradients3DMultiplierCoverageExponent = 2;
    private const int NGradients3D = (1 << NGradients3DBaseExponent) * NGradients3DMultiplier;
    private const int Gradients3DBitsConsidered = 31 - NGradients3DMultiplierCoverageExponent;
    private const int Gradients3DHashMask = (1 << Gradients3DBitsConsidered) - 1;
    private const int Gradients3DBitShift = Gradients3DBitsConsidered - NGradients3DBaseExponent;
    private const int Gradients3DSelectionMask = (1 << (NGradients3DBaseExponent + NGradients3DMultiplierCoverageExponent)) - 1;

    private const double N3Div = NormalizationDivisor3D;
    private static readonly (double X, double Y, double Z)[] Gradients3D = new (double X, double Y, double Z)[NGradients3D] {
        ( 2.22474487139 / N3Div,       2.22474487139 / N3Div,      -1.0 / N3Div                ),
        ( 2.22474487139 / N3Div,       2.22474487139 / N3Div,       1.0 / N3Div                ),
        ( 3.0862664687972017 / N3Div,  1.1721513422464978 / N3Div,  0.0 / N3Div                ),
        ( 1.1721513422464978 / N3Div,  3.0862664687972017 / N3Div,  0.0 / N3Div                ),
        (-2.22474487139 / N3Div,       2.22474487139 / N3Div,      -1.0 / N3Div                ),
        (-2.22474487139 / N3Div,       2.22474487139 / N3Div,       1.0 / N3Div                ),
        (-1.1721513422464978 / N3Div,  3.0862664687972017 / N3Div,  0.0 / N3Div                ),
        (-3.0862664687972017 / N3Div,  1.1721513422464978 / N3Div,  0.0 / N3Div                ),
        (-1.0 / N3Div,                -2.22474487139 / N3Div,      -2.22474487139 / N3Div      ),
        ( 1.0 / N3Div,                -2.22474487139 / N3Div,      -2.22474487139 / N3Div      ),
        ( 0.0 / N3Div,                -3.0862664687972017 / N3Div, -1.1721513422464978 / N3Div ),
        ( 0.0 / N3Div,                -1.1721513422464978 / N3Div, -3.0862664687972017 / N3Div ),
        (-1.0 / N3Div,                -2.22474487139 / N3Div,       2.22474487139 / N3Div      ),
        ( 1.0 / N3Div,                -2.22474487139 / N3Div,       2.22474487139 / N3Div      ),
        ( 0.0 / N3Div,                -1.1721513422464978 / N3Div,  3.0862664687972017 / N3Div ),
        ( 0.0 / N3Div,                -3.0862664687972017 / N3Div,  1.1721513422464978 / N3Div ),
        (-2.22474487139 / N3Div,      -2.22474487139 / N3Div,      -1.0 / N3Div                ),
        (-2.22474487139 / N3Div,      -2.22474487139 / N3Div,       1.0 / N3Div                ),
        (-3.0862664687972017 / N3Div, -1.1721513422464978 / N3Div,  0.0 / N3Div                ),
        (-1.1721513422464978 / N3Div, -3.0862664687972017 / N3Div,  0.0 / N3Div                ),
        (-2.22474487139 / N3Div,      -1.0 / N3Div,                -2.22474487139 / N3Div      ),
        (-2.22474487139 / N3Div,       1.0 / N3Div,                -2.22474487139 / N3Div      ),
        (-1.1721513422464978 / N3Div,  0.0 / N3Div,                -3.0862664687972017 / N3Div ),
        (-3.0862664687972017 / N3Div,  0.0 / N3Div,                -1.1721513422464978 / N3Div ),
        (-2.22474487139 / N3Div,      -1.0 / N3Div,                 2.22474487139 / N3Div      ),
        (-2.22474487139 / N3Div,       1.0 / N3Div,                 2.22474487139 / N3Div      ),
        (-3.0862664687972017 / N3Div,  0.0 / N3Div,                 1.1721513422464978 / N3Div ),
        (-1.1721513422464978 / N3Div,  0.0 / N3Div,                 3.0862664687972017 / N3Div ),
        (-1.0 / N3Div,                 2.22474487139 / N3Div,      -2.22474487139 / N3Div      ),
        ( 1.0 / N3Div,                 2.22474487139 / N3Div,      -2.22474487139 / N3Div      ),
        ( 0.0 / N3Div,                 1.1721513422464978 / N3Div, -3.0862664687972017 / N3Div ),
        ( 0.0 / N3Div,                 3.0862664687972017 / N3Div, -1.1721513422464978 / N3Div ),
        (-1.0 / N3Div,                 2.22474487139 / N3Div,       2.22474487139 / N3Div      ),
        ( 1.0 / N3Div,                 2.22474487139 / N3Div,       2.22474487139 / N3Div      ),
        ( 0.0 / N3Div,                 3.0862664687972017 / N3Div,  1.1721513422464978 / N3Div ),
        ( 0.0 / N3Div,                 1.1721513422464978 / N3Div,  3.0862664687972017 / N3Div ),
        ( 2.22474487139 / N3Div,      -2.22474487139 / N3Div,      -1.0 / N3Div                ),
        ( 2.22474487139 / N3Div,      -2.22474487139 / N3Div,       1.0 / N3Div                ),
        ( 1.1721513422464978 / N3Div, -3.0862664687972017 / N3Div,  0.0 / N3Div                ),
        ( 3.0862664687972017 / N3Div, -1.1721513422464978 / N3Div,  0.0 / N3Div                ),
        ( 2.22474487139 / N3Div,      -1.0 / N3Div,                -2.22474487139 / N3Div      ),
        ( 2.22474487139 / N3Div,       1.0 / N3Div,                -2.22474487139 / N3Div      ),
        ( 3.0862664687972017 / N3Div,  0.0 / N3Div,                -1.1721513422464978 / N3Div ),
        ( 1.1721513422464978 / N3Div,  0.0 / N3Div,                -3.0862664687972017 / N3Div ),
        ( 2.22474487139 / N3Div,      -1.0 / N3Div,                 2.22474487139 / N3Div      ),
        ( 2.22474487139 / N3Div,       1.0 / N3Div,                 2.22474487139 / N3Div      ),
        ( 1.1721513422464978 / N3Div,  0.0 / N3Div,                 3.0862664687972017 / N3Div ),
        ( 3.0862664687972017 / N3Div,  0.0 / N3Div,                 1.1721513422464978 / N3Div ),
    };
}
