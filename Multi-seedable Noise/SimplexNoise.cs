using System.Runtime.CompilerServices;

public static class SimplexNoise {
    private const long PRIME_X = 0x5205402B9270C86FL;
    private const long PRIME_Y = 0x598CD327003817B5L;
    private const long PRIME_Z = 0x5BCC226E9FA0BACBL;
    private const long HashMultiplier = 0x53A3F72DEEC546F5L;

    private const double Skew3D = 1.0 / 3.0;
    private const double Unskew3D = -1.0 / 2.0; // In place of `-1.0 / 6.0`. Performs a reflected unskew to match the gradient table.

    private const double Root3Over3 = 0.577350269189626;
    private const double RotationOrthogonalizer3D = -0.21132486540518713;

    private const double FalloffRadiusSquared = 0.6;
    private const double NormalizationDivisor3D = 0.07969837668935331;

    public static double Noise3_ImproveXZ(long seed, double x, double y, double z) {

        // Combined skew + rotation to tune noise for XZ horizontal and Y vertical (or time).
        // In place of `(xSkewed, ySkewed, zSkewed) = (x, y, z) + (x + y + z) * Skew3D`
        double xPlusZ = x + z;
        double yRescaled = y * (Root3Over3 * 2.0); // This *2.0 performs the skew.
        double xzOffset = xPlusZ * RotationOrthogonalizer3D + yRescaled;
        double xSkewed = x + xzOffset;
        double ySkewed = yRescaled - xPlusZ * Root3Over3;
        double zSkewed = z + xzOffset;

        int xSkewedBase = FastFloor(xSkewed), ySkewedBase = FastFloor(ySkewed), zSkewedBase = FastFloor(zSkewed);
        long xBasePrimed = xSkewedBase * PRIME_X, yBasePrimed = ySkewedBase * PRIME_Y, zBasePrimed = zSkewedBase * PRIME_Z;

        double xDeltaSkewed = xSkewed - xSkewedBase, yDeltaSkewed = ySkewed - ySkewedBase, zDeltaSkewed = zSkewed - zSkewedBase;
        double unskewOffset = (xDeltaSkewed + yDeltaSkewed + zDeltaSkewed) * Unskew3D;
        double xDelta000 = xDeltaSkewed + unskewOffset, yDelta000 = yDeltaSkewed + unskewOffset, zDelta000 = zDeltaSkewed + unskewOffset;
        double distanceSquared000 = xDelta000 * xDelta000 + yDelta000 * yDelta000 + zDelta000 * zDelta000;
        double xDelta111 = xDelta000 - (1 + 3 * Unskew3D), yDelta111 = yDelta000 - (1 + 3 * Unskew3D), zDelta111 = zDelta000 - (1 + 3 * Unskew3D);
        double distanceSquared111 = xDelta111 * xDelta111 + yDelta111 * yDelta111 + zDelta111 * zDelta111;

        double value = 0;
        if (distanceSquared000 < FalloffRadiusSquared) {
            value += Contribution(seed, distanceSquared000, xBasePrimed, yBasePrimed, zBasePrimed, xDelta000, yDelta000, zDelta000);
        }
        if (distanceSquared111 < FalloffRadiusSquared) {
            value += Contribution(seed, distanceSquared111, xBasePrimed + PRIME_X, yBasePrimed + PRIME_Y, zBasePrimed + PRIME_Z, xDelta111, yDelta111, zDelta111);
        }

        bool xGreaterThanY = (xDeltaSkewed > yDeltaSkewed);
        bool xGreaterThanZ = (xDeltaSkewed > zDeltaSkewed);
        bool yGreaterThanZ = (yDeltaSkewed > zDeltaSkewed);

        if (xGreaterThanY && xGreaterThanZ) {
            double xDelta100 = xDelta000 - (1 + Unskew3D), yDelta100 = yDelta000 - Unskew3D, zDelta100 = zDelta000 - Unskew3D;
            double distanceSquared100 = xDelta100 * xDelta100 + yDelta100 * yDelta100 + zDelta100 * zDelta100;
            if (distanceSquared100 < FalloffRadiusSquared) {
                value += Contribution(seed, distanceSquared100, xBasePrimed + PRIME_X, yBasePrimed, zBasePrimed, xDelta100, yDelta100, zDelta100);
            }
        } else if (!xGreaterThanY && yGreaterThanZ) {
            double xDelta010 = xDelta000 - Unskew3D, yDelta010 = yDelta000 - (1 + Unskew3D), zDelta010 = zDelta000 - Unskew3D;
            double distanceSquared010 = xDelta010 * xDelta010 + yDelta010 * yDelta010 + zDelta010 * zDelta010;
            if (distanceSquared010 < FalloffRadiusSquared) {
                value += Contribution(seed, distanceSquared010, xBasePrimed, yBasePrimed + PRIME_Y, zBasePrimed, xDelta010, yDelta010, zDelta010);
            }
        } else {
            double xDelta001 = xDelta000 - Unskew3D, yDelta001 = yDelta000 - Unskew3D, zDelta001 = zDelta000 - (1 + Unskew3D);
            double distanceSquared001 = xDelta001 * xDelta001 + yDelta001 * yDelta001 + zDelta001 * zDelta001;
            if (distanceSquared001 < FalloffRadiusSquared) {
                value += Contribution(seed, distanceSquared001, xBasePrimed, yBasePrimed, zBasePrimed + PRIME_Z, xDelta001, yDelta001, zDelta001);
            }
        }

        if (xGreaterThanZ && yGreaterThanZ) {
            double xDelta110 = xDelta000 - (1 + 2 * Unskew3D), yDelta110 = yDelta000 - (1 + 2 * Unskew3D), zDelta110 = zDelta000 - 2 * Unskew3D;
            double distanceSquared110 = xDelta110 * xDelta110 + yDelta110 * yDelta110 + zDelta110 * zDelta110;
            if (distanceSquared110 < FalloffRadiusSquared) {
                value += Contribution(seed, distanceSquared110, xBasePrimed + PRIME_X, yBasePrimed + PRIME_Y, zBasePrimed, xDelta110, yDelta110, zDelta110);
            }
        } else if (xGreaterThanY && !yGreaterThanZ) {
            double xDelta101 = xDelta000 - (1 + 2 * Unskew3D), yDelta101 = yDelta000 - 2 * Unskew3D, zDelta101 = zDelta000 - (1 + 2 * Unskew3D);
            double distanceSquared101 = xDelta101 * xDelta101 + yDelta101 * yDelta101 + zDelta101 * zDelta101;
            if (distanceSquared101 < FalloffRadiusSquared) {
                value += Contribution(seed, distanceSquared101, xBasePrimed + PRIME_X, yBasePrimed, zBasePrimed + PRIME_Z, xDelta101, yDelta101, zDelta101);
            }
        } else {
            double xDelta011 = xDelta000 - 2 * Unskew3D, yDelta011 = yDelta000 - (1 + 2 * Unskew3D), zDelta011 = zDelta000 - (1 + 2 * Unskew3D);
            double distanceSquared011 = xDelta011 * xDelta011 + yDelta011 * yDelta011 + zDelta011 * zDelta011;
            if (distanceSquared011 < FalloffRadiusSquared) {
                value += Contribution(seed, distanceSquared011, xBasePrimed, yBasePrimed + PRIME_Y, zBasePrimed + PRIME_Z, xDelta011, yDelta011, zDelta011);
            }
        }

        return value;
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static double Contribution(
            long seed, double distanceSquared, long xPrimed, long yPrimed, long zPrimed, double xDelta, double yDelta, double zDelta) {
        long hash = seed ^ xPrimed ^ yPrimed ^ zPrimed;
        hash *= HashMultiplier;
        hash ^= hash >> (64 - Gradients3DBitsConsidered);

        int index = (int)hash & Gradients3DHashMask;
        index = (index * NGradients3DMultiplier) >> Gradients3DBitShift;
        index &= Gradients3DSelectionMask;

        double falloff = Pow4(distanceSquared - FalloffRadiusSquared);
        var gradient = Gradients3D[index];
        return falloff * (gradient.X * xDelta + gradient.Y * yDelta + gradient.Z * zDelta);
    }

    public static void Noise3_ImproveXZ(Span<long> seeds, Span<double> destination, double x, double y, double z) {

        // Combined skew + rotation to tune noise for XZ horizontal and Y vertical (or time).
        // In place of `(xSkewed, ySkewed, zSkewed) = (x, y, z) + (x + y + z) * Skew3D`
        double xPlusZ = x + z;
        double yRescaled = y * (Root3Over3 * 2.0); // This *2.0 performs the skew.
        double xzOffset = xPlusZ * RotationOrthogonalizer3D + yRescaled;
        double xSkewed = x + xzOffset;
        double ySkewed = yRescaled - xPlusZ * Root3Over3;
        double zSkewed = z + xzOffset;

        int xSkewedBase = FastFloor(xSkewed), ySkewedBase = FastFloor(ySkewed), zSkewedBase = FastFloor(zSkewed);
        long xBasePrimed = xSkewedBase * PRIME_X, yBasePrimed = ySkewedBase * PRIME_Y, zBasePrimed = zSkewedBase * PRIME_Z;

        for (int i = 0; i < seeds.Length; i++) destination[i] = 0;
        double xDeltaSkewed = xSkewed - xSkewedBase, yDeltaSkewed = ySkewed - ySkewedBase, zDeltaSkewed = zSkewed - zSkewedBase;
        double unskewOffset = (xDeltaSkewed + yDeltaSkewed + zDeltaSkewed) * Unskew3D;
        double xDelta000 = xDeltaSkewed + unskewOffset, yDelta000 = yDeltaSkewed + unskewOffset, zDelta000 = zDeltaSkewed + unskewOffset;
        double distanceSquared000 = xDelta000 * xDelta000 + yDelta000 * yDelta000 + zDelta000 * zDelta000;
        double xDelta111 = xDelta000 - (1 + 3 * Unskew3D), yDelta111 = yDelta000 - (1 + 3 * Unskew3D), zDelta111 = zDelta000 - (1 + 3 * Unskew3D);
        double distanceSquared111 = xDelta111 * xDelta111 + yDelta111 * yDelta111 + zDelta111 * zDelta111;

        if (distanceSquared000 < FalloffRadiusSquared) {
            Contribution(seeds, destination, distanceSquared000, xBasePrimed, yBasePrimed, zBasePrimed, xDelta000, yDelta000, zDelta000);
        }
        if (distanceSquared111 < FalloffRadiusSquared) {
            Contribution(seeds, destination, distanceSquared111, xBasePrimed + PRIME_X, yBasePrimed + PRIME_Y, zBasePrimed + PRIME_Z, xDelta111, yDelta111, zDelta111);
        }

        bool xGreaterThanY = (xDeltaSkewed > yDeltaSkewed);
        bool xGreaterThanZ = (xDeltaSkewed > zDeltaSkewed);
        bool yGreaterThanZ = (yDeltaSkewed > zDeltaSkewed);

        if (xGreaterThanY && xGreaterThanZ) {
            double xDelta100 = xDelta000 - (1 + Unskew3D), yDelta100 = yDelta000 - Unskew3D, zDelta100 = zDelta000 - Unskew3D;
            double distanceSquared100 = xDelta100 * xDelta100 + yDelta100 * yDelta100 + zDelta100 * zDelta100;
            if (distanceSquared100 < FalloffRadiusSquared) {
                Contribution(seeds, destination, distanceSquared100, xBasePrimed + PRIME_X, yBasePrimed, zBasePrimed, xDelta100, yDelta100, zDelta100);
            }
        } else if (!xGreaterThanY && yGreaterThanZ) {
            double xDelta010 = xDelta000 - Unskew3D, yDelta010 = yDelta000 - (1 + Unskew3D), zDelta010 = zDelta000 - Unskew3D;
            double distanceSquared010 = xDelta010 * xDelta010 + yDelta010 * yDelta010 + zDelta010 * zDelta010;
            if (distanceSquared010 < FalloffRadiusSquared) {
                Contribution(seeds, destination, distanceSquared010, xBasePrimed, yBasePrimed + PRIME_Y, zBasePrimed, xDelta010, yDelta010, zDelta010);
            }
        } else {
            double xDelta001 = xDelta000 - Unskew3D, yDelta001 = yDelta000 - Unskew3D, zDelta001 = zDelta000 - (1 + Unskew3D);
            double distanceSquared001 = xDelta001 * xDelta001 + yDelta001 * yDelta001 + zDelta001 * zDelta001;
            if (distanceSquared001 < FalloffRadiusSquared) {
                Contribution(seeds, destination, distanceSquared001, xBasePrimed, yBasePrimed, zBasePrimed + PRIME_Z, xDelta001, yDelta001, zDelta001);
            }
        }

        if (xGreaterThanZ && yGreaterThanZ) {
            double xDelta110 = xDelta000 - (1 + 2 * Unskew3D), yDelta110 = yDelta000 - (1 + 2 * Unskew3D), zDelta110 = zDelta000 - 2 * Unskew3D;
            double distanceSquared110 = xDelta110 * xDelta110 + yDelta110 * yDelta110 + zDelta110 * zDelta110;
            if (distanceSquared110 < FalloffRadiusSquared) {
                Contribution(seeds, destination, distanceSquared110, xBasePrimed + PRIME_X, yBasePrimed + PRIME_Y, zBasePrimed, xDelta110, yDelta110, zDelta110);
            }
        } else if (xGreaterThanY && !yGreaterThanZ) {
            double xDelta101 = xDelta000 - (1 + 2 * Unskew3D), yDelta101 = yDelta000 - 2 * Unskew3D, zDelta101 = zDelta000 - (1 + 2 * Unskew3D);
            double distanceSquared101 = xDelta101 * xDelta101 + yDelta101 * yDelta101 + zDelta101 * zDelta101;
            if (distanceSquared101 < FalloffRadiusSquared) {
                Contribution(seeds, destination, distanceSquared101, xBasePrimed + PRIME_X, yBasePrimed, zBasePrimed + PRIME_Z, xDelta101, yDelta101, zDelta101);
            }
        } else {
            double xDelta011 = xDelta000 - 2 * Unskew3D, yDelta011 = yDelta000 - (1 + 2 * Unskew3D), zDelta011 = zDelta000 - (1 + 2 * Unskew3D);
            double distanceSquared011 = xDelta011 * xDelta011 + yDelta011 * yDelta011 + zDelta011 * zDelta011;
            if (distanceSquared011 < FalloffRadiusSquared) {
                Contribution(seeds, destination, distanceSquared011, xBasePrimed, yBasePrimed + PRIME_Y, zBasePrimed + PRIME_Z, xDelta011, yDelta011, zDelta011);
            }
        }
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static void Contribution(
            Span<long> seeds, Span<double> destination, double distanceSquared, long xPrimed, long yPrimed, long zPrimed, double xDelta, double yDelta, double zDelta) {
        long vertexHashComponent = xPrimed ^ yPrimed ^ zPrimed;
        double falloff = Pow4(distanceSquared - FalloffRadiusSquared);

        for (int i = 0; i < seeds.Length; i++) {
            long hash = seeds[i] ^ vertexHashComponent;
            hash *= HashMultiplier;
            hash ^= hash >> (64 - Gradients3DBitsConsidered);

            int index = (int)hash & Gradients3DHashMask;
            index = (index * NGradients3DMultiplier) >> Gradients3DBitShift;
            index &= Gradients3DSelectionMask;

            var gradient = Gradients3D[index];
            destination[i] += falloff * (gradient.X * xDelta + gradient.Y * yDelta + gradient.Z * zDelta);
        }
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static int FastFloor(double value) => value < (int)value ? (int)value - 1 : (int)value;

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static double Pow4(double value) => (value * value) * (value * value);

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
