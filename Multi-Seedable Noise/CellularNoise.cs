using System.Runtime.CompilerServices;

public static class CellularNoise {
    private const long PRIME_X = 0x5205402B9270C86FL;
    private const long PRIME_Y = 0x598CD327003817B5L;
    private const long PRIME_Z = 0x5BCC226E9FA0BACBL;
    private const long HashMultiplier = 0x53A3F72DEEC546F5L;

    private const double Root3Over3 = 0.577350269189626;
    private const double RotationOrthogonalizer3D = -0.21132486540518713;

    private const int DisplacementHashBitsPerAxis = 21;
    private const int DisplacementHashMask = (1 << DisplacementHashBitsPerAxis) - 1;
    private const double DisplacementHashRescale = 1.0 / (1 << DisplacementHashBitsPerAxis);

    public const double MaxDistanceSquaredToClosest = 3.0;

    private const int MaxCellAxisDelta = 2;
    private const int CellAxisDeltaSpan = 2 * MaxCellAxisDelta + 1;

    public static double Noise3_ImproveXZ(long seed, double x, double y, double z) {

        // Rotation to further conceal grid and tune noise for XZ horizontal and Y vertical (or time).
        double xPlusZ = x + z;
        double yRescaled = y * Root3Over3; // No skew here.
        double xzOffset = xPlusZ * RotationOrthogonalizer3D + yRescaled;
        x += xzOffset;
        z += xzOffset;
        y = yRescaled - xPlusZ * Root3Over3;

        int xBase = ExtraMath.FastFloor(x), yBase = ExtraMath.FastFloor(y), zBase = ExtraMath.FastFloor(z);
        double xDelta = x - xBase, yDelta = y - yBase, zDelta = z - zBase;
        long xBasePrimed = xBase * PRIME_X, yBasePrimed = yBase * PRIME_Y, zBasePrimed = zBase * PRIME_Z;

        // Axis components for calculating minimum distances to the grid cells. Used for skip conditions.
        Span<double> cellDeltasSquared = stackalloc double[] {
            (xDelta + 1) * (xDelta + 1), xDelta * xDelta, 0, (xDelta - 1) * (xDelta - 1), (xDelta - 2) * (xDelta - 2),
            (yDelta + 1) * (yDelta + 1), yDelta * yDelta, 0, (yDelta - 1) * (yDelta - 1), (yDelta - 2) * (yDelta - 2),
            (zDelta + 1) * (zDelta + 1), zDelta * zDelta, 0, (zDelta - 1) * (zDelta - 1), (zDelta - 2) * (zDelta - 2),
        };

        double distanceSquaredToClosest = MaxDistanceSquaredToClosest;
        bool groupHadCellsInRange = true;

        for (int cellIndex = 0; cellIndex < Cells.Length; cellIndex++) {
            ref readonly Cell cell = ref Cells[cellIndex]; // Using `ref` here seems to drastically improve performance.

            // If the last group didn't have any cells in range, then we know none of the rest will either.
            if (cell.IsGroupStart) {
                if (!groupHadCellsInRange) break;
                groupHadCellsInRange = false;
            }

            // Check if the cell is close enough to contain a point closer than the closest one we've found so far.
            if (distanceSquaredToClosest - cellDeltasSquared[cell.CellDeltaSquaredIndexZ] <=
                cellDeltasSquared[cell.CellDeltaSquaredIndexX] + cellDeltasSquared[cell.CellDeltaSquaredIndexY]) continue;
            groupHadCellsInRange = true;

            double distanceSquaredHere = CalculateDistanceSquared(seed,
                xBasePrimed + cell.PrimeRelativeX, yBasePrimed + cell.PrimeRelativeY, zBasePrimed + cell.PrimeRelativeZ,
                xDelta + cell.DeltaRelativeX, yDelta + cell.DeltaRelativeY, zDelta + cell.DeltaRelativeZ);

            if (distanceSquaredHere < distanceSquaredToClosest) {
                distanceSquaredToClosest = distanceSquaredHere;
            }
        }

        return distanceSquaredToClosest;
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static double CalculateDistanceSquared(long seed, long xPrimed, long yPrimed, long zPrimed, double xDelta, double yDelta, double zDelta) {
        long hash = seed ^ xPrimed ^ yPrimed ^ zPrimed;
        hash *= HashMultiplier;

        double xDeltaDisplaced = xDelta - (((hash >> 43) ^ hash) & DisplacementHashMask) * DisplacementHashRescale;
        double yDeltaDisplaced = yDelta - (((hash >> 43) ^ (hash >> 22)) & DisplacementHashMask) * DisplacementHashRescale;
        double zDeltaDisplaced = zDelta - ((hash >> 43) & DisplacementHashMask) * DisplacementHashRescale;

        return xDeltaDisplaced * xDeltaDisplaced + yDeltaDisplaced * yDeltaDisplaced + zDeltaDisplaced * zDeltaDisplaced;
    }

    public static void Noise3_ImproveXZ(Span<long> seeds, Span<double> destination, double x, double y, double z) {

        // Rotation to further conceal grid and tune noise for XZ horizontal and Y vertical (or time).
        double xPlusZ = x + z;
        double yRescaled = y * Root3Over3; // No skew here.
        double xzOffset = xPlusZ * RotationOrthogonalizer3D + yRescaled;
        x += xzOffset;
        z += xzOffset;
        y = yRescaled - xPlusZ * Root3Over3;

        int xBase = ExtraMath.FastFloor(x), yBase = ExtraMath.FastFloor(y), zBase = ExtraMath.FastFloor(z);
        double xDelta = x - xBase, yDelta = y - yBase, zDelta = z - zBase;
        long xBasePrimed = xBase * PRIME_X, yBasePrimed = yBase * PRIME_Y, zBasePrimed = zBase * PRIME_Z;

        // Axis components for calculating minimum distances to the grid cells
        Span<double> cellDeltasSquared = stackalloc double[] {
            (xDelta + 1) * (xDelta + 1), xDelta * xDelta, 0, (xDelta - 1) * (xDelta - 1), (xDelta - 2) * (xDelta - 2),
            (yDelta + 1) * (yDelta + 1), yDelta * yDelta, 0, (yDelta - 1) * (yDelta - 1), (yDelta - 2) * (yDelta - 2),
            (zDelta + 1) * (zDelta + 1), zDelta * zDelta, 0, (zDelta - 1) * (zDelta - 1), (zDelta - 2) * (zDelta - 2),
        };

        for (int i = 0; i < seeds.Length; i++) destination[i] = MaxDistanceSquaredToClosest;
        double maxDistanceSquaredToClosest = MaxDistanceSquaredToClosest;
        bool groupHadCellsInRange = true;

        for (int cellIndex = 0; cellIndex < Cells.Length; cellIndex++) {
            ref readonly Cell cell = ref Cells[cellIndex]; // Using `ref` here seems to make a huge difference.

            // If the last group didn't have any cells in range, then we know none of the rest will either.
            if (cell.IsGroupStart) {
                if (!groupHadCellsInRange) break;
                groupHadCellsInRange = false;
            }

            // Check if the cell is close enough to contain a point closer than the closest one we've found so far.
            if (maxDistanceSquaredToClosest - cellDeltasSquared[cell.CellDeltaSquaredIndexZ] <=
                cellDeltasSquared[cell.CellDeltaSquaredIndexX] + cellDeltasSquared[cell.CellDeltaSquaredIndexY]) continue;
            groupHadCellsInRange = true;

            maxDistanceSquaredToClosest = CalculateAndUpdateDistanceSquared(seeds, destination,
                xBasePrimed + cell.PrimeRelativeX, yBasePrimed + cell.PrimeRelativeY, zBasePrimed + cell.PrimeRelativeZ,
                xDelta + cell.DeltaRelativeX, yDelta + cell.DeltaRelativeY, zDelta + cell.DeltaRelativeZ);
        }
    }
    
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static double CalculateAndUpdateDistanceSquared(
            Span<long> seeds, Span<double> destination, long xPrimed, long yPrimed, long zPrimed, double xDelta, double yDelta, double zDelta) {
        long cellHashComponent = xPrimed ^ yPrimed ^ zPrimed;
        double maxDistanceSquaredToClosest = 0;

        for (int i = 0; i < seeds.Length; i++) {
            long hash = seeds[i] ^ cellHashComponent;
            hash *= HashMultiplier;

            double xDeltaDisplaced = xDelta - (((hash >> 43) ^ hash) & DisplacementHashMask) * DisplacementHashRescale;
            double yDeltaDisplaced = yDelta - (((hash >> 43) ^ (hash >> 22)) & DisplacementHashMask) * DisplacementHashRescale;
            double zDeltaDisplaced = zDelta - ((hash >> 43) & DisplacementHashMask) * DisplacementHashRescale;

            double distanceSquared = xDeltaDisplaced * xDeltaDisplaced + yDeltaDisplaced * yDeltaDisplaced + zDeltaDisplaced * zDeltaDisplaced;
            if (distanceSquared < destination[i]) {
                destination[i] = distanceSquared;
            }

            if (destination[i] > maxDistanceSquaredToClosest) {
                maxDistanceSquaredToClosest = destination[i];
            }
        }

        return maxDistanceSquaredToClosest;
    }

    private static readonly Cell[] Cells;
    
    static CellularNoise() {

        Cells = new Cell[81];
        int cellIndexCurrent = 0;

        void defineCellGroup(int firstAxisValue, int secondAxisValue, int thirdAxisValue) {
            System.Diagnostics.Debug.Assert(firstAxisValue  >= secondAxisValue);
            System.Diagnostics.Debug.Assert(secondAxisValue >= thirdAxisValue);

            bool isFirstCell = true;
            Span<int> offsetVector = stackalloc int[3];
            for (int firstAxisIndex = 0; firstAxisIndex < 3; firstAxisIndex++) {
                offsetVector[firstAxisIndex] = firstAxisValue;

                for (int secondAxisIndex = 0; secondAxisIndex < 3; secondAxisIndex++) {
                    if (secondAxisIndex == firstAxisIndex) continue;
                    if (secondAxisIndex < firstAxisIndex && secondAxisValue == firstAxisValue) continue;
                    offsetVector[secondAxisIndex] = secondAxisValue;

                    int thirdAxisIndex = 3 - (firstAxisIndex + secondAxisIndex);
                    offsetVector[thirdAxisIndex] = thirdAxisValue;
                    if (thirdAxisIndex < secondAxisIndex && thirdAxisValue == secondAxisValue) continue;

                    for (int axisPositivityBits = 0; axisPositivityBits <= 0b111; axisPositivityBits++) {
                        bool firstAxisPositive  = (axisPositivityBits & 0b001) != 0;
                        bool secondAxisPositive = (axisPositivityBits & 0b010) != 0;
                        bool thirdAxisPositive  = (axisPositivityBits & 0b100) != 0;

                        if (firstAxisPositive  && firstAxisValue  == 0) continue;
                        if (secondAxisPositive && secondAxisValue == 0) continue;
                        if (thirdAxisPositive  && thirdAxisValue  == 0) continue;

                        int axisPositivityBitsResolved = 0;
                        if (firstAxisPositive)  axisPositivityBitsResolved |= (1 << firstAxisIndex);
                        if (secondAxisPositive) axisPositivityBitsResolved |= (1 << secondAxisIndex);
                        if (thirdAxisPositive)  axisPositivityBitsResolved |= (1 << thirdAxisIndex);

                        defineCell(
                            (axisPositivityBitsResolved & 0b001) != 0 ? offsetVector[0] : -offsetVector[0],
                            (axisPositivityBitsResolved & 0b010) != 0 ? offsetVector[1] : -offsetVector[1],
                            (axisPositivityBitsResolved & 0b100) != 0 ? offsetVector[2] : -offsetVector[2],
                            isFirstCell
                        );

                        isFirstCell = false;
                    }
                }
            }
        }

        void defineCell(int x, int y, int z, bool notGroupStart) {
            Cells[cellIndexCurrent] = new Cell(
                MaxCellAxisDelta + x, CellAxisDeltaSpan + MaxCellAxisDelta + y, 2 * CellAxisDeltaSpan + MaxCellAxisDelta + z,
                x * PRIME_X,y * PRIME_Y, z * PRIME_Z,
                -x, -y, -z,
                notGroupStart
            );
            cellIndexCurrent++;
        }
        
        defineCellGroup(0, 0, 0);
        defineCellGroup(1, 0, 0);
        defineCellGroup(1, 1, 0);
        defineCellGroup(1, 1, 1);
        defineCellGroup(2, 0, 0);
        defineCellGroup(2, 1, 0);
        defineCellGroup(2, 1, 1);

        System.Diagnostics.Debug.Assert(cellIndexCurrent == Cells.Length);
    }

    private readonly struct Cell {
        public readonly int CellDeltaSquaredIndexX, CellDeltaSquaredIndexY, CellDeltaSquaredIndexZ;
        public readonly long PrimeRelativeX, PrimeRelativeY, PrimeRelativeZ;
        public readonly double DeltaRelativeX, DeltaRelativeY, DeltaRelativeZ;
        public readonly bool IsGroupStart;

        public Cell(
            int cellDeltaSquaredIndexX, int cellDeltaSquaredIndexY, int cellDeltaSquaredIndexZ,
            long primeRelativeX, long primeRelativeY, long primeRelativeZ,
            double deltaRelativeX, double deltaRelativeY, double deltaRelativeZ,
            bool isGroupStart
        ) {
            CellDeltaSquaredIndexX = cellDeltaSquaredIndexX;
            CellDeltaSquaredIndexY = cellDeltaSquaredIndexY;
            CellDeltaSquaredIndexZ = cellDeltaSquaredIndexZ;
            PrimeRelativeX = primeRelativeX;
            PrimeRelativeY = primeRelativeY;
            PrimeRelativeZ = primeRelativeZ;
            DeltaRelativeX = deltaRelativeX;
            DeltaRelativeY = deltaRelativeY;
            DeltaRelativeZ = deltaRelativeZ;
            IsGroupStart = isGroupStart;
        }
    }
}
