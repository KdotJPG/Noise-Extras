public class DomainRotatedClassicNoise {

    private static final int PERMUTATION_TABLE_SIZE_EXPONENT = 11;

    private static final int PSIZE = 1 << PERMUTATION_TABLE_SIZE_EXPONENT;
    private static final int PMASK = PSIZE - 1;

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

    /*
     * 3D
     */

    public double noise2_from3(double x, double y) {
        return noise3_ImproveXY(x, y, 0);
    }

    public double noise3_ImproveXY(double x, double y, double z) {
        double xy = x + y;
        double s2 = xy * -0.211324865405187;
        double zz = z * 0.577350269189626;
        double xr = x + (s2 + zz); double yr = y + (s2 + zz);
        double zr = xy * -0.577350269189626 + zz;

        return noise3_UnrotatedBase(xr, yr, zr);
    }

    public double noise3_ImproveXY(Vector3 derivatives, double x, double y, double z) {
        double xy = x + y;
        double s2 = xy * -0.211324865405187;
        double zz = z * 0.577350269189626;
        double xr = x + (s2 + zz); double yr = y + (s2 + zz);
        double zr = xy * -0.577350269189626 + zz;

        double value = noise3_UnrotatedBase(derivatives, xr, yr, zr);

        double xyd = derivatives.x + derivatives.y;
        double s2d = xyd * -0.211324865405187;
        double zzd = derivatives.z * 0.577350269189626;
        derivatives.x = derivatives.x + (s2d - zzd);
        derivatives.y = derivatives.y + (s2d - zzd);
        derivatives.z = xyd * 0.577350269189626 + zzd;

        return value;
    }

    public double noise3_ImproveXZ(double x, double y, double z) {
        double xz = x + z;
        double s2 = xz * -0.211324865405187;
        double yy = y * 0.577350269189626;
        double xr = x + (s2 + yy); double zr = z + (s2 + yy);
        double yr = xz * -0.577350269189626 + yy;

        return noise3_UnrotatedBase(xr, yr, zr);
    }

    public double noise3_ImproveXZ(Vector3 derivatives, double x, double y, double z) {
        double xz = x + z;
        double s2 = xz * -0.211324865405187;
        double yy = y * 0.577350269189626;
        double xr = x + (s2 + yy); double zr = z + (s2 + yy);
        double yr = xz * -0.577350269189626 + yy;

        double value = noise3_UnrotatedBase(derivatives, xr, yr, zr);

        double xzd = derivatives.x + derivatives.z;
        double s2d = xzd * -0.211324865405187;
        double yyd = derivatives.y * 0.577350269189626;
        derivatives.x = derivatives.x + (s2d - yyd);
        derivatives.z = derivatives.z + (s2d - yyd);
        derivatives.y = xzd * 0.577350269189626 + yyd;

        return value;
    }

    private double noise3_UnrotatedBase(double xr, double yr, double zr) {

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

    private double noise3_UnrotatedBase(Vector3 derivatives, double xr, double yr, double zr) {

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
        derivatives.x = dXgXYZ;
        derivatives.y = dYgXYZ;
        derivatives.z = dZgXYZ;

        return gXYZ;
    }

    /*
     * 4D
     */

    public double noise3_from4(double x, double y, double z) {
        return noise4_ImproveXYZ(x, y, z, 0);
    }

    public double noise3_from4_ImproveXY(double x, double y, double z) {
        return noise4_ImproveXYZ_ImproveXY(x, y, z, 0);
    }

    public double noise3_from4_ImproveXZ(double x, double y, double z) {
        return noise4_ImproveXYZ_ImproveXZ(x, y, z, 0);
    }

    public double noise4_ImproveXYZ_ImproveXY(double x, double y, double z, double w) {

        double xy = x + y;
        double s2 = xy * -0.21132486540518699998;
        double zz = z * 0.28867513459481294226;
        double ww = w * 0.5;
        double xr = x + (zz + ww + s2), yr = y + (zz + ww + s2);
        double zr = xy * -0.57735026918962599998 + (zz + ww);
        double wr = z * -0.866025403784439 + ww;

        return noise4_UnrotatedBase(xr, yr, zr, wr);
    }

    public double noise4_ImproveXYZ_ImproveXY(Vector4 derivatives, double x, double y, double z, double w) {

        double xy = x + y;
        double s2 = xy * -0.21132486540518699998;
        double zz = z * 0.28867513459481294226;
        double ww = w * 0.5;
        double xr = x + (zz + ww + s2), yr = y + (zz + ww + s2);
        double zr = xy * -0.57735026918962599998 + (zz + ww);
        double wr = z * -0.866025403784439 + ww;

        double value = noise4_UnrotatedBase(derivatives, xr, yr, zr, wr);

        double xyd = derivatives.x + derivatives.y;
        double xyzd = xyd + derivatives.z;
        double s2d = xyd * -0.21132486540518699998;
        double zzd = derivatives.z * 0.57735026918962599998;
        derivatives.x = derivatives.x + (s2d - zzd);
        derivatives.y = derivatives.y + (s2d - zzd);
        derivatives.z = xyzd * 0.28867513459481294226 + derivatives.w * -0.866025403784439;
        derivatives.w = (xyzd + derivatives.w) * 0.5;

        return value;
    }

    public double noise4_ImproveXYZ_ImproveXZ(double x, double y, double z, double w) {

        double xz = x + z;
        double s2 = xz * -0.21132486540518699998;
        double yy = y * 0.28867513459481294226;
        double ww = w * 0.5;
        double xr = x + (yy + ww + s2), zr = z + (yy + ww + s2);
        double yr = xz * -0.57735026918962599998 + (yy + ww);
        double wr = y * -0.866025403784439 + ww;

        return noise4_UnrotatedBase(xr, yr, zr, wr);
    }

    public double noise4_ImproveXYZ_ImproveXZ(Vector4 derivatives, double x, double y, double z, double w) {

        double xz = x + z;
        double s2 = xz * -0.21132486540518699998;
        double yy = y * 0.28867513459481294226;
        double ww = w * 0.5;
        double xr = x + (yy + ww + s2), zr = z + (yy + ww + s2);
        double yr = xz * -0.57735026918962599998 + (yy + ww);
        double wr = y * -0.866025403784439 + ww;

        double value = noise4_UnrotatedBase(derivatives, xr, yr, zr, wr);

        double xzd = derivatives.x + derivatives.z;
        double xyzd = xzd + derivatives.y;
        double s2d = xzd * -0.21132486540518699998;
        double yyd = derivatives.y * 0.57735026918962599998;
        derivatives.x = derivatives.x + (s2d - yyd);
        derivatives.z = derivatives.z + (s2d - yyd);
        derivatives.y = xyzd * 0.28867513459481294226 + derivatives.w * -0.866025403784439;
        derivatives.w = (xyzd + derivatives.w) * 0.5;

        return value;
    }

    public double noise4_ImproveXYZ(double x, double y, double z, double w) {
        double xyz = x + y + z;
        double s3 = xyz * (-1.0 / 6.0);
        double ww = w * 0.5;
        double xr = x + s3 + ww, yr = y + s3 + ww, zr = z + s3 + ww;
        double wr = xyz * -0.5 + ww;

        return noise4_UnrotatedBase(xr, yr, zr, wr);
    }

    public double noise4_ImproveXYZ(Vector4 derivatives, double x, double y, double z, double w) {
        double xyz = x + y + z;
        double s3 = xyz * (-1.0 / 6.0);
        double ww = w * 0.5;
        double xr = x + s3 + ww, yr = y + s3 + ww, zr = z + s3 + ww;
        double wr = xyz * -0.5 + ww;

        double value = noise4_UnrotatedBase(derivatives, xr, yr, zr, wr);

        double xyzd = derivatives.x + derivatives.y + derivatives.z;
        double s3d = xyzd * (-1.0 / 6.0);
        double wwd = derivatives.w * 0.5;
        derivatives.x = derivatives.x + (s3d - wwd);
        derivatives.y = derivatives.y + (s3d - wwd);
        derivatives.z = derivatives.z + (s3d - wwd);
        derivatives.w = xyzd * 0.5 + wwd;

        return value;
    }

    private double noise4_UnrotatedBase(double xr, double yr, double zr, double wr) {

        int xrb = fastFloor(xr), yrb = fastFloor(yr), zrb = fastFloor(zr), wrb = fastFloor(wr);
        double xri = xr - xrb, yri = yr - yrb, zri = zr - zrb, wri = wr - wrb;

        double g0000 = grad4(4 * (perm[perm[perm[perm[xrb & PMASK] ^ (yrb & PMASK)] ^ (zrb & PMASK)] ^ (wrb & PMASK)] & 0x1F), xri, yri, zri, wri);
        double g0001 = grad4(4 * (perm[perm[perm[perm[xrb & PMASK] ^ (yrb & PMASK)] ^ (zrb & PMASK)] ^ ((wrb + 1) & PMASK)] & 0x1F), xri, yri, zri, wri - 1);
        double g0010 = grad4(4 * (perm[perm[perm[perm[xrb & PMASK] ^ (yrb & PMASK)] ^ ((zrb + 1) & PMASK)] ^ (wrb & PMASK)] & 0x1F), xri, yri, zri - 1, wri);
        double g0011 = grad4(4 * (perm[perm[perm[perm[xrb & PMASK] ^ (yrb & PMASK)] ^ ((zrb + 1) & PMASK)] ^ ((wrb + 1) & PMASK)] & 0x1F), xri, yri, zri - 1, wri - 1);
        double g0100 = grad4(4 * (perm[perm[perm[perm[xrb & PMASK] ^ ((yrb + 1) & PMASK)] ^ (zrb & PMASK)] ^ (wrb & PMASK)] & 0x1F), xri, yri - 1, zri, wri);
        double g0101 = grad4(4 * (perm[perm[perm[perm[xrb & PMASK] ^ ((yrb + 1) & PMASK)] ^ (zrb & PMASK)] ^ ((wrb + 1) & PMASK)] & 0x1F), xri, yri - 1, zri, wri - 1);
        double g0110 = grad4(4 * (perm[perm[perm[perm[xrb & PMASK] ^ ((yrb + 1) & PMASK)] ^ ((zrb + 1) & PMASK)] ^ (wrb & PMASK)] & 0x1F), xri, yri - 1, zri - 1, wri);
        double g0111 = grad4(4 * (perm[perm[perm[perm[xrb & PMASK] ^ ((yrb + 1) & PMASK)] ^ ((zrb + 1) & PMASK)] ^ ((wrb + 1) & PMASK)] & 0x1F), xri, yri - 1, zri - 1, wri - 1);
        double g1000 = grad4(4 * (perm[perm[perm[perm[(xrb + 1) & PMASK] ^ (yrb & PMASK)] ^ (zrb & PMASK)] ^ (wrb & PMASK)] & 0x1F), xri - 1, yri, zri, wri);
        double g1001 = grad4(4 * (perm[perm[perm[perm[(xrb + 1) & PMASK] ^ (yrb & PMASK)] ^ (zrb & PMASK)] ^ ((wrb + 1) & PMASK)] & 0x1F), xri - 1, yri, zri, wri - 1);
        double g1010 = grad4(4 * (perm[perm[perm[perm[(xrb + 1) & PMASK] ^ (yrb & PMASK)] ^ ((zrb + 1) & PMASK)] ^ (wrb & PMASK)] & 0x1F), xri - 1, yri, zri - 1, wri);
        double g1011 = grad4(4 * (perm[perm[perm[perm[(xrb + 1) & PMASK] ^ (yrb & PMASK)] ^ ((zrb + 1) & PMASK)] ^ ((wrb + 1) & PMASK)] & 0x1F), xri - 1, yri, zri - 1, wri - 1);
        double g1100 = grad4(4 * (perm[perm[perm[perm[(xrb + 1) & PMASK] ^ ((yrb + 1) & PMASK)] ^ (zrb & PMASK)] ^ (wrb & PMASK)] & 0x1F), xri - 1, yri - 1, zri, wri);
        double g1101 = grad4(4 * (perm[perm[perm[perm[(xrb + 1) & PMASK] ^ ((yrb + 1) & PMASK)] ^ (zrb & PMASK)] ^ ((wrb + 1) & PMASK)] & 0x1F), xri - 1, yri - 1, zri, wri - 1);
        double g1110 = grad4(4 * (perm[perm[perm[perm[(xrb + 1) & PMASK] ^ ((yrb + 1) & PMASK)] ^ ((zrb + 1) & PMASK)] ^ (wrb & PMASK)] & 0x1F), xri - 1, yri - 1, zri - 1, wri);
        double g1111 = grad4(4 * (perm[perm[perm[perm[(xrb + 1) & PMASK] ^ ((yrb + 1) & PMASK)] ^ ((zrb + 1) & PMASK)] ^ ((wrb + 1) & PMASK)] & 0x1F), xri - 1, yri - 1, zri - 1, wri - 1);

        double fadeX = fadeCurve(xri);
        double fadeY = fadeCurve(yri);
        double fadeZ = fadeCurve(zri);
        double fadeW = fadeCurve(wri);
        double g000W = (1 - fadeW) * g0000 + fadeW * g0001;
        double g001W = (1 - fadeW) * g0010 + fadeW * g0011;
        double g010W = (1 - fadeW) * g0100 + fadeW * g0101;
        double g011W = (1 - fadeW) * g0110 + fadeW * g0111;
        double g100W = (1 - fadeW) * g1000 + fadeW * g1001;
        double g101W = (1 - fadeW) * g1010 + fadeW * g1011;
        double g110W = (1 - fadeW) * g1100 + fadeW * g1101;
        double g111W = (1 - fadeW) * g1110 + fadeW * g1111;
        double g00ZW = (1 - fadeZ) * g000W + fadeZ * g001W;
        double g01ZW = (1 - fadeZ) * g010W + fadeZ * g011W;
        double g10ZW = (1 - fadeZ) * g100W + fadeZ * g101W;
        double g11ZW = (1 - fadeZ) * g110W + fadeZ * g111W;
        double g0YZW = (1 - fadeY) * g00ZW + fadeY * g01ZW;
        double g1YZW = (1 - fadeY) * g10ZW + fadeY * g11ZW;
        double gXYZW = (1 - fadeX) * g0YZW + fadeX * g1YZW;

        return gXYZW;
    }

    private double noise4_UnrotatedBase(Vector4 derivatives, double xr, double yr, double zr, double wr) {

        int xrb = fastFloor(xr), yrb = fastFloor(yr), zrb = fastFloor(zr), wrb = fastFloor(wr);
        double xri = xr - xrb, yri = yr - yrb, zri = zr - zrb, wri = wr - wrb;

        int gi0000 = 4 * (perm[perm[perm[perm[xrb & PMASK] ^ (yrb & PMASK)] ^ (zrb & PMASK)] ^ (wrb & PMASK)] & 0x1F);
        int gi0001 = 4 * (perm[perm[perm[perm[xrb & PMASK] ^ (yrb & PMASK)] ^ (zrb & PMASK)] ^ ((wrb + 1) & PMASK)] & 0x1F);
        int gi0010 = 4 * (perm[perm[perm[perm[xrb & PMASK] ^ (yrb & PMASK)] ^ ((zrb + 1) & PMASK)] ^ (wrb & PMASK)] & 0x1F);
        int gi0011 = 4 * (perm[perm[perm[perm[xrb & PMASK] ^ (yrb & PMASK)] ^ ((zrb + 1) & PMASK)] ^ ((wrb + 1) & PMASK)] & 0x1F);
        int gi0100 = 4 * (perm[perm[perm[perm[xrb & PMASK] ^ ((yrb + 1) & PMASK)] ^ (zrb & PMASK)] ^ (wrb & PMASK)] & 0x1F);
        int gi0101 = 4 * (perm[perm[perm[perm[xrb & PMASK] ^ ((yrb + 1) & PMASK)] ^ (zrb & PMASK)] ^ ((wrb + 1) & PMASK)] & 0x1F);
        int gi0110 = 4 * (perm[perm[perm[perm[xrb & PMASK] ^ ((yrb + 1) & PMASK)] ^ ((zrb + 1) & PMASK)] ^ (wrb & PMASK)] & 0x1F);
        int gi0111 = 4 * (perm[perm[perm[perm[xrb & PMASK] ^ ((yrb + 1) & PMASK)] ^ ((zrb + 1) & PMASK)] ^ ((wrb + 1) & PMASK)] & 0x1F);
        int gi1000 = 4 * (perm[perm[perm[perm[(xrb + 1) & PMASK] ^ (yrb & PMASK)] ^ (zrb & PMASK)] ^ (wrb & PMASK)] & 0x1F);
        int gi1001 = 4 * (perm[perm[perm[perm[(xrb + 1) & PMASK] ^ (yrb & PMASK)] ^ (zrb & PMASK)] ^ ((wrb + 1) & PMASK)] & 0x1F);
        int gi1010 = 4 * (perm[perm[perm[perm[(xrb + 1) & PMASK] ^ (yrb & PMASK)] ^ ((zrb + 1) & PMASK)] ^ (wrb & PMASK)] & 0x1F);
        int gi1011 = 4 * (perm[perm[perm[perm[(xrb + 1) & PMASK] ^ (yrb & PMASK)] ^ ((zrb + 1) & PMASK)] ^ ((wrb + 1) & PMASK)] & 0x1F);
        int gi1100 = 4 * (perm[perm[perm[perm[(xrb + 1) & PMASK] ^ ((yrb + 1) & PMASK)] ^ (zrb & PMASK)] ^ (wrb & PMASK)] & 0x1F);
        int gi1101 = 4 * (perm[perm[perm[perm[(xrb + 1) & PMASK] ^ ((yrb + 1) & PMASK)] ^ (zrb & PMASK)] ^ ((wrb + 1) & PMASK)] & 0x1F);
        int gi1110 = 4 * (perm[perm[perm[perm[(xrb + 1) & PMASK] ^ ((yrb + 1) & PMASK)] ^ ((zrb + 1) & PMASK)] ^ (wrb & PMASK)] & 0x1F);
        int gi1111 = 4 * (perm[perm[perm[perm[(xrb + 1) & PMASK] ^ ((yrb + 1) & PMASK)] ^ ((zrb + 1) & PMASK)] ^ ((wrb + 1) & PMASK)] & 0x1F);

        double g0000 = grad4(gi0000, xri, yri, zri, wri);
        double g0001 = grad4(gi0001, xri, yri, zri, wri - 1);
        double g0010 = grad4(gi0010, xri, yri, zri - 1, wri);
        double g0011 = grad4(gi0011, xri, yri, zri - 1, wri - 1);
        double g0100 = grad4(gi0100, xri, yri - 1, zri, wri);
        double g0101 = grad4(gi0101, xri, yri - 1, zri, wri - 1);
        double g0110 = grad4(gi0110, xri, yri - 1, zri - 1, wri);
        double g0111 = grad4(gi0111, xri, yri - 1, zri - 1, wri - 1);
        double g1000 = grad4(gi1000, xri - 1, yri, zri, wri);
        double g1001 = grad4(gi1001, xri - 1, yri, zri, wri - 1);
        double g1010 = grad4(gi1010, xri - 1, yri, zri - 1, wri);
        double g1011 = grad4(gi1011, xri - 1, yri, zri - 1, wri - 1);
        double g1100 = grad4(gi1100, xri - 1, yri - 1, zri, wri);
        double g1101 = grad4(gi1101, xri - 1, yri - 1, zri, wri - 1);
        double g1110 = grad4(gi1110, xri - 1, yri - 1, zri - 1, wri);
        double g1111 = grad4(gi1111, xri - 1, yri - 1, zri - 1, wri - 1);

        double fadeX = fadeCurve(xri);
        double fadeY = fadeCurve(yri);
        double fadeZ = fadeCurve(zri);
        double fadeW = fadeCurve(wri);
        double g000W = (1 - fadeW) * g0000 + fadeW * g0001;
        double g001W = (1 - fadeW) * g0010 + fadeW * g0011;
        double g010W = (1 - fadeW) * g0100 + fadeW * g0101;
        double g011W = (1 - fadeW) * g0110 + fadeW * g0111;
        double g100W = (1 - fadeW) * g1000 + fadeW * g1001;
        double g101W = (1 - fadeW) * g1010 + fadeW * g1011;
        double g110W = (1 - fadeW) * g1100 + fadeW * g1101;
        double g111W = (1 - fadeW) * g1110 + fadeW * g1111;
        double g00ZW = (1 - fadeZ) * g000W + fadeZ * g001W;
        double g01ZW = (1 - fadeZ) * g010W + fadeZ * g011W;
        double g10ZW = (1 - fadeZ) * g100W + fadeZ * g101W;
        double g11ZW = (1 - fadeZ) * g110W + fadeZ * g111W;
        double g0YZW = (1 - fadeY) * g00ZW + fadeY * g01ZW;
        double g1YZW = (1 - fadeY) * g10ZW + fadeY * g11ZW;
        double gXYZW = (1 - fadeX) * g0YZW + fadeX * g1YZW;

        double dFadeX = dFadeCurve(xri);
        double dFadeY = dFadeCurve(yri);
        double dFadeZ = dFadeCurve(zri);
        double dFadeW = dFadeCurve(wri);
        double dXg000W = (1 - fadeW) * GRADIENTS_4D[gi0000 + 0] + fadeW * GRADIENTS_4D[gi0001 + 0];
        double dYg000W = (1 - fadeW) * GRADIENTS_4D[gi0000 + 1] + fadeW * GRADIENTS_4D[gi0001 + 1];
        double dZg000W = (1 - fadeW) * GRADIENTS_4D[gi0000 + 2] + fadeW * GRADIENTS_4D[gi0001 + 2];
        double dWg000W = (1 - fadeW) * GRADIENTS_4D[gi0000 + 3] + fadeW * GRADIENTS_4D[gi0001 + 3] + (-dFadeW * g0000 + dFadeW * g0001);
        double dXg001W = (1 - fadeW) * GRADIENTS_4D[gi0010 + 0] + fadeW * GRADIENTS_4D[gi0011 + 0];
        double dYg001W = (1 - fadeW) * GRADIENTS_4D[gi0010 + 1] + fadeW * GRADIENTS_4D[gi0011 + 1];
        double dZg001W = (1 - fadeW) * GRADIENTS_4D[gi0010 + 2] + fadeW * GRADIENTS_4D[gi0011 + 2];
        double dWg001W = (1 - fadeW) * GRADIENTS_4D[gi0010 + 3] + fadeW * GRADIENTS_4D[gi0011 + 3] + (-dFadeW * g0010 + dFadeW * g0011);
        double dXg010W = (1 - fadeW) * GRADIENTS_4D[gi0100 + 0] + fadeW * GRADIENTS_4D[gi0101 + 0];
        double dYg010W = (1 - fadeW) * GRADIENTS_4D[gi0100 + 1] + fadeW * GRADIENTS_4D[gi0101 + 1];
        double dZg010W = (1 - fadeW) * GRADIENTS_4D[gi0100 + 2] + fadeW * GRADIENTS_4D[gi0101 + 2];
        double dWg010W = (1 - fadeW) * GRADIENTS_4D[gi0100 + 3] + fadeW * GRADIENTS_4D[gi0101 + 3] + (-dFadeW * g0100 + dFadeW * g0101);
        double dXg011W = (1 - fadeW) * GRADIENTS_4D[gi0110 + 0] + fadeW * GRADIENTS_4D[gi0111 + 0];
        double dYg011W = (1 - fadeW) * GRADIENTS_4D[gi0110 + 1] + fadeW * GRADIENTS_4D[gi0111 + 1];
        double dZg011W = (1 - fadeW) * GRADIENTS_4D[gi0110 + 2] + fadeW * GRADIENTS_4D[gi0111 + 2];
        double dWg011W = (1 - fadeW) * GRADIENTS_4D[gi0110 + 3] + fadeW * GRADIENTS_4D[gi0111 + 3] + (-dFadeW * g0110 + dFadeW * g0111);
        double dXg100W = (1 - fadeW) * GRADIENTS_4D[gi1000 + 0] + fadeW * GRADIENTS_4D[gi1001 + 0];
        double dYg100W = (1 - fadeW) * GRADIENTS_4D[gi1000 + 1] + fadeW * GRADIENTS_4D[gi1001 + 1];
        double dZg100W = (1 - fadeW) * GRADIENTS_4D[gi1000 + 2] + fadeW * GRADIENTS_4D[gi1001 + 2];
        double dWg100W = (1 - fadeW) * GRADIENTS_4D[gi1000 + 3] + fadeW * GRADIENTS_4D[gi1001 + 3] + (-dFadeW * g1000 + dFadeW * g1001);
        double dXg101W = (1 - fadeW) * GRADIENTS_4D[gi1010 + 0] + fadeW * GRADIENTS_4D[gi1011 + 0];
        double dYg101W = (1 - fadeW) * GRADIENTS_4D[gi1010 + 1] + fadeW * GRADIENTS_4D[gi1011 + 1];
        double dZg101W = (1 - fadeW) * GRADIENTS_4D[gi1010 + 2] + fadeW * GRADIENTS_4D[gi1011 + 2];
        double dWg101W = (1 - fadeW) * GRADIENTS_4D[gi1010 + 3] + fadeW * GRADIENTS_4D[gi1011 + 3] + (-dFadeW * g1010 + dFadeW * g1011);
        double dXg110W = (1 - fadeW) * GRADIENTS_4D[gi1100 + 0] + fadeW * GRADIENTS_4D[gi1101 + 0];
        double dYg110W = (1 - fadeW) * GRADIENTS_4D[gi1100 + 1] + fadeW * GRADIENTS_4D[gi1101 + 1];
        double dZg110W = (1 - fadeW) * GRADIENTS_4D[gi1100 + 2] + fadeW * GRADIENTS_4D[gi1101 + 2];
        double dWg110W = (1 - fadeW) * GRADIENTS_4D[gi1100 + 3] + fadeW * GRADIENTS_4D[gi1101 + 3] + (-dFadeW * g1100 + dFadeW * g1101);
        double dXg111W = (1 - fadeW) * GRADIENTS_4D[gi1110 + 0] + fadeW * GRADIENTS_4D[gi1111 + 0];
        double dYg111W = (1 - fadeW) * GRADIENTS_4D[gi1110 + 1] + fadeW * GRADIENTS_4D[gi1111 + 1];
        double dZg111W = (1 - fadeW) * GRADIENTS_4D[gi1110 + 2] + fadeW * GRADIENTS_4D[gi1111 + 2];
        double dWg111W = (1 - fadeW) * GRADIENTS_4D[gi1110 + 3] + fadeW * GRADIENTS_4D[gi1111 + 3] + (-dFadeW * g1110 + dFadeW * g1111);
        double dXg00ZW = (1 - fadeZ) * dXg000W + fadeZ * dXg001W;
        double dYg00ZW = (1 - fadeZ) * dYg000W + fadeZ * dYg001W;
        double dZg00ZW = (1 - fadeZ) * dZg000W + fadeZ * dZg001W + (-dFadeZ * g000W + dFadeZ * g001W);
        double dWg00ZW = (1 - fadeZ) * dWg000W + fadeZ * dWg001W;
        double dXg01ZW = (1 - fadeZ) * dXg010W + fadeZ * dXg011W;
        double dYg01ZW = (1 - fadeZ) * dYg010W + fadeZ * dYg011W;
        double dZg01ZW = (1 - fadeZ) * dZg010W + fadeZ * dZg011W + (-dFadeZ * g010W + dFadeZ * g011W);
        double dWg01ZW = (1 - fadeZ) * dWg010W + fadeZ * dWg011W;
        double dXg10ZW = (1 - fadeZ) * dXg100W + fadeZ * dXg101W;
        double dYg10ZW = (1 - fadeZ) * dYg100W + fadeZ * dYg101W;
        double dZg10ZW = (1 - fadeZ) * dZg100W + fadeZ * dZg101W + (-dFadeZ * g100W + dFadeZ * g101W);
        double dWg10ZW = (1 - fadeZ) * dWg100W + fadeZ * dWg101W;
        double dXg11ZW = (1 - fadeZ) * dXg110W + fadeZ * dXg111W;
        double dYg11ZW = (1 - fadeZ) * dYg110W + fadeZ * dYg111W;
        double dZg11ZW = (1 - fadeZ) * dZg110W + fadeZ * dZg111W + (-dFadeZ * g110W + dFadeZ * g111W);
        double dWg11ZW = (1 - fadeZ) * dWg110W + fadeZ * dWg111W;
        double dXg0YZW = (1 - fadeY) * dXg00ZW + fadeY * dXg01ZW;
        double dYg0YZW = (1 - fadeY) * dYg00ZW + fadeY * dYg01ZW + (-dFadeY * g00ZW + dFadeY * g01ZW);
        double dZg0YZW = (1 - fadeY) * dZg00ZW + fadeY * dZg01ZW;
        double dWg0YZW = (1 - fadeY) * dWg00ZW + fadeY * dWg01ZW;
        double dXg1YZW = (1 - fadeY) * dXg10ZW + fadeY * dXg11ZW;
        double dYg1YZW = (1 - fadeY) * dYg10ZW + fadeY * dYg11ZW + (-dFadeY * g10ZW + dFadeY * g11ZW);
        double dZg1YZW = (1 - fadeY) * dZg10ZW + fadeY * dZg11ZW;
        double dWg1YZW = (1 - fadeY) * dWg10ZW + fadeY * dWg11ZW;
        double dXgXYZW = (1 - fadeX) * dXg0YZW + fadeX * dXg1YZW + (-dFadeX * g0YZW + dFadeX * g1YZW);
        double dYgXYZW = (1 - fadeX) * dYg0YZW + fadeX * dYg1YZW;
        double dZgXYZW = (1 - fadeX) * dZg0YZW + fadeX * dZg1YZW;
        double dWgXYZW = (1 - fadeX) * dWg0YZW + fadeX * dWg1YZW;
        derivatives.x = dXgXYZW;
        derivatives.y = dYgXYZW;
        derivatives.z = dZgXYZW;
        derivatives.w = dWgXYZW;

        return gXYZW;
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

    private static double grad4(int gi, double dx, double dy, double dz, double dw) {
        return GRADIENTS_4D[gi + 0] * dx + GRADIENTS_4D[gi + 1] * dy + GRADIENTS_4D[gi + 2] * dz + GRADIENTS_4D[gi + 3] * dw;
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

    private static final double[] GRADIENTS_3D;
    private static final double[] GRADIENTS_4D = {
             1,  1,  1,  0,
            -1,  1,  1,  0,
             1, -1,  1,  0,
            -1, -1,  1,  0,
             1,  1, -1,  0,
            -1,  1, -1,  0,
             1, -1, -1,  0,
            -1, -1, -1,  0,
             1,  1,  0,  1,
            -1,  1,  0,  1,
             1, -1,  0,  1,
            -1, -1,  0,  1,
             1,  1,  0, -1,
            -1,  1,  0, -1,
             1, -1,  0, -1,
            -1, -1,  0, -1,
             1,  0,  1,  1,
            -1,  0,  1,  1,
             1,  0, -1,  1,
            -1,  0, -1,  1,
             1,  0,  1, -1,
            -1,  0,  1, -1,
             1,  0, -1, -1,
            -1,  0, -1, -1,
             0,  1,  1,  1,
             0, -1,  1,  1,
             0,  1, -1,  1,
             0, -1, -1,  1,
             0,  1,  1, -1,
             0, -1,  1, -1,
             0,  1, -1, -1,
             0, -1, -1, -1
    };

    public static final double NORMALIZING_MULTIPLIER_3D = 1.0 / 1.0363538112118038;
    public static final double NORMALIZING_MULTIPLIER_4D = 1.0 / 1.5365823340468219;
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
            grad3[i] *= NORMALIZING_MULTIPLIER_3D;
        }
        for (int i = 0; i < PSIZE; i++) {
            int iSrc = i % (grad3.length / 4);
            for (int j = 0; j < 3; j++) {
                GRADIENTS_3D[i * 4 + j] = grad3[iSrc * 4 + j];
            }
        }
        for (int i = 0; i < GRADIENTS_4D.length; i++) {
            GRADIENTS_4D[i] *= NORMALIZING_MULTIPLIER_4D;
        }
    }

    /*
     * Inner Classes
     */

    public static class Vector3 {
        public double x, y, z;
    }

    public static class Vector4 {
        public double x, y, z, w;
    }
}
