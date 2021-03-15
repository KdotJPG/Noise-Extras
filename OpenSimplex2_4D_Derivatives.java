public class OpenSimplex2 {
	
	private static final int PSIZE = 2048;
	private static final int PMASK = 2047;

	private short[] perm;
	private Vector4[] permGrad4;

	public OpenSimplex2_4D_Derivatives(long seed) {
		perm = new short[PSIZE];
		permGrad4 = new Vector4[PSIZE];
		short[] source = new short[PSIZE]; 
		for (short i = 0; i < PSIZE; i++)
			source[i] = i;
		for (int i = PSIZE - 1; i >= 0; i--) {
			seed = seed * 6364136223846793005L + 1442695040888963407L;
			int r = (int)((seed + 31) % (i + 1));
			if (r < 0)
				r += (i + 1);
			perm[i] = source[r];
			permGrad4[i] = GRADIENTS_4D[perm[i]];
			source[r] = source[i];
		}
	}
	
	/*
	 * Noise Evaluators
	 */
	
	/**
	 * 4D OpenSimplex2F noise, classic lattice orientation.
	 */
	public double noise4_Unoriented(double x, double y, double z, double w) {
		
		// Get points for A4 lattice
		double s = -0.138196601125011 * (x + y + z + w);
		double xs = x + s, ys = y + s, zs = z + s, ws = w + s;
		
		return noise4_Base(xs, ys, zs, ws);
	}
    
	public double noise4_Unoriented(Vector4 derivatives, double x, double y, double z, double w) {
		
		// Get points for A4 lattice
		double s = -0.138196601125011 * (x + y + z + w);
		double xs = x + s, ys = y + s, zs = z + s, ws = w + s;
		
		return noise4_Base(derivatives, xs, ys, zs, ws);
	}
	
	/**
	 * 4D OpenSimplex2F noise, with XYZ oriented like noise3_Classic,
	 * and W for an extra degree of freedom. W repeats eventually.
	 * Recommended for time-varied animations which texture a 3D object (W=time)
	 */
	public double noise4_ImproveXYZ(double x, double y, double z, double w) {
		double xyz = x + y + z;
		double ww = w * 0.2236067977499788;
		double s2 = xyz * -0.16666666666666666 + ww;
		double xs = x + s2, ys = y + s2, zs = z + s2, ws = -0.5 * xyz + ww;
		
		return noise4_Base(xs, ys, zs, ws);
	}
    
	public double noise4_ImproveXYZ(Vector4 derivatives, double x, double y, double z, double w) {
		double xyz = x + y + z;
		double ww = w * 0.2236067977499788;
		double s2 = xyz * -0.16666666666666666 + ww;
		double xs = x + s2, ys = y + s2, zs = z + s2, ws = -0.5 * xyz + ww;
		
		double value = noise4_Base(derivatives, xs, ys, zs, ws);
        
        double xyzd = derivatives.x + derivatives.y + derivatives.z;
        double s3d = xyzd * (-1.0 / 6.0);
        double wwd = derivatives.w * 0.5;
        derivatives.x = derivatives.x + (s3d - wwd);
        derivatives.y = derivatives.y + (s3d - wwd);
        derivatives.z = derivatives.z + (s3d - wwd);
        derivatives.w = xyzd * 0.5 + wwd;
        
        return value;
	}
	
	public double noise4_ImproveXYZ_ImproveXY(double x, double y, double z, double w) {
        double xy = x + y;
        double s2 = xy * -0.21132486540518699998;
        double zz = z * 0.28867513459481294226;
        double ww = w * 0.223606797749978;
        double xs = x + (zz + ww + s2), ys = y + (zz + ww + s2);
        double zs = xy * -0.57735026918962599998 + (zz + ww);
        double ws = z * -0.866025403784439 + ww;
		
		return noise4_Base(xs, ys, zs, ws);
	}
	
	public double noise4_ImproveXYZ_ImproveXY(Vector4 derivatives, double x, double y, double z, double w) {
        double xy = x + y;
        double s2 = xy * -0.21132486540518699998;
        double zz = z * 0.28867513459481294226;
        double ww = w * 0.223606797749978;
        double xs = x + (zz + ww + s2), ys = y + (zz + ww + s2);
        double zs = xy * -0.57735026918962599998 + (zz + ww);
        double ws = z * -0.866025403784439 + ww;
		
		double value = noise4_Base(derivatives, xs, ys, zs, ws);

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
        double ww = w * 0.223606797749978;
        double xs = x + (yy + ww + s2), zs = z + (yy + ww + s2);
        double ys = xz * -0.57735026918962599998 + (yy + ww);
        double ws = y * -0.866025403784439 + ww;
		
		return noise4_Base(xs, ys, zs, ws);
	}
	
	public double noise4_ImproveXYZ_ImproveXZ(Vector4 derivatives, double x, double y, double z, double w) {
        double xz = x + z;
        double s2 = xz * -0.21132486540518699998;
        double yy = y * 0.28867513459481294226;
        double ww = w * 0.223606797749978;
        double xs = x + (yy + ww + s2), zs = z + (yy + ww + s2);
        double ys = xz * -0.57735026918962599998 + (yy + ww);
        double ws = y * -0.866025403784439 + ww;
		
		double value = noise4_Base(derivatives, xs, ys, zs, ws);

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
	
	/**
	 * 4D OpenSimplex2 noise base.
	 */
	private double noise4_Base(double xs, double ys, double zs, double ws) {
		double value = 0;
		
		// Get base points and offsets
		int xsb = fastFloor(xs), ysb = fastFloor(ys), zsb = fastFloor(zs), wsb = fastFloor(ws);
		double xsi = xs - xsb, ysi = ys - ysb, zsi = zs - zsb, wsi = ws - wsb;
		
		// If we're in the lower half, flip so we can repeat the code for the upper half. We'll flip back later.
		double siSum = xsi + ysi + zsi + wsi;
		double ssi = siSum * 0.309016994374947; // Prep for vertex contributions.
		boolean inLowerHalf = (siSum < 2);
		if (inLowerHalf) {
			xsi = 1 - xsi; ysi = 1 - ysi; zsi = 1 - zsi; wsi = 1 - wsi;
			siSum = 4 - siSum;
		}
		
		// Consider opposing vertex pairs of the octahedron formed by the central cross-section of the stretched tesseract
		double aabb = xsi + ysi - zsi - wsi, abab = xsi - ysi + zsi - wsi, abba = xsi - ysi - zsi + wsi;
		double aabbScore = Math.abs(aabb), ababScore = Math.abs(abab), abbaScore = Math.abs(abba);
		
		// Find the closest point on the stretched tesseract as if it were the upper half
		int vertexIndex, via, vib;
		double asi, bsi;
		if (aabbScore > ababScore && aabbScore > abbaScore) {
			if (aabb > 0) {
				asi = zsi; bsi = wsi; vertexIndex = 0b0011; via = 0b0111; vib = 0b1011;
			} else {
				asi = xsi; bsi = ysi; vertexIndex = 0b1100; via = 0b1101; vib = 0b1110;
			}
		} else if (ababScore > abbaScore) {
			if (abab > 0) {
				asi = ysi; bsi = wsi; vertexIndex = 0b0101; via = 0b0111; vib = 0b1101;
			} else {
				asi = xsi; bsi = zsi; vertexIndex = 0b1010; via = 0b1011; vib = 0b1110;
			}
		} else {
			if (abba > 0) {
				asi = ysi; bsi = zsi; vertexIndex = 0b1001; via = 0b1011; vib = 0b1101;
			} else {
				asi = xsi; bsi = wsi; vertexIndex = 0b0110; via = 0b0111; vib = 0b1110;
			}
		}
		if (bsi > asi) {
			via = vib;
			double temp = bsi;
			bsi = asi;
			asi = temp;
		}
		if (siSum + asi > 3) {
			vertexIndex = via;
			if (siSum + bsi > 4) {
				vertexIndex = 0b1111;
			}
		}
		
		// Now flip back if we're actually in the lower half.
		if (inLowerHalf) {
			xsi = 1 - xsi; ysi = 1 - ysi; zsi = 1 - zsi; wsi = 1 - wsi;
			vertexIndex ^= 0b1111;
		}
		
		// Five points to add, total, from five copies of the A4 lattice.
		for (int i = 0; i < 5; i++) {
		
			// Update xsb/etc. and add the lattice point's contribution.
			LatticePoint4D c = VERTICES_4D[vertexIndex];
			xsb += c.xsv; ysb += c.ysv; zsb += c.zsv; wsb += c.wsv;
			double xi = xsi + ssi, yi = ysi + ssi, zi = zsi + ssi, wi = wsi + ssi;
			double dx = xi + c.dx, dy = yi + c.dy, dz = zi + c.dz, dw = wi + c.dw;
			double attn = 0.6 - dx * dx - dy * dy - dz * dz - dw * dw;
			if (attn > 0) {
				int pxm = xsb & PMASK, pym = ysb & PMASK, pzm = zsb & PMASK, pwm = wsb & PMASK;
				Vector4 grad = permGrad4[perm[perm[perm[pxm] ^ pym] ^ pzm] ^ pwm];
				double ramped = grad.x * dx + grad.y * dy + grad.z * dz + grad.w * dw;
				
				attn *= attn;
				value += attn * attn * ramped;
			}
			
			// Maybe this helps the compiler/JVM/LLVM/etc. know we can end the loop here. Maybe not.
			if (i == 4) break;
			
			// Update the relative skewed coordinates to reference the vertex we just added.
			// Rather, reference its counterpart on the lattice copy that is shifted down by
			// the vector <-0.2, -0.2, -0.2, -0.2>
			xsi += c.xsi; ysi += c.ysi; zsi += c.zsi; wsi += c.wsi;
			ssi += c.ssiDelta;
			
			// Next point is the closest vertex on the 4-simplex whose base vertex is the aforementioned vertex.
			double score0 = 1.0 + ssi * (-1.0 / 0.309016994374947); // Seems slightly faster than 1.0-xsi-ysi-zsi-wsi
			vertexIndex = 0b0000;
			if (xsi >= ysi && xsi >= zsi && xsi >= wsi && xsi >= score0) {
				vertexIndex = 0b0001;
			}
			else if (ysi > xsi && ysi >= zsi && ysi >= wsi && ysi >= score0) {
				vertexIndex = 0b0010;
			}
			else if (zsi > xsi && zsi > ysi && zsi >= wsi && zsi >= score0) {
				vertexIndex = 0b0100;
			}
			else if (wsi > xsi && wsi > ysi && wsi > zsi && wsi >= score0) {
				vertexIndex = 0b1000;
			}
		}
		
		return value;
	}
	
	/**
	 * 4D OpenSimplex2 noise base.
	 */
	private double noise4_Base(Vector4 derivatives, double xs, double ys, double zs, double ws) {
		double value = 0;
		
		// Get base points and offsets
		int xsb = fastFloor(xs), ysb = fastFloor(ys), zsb = fastFloor(zs), wsb = fastFloor(ws);
		double xsi = xs - xsb, ysi = ys - ysb, zsi = zs - zsb, wsi = ws - wsb;
		
		// If we're in the lower half, flip so we can repeat the code for the upper half. We'll flip back later.
		double siSum = xsi + ysi + zsi + wsi;
		double ssi = siSum * 0.309016994374947; // Prep for vertex contributions.
		boolean inLowerHalf = (siSum < 2);
		if (inLowerHalf) {
			xsi = 1 - xsi; ysi = 1 - ysi; zsi = 1 - zsi; wsi = 1 - wsi;
			siSum = 4 - siSum;
		}
		
		// Consider opposing vertex pairs of the octahedron formed by the central cross-section of the stretched tesseract
		double aabb = xsi + ysi - zsi - wsi, abab = xsi - ysi + zsi - wsi, abba = xsi - ysi - zsi + wsi;
		double aabbScore = Math.abs(aabb), ababScore = Math.abs(abab), abbaScore = Math.abs(abba);
		
		// Find the closest point on the stretched tesseract as if it were the upper half
		int vertexIndex, via, vib;
		double asi, bsi;
		if (aabbScore > ababScore && aabbScore > abbaScore) {
			if (aabb > 0) {
				asi = zsi; bsi = wsi; vertexIndex = 0b0011; via = 0b0111; vib = 0b1011;
			} else {
				asi = xsi; bsi = ysi; vertexIndex = 0b1100; via = 0b1101; vib = 0b1110;
			}
		} else if (ababScore > abbaScore) {
			if (abab > 0) {
				asi = ysi; bsi = wsi; vertexIndex = 0b0101; via = 0b0111; vib = 0b1101;
			} else {
				asi = xsi; bsi = zsi; vertexIndex = 0b1010; via = 0b1011; vib = 0b1110;
			}
		} else {
			if (abba > 0) {
				asi = ysi; bsi = zsi; vertexIndex = 0b1001; via = 0b1011; vib = 0b1101;
			} else {
				asi = xsi; bsi = wsi; vertexIndex = 0b0110; via = 0b0111; vib = 0b1110;
			}
		}
		if (bsi > asi) {
			via = vib;
			double temp = bsi;
			bsi = asi;
			asi = temp;
		}
		if (siSum + asi > 3) {
			vertexIndex = via;
			if (siSum + bsi > 4) {
				vertexIndex = 0b1111;
			}
		}
		
		// Now flip back if we're actually in the lower half.
		if (inLowerHalf) {
			xsi = 1 - xsi; ysi = 1 - ysi; zsi = 1 - zsi; wsi = 1 - wsi;
			vertexIndex ^= 0b1111;
		}
		
		// Five points to add, total, from five copies of the A4 lattice.
		for (int i = 0; i < 5; i++) {
		
			// Update xsb/etc. and add the lattice point's contribution.
			LatticePoint4D c = VERTICES_4D[vertexIndex];
			xsb += c.xsv; ysb += c.ysv; zsb += c.zsv; wsb += c.wsv;
			double xi = xsi + ssi, yi = ysi + ssi, zi = zsi + ssi, wi = wsi + ssi;
			double dx = xi + c.dx, dy = yi + c.dy, dz = zi + c.dz, dw = wi + c.dw;
			double attn = 0.6 - dx * dx - dy * dy - dz * dz - dw * dw;
			if (attn > 0) {
				int pxm = xsb & PMASK, pym = ysb & PMASK, pzm = zsb & PMASK, pwm = wsb & PMASK;
				Vector4 grad = permGrad4[perm[perm[perm[pxm] ^ pym] ^ pzm] ^ pwm];
				double ramped = grad.x * dx + grad.y * dy + grad.z * dz + grad.w * dw;
				
				double attnSq = attn * attn;
				value += attnSq * attnSq * ramped;
                
				derivatives.x += (grad.x * attn - 8 * dx * ramped) * attnSq * attn;
				derivatives.y += (grad.y * attn - 8 * dy * ramped) * attnSq * attn;
				derivatives.z += (grad.z * attn - 8 * dz * ramped) * attnSq * attn;
				derivatives.w += (grad.w * attn - 8 * dw * ramped) * attnSq * attn;
			}
			
			// Maybe this helps the compiler/JVM/LLVM/etc. know we can end the loop here. Maybe not.
			if (i == 4) break;
			
			// Update the relative skewed coordinates to reference the vertex we just added.
			// Rather, reference its counterpart on the lattice copy that is shifted down by
			// the vector <-0.2, -0.2, -0.2, -0.2>
			xsi += c.xsi; ysi += c.ysi; zsi += c.zsi; wsi += c.wsi;
			ssi += c.ssiDelta;
			
			// Next point is the closest vertex on the 4-simplex whose base vertex is the aforementioned vertex.
			double score0 = 1.0 + ssi * (-1.0 / 0.309016994374947); // Seems slightly faster than 1.0-xsi-ysi-zsi-wsi
			vertexIndex = 0b0000;
			if (xsi >= ysi && xsi >= zsi && xsi >= wsi && xsi >= score0) {
				vertexIndex = 0b0001;
			}
			else if (ysi > xsi && ysi >= zsi && ysi >= wsi && ysi >= score0) {
				vertexIndex = 0b0010;
			}
			else if (zsi > xsi && zsi > ysi && zsi >= wsi && zsi >= score0) {
				vertexIndex = 0b0100;
			}
			else if (wsi > xsi && wsi > ysi && wsi > zsi && wsi >= score0) {
				vertexIndex = 0b1000;
			}
		}
		
		return value;
	}
	
	/*
	 * Utility
	 */
	
	private static int fastFloor(double x) {
		int xi = (int)x;
		return x < xi ? xi - 1 : xi;
	}
	
	/*
	 * Definitions
	 */

	private static final LatticePoint4D[] VERTICES_4D;
	static {
		VERTICES_4D = new LatticePoint4D[16];
		
		for (int i = 0; i < 16; i++) {
			VERTICES_4D[i] = new LatticePoint4D((i >> 0) & 1, (i >> 1) & 1, (i >> 2) & 1, (i >> 3) & 1);
		}
	}
	
	private static class LatticePoint4D {
		int xsv, ysv, zsv, wsv;
		double dx, dy, dz, dw;
		double xsi, ysi, zsi, wsi;
		double ssiDelta;
		public LatticePoint4D(int xsv, int ysv, int zsv, int wsv) {
			this.xsv = xsv + 409; this.ysv = ysv + 409; this.zsv = zsv + 409; this.wsv = wsv + 409;
			double ssv = (xsv + ysv + zsv + wsv) * 0.309016994374947;
			this.dx = -xsv - ssv;
			this.dy = -ysv - ssv;
			this.dz = -zsv - ssv;
			this.dw = -wsv - ssv;
			this.xsi = xsi = 0.2 - xsv;
			this.ysi = ysi = 0.2 - ysv;
			this.zsi = zsi = 0.2 - zsv;
			this.wsi = wsi = 0.2 - wsv;
			this.ssiDelta = (0.8 - xsv - ysv - zsv - wsv) * 0.309016994374947;
		}
	}
	
	/*
	 * Gradients
	 */
	
	public static class Vector4 {
		double x, y, z, w;
		public Vector4() { }
		public Vector4(double x, double y, double z,  double w) {
			this.x = x; this.y = y; this.z = z; this.w = w;
		}
	}
	
	private static final double N4 = 0.020891026027489927;
	private static final Vector4[] GRADIENTS_4D;
	static {
		GRADIENTS_4D = new Vector4[PSIZE];
		Vector4[] grad4 = {
			new Vector4(-0.753341017856078,    -0.37968289875261624,  -0.37968289875261624,  -0.37968289875261624),
			new Vector4(-0.7821684431180708,   -0.4321472685365301,   -0.4321472685365301,    0.12128480194602098),
			new Vector4(-0.7821684431180708,   -0.4321472685365301,    0.12128480194602098,  -0.4321472685365301),
			new Vector4(-0.7821684431180708,    0.12128480194602098,  -0.4321472685365301,   -0.4321472685365301),
			new Vector4(-0.8586508742123365,   -0.508629699630796,     0.044802370851755174,  0.044802370851755174),
			new Vector4(-0.8586508742123365,    0.044802370851755174, -0.508629699630796,     0.044802370851755174),
			new Vector4(-0.8586508742123365,    0.044802370851755174,  0.044802370851755174, -0.508629699630796),
			new Vector4(-0.9982828964265062,   -0.03381941603233842,  -0.03381941603233842,  -0.03381941603233842),
			new Vector4(-0.37968289875261624,  -0.753341017856078,    -0.37968289875261624,  -0.37968289875261624),
			new Vector4(-0.4321472685365301,   -0.7821684431180708,   -0.4321472685365301,    0.12128480194602098),
			new Vector4(-0.4321472685365301,   -0.7821684431180708,    0.12128480194602098,  -0.4321472685365301),
			new Vector4( 0.12128480194602098,  -0.7821684431180708,   -0.4321472685365301,   -0.4321472685365301),
			new Vector4(-0.508629699630796,    -0.8586508742123365,    0.044802370851755174,  0.044802370851755174),
			new Vector4( 0.044802370851755174, -0.8586508742123365,   -0.508629699630796,     0.044802370851755174),
			new Vector4( 0.044802370851755174, -0.8586508742123365,    0.044802370851755174, -0.508629699630796),
			new Vector4(-0.03381941603233842,  -0.9982828964265062,   -0.03381941603233842,  -0.03381941603233842),
			new Vector4(-0.37968289875261624,  -0.37968289875261624,  -0.753341017856078,    -0.37968289875261624),
			new Vector4(-0.4321472685365301,   -0.4321472685365301,   -0.7821684431180708,    0.12128480194602098),
			new Vector4(-0.4321472685365301,    0.12128480194602098,  -0.7821684431180708,   -0.4321472685365301),
			new Vector4( 0.12128480194602098,  -0.4321472685365301,   -0.7821684431180708,   -0.4321472685365301),
			new Vector4(-0.508629699630796,     0.044802370851755174, -0.8586508742123365,    0.044802370851755174),
			new Vector4( 0.044802370851755174, -0.508629699630796,    -0.8586508742123365,    0.044802370851755174),
			new Vector4( 0.044802370851755174,  0.044802370851755174, -0.8586508742123365,   -0.508629699630796),
			new Vector4(-0.03381941603233842,  -0.03381941603233842,  -0.9982828964265062,   -0.03381941603233842),
			new Vector4(-0.37968289875261624,  -0.37968289875261624,  -0.37968289875261624,  -0.753341017856078),
			new Vector4(-0.4321472685365301,   -0.4321472685365301,    0.12128480194602098,  -0.7821684431180708),
			new Vector4(-0.4321472685365301,    0.12128480194602098,  -0.4321472685365301,   -0.7821684431180708),
			new Vector4( 0.12128480194602098,  -0.4321472685365301,   -0.4321472685365301,   -0.7821684431180708),
			new Vector4(-0.508629699630796,     0.044802370851755174,  0.044802370851755174, -0.8586508742123365),
			new Vector4( 0.044802370851755174, -0.508629699630796,     0.044802370851755174, -0.8586508742123365),
			new Vector4( 0.044802370851755174,  0.044802370851755174, -0.508629699630796,    -0.8586508742123365),
			new Vector4(-0.03381941603233842,  -0.03381941603233842,  -0.03381941603233842,  -0.9982828964265062),
			new Vector4(-0.6740059517812944,   -0.3239847771997537,   -0.3239847771997537,    0.5794684678643381),
			new Vector4(-0.7504883828755602,   -0.4004672082940195,    0.15296486218853164,   0.5029860367700724),
			new Vector4(-0.7504883828755602,    0.15296486218853164,  -0.4004672082940195,    0.5029860367700724),
			new Vector4(-0.8828161875373585,    0.08164729285680945,   0.08164729285680945,   0.4553054119602712),
			new Vector4(-0.4553054119602712,   -0.08164729285680945,  -0.08164729285680945,   0.8828161875373585),
			new Vector4(-0.5029860367700724,   -0.15296486218853164,   0.4004672082940195,    0.7504883828755602),
			new Vector4(-0.5029860367700724,    0.4004672082940195,   -0.15296486218853164,   0.7504883828755602),
			new Vector4(-0.5794684678643381,    0.3239847771997537,    0.3239847771997537,    0.6740059517812944),
			new Vector4(-0.3239847771997537,   -0.6740059517812944,   -0.3239847771997537,    0.5794684678643381),
			new Vector4(-0.4004672082940195,   -0.7504883828755602,    0.15296486218853164,   0.5029860367700724),
			new Vector4( 0.15296486218853164,  -0.7504883828755602,   -0.4004672082940195,    0.5029860367700724),
			new Vector4( 0.08164729285680945,  -0.8828161875373585,    0.08164729285680945,   0.4553054119602712),
			new Vector4(-0.08164729285680945,  -0.4553054119602712,   -0.08164729285680945,   0.8828161875373585),
			new Vector4(-0.15296486218853164,  -0.5029860367700724,    0.4004672082940195,    0.7504883828755602),
			new Vector4( 0.4004672082940195,   -0.5029860367700724,   -0.15296486218853164,   0.7504883828755602),
			new Vector4( 0.3239847771997537,   -0.5794684678643381,    0.3239847771997537,    0.6740059517812944),
			new Vector4(-0.3239847771997537,   -0.3239847771997537,   -0.6740059517812944,    0.5794684678643381),
			new Vector4(-0.4004672082940195,    0.15296486218853164,  -0.7504883828755602,    0.5029860367700724),
			new Vector4( 0.15296486218853164,  -0.4004672082940195,   -0.7504883828755602,    0.5029860367700724),
			new Vector4( 0.08164729285680945,   0.08164729285680945,  -0.8828161875373585,    0.4553054119602712),
			new Vector4(-0.08164729285680945,  -0.08164729285680945,  -0.4553054119602712,    0.8828161875373585),
			new Vector4(-0.15296486218853164,   0.4004672082940195,   -0.5029860367700724,    0.7504883828755602),
			new Vector4( 0.4004672082940195,   -0.15296486218853164,  -0.5029860367700724,    0.7504883828755602),
			new Vector4( 0.3239847771997537,    0.3239847771997537,   -0.5794684678643381,    0.6740059517812944),
			new Vector4(-0.6740059517812944,   -0.3239847771997537,    0.5794684678643381,   -0.3239847771997537),
			new Vector4(-0.7504883828755602,   -0.4004672082940195,    0.5029860367700724,    0.15296486218853164),
			new Vector4(-0.7504883828755602,    0.15296486218853164,   0.5029860367700724,   -0.4004672082940195),
			new Vector4(-0.8828161875373585,    0.08164729285680945,   0.4553054119602712,    0.08164729285680945),
			new Vector4(-0.4553054119602712,   -0.08164729285680945,   0.8828161875373585,   -0.08164729285680945),
			new Vector4(-0.5029860367700724,   -0.15296486218853164,   0.7504883828755602,    0.4004672082940195),
			new Vector4(-0.5029860367700724,    0.4004672082940195,    0.7504883828755602,   -0.15296486218853164),
			new Vector4(-0.5794684678643381,    0.3239847771997537,    0.6740059517812944,    0.3239847771997537),
			new Vector4(-0.3239847771997537,   -0.6740059517812944,    0.5794684678643381,   -0.3239847771997537),
			new Vector4(-0.4004672082940195,   -0.7504883828755602,    0.5029860367700724,    0.15296486218853164),
			new Vector4( 0.15296486218853164,  -0.7504883828755602,    0.5029860367700724,   -0.4004672082940195),
			new Vector4( 0.08164729285680945,  -0.8828161875373585,    0.4553054119602712,    0.08164729285680945),
			new Vector4(-0.08164729285680945,  -0.4553054119602712,    0.8828161875373585,   -0.08164729285680945),
			new Vector4(-0.15296486218853164,  -0.5029860367700724,    0.7504883828755602,    0.4004672082940195),
			new Vector4( 0.4004672082940195,   -0.5029860367700724,    0.7504883828755602,   -0.15296486218853164),
			new Vector4( 0.3239847771997537,   -0.5794684678643381,    0.6740059517812944,    0.3239847771997537),
			new Vector4(-0.3239847771997537,   -0.3239847771997537,    0.5794684678643381,   -0.6740059517812944),
			new Vector4(-0.4004672082940195,    0.15296486218853164,   0.5029860367700724,   -0.7504883828755602),
			new Vector4( 0.15296486218853164,  -0.4004672082940195,    0.5029860367700724,   -0.7504883828755602),
			new Vector4( 0.08164729285680945,   0.08164729285680945,   0.4553054119602712,   -0.8828161875373585),
			new Vector4(-0.08164729285680945,  -0.08164729285680945,   0.8828161875373585,   -0.4553054119602712),
			new Vector4(-0.15296486218853164,   0.4004672082940195,    0.7504883828755602,   -0.5029860367700724),
			new Vector4( 0.4004672082940195,   -0.15296486218853164,   0.7504883828755602,   -0.5029860367700724),
			new Vector4( 0.3239847771997537,    0.3239847771997537,    0.6740059517812944,   -0.5794684678643381),
			new Vector4(-0.6740059517812944,    0.5794684678643381,   -0.3239847771997537,   -0.3239847771997537),
			new Vector4(-0.7504883828755602,    0.5029860367700724,   -0.4004672082940195,    0.15296486218853164),
			new Vector4(-0.7504883828755602,    0.5029860367700724,    0.15296486218853164,  -0.4004672082940195),
			new Vector4(-0.8828161875373585,    0.4553054119602712,    0.08164729285680945,   0.08164729285680945),
			new Vector4(-0.4553054119602712,    0.8828161875373585,   -0.08164729285680945,  -0.08164729285680945),
			new Vector4(-0.5029860367700724,    0.7504883828755602,   -0.15296486218853164,   0.4004672082940195),
			new Vector4(-0.5029860367700724,    0.7504883828755602,    0.4004672082940195,   -0.15296486218853164),
			new Vector4(-0.5794684678643381,    0.6740059517812944,    0.3239847771997537,    0.3239847771997537),
			new Vector4(-0.3239847771997537,    0.5794684678643381,   -0.6740059517812944,   -0.3239847771997537),
			new Vector4(-0.4004672082940195,    0.5029860367700724,   -0.7504883828755602,    0.15296486218853164),
			new Vector4( 0.15296486218853164,   0.5029860367700724,   -0.7504883828755602,   -0.4004672082940195),
			new Vector4( 0.08164729285680945,   0.4553054119602712,   -0.8828161875373585,    0.08164729285680945),
			new Vector4(-0.08164729285680945,   0.8828161875373585,   -0.4553054119602712,   -0.08164729285680945),
			new Vector4(-0.15296486218853164,   0.7504883828755602,   -0.5029860367700724,    0.4004672082940195),
			new Vector4( 0.4004672082940195,    0.7504883828755602,   -0.5029860367700724,   -0.15296486218853164),
			new Vector4( 0.3239847771997537,    0.6740059517812944,   -0.5794684678643381,    0.3239847771997537),
			new Vector4(-0.3239847771997537,    0.5794684678643381,   -0.3239847771997537,   -0.6740059517812944),
			new Vector4(-0.4004672082940195,    0.5029860367700724,    0.15296486218853164,  -0.7504883828755602),
			new Vector4( 0.15296486218853164,   0.5029860367700724,   -0.4004672082940195,   -0.7504883828755602),
			new Vector4( 0.08164729285680945,   0.4553054119602712,    0.08164729285680945,  -0.8828161875373585),
			new Vector4(-0.08164729285680945,   0.8828161875373585,   -0.08164729285680945,  -0.4553054119602712),
			new Vector4(-0.15296486218853164,   0.7504883828755602,    0.4004672082940195,   -0.5029860367700724),
			new Vector4( 0.4004672082940195,    0.7504883828755602,   -0.15296486218853164,  -0.5029860367700724),
			new Vector4( 0.3239847771997537,    0.6740059517812944,    0.3239847771997537,   -0.5794684678643381),
			new Vector4( 0.5794684678643381,   -0.6740059517812944,   -0.3239847771997537,   -0.3239847771997537),
			new Vector4( 0.5029860367700724,   -0.7504883828755602,   -0.4004672082940195,    0.15296486218853164),
			new Vector4( 0.5029860367700724,   -0.7504883828755602,    0.15296486218853164,  -0.4004672082940195),
			new Vector4( 0.4553054119602712,   -0.8828161875373585,    0.08164729285680945,   0.08164729285680945),
			new Vector4( 0.8828161875373585,   -0.4553054119602712,   -0.08164729285680945,  -0.08164729285680945),
			new Vector4( 0.7504883828755602,   -0.5029860367700724,   -0.15296486218853164,   0.4004672082940195),
			new Vector4( 0.7504883828755602,   -0.5029860367700724,    0.4004672082940195,   -0.15296486218853164),
			new Vector4( 0.6740059517812944,   -0.5794684678643381,    0.3239847771997537,    0.3239847771997537),
			new Vector4( 0.5794684678643381,   -0.3239847771997537,   -0.6740059517812944,   -0.3239847771997537),
			new Vector4( 0.5029860367700724,   -0.4004672082940195,   -0.7504883828755602,    0.15296486218853164),
			new Vector4( 0.5029860367700724,    0.15296486218853164,  -0.7504883828755602,   -0.4004672082940195),
			new Vector4( 0.4553054119602712,    0.08164729285680945,  -0.8828161875373585,    0.08164729285680945),
			new Vector4( 0.8828161875373585,   -0.08164729285680945,  -0.4553054119602712,   -0.08164729285680945),
			new Vector4( 0.7504883828755602,   -0.15296486218853164,  -0.5029860367700724,    0.4004672082940195),
			new Vector4( 0.7504883828755602,    0.4004672082940195,   -0.5029860367700724,   -0.15296486218853164),
			new Vector4( 0.6740059517812944,    0.3239847771997537,   -0.5794684678643381,    0.3239847771997537),
			new Vector4( 0.5794684678643381,   -0.3239847771997537,   -0.3239847771997537,   -0.6740059517812944),
			new Vector4( 0.5029860367700724,   -0.4004672082940195,    0.15296486218853164,  -0.7504883828755602),
			new Vector4( 0.5029860367700724,    0.15296486218853164,  -0.4004672082940195,   -0.7504883828755602),
			new Vector4( 0.4553054119602712,    0.08164729285680945,   0.08164729285680945,  -0.8828161875373585),
			new Vector4( 0.8828161875373585,   -0.08164729285680945,  -0.08164729285680945,  -0.4553054119602712),
			new Vector4( 0.7504883828755602,   -0.15296486218853164,   0.4004672082940195,   -0.5029860367700724),
			new Vector4( 0.7504883828755602,    0.4004672082940195,   -0.15296486218853164,  -0.5029860367700724),
			new Vector4( 0.6740059517812944,    0.3239847771997537,    0.3239847771997537,   -0.5794684678643381),
			new Vector4( 0.03381941603233842,   0.03381941603233842,   0.03381941603233842,   0.9982828964265062),
			new Vector4(-0.044802370851755174, -0.044802370851755174,  0.508629699630796,     0.8586508742123365),
			new Vector4(-0.044802370851755174,  0.508629699630796,    -0.044802370851755174,  0.8586508742123365),
			new Vector4(-0.12128480194602098,   0.4321472685365301,    0.4321472685365301,    0.7821684431180708),
			new Vector4( 0.508629699630796,    -0.044802370851755174, -0.044802370851755174,  0.8586508742123365),
			new Vector4( 0.4321472685365301,   -0.12128480194602098,   0.4321472685365301,    0.7821684431180708),
			new Vector4( 0.4321472685365301,    0.4321472685365301,   -0.12128480194602098,   0.7821684431180708),
			new Vector4( 0.37968289875261624,   0.37968289875261624,   0.37968289875261624,   0.753341017856078),
			new Vector4( 0.03381941603233842,   0.03381941603233842,   0.9982828964265062,    0.03381941603233842),
			new Vector4(-0.044802370851755174,  0.044802370851755174,  0.8586508742123365,    0.508629699630796),
			new Vector4(-0.044802370851755174,  0.508629699630796,     0.8586508742123365,   -0.044802370851755174),
			new Vector4(-0.12128480194602098,   0.4321472685365301,    0.7821684431180708,    0.4321472685365301),
			new Vector4( 0.508629699630796,    -0.044802370851755174,  0.8586508742123365,   -0.044802370851755174),
			new Vector4( 0.4321472685365301,   -0.12128480194602098,   0.7821684431180708,    0.4321472685365301),
			new Vector4( 0.4321472685365301,    0.4321472685365301,    0.7821684431180708,   -0.12128480194602098),
			new Vector4( 0.37968289875261624,   0.37968289875261624,   0.753341017856078,     0.37968289875261624),
			new Vector4( 0.03381941603233842,   0.9982828964265062,    0.03381941603233842,   0.03381941603233842),
			new Vector4(-0.044802370851755174,  0.8586508742123365,   -0.044802370851755174,  0.508629699630796),
			new Vector4(-0.044802370851755174,  0.8586508742123365,    0.508629699630796,    -0.044802370851755174),
			new Vector4(-0.12128480194602098,   0.7821684431180708,    0.4321472685365301,    0.4321472685365301),
			new Vector4( 0.508629699630796,     0.8586508742123365,   -0.044802370851755174, -0.044802370851755174),
			new Vector4( 0.4321472685365301,    0.7821684431180708,   -0.12128480194602098,   0.4321472685365301),
			new Vector4( 0.4321472685365301,    0.7821684431180708,    0.4321472685365301,   -0.12128480194602098),
			new Vector4( 0.37968289875261624,   0.753341017856078,     0.37968289875261624,   0.37968289875261624),
			new Vector4( 0.9982828964265062,    0.03381941603233842,   0.03381941603233842,   0.03381941603233842),
			new Vector4( 0.8586508742123365,   -0.044802370851755174, -0.044802370851755174,  0.508629699630796),
			new Vector4( 0.8586508742123365,   -0.044802370851755174,  0.508629699630796,    -0.044802370851755174),
			new Vector4( 0.7821684431180708,   -0.12128480194602098,   0.4321472685365301,    0.4321472685365301),
			new Vector4( 0.8586508742123365,    0.508629699630796,    -0.044802370851755174, -0.044802370851755174),
			new Vector4( 0.7821684431180708,    0.4321472685365301,   -0.12128480194602098,   0.4321472685365301),
			new Vector4( 0.7821684431180708,    0.4321472685365301,    0.4321472685365301,   -0.12128480194602098),
			new Vector4( 0.753341017856078,     0.37968289875261624,   0.37968289875261624,   0.37968289875261624)
		};
		for (int i = 0; i < grad4.length; i++) {
			grad4[i].x /= N4; grad4[i].y /= N4; grad4[i].z /= N4; grad4[i].w /= N4;
		}
		for (int i = 0; i < PSIZE; i++) {
			GRADIENTS_4D[i] = grad4[i % grad4.length];
		}
	}
}
