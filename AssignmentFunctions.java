import java.io.PrintStream;

public class AssignmentFunctions
{
	public static Vector hilbertB(int size)
	{
		double[] vals = new double[size];

		double b = Math.pow(0.1, size / 3.0);

		for(int i = 0; i < size; i++)
		{
			vals[i] = b;
		}

		return new Vector(vals);
	}

	public static void hilbertMatrixFunction()
	{
		PrintStream errFile = null;

		try
		{
			 errFile = new PrintStream("hilbertErrors.csv");		
		} catch (Exception exc) { }

			errFile.println("n\tLU\tHouseholder\tGivens\tSolution LU\tSolution QR Householder\tSolution QR Givens");

			for(int n = 3; n <= 75; n++)
			{
				System.out.println("n = " + n);
				Matrix a = Matrix.hilbertMatrix(n);
				Vector b = hilbertB(n);

				LUFactorizationResult lu = a.lu_fact();
				QRFactorizationResult qrH = a.qr_fact_househ();
				QRFactorizationResult qrG = a.qr_fact_givens();

				System.out.println("LU Error: " + lu.getError());
				System.out.println("QR Householder Error: " + qrH.getError());
				System.out.println("QR Givens Error: " + qrG.getError());

				// Compute a solution to the system
				Vector xlu = Solver.solve_lu_b(lu, b);
				Vector xqrh = Solver.solve_qr_b(qrH, b);
				Vector xqrg = Solver.solve_qr_b(qrG, b);

				System.out.println("Xsol:\n" + xlu);

				// Calculate b with Ax = b
				Vector calcBlu = Matrix.multiply(a, xlu);
				// ||Ax - B||inf
				double solErrorLu = Vector.add(Vector.scalarMult(-1, calcBlu), b).magnitude();

				// Calculate b with Ax = b
				Vector calcBqrh = Matrix.multiply(a, xqrh);
				// ||Ax - B||inf
				double solErrorQrh = Vector.add(Vector.scalarMult(-1, calcBqrh), b).magnitude();
				
				// Calculate b with Ax = b
				Vector calcBqrg = Matrix.multiply(a, xqrg);
				// ||Ax - B||inf
				double solErrorQrg = Vector.add(Vector.scalarMult(-1, calcBqrg), b).magnitude();

				// Write to the csv
				errFile.println(n + "\t" + lu.getError() + "\t" + qrH.getError() + "\t" + qrG.getError()  + "\t" + solErrorLu + "\t" + solErrorQrh + "\t" + solErrorQrg);

				System.out.println();
				System.out.println();
			}

			errFile.close();
	}

	public static void convolutionalEncodingFunction()
	{
		PrintStream iterFile = null;

		try
		{
			 iterFile = new PrintStream("convIter.csv");		
		} catch (Exception exc) { }

		iterFile.println("n\tj0\tj1\tg0\tg1\tej0\tej1\teg0\teg1");

		double tolerance = 0.00000001;

		// N is vector size
		for (int n = 1; n < 200; n++) {
			// We do some trials for every n
			for(int i = 0; i < 10000; i++)
			{
				Vector x = Vector.randomVector(n);

				Vector encoded = Convolution.encode(x);

				ConvolutionDecodeResult j0 = Convolution.decodeY0j(encoded, tolerance);
				ConvolutionDecodeResult j1 = Convolution.decodeY1j(encoded, tolerance);
				ConvolutionDecodeResult g0 = Convolution.decodeY0g(encoded, tolerance);
				ConvolutionDecodeResult g1 = Convolution.decodeY1g(encoded, tolerance);

				double ej0 = Vector.add(Vector.scalarMult(-1, j0.getV()), x).magnitude();
				double ej1 = Vector.add(Vector.scalarMult(-1, j1.getV()), x).magnitude();
				double eg0 = Vector.add(Vector.scalarMult(-1, g0.getV()), x).magnitude();
				double eg1 = Vector.add(Vector.scalarMult(-1, g1.getV()), x).magnitude();

				iterFile.println(n + "\t" + j0.getIterations() + "\t" + j1.getIterations() + "\t" + g0.getIterations() + "\t" + g1.getIterations() + "\t" + ej0 + "\t" + ej1 + "\t" + eg0 + "\t" + eg1);

				System.out.println("n: " + n + " i: " + i);
			}
		}

		iterFile.close();
	}

	public static void leslieFunction()
	{
		Vector x = new Vector(new double[] { 5.19, 4.45, 1.25, 1.45, 1.46, 1.66, 1.41, 1.05, 0.37 });

		double[][] ad = new double[9][9];
		ad[0][1] = 0.6;
		ad[0][2] = 1.1;
		ad[0][3] = 0.9;
		ad[0][4] = 0.1;

		ad[1][0] = 0.7;
		ad[2][1] = 0.85;
		ad[3][2] = 0.9;
		ad[4][3] = 0.9;
		ad[5][4] = 0.88;
		ad[6][5] = 0.8;
		ad[7][6] = 0.77;
		ad[8][7] = 0.4;

		Matrix a = new Matrix(ad);

		for (int n = 0; n < 6; n++) {
			System.out.print(2020 + 10 * n);

			for (double[] d : x.getBacking()) {
				System.out.print("\t" + d[0]);
			}

			System.out.println();
			x = Matrix.multiply(a, x);
		}

		PowerMethodResult res = Solver.power_method(a, x, 0.00000001);

		System.out.println(res.getEigenvalue());
	}
}
