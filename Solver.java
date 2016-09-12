
class Solver
{
	private static double[] fwdSub(Matrix a, Vector b) {
		int n = b.length();

		double[] x = new double[n];

		for(int i = 0; i < n; i++)
		{
			x[i] = b.get(i);
			
			for(int j = 0; j < i; j++)
			{
				x[i] -= a.get(i, j) * x[j];
			}

			x[i] /= a.get(i, i);
		}

		return x;
	}

	private static double[] backSub(Matrix a, Vector b)	{
		int n = b.length();

		double[] x = new double[n];

		for(int i = n - 1; i >= 0; i--)
		{
			x[i] = b.get(i);
			
			for(int j = i + 1; j < n; j++)
			{
				x[i] -= a.get(i, j) * x[j];
			}

			x[i] /= a.get(i, i);
		}

		return x;
	}

	public static Vector solve_lu_b(Matrix a, Vector b) {
		LUFactorizationResult lu = a.lu_fact();
		return solve_lu_b(lu, b);
	}

	public static Vector solve_lu_b(LUFactorizationResult lu, Vector b) {
		// The idea here is that LUx = b is easier to solve than Ax=y

		// First we forward substitute to solve Ly = b
		double[] y = fwdSub(lu.getL(), b);

		// Then we backward substitute to solve Ux = y
		double[] x = backSub(lu.getU(), new Vector(y));

		return new Vector(x);
	}

	public static Vector solve_qr_b(QRFactorizationResult qr, Vector b) {
		// Ax = b
		// A = QR
		// QRx = b
		// Rx = Q-1 b
		// Q-1 * b = y
		// Rx = y

		// First we multiply q inverse (or transpose, as it's orthonormal) by b to get y
		Vector y = Matrix.multiply(qr.getQ().transpose(), b);

		// Then we backward substitute to get x
		double[] x = backSub(qr.getR(), y);

		return new Vector(x);
	}

	private static Matrix upper(Matrix a){
		int rows = a.getRows();
        int cols = a.getColumns();
        double[][] u = new double[rows][cols];

        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
        		if (i < j) {
                    u[i][j] = a.get(i, j);
                }
            }
        }

        return new Matrix(u);
	}

	private static Matrix lower(Matrix a) {
		int rows = a.getRows();
        int cols = a.getColumns();
        double[][] l = new double[rows][cols];

        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
        		if (i > j) {
                    l[i][j] = a.get(i, j);
                }
            }
        }

        return new Matrix(l);
	}

	private static Matrix diagonal(Matrix a) {
		int rows = a.getRows();
        int cols = a.getColumns();
        double[][] d = new double[rows][cols];

        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
        		if (i == j) {
                    d[i][j] = a.get(i, j);
                }
            }
        }

        return new Matrix(d);
	}

	private static final double MAX_ITERATIONS = 50000;

	public static IterativeSolutionResult jacobi(Matrix a, Vector x0, Vector b, double tolerance) {
        int rows = a.getRows();
        int cols = a.getColumns();

		//Make mlu = -(L + U)
        Matrix mlu = Matrix.scalarMult(-1, Matrix.add(lower(a), upper(a)));

        // Make d = D
		Matrix d = diagonal(a);

        //Iterate
        for (int iter = 0; iter < MAX_ITERATIONS; iter++) {
            //Make -(L + U)xi + b
        	double[][] luxb = Matrix.add(Matrix.multiply(mlu, x0), b).getBacking();

            double[] nextX = new double[x0.length()];
            for (int i = 0; i < rows; i++) {
                nextX[i] = luxb[i][0] / d.get(i, i);
            }

            Vector lastX0 = x0;
            x0 = new Vector(nextX);

            double dx = Vector.add(Vector.scalarMult(-1, lastX0), x0).magnitude();
            if(dx < tolerance)
            {
                return new IterativeSolutionResult(x0, iter);
            }
        }

        return new IterativeSolutionResult(x0, -1);
    }

    public static IterativeSolutionResult gaussSiedel(Matrix a, Vector x0, Vector b, double tolerance) {
        int rows = a.getRows();
        int cols = a.getColumns();
        
        Matrix l = lower(a);
        Matrix u = upper(a);
        Matrix d = diagonal(a);

        //Make (L + D) called ld
        Matrix ld = Matrix.add(l, d);

        //Make -U called negU
        Matrix negU = Matrix.scalarMult(-1, u);

        //Iterate
        for (int iter = 0; iter < MAX_ITERATIONS; iter++) {
            //Make -Uxi + b
            Vector negUxb = Vector.add(Matrix.multiply(negU, x0), b);
            // S(xi1) = -Uxi + b
            Vector lastX0 = x0;
            x0 = new Vector(Solver.fwdSub(ld, negUxb));

            double dx = Vector.add(Vector.scalarMult(-1, lastX0), x0).magnitude();
            if(dx < tolerance)
            {
                return new IterativeSolutionResult(x0, iter);
            }
        }

        return new IterativeSolutionResult(x0, -1);
    }


    public static PowerMethodResult power_method(Matrix a, Vector x0, double tolerance) {
    	double eigenvalue = 0;
    	double[] eigenvector = null;

    	for (int i = 0; i < MAX_ITERATIONS; i++) {
    		Vector lastx0 = x0;

    		x0 = Matrix.multiply(a, x0);
    		eigenvalue = x0.get(0);

	    	eigenvector = new double[x0.length()];

	    	for (int j = 0; j < x0.length(); j++) {
	    		eigenvector[j] = x0.get(j) / eigenvalue;
	    	}

	    	x0 = new Vector(eigenvector);

	    	double dx = Vector.add(Vector.scalarMult(-1, lastx0), x0).magnitude();
            if(dx < tolerance)
            {
            	return new PowerMethodResult(eigenvalue, x0, i);
            }
    	}

    	return new PowerMethodResult(eigenvalue, x0, -1);
    }
}
