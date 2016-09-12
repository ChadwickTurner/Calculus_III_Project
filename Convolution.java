
class Convolution {

    public static Matrix makeY0Matrix(int rows) {
        double[][] x = new double[rows][rows];

        for (int i = 0; i < rows; i++) {
            x[i][i] = 1;
            if (i + 2 < rows) {
                x[i + 2][i] = 1;
            }
            if (i + 3 < rows) {
                x[i + 3][i] = 1;
            }
        }

        return new Matrix(x);
    }

    public static Matrix makeY1Matrix(int rows) {
        double[][] x = new double[rows][rows];

        for (int i = 0; i < rows; i++) {
            x[i][i] = 1;
            if (i + 1 < rows) {
                x[i + 1][i] = 1;
            }
            if (i + 3 < rows) {
                x[i + 3][i] = 1;
            }
        }

        return new Matrix(x);
    }

    private static Vector binaryize(Vector y) {
        double[] stuff = new double[y.length()];
        for (int i = 0; i < stuff.length; i++) {
            stuff[i] = y.get(i) % 2;
        }
        return new Vector(stuff);
    }

    public static Vector encode(Vector data) {
        int dataSize = data.length() + 3;

        double[] xd = new double[dataSize];
        double[][] dataBack = data.getBacking();

        for (int i = 0; i < dataSize - 3; i++) {
            xd[i] = dataBack[i][0];
        }

        Vector v = new Vector(xd);

        Matrix a0 = makeY0Matrix(dataSize);
        Matrix a1 = makeY1Matrix(dataSize);


        Vector y0 = binaryize(Matrix.multiply(a0, v));
        Vector y1 = binaryize(Matrix.multiply(a1, v));

        double[] output = new double[dataSize * 2];

        for (int i = 0; i < dataSize; i++) {
            output[2 * i] = y0.get(i);
            output[2 * i + 1] = y1.get(i);
        }

        return new Vector(output);
    }

    public static ConvolutionDecodeResult decodeY0j(Vector v, double tolerance) {
        Matrix a0 = makeY0Matrix(v.length() / 2);

        double[] y0d = new double[v.length() / 2];

        for (int i = 0; i < y0d.length; i++) {
            y0d[i] = v.get(i * 2);
        }

        Vector y0 = new Vector(y0d);

        double[] x0d = new double[y0d.length];

        for (int i = 0; i < x0d.length; i++) {
            x0d[i] = 0;
        }

        Vector x0 = new Vector(x0d);

        IterativeSolutionResult res = Solver.jacobi(a0, x0, y0, tolerance);
        Vector x = res.getV();
        int iterations = res.getIterations();
        double[] result = new double[x.length() - 3];

        for(int i = 0; i < result.length; i++)
        {
            result[i] = Math.abs(x.get(i)) % 2;
        }

        return new ConvolutionDecodeResult(new Vector(result), res.getIterations());
    }

    public static ConvolutionDecodeResult decodeY1j(Vector v, double tolerance) {
        Matrix a1 = makeY1Matrix(v.length() / 2);

        double[] y1d = new double[v.length() / 2];

        for (int i = 0; i < y1d.length; i++) {
            y1d[i] = v.get(i * 2 + 1);
        }

        Vector y1 = new Vector(y1d);

        double[] x0d = new double[y1d.length];

        for (int i = 0; i < x0d.length; i++) {
            x0d[i] = 0;
        }

        Vector x0 = new Vector(x0d);

        IterativeSolutionResult res = Solver.jacobi(a1, x0, y1, tolerance);
        Vector x = res.getV();
        int iterations = res.getIterations();

        double[] result = new double[x.length() - 3];

        for(int i = 0; i < result.length; i++)
        {
            result[i] = Math.abs(x.get(i)) % 2;
        }

        return new ConvolutionDecodeResult(new Vector(result), res.getIterations());
    }

    public static ConvolutionDecodeResult decodeY0g(Vector v, double tolerance) {
        Matrix a0 = makeY0Matrix(v.length() / 2);

        double[] y0d = new double[v.length() / 2];

        for (int i = 0; i < y0d.length; i++) {
            y0d[i] = v.get(i * 2);
        }

        Vector y0 = new Vector(y0d);

        double[] x0d = new double[y0d.length];

        for (int i = 0; i < x0d.length; i++) {
            x0d[i] = 0;
        }

        Vector x0 = new Vector(x0d);

        IterativeSolutionResult res = Solver.gaussSiedel(a0, x0, y0, tolerance);
        Vector x = res.getV();
        int iterations = res.getIterations();
        double[] result = new double[x.length() - 3];

        for(int i = 0; i < result.length; i++)
        {
            result[i] = Math.abs(x.get(i)) % 2;
        }

        return new ConvolutionDecodeResult(new Vector(result), res.getIterations());
    }

    public static ConvolutionDecodeResult decodeY1g(Vector v, double tolerance) {
        Matrix a1 = makeY1Matrix(v.length() / 2);

        double[] y1d = new double[v.length() / 2];

        for (int i = 0; i < y1d.length; i++) {
            y1d[i] = v.get(i * 2 + 1);
        }

        Vector y1 = new Vector(y1d);

        double[] x0d = new double[y1d.length];

        for (int i = 0; i < x0d.length; i++) {
            x0d[i] = 0;
        }

        Vector x0 = new Vector(x0d);

        IterativeSolutionResult res = Solver.gaussSiedel(a1, x0, y1, tolerance);
        Vector x = res.getV();
        int iterations = res.getIterations();

        double[] result = new double[x.length() - 3];

        for(int i = 0; i < result.length; i++)
        {
            result[i] = Math.abs(x.get(i)) % 2;
        }

        return new ConvolutionDecodeResult(new Vector(result), res.getIterations());
    }
}