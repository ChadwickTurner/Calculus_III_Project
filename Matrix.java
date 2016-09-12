import java.util.Arrays;
import java.util.Queue;
import java.util.LinkedList;
import java.util.Random;

/**
 * Base class for a matrix and matrix operations.
 *
 * @author Joshua Redding
 * @author Chadwick Turner
 * @author Matthew Kennedy
 * @version 0.8
 */
public class Matrix {
    protected double[][] entries;

    /**
     * Creates a matrix from a two-dimensional array of entries
     * @throws IllegalArgumentException if there are no entries for the matrix
     * @param entries the entries in the matrix
     */
    public Matrix(double[][] entries) {
        if (entries == null || entries.length == 0 || entries[0].length == 0) {
            throw new IllegalArgumentException("No entries for the Matrix!");
        }
        this.entries = entries;
    }

    /**
     * Creates a vector-like Matrix of one column
     * @throws IllegalArgumentException if there are no entries for the column
     * @param firstColumn the first (and only) column of the matrix
     */
    public Matrix(double[] firstColumn) {
        if (firstColumn == null || firstColumn.length == 0) {
            throw new IllegalArgumentException("No entries for the Vector!");
        }
        double[][] entries = new double[firstColumn.length][1];
        int i = 0;
        for (double d : firstColumn) {
            entries[i++][0] = d;
        }
        this.entries = entries;
    }

    /**
     * Creates a Matrix by concatenating Vectors
     * @throws IllegalArgumentException if there are no entries for the Matrix
     * @param vectors the vector columns of the Matrix
     */
    public Matrix(Vector... vectors) {
        if (vectors == null || vectors.length == 0 || vectors[0].length() == 0) {
            throw new IllegalArgumentException("No entries for the Matrix!");
        }
        double[][] entries = new double[vectors[0].length()][vectors.length];
        for (int r = 0; r < vectors[0].length(); r++) {
            for(int c = 0; c < vectors.length; c++) {
                entries[r][c] = vectors[c].get(r);
            }
        }
        this.entries = entries;
    }

    public static Matrix randomMatrix(int rows, int columns, int min, int max) {
        Random rnd = new Random();

        double[][] dat = new double[columns][rows];

        for(int r = 0; r < rows; r++)
        {
            for(int c = 0; c < columns; c++)
            {
                dat[r][c] = rnd.nextInt(max - min) + min;
            }
        }

        return new Matrix(dat);
    }

    /**
     * Instantiates a matrix from the user's input.
     * @throws IllegalArgumentException if the row or column input is invalid
     * @return a Matrix created from input
     */
    public static Matrix make() {
        //Scanner for the input
        java.util.Scanner in = new java.util.Scanner(System.in);

        //get row number and column number from user
        System.out.print("# of rows: ");
        int rows;
        try {
            rows = in.nextInt();
        } catch (Exception e) {
            throw new IllegalArgumentException("Invalid row number!");
        }
        if (rows <= 0) {
            throw new IllegalArgumentException("Invalid row number!");
        }
        System.out.print("# of columns: ");
        int columns;
        try {
            columns = in.nextInt();
        } catch (Exception e) {
            throw new IllegalArgumentException("Invalid column number!");
        }
        if (columns <= 0) {
            throw new IllegalArgumentException("Invalid column number!");
        }

        //create appropriately sized 2D array
        double[][] entries = new double[rows][columns];

        //get data to fill
        /*
            Takes data in form x x x x x
                               x x x x x
                               x x x x x
         */
        System.out.println("Type out data below (columns separated by spaces, press Enter for new row):");
        for (int i = 0; i < rows; i++) {
            //take in the row's line as a string
            String row = in.next();
            row += in.nextLine();

            //take apart the string to make column entries
            for (int j = 0; j < columns - 1; j++) {
                entries[i][j] = Double.parseDouble(row.substring(0, row.indexOf(" ")));
                row = row.substring(row.indexOf(" ") + 1);
            }
            entries[i][columns - 1] = Double.parseDouble(row);
        }

        return new Matrix(entries);
    }

    public static Matrix hilbertMatrix(int size)
    {
        double[][] dat = new double[size][size];

        for(int n = 0; n < size * 2; n++)
        {
            double frac = 1.0 / n;

            for(int i = 0; i < n; i++){
                int r = n - i - 1;
                
                if(r < size && i < size)
                {
                   dat[n - i - 1][i] = frac;
                }
            }
        }

        return new Matrix(dat);
    }

    /**
     * Gets a vector column of a Matrix
     * @param column The column of the Matrix the vector represents
     * @throws IllegalArgumentException if the columm is not in the Matrix
     * @return the Vector column
     */
    public Vector getVector(int column) {
        if (column < 0 || column >= entries[0].length) {
            throw new IllegalArgumentException("Column not in Matrix!");
        }
        double[] v = new double[getRows()];
        for (int i = 0; i < getRows(); i++) {
            v[i] = get(i,column);
        }
        return new Vector(v);
    }

    public static Matrix makeBigger(Matrix smaller, int r, int c, int x, int y) {
        Matrix bigger = identity(r, c);
        double[][] temp = bigger.getBacking();
        for (int i = x; i - x < smaller.getRows(); i++) {
            for (int j = y; j - y < smaller.getColumns(); j++) {
                temp[i][j] = smaller.get(i - x, j - y);
            }
        }
        bigger = new Matrix(temp);
        return bigger;
    }

    /**
     * Transposes the current Matrix
     * @return the transposed Matrix
     */
    public Matrix transpose() {
        double[][] a = getBacking();
        double[][] aT = new double[a[0].length][a.length];
        for (int i = 0; i < a.length; i++) {
            for (int j = 0; j < a[0].length; j++) {
                aT[j][i] = a[i][j];
            }
        }
        return new Matrix(aT);
    }

    /**
     * Adds a list of matricies.
     * @param matricies Variable number of matricies to be added
     * @return the sum of the matricies
     */
    public static Matrix add(Matrix ... matricies) {
        double[][] output = new double[matricies[0].getRows()][matricies[0].getColumns()];
        for (int i = 0; i < matricies.length; i++) {
            for (int j = 0; j < output.length; j++) {
                for (int k = 0; k < output[0].length; k++) {
                    output[j][k] += matricies[i].get(j, k);
                }
            }
        }
        return new Matrix(output);
    }


    public static Matrix subtract(Matrix minuend, Matrix subtrahend) {
        int rows = minuend.getRows();
        int cols = minuend.getColumns();

        double[][] output = new double[rows][cols];

        for(int r = 0; r < rows; r++)
        {
            for(int c = 0; c < cols; c++)
            {
                output[r][c] = minuend.get(r, c) - subtrahend.get(r, c);
            }
        }

        return new Matrix(output);
    }


    /**
     * Multiplies the matrix by a scalar
     * @param scalar the scalar by which to multiply
     * @param matrix the matrix by which to multiply
     * @return The Matrix product
     */
    public static Matrix scalarMult(double scalar, Matrix matrix) {
        double[][] entries = matrix.getBacking();
        double[][] output = new double[entries.length][entries[0].length];
        for (int i = 0; i < output.length; i++) {
            for (int j = 0; j < output[0].length; j++) {
                output[i][j] = entries[i][j] * scalar;
            }
        }
        return new Matrix(output);
    }

    /**
     * Multiplies two matricies.
     * @param matrix1 The first factor matrix
     * @param matrix2 The second factor matrix
     * @throws IllegalArgumentException if dimensions are such that multiplication is impossible
     * @return The matrix product
     */
    public static Matrix multiply(Matrix matrix1, Matrix matrix2) {
        if (matrix1.getColumns() != matrix2.getRows()) {
            throw new IllegalArgumentException("Matricies with these dimensions cannot be multiplied!");
        }
        double[][] m1 = matrix1.getBacking();
        double[][] m2 = matrix2.getBacking();
        double[][] output = new double[m1.length][m2[0].length];
        for (int i = 0; i < m1.length; i++) {
            for (int j = 0; j < m2[0].length; j++) {
                double sum = 0;
                for (int k = 0; k < m1[0].length; k++) {
                    sum += m1[i][k] * m2[k][j];
                }
                output[i][j] = sum;
            }
        }
        return new Matrix(output);
    }

    public static Vector multiply(Matrix matrix, Vector vector)
    {
        if (matrix.getColumns() != vector.getRows()) {
            throw new IllegalArgumentException("Matricies with these dimensions cannot be multiplied!");
        }

        double[][] m = matrix.getBacking();
        double[][] v = vector.getBacking();
        double[] output = new double[matrix.getRows()];

        for(int i = 0; i < output.length; i++)
        {
            double sum = 0;
            for(int j = 0; j < v.length; j++)
            {
                sum += m[i][j] * v[j][0];
            }
            output[i] = sum;
        }

        return new Vector(output);
    }

    /**
     * Gets an entry from the matrix.
     * @param row The row index of the entry
     * @param column The column index of the entry
     * @throws IllegalArgumentException if row or column is not in the Matrix
     * @return The entry
     */
    public double get(int row, int column) {
        if (row < 0 || row >= getRows()) {
            throw new IllegalArgumentException("Row does not exist in Matrix!");
        }
        if (column < 0 || column >= getColumns()) {
            throw new IllegalArgumentException("Column does not exist in Matrix!");
        }
        return entries[row][column];
    }

    /**
     * Gets the backing 2D array of entries for the matrix
     * @return the 2D entries array
     */
    public double[][] getBacking() {
        return entries;
    }

    /**
     * Gets the number of rows in the matrix.
     * @return the number of rows
     */
    public int getRows() {
        return entries.length;
    }

    /**
     * Get the number of columns of the matrix
     * @return the number of columns
     */
    public int getColumns() {
        return entries[0].length;
    }

    /**
     * Gives a string representation of the Matrix
     * @return the Matrix in String form
     */
    public String toString() {
        String output = "[";
        for (int i = 0; i < entries.length; i++) {
            for (int j = 0; j < entries[0].length - 1; j++) {
                output += entries[i][j] + ", ";
            }
            output += entries[i][entries[0].length - 1] + "]\n[";
        }
        output = output.substring(0, output.length() - 2);
        return output;
    }

    public static Matrix makeEmpty(int m, int n) {
        double[][] matrix = new double[m][n];
        return new Matrix(matrix);
    }

    /**
     * Generates an identity matrix of the given dimensions
     * @param m the rows for the identity matrix
     * @param n the columns for the identity matrix
     * @return the identity matrix
     */
    public static Matrix identity(int m, int n) {
        double[][] idty = new double[m][n];
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                if (i == j) {
                    idty[i][j] = 1;
                }
            }
        }
        return new Matrix(idty);
    }


    /**
     * Computes (without returning) the LU Decomposition of the Matrix
     */
    public LUFactorizationResult lu_fact() {
        Matrix[] gList = new Matrix[getRows() - 1];
        double[][] l = identity(getRows(), getRows()).getBacking();
        double[][] u = getBacking();
        for (int i = 0; i < gList.length; i++) {
            gList[i] = identity(getRows(), getRows());
            double[][] g = gList[i].getBacking();
            for (int j = i + 1; j < g.length; j++) {
                if (u[i][j] == 0.0) {
                    g[j][i] = 0.0;
                } else {
                    g[j][i] = -u[j][i] / u[i][i];
                }
            }

            u = Matrix.multiply(gList[i], new Matrix(u)).getBacking();
            //compute the inverse of g and add it to L
            for (int j = i + 1; j < g.length; j++) {
                l[j][i] = -g[j][i];
            }
        }

        Matrix ml = new Matrix(l);
        Matrix mu = new Matrix(u);

        // now we find error
        Matrix product = Matrix.multiply(ml, mu);
        Matrix errMatrix = Matrix.subtract(product, this);

        double error = errMatrix.infinityNorm();

        return new LUFactorizationResult(ml, mu, error);
    }

    public QRFactorizationResult qr_fact_househ() {
        double[][] m = getBacking();
        int houseCount;

        if (getRows() <= getColumns()) {
            houseCount = getRows() - 1;
        } else {
            houseCount = getColumns() - 1;
        }

        Matrix[] houseHolders = new Matrix[houseCount];
        houseCount = 0;

        for(int j = 0; j < houseHolders.length; j++) {
            Vector column = new Matrix(m).getVector(j);
            double[] cValues = new double[column.length() - j];
            for (int k = 0; k < cValues.length; k++) {
                cValues[k] = column.getBacking()[k + j][0];
            }

            Vector x = new Vector(cValues);
            double xmag = x.magnitude();
            Vector e1 = identity(getRows() - j,1).getVector(0);

            Vector xe1 = Vector.add(x, Vector.scalarMult(xmag,e1));
            Matrix u = Matrix.scalarMult(1/xe1.magnitude(),xe1);
            Matrix uT = u.transpose();
            Matrix uuT = multiply(u, uT);

            Matrix h = add(identity(uuT.getRows(), uuT.getColumns()), scalarMult(-2, uuT));
            Matrix house = makeBigger(h, getRows(), getColumns(), j, j);
            houseHolders[houseCount] = house;
            m = Matrix.multiply(houseHolders[houseCount], new Matrix(m)).getBacking();
            houseCount++;
        }

        Matrix r = new Matrix(m);
        Matrix q = identity(getRows(),getColumns());
        for (int i = 0; i < houseHolders.length; i++) {
            q = Matrix.multiply(q, houseHolders[i]);
        }

        // now we find error
        Matrix product = Matrix.multiply(q, r);
        Matrix errMatrix = Matrix.subtract(product, this);

        double error = errMatrix.infinityNorm();

        return new QRFactorizationResult(q, r, error);
    }

    
    /**
     * Calculates (without returning) the QR Factorization of the Matrix using Givens Rotation
     */
    public QRFactorizationResult qr_fact_givens() {
        //create a copy of matrix's backing to use for R
        double[][] m = getBacking();

        //find number of Givens rotations needed
        int givensCount = 0;
        int givensIncrement = entries.length - 1;
        for (int i = 0; i < entries[0].length && givensIncrement > 0; i++) {
            givensCount += givensIncrement--;
        }
        Matrix[] givens = new Matrix[givensCount];

        //perform the rotations to find r
        givensCount = 0;
        for (int i = 0; i < m[0].length; i++) {
            for (int j = i + 1; j < m.length; j++) {
                double[][] g = identity(getRows(), getRows()).getBacking();
                if (Math.sqrt(m[i][i] * m[i][i] + m[j][i] * m[j][i]) != 0.0) {
                    g[i][i] = m[i][i] / Math.sqrt(m[i][i] * m[i][i] + m[j][i] * m[j][i]);
                    
                } else {
                    g[i][i] = 0.0;
                }
                g[j][j] = g[i][i];
                if (Math.sqrt(m[i][i] * m[i][i] + m[j][i] * m[j][i]) != 0.0) {
                    g[i][j] = m[j][i] / Math.sqrt(m[i][i] * m[i][i] + m[j][i] * m[j][i]);
                } else {
                    g[i][j] = 0.0;
                }
                g[j][i] = -g[i][j];
                givens[givensCount] = new Matrix(g);
                m = Matrix.multiply(givens[givensCount], new Matrix(m)).getBacking();
                givensCount++;
            }
        }

        Matrix r = new Matrix(m);
        Matrix q = identity(getRows(), getRows());
        for (int i = 0; i < givens.length; i++) {
            q = Matrix.multiply(q, givens[i].transpose());
        }

        // now we find error
        Matrix product = Matrix.multiply(q, r);
        Matrix errMatrix = Matrix.subtract(product, this);

        double error = errMatrix.infinityNorm();

        return new QRFactorizationResult(q, r, error);
    }

    public double infinityNorm()
    {
        int rows = entries.length;
        int cols = entries[0].length;
        
        double max = -1;

        for(int i = 0; i < rows; i++)
        {
            for(int j = 0; j < cols; j++)
            {
                double val = entries[i][j];
                if(val < 0)
                {
                	val = -val;
                }

                if(val > max)
                {
                	max = val;
                }
            }            
        }

        return max;
    }
}
