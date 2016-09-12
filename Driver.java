import java.util.Scanner;
import java.util.ArrayList;
import java.io.File;

public class Driver {
    public static void main(String[] args) {
        try
        {
        Scanner scan = new Scanner(System.in);
        while(true) {
            System.out.println("Enter the number of what you would like to test:\n"
                               + "1. LU Factorization \n"
                               + "2. QR Factorization via HouseHolder Reflections \n"
                               + "3. QR Factorization via Givens Rotations \n"
                               + "4. Solve system of linear equations via QR or LU decomp\n"
                               + "5. Generate and encode a random bitstream then decode it or another stream\n"
                               + "6. Power Method \n"
                               + "0. Exit");
            String x = scan.nextLine();
            if (x.equals("1")) {
                Matrix m = getInputMatrix(scan);

                Tester.testLu(m);

            } else if (x.equals("2")) {
                Matrix m = getInputMatrix(scan);

                Tester.testQrHouseholder(m);

            } else if (x.equals("3")) {
                Matrix m = getInputMatrix(scan);
                Tester.testQrGivens(m);
            } else if (x.equals("4")) {
                System.out.println("What would you like to solve?:\n" + 
                                    "1. Import from .dat file\n" + 
                                    "2. Input the matrix by hand\n" +
                                    "3. Use a hilbert matrix");
                int n = scan.nextInt();
                scan.nextLine();
                Matrix a = null;
                Vector b = null;

                switch(n)
                {
                case 1:
                    System.out.print("Enter the file location: ");
                    String filename = scan.nextLine();

                    a = new Matrix(readInput(filename));
                    break;
                case 2:
                    a = Matrix.make();
                    break;
                case 3:
                    System.out.print("How big is your hilbert? ");
                    int hilbSize = scan.nextInt();
                    scan.nextLine();

                    a = Matrix.hilbertMatrix(hilbSize);
                    b = AssignmentFunctions.hilbertB(hilbSize);
                    break;
                default:
                    return;
                }

                if(b == null)
                {
                    ArrayList<Vector> vectors = new ArrayList<Vector>();

                    for(int i = 0; i < a.getColumns() - 1; i++)
                    {
                        Vector v = a.getVector(i);
                        vectors.add(v);
                    }

                    b = a.getVector(a.getColumns() - 1);

                    a = new Matrix(vectors.toArray(new Vector[vectors.size()]));
                }


                System.out.println("Would you like to:\n" + 
                                    "1. Solve Ax=b by LU\n" + 
                                    "2. Solve Ax=b by QR");

                n = scan.nextInt();
                scan.nextLine();

                switch(n)
                {
                case 1:
                    LUFactorizationResult lu = a.lu_fact();
                    Vector lur = Solver.solve_lu_b(lu, b);

                    Vector calcB = Matrix.multiply(a, lur);
                    double error = Vector.add(Vector.scalarMult(-1, calcB), b).magnitude();

                    System.out.println("\nResult\n" + lur);
                    System.out.println("Error: " + error);
                    break;
                case 2:
                    QRFactorizationResult qr = a.qr_fact_givens();
                    Vector qrr = Solver.solve_qr_b(qr, b);

                          calcB = Matrix.multiply(a, qrr);
                    error = Vector.add(Vector.scalarMult(-1, calcB), b).magnitude();

                    System.out.println("\nResult\n" + qrr);
                    System.out.println("Error: " + error);
                    break;
                }



            } else if (x.equals("5")) {
                System.out.print("Enter a value n for the length of the binary stream: ");
                int n = scan.nextInt();
                scan.nextLine();
                Vector y = Vector.randomVector(n);

                Matrix a0 = Convolution.makeY0Matrix(n + 3);
                Matrix a1 = Convolution.makeY1Matrix(n + 3);

                System.out.println("A0:");
                System.out.println(a0);

                System.out.println("A1:");
                System.out.println(a1);

                Vector z = Convolution.encode(y);

                System.out.println("Y:");

                for (int i = 0; i < z.length(); i++) {
                    System.out.print((int)(z.get(i)+0.1));
                }

                System.out.println("\n");

                // System.out.println("What next?\n1. Decode this bitstream\n2. Input a bitstream to decode");
                // n = scan.nextInt();
                // scan.nextLine();

                // switch(n)
                // {
                // case 1:
                //     // do nothing
                //     break;
                // case 2:

                //     break;
                // }

                System.out.println("How to decode?\n1. Jacobi on A0\n2. Jacobi on A1\n3. Gauss-Seidel on A0\n4. Gauss-Seidel on A1");
                n = scan.nextInt();
                scan.nextLine();

                ConvolutionDecodeResult result = null;

                System.out.print("Enter a tolerance. (Ex: for 10E-8, enter .00000001) : ");
                double tolerance = scan.nextDouble();
                scan.nextLine();

                switch(n)
                {
                case 1:
                    result = Convolution.decodeY0j(z, tolerance);
                    break;
                case 2:
                    result = Convolution.decodeY1j(z, tolerance);
                    break;
                case 3:
                    result = Convolution.decodeY0g(z, tolerance);
                    break;
                case 4:
                    result = Convolution.decodeY1g(z, tolerance);
                    break;
                }

                if(result.getIterations() == -1)
                {
                    System.out.println("The function did not complete within 50000 iterations. Oops.");
                }
                else
                {
                    System.out.println("It took " + result.getIterations() + " iterations to converge");
                    System.out.println("Initial code:");

                    for (int i = 0; i < result.getV().length(); i++) {
                        System.out.print((int)(result.getV().get(i)+0.1));
                    }
                }
            } else if (x.equals("6")) {
                System.out.println("Would you like to: \n"
                                + "1. Input a .dat file \n"
                                + "2. Input the matrix by hand \n"
                                + "3. Use the example Leslie Matrix");
                x = scan.nextLine();

                Matrix m;

                if (x.equals("1")) {
                    System.out.print("Enter the file location: ");
                    String filename = scan.nextLine();
                    m = new Matrix(readInput(filename));
                } else if (x.equals("2")) {
                    m = Matrix.make();
                } else {
                    double[][] ad = new double[9][9];
                    ad[0][1] = 1.2;
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

                    m = new Matrix(ad);
                }

                System.out.print("Enter a tolerance. (Ex: for 10E-8, enter .00000001) : ");
                double tol = scan.nextDouble();
                scan.nextLine();

                double[] x0d = new double[m.getRows()];

                for(int i = 0; i < x0d.length; i++)
                {
                    x0d[i] = 1;
                }

                Vector x0 = new Vector(x0d);
                Tester.testPower(m,x0,tol);
            } else {
                return;
            }
            System.out.println("\n\n\n");
        }
    } catch(Exception exc){
        System.out.println(exc);
    }
    }

    private static Matrix getInputMatrix(Scanner scan) throws Exception
    {
            System.out.println("Would you like to: \n"
                            + "1. Input a .dat file \n"
                            + "2. Input the matrix by hand");
            String x = scan.nextLine();

            Matrix m = null;

            if (x.equals("1")) {
                System.out.print("Enter the file location: ");
                String filename = scan.nextLine();

                return new Matrix(readInput(filename));
            } else {
                return Matrix.make();
            }
    }

    private static double[][] readInput(String filename) throws java.io.FileNotFoundException {
        //when you copy this over, do not forget to import java.io.File, java.util.ArrayList, and java.util.Scanner
        File input;
        try {
            input = new File(filename);
        } catch (Exception e) {
            throw new java.io.FileNotFoundException("File not found");
        }
        Scanner in = new Scanner(input);
        ArrayList<String> rows = new ArrayList<>();
        while (in.hasNext()) {
            rows.add(in.nextLine());
        }
        ArrayList<ArrayList<Double>> matrix = new ArrayList<>();
        for(int i = 0; i < rows.size(); i++) {
            matrix.add(new ArrayList<>());
            while (rows.get(i).length() > 0) {
                int spaceIndex = rows.get(i).indexOf(" ");
                int commaIndex = rows.get(i).indexOf(",");
                if (spaceIndex > 0 || commaIndex > 0) {
                    if (spaceIndex < 0 || commaIndex < 0) {
                        spaceIndex = Math.max(spaceIndex, commaIndex);
                        commaIndex = spaceIndex;
                    }
                    matrix.get(i).add(Double.parseDouble(rows.get(i).substring(0, Math.min(spaceIndex, commaIndex))));
                    rows.set(i, rows.get(i).substring(Math.max(spaceIndex, commaIndex) + 1));
                } else {
                    matrix.get(i).add(Double.parseDouble(rows.get(i)));
                    rows.set(i, "");
                }
            }
        }
        double[][] entries = new double[matrix.size()][matrix.get(0).size()];
        for (int i = 0; i < entries.length; i++) {
            for (int j = 0; j < entries[0].length; j++) {
                entries[i][j] = matrix.get(i).remove(0);
            }
        }
        return entries;
    }
}