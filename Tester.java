import java.util.*;

public class Tester
{
    public static void main(String[] args) {
        //  Matrix x = Matrix.make();
        // Matrix x = Matrix.randomMatrix(5,5,-10,10);
        // QRFactorizationResult qr = testQrHouseholder(x);
        // testSolveQR(qr, new Vector(new double[] {7, 8, 9}));
        // Vector x0 = new Vector(new double[] {1, 0});
        // testPower(x, x0, 0.00000001);

        AssignmentFunctions.leslieFunction();
    }

    public static void testPower(Matrix a, Vector x0, double tolerance){
        System.out.println("Power Method:\n");
        PowerMethodResult pmr = Solver.power_method(a, x0, tolerance);
        
        if(pmr.getIterations() == -1)
        {
            System.out.println("The function did not complete within 50000 iterations. Oops.");
            return;
        }
        
        System.out.println("Eigenvalue: " + pmr.getEigenvalue());
        System.out.println();
        System.out.println("Eigenvector: \n" + pmr.getEigenvector());
        System.out.println();
        System.out.println("Iterations: \n" + pmr.getIterations());
    }

    public static void testJacobi(Matrix a, Vector b, Vector x0, double tolerance) {
        System.out.println("Jacobi\n");
        System.out.println(Solver.jacobi(a, x0, b, tolerance));
    }

    public static void testGaussSeidel(Matrix a, Vector b, Vector x0, double tolerance) {
        System.out.println("Gauss-Siedel:\n");
        System.out.println(Solver.gaussSiedel(a, x0, b, tolerance));
    }

    public static QRFactorizationResult testQrGivens(Matrix a) {
        // System.out.println("QR Factorization via Givens");
        QRFactorizationResult resG = a.qr_fact_givens();
        System.out.println("Q:\n" + resG.getQ());
        System.out.println();
        System.out.println("R:\n" + resG.getR());
        System.out.println();
        System.out.println("Error: " + resG.getError());
        System.out.println();

        return resG;
    }

    public static QRFactorizationResult testQrHouseholder(Matrix a) {
        System.out.println("\n\n\n\nQR Factorization via Householder");
        QRFactorizationResult resH = a.qr_fact_househ();
        System.out.println("Q:\n" + resH.getQ());
        System.out.println();
        System.out.println("R:\n" + resH.getR());
        System.out.println();
        System.out.println("Error: " + resH.getError());
        System.out.println();

        return resH;
    }

    public static LUFactorizationResult testLu(Matrix a) {
        System.out.println("\n\n\n\nLU Factorization");
        LUFactorizationResult resLU = a.lu_fact();
        System.out.println("L:\n" + resLU.getL());
        System.out.println();
        System.out.println("U:\n" + resLU.getU());
        System.out.println();
        System.out.println("Error: " + resLU.getError());
        System.out.println();

        return resLU;
    }

    public static void testSolveLU(LUFactorizationResult a, Vector b) {
        System.out.println("\n\n\n\nSolve via LU Decomp, b = ");
        System.out.println(b);
        Vector result = Solver.solve_lu_b(a, b);
        System.out.println("\nResult\n" + result);
        System.out.println();
    }

    public static void testSolveQR(QRFactorizationResult a, Vector b) {
        System.out.println("\n\n\n\nSolve via QR decomp, b = ");
        System.out.println(b);
        Vector result = Solver.solve_qr_b(a, b);
        System.out.println(result);
        System.out.println();
    }
}
