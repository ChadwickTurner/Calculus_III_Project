
class PowerMethodResult
{
    private final double eigenvalue;
    private final Vector eigenvector;
    private final int iterations;

    public PowerMethodResult(double eigenvalue, Vector eigenvector, int iterations)
    {
        this.eigenvalue = eigenvalue;
        this.eigenvector = eigenvector;
        this.iterations = iterations;
    }

    public double getEigenvalue()
    {
        return eigenvalue;
    }

    public Vector getEigenvector()
    {
        return eigenvector;
    }

    public int getIterations()
    {
        return iterations;
    }
}