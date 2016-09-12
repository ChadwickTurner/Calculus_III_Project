
class LUFactorizationResult
{
	private final Matrix l, u;
	private final double e;

	public LUFactorizationResult(Matrix l, Matrix u, double error)
	{
		this.l = l;
		this.u = u;
		this.e = error;
	}

	public Matrix getL()
	{
		return l;
	}

	public Matrix getU()
	{
		return u;
	}

	public double getError()
	{
		return e;
	}
}