
class IterativeSolutionResult
{
	final Vector v;
	final int i;

	public IterativeSolutionResult(Vector v, int iter)
	{
		this.v = v;
		this.i = iter;
	}

	public Vector getV()
	{
		return v;
	}

	public int getIterations()
	{
		return i;
	}
}