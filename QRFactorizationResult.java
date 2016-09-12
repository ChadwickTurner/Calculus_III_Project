
class QRFactorizationResult
{
	final Matrix q, r;
	final double e;

	public QRFactorizationResult(Matrix q, Matrix r, double error)
	{
		this.q = q;
		this.r = r;
		this.e = error;
	}

	public Matrix getQ()
	{
		return q;
	}

	public Matrix getR()
	{
		return r;
	}

	public double getError()
	{
		return e;
	}
}