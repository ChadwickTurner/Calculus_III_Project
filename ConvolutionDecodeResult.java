public class ConvolutionDecodeResult
{
	final Vector v;
	final int i;

	public ConvolutionDecodeResult(Vector v, int iter)
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