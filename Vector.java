import java.util.Random;

/**
 * Base class for a vector and vector operations
 *
 * @author Joshua Redding
 * @version 0.6
 */
public final class Vector extends Matrix {
    /**
     * Creates a vector from an array of entries
     *
     * @param entries The array of entries for the vector
     */
    public Vector(double[] values) {
    	super(values);
    }

    public static Vector scalarMult(double scalar, Vector v)
    {
    	double[][] entries = v.getBacking();
    	double[] output = new double[entries.length];
    	for(int i = 0; i < entries.length; i++)
    	{
    		output[i] = entries[i][0] * scalar;
    	}
    	return new Vector(output);
    }

    public static Vector randomVector(int n)
    {
        double[] data = new double[n];
        Random rng = new Random();
        for (int i =0; i < n; i++) {
            data[i] = rng.nextInt(2);
        }
        return new Vector(data);
    }

    /**
     * Instantiates a vector from the user's input.
     *
     * @return a Vector created from input
     */
    public static Vector make() {
        //Scanner for the input
        java.util.Scanner in = new java.util.Scanner(System.in);

        //get entry number from user
        System.out.print("# of entries: ");
        int length = in.nextInt();

        //create the appropriately sized array
        double[] entries = new double[length];

        //get data to fill
        //Takes data in form x x x x x
        System.out.println("Type out data below:");

        //take in the vector's line as a string
        String entry = in.next();
        entry += in.nextLine();

        //take apart the string to make entries
        for (int i = 0; i < length - 1; i++) {
            entries[i] = Double.parseDouble(entry.substring(0, entry.indexOf(" ")));
            entry = entry.substring(entry.indexOf(" ") + 1);
        }
        entries[length - 1] = Double.parseDouble(entry);

        return new Vector(entries);
    }

    /**
     * Adds the list of vectors.
     * @param vectors the vectors to be added
     * @return The sum of the vectors
     */
    public static Vector add(Vector ... vectors) {
        double[] output = new double[vectors[0].length()];
        for (int i = 0; i < vectors.length; i++) {
            for (int j = 0; j < output.length; j++) {
                output[j] += vectors[i].get(j);
            }
        }
        return new Vector(output);
    }

    /**
     * Computes the dot product of two vectors
     * @param v1 The first vector
     * @param v2 The second vectogr
     * @return The dot product
     */
    public static double dot(Vector vector1, Vector vector2) {
        Matrix product = multiply(vector1, vector2);

     	double sum = 0.0;
        for (double d : product.getBacking()[0]) {
            sum += Math.pow(d,2);
        }
		
        return sum;
    }

    public static Vector add(Vector v1, Vector v2)
    {
        if(v1.length() != v2.length())
        {
            throw new IllegalArgumentException("One does not simply add vectors with mismatched sizes.");
        }

        int len = v1.length();

        double[] arr = new double[len];

        for(int i = 0; i < len; i++)
        {
            arr[i] = v1.get(i) + v2.get(i);
        }

        return new Vector(arr);
    }

    /**
     * Gets an entry from the vector
     * @param entry the index of the entry
     * @return the entry
     */
    public double get(int entry) {
        return entries[entry][0];
    }

    /**
     * Gets the length of the vector
     * @return the length
     */
    public int length() {
        return entries.length;
    }

    /**
     * Gives a string representation of the Matrix
     *
     * @return the Matrix in String form
     */
    public String toString() {
        String output = "";
        for (int i = 0; i < entries.length - 1; i++) {
            output += "[" + entries[i][0] + "] \n";
        }
        output += "[" +entries[entries.length - 1][0] + "]";
        return output;
    }

    public double magnitudeSquared() {
		double sum = 0.0;
        for (int i = 0; i < getRows(); i++) {
            double d = getBacking()[i][0];
            sum += Math.pow(d,2);
        }
		
		return sum;
    }

    public double magnitude() {
        return Math.sqrt(magnitudeSquared());
    }

    public boolean equals(Object other) {
        if (!(other instanceof Vector)) {
            System.out.println("fail on type");
            return false;
        }
        for (int i = 0; i < this.length(); i++) {
            if ((this.get(i) - ((Vector) other).get(i))> 0.01) {
                System.out.println("x" + this.get(i));
                System.out.println("y" + ((Vector) other).get(i));
                System.out.println("faile at " + i);
                return false;
            }
        }
        return true;
    }

    // public static void main(String[] args) {
    //     // double[] d = {1,2,2};
    //     Vector v = new Vector(d);
    //     System.out.println(v.magnitude());
    // }
}