package nutils;

/**
 * A convenience class to cast primtives to one another.  The intent is that developers
 * can track the casting than 
 * @author Ninad Dewal
 *
 */
public class Cast {

	public static int toInt(double d) { return (int) d; }
	public static int toInt(long l)   { return (int) l; }
	public static int toInt(float f)  { return (int) f; }	
	
	public static long toLong(double d) { return Math.round(d); }
	public static long toLong(float  f) { return Math.round(f); }
	
	public static float toFloat(double d) { return (float) d; }
	public static float toFloat(long l)   { return (float) l; }
	
	public static<T, S extends T> T upcast(S source)   { return (T) source; }
	public static<S, T extends S> T downcast(S source) { return (T) source; }

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub

	}

}
