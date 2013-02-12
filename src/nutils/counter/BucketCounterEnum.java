package nutils.counter;

/**
 * This class provides a counter for all the constants of an enum type
 * @author vegit_000
 *
 * @param <E>
 */
public class BucketCounterEnum<E extends Enum<E>> extends BucketCounter {
	
	public BucketCounterEnum(Class<E> theClass) {		
		super(theClass.getEnumConstants().length, 0);		
	}
	
	public void increment(E e) {
		increment(e, 1);
	}
	
	public void increment(E e, int numToIncrement) {
		super.increment(e.ordinal(), numToIncrement);
	}
	
	public int    getCount(E e)      { return getCount(e.ordinal()); }
	public double getProportion(E e) { return getProportion(e.ordinal()); }
	
	/** @override Override this method to throw an exception. 
	 *  Will throw exception and exit the program.  
	 */
	public void increment(int index) {
		this.increment(index, 1);
	}
	
	/** @override Overrides the method inherited from superclass. 
	 *  Will throw exception and exit the program. 
	 */
	public void increment(int index, int numToIncrement) {
		(new Exception("ERROR: Cannot increment a numeric index for BucketCounterEnum")).printStackTrace();
		System.exit(-1);		
	}
}
