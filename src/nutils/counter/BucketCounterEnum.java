package nutils.counter;

/**
 * This class provides a counter for all the constants of an enum type
 * @author vegit_000
 *
 * @param <E>
 */
public class BucketCounterEnum<E extends Enum<E>> extends BucketCounterCore {
	
	// ========================================================================
	public BucketCounterEnum(Class<E> theClass) {		
		super(theClass.getEnumConstants().length, 0);		
	}
		
	// ========================================================================
	public BucketCounterEnum(BucketCounterEnum<E> rhs) {
		super(rhs);
	}
	
	// ========================================================================
	public BucketCounterEnum<E> getCopy() { return new BucketCounterEnum<E>(this); }
	
	// ========================================================================	
	public void increment(E e) { 
		increment(e, 1); 
	}
	
	// ========================================================================
	public void increment(E e, int numToIncrement) {
		super.incrementBucket(e.ordinal(), numToIncrement);
	}
	
	// ========================================================================
	public int getCount(E e)         { return super.getCountBucket(e.ordinal()); }
	
	// ========================================================================
	public double getProportion(E e) { return super.getProportionBucket(e.ordinal()); }
		
	// ========================================================================
	// A private static enum used for unit testing purposes
	private static enum Alphabet {
		A, B, C, D;
	}
	
	// ========================================================================
	private static void TestRobust() {
		BucketCounterEnum<Alphabet> alphabetCounter  = new BucketCounterEnum<Alphabet>(Alphabet.class);
		BucketCounterEnum<Alphabet> alphabetCounter2 = new BucketCounterEnum<Alphabet>(Alphabet.class);
		
		for (Alphabet letter : Alphabet.values()) {
			for (int i = 0; i < letter.ordinal(); i++) {
				alphabetCounter.increment(letter);
			}
		}
		
		for (Alphabet letter : Alphabet.values()) {
			System.out.printf("%s\t%d\n", letter.name(), alphabetCounter.getCount(letter));
		}
	}
	
	// ========================================================================
	public static void main(String[] args) {
		TestRobust();
	}
}
