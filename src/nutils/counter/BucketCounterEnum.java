package nutils.counter;

import nutils.CloneInf;

/**
 * This class provides a counter for all the constants of an enum type
 * @author vegit_000
 *
 * @param <E>
 */
public class BucketCounterEnum<E extends Enum<E>> extends BucketCounterCore implements CloneInf<BucketCounterEnum<E>> {
	
	// ========================================================================
	public BucketCounterEnum(Class<E> theClass) {		
		super(theClass.getEnumConstants().length, 0);
	}
	
	// ========================================================================
	public BucketCounterEnum(Class<E> theClass, int initialCountForEachBucket) {
		this(theClass);
		
		E[] enumConstants = theClass.getEnumConstants();
		for (E enumConstant : enumConstants) {
			increment(enumConstant, initialCountForEachBucket);
		}
	}
	
	// ========================================================================
	public BucketCounterEnum(Class<E> theClass, int[] countsForEnums) {
		this(theClass);

		E[] enumConstants = theClass.getEnumConstants();
		int index = -1;
		for (E enumConstant : enumConstants) {
			increment(enumConstant, countsForEnums[++index]);
		}
	}	
		
	// ========================================================================
	public BucketCounterEnum(BucketCounterEnum<E> rhs) {
		super(rhs);
	}
	
	// ========================================================================
	@Override
	public BucketCounterEnum<E> makeClone() { return makeClone(true); }

	// ========================================================================
	@Override
	public BucketCounterEnum<E> makeClone(boolean deepCopy) {
		return new BucketCounterEnum<E>(this);
	}
	
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
	protected static enum Alphabet {
		A, B, C, D;
	}
	
	// ========================================================================
	private static void TestRobust() {
		BucketCounterEnum<Alphabet> alphabetCounter  = new BucketCounterEnum<Alphabet>(Alphabet.class, 5);
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
