package nutils.counter;

import nutils.ArrayUtils;
import nutils.counter.BucketCounterEnum.Alphabet;

/**
 * This class creates a matrix for an enumerated type.  Each matrix cell stores a count.
 * This matrix is not symmetric across the diagonal, as being below or above the diagonal
 * can represent different meanings.  
 * 
 * @author Ninad Dewal
 *
 */

public class BucketCounterEnumMatrix<E extends Enum<E>> {

	public static int DefaultInitialValue = 0;
	
	// ========================================================================
	// MEMBER VARIABLES
	int[][] mCounts;
	
	// ========================================================================
	public BucketCounterEnumMatrix(Class<E> theClass, int initialValue) {
		E[] enumConstants = theClass.getEnumConstants();
		
		mCounts = new int[enumConstants.length][enumConstants.length];
		clear(initialValue);
	}
	
	// ========================================================================
	public BucketCounterEnumMatrix(Class<E> theClass) {
		this(theClass, DefaultInitialValue);
	}
		
	// ========================================================================
	public void clear(int initialValue) { ArrayUtils.arrayFill(mCounts, initialValue); }
	
	// ========================================================================
	public void clear()                 { ArrayUtils.arrayFill(mCounts, DefaultInitialValue); }
	
	// ========================================================================
	public int getCount(E row, E column) { return mCounts[row.ordinal()][column.ordinal()]; }
	
	// ========================================================================
	public void setCount(E row, E column, int value) { mCounts[row.ordinal()][column.ordinal()] = value; }	
	
	// ========================================================================
	public int increment(E row, E column) { return increment(row, column, 1); }
	
	// ========================================================================
	public int increment(E row, E column, int incrementValue) {		
		mCounts[row.ordinal()][column.ordinal()] += incrementValue;
		return mCounts[row.ordinal()][column.ordinal()];
	}
		
	// ========================================================================
	private static void TestMatrix() {
		BucketCounterEnumMatrix<Alphabet> matrix = new BucketCounterEnumMatrix<>(Alphabet.class);
		int counter = 0;
		int increment = 2;
		for (Alphabet letter1 : Alphabet.values()) {
			for (Alphabet letter2 : Alphabet.values()) {
				matrix.increment(letter1, letter2, ++counter);
			}
		}
		
		for (Alphabet letter1 : Alphabet.values()) {
			for (Alphabet letter2 : Alphabet.values()) {
				System.out.printf("%d\t", matrix.getCount(letter1, letter2));
			}
			System.out.println("");
		}

	}
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		TestMatrix();
	}

}
