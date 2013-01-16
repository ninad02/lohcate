package shared;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.Iterator;

import com.carrotsearch.hppc.DoubleArrayList;
import com.carrotsearch.hppc.LongArrayList;

public class ArrayUtils {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		DynamicBucketCounter.TestDynamicBucketCounter();
	}

	/** Copies arr2 into arr1. */
	public static void arrayCopy(boolean[] arr1, final boolean[] arr2, int length) {
		for (int i = 0; i < length; i++) {
			arr1[i] = arr2[i];
		}
	}

	/** Copies arr2 into arr1. */
	public static void arrayCopy(char[] arr1, final char[] arr2, int length) {
		for (int i = 0; i < length; i++) {
			arr1[i] = arr2[i];
		}
	}

	/** Copies arr2 into arr1. */
	public static void arrayCopy(long[] arr1, final long[] arr2, int length) {
		for (int i = 0; i < length; i++) {
			arr1[i] = arr2[i];
		}
	}

	/** Copies arr2 into arr1. */
	public static void arrayCopy(byte[] arr1, final byte[] arr2, int length) {
		for (int i = 0; i < length; i++) {
			arr1[i] = arr2[i];
		}
	}

	/** Copies arr2 into arr1. */
	public static void arrayCopy(int[] arr1, final int[] arr2, int length) {
		for (int i = 0; i < length; i++) {
			arr1[i] = arr2[i];
		}
	}
	
	/** Copies arr2 into arr1. */
	public static void arrayCopy(short[] arr1, final short[] arr2, int length) {
		for (int i = 0; i < length; i++) {
			arr1[i] = arr2[i];
		}
	}

	/** Copies arr2 into arr1. */
	public static<T> void arrayCopy(T[] arr1, final T[] arr2, int length) {
		for (int i = 0; i < length; i++) {
			arr1[i] = arr2[i];
		}
	}

	/** Copies arr2 into arr1. */
	public static void arrayCopy(float[] arr1, final float[] arr2, int length) {
		for (int i = 0; i < length; i++) {
			arr1[i] = arr2[i];
		}
	}
	
	/** Copies arr2 into arr1. */
	public static void arrayCopy(double[] arr1, final double[] arr2, int length) {
		for (int i = 0; i < length; i++) {
			arr1[i] = arr2[i];
		}
	}

	public static<T> void reverseArray(T[] theArray) {
		T tempVal;
		int oppositeIndex;
		for (int i = 0; i < (theArray.length >>> 2); i++) {
			oppositeIndex = theArray.length - i - 1;
			tempVal = theArray[i];
			theArray[i] = theArray[oppositeIndex];
			theArray[oppositeIndex] = tempVal;
		}
	}
	
	public static void reverseArray(int[] theArray) {
		int tempVal;
		int oppositeIndex;
		for (int i = 0; i < theArray.length / 2; i++) {
			oppositeIndex = theArray.length - i - 1;
			tempVal = theArray[i];
			theArray[i] = theArray[oppositeIndex];
			theArray[oppositeIndex] = tempVal;
		}
	}

	public static void reverseArray(double[] theArray) {
		double tempVal;
		int oppositeIndex;
		for (int i = 0; i < theArray.length / 2; i++) {
			oppositeIndex = theArray.length - i - 1;
			tempVal = theArray[i];
			theArray[i] = theArray[oppositeIndex];
			theArray[oppositeIndex] = tempVal;
		}
	}

	public static void reverseArray(char[] theArray) {
		char tempVal;
		int oppositeIndex;
		for (int i = 0; i < theArray.length / 2; i++) {
			oppositeIndex = theArray.length - i - 1;
			tempVal = theArray[i];
			theArray[i] = theArray[oppositeIndex];
			theArray[oppositeIndex] = tempVal;
		}
	}

	/** Given an array, this reverses the array. */
	public static void reverseArray(byte[] theArray) {
		byte tempVal;
		int oppositeIndex;
		for (int i = 0; i < theArray.length / 2; i++) {
			oppositeIndex = theArray.length - i - 1;
			tempVal = theArray[i];
			theArray[i] = theArray[oppositeIndex];
			theArray[oppositeIndex] = tempVal;
		}
	}

	public static void convertArray(String[] toArray, final double[] fromArray) {
		for (int i = 0; i < toArray.length; i++) {
			toArray[i] = "" + fromArray[i];
		}
	}

	public static int[] convertArrayToInt(String[] theArray) {
		int[] rV = new int[theArray.length];
		for (int i = 0; i < theArray.length; i++) {
			rV[i] = Integer.parseInt(theArray[i]);
		}
		return rV;
	}

	public static double[] convertArray(String[] theArray) {
		double[] rV = new double[theArray.length];
		for (int i = 0; i < theArray.length; i++) {
			rV[i] = Double.parseDouble(theArray[i]);
		}
		return rV;
	}

	public static double[] convertArray(String[] theArray, int startIndex) {
		int newArrayLength = theArray.length - startIndex;
		double[] rV = new double[newArrayLength];
		for (int i = 0; i < newArrayLength; i++) {
			rV[i] = Double.parseDouble(theArray[i + startIndex]);
		}
		return rV;
	}

	public static double[] convertArray(int[] theArray) {
		double[] rV = new double[theArray.length];
		for (int i = 0; i < theArray.length; i++) {
			rV[i] = theArray[i];
		}
		return rV;
	}

	/** Converts a collection of Integers to a primitive list. */
	public static int[] getPrimitiveIntegerArray(Collection<Integer> theList) {
		int[] rV = new int[theList.size()];
		int index = 0;
		for (Iterator<Integer> iter = theList.iterator(); iter.hasNext(); ) {
			rV[index++] = iter.next().intValue();
		}
		return rV;
	}

	/** Converts a collection of Double to a primitive list. */
	public static double[] getPrimitiveDoubleArray(Collection<Double> theList) {
		double[] rV = new double[theList.size()];
		int index = 0;
		for (Iterator<Double> iter = theList.iterator(); iter.hasNext(); ) {
			rV[index++] = iter.next().doubleValue();
		}
		return rV;
	}

	/** Given an array, this returns the sum of the array elements. */
	public static int arraySum(int[] theArray) {
		int sum = 0;
		for (int i = 0; i < theArray.length; i++) { sum += theArray[i]; }
		return sum;
	}

	/** Given an array, this returns the sum of the array elements. */
	public static int arraySum(byte[] theArray) {
		int sum = 0;
		for (int i = 0; i < theArray.length; i++) { sum += theArray[i]; }
		return sum;
	}

	
	/** Given an array of doubles, this returns the index of the maximum element.  If there
	 *  are multiple maximal elements, this returns the index of the first one.  The second
	 *  argument is the index at which searching should begin. */
	public static int getIndexOfMaxElement(double[] theArray, int indexToStartSearch) {
		int maxIndex = 0;
		double maxElement = Double.MAX_VALUE * -1; 
		for (int i = indexToStartSearch; i < theArray.length; i++) {
			if (theArray[i] > maxElement) {
				maxElement = theArray[i];
				maxIndex = i;
			}
		}
		return maxIndex;
	}
	
	/** Given an array of doubles, this returns the index of the maximum element.  If there
	 *  are multiple maximal elements, this returns the index of the first one.  The second
	 *  argument is the index at which searching should begin. */
	public static int getIndexOfMaxElement(int[] theArray, int indexToStartSearch) {
		int maxIndex = -1;
		int maxElement = Integer.MIN_VALUE; 			
		for (int i = indexToStartSearch; i < theArray.length; i++) {
			if (theArray[i] > maxElement) {
				maxElement = theArray[i];
				maxIndex = i;
			}
		}
		return maxIndex;
	}

	/** Given an array of doubles, this returns the index of the minimum element.  If there
	 *  are multiple minimum elements, this returns the index of the first one.  The second
	 *  argument is the index at which searching should begin. */
	public static int getIndexOfMinElement(double[] theArray, int indexToStartSearch) {
		int minIndex = -1;
		double minElement = Double.MAX_VALUE; 
		for (int i = indexToStartSearch; i < theArray.length; i++) {
			if (theArray[i] < minElement) {
				minElement = theArray[i];
				minIndex = i;
			}
		}
		return minIndex;
	}

	/** Given an array, this normalizes all the values in the array with the given normalizing constant. */
	public static double[] normalizeArray(double[] theArray, double normalizingFactor) {
		for (int i = 0; i < theArray.length; i++) {
			theArray[i] /= normalizingFactor;
		}
		return theArray;
	}

	/** Prints the double array to screen. */ 
	public static void printDoubleArray(double[] theArray) {
		for (int i = 0; i < theArray.length; i++) {
			if (i > 0) {
				System.out.print("\t");
			}
			System.out.print(theArray[i]);			
		}
		System.out.println("");
	}

	/** Prints the int array to screen. */ 
	public static void printIntArray(int[] theArray) {
		for (int i = 0; i < theArray.length; i++) {
			if (i > 0) {
				System.out.print("\t");
			}
			System.out.print(theArray[i]);			
		}
		System.out.println("");
	}
	
	/** Prints the long array to screen. */ 
	public static void printLongArray(long[] theArray) {
		for (int i = 0; i < theArray.length; i++) {
			if (i > 0) {
				System.out.print("\t");
			}
			System.out.print(theArray[i]);			
		}
		System.out.println("");
	}

	public static void printCollectionToScreen(Collection<?> c) {
		for (Object o : c) {
			System.out.println(o.toString());
		}
	}

	/** Given a collection, a dummy object, and a count, this adds "count" instances of 
	 *  the dummy object to the collection (only a shallow adding, not a deep adding).
	 *  
	 *  This method is unsafe in that in can produce a run-time exception if the 
	 *  collection is defined to hold a particular type and the object is not that type. 
	 */
	public static <T> void fillCollection(Collection<T> c, T o, int count) {
		for (int i = 0; i < count; i++) {
			c.add(o);
		}
	}

	// ========================================================================
	// ===== INNER CLASS =====
	// ========================================================================
	public static abstract class NumArray {
		public int mArrayLength;
		
		public NumArray() { reset(); }		
		public void reset() { mArrayLength = 0; }
	}

	public static class IntArray extends NumArray {
		public int[] mArray;
					
		public IntArray(int arrayLength) { 
			super();
			mArray = new int[arrayLength];
		}
		
		public int add(int num) { 
			mArray[mArrayLength] = num; 
			return ++mArrayLength; 
		}
	}

	public static class DoubleArray extends NumArray {
		public double[] mArray;
		
		public DoubleArray(int arrayLength) { 
			super();
			mArray = new double[arrayLength];
		}
		
		public int add(double num) { 
			mArray[mArrayLength] = num; 
			return ++mArrayLength; 
		}
	}
	
	// ========================================================================
	// ===== INNER CLASS =====
	// ========================================================================
	public static class BucketCounter {
		private int[] mArray;		
		private int mIndexStart;
		private int mSumOfCounts;
		private int mTotalSum;
		
		public BucketCounter(int arrayLength, int valueStart) {
			mIndexStart = valueStart;
			mArray = new int[arrayLength];
			clear();
		}
		
		public void clear() { 
			Arrays.fill(mArray, 0); 
			mSumOfCounts = 0; 	
			mTotalSum = 0;
		}
		
		public void increment(int bucket) { 
			increment(bucket, 1);
		}
		
		public void increment(int bucket, int numTimesToIncrement) {
			mArray[bucket - mIndexStart] += numTimesToIncrement;
			mSumOfCounts                 += numTimesToIncrement;
			mTotalSum                    += (numTimesToIncrement * bucket);
		}
		
		public int getLength() { return mArray.length; }
		public int getCount(int bucket) { return mArray[bucket - mIndexStart]; }
		public double getProportion(int bucket) { return (double) getCount(bucket) / (double) mSumOfCounts; }
		public double calcMean() { return (double) mTotalSum / (double) mSumOfCounts; }
		
		/*
		public ArrayUtils.BucketCounter downsampleCounts(double probabilityToRetainCount) {
			int thisLength = getLength();
			ArrayUtils.BucketCounter bcSubsampled = new ArrayUtils.BucketCounter(thisLength, mIndexStart);
			
			for (int b = 0; b < thisLength; b++) {
				int trueBucketValue = b + mIndexStart;
				int numSitesInBucket = getCount(trueBucketValue);
				for (int j = 0; j < numSitesInBucket; j++) {
					int downsampledNumReads = NumberUtils.numSuccessesInTrials(trueBucketValue, probabilityToRetainCount, 0);
					bcSubsampled.increment(downsampledNumReads);
				}
			}	
			return bcSubsampled;
		}*/
		
		public boolean add(BucketCounter bc) {
			if (this.mArray.length != bc.mArray.length) return false;
			if (this.mIndexStart != bc.mIndexStart) return false;
			
			int bcLength = bc.getLength();
			for (int b = 0; b < bcLength; b++) {
				int trueOtherBucketValue = b + bc.mIndexStart;
				int otherCount = bc.getCount(trueOtherBucketValue);
				this.increment(trueOtherBucketValue, otherCount);
			}
			return true;
		}
		
		public void print(PrintStream out) {			
			print(out, true);
		}
		
		public void print(PrintStream out, boolean printZeroCountElements) {
			for (int i = 0; i < getLength(); i++) {
				if (printZeroCountElements || (getCount(i) > 0)) {
					out.println(i + "\t" + getCount(i) + "\t" + getProportion(i));
				}
			}
		}
	}

	// ========================================================================
	
	
	// ========================================================================
	public static class DynamicBucketCounter {

		LongArrayList mArray;
		private int mSumOfCounts;
		private int mTotalSum;
		
		public DynamicBucketCounter() {
			mArray = new LongArrayList();
			clear();
		}
		
		public void clear() { 
			mArray.clear();
			mSumOfCounts = 0; 	
			mTotalSum = 0;
		}
		
		public int getCount(int key) {
			int indexOfKey = getIndexOfKey(key);
			if (indexOfKey < 0) return indexOfKey;
			return getCountAtIndex(indexOfKey);			
		}
		
		private int getCountAtIndex(int indexOfKey) {
			return (int) IntExtractorLSB.extractValue(mArray.get(indexOfKey));
		}
		
		private int getKeyAtIndex(int indexOfKey) {
			return (int) IntExtractorMSB.extractValue(mArray.get(indexOfKey));
		}
		
		private int getIndexOfKey(int key) {
			return binarySearchValue(key, mArray, IntExtractorMSB);
		}
		
		public int getKeyLast() { return (int) IntExtractorMSB.extractValue(mArray.get(mArray.size() - 1)); }
						
		public int incrementCount(int key) { return incrementCount(key, 1); }
		
		public int incrementCount(int key, int numToIncrement) {
			int indexOfKey = getIndexOfKey(key);
			if (indexOfKey < 0) {
				mArray.insert(getInsertPoint(indexOfKey), getKeyValueAsLong(key, numToIncrement));
				return numToIncrement;
			} else {
				int newCount = numToIncrement + getCountAtIndex(indexOfKey);
				mArray.set(indexOfKey, getKeyValueAsLong(key, newCount));		
				return newCount;
			}
		}
		
		private static long getKeyValueAsLong(int key, int value) {
			return ( 0L | (((long) key) << Integer.SIZE) | value );
		}
		
		
		public int getKeyWithMaxCount() {
			long keyValueMax = getKeyValueWithMaxCount();
			return (int) IntExtractorMSB.extractValue(keyValueMax);
		}
		
		public long getKeyValueWithMaxCount() {
			return searchMaxValue(mArray, IntExtractorLSB);
		}
		
		public static void TestDynamicBucketCounter() {
			int numVals = 100;		
			DynamicBucketCounter dbc = new DynamicBucketCounter();
			
			for (int i = 1; i <= numVals; i++) {				
				for (int j = 0; j < i; j++) {
					dbc.incrementCount(i);
				}
			}
			
			for (int i = 1; i <= numVals; i++) {
				System.out.printf("%d\t%d\n", i, dbc.getCount(i));
			}
			
			System.out.println("Key with max count: " + dbc.getKeyWithMaxCount());
		}
		
		public DoubleArrayList[] toArrayListDouble() {
			DoubleArrayList[] rV = new DoubleArrayList[2];
			rV[0] = new DoubleArrayList();
			rV[1] = new DoubleArrayList();
			
			for (int i = 0; i < mArray.size(); i++) {				
				rV[0].add(getKeyAtIndex(i));
				rV[1].add(getCountAtIndex(i));
			}
			return rV;
		}
		
	}
	
	// ========================================================================
	public static int getInsertPoint(int negativeArrayIndex) { return -(negativeArrayIndex + 1); }
	
	// ========================================================================
	/** Assumes that both lists are of the same size.  If not, it returns null. */ 
	public static double[][] convertToTwoDimensionalArrays(Collection<Double> list0, Collection<Double> list1) {
		if (list0.size() != list1.size()) return null;
			
		int listSize = list0.size();
		double[][] rV = new double[2][listSize];
		Iterator<Double> iter0 = list0.iterator();
		Iterator<Double> iter1 = list1.iterator();
		for (int i = 0; i < listSize; i++) {
			rV[0][i] = iter0.next().doubleValue(); 
			rV[1][i] = iter1.next().doubleValue();	
		}
		return rV;
	}

	// ========================================================================
	/** Given a String list in braces {}, this returns the list of numbers contained
	 *  within those braces.  This also allows for ranges to be included.  For example,
	 *  an input of {1, 2, 3-5, 7} would return the list {1, 2, 3, 4, 5, 7}
	 */ 
	public static ArrayList<Integer> getIntListFromStringForm(String listAsStr, boolean performSort) {
		char braceLeft  = '{';
		char braceRight = '}';
		char comma = ',';
		char hyphen = '-';
		
		listAsStr = listAsStr.replace(braceLeft, comma);
		listAsStr = listAsStr.replace(braceRight, comma);		
		String components[] = listAsStr.split("" + comma);
		ArrayList<Integer> finalList = new ArrayList<Integer>();
		
		for (int i = 0; i < components.length; i++) {
			String item = components[i].trim();
			if (item.equals("")) continue;
			
			int hyphenIndex = item.indexOf(hyphen);
			if (hyphenIndex < 0) {
				finalList.add(new Integer(item));
			} else {
				String numbers[] = item.split("" + hyphen);
				if (numbers.length == 2) {
					int firstNum = Integer.parseInt(numbers[0]);
					int lastNum = Integer.parseInt(numbers[numbers.length - 1]);
					if (firstNum > lastNum) {
						int tempNum = firstNum;
						firstNum = lastNum;
						lastNum = tempNum;
					}
					for (int j = firstNum; j <= lastNum; j++) {
						finalList.add(new Integer(j));
					}
				}
			}			
		}
		
		if (performSort) { Collections.sort(finalList); }		
		return finalList;
	}

	public static ArrayList<Integer> getIntListFromStringForm(String listAsStr) {
		return getIntListFromStringForm(listAsStr, true);
	}

	// ========================================================================
	/** Given a String list of double values in braces {}, this returns the list
	 *  of doubles within those braces.  
	 */
	public static double[] getDoubleListFromStringForm(String listAsStr, boolean performSort) {	
		listAsStr = stripBraces(listAsStr); 				
		String components[] = listAsStr.split(",");
		double[] rv = new double[components.length];
		for (int i = 0; i < components.length; i++) {
			rv[i] = Double.parseDouble(components[i].trim());
		}
		
		if (performSort) { Arrays.sort(rv); }
		return rv;
	}

	// ========================================================================
	/** Given a String list in the format: (list1;list2;...;list_n), where each list
	 *  is of the form {3,2.6,9} (for example), this returns a list of lists represented
	 *  by the string
	 */
	public static ArrayList<double[]> getListOfDoubleListsFromStringForm(String listAsStr, boolean performSortInEachList) {
		char parenLeft = '(';
		char parenRight = ')';
		char separator = ';';
		
		listAsStr = listAsStr.replace(parenLeft, separator);
		listAsStr = listAsStr.replace(parenRight, separator);		
		String components[] = listAsStr.split("" + separator);
		ArrayList<double[]> finalList = new ArrayList<double[]>(components.length);
				
		for (int i = 0; i < components.length; i++) {
			String item = components[i].trim();
			if (item.equals("")) continue;
			double[] subList = getDoubleListFromStringForm(components[i], performSortInEachList);
			finalList.add(subList);
		}
		return finalList;
	}

	// ========================================================================
	/** Given a String list in braces {}, this returns the list of numbers contained
	 *  within those braces.  This also allows for ranges to be included.  For example,
	 *  an input of {1, 2, 3-5, 7} would return the list {1, 2, 3, 4, 5, 7}
	 */ 
	public static ArrayList<ArrayList<Integer>> getListOfIntListsFromStringForm(String listAsStr, boolean performSort) {
		char braceLeft  = '{';
		char braceRight = '}';
		char comma = ',';
		
		listAsStr = listAsStr.replace(braceLeft, comma);
		listAsStr = listAsStr.replace(braceRight, comma);		
		String components[] = listAsStr.split("" + comma);
		ArrayList<ArrayList<Integer>> finalList = new ArrayList<ArrayList<Integer>>(components.length);
		
		for (int i = 0; i < components.length; i++) {
			String item = components[i].trim();
			if (item.equals("")) continue;			
			ArrayList<Integer> subList = getIntListFromStringForm(components[i], performSort);
			finalList.add(subList);
		}
				
		return finalList;
	}
	
	// ========================================================================
	/** Given a String list of comma-separated numbers such that the head and tail
	 *  of the list are surrounded by braces, this returns the string without the
	 *  braces.  For example { 2, 3.5, 5, 10 } is returned as: 2, 3.5, 5, 10 
	 */
	public static String stripBraces(String listAsStr) {
		int indexBraceLeft  = listAsStr.indexOf('{');		
		int indexBraceRight = listAsStr.lastIndexOf('}');
		
		int substringStartIndex = (indexBraceLeft  < 0) ? 0 : indexBraceLeft + 1;
		int substringEndIndex   = (indexBraceRight < 0) ? listAsStr.length() : indexBraceRight;
		return listAsStr.substring(substringStartIndex, substringEndIndex);
	}
	
	
	// ========================================================================
	public interface ValueExtractor {
		public long extractValue(long compactUnit);		
	}
	
	public static ValueExtractor IntExtractorLSB = new ValueExtractor() { 
		public long extractValue(long compactUnit) { return (compactUnit & 0xFFFFFFFFL); }
	};
		
	public static ValueExtractor IntExtractorMSB = new ValueExtractor() { 
		public long extractValue(long compactUnit) { return (compactUnit >>> Integer.SIZE); }
	};
	

	// ========================================================================
	public static long searchMaxValue(LongArrayList longArray, ValueExtractor extractor) {
		long maxValue = Long.MIN_VALUE;
		int maxValueIndex = -1;
		
		for (int i = 0; i < longArray.size(); i++) {
			long value = extractor.extractValue(longArray.get(i));
			//System.out.println("" + i + "\t" + value);
			if (value > maxValue) {
				maxValue = value;
				maxValueIndex = i;
			}
		}
		
		return longArray.get(maxValueIndex);
	}
	
	// ========================================================================
	/** Given a sample index and integer array, this performs a binary search
	 *  for the sample index on the integer array.  If the sample index is 
	 *  found, an index >= 0 is returned, else (-(insertion point) - 1) is returned. */
	public static int binarySearchValue(final long value, LongArrayList longArray, ValueExtractor extractor) {
		int lowerIndex = 0;
		int upperIndex = longArray.size() - 1;
		int midIndex = 0;
		
		// If empty array, we return with the insertion point at index 0
		if (upperIndex < 0) return -1;
		
		// Shortcut Lower: Test if the value is <= than the lowest value the array.
		// If so, then we don't perform the binary search, and we simply return.
		long valueLowerIndex = extractor.extractValue(longArray.get(lowerIndex));
		if (value < valueLowerIndex) {
			return ((-lowerIndex) - 1);			
		} else if (value == valueLowerIndex) {
			return lowerIndex;
		} 

		// Shortcut Upper: Test if the value is >= the high value in the array.
		// If so, then we don't perform the binary search, and we simply return.
		long valueUpperIndex = extractor.extractValue(longArray.get(upperIndex));
		if (value > valueUpperIndex) {
			return (-(upperIndex + 1) - 1);
		} else if (value == valueUpperIndex) {
			return upperIndex;
		}
				
		long valueAtMidIndex;
		while (lowerIndex <= upperIndex) {
			midIndex = (lowerIndex + upperIndex) >>> 1;  // right-shift by 1 to divide by 2; 
			valueAtMidIndex = extractor.extractValue(longArray.get(midIndex));
			
			if (value == valueAtMidIndex) {
				return midIndex;
			} else if (value > valueAtMidIndex) {
				lowerIndex = midIndex + 1;
			} else {
				upperIndex = midIndex - 1;
			}
		}
		
		return -(lowerIndex + 1);
	}
}
