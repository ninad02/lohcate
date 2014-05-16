package nutils;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.Random;

import nutils.BitUtils.ValueExtractor;
import nutils.counter.BucketCounterEnum;

import com.carrotsearch.hppc.DoubleArrayList;
import com.carrotsearch.hppc.LongArrayList;

public class ArrayUtils {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		//Test_removeNullElements();
		Test_removeNullElementsSimple();
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
		for (int i = 0; i < (theArray.length >>> 1); i++) {
			oppositeIndex = theArray.length - i - 1;
			tempVal = theArray[i];
			theArray[i] = theArray[oppositeIndex];
			theArray[oppositeIndex] = tempVal;
		}
	}
	
	public static void reverseArray(int[] theArray) {
		int tempVal;
		int oppositeIndex;
		for (int i = 0; i < theArray.length >>> 1; i++) {
			oppositeIndex = theArray.length - i - 1;
			tempVal = theArray[i];
			theArray[i] = theArray[oppositeIndex];
			theArray[oppositeIndex] = tempVal;
		}
	}

	public static void reverseArray(double[] theArray) {
		double tempVal;
		int oppositeIndex;
		for (int i = 0; i < theArray.length >>> 1; i++) {
			oppositeIndex = theArray.length - i - 1;
			tempVal = theArray[i];
			theArray[i] = theArray[oppositeIndex];
			theArray[oppositeIndex] = tempVal;
		}
	}

	public static void reverseArray(char[] theArray) {
		char tempVal;
		int oppositeIndex;
		for (int i = 0; i < theArray.length >>> 1; i++) {
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
		for (int i = 0; i < theArray.length >>> 1; i++) {
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
	public static int[] getPrimitiveArrayInt(Collection<Integer> theList) {
		int[] rV = new int[theList.size()];
		int index = 0;
		for (Integer num : theList) {
			rV[index++] = num.intValue();
		}		
		return rV;
	}

	/** Converts a collection of Double to a primitive list. */
	public static double[] getPrimitiveArrayDouble(Collection<Double> theList) {
		double[] rV = new double[theList.size()];
		int index = 0;
		for (Double d : theList) {		
			rV[index++] = d.doubleValue();
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
	
	/** Given an array, this returns the sum of the array elements. */
	public static double arraySum(double[] theArray) {
		double sum = 0;
		for (int i = 0; i < theArray.length; i++) { sum += theArray[i]; }
		return sum;
	}
	
	/** Given a primitive array list, this sorts it. */
	public static void sort(com.carrotsearch.hppc.DoubleArrayList theArray) { Arrays.sort(theArray.buffer, 0, theArray.size()); }
	public static void sort(com.carrotsearch.hppc.IntArrayList theArray)    { Arrays.sort(theArray.buffer, 0, theArray.size()); }
	public static void sort(com.carrotsearch.hppc.LongArrayList theArray)   { Arrays.sort(theArray.buffer, 0, theArray.size()); }
	public static void sort(com.carrotsearch.hppc.ByteArrayList theArray)   { Arrays.sort(theArray.buffer, 0, theArray.size()); }
	public static void sort(com.carrotsearch.hppc.ShortArrayList theArray)  { Arrays.sort(theArray.buffer, 0, theArray.size()); }
	public static void sort(com.carrotsearch.hppc.FloatArrayList theArray)  { Arrays.sort(theArray.buffer, 0, theArray.size()); }


	
	/** Given an array of doubles, this returns the index of the maximum element.  If there
	 *  are multiple maximal elements, this returns the index of the first one.  The second
	 *  argument is the index at which searching should begin. */
	public static int getIndexOfMaxElement(double[] theArray, int indexToStartSearch) {
		int maxIndex = -1;
		double maxElement = -Double.MAX_VALUE; 
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
	
	public static int getIndexOfMinElement(int[] theArray, int indexToStartSearch) {
		int minIndex = -1;
		int minElement = Integer.MAX_VALUE; 
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
	 */
	public static <T, C extends Collection<T>> C addToCollection(C c, T o, int count, boolean clearCollection) {
		if (clearCollection) {
			c.clear();
		}
				
		for (int i = 0; i < count; i++) {
			c.add(o);
		}
		
		return c;
	}

	// ========================================================================
	// ===== INNER CLASS =====
	// ========================================================================
	public static abstract class NumArray {
		public int mArraySize;
		
		public NumArray() { 
			reset(); 
		}		
		
		public void reset() { mArraySize = 0; }
		public int size() { return mArraySize; }
	}

	// ========================================================================
	public static class IntArray extends NumArray {
		public int[] mArray;
					
		public IntArray(int numElements) { 
			super();
			mArray = new int[numElements];
		}
		
		public int add(int num) { 
			mArray[mArraySize] = num; 
			return ++mArraySize; 
		}
	}

	// ========================================================================
	public static class DoubleArray extends NumArray {
		public double[] mArray;
		
		public DoubleArray(int numElements) { 
			super();
			mArray = new double[numElements];
		}
		
		public int add(double num) { 
			mArray[mArraySize] = num; 
			return ++mArraySize; 
		}
	}
	
	// ========================================================================
	public static class ParallelArrayDouble extends NumArray {
		public double[][] mArray;		
		
		public ParallelArrayDouble(int numElements) { 
			super();
			mArray = new double[2][numElements];
			ArrayUtils.arrayFill(mArray, 0);
		}
		
		public int add(double numX, double numY) {
			mArray[0][mArraySize] = numX;
			mArray[1][mArraySize] = numY;
			return  ++mArraySize; 
		}
		
		public double[][] toArrays() {
			return mArray;
		}
	}
	
	// ========================================================================
	public static class ParallelArrayDoubleDynamic {
		public DoubleArrayList[] mLists;
		
		public ParallelArrayDoubleDynamic() {
			mLists = new DoubleArrayList[] { new DoubleArrayList(), new DoubleArrayList() } ;			
		}
		
		public ParallelArrayDoubleDynamic(int numEstimatedElements) {
			mLists = new DoubleArrayList[] { new DoubleArrayList(numEstimatedElements), new DoubleArrayList(numEstimatedElements) } ;
		}
		
		public int add(double numX, double numY) {
			mLists[0].add(numX);
			mLists[1].add(numY);
			return mLists[0].size();
		}
		
		public double[][] toArrays() {
			double[][] rV = new double[2][];
			for (int i = 0; i < mLists.length; i++) {
				rV[i] = mLists[i].toArray();
			}
			return rV;
		}
	}
	
	
	// ========================================================================
	// ===== INNER CLASS =====
	// ========================================================================
	
	// ========================================================================
	public static int getInsertPoint(int negativeArrayIndex) { 
		return ((negativeArrayIndex < 0) ? -(negativeArrayIndex + 1) : negativeArrayIndex); 
	}
	
	// ========================================================================
	public static int getFailPoint(int arrayIndex) {
		if (arrayIndex < 0) {
			CompareUtils.ensureTrue(false, "ERROR: Failing point index must be >= 0 instead of: " + arrayIndex);
		}
		return (-arrayIndex - 1);
	}
	
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
	public static DoubleArrayList getDoubleListFromStringForm(String listAsStr, boolean performSort) {	
		listAsStr = stripBraces(listAsStr); 				
		String components[] = listAsStr.split(",");
		DoubleArrayList rV = new DoubleArrayList(components.length);
		for (String component : components) {
			component = component.trim();			
			if (component.length() > 0) {				
				rV.add(Double.parseDouble(component.trim()));				
			}
		}				
		if (performSort) { ArrayUtils.sort(rV); }
		return rV;
	}

	// ========================================================================
	/** Given a String list in the format: (list1;list2;...;list_n), where each list
	 *  is of the form {3,2.6,9} (for example), this returns a list of lists represented
	 *  by the string
	 */
	public static ArrayList<DoubleArrayList> getListOfDoubleListsFromStringForm(String listAsStr, boolean performSortInEachList) {
		char parenLeft = '(';
		char parenRight = ')';
		char separator = ';';
		
		listAsStr = listAsStr.replace(parenLeft, separator);
		listAsStr = listAsStr.replace(parenRight, separator);		
		String components[] = listAsStr.split("" + separator);
		ArrayList<DoubleArrayList> finalList = new ArrayList<DoubleArrayList>(components.length);
				
		for (int i = 0; i < components.length; i++) {
			String item = components[i].trim();
			if (item.equals("")) continue;			
			DoubleArrayList subList = getDoubleListFromStringForm(components[i], performSortInEachList);
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
	
	
	// ========================================================================
	/** Given a sample index and integer array, this performs a binary search
	 *  for the sample index on the integer array.  If the sample index is 
	 *  found, an index >= 0 is returned, else (-(insertion point) - 1) is returned. */
	public static int binarySearchValue(final long value, com.carrotsearch.hppc.LongArrayList longArray, ValueExtractor extractor) {
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

	// ========================================================================
	/** Given a sample index and integer array, this performs a binary search
	 *  for the sample index on the integer array.  If the sample index is 
	 *  found, an index >= 0 is returned, else (-(insertion point) - 1) is returned. */
	public static int binarySearchValue(final long value, cern.colt.list.LongArrayList longArray, ValueExtractor extractor) {
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
	
	// ========================================================================
	/** Given a target collection, this adds elements from the source collection at the end
	 *  of the target collection.  We create this function because it is ironically more 
	 *  efficient than the addAll() method in some collections, which stupidly allocates
	 *  memory to create an extra and needless temporary array in its implementation.
	 * @param listToWhichToAdd
	 * @param elementsToAdd
	 */
	// We create this function because it is ironically more efficient than the ArrayList.addAll() method,
	// which stupidly allocates memory to create an extra and needless temporary array in its implementaion.
	public static<T> void addAll(ArrayList<T> listToWhichToAdd, ArrayList<T> elementsToAdd) {
		listToWhichToAdd.ensureCapacity( listToWhichToAdd.size() + elementsToAdd.size() );  // Ensure capacity so no memory reallocation later		
		for (T element : elementsToAdd) {
			listToWhichToAdd.add(element);
		}
	}

	// ========================================================================
	/** Given a two dimensional double array, this fills all elements with the specified value. */
	public static void arrayFill(double[][] theArray, double value) {
		for (double[] subArray : theArray) {
			Arrays.fill(subArray, value);
		}
	}	
	
	// ========================================================================
	/** Given a two dimensional object array, this fills all elements with the specified value. */
	public static<E> void arrayFill(E[][] theMatrix, E value) {
		for (E[] theArray : theMatrix) {
			Arrays.fill(theArray, value);
		}
	}
	
	// ========================================================================
	/** Given a two dimensional integer array, this fills all elements with the specified value. */
	public static void arrayFill(int[][] theArray, int value) {
		for (int[] subArray : theArray) {
			Arrays.fill(subArray, value);
		}
	}

	// ========================================================================
	/** Given a three dimensional integer array, this fills all the elements with the specified value. */
	public static void arrayFill(int[][][] theArray, int value) {
		for (int[][] subArray : theArray) {
			arrayFill(subArray, value);
		}
	}

	// ========================================================================
	/** Given a four dimensional integer array, thsi fills all the elements with the specified value. */
	public static void arrayFill(int[][][][] theArray, int value) {
		for (int[][][] subArray : theArray) {
			arrayFill(subArray, value);
		}
	}

	// ========================================================================
	/** This returns a new int[] array with the specified length and initialized to the specified values. */
	public static int[] newIntArray(int numElements, int fillValue) {
		int[] newArray = new int[numElements];
		Arrays.fill(newArray, fillValue);
		return newArray;
	}

	// ========================================================================
	/** This returns a new int[][] array with the specified dimensions and with all elements
	 *  initialized to the specified value.
	 */
	public static int[][] newIntArray2D(int numElementsX, int numElementsY, int fillValue) {
		int[][] newArray = new int[numElementsX][numElementsY];
		arrayFill(newArray, fillValue);
		return newArray;
	}

	// ========================================================================
	/** This returns a new double[][] array with the specified dimensions and with all elements
	 *  initialized to the specified value.
	 */
	public static double[][] newDoubleArray2D(int numElementsX, int numElementsY, double fillValue) {
		double[][] newArray = new double[numElementsX][numElementsY];
		arrayFill(newArray, fillValue);
		return newArray;
	}
	
	// ========================================================================
	/** Given an arraylist of arraylists of a certain type, this adds numNewArraysToAdd new ArrayLists of
	 *  that type to the target list.
	 * @param targetList
	 * @param numNewArraysToAdd
	 */
	public static<T> ArrayList<ArrayList<T>> addNewEmptyArrayLists(ArrayList<ArrayList<T>> targetList, int numNewArraysToAdd) {
		targetList = (targetList == null) ? new ArrayList<ArrayList<T>>() : targetList;
		
		for (int i = 0; i < numNewArraysToAdd; i++) {
			targetList.add(new ArrayList<T>());
		}
		
		return targetList;
	}

	// ========================================================================
	// removes last n elements from the list
	public static<T> ArrayList<T> removeLastNElements(ArrayList<T> theList, int numLastElementsToRemove) {
		numLastElementsToRemove = Math.min(numLastElementsToRemove, theList.size());
		int listSize = theList.size();
		for (int i = 0; i < numLastElementsToRemove; i++) {
			theList.remove(listSize - 1);
			listSize--;
		}
		return theList;
	}

	// ========================================================================
	// removes any elements set to null in the list
	public static<T> ArrayList<T> removeNullElements(ArrayList<T> theList) {
		
		int indexNow = 0;
		int indexAhead = indexNow + 1;
		int numNullValues = 0;
		
		boolean remainingElementsNonNull = true;
		int size = theList.size();
		for ( ; (indexNow < size) && remainingElementsNonNull; indexNow++) {
			if (theList.get(indexNow) == null) {
				// Currently we have a null element.
				numNullValues++;
				
				// Now we move indexAhead forward until we 
				// reach a non-null value or the end of the list				
				for (indexAhead = Math.max(indexNow + 1, indexAhead); 
						(indexAhead < size) && (theList.get(indexAhead) == null); 
						indexAhead++) {}				
				
				// At this point, indexAhead has either crossed its bounds
				// or has found an element
				if (indexAhead >= size) {    
					// crossed bounds.  Means that there were no non-null elements past 
					// the current index.  All we do now is remove the current and 
					// remaining (null) elements from the list.
					int numLastElementsToRemove = size - indexNow;  					
					removeLastNElements(theList, numLastElementsToRemove);
					numNullValues += numLastElementsToRemove - 1;
					remainingElementsNonNull = false;
				} else {  
					// found a non-null element, so swap them
					theList.set(indexNow, theList.get(indexAhead));
					theList.set(indexAhead, null);
				}
				
			} 
		}
				
		return theList;		
	}

	// ========================================================================
	private static void Test_removeNullElementsSimple() {
		ArrayList<Integer> list1 = new ArrayList<Integer>();
		for (int i = 0; i < 15; i++) { list1.add(i); }
		list1.set(14, null);
		list1.set(13, null);
		list1.set(12, null);
		list1.set(11, null);
		
		removeNullElements(list1);
		System.out.println(list1);
	}
	
	// ========================================================================
	private static void Test_removeNullElements() {
		
		ArrayList<Integer> list1 = new ArrayList<Integer>();
		ArrayList<Integer> list2 = new ArrayList<Integer>();
		ArrayList<Integer> listM = new ArrayList<Integer>();  // a master list
		
		ArrayList<Integer> list2IndicesToRemove = new ArrayList<Integer>();
		
		int numTrials = 10000;
		int numElementsInArray = 100;
		int numElementsToDelete = 30;
		Random randomGen = new Random();
		
		for (int i = 0; i < numElementsInArray; i++) {			
			listM.add(new Integer(i));
		}
		
		for (int trial = 0; trial < numTrials; trial++) {
			list1.clear();
			list2.clear();
			list1.addAll(listM);
			list2.addAll(listM);
			list2IndicesToRemove.clear();
			
			// Now select random elements to delete from list1 
			for (int j = 0; j < numElementsToDelete; ) {
				int indexToDelete = randomGen.nextInt(numElementsInArray);
				if (list1.get(indexToDelete) != null) {
					list1.set(indexToDelete, null);
					list2IndicesToRemove.add(indexToDelete);
					j++;						
				}
			}
			
			// Now we must delete from list2, though in reverse order
			Collections.sort(list2IndicesToRemove);
			Collections.reverse(list2IndicesToRemove);
			for (Integer theIndex : list2IndicesToRemove) {
				list2.remove(theIndex);
			}
			
			// Now perform the deletion on the list 1
			removeNullElements(list1);
	
			// Now compare
			boolean areListsEqual = list1.equals(list2);
			if (areListsEqual) {
				//System.out.println(list1);
				//System.out.println(list2);
				//System.out.println("");
				//break;
			} else {
				System.out.println("Lists Unequal:");
				System.out.println(list1);
				System.out.println(list2);
				System.out.println("");

			}
		}
	}

	// ========================================================================
	/** This assumes that the list is sorted.  It tests whether the element already exists in
	 *  the list and inserts the element only if it doesn't already exist in the list. 
	 * @param targetList The list that is tested for the presence of the element in the list  
	 * @param element The element to be added if it doesn't already exist in the list
	 * @return true if the element already existed in the list, false otherwise
	 */
	public static<T extends Comparable<? super T>> boolean checkInListAndInsertIfMissing(List<T> list, T key) {
		int resultIndex = Collections.binarySearch(list, key);
		if (resultIndex < 0) {
			int insertionIndex = getInsertPoint(resultIndex);
			list.add(insertionIndex, key);			
			return false;
		} else {
			return true;
		}		
	}

	// ========================================================================
	/** This linearly searches the aray for the integer value.  
	 * 
	 * @param theArray The input array
	 * @param targetValue The target value
	 * @return The index of the target value in the input array, or -1 if not found.
	 */
	public static int linearSearch(ArrayList<Integer> theArray, final int targetValue) {
		for (Integer i : theArray) {
			if (targetValue == i) return i;
		}
		return -1;
	}
	
	// ========================================================================
	/** Creates an arrayList of new objects given the desired array size and class factory. */
	
	
	// ========================================================================
	/** Given a list of enums, this returns the counts of those enum. */ 
	public static<T extends Enum<T>> BucketCounterEnum<T> getEnumTypeCounts(Collection<T> enumTypeList, Class<T> enumClass) {		
		BucketCounterEnum<T> enumCounts = new BucketCounterEnum<T>(enumClass);				
		for (T enumType : enumTypeList) {
			enumCounts.increment(enumType);			
		}		
		return enumCounts;
	}
	
	// ========================================================================
	/** This linearly searches the aray for the integer value.  
	 * 
	 * @param theArray The input array
	 * @param elem The target object
	 * @return The index of the target value in the input array, or -1 if not found.
	 */	
	public static<E> int linearSearch(ArrayList<E> theArray, E elem) {
		return theArray.indexOf(elem);
	}

	// ========================================================================
	public static double[][] combineTwoDynamicArraysIntoOneStatic(DoubleArrayList list1, DoubleArrayList list2) {
		return new double[][] { list1.toArray(), list2.toArray() };
	}
}
