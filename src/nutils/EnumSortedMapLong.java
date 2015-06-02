package nutils;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;

import nutils.BitUtils.ValueExtractor;
import nutils.array.FastArrayListLong.FastArrayListIndex;
import nutils.primitives.wrapper.PrimitiveWrapper;

import cern.colt.function.LongComparator;
import cern.colt.list.LongArrayList;

/**
 * This data structure is an EnumMap with a twist -- for each enum value, there
 * is a sorted dynamic array list of long primitives.  A special comparator must be
 * passed in in order to sort.  
 * @author Ninad Dewal
 *
 */
public class EnumSortedMapLong<E extends Enum<E>> {

	private static NullaryClassFactory<LongArrayList> LongArrayListClassFactory = new NullaryClassFactory<>(LongArrayList.class);
	
	protected EnumMapSafe<E, LongArrayList> mMap;
	protected EnumMapSafe<E, Boolean> mIsSorted;
	protected boolean mEntireMapSorted;  // used for speed optimization
	protected final Class<E> mEnumClass;
	protected ValueExtractor mExtractorForSorting;
	protected LongComparator mComparator;
	
	// ========================================================================
	public EnumSortedMapLong(Class<E> enumClass, ValueExtractor extractorForSorting) {
		this.mEnumClass = enumClass;		
		this.mExtractorForSorting = extractorForSorting;
		createMap(enumClass);		
		setComparator();		
	}
	
	// ========================================================================
	private void setComparator() {
		mComparator = new LongComparator() {			
			@Override
			public int compare(long arg0, long arg1) {
				return Long.compare(mExtractorForSorting.extractValue(arg0), mExtractorForSorting.extractValue(arg1));				
			}
		};
	}
			
	// ========================================================================
	private void createMap(Class<E> enumClass) {
		mMap = LongArrayListClassFactory.newEnumMap(enumClass);		
		mIsSorted = new EnumMapSafe<E, Boolean>(enumClass); 
		E[] enumValues = enumClass.getEnumConstants();

		for (E enumValue : enumValues) {			
			mIsSorted.put(enumValue, Boolean.TRUE);
		}
		mEntireMapSorted = true;
	}

	// ========================================================================
	public void clear(final E theEnum) {
		mMap.get(theEnum).clear();
	}
	
	// ========================================================================
	public void clear() {
		for (E enumValue : mEnumClass.getEnumConstants()) {
			mMap.get(enumValue).clear();
		}
	}
	
	// ========================================================================
	public int getIndex(final E theEnum, final long key) {
		ifNotSortedThenSort();
		return ArrayUtils.binarySearchValue(mExtractorForSorting.extractValue(key), mMap.get(theEnum), mExtractorForSorting);		
	}

	// ========================================================================
	public int size(final E theEnum) { return mMap.get(theEnum).size(); }
	
	// ========================================================================
	public int size() {
		int sum = 0;
		for (E theEnum : mEnumClass.getEnumConstants()) {
			sum += size(theEnum);
		}
		return sum;
	}	
	
	// ========================================================================
	public long get(final E theEnum, final int index) {		
		return mMap.get(theEnum).get(index);
	}
	
	// ========================================================================
	public boolean get(final E theEnum, final long key, PrimitiveWrapper.WLong unit) {		
		int resultIndex = getIndex(theEnum, key);
		if (resultIndex < 0) {
			return false;
		} else {
			unit.mLong = get(theEnum, resultIndex);
			return true;
		}  
	}

	// ========================================================================
	/** Replaces the value at the given index with a new value. */
	public boolean replace(final E theEnum, final int index, final long newValue) {
		LongArrayList targetList = mMap.get(theEnum); 				
		
		if (targetList.isEmpty()) {
			CompareUtils.ensureTrue(false, "Cannot set element on empty list!");
		}
		
		int nextIndex = index + 1;
		int prevIndex = index - 1;
		
		// Test only the relevant portion
		ValueExtractor extractor = mExtractorForSorting;
		final long value = extractor.extractValue(newValue);
		
		if (index == 0) {
			if ( (targetList.size() > 1) && (value > extractor.extractValue(targetList.get(nextIndex))) ) return false;  // return because value to be set will break sorted state	
			
		} else if (index == targetList.size() - 1) {
			if (value < extractor.extractValue(targetList.get(prevIndex))) return false;  // return because value to be set will break sorted state
				
		} else {
			if ((value > extractor.extractValue(targetList.get(nextIndex))) || 
				(value < extractor.extractValue(targetList.get(prevIndex)))) return false; 
		}
		
		// Now we are guaranteed that we will not break the ordering
		targetList.set(index, newValue);
		return true;
	}
	
	// ========================================================================
	/** Returns the index >= 0 of the unit if it already exists, else returns a negative index.
	 */
	public int addSorted(final E theEnum, final long theElement, PrimitiveWrapper.WLong value) {
		ifNotSortedThenSort();
		
		int resultIndex = getIndex(theEnum, theElement);
		if (resultIndex < 0) {						
			mMap.get(theEnum).beforeInsert(ArrayUtils.getInsertPoint(resultIndex), theElement);			
		} else {
			value.mLong = get(theEnum, resultIndex);			
		}
		return resultIndex;
	}

	// ========================================================================
	public void ensureCapacity(final E theEnum, int theSize) {
		mMap.get(theEnum).ensureCapacity(theSize);
	}
	
	// ========================================================================
	/** Adds an object at the specified enum.  This does NOT check whether the 
	 *  element already exists.  The purpose of this function is to be faster than
	 *  the other add function. 
	 */
	public boolean addToTail(final E theEnum, final long theElement) {
		LongArrayList theList = mMap.get(theEnum); 
		theList.add(theElement);
		
		// Check if we need to set the unsorted flag
		if (theList.size() > 1) {
			long nextToLastElement = theList.get(theList.size() - 2);
			//System.out.println(theElement + "\t" + nextToLastElement + "\t" + mComparator.compare(theElement, nextToLastElement));
			if (mComparator.compare(theElement, nextToLastElement) < 0) {					
				setUnsorted(theEnum);
			}
		}

		return true;
	}
	
	// ========================================================================
	protected void ifNotSortedThenSort() {
		if (mEntireMapSorted) return;  // As a speed up, check if any changes to overall map.
		
		for (E enumValue : mEnumClass.getEnumConstants()) {
			if (!mIsSorted.get(enumValue).booleanValue()) {
				LongArrayList theList = mMap.get(enumValue); 
				theList.quickSortFromTo(0, theList.size() - 1, mComparator);
				mIsSorted.put(enumValue, Boolean.TRUE);
				
				System.out.println("SORTING! " + enumValue);
			}
		}		
		
		
		mEntireMapSorted = true; // Set that the map is sorted now
	}

	// ========================================================================
	public void sortTable() { ifNotSortedThenSort(); }
	
	// ========================================================================
	private void setUnsorted(E theEnum) {
		mEntireMapSorted = false;
		mIsSorted.put(theEnum, Boolean.FALSE);
	}
	
	// ========================================================================
	public void print(PrintStream outStream) {
		ifNotSortedThenSort();
		for (E et : mEnumClass.getEnumConstants()) {
			LongArrayList theList = mMap.get(et);
			for (int i = 0; i < theList.size(); i++) {
				outStream.printf("%s\t%d\t%d\n", et, mExtractorForSorting.extractValue(theList.get(i)), theList.get(i));
			}
		}
	}

	// ========================================================================
	private static enum TestEnum { A, B, C, D, E; }
	
	private static void Test() {
		EnumSortedMapLong<TestEnum> testMap = new EnumSortedMapLong<>(TestEnum.class, ValueExtractor.IntExtractorMSB);
		int numIter = 100000;
		
		for (TestEnum te : TestEnum.values()) {
			
			for (int i = numIter; i >= 0; i--) {
				long key = ((long) i) << Integer.SIZE;
				testMap.addToTail(te, key);
			}
		}
		//testMap.print(System.out);
		
		PrimitiveWrapper.WLong longHolder = new PrimitiveWrapper.WLong(0); 
		for (TestEnum te : TestEnum.values()) {
			for (int i = 0; i <= numIter; i++) {
				long key = ((long) i) << Integer.SIZE;
				boolean result = testMap.get(te, key, longHolder);
				CompareUtils.ensureTrue(result, "ERROR: Could not find inserted value!");
			}
		}
		
		// Now add More Elements
		for (TestEnum te : TestEnum.values()) {
			for (int i = numIter + numIter; i > numIter; i--) {
				long key = ((long) i) << Integer.SIZE;
				int resultIndex = testMap.addSorted(te, key, longHolder);
				if (resultIndex >= 0) {
					CompareUtils.ensureTrue(false, "ERROR: Should not be in list!\t" + i + "\t" + te);
				}
			}
		}
		//testMap.print(System.out);
		
		for (TestEnum te : TestEnum.values()) {
			for (int i = 0; i <= numIter + numIter; i++) {
				long key = ((long) i) << Integer.SIZE;
				boolean result = testMap.get(te, key, longHolder);
				CompareUtils.ensureTrue(result, "ERROR: Could not find inserted value!");
			}
		}
		
		
	}
	
	// ========================================================================
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		Test();
	}

}
