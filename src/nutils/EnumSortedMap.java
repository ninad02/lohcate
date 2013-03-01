package nutils;

import java.util.ArrayList;
import java.util.Collections;

/**
 * This data structure is an EnumMap with a twist -- for each enum value, there
 * is a sorted ArrayList of element objects.  The element objects themselves
 * serve as their own keys via the compareTo() function.  The element object V
 * must also either: support a null constructor, or a dummy object for use by the 
 * class must be passed in.
 *  
 * @author Ninad Dewal
 *
 */
public class EnumSortedMap<E extends Enum<E>, V extends Comparable<V>> {
	
	protected EnumMapSafe<E, ArrayList<V>> mMap;
	protected EnumMapSafe<E, Boolean> mIsSorted;
	protected boolean mEntireMapSorted;  // used for speed optimization
	protected final Class<E> mEnumClass;
	
	
	// ========================================================================
	public EnumSortedMap(Class<E> enumClass) {
		mEnumClass = enumClass;		
		createMap(enumClass);
	}
			
	// ========================================================================
	private void createMap(Class<E> enumClass) {
		mMap = new EnumMapSafe<E, ArrayList<V>>(enumClass);
		mIsSorted = new EnumMapSafe<E, Boolean>(enumClass); 
		E[] enumValues = enumClass.getEnumConstants();

		for (E enumValue : enumValues) {
			mMap.put(enumValue, new ArrayList<V>());
			mIsSorted.put(enumValue, Boolean.TRUE);
		}
		mEntireMapSorted = true;
	}
	
	// ========================================================================
	public int getIndex(final E theEnum, final V key) {
		ifNotSortedThenSort();
		return Collections.binarySearch(mMap.get(theEnum), key);
	}
	
	// ========================================================================
	public V get(final E theEnum, final int index) {		
		return mMap.get(theEnum).get(index);
	}
	
	// ========================================================================
	public V get(final E theEnum, final V key) {		
		int resultIndex = getIndex(theEnum, key);
		return ((resultIndex < 0) ? null : get(theEnum, resultIndex)); 
	}
	
	// ========================================================================
	/** Returns the object if it already exists, or null if it didn't (and it adds this one).
	 */
	public V addSorted(final E theEnum, final V theElement) {
		ifNotSortedThenSort();
		
		int resultIndex = getIndex(theEnum, theElement);
		if (resultIndex < 0) {						
			mMap.get(theEnum).add(ArrayUtils.getInsertPoint(resultIndex), theElement);
			return null;
		} else {
			return mMap.get(theEnum).get(resultIndex);
		}
	}
	
	// ========================================================================
	/** Adds an object at the specified enum.  This does NOT check whether the 
	 *  element already exists.  The purpose of this function is to be faster than
	 *  the other add function. 
	 */
	public boolean addToTail(final E theEnum, final V theElement) {
		ArrayList<V> theList = mMap.get(theEnum); 
		boolean result = theList.add(theElement);
		
		// Check if we need to set the unsorted flag
		if (result && theList.size() > 1) {
			V nextToLastElement = theList.get(theList.size() - 2);
			if (theElement.compareTo(nextToLastElement) < 0) {
				setUnsorted(theEnum);
			}
		}

		return result;
	}

	// ========================================================================
	protected void ifNotSortedThenSort() {
		if (mEntireMapSorted) return;  // As a speed up, check if any changes to overall map.
		
		for (E enumValue :  mEnumClass.getEnumConstants()) {
			if (!mIsSorted.get(enumValue).booleanValue()) {
				Collections.sort(mMap.get(enumValue));
				mIsSorted.put(enumValue, Boolean.TRUE);
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
	
	// ========================================================================

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub

	}

}
