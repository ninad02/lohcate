package nutils;

import nutils.BitUtils.ValueExtractor;
import nutils.array.FastArrayListLong;
import nutils.array.FastArrayListLong.FastArrayListIndex;
import nutils.primitives.wrapper.PrimitiveWrapper;

public class EnumSortedMapLongFast<E extends Enum<E>> {
	
	private static final NullaryClassFactory<FastArrayListLong> LongArrayListClassFactory = new NullaryClassFactory<>(FastArrayListLong.class);
	
	protected EnumMapSafe<E, FastArrayListLong> mMap;
	protected final Class<E> mEnumClass;
	protected ValueExtractor mExtractorForSorting;
	private FastArrayListIndex mDummyIndex; 	
	
	// ========================================================================
	public EnumSortedMapLongFast(Class<E> enumClass, ValueExtractor extractorForSorting) {
		this.mEnumClass = enumClass;		
		this.mExtractorForSorting = extractorForSorting;
		this.mDummyIndex = FastArrayListLong.getNewIndexToken();
		createMap(enumClass);				
	}
	
	// ========================================================================
	private void createMap(Class<E> enumClass) {
		mMap = LongArrayListClassFactory.newEnumMap(enumClass);				
		E[] enumValues = enumClass.getEnumConstants();

		// Set the long-extractor for each element
		for (E enumValue : enumValues) {	
			mMap.get(enumValue).setValueExtractor(mExtractorForSorting);			
		}
	}

	// ========================================================================
	public void clear(final E theEnum) {
		mMap.get(theEnum).clear();
	}
	
	// ========================================================================
	public void clear() {
		for (E enumValue : mEnumClass.getEnumConstants()) {
			this.clear(enumValue);			
		}
	}
	
	// ========================================================================
	public FastArrayListIndex getIndex(final E theEnum, final long key, FastArrayListIndex theIndex) {		
		return mMap.get(theEnum).getIndex(key, theIndex);
	}
	
	// ========================================================================
	public long get(final E theEnum, final FastArrayListIndex index) {		
		return mMap.get(theEnum).get(index);
	}
	
	// ========================================================================
	public synchronized boolean get(final E theEnum, final long key, PrimitiveWrapper.WLong unit) {

		FastArrayListIndex resultIndex = getIndex(theEnum, key, mDummyIndex);
		if (resultIndex == null) {
			return false;
		} else {
			unit.mLong = get(theEnum, resultIndex);
			return true;
		}  
	}
	
	// ========================================================================
	/** Replaces the value at the given index with a new value. */
	public synchronized boolean replace(final E theEnum, final long newValue) {
		FastArrayListIndex resultIndex = getIndex(theEnum, newValue, mDummyIndex);
		if (resultIndex == null) {
			return false;
		} else {
			mMap.get(theEnum).replace(resultIndex, newValue);
			return true;
		}
	}
	
	// ========================================================================
	/** Returns the index >= 0 of the unit if it already exists, else returns a negative index.
	 */
	public synchronized boolean addIfUnique(final E theEnum, final long theElement, PrimitiveWrapper.WLong value) {
		return mMap.get(theEnum).addIfUnique(theElement, value);		
	}
	
	// ========================================================================
	/** Adds an object at the specified enum.  This does NOT check whether the 
	 *  element already exists.  The purpose of this function is to be faster than
	 *  the other add function. 
	 */
	public boolean addToTail(final E theEnum, final long theElement) {
		mMap.get(theEnum).add(theElement);
		return true;
	}
	
	// ========================================================================
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub

	}

}
