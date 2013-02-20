package nutils.counter;

import nutils.ArrayUtils;
import nutils.BitSetUtils;
import nutils.PrimitiveWrapper;
import nutils.NullaryClassFactory;

import com.carrotsearch.hppc.DoubleArrayList;
import com.carrotsearch.hppc.IntArrayList;
import com.carrotsearch.hppc.LongArrayList;

// ========================================================================
public class DynamicBucketCounter {

	public static NullaryClassFactory<DynamicBucketCounter> ClassFactory = new NullaryClassFactory<DynamicBucketCounter>(DynamicBucketCounter.class);
	
	LongArrayList mArray;
	private int mSumOfCounts;
	private int mTotalSum;
	private int mIndexOfMaxCount;
	
	public DynamicBucketCounter() {
		mArray = new LongArrayList();
		clear();
	}
	
	// ========================================================================
	public void clear() { 
		mArray.clear();
		mSumOfCounts = 0; 	
		mTotalSum = 0;
		mIndexOfMaxCount = -1;
	}
	
	// ========================================================================
	public int getCount(int key) {
		int indexOfKey = getIndexOfKey(key);
		if (indexOfKey < 0) return indexOfKey;
		return getCountAtIndex(indexOfKey);			
	}
	
	// ========================================================================
	public int getCountMax() {
		return (mIndexOfMaxCount < 0) ? 0 : getCountAtIndex(mIndexOfMaxCount);
	}
	
	// ========================================================================
	private int getCountAtIndex(int indexOfKey) {
		return (int) BitSetUtils.IntExtractorLSB.extractValue(mArray.get(indexOfKey));
	}
	
	// ========================================================================
	private int getKeyAtIndex(int indexOfKey) {
		return (int) BitSetUtils.IntExtractorMSB.extractValue(mArray.get(indexOfKey));
	}
	
	// ========================================================================
	private int getIndexOfKey(int key) {
		return ArrayUtils.binarySearchValue(key, mArray, BitSetUtils.IntExtractorMSB);
	}
	
	// ========================================================================
	public boolean getKeyLast(PrimitiveWrapper.WInteger key) {
		if (mArray.isEmpty()) return false;
		key.mInt = (int) BitSetUtils.IntExtractorMSB.extractValue(mArray.get(mArray.size() - 1));
		return true;
	}
	
	// ========================================================================
	public int incrementCount(int key) { return incrementCount(key, 1); }
	
	// ========================================================================
	public int incrementCount(int key, int numToIncrement) {
		mSumOfCounts += numToIncrement;
		
		int indexOfKey = getIndexOfKey(key);
		if (indexOfKey < 0) {
			int insertPoint = ArrayUtils.getInsertPoint(indexOfKey); 			
			mArray.insert(insertPoint, getKeyValueAsLong(key, numToIncrement));
			determineIndexOfMaxCountWithIncrementedIndex(insertPoint, numToIncrement);
			return numToIncrement;
		} else {
			int newCount = numToIncrement + getCountAtIndex(indexOfKey);
			mArray.set(indexOfKey, getKeyValueAsLong(key, newCount));		
			determineIndexOfMaxCountWithIncrementedIndex(indexOfKey, newCount);
			return newCount;
		}
	}

	// ========================================================================
	private void determineIndexOfMaxCountWithIncrementedIndex(int index, int count) {
		if ((mIndexOfMaxCount < 0) || (count > getCountAtIndex(mIndexOfMaxCount))) {
			mIndexOfMaxCount = index;
		}
	}
	
	// ========================================================================
	private static long getKeyValueAsLong(int key, int value) {
		return ( 0L | (((long) key) << Integer.SIZE) | value );
	}
	
	// ========================================================================
	public int getKeyWithMaxCount() {
		long keyValueMax = getKeyValueWithMaxCount();
		return (int) BitSetUtils.IntExtractorMSB.extractValue(keyValueMax);
	}
	
	// ========================================================================
	public long getKeyValueWithMaxCount() {
		return ArrayUtils.searchMaxValue(mArray, BitSetUtils.IntExtractorLSB);
	}
	
	// ========================================================================
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

	// ========================================================================
	public IntArrayList[] toArrayListInt() {
		IntArrayList[] rV = new IntArrayList[2];
		rV[0] = new IntArrayList();
		rV[1] = new IntArrayList();
		
		for (int i = 0; i < mArray.size(); i++) {				
			rV[0].add(getKeyAtIndex(i));
			rV[1].add(getCountAtIndex(i));
		}
		return rV;
	}
	
	// ========================================================================
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