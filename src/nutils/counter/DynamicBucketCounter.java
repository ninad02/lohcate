package nutils.counter;

import nutils.ArrayUtils;
import nutils.BitSetUtils;

import com.carrotsearch.hppc.DoubleArrayList;
import com.carrotsearch.hppc.LongArrayList;

// ========================================================================
public class DynamicBucketCounter {

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
		return (int) BitSetUtils.IntExtractorLSB.extractValue(mArray.get(indexOfKey));
	}
	
	private int getKeyAtIndex(int indexOfKey) {
		return (int) BitSetUtils.IntExtractorMSB.extractValue(mArray.get(indexOfKey));
	}
	
	private int getIndexOfKey(int key) {
		return ArrayUtils.binarySearchValue(key, mArray, BitSetUtils.IntExtractorMSB);
	}
	
	public int getKeyLast() { return (int) BitSetUtils.IntExtractorMSB.extractValue(mArray.get(mArray.size() - 1)); }
					
	public int incrementCount(int key) { return incrementCount(key, 1); }
	
	public int incrementCount(int key, int numToIncrement) {
		int indexOfKey = getIndexOfKey(key);
		if (indexOfKey < 0) {
			mArray.insert(ArrayUtils.getInsertPoint(indexOfKey), getKeyValueAsLong(key, numToIncrement));
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
		return (int) BitSetUtils.IntExtractorMSB.extractValue(keyValueMax);
	}
	
	public long getKeyValueWithMaxCount() {
		return ArrayUtils.searchMaxValue(mArray, BitSetUtils.IntExtractorLSB);
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