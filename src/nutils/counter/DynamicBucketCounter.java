package nutils.counter;

import nutils.ArrayUtils;
import nutils.BitSetUtils;
import nutils.PrimitiveWrapper;
import nutils.NullaryClassFactory;

import com.carrotsearch.hppc.FloatArrayList;
import com.carrotsearch.hppc.DoubleArrayList;
import com.carrotsearch.hppc.IntArrayList;
import com.carrotsearch.hppc.LongArrayList;

// ========================================================================
public class DynamicBucketCounter {

	public static NullaryClassFactory<DynamicBucketCounter> ClassFactory = new NullaryClassFactory<DynamicBucketCounter>(DynamicBucketCounter.class);
	
	LongArrayList mArray;
	private int mSumOfCounts;
	private int mWeightedTotalSum;
	private int mIndexOfMaxCount;
	
	// ========================================================================
	public DynamicBucketCounter() {
		this(10);
	}
	
	// ========================================================================
	public DynamicBucketCounter(int initialCapacity) {
		mArray = new LongArrayList(initialCapacity);
		clear();		
	}
	
	// ========================================================================
	public void clear() { 
		mArray.clear();
		mSumOfCounts = 0; 	
		mWeightedTotalSum = 0;
		mIndexOfMaxCount = -1;
	}
	
	// ========================================================================
	public int getCount(int key) {
		int indexOfKey = getIndexOfKey(key);
		if (indexOfKey < 0) return indexOfKey;
		return getCountAtIndex(indexOfKey);			
	}
	
	// ========================================================================
	public void getCounts(IntArrayList countsBuffer) {
		countsBuffer.clear();
		
		int numKeys = getNumKeys();
		countsBuffer.ensureCapacity(numKeys);
		
		for (int i = 0; i < numKeys; i++) {
			countsBuffer.add(getCountAtIndex(i));
		}
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
	public int getNumKeys() { return mArray.size(); }
	
	// ========================================================================
	public int incrementCount(int key) { return incrementCount(key, 1); }
	
	// ========================================================================
	public int incrementCount(int key, int numToIncrement) {
		mSumOfCounts += numToIncrement;
		mWeightedTotalSum += (key * numToIncrement);
		
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
	public void getProportions(FloatArrayList proportionsBuffer) {
		proportionsBuffer.clear();
		
		int numKeys = getNumKeys();
		proportionsBuffer.ensureCapacity(numKeys);
		for (int i = 0; i < numKeys; i++) {
			float fraction = (float) getCountAtIndex(i) / (float) mSumOfCounts;
			proportionsBuffer.add(fraction);
		}
	}
	
	// ========================================================================
	public void getProportionsCumulativeForward(FloatArrayList proportionsBuffer) {
		getProportions(proportionsBuffer);
		
		// Return early if empty list		
		int numKeys = getNumKeys();
		if (numKeys == 0) return;
		
		for (int i = 1; i < numKeys; i++) {
			proportionsBuffer.set(i, proportionsBuffer.get(i) + proportionsBuffer.get(i - 1));			
		}		
	}
	
	// ========================================================================
	public void getCumulativeSumsForward(IntArrayList cumulativeCountsBuffer) {		
		cumulativeCountsBuffer.clear();		
		
		// Return early if empty list		
		int numKeys = getNumKeys();
		if (numKeys == 0) return;
		cumulativeCountsBuffer.ensureCapacity(numKeys);
		
		// First initialize first element
		cumulativeCountsBuffer.add(getCountAtIndex(0));
		
		for (int i = 1; i < numKeys; i++) {
			cumulativeCountsBuffer.add(cumulativeCountsBuffer.get(i - 1) + getCountAtIndex(i));			
		}
	}
	
	// ========================================================================
	public void getCumulativeSumsBackward(IntArrayList cumulativeCountsBuffer) {
		cumulativeCountsBuffer.clear();		
		
		// Return early if empty list		
		int numKeys = getNumKeys();
		if (numKeys == 0) return;
		cumulativeCountsBuffer.ensureCapacity(numKeys);
		
		// Initialize
		getCounts(cumulativeCountsBuffer);
		
		// Now go through and add
		for (int i = numKeys - 2; i >= 0; i--) {
			cumulativeCountsBuffer.set(i, cumulativeCountsBuffer.get(i) + cumulativeCountsBuffer.get(i + 1));			 
		 }
	}
	
	// ========================================================================
	public static void TestDynamicBucketCounter() {
		int numVals = 100;		
		DynamicBucketCounter dbc = new DynamicBucketCounter();
		
		for (int i = 1; i <= numVals; i++) {				
			for (int j = 0; j < i; j++) {
				dbc.incrementCount(i * 2);
			}
		}
		
		for (int i = 1; i <= numVals; i++) {
			System.out.printf("%d\t%d\n", i, dbc.getCount(i));
		}
		
		IntArrayList cumulativeForward = new IntArrayList();
		IntArrayList cumulativeBackward = new IntArrayList();
		FloatArrayList proportionsBuffer =  new FloatArrayList();
		FloatArrayList proportionsBufferCumulativeForward =  new FloatArrayList();
		dbc.getCumulativeSumsForward(cumulativeForward);
		dbc.getCumulativeSumsBackward(cumulativeBackward);
		dbc.getProportions(proportionsBuffer);
		dbc.getProportionsCumulativeForward(proportionsBufferCumulativeForward);
		
		for (int i = 0; i < dbc.getNumKeys(); i++) {
			System.out.printf("%d\t%d\t%d\t%d\t%d\t%g\t%g\n", i, dbc.getKeyAtIndex(i), dbc.getCountAtIndex(i), cumulativeForward.get(i), cumulativeBackward.get(i), proportionsBuffer.get(i), proportionsBufferCumulativeForward.get(i));
		}
		
		System.out.println("Key with max count: " + dbc.getKeyWithMaxCount());
		System.out.println("Sum of counts: " + dbc.mSumOfCounts);
		System.out.println("Total counts: " + dbc.mWeightedTotalSum);
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
	
	// ========================================================================
	public String toString() {
		StringBuilder sb = new StringBuilder(1024);
		int numElements = mArray.size();		
		for (int i = 0; i < numElements; i++) {
			if (i > 0) {
				sb.append(',');
			}
			sb.append(getKeyAtIndex(i));
			sb.append(':');
			sb.append(getCountAtIndex(i));
		}
		return sb.toString();
	}
	
	// ========================================================================
	public static void main(String[] args) {
		TestDynamicBucketCounter();
	}
	
}