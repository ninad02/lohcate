package nutils.counter;

import genomeUtils.ChromPosition;

import java.util.HashMap;

import nutils.ArrayUtils;
import nutils.NumberUtils;
import nutils.PrimitiveWrapper;
import nutils.NullaryClassFactory;
import nutils.BitUtils.ValueExtractor;

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
		if (indexOfKey < 0) return 0;
		return getCountAtIndex(indexOfKey);			
	}
	
	// ========================================================================
	public IntArrayList getCounts(IntArrayList countsBuffer) {
		countsBuffer = (countsBuffer == null) ? new IntArrayList() : countsBuffer; 
		countsBuffer.clear();
		
		int numKeys = getNumKeys();
		countsBuffer.ensureCapacity(numKeys);
		
		for (int i = 0; i < numKeys; i++) {
			countsBuffer.add(getCountAtIndex(i));
		}
		return countsBuffer;
	}
	
	// ========================================================================
	public int getCountMax() {			
		return (mIndexOfMaxCount < 0) ? 0 : getCountAtIndex(mIndexOfMaxCount);
	}
	
	// ========================================================================
	private int getCountAtIndex(int indexOfKey) {
		return (int) ValueExtractor.IntExtractorLSB.extractValue(mArray.get(indexOfKey));
	}
	
	// ========================================================================
	private int getKeyAtIndex(int indexOfKey) {
		return (int) ValueExtractor.IntExtractorMSB.extractValue(mArray.get(indexOfKey));
	}
	
	// ========================================================================
	private int getIndexOfKey(int key) {
		return ArrayUtils.binarySearchValue(key, mArray, ValueExtractor.IntExtractorMSB);
	}
	
	// ========================================================================
	private float calcProportionAtIndex(int indexOfKey) {
		return (float) getCountAtIndex(indexOfKey) / (float) mSumOfCounts;
	}
	
	// ========================================================================
	public boolean getKeyLast(PrimitiveWrapper.WInteger key) {
		if (mArray.isEmpty()) return false;
		key.mInt = (int) ValueExtractor.IntExtractorMSB.extractValue(mArray.get(mArray.size() - 1));
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
			if (mIndexOfMaxCount >= insertPoint) { ++mIndexOfMaxCount; }
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
		int index = getIndexOfKeyWithMaxCount();
		return (int) ValueExtractor.IntExtractorMSB.extractValue(mArray.get(index));
	}
	
	// ========================================================================
	public int getIndexOfKeyWithMaxCount() {	
		return ArrayUtils.searchIndexOfMaxValue(mArray, ValueExtractor.IntExtractorLSB);		
	}

	// ========================================================================
	public IntArrayList getKeys(IntArrayList keysBuffer) {
		keysBuffer = (keysBuffer == null) ? new IntArrayList() : keysBuffer;
		
		keysBuffer.clear();
		int numKeys = getNumKeys();
		keysBuffer.ensureCapacity(numKeys);
		for (int i = 0; i < numKeys; i++) {
			keysBuffer.add(getKeyAtIndex(i));
		}
		return keysBuffer;
	}
	
	// ========================================================================
	public FloatArrayList getProportions(FloatArrayList proportionsBuffer) {
		
		proportionsBuffer = (proportionsBuffer == null) ? new FloatArrayList() : proportionsBuffer;		
		proportionsBuffer.clear();
		
		int numKeys = getNumKeys();
		proportionsBuffer.ensureCapacity(numKeys);
		for (int i = 0; i < numKeys; i++) {			
			proportionsBuffer.add(calcProportionAtIndex(i));
		}
		return proportionsBuffer;
	}
	
	// ========================================================================
	public FloatArrayList getCumulativeProportionsForward(FloatArrayList proportionsBuffer) {
		proportionsBuffer = getProportions(proportionsBuffer);
		
		// Return early if empty list		
		int numKeys = getNumKeys();
		if (numKeys == 0) return proportionsBuffer;
		
		for (int i = 1; i < numKeys; i++) {
			proportionsBuffer.set(i, proportionsBuffer.get(i) + proportionsBuffer.get(i - 1));			
		}	
		
		return proportionsBuffer;
	}
	
	// ========================================================================
	public IntArrayList getCumulativeSumsForward(IntArrayList cumulativeCountsBuffer) {
		cumulativeCountsBuffer = (cumulativeCountsBuffer == null) ? new IntArrayList() : cumulativeCountsBuffer;
		cumulativeCountsBuffer.clear();		
		
		// Return early if empty list		
		int numKeys = getNumKeys();
		if (numKeys == 0) return cumulativeCountsBuffer;
		cumulativeCountsBuffer.ensureCapacity(numKeys);
		
		// First initialize first element
		cumulativeCountsBuffer.add(getCountAtIndex(0));
		
		for (int i = 1; i < numKeys; i++) {
			cumulativeCountsBuffer.add(cumulativeCountsBuffer.get(i - 1) + getCountAtIndex(i));			
		}
		return cumulativeCountsBuffer;
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
		dbc.getCumulativeProportionsForward(proportionsBufferCumulativeForward);
		
		for (int i = 0; i < dbc.getNumKeys(); i++) {
			System.out.printf("%d\t%d\t%d\t%d\t%d\t%g\t%g\n", i, dbc.getKeyAtIndex(i), dbc.getCountAtIndex(i), cumulativeForward.get(i), cumulativeBackward.get(i), proportionsBuffer.get(i), proportionsBufferCumulativeForward.get(i));
		}
		
		System.out.println("Key with max count: " + dbc.getKeyWithMaxCount());
		System.out.println("Sum of counts: " + dbc.mSumOfCounts);
		System.out.println("Total counts: " + dbc.mWeightedTotalSum);
	}
	
	// ========================================================================
	private static void TestDBC2() {
		DynamicBucketCounter dbc = new DynamicBucketCounter();		
		HashMap<PrimitiveWrapper.WInteger, PrimitiveWrapper.WInteger> truthMap = new HashMap<>();
		PrimitiveWrapper.WInteger dummy = new PrimitiveWrapper.WInteger(0);
		
		
		int numTrials = 100_000_000;		
		int printIncr = 100_000;
		for (int trial = 0; trial < numTrials; trial++) {
			if (trial % printIncr == 0) System.out.printf("Trial: %d\n", trial);
			int key = NumberUtils.getRandomInteger(0, 25);			
			dbc.incrementCount(key);
			
			// Now increment in truth table
			dummy.mInt = key;
			PrimitiveWrapper.WInteger result = truthMap.get(dummy);
			if (result == null) {
				truthMap.put(new PrimitiveWrapper.WInteger(key), new PrimitiveWrapper.WInteger(1));
			} else {
				result.mInt++;
			}
		}
		
		// Now check
		int maxTruthValue = -1;
		int keyWithMaxTruthValue = -1;
		
		for (PrimitiveWrapper.WInteger theKey : truthMap.keySet()) {
			PrimitiveWrapper.WInteger truthValue = truthMap.get(theKey);
			int testValue = dbc.getCount(theKey.mInt);
			
			if (truthValue.mInt != testValue) {
				System.err.printf("ERROR: Counts do not match: %d\t%d\t%d\n", theKey.mInt, truthValue.mInt, testValue);
				
			}
			
			if (truthValue.mInt > maxTruthValue) {
				maxTruthValue = truthValue.mInt;
				keyWithMaxTruthValue = theKey.mInt;
			}
		}
		
		// Now ensure the maximum tests are true		
		if (keyWithMaxTruthValue != dbc.getKeyWithMaxCount()) {
			System.err.printf("Max keys don't match: %d\t%d\n", keyWithMaxTruthValue, dbc.getKeyWithMaxCount());
		}
		System.out.println("Key Max: " + keyWithMaxTruthValue);
		
		if (maxTruthValue != dbc.getCountMax()) {
			System.err.printf("Max counts don't match: %d\t%d\t%d\n", keyWithMaxTruthValue, maxTruthValue, dbc.getCountMax());
		}		
		System.out.println("Count Max: " + maxTruthValue);
		
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
		//TestDynamicBucketCounter();
		TestDBC2();
	}
	
}