package nutils.counter;

import java.io.PrintStream;
import java.util.Arrays;

import nutils.ArrayUtils;
import nutils.NumberUtils;

public class BucketCounter {
	private int[] mArray;		
	private int mIndexStart;
	private int mSumOfCounts;
	private int mTotalSum;
	
	public BucketCounter(int arrayLength, int valueStart) {
		mIndexStart = valueStart;
		mArray = new int[arrayLength];
		clear();
	}
	
	public BucketCounter(BucketCounter rhs) {
		this(rhs.mArray.length, rhs.mIndexStart);		
		add(rhs);
	}
	
	public BucketCounter getCopy() { return new BucketCounter(this); }
	
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
	
	public BucketCounter downsampleCounts(double probabilityToRetainCount) {
		int thisLength = getLength();
		BucketCounter bcSubsampled = new BucketCounter(thisLength, mIndexStart);
		
		for (int b = 0; b < thisLength; b++) {
			int trueBucketValue = b + mIndexStart;
			int numSitesInBucket = getCount(trueBucketValue);
			for (int j = 0; j < numSitesInBucket; j++) {
				int downsampledNumReads = NumberUtils.numSuccessesInTrials(trueBucketValue, probabilityToRetainCount, 0);
				bcSubsampled.increment(downsampledNumReads);
			}
		}	
		return bcSubsampled;
	}
	
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
		for (int i = 0; i < mArray.length; i++) {
			if (printZeroCountElements || (mArray[i] > 0)) {
				int bucketValue = i + mIndexStart;
				out.println(bucketValue + "\t" + getCount(bucketValue) + "\t" + getProportion(bucketValue));
			}
		}		
	}
	
	/** Writes all the counts to the string builder passed in. */
	public void constructString(StringBuilder sb, boolean clearStringBuilder, String delimiter) {
		if (clearStringBuilder) {
			sb.setLength(0);
		} 
		
		for (int i = 0; i < mArray.length; i++) {
			if ((i > 0) || ((i == 0) && !clearStringBuilder)) {
				sb.append(delimiter);
			}
			sb.append(mArray[i]);
		}	
	}
	
	public static void TestBucketCounter() {
		int numBuckets = 8;
		int valueStart = 1;
		BucketCounter bc = new BucketCounter(numBuckets, valueStart);
		for (int i = 0; i < numBuckets; i++) {
			int bucketValue = i + valueStart;
			bc.increment(bucketValue, bucketValue * bucketValue);
		}
		BucketCounter bcCopy = bc.getCopy();
		bc.print(System.out);
		bcCopy.print(System.out);
	}
	
	public static void main(String[] args) {
		TestBucketCounter();
	}
}
