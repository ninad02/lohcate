package nutils.counter;

import java.io.PrintStream;
import java.util.Arrays;

import nutils.CompareUtils;
import nutils.NumberUtils;

public abstract class BucketCounterCore {
	
	// ========================================================================
	protected int[] mArray;		
	protected int mIndexStart;
	protected int mSumOfCounts;
	protected int mTotalSum;
	
	// ========================================================================
	protected BucketCounterCore(int arrayLength, int valueStart) {
		mIndexStart = valueStart;
		mArray = new int[arrayLength];
		clear();
	}
	
	// ========================================================================
	protected BucketCounterCore(BucketCounterCore rhs) {
		this(rhs.mArray.length, rhs.mIndexStart);		
		add(rhs);
	}

	// ========================================================================
	public void clear() { 
		Arrays.fill(mArray, 0); 
		mSumOfCounts = 0; 	
		mTotalSum = 0;
	}
	
	// ========================================================================
	public int getLength() { return mArray.length; }

	// ========================================================================
	public double calcMean() { return (double) mTotalSum / (double) mSumOfCounts; }
	
	// ========================================================================
	public void print(PrintStream out) {			
		print(out, true);
	}
	
	// ========================================================================
	public void print(PrintStream out, boolean printZeroCountElements) {
		for (int i = 0; i < mArray.length; i++) {
			if (printZeroCountElements || (mArray[i] > 0)) {
				int bucketValue = i + mIndexStart;
				out.printf("%d\t%d\t%g\n", bucketValue, getCountBucket(bucketValue), getProportionBucket(bucketValue));
			}
		}		
	}
	
	// ========================================================================
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
	
	// ========================================================================
	public boolean add(BucketCounterCore bc) {
		if (!this.sameLengthAs(bc))        return false;
		if (!this.sameStartingIndexAs(bc)) return false;
		
		int bcLength = bc.getLength();
		for (int b = 0; b < bcLength; b++) {
			int trueOtherBucketValue = b + bc.mIndexStart;
			int otherCount = bc.getCountBucket(trueOtherBucketValue);
			this.incrementBucket(trueOtherBucketValue, otherCount);
		}
		return true;
	}
	
	// ========================================================================
	public boolean sameLengthAs(BucketCounterCore bc) {
		return this.mArray.length == bc.mArray.length;
	}
	
	// ========================================================================
	public boolean sameStartingIndexAs(BucketCounterCore bc) {
		return this.mIndexStart == bc.mIndexStart;
	}
	
	// ========================================================================
	protected void incrementBucket(int bucket) { 
		incrementBucket(bucket, 1);
	}
	
	// ========================================================================
	protected void incrementBucket(int bucket, int numTimesToIncrement) {
		mArray[bucket - mIndexStart] += numTimesToIncrement;
		mSumOfCounts                 += numTimesToIncrement;
		mTotalSum                    += (numTimesToIncrement * bucket);
	}

	// ========================================================================
	protected int getCountBucket(int bucket) { return mArray[bucket - mIndexStart]; }
	
	// ========================================================================
	protected double getProportionBucket(int bucket) { return (double) getCountBucket(bucket) / (double) mSumOfCounts; }
	
	
	// ========================================================================
	protected void downsampleCounts(double probabilityToRetainCount, BucketCounterCore bcSubsampled) {
		int thisLength = getLength();
		
		CompareUtils.ensureTrue(sameLengthAs(bcSubsampled),        "ERROR: BucketCounter.downsampleCounts(): Buffer BucketCounter must be of same length as original!");
		CompareUtils.ensureTrue(sameStartingIndexAs(bcSubsampled), "ERROR: BucketCounter.downsampleCounts(): Buffer BucketCounter must have same starting index as that of the original!");
		
		for (int b = 0; b < thisLength; b++) {
			
			int trueBucketValue = b + mIndexStart;
			int numSitesInBucket = getCountBucket(trueBucketValue);
			
			for (int j = 0; j < numSitesInBucket; j++) {
				int downsampledNumReads = NumberUtils.numSuccessesInTrials(trueBucketValue, probabilityToRetainCount, 0);
				bcSubsampled.incrementBucket(downsampledNumReads);
			}
		}			
	}

	// ========================================================================
	public static void main(String[] args) {		
	}
}
