package shared;

import java.io.PrintStream;
import java.util.Arrays;

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
		for (int i = 0; i < getLength(); i++) {
			if (printZeroCountElements || (getCount(i) > 0)) {
				out.println(i + "\t" + getCount(i) + "\t" + getProportion(i));
			}
		}
	}
}
