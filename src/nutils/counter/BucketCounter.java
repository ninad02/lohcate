/**
 * 
 */
package nutils.counter;

import java.util.Arrays;

/**
 * @author Ninad Dewal
 *
 */
public class BucketCounter extends BucketCounterCore {
	
	// ========================================================================
	public BucketCounter(int arrayLength, int valueStart) {
		super(arrayLength, valueStart);
	}
	
	// ========================================================================
	public BucketCounter(BucketCounter rhs) {
		super(rhs);
	}
	
	// ========================================================================
	public BucketCounter getCopy() {
		return new BucketCounter(this);
	}

	// ========================================================================
	public void increment(int bucket) { 
		increment(bucket, 1);
	}
	
	// ========================================================================
	public void increment(int bucket, int numTimesToIncrement) {
		super.incrementBucket(bucket, numTimesToIncrement);
	}

	// ========================================================================
	public int getCount(int bucket) { return super.getCountBucket(bucket); }
	
	// ========================================================================
	public double getProportion(int bucket) { return super.getProportionBucket(bucket); }
 
	// ========================================================================
	public static void TestBucketCounter() {
		int numBuckets = 8;
		int valueStart = 1;
		BucketCounter bc = new BucketCounter(numBuckets, valueStart);
		for (int i = 0; i < numBuckets; i++) {
			int bucketValue = i + valueStart;
			bc.increment(bucketValue, bucketValue * bucketValue);
		}
		BucketCounterCore bcCopy = bc.getCopy();
		bc.print(System.out);
		bcCopy.print(System.out);
	}
	
	// ========================================================================	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		TestBucketCounter();
	}

}
