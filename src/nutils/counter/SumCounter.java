package nutils.counter;

/**
 * A simple convenience class for keeping track of a sum and reporting an average.
 * @author Ninad Dewal
 *
 */
public abstract class SumCounter<E extends SumCounter<E>> implements Comparable<E> {
	protected long mNumElements;	
	
	// ========================================================================
	public SumCounter() {
		clear();
	}

	// ========================================================================
	public abstract void clear();
	
	// ========================================================================
	public abstract double getMean(); 

	// ========================================================================
	public static class DoubleC extends SumCounter<DoubleC> {
		
		protected double mSum;

		// ========================================================================
		public DoubleC() { clear(); }
				
		// ========================================================================
		@Override
		public int compareTo(DoubleC o) {
			int result = Long.compare(this.mNumElements, o.mNumElements);
			if (result == 0) {
				result = Double.compare(this.mSum, o.mSum);
			} 
			return result;
		}

		// ========================================================================
		@Override
		public void clear() { mSum = mNumElements = 0; }

		// ========================================================================
		@Override
		public double getMean() { return mSum / (double) mNumElements; }
	}
		
	// ========================================================================	
	public static class LongC extends SumCounter<LongC> {
		
		protected long mSum;

		// ========================================================================
		public LongC() { clear(); }
		
		// ========================================================================
		@Override
		public int compareTo(LongC o) {
			int result = Long.compare(this.mNumElements, o.mNumElements);
			if (result == 0) {
				result = Long.compare(this.mSum, o.mSum);
			} 
			return result;
		}

		// ========================================================================
		@Override
		public void clear() { this.mNumElements = this.mSum = 0; }

		// ========================================================================
		@Override
		public double getMean() { return (double) mSum / (double) mNumElements; }
		
	}
	
}
