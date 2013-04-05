package nutils.counter;

// ========================================================================
// ===== INNER CLASS =====
// ========================================================================
public class DynamicRoundedDoubleCounter {
	
	// ========================================================================
	private static final int SpaceMultiplier = 10;

	// ========================================================================
	protected int mDecimalMultiplier;
	protected DynamicBucketCounter mRoundedDoubleList; 		
	
	// ========================================================================
	public DynamicRoundedDoubleCounter(int numDecimalsAllowed) {
		numDecimalsAllowed = Math.max(0, numDecimalsAllowed);
		mDecimalMultiplier = (int) Math.round(Math.pow(10, numDecimalsAllowed));
		mRoundedDoubleList = new DynamicBucketCounter(mDecimalMultiplier * SpaceMultiplier);				
	}
	
	// ========================================================================
	public void register(double value) {
		int valueAdjusted = getAdjustedDouble(value);
		mRoundedDoubleList.incrementCount(valueAdjusted);
	}
		
	// ========================================================================
	public double getValueWithMaxCount() {
		int value = mRoundedDoubleList.getKeyWithMaxCount();
		return getDeadjustedDouble(value);
	}
	
	// ========================================================================
	protected int getAdjustedDouble(double value) {
		return (int) Math.round(value * mDecimalMultiplier);
	}
	
	// ========================================================================
	protected double getDeadjustedDouble(int value) {
		return (value / ((double) mDecimalMultiplier));
	}
	
	// ========================================================================
	public static void Test() {
		
		DynamicRoundedDoubleCounter lodCounter = new DynamicRoundedDoubleCounter(2);
		
		lodCounter.register(7.135);
		lodCounter.register(7.135);
		lodCounter.register(7.135);
		lodCounter.register(7.1399);
		lodCounter.register(5.1);
		lodCounter.register(.91);
		lodCounter.register(.913);
		lodCounter.register(15.3);
		lodCounter.register(0.02);
		lodCounter.register(0.017);
		lodCounter.register(0.005);
		lodCounter.register(0.004);
		lodCounter.register(0.010006);
		
		System.out.println(lodCounter.getValueWithMaxCount());
		
		/*
		ValueDoubleCounterPackage llwc = ValueDoubleCounterPackage();
		
		
		for (int i = 0; i < llwc.mLodScores.length; i++) {
			System.out.println("" + llwc.mLodScores[i] + "\t" + llwc.mCounts[i] + "\t" + llwc.mCumulativeCountsForward[i] + "\t" + llwc.mCumulativeCountsBackward[i]);
		}*/
	}

	// ========================================================================
	public static void main(String[] args) {
		Test();
	}
			
	
}