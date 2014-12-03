package nutils.counter;

import com.carrotsearch.hppc.DoubleArrayList;
import com.carrotsearch.hppc.FloatArrayList;
import com.carrotsearch.hppc.IntArrayList;

import nutils.Cast;
import nutils.BitUtils.Compactor.CompactorInf;
import nutils.BitUtils.Compactor.CompactorIntoInt;
import nutils.BitUtils.Compactor.CompactorIntoLong;

// ========================================================================
// ===== INNER CLASS =====
// ========================================================================
public class DynamicRoundedDoubleCounter {
	
	// ========================================================================
	private static final int SpaceMultiplier = 10;

	// ========================================================================
	protected int mNumSignificantFigures;
	protected DynamicBucketCounter mRoundedDoubleList; 		
	
	// ========================================================================
	public DynamicRoundedDoubleCounter(int numSignificantFigures) {
		mNumSignificantFigures = Math.max(0, numSignificantFigures);		
		mRoundedDoubleList = new DynamicBucketCounter();				
	}
	
	// ========================================================================
	public void register(double value) {
		int roundedNumCompact = roundToSignificantFiguresCompact(value, mNumSignificantFigures);		
		mRoundedDoubleList.incrementCount(roundedNumCompact);
	}

	// ========================================================================
	private static double getDoubleValueFromCompactForm(int compactUnit) {
		int belowZero        = CompactDecimal.Compactor.getValue(CompactDecimal.BelowZero, compactUnit);
		int shifted          = CompactDecimal.Compactor.getValue(CompactDecimal.Base,      compactUnit);
		int exponentNegative = CompactDecimal.Compactor.getValue(CompactDecimal.ExponentNegative, compactUnit);
		int exponentBase10   = CompactDecimal.Compactor.getValue(CompactDecimal.Exponent_Base10, compactUnit);
				
		double magnitude = Math.pow(10, ((exponentNegative > 0) ? -exponentBase10 : exponentBase10));
		double value = shifted / magnitude;
		if (belowZero > 0) {
			value = -value;
		}
		return value;
	}
	
	// ========================================================================
	public double getValueWithMaxCount() {		
		return getDoubleValueFromCompactForm(mRoundedDoubleList.getKeyWithMaxCount());
	}
	
	// ========================================================================
	public void getValuesWithProportions(DoubleArrayList values, FloatArrayList proportions) {
		IntArrayList keysBuffer = new IntArrayList();
		mRoundedDoubleList.getKeys(keysBuffer);
		
		values.clear();
		values.ensureCapacity(keysBuffer.size());
		for (int i = 0; i < keysBuffer.size(); i++) {
			values.add(getDoubleValueFromCompactForm(keysBuffer.get(i)));
		}		
		mRoundedDoubleList.getProportions(proportions);
		
	}	

	// ========================================================================
	// SOURCE: StackOverflow
	// http://stackoverflow.com/questions/202302/rounding-to-an-arbitrary-number-of-significant-digits
	
	// ========================================================================
	/** Returns a rounded number in a compact form. */
	private static int roundToSignificantFiguresCompact(double num, int n) {
	    if (num == 0) {
	        return 0;
	    }
	    
	    int belowZero = (num < 0) ? 1 : 0;
	    final long d = Math.round(Math.ceil(Math.log10(num < 0 ? -num: num)));
	    final int power = n - (int) d;

	    final double magnitude = Math.pow(10, power);
	    final int shifted = Cast.toInt(Math.round(num * magnitude));
	    
	    int compactUnit = 0;
	    int exponentNegative = (power < 0) ? 1 : 0;
	    compactUnit = CompactDecimal.Compactor.setValue(CompactDecimal.BelowZero,               belowZero, compactUnit);
	    compactUnit = CompactDecimal.Compactor.setValue(CompactDecimal.Base,                      shifted, compactUnit);
	    compactUnit = CompactDecimal.Compactor.setValue(CompactDecimal.ExponentNegative, exponentNegative, compactUnit);
	    compactUnit = CompactDecimal.Compactor.setValue(CompactDecimal.Exponent_Base10,   Math.abs(power), compactUnit);
	    return compactUnit;
	}
	
	// ========================================================================	
	public static enum CompactDecimal implements CompactorInf<CompactDecimal> {
		BelowZero(1),		
		ExponentNegative(1),
		Exponent_Base10(14),
		Base(15),
		;

		public static CompactorIntoInt<CompactDecimal> Compactor = new CompactorIntoInt<>(CompactDecimal.class, false);
		
		private int mNumBits;
		private CompactDecimal(int numBits) { mNumBits = numBits; }
		public int getNumBits() { return mNumBits; }
	}

	// ========================================================================
	public static double roundToNumSigFigs(double num, int numSigFigs) {
		return getDoubleValueFromCompactForm(roundToSignificantFiguresCompact(num, numSigFigs));
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
		
		System.out.println(roundToNumSigFigs(2, 2));
		System.out.println(roundToNumSigFigs(2.1, 2));
		System.out.println(roundToNumSigFigs(2.14, 2));
		System.out.println(roundToNumSigFigs(2.16, 2));
		System.out.println(roundToNumSigFigs(2.163, 2));
		System.out.println(roundToNumSigFigs(21, 2));
		System.out.println(roundToNumSigFigs(21.4, 2));
		System.out.println(roundToNumSigFigs(21.6, 2));
		System.out.println(roundToNumSigFigs(0.2, 2));
		System.out.println(roundToNumSigFigs(0.21, 2));
		System.out.println(roundToNumSigFigs(0.214, 2));
		System.out.println(roundToNumSigFigs(0.216, 2));
		System.out.println(roundToNumSigFigs(0.02, 2));
		System.out.println(roundToNumSigFigs(0.021, 2));
		System.out.println(roundToNumSigFigs(0.0214, 2));
		System.out.println(roundToNumSigFigs(0.0216, 2));
		System.out.println(roundToNumSigFigs(0.002, 2));
		System.out.println(roundToNumSigFigs(0.0021, 2));
		System.out.println(roundToNumSigFigs(0.00214, 2));
		System.out.println(roundToNumSigFigs(0.00216, 2));
		System.out.println(roundToNumSigFigs(0.00021, 2));
		
		

		
		/*
		ValueDoubleCounterPackage llwc = ValueDoubleCounterPackage();
		
		
		for (int i = 0; i < llwc.mLodScores.length; i++) {
			System.out.println("" + llwc.mLodScores[i] + "\t" + llwc.mCounts[i] + "\t" + llwc.mCumulativeCountsForward[i] + "\t" + llwc.mCumulativeCountsBackward[i]);
		}*/
	}
	
	// ========================================================================
	private static void Test2() {
		int numSigFig = 1;
		DynamicRoundedDoubleCounter theCounter = new DynamicRoundedDoubleCounter(numSigFig);
		for (double d = 0; d <= 1.0; d += 0.05) {
			int repeatCount = (int) (Math.round(d * 10) + 1);
			System.out.println(roundToNumSigFigs(d, numSigFig));
			for (int j = 0; j <= repeatCount; j++) {				
				theCounter.register(d);
			}
		}
		
		DoubleArrayList values = new DoubleArrayList();
		FloatArrayList proportions = new FloatArrayList();
		
		theCounter.getValuesWithProportions(values, proportions);
		
		double totalProportion = 0;
		for (int i = 0; i < values.size(); i++) {
			System.out.printf("Value: %g\t%g\n", values.get(i), proportions.get(i));
			totalProportion += proportions.get(i);
		}
		System.out.println(totalProportion);
	}

	// ========================================================================
	public static void main(String[] args) {
		Test2();
	}
			
	
}