package nutils.counter;

import java.io.BufferedReader;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Enumeration;
import java.util.Hashtable;
import java.util.LinkedList;

import org.apache.commons.math3.stat.inference.ChiSquareTest;
import org.apache.commons.math3.stat.inference.GTest;

import nutils.ArrayUtils;
import nutils.CompareUtils;
import nutils.IOUtils;
import nutils.NumberUtils;

public class DecimalCounter {


	// ========================================================================
	
	
	// ========================================================================
	// ===== INNER CLASS =====
	// ========================================================================
	public static class DoubleValuesWithCounts {
		public double[] mValues;
		public int[] mCounts;
		public int[] mCumulativeCountsForward;
		public int[] mCumulativeCountsBackward;
		
		public DoubleValuesWithCounts(int arraySize) {
			mValues = new double[arraySize];
			mCounts = new int[arraySize];
		}
		
		// Given a filename, this reads in the lod list
		public static DoubleValuesWithCounts getLodListsWithCountsFromFile(String inFilename, double minLodOfInterest, double maxLodOfInterest) {				
			ArrayList<String> allLines = IOUtils.readAllLinesFromFile(inFilename);
			LinkedList<Double> lodScoresList = new LinkedList<Double>();
			LinkedList<Integer> countsList = new LinkedList<Integer>();
			
			for (int i = 0; i < allLines.size(); i++) {
				String[] components = allLines.get(i).split("\\s");
				double lodScore = Double.parseDouble(components[0]);
				if ((lodScore >= minLodOfInterest) && (lodScore <= maxLodOfInterest)){
					int count = Integer.parseInt(components[1]);
					lodScoresList.add(new Double(lodScore));
					countsList.add(new Integer(count));
				}
			}
	
			DoubleValuesWithCounts dvwc = new DoubleValuesWithCounts(lodScoresList.size());
			dvwc.mValues = ArrayUtils.getPrimitiveArrayDouble(lodScoresList);
			dvwc.mCounts    = ArrayUtils.getPrimitiveArrayInt(countsList);				
		//	dvwc.mCumulativeCountsForward = DynamicRoundedDoubleCounter.getCumulativeCounts(dvwc.mCounts, true);
		//	dvwc.mCumulativeCountsBackward = DynamicRoundedDoubleCounter.getCumulativeCounts(dvwc.mCounts, false);
			return dvwc;				
		}			
	}
	
	// ========================================================================
	// ===== INNER CLASS =====
	// ========================================================================
	static class RoundedDoubleCountPackage implements Comparable<RoundedDoubleCountPackage> {
		int mLodValue = 0;
		int mCount = 0;

		public int compareTo(RoundedDoubleCountPackage rhs) {
			return CompareUtils.compareInt(rhs.mLodValue, mLodValue);
		}
	}

	// ========================================================================
	// ===== INNER CLASS =====
	// ========================================================================
	public static class PValueCounter {
		private static final int GranularityDefault = 10;		
		int[] pValueCounts;  // there should be 10001, thus index 0 represents pvalue of 0, 1 represetns pvalue of 0.01
		int mGranularity;
		
		public PValueCounter() {
			constructorHelper(GranularityDefault);
		}
		
		public PValueCounter(int granularityMagnitude) {
			int granularity = (int) Math.round(Math.pow(10, granularityMagnitude));
			constructorHelper(granularity);
		}
		
		private void constructorHelper(int granularity) {
			mGranularity = granularity;
			pValueCounts = new int[mGranularity + 1];
			Arrays.fill(pValueCounts, 0);
		}
		
		/** Returns the counts of the p values. */
		public int[] getPValueCounts() { return pValueCounts; }
		
		/** Returns the number of elements. */
		public int getNumBins() { return pValueCounts.length; }
		
		public int submitPValue(float p, boolean shouldRound) {
			// First, we round this to the correct decimal position to get the index
			int index = shouldRound ? Math.round(p * mGranularity) : (int) (p * mGranularity);
			return (++(pValueCounts[index]));
		}						
		
		/** Returns the sum of the counts. */
		public int getTotalCount() {
			int sum = 0;
			for (int i = 0; i < pValueCounts.length; i++) {
				sum += pValueCounts[i];
			}
			return sum;
		}
		
		
		/** Returns the cumulative counts of the p values over the buckets. */
		public int[] getPValueCountsCumulative() {
			int[] cumCounts = (int[]) pValueCounts.clone();
			
			for (int i = 1; i < cumCounts.length; i++) {
				cumCounts[i] += cumCounts[i - 1];
			}
			
			return cumCounts;
		}
		
		/** Returns the cumulative counts of the p values over the buckets and tests against a uniform distribution. */
		public double performChiSquareWithUniformExpected(PrintStream out) {
			double[] observed = getProportions();
			int multiplier = 100;					
			
			double expectedValue = Math.max(1, (1.0 / (double) getNumBins()) * multiplier);
			double[] expected = new double[observed.length];
			Arrays.fill(expected, expectedValue);
			
			long[] observedLong = new long[observed.length];
			for (int i = 0; i < observed.length; i++) {
				observedLong[i] = Math.round((observed[i] * multiplier));
			}
			
			if (nutils.CompareUtils.isNotNull(out)) {
				out.println("-----\nChi Square or G Test\n-----\n");
				for (int i = 0; i < expected.length; i++) {
					out.printf("%d\t%g\t%g\n", i, expected[i], observed[i]);
				}
			}
			
			GTest gTest = new GTest();			
			//org.apache.commons.math3.stat.inference.ChiSquareTest
			ChiSquareTest csTest = new ChiSquareTest();			
			double pVal = csTest.chiSquareTest(expected, observedLong);
			//double pVal = gTest.gTest(expected, observedLong);
			return pVal;
		}
		
		/** Returns the proportions of the counts. */
		public double[] getProportions() {
			double[] proportions = new double[pValueCounts.length];
			int sumOfCounts = getTotalCount();
			
			for (int i = 0; i < proportions.length; i++) {
				proportions[i] = pValueCounts[i] / (double) sumOfCounts;
			}
			
			return proportions;
		}
		
		/** Returns the cumulative log counts of the p values over the buckets. */
		public double[] getPValueCountsCumulativeLog(boolean forward) {
			int[] results = forward ? getPValueCountsCumulative() : getPValueCountsCumulativeReverse();
			double[] rV = new double[results.length];
			
			for (int i = 0; i < results.length; i++) {
				rV[i] = NumberUtils.MathLog10Safe(results[i]);
			}
			return rV;
		}
		
		/** Returns the cumulative counts of hte p values over the buckets in the reverse direction. */
		public int[] getPValueCountsCumulativeReverse() {
			int[] cumCounts = (int[]) pValueCounts.clone();
			
			for (int i = cumCounts.length - 2; i >= 0; i--) {
				cumCounts[i] += cumCounts[i + 1];
			}
			
			return cumCounts;
		}
		
		
		/** Returns the proportion of each count (over the cumulative count). */
		public float[] getPValueCountsCumulativeProportion(boolean calcInLog) {
			int[] cumCounts = getPValueCountsCumulative();
			float totalCount = (float) cumCounts[cumCounts.length - 1];
			
			float[] countProp = new float[cumCounts.length];
			for (int i = 0; i < countProp.length; i++) {
				countProp[i] = cumCounts[i] / totalCount; 
				countProp[i] = (float) (calcInLog ? NumberUtils.MathLog10Safe(countProp[i]) : countProp[i]);
			}
			return countProp;
		}
		
		public static void Test() {
			PValueCounter pvc = new PValueCounter(2);
			for (int i = 0; i <= 100; i++) {
				pvc.submitPValue((float) i / (float) 100.0, true);
				pvc.submitPValue((float) i / (float) 100.0, true);
			}
			
			int[] counts = pvc.getPValueCounts();
			int[] cumCounts = pvc.getPValueCountsCumulative();
			
			for (int i = 0; i < counts.length; i++) {
				System.out.println(counts[i] + "\t" + cumCounts[i]);
			}
		}
		
		private static void Test2() {
			PValueCounter pvc = new PValueCounter(1);
			
			for (double d = 0; d <= 1.0; d += 0.1) {
				int maxIter = (int) Math.round(d * 10);
				System.out.printf("Val: %g\t%d\n", d, maxIter);
				for (int j = 0; j < maxIter; j++) {					
					pvc.submitPValue((float) d, true);
				}
			}
			
			System.out.println("Total count: " + pvc.getTotalCount());
			
			double[] proportions = pvc.getProportions();
			int[] counts = pvc.getPValueCounts();
			
			for (int i = 0; i < counts.length; i++) {
				System.out.println(counts[i] + "\t" + proportions[i]);
			}
		}
	}

	// ========================================================================
	// ===== INNER CLASS =====
	// ========================================================================
	public static class ValueDoubleCounterPackage implements Comparable<ValueDoubleCounterPackage> {
		public double mValue;
		public double mCounts;
		public double mCountsCumulative;
		public double mProportions;
		public double mProportionsCumulative;
		
		public ValueDoubleCounterPackage(double value) {
			mValue = value;
			mCounts = mCountsCumulative = mProportions = mProportionsCumulative = 0;
		}
		
		public int compareTo(ValueDoubleCounterPackage rhs) {
			return Double.compare(mValue, rhs.mValue);
		}
	}

	// ========================================================================
	// ===== INNER CLASS =====
	// ========================================================================
	// UNFINISHED
	public static class PValueCounterSparse {
		Hashtable<Double, ValueDoubleCounterPackage> mCounterTable;
		
		public PValueCounterSparse() {
			mCounterTable = new Hashtable<Double, ValueDoubleCounterPackage>();
		}
		
		public void clear() {
			mCounterTable.clear();
		}
		
		/** Rounds the p-value to the 2nd non-zero digit and increments the count. */
		public int submitPValue(double p, int nthNonZeroDigit) {
			nthNonZeroDigit = Math.max(1, nthNonZeroDigit);
			int exponent = ((int) (-1 * NumberUtils.MathLog10Safe(p))) + nthNonZeroDigit;
			long multiplier = Math.round(Math.pow(10, exponent));
			long pShifted = Math.round(p * multiplier);
			double pNew = (double) pShifted / (double) multiplier;
			
			ValueDoubleCounterPackage pvp = mCounterTable.get(pNew);
			if (pvp == null) {
				pvp = new ValueDoubleCounterPackage(pNew);
				mCounterTable.put(pNew, pvp);
			} 
			pvp.mCounts++;
			return (int) Math.round(pvp.mCounts);
		}
		
		/** Returns the list of p-values, counts, cumulative counts, and proportions. */ 
		public ArrayList<ValueDoubleCounterPackage> getCountsAndProportions(boolean reverseList) {
			ArrayList<ValueDoubleCounterPackage> pvpList = new ArrayList<ValueDoubleCounterPackage>(mCounterTable.size());			
			for (Enumeration<ValueDoubleCounterPackage> e = mCounterTable.elements(); e.hasMoreElements(); ) {
				ValueDoubleCounterPackage pvp = e.nextElement();
				pvpList.add(pvp);
			}
			
			// Next, sort the list
			Collections.sort(pvpList);
			if (reverseList) {
				Collections.reverse(pvpList);
			}
			
			// Now that the list is sorted, we can find the cumulative counts and proportions
			if (pvpList.size() > 0) {
				pvpList.get(0).mCountsCumulative = pvpList.get(0).mCounts;
				for (int i = 1; i < pvpList.size(); i++) {
					ValueDoubleCounterPackage pvp = pvpList.get(i); 
					pvp.mCountsCumulative = pvpList.get(i - 1).mCountsCumulative + pvp.mCounts;
				}
				
				// Now, calculate the proportions
				double totalCount = pvpList.get(pvpList.size() - 1).mCountsCumulative;
				pvpList.get(0).mProportions = pvpList.get(0).mCounts / totalCount;
				pvpList.get(0).mProportionsCumulative = pvpList.get(0).mProportions; 
				for (int i = 1; i < pvpList.size(); i++) {
					ValueDoubleCounterPackage pvp = pvpList.get(i);
					pvp.mProportions = pvp.mCounts / totalCount;
					pvp.mProportionsCumulative = pvpList.get(i - 1).mProportionsCumulative + pvp.mProportions;
				}				
			}
			
			return pvpList;
		}
	
		/** Given a filename, this reads the values of a column in that filename into a pvc. If
		 *  an existing p-value counter object is provided, it is used.  Else if null is passed
		 *  in for this, a new one is created. */
		public static PValueCounterSparse readColumnValuesIntoPVC(PValueCounterSparse existingPVC, String inFilename, int columnOfInterest, boolean hasHeader, int numSignificantDigits) {
			DecimalCounter.PValueCounterSparse pvc = (existingPVC == null) ? new DecimalCounter.PValueCounterSparse() : existingPVC;			
			
			BufferedReader in = IOUtils.getBufferedReader(inFilename);
			int lineCounter = -1;
			String line;			
			while ((line = IOUtils.getNextLineInBufferedReader(in)) != null) {
				lineCounter++;
				if (hasHeader && (lineCounter == 0)) continue;
				
				String components[] = line.split("\\s");
				double pvalue = Double.parseDouble(components[columnOfInterest]);
				pvc.submitPValue(pvalue, numSignificantDigits);			
			}
			IOUtils.closeBufferedReader(in);
			return pvc;
		}
	}

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		PValueCounter.Test2();
	}

}
