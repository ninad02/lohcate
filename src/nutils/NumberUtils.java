package nutils;

import genomeEnums.Nuc;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.EnumSet;
import java.util.Random;

import org.apache.commons.math3.distribution.BinomialDistribution;
import org.apache.commons.math3.distribution.PoissonDistribution;
import org.apache.commons.math3.stat.Frequency;

import com.carrotsearch.hppc.DoubleArrayList;
import com.carrotsearch.hppc.cursors.DoubleCursor;

import cern.jet.random.Binomial;
import cern.jet.random.engine.DRand;

public class NumberUtils {
	
	private static Binomial BinomialDistGlobalInstance = new Binomial(100, 0.5, new DRand());
	public static Random RandomGen = new Random();
	
	// ========================================================================
	public static final double Log10OfZero = -50;
	
	// ========================================================================
	/** If n is 0, this approximates n to a very small value. */
	public static double MathLog10Safe(double n) {
		return (n == 0.0) ? Log10OfZero : Math.log10(n); 
	}

	// ========================================================================
	public static double MathLog2(double n) { return MathLogBase(n, 2); }
	
	// ========================================================================
	public static double MathLogBase(double n, double base) {
		return (n == 0.0) ? Double.NaN : (Math.log10(n) / Math.log10(base));
	}

	/** Parses an integer written possibly and scientific format (with exponents) and returns the integer. 
	 * For example, 6.5e+07 is returned as 65000000 */
	public static int parseIntScientific(String token) {
		String sentinel = "e+";
		int result = token.indexOf(sentinel);
		
		if (result >= 0) {
			double thePrecision = Double.parseDouble(token.substring(0, result));
			int theScale = Integer.parseInt(token.substring(result + sentinel.length()));
			return ((int) Math.round(thePrecision * Math.pow(10, theScale)));
		} else {
			return Integer.parseInt(token);
		}
	}

	// ========================================================================
	// RANDOM FUNCTIONS
	// ========================================================================	
	public static void setSeed(long seed) { RandomGen.setSeed(seed); }
	
	// ========================================================================	
	public static boolean getRandomBit() { return RandomGen.nextBoolean(); }

	// ========================================================================	
	public static int getRandomInteger(int minValue, int maxValue) {
		return ((int) (RandomGen.nextDouble() * ((long) maxValue - (long) minValue + 1))) + minValue;
	}

	// ========================================================================	
	public static long getRandomLong(long minValue, long maxValue) {
		return ((long) (RandomGen.nextDouble() * (maxValue - minValue + 1))) + minValue;
	}


	// ========================================================================
	/** Given an integer, this returns the number of bits needed to represent that integer. */
	public static int calcNumBitsToRepresent(int num) {
		return Integer.SIZE - Integer.numberOfLeadingZeros(num);
	}
	
	// ========================================================================
	/** Adds two numbers.  Assumes the two numbers are in log format and that
	 *  the true (non-log) values for the numbers are very low. */
	public static final int LowLogSentinelDefault = ((int) Math.log(Double.MIN_VALUE)) + 1;
	
	public static double addLogNumbers(double logNum1, double logNum2) {
		int sentinel = LowLogSentinelDefault;  // do not change this.  
		//System.out.println(sentinel);

		// Break down the numbers if they are too large (to cause underflow)
		if ((logNum1 < sentinel) || (logNum2 < sentinel)) {
			if ((logNum1 < sentinel) && (logNum2 < sentinel)) {
				return (addLogNumbers(logNum1 - sentinel, logNum2 - sentinel) + sentinel);
			} else if (logNum1 < sentinel) {
				return logNum2;
			} else {
				return logNum1;
			}
		} else {
			return Math.log(Math.exp(logNum1) + Math.exp(logNum2));
		}
	}
	
	public static final double LowLog10SentinelDefault = ((int) Math.log10(Double.MIN_VALUE)) + 1;

	public static double addLog10NumbersOld(double logNum1, double logNum2) {
		final double sentinel = LowLog10SentinelDefault;  // do not change this.  
		
		// Break down the numbers if they are too large (to cause underflow)
		if ((logNum1 < sentinel) || (logNum2 < sentinel)) {
			if ((logNum1 < sentinel) && (logNum2 < sentinel)) {
				return (addLog10Numbers(logNum1 - sentinel, logNum2 - sentinel) + sentinel);
			} else if (logNum1 < sentinel) {
				return logNum2;
			} else {
				return logNum1;
			}
		} else {
			return Math.log10(Math.pow(10, logNum1) + Math.pow(10, logNum2));
		}
	}
	
	// ========================================================================
	/** Divides the divisor by the dividend.  If the divident is zero, then
	 *  the function returns the divBy0ReturnVal value. */
	public static double divideSafeByZero(double divisor, double dividend, double divBy0ReturnVal) {
		return (dividend == 0) ? divBy0ReturnVal : (divisor / dividend);
	}
	
	// ========================================================================
	public static final double InvalidLogProbabilityThreshold = 0.01;
	public boolean isInvalidLogProbability(double logProb) { return logProb > InvalidLogProbabilityThreshold; } 

	// ========================================================================
	public static final int getIndexOfPoissonWithHighestProbability(int k, PoissonDistribution[] dists) {
		double maxProb = -Double.MAX_VALUE;
		int    maxProbIndex = -1;
				
		for (int i = 0; i < dists.length; i++) {
			double prob = dists[i].probability(k);
			if (prob > maxProb) {
				maxProb = prob;
				maxProbIndex = i;
			}
		}
		
		return maxProbIndex;
	}

	// ========================================================================
	public static final int roundToNearest5(int num) {
		return (int) (Math.round(num / 5.0) * 5);
	}
	
	// ========================================================================
	public static double getCumProbInterval(Frequency freqDist, int center, int leeway) {
		return freqDist.getCumPct(center + leeway) - freqDist.getCumPct(center - leeway);
	}

	// ========================================================================
	/** Calculates the FDR via the Benjamini-Hochberg correction and returns the p-value
	 *  threshold.
	 *  
	 *  BH correction works as follows.  Sort the p-values from least to greatest.  Then,
	 *  starting from greatest p-value (end of array), move to least p-value (beginning 
	 *  of array) such that you reach a p-value that meets the condition:
	 *  
	 *  p-value[k] <= (k / array-length) * alpha
	 *  
	 *  NOTE: Since this is mathematics, k starts at 1, not 0
	 *  
	 */
	public static double getFDR_BenjaminiHochberg(double[] pValues, double fdrAlpha, boolean pointsIndependent) {
		double[] pValuesClone = pValues.clone();
		Arrays.sort(pValuesClone);
		
		double c_m = 1.0;  // Default under condition of independence
		
		if (!pointsIndependent) {
			double eulerMascheroniConstant = 0.577215664901532;
			double e_m = 1.0 / (2 * pValuesClone.length);
			c_m = Math.log(pValuesClone.length) + eulerMascheroniConstant + e_m;
		}
		
		for (int k = pValuesClone.length - 1; k >= 0; k--) {
			double threshold = ((double) (k + 1) / ((double) pValuesClone.length * c_m)) * fdrAlpha;
			//System.out.println("" + k + "\t" + pValuesClone[k] + "\t" + threshold);
			if (pValuesClone[k] <= threshold) {
				return pValuesClone[k];
			}
		}
		
		return 0;
	}

	public static int incrementNum(int x) { return x + 1; }

	/** This returns the probability of k successes in an experiment of n trials, with
	 *  the probability of success for a trial represented as p
	 * @param n The number of trials
	 * @param k The number of successes
	 * @param p The probability of one success
	 * @return The probability P(X >= k), where X is a random variable of number of successes
	 */
	public static double cumulativeProbabilitySuccessBinomial(int n, int k, double p) throws IllegalArgumentException {
		
		if (k > n) {
			throw new IllegalArgumentException("Argument k cannot be greater than argument n!");
		}
		
		if ((p < 0) || (p > 1)) {
			throw new IllegalArgumentException("Argument p must be between 0 and 1, inclusive!");
		}
		
		if (n <= 0) {
			return 1.0;
		}
		
		if (k == 0) {
			return 1.0;  // since probability Pr(X >= 0) is always 1.0 with positive number of trials
		}
	
		//BinomialDistribution bd = new BinomialDistribution(n, p);				
		//double result = 1 - bd.cumulativeProbability(k - 1);
		BinomialDistGlobalInstance.setNandP(n, p);
		double result = 1 - BinomialDistGlobalInstance.cdf(k - 1);
		
		return result;
	}
	
	/*
	public static double binomDist_PMF(int n, int k, double p) throws IllegalArgumentException {
		
		if (k > n) {
			throw new IllegalArgumentException("Argument k cannot be greater than argument n!");
		}
		
		if ((p < 0) || (p > 1)) {
			throw new IllegalArgumentException("Argument p must be between 0 and 1, inclusive!");
		}
		
		if (n <= 0) {
			return 1.0;
		}
		
		double nChooseKResultLog = nChooseK(n, k, true);
		
	}*/

	/** Given an integer array of counts, this picks a random number, and based
	 *  on the proportion of the counts, this returns the index of the chosen count.
	 */
	public static int pickRandomCountFromArray(int[] counts) {
		int arraySum = ArrayUtils.arraySum(counts);		
		if (arraySum == 0) {
			return -1;
		}
		
		int randomInt = (int) (Math.random() * arraySum);
		int runningSum = 0;
		for (int i = 0; i < counts.length; i++) {
			if (counts[i] != 0) {
				if ((randomInt >= runningSum) && (randomInt < runningSum + counts[i])) {
					return i;
				}
				runningSum += counts[i];
			}
		}
		
		return counts.length - 1;		
	}

	/** Given a base and the exponent, this raises the base to the exponent.  If the
	 *  
	 * @param base The base that is to be raised
	 * @param exponent The number by which to raise the base
	 * @param resultInLog Whether the return value should be in natural log form or not
	 * @return The result in regular or natural log form, based on the @param resultInLog argument
	 */
	public static double power(double base, double exponent, boolean resultInLog) {
		if (exponent == 0) {
			return 1;
		}
		
		double result = Math.log(base) * exponent;		
		if (resultInLog) {
			return result;
		} else {
			return Math.exp(result);
		}
	}

	public static double getPvalueFromLODBase10(double lodValue) {
		return Math.pow(10, -1 * lodValue);
	}

	/** Given a collection of doubles, this returns the mean and standard dev. The
	 *  isPopulation argument is needed to determine the degrees of freedom. */
	public static double calculateMeanAndSD(Collection<Double> c, PrimitiveWrapper.WDouble sDev, boolean isPopulation) {
		double total = 0;
		int numElem = 0;		
		for (Double d : c) {						
			total += d.doubleValue();
			numElem++;
		}		
		double mean = total / numElem;
		
		// Now calculate the SD
		double totalDeviation = 0;
		for (Double d : c) {
			double n = d.doubleValue();
			totalDeviation += ((n - mean) * (n - mean));
		}
		int degFreedom = (isPopulation) ? numElem : numElem - 1;
		if (sDev != null) {
			sDev.mDouble = Math.sqrt(totalDeviation / (double) degFreedom);
		}
		return mean;		
	}

	/** This returns the probability of exactly k successes in an experiment of n trials, 
	 *  with the probability of success for a trial represented as p
	 * @param n The number of trials
	 * @param k The number of successes
	 * @param p The probability of one success
	 * @return The probability P(X == k), where X is a random variable of number of successes
	 */
	public static double inidividualProbabilitySuccess(int n, int k, double p) {
		CompareUtils.ensureTrue(k <= n, "ERROR: CopyNumberMath.individualProbabilitySuccess(): Argument k cannot be greater than argument n!");
		CompareUtils.ensureTrue(k >= 0, "ERROR: CopyNumberMath.individualProbabilitySuccess(): Argument k must be 0 or greater!");		
		CompareUtils.ensureTrue((p >= 0) && (p <= 1), "Argument p must be between 0 and 1, inclusive!");
				
		if (n <= 0) {
			return 1.0;
		}
				
		BinomialDistribution bd = new BinomialDistribution(n, p);
		return bd.probability(k);
	}

	/** This returns the probability of  at least k successes in a poisson process of mean lambda. */
	public static double cumulativeProbabilitySuccessPoisson(int k, double lambda) {
		if (k < 0) {
			throw new IllegalArgumentException("Argument k must be >= 0!");
		}
		
		if (k == 0) {
			return 1.0;
		}
				
		PoissonDistribution pd = new PoissonDistribution(lambda); 				
		int kNew = Math.max(k - 1, 0);
		
		double result = 1 - pd.cumulativeProbability(kNew);		
		return result;
	}

	public static double nChooseK(int n, int k, boolean resultInLog) throws IllegalArgumentException {
		if (k > n) {
			throw new IllegalArgumentException("The argument k cannot be larger than the argument n!");
		}
		
		// Here, we are going to test k in relation to n.  If k is greater
		// than (n - k), we will switch the value of k and (n - k), as this
		// will help optimize the function, as seen later.
		k = Math.min(k, n - k);
		
		// First, add the numerator up (in log)
		double result = Math.log(n);
		int nMinusK = n - k;		
		for (int i = n - 1; i > nMinusK; i--) {
			result += Math.log(i);
		}
		
		// Next, subtract the demoninator out (in log)
		for (int i = k; i > 0; i--) {
			result -= Math.log(i);
		}
		
		if (resultInLog) {
			return result;
		} else {
			return Math.exp(result);
		}
	}

	// ========================================================================
	/** Given a number of trials and a probability of success, this returns the number of trials that succeed via random process. */
	public static int numSuccessesInTrials(int numTrials, double probabilityOfSuccess, int minNumberOfSuccesses) {
		int sum = 0;
		for (int i = 0; i < numTrials; i++) {
			sum += ((Math.random() <= probabilityOfSuccess) ? 1 : 0);
		}		
		return Math.max(sum, minNumberOfSuccesses);
	}

	// ========================================================================
	public static double roundToNthDecimal(double num, int numDecimalPlaces) {
		// We need to ensure we don't perform any buffer overruns on the double,
		// therefore we will need to use a better method than just multiplying by
		// 10^n, rounding, then dividing by 10^n.
		double theDecimalPortion = num - ((int) num);
		double multiplier = Math.pow(10, numDecimalPlaces);
		theDecimalPortion = Math.round(theDecimalPortion * multiplier) / multiplier; 
		return ((int) num) + theDecimalPortion;
	}
	
	// ========================================================================
	public static double addLog10Numbers(double logNum1, double logNum2) {
		double greater = Math.max(logNum1, logNum2);
		double smaller = Math.min(logNum1, logNum2);
		return (greater + Math.log10(1 + Math.pow(10, smaller - greater)));
	}
	
	// ========================================================================
	public static boolean isOdd(long num) { return ((num & 0x01L) == 1); }
	
	// ========================================================================
	public static boolean isEven(long num) { return ((num & 0x01L) == 0); }
	
	// ========================================================================
	public static <E extends Enum<E>> ArrayList<EnumSet<E>> powerSetEnum(Class<E> enumClass, boolean allowEmptySet) {
		EnumSet<E> emptySet = EnumSet.noneOf(enumClass);
		E[] enumConstants = enumClass.getEnumConstants();		
		int maxValue = ((int) Math.round(Math.pow(2, enumConstants.length))) - 1;
		
		// Set the bits for each enum		
		int[] bitRepresentationEnumOrdinal = new int[enumConstants.length];		
		for (int i = 0; i < enumConstants.length; i++) {
			bitRepresentationEnumOrdinal[i] = (1 << enumConstants[i].ordinal());			
		}

		// The return value
		ArrayList<EnumSet<E>> allSets = new ArrayList<>(maxValue);
		for (int value = 0; value <= maxValue; value++) {
			EnumSet<E> tempSet = EnumSet.copyOf(emptySet);
			for (int i = 0; i < bitRepresentationEnumOrdinal.length; i++) {
				if ((value & bitRepresentationEnumOrdinal[i]) == bitRepresentationEnumOrdinal[i]) {
					tempSet.add(enumConstants[i]);
				}
			}
			
			if (allowEmptySet || !tempSet.isEmpty()) {
				allSets.add(tempSet);
			}
		}
		
		return allSets;		
	}
	
	// ========================================================================
	/** This function adds two log-base-10 numbers, with the added boolean condition.  
	 *  The boolean parameter, if true, makes the function just return the first 
	 *  number.  If false, the function actually returns the sum of both numbers
	 *  This convenience function helps when tallying log numbers in a loop.
	 */
	public static double addLog10Numbers(boolean useOnlyFirstTerm, double logNum1, double logNum2) {
		return (useOnlyFirstTerm ? logNum1 : addLog10Numbers(logNum1, logNum2));
	}

	// ========================================================================
	private static void TestAddLog10Numbers() {
		System.out.println(Double.MIN_VALUE);
		System.out.println(Double.MIN_VALUE - 1);
		System.out.println(Double.MAX_VALUE);
		System.out.println(Math.log10(Double.MAX_VALUE));
		System.out.println(LowLog10SentinelDefault);
		
		double logNum1 = -5000;
		double logNum2 = -499;
		System.out.println(addLog10Numbers(logNum1, logNum2));
		System.out.println(addLog10NumbersOld(logNum1, logNum2));
	}
	
	// ========================================================================
	private static void TestDoubleListFromString() {
		String s = "{ 2.5, 3, 3.6, 4.17e-10, 1001 }";
		DoubleArrayList dArr = ArrayUtils.getDoubleListFromStringForm(s, true);
		for (DoubleCursor d : dArr) {
			System.out.print(d.value + " ");
		}
		System.out.println("");
	}
	
	// ========================================================================
	private static void TestIntListFromString() {
		String s = "{60-60}";
		ArrayList<Integer> iArr = ArrayUtils.getIntListFromStringForm(s);
		for (Integer i : iArr) {
			System.out.print(i + " ");
		}
		System.out.println("");
	}
	
	

	// ========================================================================
	private static enum TestEnum { A, B, C, D };
	
	
	// ========================================================================
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		//TestAddLog10Numbers();
		//TestDoubleListFromString();
		//TestIntListFromString();;
		
		powerSetEnum(TestEnum.class, false);
		
		/*
		String dd = "2,3,4,5";
		System.out.println(ArrayUtils.stripBraces(dd));
		
		String listAsStr = "({2.5,3,2.8};{1.1,1.6,0.5};{0.2,10.3,7.5})";
		ArrayList<DoubleArrayList> theList = ArrayUtils.getListOfDoubleListsFromStringForm(listAsStr, false);
		for (DoubleArrayList d : theList) {
			System.out.println(d);
		}
		
		for (int i = 0; i < 10000; i++) {
		//	System.out.println(getRandomInteger(0, 1));
		}
		*/
	}
	
}
