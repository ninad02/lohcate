package shared;

import org.apache.commons.math3.distribution.BinomialDistribution;

public class NumberUtils {

	
	/** This returns the probability of k successes in an experiment of n trials, with
	 *  the probability of success for a trial represented as p
	 * @param n The number of trials
	 * @param k The number of successes
	 * @param p The probability of one success
	 * @return The probability P(X >= k), where X is a random variable of number of successes
	 */
	public static double cumulativeProbabilitySuccess(int n, int k, double p) throws IllegalArgumentException {
		
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

		BinomialDistribution bd = new BinomialDistribution(n, p);
		double result = 1 - bd.cumulativeProbability(k - 1);					
		return result;
	}
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub

	}

}
