package nutils.math;

import java.io.PrintStream;
import java.util.ArrayList;

import nutils.PrimitiveWrapper;
import nutils.collectionsSorted.ArrayListSorted;
import nutils.collectionsSorted.ArrayListSortedComparable;


/**
 * A class built to house multiple Poisson Distributions.  The primary purpose
 * of this class is to allow the user to determine the maximum likelihood Poisson
 * distribution for a given value.
 * @author Ninad Dewal
 *
 */

public class PoissonDistributionList {

	// //========================================================================
	ArrayListSortedComparable<PoissonDistributionSortable> mPdList;
		
	// ========================================================================
	public PoissonDistributionList() {
		mPdList = new ArrayListSortedComparable<PoissonDistributionSortable>();
	}
	
	// ========================================================================
	/** Returns whether this list already possessed a Poisson Distribution with the given mean. */
	public boolean registerMean(double mean) {
		PoissonDistributionSortable pds = new PoissonDistributionSortable(mean);
		return mPdList.add(pds);		
	}
	
	// ========================================================================
	/** Returns the mean of the Poisson distrubtion with the highest PDF for the value. */
	public double getMeanOfMostLikelyDistribution(int value, PrimitiveWrapper.WDouble probability) {
		if (mPdList.isEmpty()) return 0;
		
		int indexOfMax = 0;
		double maxPDF  = mPdList.get(indexOfMax).probability(value);
		
		for (int i = 1; i < mPdList.size(); i++) {
			double pdf = mPdList.get(i).probability(value);
			if (pdf > maxPDF) {
				maxPDF = pdf;
				indexOfMax = i;
			} else {
				break;  // since the distributions are sorted by mean
			}
		}
		
		probability.mDouble = maxPDF;
		return mPdList.get(indexOfMax).getMean(); 
	}
	
	// ========================================================================
	public void print(PrintStream out) {
		for (PoissonDistributionSortable pds : mPdList) {
			out.printf("%g\t", pds.getMean());
		}
		out.println();
	}
	
	// ========================================================================
	
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		PoissonDistributionList pdList = new PoissonDistributionList();
		pdList.registerMean(2);		
		pdList.registerMean(100);
		pdList.registerMean(10);
		pdList.registerMean(100);
		pdList.registerMean(15);
		pdList.print(System.out);
		
		PrimitiveWrapper.WDouble highestProbability = new PrimitiveWrapper.WDouble(0);
		
		for (int i = 0; i < 100; i++) {
			System.out.println(i + "\t" + pdList.getMeanOfMostLikelyDistribution(i, highestProbability) + "\t" + highestProbability.mDouble);
		}
	}

}
