package nutils;

import java.util.Arrays;

/** This class represents the contingency table. */

public class ContingencyTable {

	// ============================================================================================
	public static final NullaryClassFactory<ContingencyTable> ClassFactory = new NullaryClassFactory<ContingencyTable>(ContingencyTable.class);
	
	// ============================================================================================
	/** Define the enumerated type for a contingency table. */
	public static enum ContingencyTableValue {
		TruePositive(0),
		FalseNegative(1),
		FalsePositive(2),
		TrueNegative(3);
		
		int mIdentifier;
		
		ContingencyTableValue(int identifier) {
			mIdentifier = identifier;
		}
		
		public int getIdentifier() { return mIdentifier; }
		
		public boolean equalsContingencyTableValue(ContingencyTableValue rhs) {
			return (mIdentifier == rhs.mIdentifier);
		}
	}
	
	// ============================================================================================
	// Now define the variables for the class. 	
	int[] mCounts;
	double mPrevalence;
	
	// ============================================================================================
	/** Constructor. */
	public ContingencyTable() {
		mCounts = new int[ContingencyTableValue.values().length];
		clear();
	}

	// ============================================================================================
	/** Clear the table. */
	public void clear() {
		Arrays.fill(mCounts, 0);
	}
	
	// ============================================================================================
	/** Returns the value of an element in the table. */
	public int getCount(ContingencyTableValue ctv) {
		return mCounts[ctv.getIdentifier()];
	}
	
	// ============================================================================================
	/** Increments the count of the element. */
	public void incrCount(ContingencyTableValue ctv) {
		mCounts[ctv.getIdentifier()]++;
	}
	
	// ============================================================================================
	/** Increments the count of the element. */
	public void incrCount(ContingencyTableValue ctv, int incrAmount) {
		mCounts[ctv.getIdentifier()] += incrAmount;
	}

	
	// ============================================================================================
	/** Sets the count of the element. */
	public void setCount(ContingencyTableValue ctv, int theCount) {
		mCounts[ctv.getIdentifier()] = theCount;
	}

	// ============================================================================================
	/** Returns the truth positive total. */
	public int getTotalTruthPositive() {
		return getCount(ContingencyTableValue.TruePositive) + getCount(ContingencyTableValue.FalseNegative);
	}
	
	// ============================================================================================
	/** Returns the truth negative total. */
	public int getTotalTruthNegative() {
		return getCount(ContingencyTableValue.FalsePositive) + getCount(ContingencyTableValue.TrueNegative);
	}
	
	// ============================================================================================
	/** Returns the sensitivity. */
	public double getSensitivity() {
		return (double) getCount(ContingencyTableValue.TruePositive) / ((double) getTotalTruthPositive());
	}

	// ============================================================================================
	/** Returns the specificity. */
	public double getSpecificity() { 
		return (double) getCount(ContingencyTableValue.TrueNegative) / ((double) getTotalTruthNegative());
	}
	
	// ============================================================================================
	/** Returns the test positive total. */
	public int getTotalTestPositive() {
		return getCount(ContingencyTableValue.TruePositive) + getCount(ContingencyTableValue.FalsePositive);
	}
	
	// ============================================================================================
	/** Returns the test negative total. */
	public int getTotalTestNegative() {
		return getCount(ContingencyTableValue.FalseNegative) + getCount(ContingencyTableValue.TrueNegative);
	}
	
	// ============================================================================================
	/** Returns the positive predictive value. */
	public double getPPV() {
		return (double) getCount(ContingencyTableValue.TruePositive) / ((double) getTotalTestPositive());
	}
	
	// ============================================================================================
	/** Returns the negative predictive value. */
	public double getNPV() {
		return (double) getCount(ContingencyTableValue.TrueNegative) / ((double) getTotalTestNegative());
	}
	
	// ============================================================================================
	/** REturns the F-measure. */
	public double getFMeasure() {
		double precision = getPPV();
		double recall    = getSensitivity();
		return (2 * ((precision * recall) / (precision + recall)));
	}
	
	// ============================================================================================
	// Tests the contingency table
	private static void TestTable() {
		ContingencyTable ct = new ContingencyTable();
		ct.setCount(ContingencyTableValue.TruePositive, 171);
		ct.setCount(ContingencyTableValue.TrueNegative, 243);
		ct.setCount(ContingencyTableValue.FalsePositive, 23);
		ct.setCount(ContingencyTableValue.FalseNegative, 45);
		System.out.println("Sens: " + ct.getSensitivity());
		System.out.println("Spec: " + ct.getSpecificity());
		System.out.println("PPV:  " + ct.getPPV());
		System.out.println("NPV:  " + ct.getNPV());				
	}
	
	// ============================================================================================
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		TestTable();
	}

}


