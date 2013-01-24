package genomeUtils;

public class GenomeConstants {

	// ========================================================================
	/** The default (copy neutral) copy number of a diploid human individual. */
	public static final int DefaultDiploidCopyNumber = 2;
	public static final int DefaultHaploidCopyNumber = 1;
	
	// ========================================================================
	public static boolean isAmplification(double copyNumberTotal) {
		return (copyNumberTotal > DefaultDiploidCopyNumber);
	}
	
	// ========================================================================
	public static boolean isDoubleAmplification(double copyNumberA, double copyNumberB) {
		return ((copyNumberA > DefaultHaploidCopyNumber) && (copyNumberB > DefaultHaploidCopyNumber));
	}

	// ========================================================================
	public static boolean isDeletion(double copyNumberTotal) {
		return (copyNumberTotal < DefaultDiploidCopyNumber);
	}
	
	// ========================================================================
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub

	}

}