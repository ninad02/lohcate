package genomeUtils;

public class GenotypeUtils {

	// ========================================================================
	/** Given an rsID string, this returns the integer from the rs ID. */
	public static int getNumberFromRsId(String rsIdString) {
		rsIdString = rsIdString.trim();
		int index = rsIdString.indexOf(GenotypeUtils.RsPrefix);
		if (index >= 0) {
			return Integer.parseInt(rsIdString.substring(index + GenotypeUtils.RsPrefix.length())); 
		} else {
			return -1;
		}
	}

	public static final String RsPrefix = "rs";

	/** Given an rs number, this returns the rs# string form. */
	public static String getRsIdFromNumber(int rsNum) { return "rs" + rsNum; }
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub

	}

}
