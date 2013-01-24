package genomeEnums;

public enum Phase {
	p0(0),
	p1(1);
	
	// ============ MEMBER VARIABLES ==============
	public final int mCode;
	final String mString;
	
	// ========================================================================
	/** Returns the number of valid allele types. */
	public static int numValid() { return values().length; }
	
	// ========================================================================
	private Phase(int value) {
		mCode = value;	
		mString = "" + mCode;
	}
	
	// ========================================================================
	public int getCode() { return mCode; }
	
	// ========================================================================
	public String toString() { return mString; }
	
	// ========================================================================
	public static Phase getPhase(int code) { return values()[code]; }

}
