package genomeEnums;

public enum TissueType {
	Normal,
	Tumor;
	
	// ============ MEMBER VARIABLES ==============
	public final int mCode;
	
	// ========================================================================
	/** Returns the number of valid tissue types. */
	public static int numValid() { return values().length; }
	
	// ========================================================================
	private TissueType() {
		mCode = ordinal();	
	}
	
	// ========================================================================
	public int getCode() { return mCode; }
}
