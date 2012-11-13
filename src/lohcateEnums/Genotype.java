package lohcateEnums;

/** This is an enumerated type to represent the genotypes. */

public enum Genotype {
	EnumHomozygous00(0, "Hom00"),
	EnumHeterozygous(1, "Het"),
	EnumHomozygous11(2, "Hom11"),
	EnumHemizygous10(3, "Hem10"),
	EnumHemizygous01(4, "Hem01"),
	EnumHomoDeletion(5, "Del"),	
	EnumInvalidGenotype(-1, "Invalid");

	// ============ MEMBER VARIABLES ==============
	final int mCode;
	final String mString;

	// ========================================================================
	/** Returns whether genotype is either (Hom00, Het, or Hom11). */
	public boolean isFullGenotype() {
		return ((getCode() >= EnumHomozygous00.getCode()) && (getCode() <= EnumHomozygous11.getCode()));
	}
	
	// ========================================================================
	/** Returns the number of valid genotypes. */
	public static int numValid() { return values().length - 1; }

	// ========================================================================
	private Genotype(int value, String commonName) {
		mCode = value;	
		mString = commonName;
	}
	
	// ========================================================================
	public int getCode() { return mCode; }

	// ========================================================================
	public static Genotype getGenotype(int code) {
		if (code >= 0 && code < values().length) {
			return values()[code];
		} else {
			return EnumInvalidGenotype;
		}
	}
	
	// ========================================================================
	public String toString() { return mString; }
}
