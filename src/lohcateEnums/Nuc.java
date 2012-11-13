package lohcateEnums;

public enum Nuc implements Comparable<Nuc> {
	A(0, 'A'),
	C(1, 'C'),
	G(2, 'G'),
	T(3, 'T'),
	N(4, 'N');
	
	// The following arrays are created for speed.
	public static final Nuc[] ValidNucs = new Nuc[] { A, C, G, T}; 
	public static final Nuc[] Nucs = values();
	
	// ============ MEMBER VARIABLES ==============
	public final byte mCode;
	final char mChar;
	final String mString;
	
	// ========================================================================
	/** Returns the number of valid allele types. */
	public static int numValid() { return values().length - 1; }
	
	// ========================================================================
	/** Returns the total number of allele types. */
	public static int numTotal() { return values().length; }
	
	// ========================================================================
	/** Given a code, this returns the allele. */
	public static Nuc getAllele(long code) { return Nucs[(int) code]; }
	
	// ========================================================================
	public static Nuc getNuc(char ch) {
		Nuc theNuc = getNucUnsafe(ch);
		if (theNuc == null) {
			System.err.println("ERROR: Nuc.getNuc(): Invalid nucleotide!: (" + ch + ")");
			System.exit(-1);
		}
		return theNuc;
	}
	
	// ========================================================================
	public static Nuc getNucUnsafe(char ch) {
		switch(ch) {
		case 'A': case 'a': return Nuc.A;
		case 'C': case 'c': return Nuc.C;
		case 'G': case 'g': return Nuc.G;
		case 'T': case 't': return Nuc.T;
		case 'N': case 'n': return Nuc.N;
		default: return null;
		}
	}
	
	// ========================================================================
	private Nuc(int code, char charValue) {
		mCode = (byte) code;
		mChar = charValue;
		mString = "" + charValue;
	}
	
	// ========================================================================
	public byte getCode() { return mCode; }
	
	// ========================================================================
	public char getChar() { return mChar; }
	
	// ========================================================================
	public String getString() { return mString; } 
	
	// ========================================================================
	/** Returns whether this is the null allele or not. */
	public boolean isNull() { return (this == N); }
	
	// ========================================================================
	/** Returns whether this allele is not the null allele */
	public boolean isValid() { return (this != N); }

	// ========================================================================
	/** Returns whether a character is a valid nucleotide. */
	public static boolean isValid(char ch) {
		Nuc theNuc = getNucUnsafe(ch);
		return ((theNuc != null) && (theNuc.isValid()));
	}
	
	// ========================================================================
	/** Returns whether this allele is a G or C allele. */
	public boolean isGC() { return (this == Nuc.G || this == Nuc.C); }
	
	// ========================================================================
	/** Returns the complement of this nucleotide. */
	public Nuc getComplement() {
		if (isNull()) return this;
		
		return getAllele(Nuc.numValid() - 1 - getCode());		
	}

	// ========================================================================
	public static char getComplementChar(char theNucChar) {
		Nuc theNuc = getNuc(theNucChar);
		return theNuc.getComplement().getChar();
	}
	
	// ========================================================================
	/** Returns whether this allele is identical to either allele1 or allele2. */
	public boolean matchesTwoAlleles(Nuc allele1, Nuc allele2) {
		return ((this == allele1) || (this == allele2));
	}
}

