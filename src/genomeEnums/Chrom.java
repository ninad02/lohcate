package genomeEnums;

import java.util.ArrayList;

public enum Chrom {
	c0,
	c1,
	c2,
	c3,
	c4,
	c5,
	c6,
	c7,
	c8,
	c9,
	c10,
	c11,
	c12,
	c13,
	c14,
	c15,
	c16,
	c17,
	c18,
	c19,
	c20,
	c21,
	c22,
	cX,
	cY,
	cM;
	
	public static final int IndexAutosomalStart = 1;
	public static final int IndexAutosomalEnd   = 22;
	public static final int IndexSexChromStart  = 23;
	public static final int IndexSexChromEnd    = 24;
	public static final int NumAutosomes = 22;
	public static final ArrayList<Chrom> Autosomes = createAutosomeArray();

	// ========================================================================
	private static ArrayList<Chrom> createAutosomeArray() {
		ArrayList<Chrom> autosomes = new ArrayList<Chrom>(NumAutosomes);
		for (int i = 0; i < NumAutosomes; i++) {
			autosomes.add(Chrom.values()[i + IndexAutosomalStart]);
		}
		return autosomes;
	}
	
	// ========================================================================
	private Chrom() {}
	
	// ========================================================================
	public byte getCode() { return (byte) ordinal(); }

	// ========================================================================
	public static Chrom getChrom(String chromStr) {
		chromStr = chromStr.trim().toLowerCase();
		String subStr = chromStr;		
		if (chromStr.indexOf("chr_") == 0) {
			subStr = chromStr.substring("chr_".length());			
		} else if (chromStr.indexOf("chr") == 0) {
			subStr = chromStr.substring("chr".length());
		} else if (chromStr.indexOf("c") == 0) {
			subStr = chromStr.substring("c".length());
		}
		
		if (subStr.charAt(0) == 'x') {
			return cX;
		} else if (subStr.charAt(0) == 'y') {
			return cY;		
		} else if (subStr.charAt(0) == 'm') {
			return cM;
		} else {
			return getChrom(Byte.parseByte(subStr));
		}
	}
	
	// ========================================================================
	public static Chrom getChrom(byte chromNum) { return values()[chromNum]; }

	// ========================================================================
	public boolean isAutosomal() { 
		int code = getCode();
		return (code >= IndexAutosomalStart) && (code <= IndexAutosomalEnd); 
	}
	
	// ========================================================================
	public boolean isSexChromosome() {
		int code = getCode();
		return (code >= IndexSexChromStart) && (code <= IndexSexChromEnd);
	}
	
	// ========================================================================
	public boolean isInvalid() { return (this == c0); }

	// =======================================================================
	public static Chrom getInvalidChrom() { return c0; }
	
	// ========================================================================
	public static int getNumAutosomal() { return IndexAutosomalEnd - IndexAutosomalStart + 1; }
	
	// ========================================================================
	public static int getNumSexChromosomes() { return IndexSexChromEnd - IndexSexChromStart + 1; }
	
	// ========================================================================
	public boolean lessThan(Chrom rhs) {
		return (this.getCode() < rhs.getCode());
	}
	
	// ========================================================================
	public boolean greaterThan(Chrom rhs) {
		return (this.getCode() > rhs.getCode());
	}
}
