package genomeEnums;

import java.util.ArrayList;


public enum Chrom {
	c0(0, 0),
	c1(248_956_422, 125000000),
	c2(242_193_529, 93300000),
	c3(198_295_559, 91000000),
	c4(190_214_555, 50400000),
	c5(181_538_259, 48400000),
	c6(170_805_979, 61000000),
	c7(159_345_973, 59900000),
	c8(145_138_636, 45600000),
	c9(138_394_717, 49000000),
	c10(133_797_422, 40200000),
	c11(135_086_622, 53700000),
	c12(133_275_309, 35800000),
	c13(114_364_328, 17900000),
	c14(107_043_718, 17600000),
	c15(101_991_189, 19000000),
	c16(90_338_345, 36600000),
	c17(83_257_441, 24000000),
	c18(80_373_285, 17200000),
	c19(58_617_616, 26500000),
	c20(64_444_167, 27500000),
	c21(46_709_983, 13200000),
	c22(50_818_468, 14700000),
	cX(156_040_895, 60600000),
	cY(57_227_415, 12_500_000),
	cM(0, 0);
	
	public static final int IndexAutosomalStart = 1;
	public static final int IndexAutosomalEnd   = 22;
	public static final int IndexSexChromStart  = 23;
	public static final int IndexSexChromEnd    = 24;
	public static final int NumAutosomes = 22;
	public static final ArrayList<Chrom> Autosomes = createAutosomeArray();
	public static final int MaxPositionOnLongestChrom = 0x0FFFFFFF;  // 28 bit mask

	public static final String ChromPrefix_chr_ = "chr_";
	public static final String ChromPrefix_chr  = "chr";
	public static final String ChromPrefix_c    = "c";
	
	// ========================================================================
	private static ArrayList<Chrom> createAutosomeArray() {
		ArrayList<Chrom> autosomes = new ArrayList<Chrom>(NumAutosomes);
		for (int i = 0; i < NumAutosomes; i++) {
			autosomes.add(Chrom.values()[i + IndexAutosomalStart]);
		}
		return autosomes;
	}

	// ========================================================================
	// MEMBER VARIABLE
	protected String mName;
	protected int mLength;
	protected int mLengthArmP; 
	protected long mGenomeWidePositionStart;
	
	// ========================================================================
	private Chrom(int length, int lengthPArm) {
		mLength = length;
		mLengthArmP = lengthPArm;
		mGenomeWidePositionStart = -1;
		mName = name().substring(1);
	}

	// ========================================================================
	public String getName() { return mName; }
	
	// ========================================================================
	public byte getCode() { return (byte) ordinal(); }

	// ========================================================================
	public int getLength() { return mLength; }

	// ========================================================================
	public int getLengthOfPArm() { return mLengthArmP; }
	
	// ========================================================================
	public Arm getArm(int positionOnChrom) {
		return (positionOnChrom >= getLengthOfPArm() ? Arm.q : Arm.p); 
	}

	// ========================================================================
	private long calculateGenomeWidePositionStart() {
		long genomeWidePositionStart = 0;
		for (Chrom chrom : Chrom.values()) {
			if (chrom == this) break;			
			genomeWidePositionStart += chrom.getLength();
		}
		return genomeWidePositionStart;
	}
	
	// ========================================================================
	public long getGenomeWidePositionStart() {
		if (mGenomeWidePositionStart < 0) {
			mGenomeWidePositionStart = calculateGenomeWidePositionStart();
		}
		return mGenomeWidePositionStart;
	}

	// ========================================================================
	public long calculateGenomeWidePositionEndOfArmP() {
		return getGenomeWidePositionStart() + mLengthArmP;
	}
	
	// ========================================================================
	public long calculateGenomeWidePositionEnd() {
		return getGenomeWidePositionStart() + getLength();
	}
	
	// ========================================================================
	public long calculateGenomeWidePositionMidpoint() {
		return getGenomeWidePositionStart() + (mLength >> 1);
	}
	
	// ========================================================================
	public static Chrom getChrom(String chromStr) {
		chromStr = chromStr.trim().toLowerCase();
		String subStr = chromStr;		
		
		if (chromStr.indexOf(ChromPrefix_chr_) == 0) {
			subStr = chromStr.substring(ChromPrefix_chr_.length());	
			
		} else if (chromStr.indexOf(ChromPrefix_chr) == 0) {
			subStr = chromStr.substring(ChromPrefix_chr.length());
			
		} else if (chromStr.indexOf(ChromPrefix_c) == 0) {
			subStr = chromStr.substring(ChromPrefix_c.length());
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
	public static Chrom getChrom(byte chromNum) { 
		return ((chromNum < 0 || chromNum >= values().length) ? null : values()[chromNum]); 
	}

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
	
	// ========================================================================
	public Chrom getNextChrom() {
		return Chrom.getChrom((byte) (ordinal() + 1));
	}
	
	// ========================================================================
	public static enum Arm { p, q; }
}
