package genomeUtils;

import genomeEnums.Chrom;
import genomeEnums.Genotype;
import genomeEnums.Nuc;
import genomeUtils.GenotypeUtils;

import java.io.BufferedWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.Iterator;

import nutils.ArrayUtils;
import nutils.CompareUtils;
import nutils.IOUtils;


import com.carrotsearch.hppc.IntArrayList;
import com.carrotsearch.hppc.LongArrayList;

import shared.Utils;

import nutils.NumberUtils;

public class SNVMap {

	public static final int NucleotideMapUnknownAlleleIndex = 0;
	public static final int IndexAlleleA = 0;
	public static final int IndexAlleleB = 1;
	public static final int IndexAlleleBoth = 2;
	public static final int NumAllelesInCallSet = 2;
	
	ArrayList<LongArrayList> mSNVsByPosition;
	ArrayList<LongArrayList> mSNVsByRsId;	
	
	public static final int DefaultRsId = 0;
	public static final boolean DefaultStrand = true;
	
	public SNVMap() {
		mSNVsByPosition = new ArrayList<LongArrayList>();
		mSNVsByRsId     = new ArrayList<LongArrayList>();	
		
		for (Chrom chrom : Chrom.values()) {
			System.out.println(chrom);
			mSNVsByPosition.add(new LongArrayList());
			mSNVsByRsId.add(new LongArrayList());
		}
	}
	
	public static boolean isMissingAllele(char allele) { return Nuc.isValid(allele); }			
	
	public void printMe(boolean orderedByPosition, String outFilename) {
		BufferedWriter out = IOUtils.getBufferedWriter(outFilename);
		ArrayList<LongArrayList> tableToPrint = orderedByPosition ? mSNVsByPosition : mSNVsByRsId;
		printMe(tableToPrint, out, orderedByPosition);
		IOUtils.closeBufferedWriter(out);
	}
	
	private static void printMe(ArrayList<LongArrayList> tableToPrint, BufferedWriter out, boolean orderedByPosition) {
		Nuc[] nucs = new Nuc[NumAllelesInCallSet];
		StringBuilder sb = new StringBuilder(2048);
		int novelCounter = 0;
		for (Chrom chrom : Chrom.values()) {
			LongArrayList listOfSites = tableToPrint.get(chrom.ordinal());
			for (int i = 0; i < listOfSites.size(); i++) {
				long compactUnit = listOfSites.get(i);
				sb.setLength(0);
				
				int rsID = extractRsIdFromCompactForm(compactUnit, !orderedByPosition);
				if (rsID == 0) {
					sb.append("novel_").append(++novelCounter);
				} else {
					sb.append(GenotypeUtils.RsPrefix).append(rsID);
				}
				
				sb.append(Utils.TabStr).append(chrom.ordinal());
				sb.append(Utils.TabStr).append(extractPositionFromCompactForm(compactUnit, !orderedByPosition));
				
				extractNucleotideIDsFromCompactForm(compactUnit, nucs);
				sb.append(Utils.TabStr).append(nucs[0]).append(Utils.TabStr).append(nucs[1]);
				sb.append(Utils.TabStr).append(Long.toBinaryString(compactUnit));

				IOUtils.writeToBufferedWriter(out, sb.toString(), true);
			}
		}
	}
	
	// ========================================================================	
	/** Given a chromosome number, position, rsId, two nucleotides, and the
	 *  strand, all of which represent a SNP, this stores the SNP in memory.
	 *  Returns whether the SNV was inserted or not (e.g. may return false
	 *  if site alreay existed in list).
	 */
	public Boolean registerSNV(Chrom chrom, int position, int rsId, Nuc nuc1, Nuc nuc2, boolean strand, boolean uniqueOnly) {
		CompareUtils.ensureTrue(rsId >= 0, "ERROR: AffySNPMap.registerSNP(): rsId cannot be < 0!");
		if (!strand) {
			nuc1 = nuc1.getComplement();
			nuc2 = nuc2.getComplement();
		}
		
		int resultIndex = getIndexOfPositionInMap(chrom, position);
		if (!uniqueOnly || (resultIndex < 0)) {
			int insertionIndex = (resultIndex < 0) ? -(resultIndex + 1) : resultIndex;
			long compactUnit = compactSNPInfo(position, rsId, nuc1, nuc2, false);
			mSNVsByPosition.get(chrom.ordinal()).insert(insertionIndex, compactUnit);
			return Boolean.TRUE;
		} else {
			return Boolean.FALSE;
		}
	}
	
	/** Given a chromosome number, position, rsId, two nucleotides, and the
	 *  indicator of whether it is Illumina Top or Bottom -- all of which 
	 *  represent a SNP, this stores the SNP in memory.
	 *  TODO - UNFINISHED FUNCTION
	 */
	public void registerSNPinIllumina(byte chromNum, int position, int rsId, Nuc nuc1, Nuc nuc2, boolean isIlluminaTop) {
		// These are the rules.  
		if (isIlluminaTop) {
			CompareUtils.ensureTrue(nuc1  == Nuc.A, "ERROR: SNVMap.registerSNP(): nuc1 is " + nuc1 + "instead of A");
			CompareUtils.ensureTrue((nuc2 == Nuc.C || nuc2 == Nuc.G), "ERROR: AffySNPMap.registerSNP(): nuc2 is " + nuc2 + "instead of C or G");
			// do the actual SNP registration
		} else {
			CompareUtils.ensureTrue(nuc2  == Nuc.T, "ERROR: SNVMap.registerSNP(): nuc2 is " + nuc2 + "instead of T");
			CompareUtils.ensureTrue((nuc1 == Nuc.C || nuc1 == Nuc.G), "ERROR: AffySNPMap.registerSNP(): nuc1 is " + nuc1 + "instead of C or G");		
			// do the actual SNP registration
		}
	}
	
	/** Compacts the information into a single unit. */
	private static final long ShiftNuc = 3;
	private static final int MaskNuc = 0x07;
	private static final int MaskPosition = 0x7FFFFFFF;
	private static final int MaskRsId     = 0x07FFFFFF;
	
	private static long compactSNPInfo(int position, int rsId, Nuc nuc1, Nuc nuc2, boolean makeRsIDMostSignificant) {		
		int bitsToShiftForPosition = makeRsIDMostSignificant ? 6 : 33;
		int bitsToShiftForRsId = makeRsIDMostSignificant ? 37 : 6;
		
		long rV = 0L | 				   
				  (long) (((long) (position & MaskPosition)) << bitsToShiftForPosition) |				  
				  (long) (((long) (rsId     & MaskRsId))     << bitsToShiftForRsId) |
				  (long) (((long) (nuc1.getCode() & MaskNuc)) << ShiftNuc) |
				  (long)          (nuc2.getCode() & MaskNuc);
		return rV;
	}

	
	/** Extracts the position from the compact form. */
	private static int extractPositionFromCompactForm(long compactUnit, boolean makeRsIDMostSignificant) {
		return ((int) (compactUnit >>> (makeRsIDMostSignificant ? 6 : 33))) & MaskPosition;
	}
	
	
	/** Extracts the RS ID from the compact form. */
	private static int extractRsIdFromCompactForm(long compactUnit, boolean makeRsIDMostSignificant) {
		return ((int) (compactUnit >>> (makeRsIDMostSignificant ? 37 : 6))) & MaskRsId;
	}
	
	/** Extracts the nucleotide IDs as integers from the compact form. */
	private static Nuc[] extractNucleotideIDsFromCompactForm(long compactUnit, Nuc[] nucleotideAlleles) {		
		for (int i = NumAllelesInCallSet - 1; i >= 0; i--) {
			nucleotideAlleles[i] = Nuc.getAllele((byte) (compactUnit & MaskNuc));
			compactUnit = (long) (compactUnit >>> ShiftNuc);
		}
		return nucleotideAlleles;
	}
	
	private static final ExtractorFromCompactUnitRsID_rsIDFirst extractorRsID_rsIDFirst = new ExtractorFromCompactUnitRsID_rsIDFirst();
	private static final ExtractorFromCompactUnitPosition_PositionFirst extractorPosition_PositionFirst = new ExtractorFromCompactUnitPosition_PositionFirst();
	
	private static interface ExtractorFromCompactUnit {
		public long extractValue(long compactUnit);
	}
	
	private static class ExtractorFromCompactUnitRsID_rsIDFirst implements ExtractorFromCompactUnit {
		public long extractValue(long compactUnit) { return extractRsIdFromCompactForm(compactUnit, true); }
	}
	
	private static class ExtractorFromCompactUnitPosition_PositionFirst implements ExtractorFromCompactUnit {
		public long extractValue(long compactUnit) { return extractPositionFromCompactForm(compactUnit, false); }
	}
	
	// ========================================================================
	// Returns the relevant array of longs.  Parameters specify the chromosome number
	// in question as well as whether the array should be ordered by position (true)
	// or rsId.
	private LongArrayList getArrayForChrom(Chrom chrom, boolean orderedByPosition) {
		return (orderedByPosition ? mSNVsByPosition.get(chrom.ordinal()) : mSNVsByRsId.get(chrom.ordinal()));				
	}
	
	// ========================================================================
	
	/** Retrieves the nucleotides as byte variables for a particular chromosome and position. */
	public Nuc[] getNucleotidesInMap(Chrom chrom, int position) {
		Nuc[] rV = new Nuc[NumAllelesInCallSet];
		return getNucleotidesInMap(chrom, position, rV);
	}

	/** Retrieves the nucleotides as byte variables for a particular chromosome and position. */
	public Nuc[] getNucleotidesInMap(Chrom chrom, int position, Nuc[] nucleotideAlleles) {
		LongArrayList arrayForChrom = getArrayForChrom(chrom, true);
		int resultIndex = binarySearchValue(position, arrayForChrom, extractorPosition_PositionFirst);
		if (resultIndex >= 0) {
			return extractNucleotideIDsFromCompactForm(arrayForChrom.get(resultIndex), nucleotideAlleles);
		} else { 
			return null;
		}
	}
	
	/** Retrives the nucleotides as byte variables for a particular chromosome and rsID number. 
	 *  This runs in linear time. */
	public Nuc[] getNucleotidesInMapByRsID(final Chrom chrom, final int rsId) {		
		return getNucleotidesInMapByRsID(chrom, rsId, new Nuc[NumAllelesInCallSet]);
	}
	
	public Nuc[] getNucleotidesInMapByRsID(final Chrom chrom, final int rsId, Nuc[] nucleotideAlleles) {
		int resultIndex = getIndexOfRsIdInMapHelper(chrom, rsId);
		if (resultIndex >= 0) {
			LongArrayList arrayForChrom = getArrayForChrom(chrom, false);
			return extractNucleotideIDsFromCompactForm(arrayForChrom.get(resultIndex), nucleotideAlleles);
		} else {
			return null;
		}
	}
	
	public Nuc[] getNucleotidesInMapByRsID(final int rsId, Nuc[] nucleotideAlleles) {
		for (Chrom chrom : Chrom.values()) {
			Nuc[] result = getNucleotidesInMapByRsID(chrom, rsId, nucleotideAlleles);
			if (result != null) {
				return result;
			}
		}
		return null;
	}
	
	/** First tries to get the bytes by rsID.  If that doesn't work, then it tries for
	 *  the position
	 */
	public Nuc[] getNucleotidesInMapByRsIDOrPosition(final Chrom chrom, final int rsId, final int position, Nuc[] nucleotideAlleles) {
		Nuc[] result = getNucleotidesInMapByRsID(chrom, rsId, nucleotideAlleles);
		if (result == null) {
			result = getNucleotidesInMap(chrom, position, nucleotideAlleles);
		}
		return result;
	}
		
	public int getIndexOfPositionInMap(final Chrom chrom, int position) {
		LongArrayList arrayForChrom = getArrayForChrom(chrom, true);
		return binarySearchValue(position, arrayForChrom, extractorPosition_PositionFirst);
	}
	
	/** Returns the position at the specified index along the chromosomal positions. 
	 * @param chrom
	 * @param index
	 * @return the position if a valid index is specified.  -1 for an invalid index
	 */
	public int getPosition(final Chrom chrom, final int index) {
		LongArrayList arrayForChrom = getArrayForChrom(chrom, true);
		return ((index >= 0 && index < arrayForChrom.size()) ? extractPositionFromCompactForm(arrayForChrom.get(index), false) : -1); 
	}
	
	/** Returns the number of sites on a particular chromosome. */
	public int getNumSitesOnChromosome(final Chrom chrom) {
		return getArrayForChrom(chrom, true).size();
	}
	
	public int getIndexOfRsIdInMap(final Chrom chrom, final int rsId) {
		int resultIndex = getIndexOfRsIdInMapHelper(chrom, rsId);
		if (resultIndex >= 0) {
			long compactUnit = getArrayForChrom(chrom, false).get(resultIndex);
			int position = extractPositionFromCompactForm(compactUnit, true);
			return getIndexOfPositionInMap(chrom, position);
		}
		return resultIndex;
	}
	
	private int getIndexOfRsIdInMapHelper(final Chrom chrom, final int rsId) {
		LongArrayList arrayForChrom = getArrayForChrom(chrom, false);
		return binarySearchValue(rsId, arrayForChrom, extractorRsID_rsIDFirst);		
	}
	
	public int[] getPositionsBetweenPositionsInMap(final Chrom chrom, int positionStart, int positionEnd, boolean posOrRsID) {		
		int indexStart = getIndexOfPositionInMap(chrom, positionStart);
		int indexEnd =  getIndexOfPositionInMap(chrom, positionEnd);
	
		if ((indexStart >= 0) && (indexEnd >= indexStart)) {
			LongArrayList arrayForChrom = getArrayForChrom(chrom, true);
			int numPositions = indexEnd - indexStart + 1;
			int[] positions = new int[numPositions];
			for (int i = indexStart; i <= indexEnd; i++) {
				if (posOrRsID) {
					positions[i - indexStart] = extractPositionFromCompactForm(arrayForChrom.get(i), false);
				} else {
					positions[i - indexStart] = extractRsIdFromCompactForm(arrayForChrom.get(i), false);
				}
			}			
			return positions;
		} else {
			return null;
		}		
	}
	
	/** Retrieves the rsId as an integer for a particular chromosome and position. */
	public int getRsIdInMap(final Chrom chrom, int position) {
		LongArrayList arrayForChrom = getArrayForChrom(chrom, true);
		int resultIndex = binarySearchValue(position, arrayForChrom, extractorPosition_PositionFirst);
		if (resultIndex >= 0) {
			return extractRsIdFromCompactForm(arrayForChrom.get(resultIndex), false);
		} else { 
			return -1; //DefaultRsId;
		}
	}
	
	// ========================================================================	
	
	/** This method compares the allele-pairs at a SNP between two platforms.  The
	 *  first argument is a byte[] that represents the allele-pair from the first 
	 *  platform, and the second byte[] argument represents the allele-pair from
	 *  the second platform.  The third byte[] argument result[] is provided by the 
	 *  caller as an output buffer; its length should be 2.  
	 *  
	 *  This method compares the two allele-pairs to test whether they are the same.
	 *  It goes even as far as to take the complement of the alleles in the second
	 *  pair or flip them if need be.  If it takes the complement, then result[0] is
	 *  set to true.  If it flips them, then result[1] is set to true.  Otherwise,
	 *  the result[] array values are set to false.  
	 *  
	 *  The method returns true if the if the allele-pairs are the same, false otherwise.
	 */
	public static boolean compareAllelePairs(final Nuc[] pair1, final Nuc[] pair2, boolean[] result) {
		Arrays.fill(result, false);		
		
		// First, test whether they are the same
		if (Arrays.equals(pair1, pair2)) {
			return true;
		} else {
			// They are not the same just by themselves.  We will need a temporary 
			// array that will hold the values of pair2.  
			Nuc[] pair2Temp = new Nuc[pair2.length];
			
			// We first take the complement of pair2 and see whether that works.			
			for (int i = 0; i < pair2Temp.length; i++) {
				pair2Temp[i] = pair2[i].getComplement(); 						
			}
			
			if (Arrays.equals(pair1, pair2Temp)) {				
				result[0] = true;  // It means we needed to only take the complement
				return true;
			}
									
			// The complement alone did not work.  Now try flipping the complemented array.			
			ArrayUtils.reverseArray(pair2Temp);
			if (Arrays.equals(pair1, pair2Temp)) {
				result[0] = result[1] = true; // We needed to complement and flip
				return true;
			}
			
			// Complement and flipping did not work.  We now undo the complementing and 
			// see whether just flipping alone is sufficient.
			for (int i = 0; i < pair2Temp.length; i++) {
				pair2Temp[i] = pair2Temp[i].getComplement();
			}
			
			if (Arrays.equals(pair1, pair2Temp)) {				
				result[1] = true;	// We only needed to flip
				return true;
			}
		}
		
		return false;
	}
	
	/** Given a genotype code, this fills in the characters for that genotype 
	 *  given the mapping characters.  For heterozygous cases, the first mapping
	 *  char is written first, then the second (unless flipForHet is true)
	 * @param genotypeCode - A genotype code represented by HardyWeinbergCalculator
	 * @param nucsChar - The array that the genotype chars will be written to
	 * @param referenceVariantAlleles - The array that houses the reference nucleotide alleles for a marker
	 * @return - false if an invalid genotype code, true otherwise
	 */
	public static boolean getGenotypeAsChars(Genotype genotypeCode, Nuc genotypeAlleles[], final Nuc referenceVariantAlleles[], boolean flipForHet) {
		switch (genotypeCode) {
		case EnumHomozygous00:
			genotypeAlleles[0] = referenceVariantAlleles[0];
			genotypeAlleles[1] = referenceVariantAlleles[0];
			break;
		case EnumHomozygous11:
			genotypeAlleles[0] = referenceVariantAlleles[1];
			genotypeAlleles[1] = referenceVariantAlleles[1];
			break;
		case EnumHeterozygous:
			genotypeAlleles[0] = flipForHet ? referenceVariantAlleles[1] : referenceVariantAlleles[0];
			genotypeAlleles[1] = flipForHet ? referenceVariantAlleles[0] : referenceVariantAlleles[1];
			break;
		case EnumHemizygous10:
			genotypeAlleles[0] = referenceVariantAlleles[0];
			genotypeAlleles[1] = Nuc.N;
			break;
		case EnumHemizygous01:
			genotypeAlleles[0] = referenceVariantAlleles[1];
			genotypeAlleles[1] = Nuc.N;
			break;
		default: 
			genotypeAlleles[0] = Nuc.N;
			genotypeAlleles[1] = Nuc.N;
			return false;
		}
		return true;
	}
	
	/** Checks whether the given genotype concords with the alleles in the map.  If one
	 *  of the alleles is invalid, then null is returned.  Else a boolean value is returned.
	 */
	public static Boolean testForConcordanceWithMap(final Nuc[] alleles, Nuc[] referenceVariantAlleles) {
		for (Nuc allele : alleles) {		
			final Boolean result = testForConcordanceWithMap(allele, referenceVariantAlleles);
			if (result == null) {
				return null;
			} else if (result == Boolean.FALSE) {
				return Boolean.FALSE;
			}			
		}
		return Boolean.TRUE;
	}
	
	/** Checks whether a given char concords with one of the alleles in the map.  If the
	 *  allele is an invalid allele, null is returned.  Else a boolean value is returned.
	 */
	public static Boolean testForConcordanceWithMap(final Nuc allele, Nuc[] referenceVariantAlleles) {
		if (allele.isNull()) {				
			return null;
		} else {
			for (Nuc otherAllele : referenceVariantAlleles) {
				if (allele == otherAllele) return Boolean.TRUE;
			}
		}
		return Boolean.FALSE;
	}
	
	/** Given a sample index and integer array, this performs a binary search
	 *  for the sample index on the integer array.  If the sample index is 
	 *  found, an index >= 0 is returned, else (-(insertion point) - 1) is returned. */
	private static int binarySearchValue(final int value, LongArrayList longArray, ExtractorFromCompactUnit extractor) {
		int lowerIndex = 0;
		int upperIndex = longArray.size() - 1;
		int midIndex = 0;
		
		// If empty array, we return with the insertion point at index 0
		if (upperIndex < 0) return -1;
		
		// Shortcut Lower: Test if the value is <= than the lowest value the array.
		// If so, then we don't perform the binary search, and we simply return.
		long valueLowerIndex = (int) extractor.extractValue(longArray.get(lowerIndex));
		if (value < valueLowerIndex) {
			return ((-lowerIndex) - 1);			
		} else if (value == valueLowerIndex) {
			return lowerIndex;
		} 

		// Shortcut Upper: Test if the value is >= the high value in the array.
		// If so, then we don't perform the binary search, and we simply return.
		long valueUpperIndex = (int) extractor.extractValue(longArray.get(upperIndex));
		if (value > valueUpperIndex) {
			return (-(upperIndex + 1) - 1);
		} else if (value == valueUpperIndex) {
			return upperIndex;
		}
				
		int valueAtMidIndex;
		while (lowerIndex <= upperIndex) {
			midIndex = (lowerIndex + upperIndex) >>> 1;  // right-shift by 1 to divide by 2; 
			valueAtMidIndex = (int) extractor.extractValue(longArray.get(midIndex));
			
			if (value == valueAtMidIndex) {
				return midIndex;
			} else if (value > valueAtMidIndex) {
				lowerIndex = midIndex + 1;
			} else {
				upperIndex = midIndex - 1;
			}
		}
		
		return -(lowerIndex + 1);
	}
	
	/** Given an array of bytes and a nucleotide, this returns the index of the 
	 *  nucleotide in the byte array. */
	public static int getIndexOfNucleotideInByteArray(final Nuc[] theArray, final Nuc nuc) {
		for (int i = 0; i < theArray.length; i++) {
			if (nuc == theArray[i]) {
				return i;
			}
		}
		return -1;
	}
	
	// =========================== TESTING FUNCTIONS ==========================
	
	/*
	public static void TestAffySNPMapWhole() {
		String inFilename = IOUtils.pathConcat(new String[]{"..", "..", "..", "..", "TCGA", "GenomeWideSNP_6.na24.annot.csv.All2.autosomes.csv"});
		AffySNPMap asm = CopyNumberParse.parseAffySNPFileHelper(inFilename);
		byte[] nucs = asm.getNucleotidesInMapByRsID((byte) 1, 34408665);
		asm.printMe(false, "AffyMapByRsID.txt");
		asm.printMe(true, "AffyMapByPosition.txt");
	}
	*/
	
	public static void TestSNVMap_Robust() {
		int numTrials = 2000000;
		
		ArrayList<Chrom> chromList = new ArrayList<Chrom>(numTrials);
		IntArrayList positionList = new IntArrayList(numTrials);
		IntArrayList rsIDList     = new IntArrayList(numTrials);
		ArrayList<Boolean> strandList = new ArrayList<Boolean>(numTrials);
		ArrayList<Nuc> nuc1List = new ArrayList<Nuc>(numTrials);
		ArrayList<Nuc> nuc2List = new ArrayList<Nuc>(numTrials);
		
		SNVMap snvMap = new SNVMap();
		
		for (int trial = 0; trial < numTrials; trial++) {
			
			if (trial % 100000 == 0) System.out.println("Trial:\t" + trial);
			
			Chrom chrom = Chrom.getChrom((byte) NumberUtils.getRandomInteger(0, Chrom.values().length - 1));
			int rsID = NumberUtils.getRandomInteger(0, 67108863);
			int position = trial + 1; //NumberUtils.getRandomInteger(0, 2147483647);
			Boolean strand = NumberUtils.getRandomBit() ? Boolean.TRUE : Boolean.FALSE;			
			Nuc nuc1 = Nuc.getAllele(NumberUtils.getRandomInteger(Nuc.A.getCode(), Nuc.N.getCode()));
			Nuc nuc2 = Nuc.getAllele(NumberUtils.getRandomInteger(Nuc.A.getCode(), Nuc.N.getCode()));
					
			Boolean addResult = snvMap.registerSNV(chrom, position, rsID, nuc1, nuc2, strand, true);
			
			if (addResult) {
				chromList.add(chrom);
				positionList.add(position);
				strandList.add(strand);
				rsIDList.add(rsID);
				nuc1List.add(nuc1);
				nuc2List.add(nuc2);
			}
		}
		
		Nuc[] extractedNucs = new Nuc[2];
		
		// Now test whether stored properly
		for (int trial = 0; trial < chromList.size(); trial++) {			
			int extractedRsId = snvMap.getRsIdInMap(chromList.get(trial), positionList.get(trial));
			if (extractedRsId != rsIDList.get(trial)) {
				String errorString = "ERROR: Correct rs-id not pulled!\n" + rsIDList.get(trial) + "\t" + extractedRsId;
				CompareUtils.ensureTrue(false, errorString);
			}
			
			snvMap.getNucleotidesInMap(chromList.get(trial), positionList.get(trial), extractedNucs);
			if (strandList.get(trial)) {
				if (extractedNucs[0] != nuc1List.get(trial) || extractedNucs[1] != nuc2List.get(trial)) {
					String errorString = "ERROR: Correct nuc alleles not pulled on positive strand!";
					CompareUtils.ensureTrue(false, errorString);
				}
			} else {
				if (extractedNucs[0] != nuc1List.get(trial).getComplement() || extractedNucs[1] != nuc2List.get(trial).getComplement()) {
					String errorString = "ERROR: Correct nuc alleles not pulled on negative strand!";
					CompareUtils.ensureTrue(false, errorString);
				}
			}
		}
	}
	
	public static void TestAffySNPMap() {
		int numChromUsed = 4;
		int dummyBuffer = 3;
		int[] chromCount = new int[numChromUsed + dummyBuffer];
		Arrays.fill(chromCount, 0);
		
		// Fill the count arrays
		for (int i = 0; i < chromCount.length; i++) {
			chromCount[i] = i + 1;
		}
		
		// create the map
		SNVMap aMap = new SNVMap();
		
		// Get the bytes for the nucleotides
		// Now, register accordingly
		aMap.registerSNV(Chrom.c1, 100, 10000, Nuc.A, Nuc.G, false, true);
		
		aMap.registerSNV(Chrom.c2, 200, 11000, Nuc.C, Nuc.G, true, true);
		aMap.registerSNV(Chrom.c2, 100, 12000, Nuc.A, Nuc.G, true, true);
		
		aMap.registerSNV(Chrom.c3, 150, 13000, Nuc.C, Nuc.G, true, true);		
		aMap.registerSNV(Chrom.c3, 100, 14000, Nuc.A, Nuc.G, true, true);
		aMap.registerSNV(Chrom.c3, 200, 15000, Nuc.G, Nuc.G, true, true);
		
		aMap.registerSNV(Chrom.c4, 133, 16000, Nuc.C, Nuc.G, true, true);
		aMap.registerSNV(Chrom.c4, 100, 17000, Nuc.A, Nuc.G, true, true);
		aMap.registerSNV(Chrom.c4, 200, 18000, Nuc.T, Nuc.G, true, true);
		aMap.registerSNV(Chrom.c4, 166, 19000, Nuc.G, Nuc.G, false, true);		
		
		
		// Print the Affy structure
		aMap.printMe(true, "AffyMap.txt");
		
		// Now see that the elements exist
		//byte[] result;
		System.out.println(getStringForBytes(aMap.getNucleotidesInMap(Chrom.c4, 166)));
		System.out.println(getStringForBytes(aMap.getNucleotidesInMap(Chrom.c3, 150)));
		System.out.println(getStringForBytes(aMap.getNucleotidesInMap(Chrom.c4, 133)));
		System.out.println(getStringForBytes(aMap.getNucleotidesInMap(Chrom.c2, 200)));
		System.out.println(getStringForBytes(aMap.getNucleotidesInMap(Chrom.c1, 100)));
		System.out.println(getStringForBytes(aMap.getNucleotidesInMap(Chrom.c3, 100)));
		System.out.println(getStringForBytes(aMap.getNucleotidesInMap(Chrom.c3, 125)));
		
	}
	
	public static String getStringForBytes(Nuc[] theNucs) {
		if (theNucs == null) return "null";
		
		StringBuilder sb = new StringBuilder(1024);
		sb.append('(');
		for (int i = 0; i < theNucs.length; i++) {
			sb.append(theNucs[i]);
			if (i < (theNucs.length - 1)) {
				sb.append(", ");
			}
		}
		sb.append(')');
		return sb.toString();
	}
	
	public static void TestSNPCompact() {
		
		boolean makeRsIDMostSignificant = false;
		Nuc[] extractedNucs = new Nuc[2];
		Nuc[] assignedNucs = new Nuc[2];
		
		do {
			int numTrials = 10000000;
			for (int trial = 0; trial < numTrials; trial++) {
				int rsID = NumberUtils.getRandomInteger(0, 67108863);
				int position = NumberUtils.getRandomInteger(0, 2147483647);
				boolean strand = NumberUtils.getRandomBit();

				assignedNucs[0] = Nuc.getAllele(NumberUtils.getRandomInteger(Nuc.A.getCode(), Nuc.N.getCode()));
				assignedNucs[1] = Nuc.getAllele(NumberUtils.getRandomInteger(Nuc.A.getCode(), Nuc.N.getCode()));

				long compactUnit = compactSNPInfo(position, rsID, assignedNucs[0], assignedNucs[1], makeRsIDMostSignificant);
				
				int extractedPosition = extractPositionFromCompactForm(compactUnit, makeRsIDMostSignificant);
				int extractedRsID = extractRsIdFromCompactForm(compactUnit, makeRsIDMostSignificant);
				extractNucleotideIDsFromCompactForm(compactUnit, extractedNucs);
				
				if (extractedPosition != position) {
					String errorString = "ERROR: Mismatched positions!\nTrial:\t" + trial + "\n" + position + "\t" + extractedPosition;
					System.out.println("BinString = " + Long.toBinaryString(compactUnit));
					CompareUtils.ensureTrue(false, errorString);
				}
				
				CompareUtils.ensureTrue(extractedRsID == rsID, "ERROR: Mismatched rsIDs!");
				
				if (!Arrays.equals(extractedNucs, assignedNucs)) {
					String errorString = "ERROR: Mismatched Nucs!\n" + getStringForBytes(extractedNucs) + "\n" + getStringForBytes(assignedNucs);
					CompareUtils.ensureTrue(false, errorString);
				}
			}
			makeRsIDMostSignificant = !makeRsIDMostSignificant;
		} while (makeRsIDMostSignificant);
		
		/*
		int rsId = 67108863;
		int position = 2147483647;
		boolean strand = true;		
		System.out.println("Strand: " + strand);
		Nuc nuc1 = Nuc.A;
		Nuc nuc2 = Nuc.T;
				
		long result = compactSNPInfo(position, rsId, nuc1, nuc2, false);
		System.out.println("Result = " + result);
		System.out.println("BinString = " + Long.toBinaryString(result));
		System.out.println("BinString Length = " + Long.toBinaryString(result).length());
		System.out.println("Extracted position = " + extractPositionFromCompactForm(result, false));
		System.out.println("Extracted RS ID = " + extractRsIdFromCompactForm(result, false));
		
		Nuc[] nucs = new Nuc[NumAllelesInCallSet];
		extractNucleotideIDsFromCompactForm(result, nucs); 
		System.out.println("Extracted Nuc1, Nuc2: " + getStringForBytes(nucs));
		*/
	}
	
	public static void TestCompareAllelePairs() {
		Nuc[] pair1 = new Nuc[] { Nuc.C, Nuc.T };
		Nuc[] pair2 = new Nuc[] { Nuc.C, Nuc.T };
		boolean[] result = new boolean[pair2.length];
		
		boolean testResult = SNVMap.compareAllelePairs(pair1, pair2, result);
		System.out.println(testResult + "\tComp:" + result[0] + "\tFlip:" + result[1]);
	}
	
	/*
	public static void TestCompareAllelePairsAffyAndIllumina() {
		String affyInFilename = "GenomeWideSNP_6.na24.annot.csv.All2.autosomes.csv";
		String illuminaInFilename = "Illumina_Hap550.bpm.edited.autosomes.strandFixed.csv";
		String commonRsIdInFilename = "affy.gbm.common.illumina550K.fromAscnCompare.rsIdList.txt";
		byte[] pairAffy = new byte[2];
		byte[] pairIllumina = new byte[pairAffy.length];
		boolean[] result = new boolean[pairAffy.length];
		
		ArrayList<String> rsIdList = IOUtils.readAllLinesFromFile(commonRsIdInFilename);
		AffySNPMap asmAffy = CopyNumberParse.parseAffySNPFileHelper(affyInFilename);
		AffySNPMap asmIllumina = CopyNumberParse.parseAffySNPFileHelper(illuminaInFilename);
		
		for (Iterator<String> iter = rsIdList.iterator(); iter.hasNext(); ) {
			int rsId = GenotypeUtils.getNumberFromRsId(iter.next());
			byte[] resultAffy = asmAffy.getNucleotidesInMapByRsID(rsId, pairAffy);
			byte[] resultIllumina = asmIllumina.getNucleotidesInMapByRsID(rsId, pairIllumina);
			
			Utils.ensureTrue(resultAffy != null, "ERROR: RsID does not exist in Affy!: " + rsId);
			Utils.ensureTrue(resultIllumina != null, "ERROR: RsID does not exist in Illumina!: " + rsId);
			
			boolean compResult = compareAllelePairs(pairAffy, pairIllumina, result);
			if (!compResult) {
				String errString = "ERROR: SNPs do not match between platforms!\trs" + rsId 
					+ "Affy: {" + pairAffy[0] + "," + pairAffy[1]
					+ "}\tIllumina: {" + pairIllumina[0] + "," + pairIllumina[1] + "}";
				Utils.ensureTrue(false, errString);
			} else {
				String outString = "Match:\trs" + rsId
					+ "\tResult: {"+ result[0] + "," + result[1] + "}";
				System.out.println(outString);
			}
		}
		System.out.println("All matched!");
	}
	*/
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		//TestSNPCompact();
		TestSNVMap_Robust();
		//TestAffySNPMap();
		//TestAffySNPMapWhole();
		//TestCompareAllelePairs();
		//TestCompareAllelePairsAffyAndIllumina();
	}

}
