package genomeUtils;

import java.io.BufferedWriter;
import java.util.Arrays;
import java.util.Comparator;

import org.apache.commons.math3.random.RandomDataGenerator;

import nutils.CompareUtils;
import nutils.IOUtils;
import nutils.NumberUtils;
import nutils.StringUtils;
import genomeEnums.Chrom;
import genomeEnums.Genotype;
import genomeEnums.Nuc;
import shared.Utils;

public class GenotypeUtils {

	// ========================================================================
	// INNER COMPARATOR CLASS
	// ========================================================================
	/** We define a comparator for sorting rows.  This assumes that the columns are tab
	 *  delimited and that the chromosome resides in the first column (in the form of an 
	 *  integer or chrN, where N represents an integer), while an integer genomic coordinate 
	 *  resides in the second column. 
	 */
	public static class GenomicCoordinateComparatorInTextFileLine implements Comparator<String> {
		public int compare(String line1, String line2) {
			long line1Coordinates = readChromNumAndBasePairPosition(line1);
			long line2Coordinates = readChromNumAndBasePairPosition(line2);
			return Long.compare(line1Coordinates, line2Coordinates);			
		}
			
		/** This function reads the chromosome number and base pair position and packs
		 *  them into a long variable, with the chromosome number taking up the first 32
		 *  most significant bits and the base pair position taking up the last 32 bits.
		 *  This allows for returning multiple integer values and for easy sorting therafter
		 */
		private long readChromNumAndBasePairPosition(String line) {
			int firstTabIndexLine = line.indexOf(StringUtils.TabStr);			
			if (firstTabIndexLine < 0) CompareUtils.throwErrorAndExit("ERROR: Columns are not tab delimited in input line 1!");
	
			int secondTabIndexLine = line.indexOf(StringUtils.TabStr, firstTabIndexLine + 1);
			if (secondTabIndexLine < 0) CompareUtils.throwErrorAndExit("ERROR: Columns are not tab delimited in input line 1!");
	
			String chromString = "";
			int indexChromPrefix = line.indexOf(Chrom.ChromPrefix_chr);
			if (indexChromPrefix < 0) { // it does not exist
				chromString = line.substring(0, firstTabIndexLine);
			} else if (indexChromPrefix < firstTabIndexLine) {
				chromString = line.substring(indexChromPrefix + Chrom.ChromPrefix_chr.length(), firstTabIndexLine);				
			} else if (indexChromPrefix < secondTabIndexLine) {
				CompareUtils.throwErrorAndExit("ERROR: " + Chrom.ChromPrefix_chr + " string in an incorrect place on line!");
			}  // We don't care if the string exists elsewhere
			
			int chromNum = Chrom.getChrom(chromString).getCode();					
			int basePairPosition = Integer.parseInt(line.substring(firstTabIndexLine + 1, secondTabIndexLine));
			//System.out.println(chromNum + "\t" + basePairPosition);
			return ( 0L | (((long) chromNum) << 32) | ((long) basePairPosition) );   // pack two values into one long	
		}
	}

	public static final int RsID_Unknown = -1;
	public static final int RsID_Novel   = -5;
	public static final String RsPrefix = "rs";

	// ========================================================================
	public static double adjustHaploidCopyNumber(double copyNumberPure, double fractionPurity) {
		return adjustCopyNumberHelper(copyNumberPure, fractionPurity, GenomeConstants.DefaultHaploidCopyNumber);
	}
	
	// ========================================================================
	public static double adjustDiploidCopyNumber(double copyNumberPure, double fractionPurity) {
		return adjustCopyNumberHelper(copyNumberPure, fractionPurity, GenomeConstants.DefaultDiploidCopyNumber);
	}

	// ========================================================================
	protected static double adjustCopyNumberHelper(double copyNumberPure, double fractionPurity, double copyNumberStromal) {
		return ((fractionPurity * copyNumberPure) + ((1 - fractionPurity) * copyNumberStromal));
	}
	
	// ========================================================================
	public static void defineGenotypeAlleles(Genotype genotypeCode, Nuc referenceAllele, Nuc otherAllele, Nuc[] allelePhase) {
		switch(genotypeCode) {
		case EnumHomozygous00: 
			allelePhase[0] = allelePhase[1] = referenceAllele; 
			break;
		case EnumHeterozygous:
			// Doesn't matter the order we put them in since genotype isn't supposed to be phased
			allelePhase[0] = referenceAllele;
			allelePhase[1] = otherAllele;
			break;
		case EnumHomozygous11: 
			allelePhase[0] = allelePhase[1] = otherAllele; 
			break;
		case EnumInvalidGenotype: case EnumHomoDeletion:
			allelePhase[0] = allelePhase[1] = Nuc.N;
			break;
		default: 
			CompareUtils.ensureTrue(false, "ERROR: defineGenotypeAlleles(): Invalid Genotype!");
			break;
		}
	}
	
	// ========================================================================
	/** Given a pair of genotype alleles and a reference, this returns the variant allele. 
	 * @return the variant allele, or null if both match the ref, or N if neither match the ref.
	 */
	public static Nuc deduceVariant(Nuc genoPhase0, Nuc genoPhase1, Nuc referenceAllele) {
		if ((genoPhase0 == referenceAllele) && (genoPhase1 == referenceAllele)) {
			return null;
		} else if (genoPhase0 == referenceAllele) {
			return genoPhase1;
		} else if (genoPhase1 == referenceAllele) {
			return genoPhase0;
		} else {
			return Nuc.N;
		}
	}
	
	// =======================================================================
	public static Chrom getRandomAutosomalChromosome() { 
		int index = NumberUtils.getRandomInteger(Chrom.IndexAutosomalStart, Chrom.IndexAutosomalEnd);
		return Chrom.getChrom((byte) index);
	}
	
	// ========================================================================
	/** Given an rsID string, this returns the integer from the rs ID. */
	public static int getNumberFromRsId(String rsIdString) {
		rsIdString = rsIdString.trim();
		int index = rsIdString.indexOf(GenotypeUtils.RsPrefix);
		if (index >= 0) {
			return Integer.parseInt(rsIdString.substring(index + GenotypeUtils.RsPrefix.length())); 
		} else {
			return RsID_Unknown;
		}
	}

	/** Given an rs number, this returns the rs# string form. */
	public static String getRsIdFromNumber(int rsNum) { return "rs" + rsNum; }
	
	/** @return -1, if no valid nucleotides are present in the string; or if valid nucleotides are present,
	 *  the fraction of G/C nucleotides
	 */
	public static synchronized double calcFractionGC(String nucleotideString) {		
		int strLen = nucleotideString.length();
		int validLength = 0;
		int numGC = 0;
		nucleotideString.getChars(0, strLen, charBuffer, 0);  // more efficient than calling charAt() each character
		
		for (int i = 0; i < strLen; i++) {
			Nuc theNuc = Nuc.getNucUnsafe(charBuffer[i]);
			if (theNuc != null) {
				validLength++;
				if (theNuc.isGC()) {
					numGC++;
				}
			}
		}
		
		// Now divide
		if (validLength == 0) return -1;
		return (double) numGC / (double) validLength;
	}

	/** Given a string of nucleotides, this returns the fraction of characters that are G or C. */	
	private static final char[] charBuffer = new char[65536];

	/** Given a string that contains an rs#, this returns an rs#.  
	 *  If the string is "rs2342,abcd", this will return: "rs2342"
	 *  @return the number after rs and before a non-numeric character, if all valid
	 *  @return "rs-1", if the rs# is -1 (rs-1), as listed in some files
	 *  @return "rs", if there is no valid rs# after rs
	 *  @return null, if there is no "rs" substring in this line.  
	 */
	public static String extractRsNumberFromLine(String line) {
		int rsIndex = line.indexOf(Utils.rsPrefix);
		if (rsIndex >= 0) {
			int startIndex = rsIndex + Utils.rsPrefix.length(); 
			for (int i = startIndex; i < line.length(); i++) {
				char ch = line.charAt(i);
				if (!Character.isDigit(ch)) {
					// Two things can happen.  We are either at the end of the rs#, 
					// or we have rs-1, which we must be able to parse					
					if ((ch == '-') 
						&& (i == startIndex) 
						&& (i < line.length() - 1) 
						&& (line.charAt(i + 1) == '1')) {						
						return Utils.rsNegative;					
					} else {
						return line.substring(rsIndex, i);
					}
				}
			}
			return line.substring(rsIndex);
		} else {
			return null;
		}
	}

	// ========================================================================
	/** This method count the number of each nucleotide in the given string.
	 *  @param seq The string to parse
	 *  @param nc  The counter that will store the counts 
	 *  @return Will return the invalid nucleotide in the string if there is one, else returns 'A'. */
	/*
	public static char countConstituentAlleles(CharSequence seq, NucCounter nc) {
		nc.clear();
	    for (int i = 0; i < seq.length(); i++) {
	    	Nuc theNuc = Nuc.getNucUnsafe(seq.charAt(i));
	    	if (theNuc == null) { return seq.charAt(i); }
	    	nc.incrementCount(theNuc);
	    }
	    return 'A';
	}*/
	
	public static void Test_extractRsNumberFromLine() {
		String[] toTest = new String[] {
				"rs-1",
				",rs-1",
				"rs-1,",
				";rs-1;",
	
				"rs-2",
				",rs-0",
				"rs-0,",
				";rs-0;",
				";rs;-1",
				";rs'123",
	
				"rs12345",
				"rs12345a",
				",rs12345",
				",rs12345,",
				
				",rs12$345,"
		};
		
		for (int i = 0; i < toTest.length; i++) {
			System.out.println(extractRsNumberFromLine(toTest[i]));
		}
	}

	// ------------------------------------------------------------------------
	public static double getPValuesImbalanceTissue(int coverageTotal, int coverageVariant) {
		int maxRefOrVarCovg = Math.max(coverageVariant, coverageTotal - coverageVariant);
		return nutils.NumberUtils.cumulativeProbabilitySuccessBinomial(coverageTotal, maxRefOrVarCovg, 0.5);
	}

	// ========================================================================
	public static double calcSigPValueByPermutation(int readCountTotal, double significanceLevel, double[] pValueBuffer) {		
		RandomDataGenerator randomGen = new RandomDataGenerator();
		
		int numIter = pValueBuffer.length;		
		for (int iter = 0; iter < numIter; iter++) {			
			int numReadsA = randomGen.nextBinomial(readCountTotal, 0.5); //(int) randomGen.nextPoisson(readCountHaplotype);
			numReadsA = Math.min(numReadsA, readCountTotal);								
			pValueBuffer[iter] = getPValuesImbalanceTissue(readCountTotal, numReadsA);
		}
		
		Arrays.sort(pValueBuffer);
		int indexOfSignificance = (int) (pValueBuffer.length * significanceLevel);
		return pValueBuffer[indexOfSignificance];
	}
	
	// ========================================================================
	public static double[] PrecalculatePerSitePermutationThresholds() {
		int numMaxTrials = 10000;
		double[] rV = new double[numMaxTrials + 1];
		
		int numIter = 10000;
		double[] buffer = new double[numIter];

		for (int trialNum = 1; trialNum <= numMaxTrials; trialNum++) {
			rV[trialNum] = calcSigPValueByPermutation(trialNum, 0.05, buffer);
			System.out.println(trialNum + "\t" + rV[trialNum]);
		}
		
		return rV;
	}

	// ========================================================================
	public static void PrintPrecalculatedPermutationSiteThresholds() {
		double[] sigPVals = PrecalculatePerSitePermutationThresholds();
		BufferedWriter out = IOUtils.getBufferedWriter("CalcPermValues.txt");
		for (int i = 0; i < sigPVals.length; i++) {		
			IOUtils.writeToBufferedWriter(out, i + "\t" + sigPVals[i], true);			
		}
		IOUtils.closeBufferedWriter(out);
	}
	
	// ========================================================================
	public static void main(String[] args) {
		PrintPrecalculatedPermutationSiteThresholds();
	}
}
