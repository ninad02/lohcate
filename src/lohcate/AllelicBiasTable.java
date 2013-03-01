package lohcate;
import genomeEnums.Chrom;

import java.io.BufferedReader;
import java.util.ArrayList;
import java.util.Collections;

import nutils.ArrayUtils;
import nutils.CompareUtils;
import nutils.EnumMapSafe;
import nutils.EnumSortedMap;
import nutils.IOUtils;
import nutils.StringUtils;



import com.carrotsearch.hppc.LongArrayList;


/** This table will contain values of bias in an optimized fashion. */

public class AllelicBiasTable {
	
	// An arraylist for each chromosome
	EnumSortedMap<Chrom, PositionAndPayload> mPositionsAndVAFs;
	PositionAndPayload mDummyPayload;

	// ========================================================================
	public AllelicBiasTable() {
		mDummyPayload = new PositionAndPayload(0, 0, 0);
		mPositionsAndVAFs = new EnumSortedMap<Chrom, PositionAndPayload>(Chrom.class); 
	}

	// ========================================================================
	/** Given a chromosome and position, this returns the average VAF reported for that position. 
	 *  @return the average vaf, or -1 if it doesn't exist. */
	public float getAvgVAF(Chrom chrom, int position) {
		return getAvgVAF(chrom, position, 0);
	}
	
	// ========================================================================
	/** Given a chromosome, position, and site threshold, this returns the average VAF reported for that position
	 *  so long as the number of sites that contribute to the average VAF meet or exceed the threshold (inclusive).
	 *  @return the average VAF, -1 if the position does not exist, or -2 if not enough sites meet the threshold.
	 */
	public float getAvgVAF(Chrom chrom, int position, int siteThresholdInclusive) {
		mDummyPayload.mPosition = position;
		final PositionAndPayload pap = mPositionsAndVAFs.get(chrom, mDummyPayload);
		return ( (pap == null) ? -1 : (pap.mNumSamplesRepresented >= siteThresholdInclusive ? pap.mVAFNormal : -2) );
	}

	// ========================================================================
	/** Registers a site with a normal vaf value. If a site already exists, the vaf value is weighted-averaged with the existing vaf value. */
	public void registerSite(Chrom chrom, int position, float vafNormal) {
		mDummyPayload.set(position, 1, vafNormal);
		PositionAndPayload papExists = mPositionsAndVAFs.addSorted(chrom, mDummyPayload);
		if (papExists == null) {
			mDummyPayload = new PositionAndPayload(0, 0, 0);  // Create new dummy object since previous was put into map			
		} else {
			papExists.registerVAFNormalForAnotherSite(vafNormal);
		}
	}
	
	// ========================================================================
	/** Given a filename that represents the vaf allelic averages, this stores the 
	 *  contents of the file. 
	 */	
	public static AllelicBiasTable readFileAndConstructTable(String inFilename, int colNumSamples, int colAvgVAFNormal) {
		BufferedReader in = IOUtils.getBufferedReader(inFilename);
		String line;
		String delim = StringUtils.FileExtensionTSV.mDelimiter;
		
		AllelicBiasTable allelicBiasTable = new AllelicBiasTable();
		while ((line = IOUtils.getNextLineInBufferedReader(in)) != null) {
			Chrom chrom    = Chrom.getChrom(StringUtils.extractNthColumnValue(line, 0, delim));
			int position = Integer.parseInt(StringUtils.extractNthColumnValue(line, 1, delim));
			int numSamplesRepresented = Integer.parseInt(StringUtils.extractNthColumnValue(line, colNumSamples,   delim));
			float avgVAFNormal        = Float.parseFloat(StringUtils.extractNthColumnValue(line, colAvgVAFNormal, delim));
			PositionAndPayload posAndPay = new PositionAndPayload(position, numSamplesRepresented, avgVAFNormal);
			
			allelicBiasTable.mPositionsAndVAFs.addToTail(chrom, posAndPay);			
		}		
		IOUtils.closeBufferedReader(in);

		allelicBiasTable.mPositionsAndVAFs.sortTable();
		return allelicBiasTable;
	}
	
	
	// ========================================================================
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub

	}

	// ========================================================================
	private static class PositionAndPayload implements Comparable<PositionAndPayload> {
		int mPosition;
		int mNumSamplesRepresented;
		float mVAFNormal;
		
		public PositionAndPayload(int position, int numSamplesRepresented, float avgVafNormal) {
			set(position, numSamplesRepresented, avgVafNormal);
		}
		
		public void set(int position, int numSamplesRepresented, float avgVafNormal) {
			mPosition = position;
			mNumSamplesRepresented = numSamplesRepresented;
			mVAFNormal = avgVafNormal;			
		}
		
		public int compareTo(PositionAndPayload rhs) {
			return Integer.compare(mPosition, rhs.mPosition);
		}
		
		public float registerVAFNormalForAnotherSite(double vafNormal) {
			CompareUtils.ensureTrue(vafNormal >= 0, "ERROR: vafNormal cannot be less than 0!");
			CompareUtils.ensureTrue(vafNormal <= 1, "ERROR: vafNormal cannot be greater than 1!");
						
			mVAFNormal *= mNumSamplesRepresented;
			mVAFNormal += vafNormal;
			mVAFNormal /= (++mNumSamplesRepresented);
			return mVAFNormal;
		}		
	}
	
	// ========================================================================
	/** Given a base pair position and a vafNormal value, and the number of samples
	 *  involved in the averaging of the vafNormal value, this packs them into a
	 *  long variable.  The vafNormal is rounded to the nearing 0.0001 and coverted
	 *  to a format that can be stored within a long.
	 */
	/* DON'T IMPLEMENT EFFICIENT VERSION FOR NOW
	private static final int ShiftPosition = 35;
	private static final long MaxValuePosition = Utils.getMask(28);
	private static final int MaskPosition = 0x0FFFFFFF;
	
	private static final int 
	public static long packPositionAndVAFAndNumSamples(int position, float avgVafNormal, int numSamplesRepresented) {
		position = Math.min(position, MaxValuePosition);
	}
	*/
	
}
