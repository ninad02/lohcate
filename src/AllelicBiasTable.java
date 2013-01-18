import genomeEnums.Chrom;

import java.io.BufferedReader;
import java.util.ArrayList;
import java.util.Collections;

import nutils.ArrayUtils;
import nutils.IOUtils;
import nutils.StringUtils;



import com.carrotsearch.hppc.LongArrayList;


/** This table will contain values of bias in an optimized fashion. */

public class AllelicBiasTable {
	
	// An arraylist for each chromosome
	ArrayList< ArrayList<PositionAndPayload> > mPositionsAndVAFs;
	PositionAndPayload mDummyPayload;
	
	public AllelicBiasTable() {
		mDummyPayload = new PositionAndPayload(0, 0, 0);
		mPositionsAndVAFs = ArrayUtils.addNewEmptyArrayLists(null, Chrom.values().length);
	}

	// ========================================================================
	/** Given a chromosome and position, this returns the average VAF reported for that position. */
	public float getAvgVAF(Chrom chrom, int position) {
		ArrayList<PositionAndPayload> listForChrom = mPositionsAndVAFs.get(chrom.ordinal());
		mDummyPayload.mPosition = position;
		int resultIndex = Collections.binarySearch(listForChrom, mDummyPayload);
		if (resultIndex >= 0) {
			PositionAndPayload pap = listForChrom.get(resultIndex);
			return pap.mVAFNormal;
		}
		
		return -1;
	}
	
	
	// ========================================================================
	/** Sorts the table by position per chromosome. */
	private void sortTable() {
		for (ArrayList<PositionAndPayload> listPerChrom : mPositionsAndVAFs) {
			Collections.sort(listPerChrom);
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
			
			allelicBiasTable.mPositionsAndVAFs.get(chrom.ordinal()).add(posAndPay);
		}		
		IOUtils.closeBufferedReader(in);
		
		allelicBiasTable.sortTable();
		return allelicBiasTable;
	}
	
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub

	}

	private static class PositionAndPayload implements Comparable<PositionAndPayload> {
		int mPosition;
		int mNumSamplesRepresented;
		float mVAFNormal;
		
		public PositionAndPayload(int position, int numSamplesRepresented, float avgVafNormal) {
			mPosition = position;
			mNumSamplesRepresented = numSamplesRepresented;
			mVAFNormal = avgVafNormal;
		}
		
		public int compareTo(PositionAndPayload rhs) {
			return Integer.compare(mPosition, rhs.mPosition);
		}
		
	}
	
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
