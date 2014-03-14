package lohcate;
import genomeEnums.Chrom;

import java.io.BufferedReader;

import nutils.CompareUtils;
import nutils.EnumSortedMapLong;
import nutils.IOUtils;
import nutils.PrimitiveWrapper;
import nutils.StringUtils;
import nutils.BitUtils.Compactor.CompactorInf;
import nutils.BitUtils.Compactor.CompactorIntoLong;

/** This table will contain values of bias in an optimized fashion. */

public class AllelicBiasTableCompact implements AllelicBiasTable {
	
	protected final static float VAFNormalMultiplier = 10000;
	
	// An arraylist for each chromosome
	EnumSortedMapLong<Chrom> mPositionsAndVAF;		
	PrimitiveWrapper.WLong mLongHolder;

	// ========================================================================
	public AllelicBiasTableCompact() {		
		mPositionsAndVAF = new EnumSortedMapLong<>(Chrom.class, PositionAndVAF.Compactor.getValueExtractor(PositionAndVAF.Position));
		mLongHolder = new PrimitiveWrapper.WLong(0);
	}

	// ========================================================================
	/* (non-Javadoc)
	 * @see lohcate.AllelicBiasTable#getAvgVAF(genomeEnums.Chrom, int)
	 */
	@Override
	public float getAvgVAF(Chrom chrom, int position) {
		return getAvgVAF(chrom, position, 0);
	}

	// ========================================================================
	private static float storedVAFtoVAF(long unit) {
		return PositionAndVAF.Compactor.getValue(PositionAndVAF.VAFNormal, unit) / VAFNormalMultiplier;
	}
	
	// ========================================================================
	private static long vafToStoredVAF(float vaf) {
		return Math.round(vaf * VAFNormalMultiplier);		
	}

	// ========================================================================
	/* (non-Javadoc)
	 * @see lohcate.AllelicBiasTable#getNumSitesRegisteredAtPosition(genomeEnums.Chrom, int)
	 */
	@Override
	public long getNumSitesRegisteredAtPosition(Chrom chrom, int position) {
		long key = PositionAndVAF.Compactor.setValue(PositionAndVAF.Position, position, 0L);
		boolean exists = mPositionsAndVAF.get(chrom, key, mLongHolder);
		return exists ? 
				PositionAndVAF.Compactor.getValue(PositionAndVAF.NumSamples, mLongHolder.mLong) :
				-1;
	}
	
	// ========================================================================
	/* (non-Javadoc)
	 * @see lohcate.AllelicBiasTable#getAvgVAF(genomeEnums.Chrom, int, int)
	 */
	@Override
	public float getAvgVAF(Chrom chrom, int position, int siteThresholdInclusive) {
		long key = PositionAndVAF.Compactor.setValue(PositionAndVAF.Position, position, 0L);
		boolean exists = mPositionsAndVAF.get(chrom, key, mLongHolder);
		if (exists) {
			long numSitesRepresented = PositionAndVAF.Compactor.getValue(PositionAndVAF.NumSamples, mLongHolder.mLong);
			if (numSitesRepresented >= siteThresholdInclusive) {
				return storedVAFtoVAF(mLongHolder.mLong);
			} else {
				return -2;
			}
		} else {
			return -1;
		}
	}

	// ========================================================================
	private static long makePositionVAFCompactUnit(int position, int numSamplesRepresented, float vaf) {
		long key = PositionAndVAF.Compactor.setValue(PositionAndVAF.Position,   position,                   0L);
		key      = PositionAndVAF.Compactor.setValue(PositionAndVAF.NumSamples, numSamplesRepresented,     key);
		key      = PositionAndVAF.Compactor.setValue(PositionAndVAF.VAFNormal,  vafToStoredVAF(vaf),       key);
		return key;
	}
	
	// ========================================================================
	/* (non-Javadoc)
	 * @see lohcate.AllelicBiasTable#registerSite(genomeEnums.Chrom, int, float)
	 */
	@Override
	public void registerSite(Chrom chrom, int position, float vafNormal) {
		// Make the key
		long key = makePositionVAFCompactUnit(position, 1, vafNormal);

		int resultIndex = mPositionsAndVAF.addSorted(chrom, key, mLongHolder);
		if (resultIndex >= 0) {
			final long newUnit = registerVAFNormalForAnotherSite(vafNormal, mLongHolder.mLong);
			mPositionsAndVAF.replace(chrom, resultIndex, newUnit);
		} 
	}

	// ========================================================================
	/** Updates the unit to incorporate the new vaf value and returns the new unit */
	private long registerVAFNormalForAnotherSite(final double vafNormaNew, final long compactUnit) {
		CompareUtils.ensureTrue(vafNormaNew >= 0, "ERROR: vafNormal cannot be less than 0!");
		CompareUtils.ensureTrue(vafNormaNew <= 1, "ERROR: vafNormal cannot be greater than 1!");
		
		float vafNormal = storedVAFtoVAF(compactUnit);
		long numSamplesRepresented = PositionAndVAF.Compactor.getValue(PositionAndVAF.NumSamples, compactUnit);
		
		vafNormal *= numSamplesRepresented;
		vafNormal += vafNormaNew;		
		vafNormal /= (++numSamplesRepresented);
		
		long compactUnitNew = compactUnit;  // Make a copy 
		compactUnitNew = PositionAndVAF.Compactor.setValue(PositionAndVAF.NumSamples, numSamplesRepresented, compactUnitNew);
		compactUnitNew = PositionAndVAF.Compactor.setValue(PositionAndVAF.VAFNormal, vafToStoredVAF(vafNormal), compactUnitNew);
		return compactUnitNew;
	}
	
	// ========================================================================
	/** Given a filename that represents the vaf allelic averages, this stores the 
	 *  contents of the file. 
	 */	
	public static AllelicBiasTable readFileAndConstructTable(String inFilename, int colNumSamples, int colAvgVAFNormal) {
		BufferedReader in = IOUtils.getBufferedReader(inFilename);
		String line;
		String delim = StringUtils.FileExtensionTSV.mDelimiter;
		
		AllelicBiasTableCompact allelicBiasTable = new AllelicBiasTableCompact();
		while ((line = IOUtils.getNextLineInBufferedReader(in)) != null) {
			Chrom chrom    = Chrom.getChrom(StringUtils.extractNthColumnValue(line, 0, delim));
			int position = Integer.parseInt(StringUtils.extractNthColumnValue(line, 1, delim));
			int numSamplesRepresented = Integer.parseInt(StringUtils.extractNthColumnValue(line, colNumSamples,   delim));
			float avgVAFNormal        = Float.parseFloat(StringUtils.extractNthColumnValue(line, colAvgVAFNormal, delim));
			
			long compactUnit = makePositionVAFCompactUnit(position, numSamplesRepresented, avgVAFNormal);
			allelicBiasTable.mPositionsAndVAF.addToTail(chrom, compactUnit);			
		}		
		IOUtils.closeBufferedReader(in);

		allelicBiasTable.mPositionsAndVAF.sortTable();
		return allelicBiasTable;
	}

	// ========================================================================
	private static enum PositionAndVAF implements CompactorInf<PositionAndVAF> {
		Position(29),
		NumSamples(16),
		VAFNormal(16)
		;

		// ========================================================================
		public static final CompactorIntoLong<PositionAndVAF> Compactor = new CompactorIntoLong<>(PositionAndVAF.class, false);
		
		// ========================================================================		
		int mNumBits;
		private PositionAndVAF(int numBits) {
			mNumBits = numBits;
		}
		// ========================================================================
		@Override
		public int getNumBits() { return mNumBits; }		
	}
	
	// ========================================================================
	private static void Test() {
		AllelicBiasTable abt = new AllelicBiasTableCompact();
		AllelicBiasTableFull abtOld = new AllelicBiasTableFull();
		
		int numIter = 100;
		for (Chrom chrom : Chrom.Autosomes) {
			for (int i = 0; i < numIter; i++) {
				float vafNormal = (float) i / (float) numIter;
				for (int j = 0; j < i; j++) {					
					abt.registerSite(chrom, i, vafNormal);
					abtOld.registerSite(chrom, i, vafNormal);
				}
			}
		}
		
		for (Chrom chrom : Chrom.Autosomes) {
			for (int i = 0; i < numIter; i++) {
				float avgVAF = abt.getAvgVAF(chrom, i);
				float avgVAFOld = abtOld.getAvgVAF(chrom, i);
				
				if (Math.abs(avgVAF - avgVAFOld) > 0.00001) {
					System.out.printf("ERROR: Mismatch: %s\t%d\t%f\t%f\n", chrom, i, avgVAF, avgVAFOld);
				}
				//System.out.printf("%s\t%d\t%d\t%g\n", chrom, i, abt.getNumSitesRegisteredAtPosition(chrom, i), avgVAF);
			}
		}
		
		System.out.println("Success!");
		//abt.mPositionsAndVAF.print(System.out);
	}
	
	// ========================================================================
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		Test();
	}
	
}
