package lohcate;
import genomeEnums.Chrom;

import java.io.BufferedReader;

import nutils.Cast;
import nutils.CompareUtils;
import nutils.EnumSortedMapLong;
import nutils.EnumSortedMapLongFast;
import nutils.IOUtils;
import nutils.StringUtils;
import nutils.BitUtils.Compactor.CompactorInf;
import nutils.BitUtils.Compactor.CompactorIntoLong;
import nutils.primitives.wrapper.PrimitiveWrapper;

/** This table will contain values of bias in an optimized fashion. */

public class AllelicBiasTableCompact implements AllelicBiasTable {
	
	protected final static float VAFNormalMultiplier = 10000;
	
	// An arraylist for each chromosome
	EnumSortedMapLong<Chrom> mPositionsAndVAF;		
	PrimitiveWrapper.WLong mLongHolder;
	EnumSortedMapLong<Chrom> mTempMap;

	// ========================================================================
	public AllelicBiasTableCompact() {		
		mPositionsAndVAF = new EnumSortedMapLong<>(Chrom.class, PositionAndVAF.Compactor.getValueExtractor(PositionAndVAF.Position));
		mTempMap         = new EnumSortedMapLong<>(Chrom.class, PositionAndVAF.Compactor.getValueExtractor(PositionAndVAF.Position));
		mLongHolder = new PrimitiveWrapper.WLong(0);		
	}

	// ========================================================================
	@Override
	public float getAvgAbsCopyNumber(Chrom chrom, int position) {
		return getAvgAbsCopyNumber(chrom, position, 0);
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
				AllelicBiasTable.ChromPositionDoesNotExist;
	}
	
	// ========================================================================
	/* (non-Javadoc)
	 * @see lohcate.AllelicBiasTable#getAvgVAF(genomeEnums.Chrom, int, int)
	 */
	@Override
	public float getAvgVAF(Chrom chrom, int position, int siteThresholdInclusive) {		
		return getAvg_VAF_CopyNum_Helper(chrom, position, siteThresholdInclusive, PositionAndVAF.VAFNormal);
	}
	
	// ========================================================================
	@Override
	public float getAvgAbsCopyNumber(Chrom chrom, int position, int siteThresholdInclusive) {
		return getAvg_VAF_CopyNum_Helper(chrom, position, siteThresholdInclusive, PositionAndVAF.AbsoluteCopyNumber);
	}

	// ========================================================================
	private float getAvg_VAF_CopyNum_Helper(Chrom chrom, int position, int siteThresholdInclusive, PositionAndVAF target) {
		long key = PositionAndVAF.Compactor.setValue(PositionAndVAF.Position, position, 0L);
		boolean exists = mPositionsAndVAF.get(chrom, key, mLongHolder);
		if (exists) {
			long numSitesRepresented = PositionAndVAF.Compactor.getValue(PositionAndVAF.NumSamples, mLongHolder.mLong);
			if (numSitesRepresented >= siteThresholdInclusive) {
				switch(target) {
				case VAFNormal:          return storedVAFtoVAF(mLongHolder.mLong);
				case AbsoluteCopyNumber: return PositionAndVAF.Compactor.getValue(PositionAndVAF.AbsoluteCopyNumber, mLongHolder.mLong);
				default: CompareUtils.ensureTrue(false, "Invalid target!");
				}				
				return 0;
			} else {
				return AllelicBiasTable.InsufficientSites;
			}
		} else {
			return AllelicBiasTable.ChromPositionDoesNotExist;
		}
	}
		
	// ========================================================================
	private static long makePositionVAFCompactUnit(int position, int numSamplesRepresented, float vaf, float absoluteCopyNumber) {
		long key = PositionAndVAF.Compactor.setValue(PositionAndVAF.Position,   position,                   0L);
		key      = PositionAndVAF.Compactor.setValue(PositionAndVAF.NumSamples, numSamplesRepresented,      key);
		key      = PositionAndVAF.Compactor.setValue(PositionAndVAF.VAFNormal,  vafToStoredVAF(vaf),        key);
		key      = PositionAndVAF.Compactor.setValue(PositionAndVAF.AbsoluteCopyNumber, Math.round(absoluteCopyNumber), key);
		return key;
	}

	// ========================================================================
	public int size() { return mPositionsAndVAF.size(); }
	
	// ========================================================================
	public void flushSitesOptimized() {
		
		// Each site in the temp map is guaranteed to be unique.  We just 
		// add these to the tail of the real map, which will sort within itself.		
		for (Chrom theEnum : Chrom.class.getEnumConstants()) {		
			mPositionsAndVAF.ensureCapacity(theEnum, mPositionsAndVAF.size(theEnum) + mTempMap.size(theEnum));
			for (int i = 0; i < mTempMap.size(theEnum); i++) {
				mPositionsAndVAF.addToTail(theEnum, mTempMap.get(theEnum, i));
			}			
			System.out.printf("\t\tChrom %s:\t%d\n", theEnum, mPositionsAndVAF.size(theEnum));
		}
		mTempMap.clear();
	}
	
	// ========================================================================
	public void registerSiteOptimized(Chrom chrom, int position, float vafNormal, float absoluteCopyNumber) {
		
		long key = makePositionVAFCompactUnit(position, 1, vafNormal, absoluteCopyNumber);
		
		int resultIndex = mPositionsAndVAF.getIndex(chrom, key);
		if (resultIndex >= 0) {
			// The key exists in the table.  We simply update			
			final long newUnit = registerVAFNormalForAnotherSite(vafNormal, absoluteCopyNumber, mPositionsAndVAF.get(chrom, resultIndex));
			boolean result = mPositionsAndVAF.replace(chrom, resultIndex, newUnit);
			CompareUtils.ensureTrue(result, "ERROR: Replacing a site must always succeed");
		} else {
			// Add sites in a unique manner to a spillover buffer 
			mTempMap.addSorted(chrom, key, mLongHolder);			
		}
	}
	
	
	// ========================================================================
	/* (non-Javadoc)
	 * @see lohcate.AllelicBiasTable#registerSite(genomeEnums.Chrom, int, float)
	 */
	@Override
	public void registerSite(Chrom chrom, int position, float vafNormal, float absoluteCopyNumber) {
		// Make the key
		long key = makePositionVAFCompactUnit(position, 1, vafNormal, absoluteCopyNumber);

		/*
		boolean result = mPositionsAndVAF.addIfUnique(chrom, key, mLongHolder);
		if (!result) {
			final long newUnit = registerVAFNormalForAnotherSite(vafNormal, absoluteCopyNumber, mLongHolder.mLong);
			mPositionsAndVAF.replace(chrom, newUnit);
		}*/
		
		int resultIndex = mPositionsAndVAF.getIndex(chrom, key);
		if (resultIndex >= 0) {
			// The key exists in the table.  We simply update			
			final long newUnit = registerVAFNormalForAnotherSite(vafNormal, absoluteCopyNumber, mPositionsAndVAF.get(chrom, resultIndex));
			mPositionsAndVAF.replace(chrom, resultIndex, newUnit);
		}
		boolean result = mPositionsAndVAF.get(chrom, key, mLongHolder);
		if (result) {
			
		}
		
/*		int resultIndex = mPositionsAndVAF.addSorted(chrom, key, mLongHolder);
		if (resultIndex >= 0) {
			final long newUnit = registerVAFNormalForAnotherSite(vafNormal, absoluteCopyNumber, mLongHolder.mLong);
			mPositionsAndVAF.replace(chrom, resultIndex, newUnit);
		} */
	}

	// ========================================================================
	/** Updates the unit to incorporate the new vaf value and returns the new unit */
	private static long registerVAFNormalForAnotherSite(final double vafNormaNew, final double absCopyNumberNew, final long compactUnit) {
		CompareUtils.ensureTrue(vafNormaNew >= 0, "ERROR: vafNormal cannot be less than 0!");
		CompareUtils.ensureTrue(vafNormaNew <= 1, "ERROR: vafNormal cannot be greater than 1!");		
		CompareUtils.ensureTrue(absCopyNumberNew >= 0, "ERROR: Absolute Copy Number cannot be less than 0!");
		
		float vafNormal = storedVAFtoVAF(compactUnit);
		long numSamplesRepresented = PositionAndVAF.Compactor.getValue(PositionAndVAF.NumSamples,         compactUnit);
		float avgCopyNum           = PositionAndVAF.Compactor.getValue(PositionAndVAF.AbsoluteCopyNumber, compactUnit);
		
		vafNormal  = Cast.toFloat((vafNormal  * numSamplesRepresented) + vafNormaNew);
		avgCopyNum = Cast.toFloat((avgCopyNum * numSamplesRepresented) + absCopyNumberNew);	
		
		// Increment the numSamplesRepresented, and update any other values
		++numSamplesRepresented;
		vafNormal  /= numSamplesRepresented;
		avgCopyNum /= numSamplesRepresented;

		// Now repack
		long compactUnitNew = compactUnit;  // Make a copy 
		compactUnitNew = PositionAndVAF.Compactor.setValue(PositionAndVAF.NumSamples, numSamplesRepresented, compactUnitNew);
		compactUnitNew = PositionAndVAF.Compactor.setValue(PositionAndVAF.VAFNormal, vafToStoredVAF(vafNormal), compactUnitNew);
		compactUnitNew = PositionAndVAF.Compactor.setValue(PositionAndVAF.AbsoluteCopyNumber, Math.round(avgCopyNum), compactUnitNew);
		return compactUnitNew;
	}
	
	// ========================================================================
	/** Given a filename that represents the vaf allelic averages, this stores the 
	 *  contents of the file. 
	 */	
	public static AllelicBiasTable readFileAndConstructTable(String inFilename, int colNumSamples, int colAvgVAFNormal, int colAvgCopyNumber) {
		BufferedReader in = IOUtils.getBufferedReader(inFilename);
		String line;
		String delim = StringUtils.FileExtensionTSV.mDelimiter;
		
		AllelicBiasTableCompact allelicBiasTable = new AllelicBiasTableCompact();
		while ((line = IOUtils.getNextLineInBufferedReader(in)) != null) {
			Chrom chrom    = Chrom.getChrom(StringUtils.extractNthColumnValue(line, 0, delim));
			int position = Integer.parseInt(StringUtils.extractNthColumnValue(line, 1, delim));
			int numSamplesRepresented = Integer.parseInt(StringUtils.extractNthColumnValue(line, colNumSamples,   delim));
			float avgVAFNormal        = Float.parseFloat(StringUtils.extractNthColumnValue(line, colAvgVAFNormal, delim));
			float avgAbsCopyNumber    = Float.parseFloat(StringUtils.extractNthColumnValue(line, colAvgCopyNumber, delim));
			
			long compactUnit = makePositionVAFCompactUnit(position, numSamplesRepresented, avgVAFNormal, avgAbsCopyNumber);
			allelicBiasTable.mPositionsAndVAF.addToTail(chrom, compactUnit);			
		}		
		IOUtils.closeBufferedReader(in);

		allelicBiasTable.mPositionsAndVAF.sortTable();
		return allelicBiasTable;
	}

	// ========================================================================
	private static enum PositionAndVAF implements CompactorInf<PositionAndVAF> {
		Position(29),
		NumSamples(15),
		VAFNormal(14),
		AbsoluteCopyNumber(4)
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
	// ========================================================================
	// ========================================================================
	// ========================================================================
	// ========================================================================
	private static void Test() {
		AllelicBiasTable abt = new AllelicBiasTableCompact();
		AllelicBiasTableFull abtOld = new AllelicBiasTableFull();
		
		int numIter = 100;
		for (Chrom chrom : Chrom.Autosomes) {
			for (int i = 0; i < numIter; i++) {
				float vafNormal = (float) i / (float) numIter;
				for (int j = 0; j < i; j++) {					
					abt.registerSite(chrom, i, vafNormal, 3);
					abtOld.registerSite(chrom, i, vafNormal);
				}
			}
		}
		
		for (Chrom chrom : Chrom.Autosomes) {
			for (int i = 0; i < numIter; i++) {
				float avgVAF = abt.getAvgVAF(chrom, i);
				float avgVAFOld = abtOld.getAvgVAF(chrom, i);
				
				System.out.printf("%g\t%g\n", avgVAF, avgVAFOld);
				
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
