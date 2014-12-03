package lohcate;

import genomeEnums.Chrom;
import genomeUtils.GenomeConstants;

public class AllelicBiasTableDummy implements AllelicBiasTable {

	private static long NumSitesRegisteredAtPosition = 2;
	private static float DefaultAverageVAF = 0.5f; 
	private static float DefaultAverageAbsoluteCopyNumber = GenomeConstants.DefaultDiploidCopyNumber;
	
	// ========================================================================
	@Override
	public float getAvgVAF(Chrom chrom, int position) {
		return getAvgVAF(chrom, position, 0);
	}

	// ========================================================================
	@Override
	public long getNumSitesRegisteredAtPosition(Chrom chrom, int position) {
		return NumSitesRegisteredAtPosition;
	}

	// ========================================================================
	@Override
	public float getAvgVAF(Chrom chrom, int position, int siteThresholdInclusive) {
		return DefaultAverageVAF;
	}

	// ========================================================================
	@Override
	public void registerSite(Chrom chrom, int position, float vafNormal, float absoluteCopyNumber) {
		// Do nothing
	}

	// ========================================================================
	@Override
	public float getAvgAbsCopyNumber(Chrom chrom, int position) { 		
		return getAvgAbsCopyNumber(chrom, position, 0);
	}

	// ========================================================================
	@Override
	public float getAvgAbsCopyNumber(Chrom chrom, int position, int siteThresholdInclusive) {
		return DefaultAverageAbsoluteCopyNumber;
	}

}
