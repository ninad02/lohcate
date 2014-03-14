package lohcate;

import genomeEnums.Chrom;

public class AllelicBiasTableDummy implements AllelicBiasTable {

	private static long NumSitesRegisteredAtPosition = 2;
	private static float DefaultAverageVAF = 0.5f; 
	
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
	public void registerSite(Chrom chrom, int position, float vafNormal) {
		// Do nothing
	}

}
