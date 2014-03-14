package lohcate;
import genomeEnums.Chrom;

public interface AllelicBiasTable {
	
	// ========================================================================
	/** Given a chromosome and position, this returns the average VAF reported for that position. 
	 *  @return the average vaf, or -1 if it doesn't exist. */
	public float getAvgVAF(Chrom chrom, int position);

	// ========================================================================
	public long getNumSitesRegisteredAtPosition(Chrom chrom, int position);

	// ========================================================================
	/** Given a chromosome, position, and site threshold, this returns the average VAF reported for that position
	 *  so long as the number of sites that contribute to the average VAF meet or exceed the threshold (inclusive).
	 *  @return the average VAF, -1 if the position does not exist, or -2 if not enough sites meet the threshold.
	 */
	public float getAvgVAF(Chrom chrom, int position, int siteThresholdInclusive);

	// ========================================================================
	/** Registers a site with a normal vaf value. If a site already exists, the vaf value is weighted-averaged with the existing vaf value. */
	public void registerSite(Chrom chrom, int position, float vafNormal);
	
}
