package lohcate;
import genomeEnums.Chrom;

public interface AllelicBiasTable {
	
	public static final int ChromPositionDoesNotExist = -1;
	public static final int InsufficientSites         = -2;
	
	// ========================================================================
	/** Given a chromosome and position, this returns the average VAF reported for that position. 
	 *  @return the average vaf, or -1 if it doesn't exist. */
	public float getAvgVAF(Chrom chrom, int position);

	// ========================================================================
	public long getNumSitesRegisteredAtPosition(Chrom chrom, int position);

	// ========================================================================
	/** Given a chromosome and a position, this returns the average absolute copy number 
	 *  across the matched normals reported for that position
	 *  @return the average absolute copy number of matched normals, or -1 if it doesn't exist. 
	 */
	public float getAvgAbsCopyNumber(Chrom chrom, int position);
	
	// ========================================================================
	/** Given a chromosome, position, and site threshold, this returns the average absolute copy number reported
	 *  for that position across the matched normals, so long as the number of sites that contribute to the 
	 *  absolute average copy number meets or exceeds the threshold (inclusive).
	 *  @return the average absolute copy number, ChromPositionDoesNotExist if the position does not exist, or 
	 *  InsufficientSites if not enough sites meet the threshold
	 */
	public float getAvgAbsCopyNumber(Chrom chrom, int position, int siteThresholdInclusive);
	
	// ========================================================================
	/** Given a chromosome, position, and site threshold, this returns the average VAF reported for that position
	 *  so long as the number of sites that contribute to the average VAF meet or exceed the threshold (inclusive).
	 *  @return the average VAF, -1 if the position does not exist, or -2 if not enough sites meet the threshold.
	 */
	public float getAvgVAF(Chrom chrom, int position, int siteThresholdInclusive);

	// ========================================================================
	/** Registers a site with a normal vaf value. If a site already exists, the vaf value is weighted-averaged with the existing vaf value. */
	public void registerSite(Chrom chrom, int position, float vafNormal, float absoluteCopyNumber);
	
}
