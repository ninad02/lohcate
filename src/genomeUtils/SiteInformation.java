package genomeUtils;

import genomeEnums.Chrom;

// ========================================================================
public interface SiteInformation {
	public int getPosition();
	public Chrom getChrom();
	
	// ========================================================================
	public static interface Writeable extends SiteInformation {
		public void set(Chrom chrom, int position);
	}
}