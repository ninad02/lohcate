package lohcate;
import genomeEnums.Chrom;
import genomeUtils.RegionRange;
import lohcateEnums.ClusterType;
import nutils.counter.BucketCounter;


// ========================================================================
/** Given an individual's copy number regions by chromosome, this partitions the
 *  regions into three such objects based on the clustering type (Dup, LOH, HET)
 */
public class CopyNumberRegionRange extends RegionRange {
	public ClusterType mCopyNumberClusterType;
	public float mRecurrenceScore;
	protected BucketCounter mClusterTypeCounts;
	
	public CopyNumberRegionRange(ClusterType copyNumberClusterType, Chrom chrom, int regionStart) {
		super(chrom, regionStart);
		mCopyNumberClusterType = copyNumberClusterType;
		mRecurrenceScore = 1.0f;
		mClusterTypeCounts = new BucketCounter(ClusterType.getNumClusterTypes(), 0);
	}		
	
	public CopyNumberRegionRange(ClusterType copyNumberClusterType, Chrom chrom, int regionStart, int regionEnd) {
		super(chrom, regionStart, regionEnd);
		mCopyNumberClusterType = copyNumberClusterType;
		mRecurrenceScore = 1.0f;
		mClusterTypeCounts = new BucketCounter(ClusterType.getNumClusterTypes(), 0);
	}
	
	public CopyNumberRegionRange(CopyNumberRegionRange rhs) {
		super(rhs);
		mCopyNumberClusterType = rhs.mCopyNumberClusterType;
		mRecurrenceScore       = rhs.mRecurrenceScore;
		mClusterTypeCounts     = rhs.mClusterTypeCounts.getCopy();
	}
	
	public CopyNumberRegionRange getCopy() { return new CopyNumberRegionRange(this); }
}