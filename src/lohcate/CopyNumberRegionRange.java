package lohcate;
import genomeEnums.Chrom;
import genomeUtils.RegionRange;
import lohcateEnums.ClusterType;
import nutils.counter.BucketCounterEnum;


// ========================================================================
/** Given an individual's copy number regions by chromosome, this partitions the
 *  regions into three such objects based on the clustering type (Dup, LOH, HET)
 */
public class CopyNumberRegionRange extends RegionRange {
	public ClusterType mCopyNumberClusterType;
	public float mRecurrenceScore;
	protected BucketCounterEnum<ClusterType> mClusterTypeCounts;
	
	public CopyNumberRegionRange(ClusterType eventType, Chrom chrom, int regionStart) {
		super(chrom, regionStart);
		constructorCommon(eventType);
	}		
	
	public CopyNumberRegionRange(ClusterType eventType, Chrom chrom, int regionStart, int regionEnd) {
		super(chrom, regionStart, regionEnd);
		constructorCommon(eventType);
	}
	
	private void constructorCommon(ClusterType eventType) {
		mCopyNumberClusterType = eventType;
		mRecurrenceScore = 1.0f;
		mClusterTypeCounts = new BucketCounterEnum<ClusterType>(ClusterType.class);				
	}
	
	public CopyNumberRegionRange(CopyNumberRegionRange rhs) {
		super(rhs);
		mCopyNumberClusterType = rhs.mCopyNumberClusterType;
		mRecurrenceScore       = rhs.mRecurrenceScore;
		mClusterTypeCounts     = rhs.mClusterTypeCounts.getCopy();
	}
	
	public CopyNumberRegionRange getCopy() { return new CopyNumberRegionRange(this); }
}