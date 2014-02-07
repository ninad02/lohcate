package lohcate;
import genomeEnums.Chrom;
import genomeUtils.GenomeConstants;
import genomeUtils.RegionRange;
import lohcateEnums.EventType;
import nutils.NullaryClassFactory;
import nutils.counter.BucketCounterEnum;


// ========================================================================
/** Given an individual's copy number regions by chromosome, this partitions the
 *  regions into three such objects based on the clustering type (Dup, LOH, HET)
 */
public class CopyNumberRegionRange extends RegionRange {
	
	// Nullary
	public static NullaryClassFactory<CopyNumberRegionRange> ClassFactory = new NullaryClassFactory<CopyNumberRegionRange>(CopyNumberRegionRange.class);
	
	// Member variables
	public EventType mCopyNumberEventType;
	public float mRecurrenceScore;
	public double mCopyNumber;
	protected BucketCounterEnum<EventType> mClusterTypeCounts;
	
	public CopyNumberRegionRange(EventType eventType, Chrom chrom, int regionStart) {
		super(chrom, regionStart);
		constructorCommon(eventType);
	}		
	
	public CopyNumberRegionRange(EventType eventType, Chrom chrom, int regionStart, int regionEnd) {
		super(chrom, regionStart, regionEnd);
		constructorCommon(eventType);
	}
	
	private void constructorCommon(EventType eventType) {
		mCopyNumberEventType = eventType;
		mRecurrenceScore = 1.0f;
		mCopyNumber = GenomeConstants.DefaultDiploidCopyNumber;
		mClusterTypeCounts = new BucketCounterEnum<EventType>(EventType.class);				
	}
	
	public CopyNumberRegionRange(CopyNumberRegionRange rhs) {
		super(rhs);
		mCopyNumberEventType = rhs.mCopyNumberEventType;
		mRecurrenceScore       = rhs.mRecurrenceScore;
		mCopyNumber            = rhs.mCopyNumber;
		mClusterTypeCounts     = rhs.mClusterTypeCounts.getCopy();
	}
	
	public CopyNumberRegionRange() {
		super();
		constructorCommon(EventType.Ignored);
	}
	
	public CopyNumberRegionRange getCopy() { return new CopyNumberRegionRange(this); }
}