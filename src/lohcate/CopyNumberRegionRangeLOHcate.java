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
public class CopyNumberRegionRangeLOHcate extends RegionRange {
	
	// Nullary
	public static NullaryClassFactory<CopyNumberRegionRangeLOHcate> ClassFactory = new NullaryClassFactory<CopyNumberRegionRangeLOHcate>(CopyNumberRegionRangeLOHcate.class);
	
	// Member variables
	public EventType mCopyNumberEventType;
	public float mRecurrenceScore;
	public double mCopyNumber;
	protected BucketCounterEnum<EventType> mEventTypeCounts;
	
	public CopyNumberRegionRangeLOHcate(EventType eventType, Chrom chrom, int regionStart) {
		super(chrom, regionStart);
		constructorCommon(eventType);
	}		
	
	public CopyNumberRegionRangeLOHcate(EventType eventType, Chrom chrom, int regionStart, int regionEnd) {
		super(chrom, regionStart, regionEnd);
		constructorCommon(eventType);
	}
	
	private void constructorCommon(EventType eventType) {
		mCopyNumberEventType = eventType;
		mRecurrenceScore = 1.0f;
		mCopyNumber = GenomeConstants.DefaultDiploidCopyNumber;
		mEventTypeCounts = new BucketCounterEnum<EventType>(EventType.class);				
	}
	
	public CopyNumberRegionRangeLOHcate(CopyNumberRegionRangeLOHcate rhs) {
		super(rhs);
		mCopyNumberEventType = rhs.mCopyNumberEventType;
		mRecurrenceScore       = rhs.mRecurrenceScore;
		mCopyNumber            = rhs.mCopyNumber;
		mEventTypeCounts     = rhs.mEventTypeCounts.getCopy();
	}
	
	public CopyNumberRegionRangeLOHcate() {
		super();
		constructorCommon(EventType.Ignored);
	}
	
	public boolean shareSameEvent(CopyNumberRegionRangeLOHcate rhs) { return this.mCopyNumberEventType == rhs.mCopyNumberEventType; }
	
	public CopyNumberRegionRangeLOHcate getCopy() { return new CopyNumberRegionRangeLOHcate(this); }
}