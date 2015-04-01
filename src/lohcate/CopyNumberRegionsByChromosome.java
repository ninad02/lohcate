package lohcate;

import genomeEnums.Chrom;
import genomeUtils.RegionRangesOverGenome;

import java.io.PrintStream;
import java.util.ListIterator;

import lohcateEnums.EventType;
import nutils.StringUtils.StringClone;

// ========================================================================
// INNER CLASS
// ========================================================================
/** Stores the regions for a given sample, split by chromosome.  
 *  Chromosomes are indexed starting at 1 */
public class CopyNumberRegionsByChromosome extends RegionRangesOverGenome<CopyNumberRegionRangeLOHcate, StringClone> {
	
	// ========================================================================
	public CopyNumberRegionsByChromosome(String sampleName) {
		super(new StringClone(sampleName));
	}
	
	// ========================================================================
	public CopyNumberRegionsByChromosome(CopyNumberRegionsByChromosome rhs) {
		super(rhs);
	}
	
	// ========================================================================
	public CopyNumberRegionsByChromosome(CopyNumberRegionsByChromosome rhs, boolean deepCopy) {
		super(rhs, deepCopy);
	}
	
	// ========================================================================
	public String getSampleName() { return this.mMetaInfo.mStr; }

	// ========================================================================
	public void clearClusterCounts() {
		for (Chrom chrom : Chrom.values()) {			
			for (CopyNumberRegionRangeLOHcate oneRegionInChrom : getRegions(chrom)) {
				oneRegionInChrom.mEventTypeCounts.clear();
			}
		}
	}

	// ========================================================================
	@Override
	/** Returns an exact replica of this entire object, including copies of any contained regions. */ 
	public CopyNumberRegionsByChromosome makeClone() { return this.makeClone(true); }

	// ========================================================================
	@Override
	public CopyNumberRegionsByChromosome makeClone(boolean deepCopy) { return new CopyNumberRegionsByChromosome(this, deepCopy); }

	// ========================================================================
	/** Prints the contents to the PrintStream. */
	public void print(PrintStream out, String delimiter) {
		StringBuilder sb = new StringBuilder(65536);

		for (Chrom chrom : Chrom.values()) {
			for (CopyNumberRegionRangeLOHcate cnrr : mRegionsByChrom.get(chrom)) {
				sb.setLength(0);
				sb.append(cnrr.mCopyNumberEventType);
				sb.append(delimiter).append(cnrr.mRecurrenceScore);
				sb.append(delimiter).append(cnrr.toString());
				sb.append(delimiter).append(cnrr.getChromosome().ordinal());
				sb.append(delimiter).append(cnrr.getRangeStart());
				sb.append(delimiter).append(cnrr.getRangeEnd());
				sb.append(delimiter).append(cnrr.getRangeLength());
				
				double densityClusterType =        ((double) cnrr.mEventTypeCounts.getCount(cnrr.mCopyNumberEventType) / (double) cnrr.getRangeLength());
				double densityClusterHetGermline = ((double) cnrr.mEventTypeCounts.getCount(EventType.HETGermline)     / (double) cnrr.getRangeLength());
				sb.append(delimiter).append(densityClusterType);
				sb.append(delimiter).append(densityClusterHetGermline);
				
				cnrr.mEventTypeCounts.constructString(sb, false, delimiter);		
				out.println(sb.toString());
			}
		}
	}		
	
	// ========================================================================
	public void removeSingletonRegions() {
		for (Chrom chrom : Chrom.values()) {
			for (ListIterator<CopyNumberRegionRangeLOHcate> iter = getIteratorForChromosome(chrom); iter.hasNext(); ) {			
				if (iter.next().spansOneSite()) {
					iter.remove();
				}
			}
		}
	}
}

// ========================================================================