package lohcate;

import genomeEnums.Chrom;
import genomeUtils.RegionRangesOverGenome;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.ListIterator;

import lohcateEnums.EventType;
import nutils.EnumMapSafe;

// ========================================================================
// INNER CLASS
// ========================================================================
/** Stores the regions for a given sample, split by chromosome.  
 *  Chromosomes are indexed starting at 1 */
public class CopyNumberRegionsByChromosome {
	
	// ========================================================================
	// MEMBER VARIABLES
	private String mSampleName;
	private EnumMapSafe<Chrom, ArrayList<CopyNumberRegionRangeLOHcate> > mRegionsByChrom;		

	// ========================================================================
	public CopyNumberRegionsByChromosome(String sampleName) {
		mRegionsByChrom = EnumMapSafe.createEnumMapOfArrayLists(Chrom.class, CopyNumberRegionRangeLOHcate.class);
		mSampleName = sampleName;
	}
	
	// ========================================================================
	public CopyNumberRegionsByChromosome(CopyNumberRegionsByChromosome rhs) {
		this(rhs.mSampleName);
		
		// Now deep copy the regions
		for (Chrom chrom : Chrom.values()) {
			ArrayList<CopyNumberRegionRangeLOHcate> cnrrListRhs = rhs.mRegionsByChrom.get(chrom);
			ArrayList<CopyNumberRegionRangeLOHcate> cnrrList    =     mRegionsByChrom.get(chrom); 
			
			for (CopyNumberRegionRangeLOHcate cnrrRhs : cnrrListRhs) {
				cnrrList.add(cnrrRhs.getCopy());					
			}
		}
	}
	
	// ========================================================================
	public String getSampleName() { return mSampleName; }
	
	// ========================================================================
	public void addRegion(Chrom chrom, CopyNumberRegionRangeLOHcate region) {
		mRegionsByChrom.get(chrom).add(region);
	}
	
	// ========================================================================
	public void addRegion(CopyNumberRegionRangeLOHcate region) {
		addRegion(region.getChromosome(), region);
	}
	
	// ========================================================================
	public void clearClusterCounts() {
		for (Chrom chrom : Chrom.values()) {			
			for (CopyNumberRegionRangeLOHcate oneRegionInChrom : mRegionsByChrom.get(chrom)) {
				oneRegionInChrom.mEventTypeCounts.clear();
			}
		}
	}

	// ========================================================================
	/** Returns the list of regions for a chromosome. */
	public ArrayList<CopyNumberRegionRangeLOHcate> getRegions(Chrom chrom) { return mRegionsByChrom.get(chrom); }
	
	// ========================================================================
	/** Returns an exact replica of this entire object, including copies of any contained regions. */ 
	public CopyNumberRegionsByChromosome getCopy() { return new CopyNumberRegionsByChromosome(this); }
	
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
	/** Given a chromosome and position, returns the region that includes the coordinate, or 
	 *  null if no region includes the coordinate.
	 */
	public CopyNumberRegionRangeLOHcate getRegion(Chrom chrom, int position) {
		int resultIndex = getIndexOfRegion(chrom, position);
		return ((resultIndex < 0) ? null : mRegionsByChrom.get(chrom).get(resultIndex));
	}
	
	
	// ========================================================================
	/** Returns an iterator for the regions for a chromosome. */
	public ListIterator<CopyNumberRegionRangeLOHcate> getIteratorForChromosome(Chrom chrom) {
		ArrayList<CopyNumberRegionRangeLOHcate> regionsForChrom = mRegionsByChrom.get(chrom);
		return regionsForChrom.listIterator();		
	}
	
	// ========================================================================
	/** Given a chromosome and position, returns the index of the region that includes the coordinate, or 
	 *  -1 if no region includes the coordinate.
	 */
	public int getIndexOfRegion(Chrom chrom, int position) {
		ArrayList<CopyNumberRegionRangeLOHcate> regions = mRegionsByChrom.get(chrom);
		for (int i = 0; i < regions.size(); i++) {			
			if (regions.get(i).inRange(chrom, position)) {
				return i; 
			} else if (regions.get(i).beforeRange(chrom, position)) {
				break;
			}
		}
		
		return -1 ;
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
	
	// ========================================================================
	public static class StringClone implements Cloneable {
		String mStr;
		public StringClone(String s) { mStr = s; }
		public Object clone() { return new String(mStr); }
	}
}

// ========================================================================