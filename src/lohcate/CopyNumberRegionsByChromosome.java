package lohcate;

import genomeEnums.Chrom;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.ListIterator;

import lohcateEnums.ClusterType;
import nutils.ArrayUtils;
import nutils.EnumMapSafe;

// ========================================================================
// INNER CLASS
// ========================================================================
/** Stores the regions for a given sample, split by chromosome.  
 *  Chromosomes are indexed starting at 1 */
public class CopyNumberRegionsByChromosome {
	private static final ArrayList<CopyNumberRegionRange> dummyListForChrom0 = new ArrayList<CopyNumberRegionRange>();
	
	// ========================================================================
	// MEMBER VARIABLES
	String mSampleName;
	EnumMapSafe<Chrom, ArrayList<CopyNumberRegionRange> > mRegionsByChrom;		

	// ========================================================================
	public CopyNumberRegionsByChromosome(String sampleName) {
		mRegionsByChrom = ArrayUtils.createEnumMapOfArrayLists(Chrom.class, CopyNumberRegionRange.class);
		mSampleName = sampleName;
	}
	
	// ========================================================================
	public CopyNumberRegionsByChromosome(CopyNumberRegionsByChromosome rhs) {
		this(rhs.mSampleName);
		
		// Now deep copy the regions
		for (Chrom chrom : Chrom.values()) {
			ArrayList<CopyNumberRegionRange> cnrrListRhs = rhs.mRegionsByChrom.get(chrom);
			ArrayList<CopyNumberRegionRange> cnrrList    =     mRegionsByChrom.get(chrom); 
			
			for (CopyNumberRegionRange cnrrRhs : cnrrListRhs) {
				cnrrList.add(cnrrRhs.getCopy());					
			}
		}
	}
	
	// ========================================================================
	public void addRegion(Chrom chrom, CopyNumberRegionRange region) {
		mRegionsByChrom.get(chrom).add(region);
	}
	
	// ========================================================================
	public void addRegion(CopyNumberRegionRange region) {
		addRegion(region.getChromosome(), region);
	}
	
	// ========================================================================
	public void clearClusterCounts() {
		for (Chrom chrom : Chrom.values()) {			
			for (CopyNumberRegionRange oneRegionInChrom : mRegionsByChrom.get(chrom)) {
				oneRegionInChrom.mClusterTypeCounts.clear();
			}
		}
	}

	// ========================================================================
	/** Returns the list of regions for a chromosome. */
	public ArrayList<CopyNumberRegionRange> getRegions(Chrom chrom) { return mRegionsByChrom.get(chrom); }
	
	// ========================================================================
	/** Returns an exact replica of this entire object, including copies of any contained regions. */ 
	public CopyNumberRegionsByChromosome getCopy() { return new CopyNumberRegionsByChromosome(this); }
	
	// ========================================================================
	/** Prints the contents to the PrintStream. */
	public void print(PrintStream out, String delimiter) {
		StringBuilder sb = new StringBuilder(65536);

		for (Chrom chrom : Chrom.values()) {
			for (CopyNumberRegionRange cnrr : mRegionsByChrom.get(chrom)) {
				sb.setLength(0);
				sb.append(cnrr.mCopyNumberClusterType);
				sb.append(delimiter).append(cnrr.mRecurrenceScore);
				sb.append(delimiter).append(cnrr.toString());
				sb.append(delimiter).append(cnrr.getChromosome().ordinal());
				sb.append(delimiter).append(cnrr.getRangeStart());
				sb.append(delimiter).append(cnrr.getRangeEnd());
				sb.append(delimiter).append(cnrr.getRangeLength());
				
				double densityClusterType =        ((double) cnrr.mClusterTypeCounts.getCount(cnrr.mCopyNumberClusterType) / (double) cnrr.getRangeLength());
				double densityClusterHetGermline = ((double) cnrr.mClusterTypeCounts.getCount(ClusterType.HETGermline)     / (double) cnrr.getRangeLength());
				sb.append(delimiter).append(densityClusterType);
				sb.append(delimiter).append(densityClusterHetGermline);
				
				cnrr.mClusterTypeCounts.constructString(sb, false, delimiter);		
				out.println(sb.toString());
			}
		}
	}		
	
	// ========================================================================
	/** Given a chromosome and position, returns the region that includes the coordinate, or 
	 *  null if no region includes the coordinate.
	 */
	public CopyNumberRegionRange getRegion(Chrom chrom, int position) {
		int resultIndex = getIndexOfRegion(chrom, position);
		return ((resultIndex < 0) ? null : mRegionsByChrom.get(chrom).get(resultIndex));
	}
	
	
	// ========================================================================
	/** Returns an iterator for the regions for a chromosome. */
	public ListIterator<CopyNumberRegionRange> getIteratorForChromosome(Chrom chrom) {
		ArrayList<CopyNumberRegionRange> regionsForChrom = mRegionsByChrom.get(chrom);
		return regionsForChrom.listIterator();		
	}
	
	// ========================================================================
	/** Given a chromosome and position, returns the index of the region that includes the coordinate, or 
	 *  -1 if no region includes the coordinate.
	 */
	public int getIndexOfRegion(Chrom chrom, int position) {
		ArrayList<CopyNumberRegionRange> regions = mRegionsByChrom.get(chrom);
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
			for (ListIterator<CopyNumberRegionRange> iter = getIteratorForChromosome(chrom); iter.hasNext(); ) {			
				if (iter.next().spansOneSite()) {
					iter.remove();
				}
			}
		}
	}
}

// ========================================================================