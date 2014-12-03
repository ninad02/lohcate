package genomeUtils;

import genomeEnums.Chrom;

import java.util.Comparator;


public class RegionRange {
	
	// ========================================================================
	// MAIN CLASS
	// ========================================================================
	Chrom mChrom;
	boolean mRangeFinalized;
	int mRangeStart;
	int mRangeEnd;		
	int mNumSitesInterrogated;
	
	public RegionRange() {
		set(Chrom.c0, 0);
	}
	
	public RegionRange(Chrom chrom, int rangeStart) {
		set(chrom, rangeStart);		
	}
	
	public RegionRange(Chrom chrom, int rangeStart, int rangeEnd) {
		this(chrom, rangeStart);		
		extendRange(chrom, rangeEnd);
	}
	
	public RegionRange(RegionRange rhs) {
		set(rhs.mChrom, rhs.mRangeStart, rhs.mRangeEnd, rhs.mRangeFinalized, rhs.mNumSitesInterrogated);
	}
	
	public RegionRange getCopy() { return new RegionRange(this); }
	
	/** Returns the range start. */
	public int getRangeStart() { return mRangeStart; }
	
	/** Returns the range end. */
	public int getRangeEnd() { return mRangeEnd; }
	
	/** Returns the range length, inclusive of start and end. */
	public int getRangeLength() { return mRangeEnd - mRangeStart + 1; }
	
	/** Returns whether the range spans only one site (in other words, has range length of 1. */
	public boolean spansOneSite() { return (getRangeLength() == 1); }
	
	/** Returns the chromosome. */
	public Chrom getChromosome() { return mChrom; }
	
	/** Increments the number of sites interrogated by the given amount. */
	public void incrementSitesInterrogated(int numSitesToAdd) { mNumSitesInterrogated += numSitesToAdd; }
	
	/** Returns the number of sites interrogated in the region. */
	public int getNumSitesInterrogated() { return mNumSitesInterrogated; }
	
	/** Makes the range finalized. */
	public void makeFinalized() { mRangeFinalized = true; }
	
	/** Returns whether this range is finalized or not. */
	public boolean isFinalized() { return mRangeFinalized; }
	
	/** Sets the range to the given chrom and position.  Automatically sets the range end. */
	public void set(Chrom chrom, int rangeStart) {
		set(chrom, rangeStart, rangeStart, false, 1);
	}
	
	/** Sets the range. */
	public void set(Chrom chrom, int rangeStart, int rangeEnd, boolean makeFinalized, int numSitesInterrogated) {
		if (rangeEnd < rangeStart) {
			throwErrorAndExit("ERROR: RegionRange.set(): Range end must follow range start!\tStart:\t" + rangeStart + "\tEnd:\t" + rangeEnd);		
		}
		
		mChrom = chrom;		
		mRangeStart = rangeStart;
		mRangeEnd   = rangeEnd;  // need to do this before extending the rnage
		mRangeFinalized = makeFinalized;
		mNumSitesInterrogated = numSitesInterrogated;
	}
	
	public void set(RegionRange rhs) {
		set(mChrom, mRangeStart, mRangeEnd, mRangeFinalized, mNumSitesInterrogated);
	}
	
	/** Sets the range start. */
	public void setRangeStart(int rangeStart) {
		if (rangeStart > mRangeEnd) throwErrorAndExit("ERROR: RegionRange.setRangeStart() Range start must precede range end!");
		mRangeStart = rangeStart;
	}
	
	/** Sets the range end.  */
	public void setRangeEnd(int rangeEnd) {
		if (rangeEnd < mRangeStart) throwErrorAndExit("ERROR: RegionRange.setRangeEnd(): Range end must follow range start!");
		mRangeEnd = rangeEnd;		
	}
	
	/** Given a chromNum and a position, this checks whether the chrom
	 *  number and the position can be used to extend the range.  First, 
	 *  this checks that the chrom number and position do not already 
	 *  exist in the range.  If so, and the chrom number is appropriate,
	 *  then it adds it to the range, and true is returned.  Else, false. 
	 */
	public boolean extendRange(Chrom chrom, int position) {
		//System.out.println(mRangeFinalized + "\t" + chrom + "\t" + mChrom + "\t" + position + "\t" + mRangeEnd);
		if (mRangeFinalized || (chrom != mChrom) || (position <= mRangeEnd)) {
			makeFinalized();
			return false;
		}
				
		mRangeEnd = position;
		mNumSitesInterrogated++;
		return true;
	}
	
	/** Given a chromNum and a position, this returns whether in range. */
	public boolean inRange(Chrom chrom, int position) {
		return (mChrom.equals(chrom) && ((position >= mRangeStart) && (position <= mRangeEnd)));
	}
	
	/** Given a chromNum and a position, this returns whether it is after the range. */
	public boolean afterRange(Chrom chrom, int position) {
		return (chrom.greaterThan(mChrom) || (mChrom.equals(chrom) && (position > mRangeEnd)));
	}
	
	/** Given a chromNum and a position, this returns whether it is before the range. */
	public boolean beforeRange(Chrom chrom, int position) {
		return (chrom.lessThan(mChrom) || (mChrom.equals(chrom) && (position < mRangeStart)));
	}
	
	/** Given another region, this returns whether there is any overlap.  */
	public boolean overlapRange(RegionRange rhs) {
		return testAndCharacterizeOverlap(rhs).isOverlapType();		
	}
			
	/** Given another region, this returns the type of overlap, if any.  
	 *  The return value is with respect to "this" object, not the argument object passed in.
	 */
	public RegionRangeOverlap testAndCharacterizeOverlap(RegionRange rhs) {
		if (mChrom.lessThan   (rhs.mChrom)) return RegionRangeOverlap.BeforeViaDiffChromosome;
		if (mChrom.greaterThan(rhs.mChrom)) return RegionRangeOverlap.AfterViaDiffChromosome;
		
		// Now we are guaranteeing we are on the same chromosome
		if (this.mRangeStart == rhs.mRangeStart) {
			// There are three possibilities: Equals, SubsumesAlignedLeft, or ConsumedByAlignedLeft
			if        (this.mRangeEnd == rhs.mRangeEnd) {
				return RegionRangeOverlap.Equals;
			} else if (this.mRangeEnd <  rhs.mRangeEnd) {
				return RegionRangeOverlap.ConsumedByAlignedLeft;
			} else {
				return RegionRangeOverlap.SubsumesAlignedLeft;
			}
			
		} else if (this.mRangeStart < rhs.mRangeStart) {
			// Possibiilites are: SubsumesTotal, SubsumesAlignedRight, BeforeWithOverlap, AdjacentBefore, BeforeWithoutOverlap
			if        (this.mRangeEnd == rhs.mRangeEnd) {
				return RegionRangeOverlap.SubsumesAlignedRight;
			} else if (this.mRangeEnd >  rhs.mRangeEnd) {
				return RegionRangeOverlap.SubsumesTotal;
			} else {
				if (this.mRangeEnd >= rhs.mRangeStart) {
					return RegionRangeOverlap.BeforeWithOverlap;					
				} else {
					int distance = rhs.mRangeStart - this.mRangeEnd;
					return (distance == 1) ? RegionRangeOverlap.AdjacentBefore : RegionRangeOverlap.BeforeWithoutOverlap;
				}				
			}
		} else {
			// Possibilities are: ConsumedByTotal, ConsumedByAlignedRight, AfterWithOverlap, AdjacentAfter, AfterWithoutOverlap
			if        (this.mRangeEnd == rhs.mRangeEnd) {
				return RegionRangeOverlap.ConsumedByAlignedRight;
			} else if (this.mRangeEnd <  rhs.mRangeEnd) {
				return RegionRangeOverlap.ConsumedByTotal;
			} else {
				if (this.mRangeStart <= rhs.mRangeEnd) {
					return RegionRangeOverlap.AfterWithOverlap;
				} else {
					int distance = this.mRangeStart - rhs.mRangeEnd;
					return (distance == 1) ? RegionRangeOverlap.AdjacentAfter : RegionRangeOverlap.AfterWithoutOverlap;
				}
			}
		}		
	}
	
	// ========================================================================
	//@Override
	public int compareTo(RegionRange rhs) {
		RegionRangeOverlap result = testAndCharacterizeOverlap(rhs);
		
		switch(result) {
		case Equals: 
			return 0;
		case BeforeWithoutOverlap: case BeforeViaDiffChromosome: case AdjacentBefore: 
			return -1;
		case SubsumesTotal: case SubsumesAlignedRight: case BeforeWithOverlap: case ConsumedByAlignedLeft:
			return -1;		
		case AfterWithoutOverlap:  case AfterViaDiffChromosome:  case AdjacentAfter:  
			return 1;
		case ConsumedByTotal: case ConsumedByAlignedRight: case AfterWithOverlap: case SubsumesAlignedLeft:
			return 1;
		default:
			throwErrorAndExit("ERROR: Impossible state!");
			return 0;
		}		
	}

	// ========================================================================
	public static void throwErrorAndExit(String errorString) {
		Exception e = new Exception(errorString);
		e.printStackTrace();
		System.exit(-1);	
	}
	
	// ========================================================================
	// INNER CLASS
	// ========================================================================
	/**
	 * We go through various examples, with a region represented by a series of x's: xxxxxxx and gaps by _
	 * 
	 * Equals:                    xxxxxxx
	 *                            xxxxxxx
	 * 
	 * SubsumesTotal              xxxxxxx
	 *                             xxxxx
	 *                             
	 * SubsumesAlignedLeft        xxxxxxx
	 *                            xxxxxx
	 * 
	 * SubsumesAlignedRight       xxxxxxx
	 *                             xxxxxx
	 *                             
	 * Consumed By Total:          xxxxx
	 *                            xxxxxxx
	 *                            
	 * ConsumedByAlignedLeft      xxxxxx
	 *                            xxxxxxx
	 *                            
	 * ConsumedByAlignedRight      xxxxxx
	 *                            xxxxxxx
	 *                            
	 * Before With Overlap:       xxxxxxx            xxxxxxx
	 *                             xxxxxxx                 xxxxxxx
	 *                                  
	 * After With Overlap:         xxxxxxx                 xxxxxxx
	 *                            xxxxxxx            xxxxxxx
	 *
	 * Adjacent Before:           xxxxxxx
	 *                                   xxxxxxx
	 *
	 * Adjacent After:                   xxxxxxx
	 *                            xxxxxxx
	 *                            
	 * Before Without Overlap:    xxxxxxx_
	 *                                   _xxxxxxx
	 *                                   
	 * After Without Overlap:            _xxxxxxx
	 *                            xxxxxxx_
	 *                             
	 * Before Via Diff Chromosome: chromosome # of "this" region less    than that of other region
	 * After  Via Diff Chromosome: chromosome # of "this" region greater than that of other region
	 *   
	 * @author Ninad Dewal
	 *
	 */
	public static enum RegionRangeOverlap {
		// Overlap Types
		Equals,                    
		SubsumesTotal,           
		SubsumesAlignedLeft,     
		SubsumesAlignedRight,    
		ConsumedByTotal,
		ConsumedByAlignedLeft,   
		ConsumedByAlignedRight,
		BeforeWithOverlap,       
		AfterWithOverlap,
		
		// Non-overlap types
		BeforeWithoutOverlap,    
		AfterWithoutOverlap,
		BeforeViaDiffChromosome, 
		AfterViaDiffChromosome,  
		AdjacentBefore,          
		AdjacentAfter;
		
		protected static final RegionRangeOverlap[] RangeEndComparisonWithAlignedLeft = new RegionRangeOverlap[] {  
			ConsumedByAlignedLeft, Equals, SubsumesAlignedLeft
			// to correspond to -1, 0, 1 from a right end alignment comparison
		};
		
		public boolean isOverlapType() { return (ordinal() < BeforeWithoutOverlap.ordinal()); }
		
	}
	
	// ========================================================================
	
	
	/** Returns this as a string. */
	public String toString() {
		StringBuilder sb = new StringBuilder(512);
		sb.append(mChrom.getCode());
		sb.append(":");
		sb.append(mRangeStart);
		sb.append("-");
		sb.append(mRangeEnd);		
		return sb.toString();
	}
	
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
	}
	
	// ========================================================================
	// INNER CLASS
	// ========================================================================
	
	public static class RegionRangeComparator implements Comparator<RegionRange> {
		
		public static final RegionRangeComparator TheComparator = new RegionRangeComparator();
			
		public int compare(RegionRange ar1, RegionRange ar2) {
			if (ar1.mRangeStart == ar1.mRangeEnd) {
				return compareHelper(ar1, ar2);
			} else {
				return (-1 * compareHelper(ar2, ar1));
			}
		}
		
		private int compareHelper(RegionRange ar1, RegionRange ar2) {
			//System.out.println(ar1.toString() + ar2.toString());
			
			// First, compare the chromosomes. 
			if (ar1.mChrom.getCode() < ar2.mChrom.getCode()) {
				return -1;
			} else if (ar1.mChrom.getCode() > ar2.mChrom.getCode()) {
				return 1;
			} else {				
				if (ar1.mRangeStart < ar2.mRangeStart) {
					return -1;
				} else {
					// means ar1.mRangeStart >= ar2.mRangeStart
					// However, we do not match if it is so large that it crosses the end					
					if (ar1.mRangeStart > ar2.mRangeEnd) {
						return 1;
					} else {
						// Means ar1.mRangeStart <= ar2.mRangeEnd
						// If we are using this class, then it means it is for
						// searching in a list of ranges.  This class was made
						// with the assumption that ar1.mRangeStart == 
						// ar1.mRangeEnd (its range was of size 1).  Thus, we
						// return 0 here
						return 0;						
					}
				}
			}
		}
	}

	
	// ========================================================================
	// END INNER CLASS
	// ========================================================================

}


