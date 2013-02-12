package genomeUtils;

import genomeEnums.Chrom;

/**
 * Oftentimes, during a traversal of a list of genomic sites, the chromosome
 * boundaries need to be identified.  This convenience class facilities the
 * routine pattern of checking whether a chromosome boundary has been crossed
 * during the iteration.  
 * 
 * @author Ninad Dewal
 *
 */
public class ChromPositionTracker {
	
	protected Chrom mChromPrev;
	protected int   mPositionPrev;
	protected long  mPositionGenomeWide;
	
	// ========================================================================
	public ChromPositionTracker() {
		clear();
	}
	
	// ========================================================================
	public void clear() {
		mChromPrev = Chrom.c0;;
		mPositionPrev = 0;
		mPositionGenomeWide = 0;
	}
	
	// ========================================================================
	/** Returns whether a chromosome boundary has been crossed. */
	public boolean chromCrossedWithCurrentCoordinates(Chrom chromCurrent, int positionCurrent) {
		boolean chromCrossed = false;
		
		if (chromCurrent == mChromPrev) {			
			mPositionGenomeWide += (positionCurrent - mPositionPrev);
		} else {
			chromCrossed = true;
			mPositionGenomeWide += positionCurrent;
		}
		
		mPositionPrev = positionCurrent;
		mChromPrev = chromCurrent;
		return chromCrossed;
	}

	// ========================================================================
	public long getPositionGenomeWide() { return mPositionGenomeWide; }
	
}
