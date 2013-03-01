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
public class ChromPositionTracker extends ObjectWalkerTracker<Chrom> {
		
	protected int   mPositionPrev;
	protected long  mPositionGenomeWide;
	
	// ========================================================================
	public ChromPositionTracker() {
		super(Chrom.c0, null);
		clear();
	}
	
	// ========================================================================
	public void clear() {
		super.clear();
		mPositionPrev = 0;
		mPositionGenomeWide = 0;
	}
	
	// ========================================================================
	/** Returns whether a chromosome boundary has been crossed. */
	public boolean chromCrossedWithCurrentCoordinates(Chrom chromCurrent, int positionCurrent) {
		
		ChangeType chromChangeType = hasChanged(chromCurrent, true);
		
		if (chromChangeType.hasChanged()) {
			mPositionGenomeWide += positionCurrent;
		} else {
			mPositionGenomeWide += (positionCurrent - mPositionPrev);
		}
		
		mPositionPrev = positionCurrent;		
		return chromChangeType.hasChanged();
	}

	// ========================================================================
	public long getPositionGenomeWide() { return mPositionGenomeWide; }
	
}
