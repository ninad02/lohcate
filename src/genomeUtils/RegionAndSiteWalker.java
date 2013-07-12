package genomeUtils;

import genomeEnums.Chrom;

import java.io.BufferedWriter;
import java.util.ArrayList;
import java.util.ListIterator;

import lohcate.CopyNumberRegionRange;
import lohcateEnums.EventType;
import nutils.CompareUtils;
import nutils.IOUtils;

/**
 * 
 * @author Ninad Dewal
 *
 */
public class RegionAndSiteWalker {

	// ========================================================================
	// ========================================================================
	public static interface Actioner {
		
		public void takeAction(boolean rangeInTargetSet, RegionRange region);
	}

	// ========================================================================
	// ========================================================================
	RegionRange dummyRange = new RegionRange(Chrom.c0, 0);

	// ========================================================================
	public void walk(ArrayList<RegionRange> regionsInChrom, SNVMap snvMap, Chrom chrom, Actioner actioner) {

		clear();
		int numSitesOnChrom = snvMap.getNumSitesOnChromosome(chrom);
		if (numSitesOnChrom <= 0) return;

		ListIterator<RegionRange> regionsInChromIter = regionsInChrom.listIterator();
		RegionRange currentRegion = null;		

		int indexInMap = 0;
		boolean dummyRangeValid = false;

		// First scan and see if we have any regions preceding the first chromosomal position in the map
		while (regionsInChromIter.hasNext()) {
			currentRegion = regionsInChromIter.next();
			if (currentRegion.afterRange(chrom, snvMap.getPosition(chrom, indexInMap))) {
				actioner.takeAction(false, currentRegion);					
				currentRegion = null;
			} else {
				break;
			}
		}

		// We could have broken the loop because the current region encloses the

		for (indexInMap = 0; indexInMap < numSitesOnChrom; indexInMap++) {
			int mapPosition = snvMap.getPosition(chrom, indexInMap);

			if ((currentRegion == null) || currentRegion.beforeRange(chrom, mapPosition)) {
				if (dummyRangeValid) {
					if (dummyRange.isFinalized()) {
						actioner.takeAction(false, dummyRange);							
						dummyRange.set(chrom, mapPosition);
					} else {
						boolean extendResult = dummyRange.extendRange(chrom, mapPosition);
						if (!extendResult) {
							String errorString = "ERROR: printSegmentedRegionsToFile(): Could not extend range! " + chrom + "\t" + mapPosition + "\t" + dummyRange;
							CompareUtils.ensureTrue(extendResult, errorString);	
						}

					}
				} else {
					dummyRange.set(chrom, mapPosition);
					dummyRangeValid = true;
				}

			} else if (currentRegion.inRange(chrom, mapPosition)) {
				if (dummyRangeValid) {
					actioner.takeAction(false, dummyRange);						
					dummyRangeValid = false;
				}

				// Print this region to file 
				int indexOfRegionStartInMap = snvMap.getIndexOfPositionInMap(chrom, currentRegion.getRangeStart());
				if (indexOfRegionStartInMap < 0) {
					String errorString = "ERROR: printSegmentedRegionsToFile(): Region start (" + currentRegion.getRangeStart() + "must exist in map!";
					CompareUtils.ensureTrue(false, errorString);
				}

				// Now change the loop index to move to the map position just after the end of this region
				int indexOfRegionEndInMap = snvMap.getIndexOfPositionInMap(chrom, currentRegion.getRangeEnd());
				if (indexOfRegionEndInMap < 0) {
					String errorString = "ERROR: printSegmentedRegionsToFile(): Region end (" + currentRegion.getRangeEnd() + "must exist in map!";
					CompareUtils.ensureTrue(false, errorString);						
				}
				indexInMap = indexOfRegionEndInMap;  // indexInMap will be incremented at loop end

				currentRegion.set(chrom, currentRegion.getRangeStart(), currentRegion.getRangeEnd(), true, indexOfRegionEndInMap - indexOfRegionStartInMap + 1);
				actioner.takeAction(true, currentRegion);								

				// Now move to the next cna-affected region
				currentRegion = regionsInChromIter.hasNext() ? regionsInChromIter.next() : null;

			} else if (currentRegion.afterRange(chrom, mapPosition)) {
				CompareUtils.ensureTrue(false, "ERROR: printSegmentedRegionsToFile(): Impossible for region to precede position!");
			} else {
				CompareUtils.ensureTrue(false, "ERROR: printSegmentedRegionsToFile(): Impossible state!");
			}
		}

		if (dummyRangeValid) {
			actioner.takeAction(false, dummyRange);				
		}
	}

	// ========================================================================
	public void clear() {
		dummyRange.set(Chrom.c0, 0, 0, false, 1);		
	}

	// ========================================================================
	public RegionAndSiteWalker() {
		dummyRange = new RegionRange(Chrom.c0, 0);
		// TODO Auto-generated constructor stub
	}

	// ========================================================================
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub

	}

}
