package genomeUtils;

import genomeEnums.Chrom;

import java.util.ArrayList;
import java.util.ListIterator;

import nutils.CompareUtils;
import nutils.NullaryClassFactory;

/**
 * 
 * @author Ninad Dewal
 *
 */
public class RegionAndSiteWalker<T extends RegionRange<T>> {

	// ========================================================================
	// ========================================================================
	public static interface Actioner<T extends RegionRange<T>> {
		
		public void takeAction(ArrayList<Boolean> rangeInTargetSet, ArrayList<T> regions, T regionLatest);
	}

	// ========================================================================
	// ========================================================================
	T mDummyRange;
	NullaryClassFactory<T> mClassFactory;
	ArrayList<Boolean> mRangeInTargetSet;

	// ========================================================================
	public void walk(ArrayList<T> regionsInChrom, SNVMap snvMap, Chrom chrom, Actioner<T> actioner) {

		clear();
		int numSitesOnChrom = snvMap.getNumSitesOnChromosome(chrom);
		if (numSitesOnChrom <= 0) return;

		ListIterator<T> regionsInChromIter = regionsInChrom.listIterator();
		T currentRegion = null;		
		ArrayList<T> mergedRegionList = new ArrayList<T>(regionsInChrom.size());

		int indexInMap = 0;
		boolean dummyRangeValid = false;

		// First scan and see if we have any regions preceding the first chromosomal position in the map
		while (regionsInChromIter.hasNext()) {
			currentRegion = regionsInChromIter.next();
			if (currentRegion.afterRange(chrom, snvMap.getPosition(chrom, indexInMap))) {
				mRangeInTargetSet.add(Boolean.FALSE);
				mergedRegionList.add(currentRegion);
				actioner.takeAction(mRangeInTargetSet, mergedRegionList, currentRegion);					
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
					if (mDummyRange.isFinalized()) {
						mRangeInTargetSet.add(Boolean.FALSE);
						mergedRegionList.add(mDummyRange);
						actioner.takeAction(mRangeInTargetSet, mergedRegionList, mDummyRange);							
						createAndSetNewBufferRegion(chrom, mapPosition);
					} else {
						boolean extendResult = mDummyRange.extendRange(chrom, mapPosition);
						if (!extendResult) {
							String errorString = "ERROR:  Could not extend range! " + chrom + "\t" + mapPosition + "\t" + mDummyRange;
							CompareUtils.ensureTrue(extendResult, errorString);	
						}

					}
				} else {
					createAndSetNewBufferRegion(chrom, mapPosition);
					dummyRangeValid = true;
				}

			} else if (currentRegion.inRange(chrom, mapPosition)) {
				if (dummyRangeValid) {
					mRangeInTargetSet.add(Boolean.FALSE);
					mergedRegionList.add(mDummyRange);
					actioner.takeAction(mRangeInTargetSet, mergedRegionList, mDummyRange);						
					dummyRangeValid = false;
				}

				// Print this region to file 
				int indexOfRegionStartInMap = snvMap.getIndexOfPositionInMap(chrom, currentRegion.getRangeStart());
				if (indexOfRegionStartInMap < 0) {
					String errorString = "ERROR:  Region start (" + currentRegion.getRangeStart() + "must exist in map!";
					CompareUtils.ensureTrue(false, errorString);
				}

				// Now change the loop index to move to the map position just after the end of this region
				int indexOfRegionEndInMap = snvMap.getIndexOfPositionInMap(chrom, currentRegion.getRangeEnd());
				if (indexOfRegionEndInMap < 0) {
					String errorString = "ERROR:  Region end (" + currentRegion.getRangeEnd() + "must exist in map!";
					CompareUtils.ensureTrue(false, errorString);						
				}
				indexInMap = indexOfRegionEndInMap;  // indexInMap will be incremented at loop end

				currentRegion.set(chrom, currentRegion.getRangeStart(), currentRegion.getRangeEnd(), true, indexOfRegionEndInMap - indexOfRegionStartInMap + 1);
				mRangeInTargetSet.add(Boolean.TRUE);
				mergedRegionList.add(currentRegion);
				actioner.takeAction(mRangeInTargetSet, mergedRegionList, currentRegion);								

				// Now move to the next cna-affected region
				currentRegion = regionsInChromIter.hasNext() ? regionsInChromIter.next() : null;

			} else if (currentRegion.afterRange(chrom, mapPosition)) {
				CompareUtils.ensureTrue(false, "ERROR:  Impossible for region to precede position!");
			} else {
				CompareUtils.ensureTrue(false, "ERROR:  Impossible state!");
			}
		}

		if (dummyRangeValid) {
			mRangeInTargetSet.add(Boolean.FALSE);
			mergedRegionList.add(mDummyRange);
			actioner.takeAction(mRangeInTargetSet, mergedRegionList, mDummyRange);	
			createAndSetNewBufferRegion(Chrom.c0, 0);
		}
		
		// Finally, we perform an internal test to make sure each site in the chromosome in the map was represented in a region
		testConstructedRegions(mergedRegionList, snvMap, chrom);
	}

	// ========================================================================
	private void testConstructedRegions(ArrayList<T> regions, SNVMap snvMap, Chrom chrom) {

		int numSitesOnChrom = snvMap.getNumSitesOnChromosome(chrom);
		if (numSitesOnChrom <= 0) return;
		
		if (regions.isEmpty()) {
			CompareUtils.ensureTrue(false, "ERROR: Must have >= 1 region on Chrom: " + chrom.getCode());
		}
		int regionIndex = 0;  
		RegionRange<T> currentRegion = regions.get(regionIndex);
		
		for (int indexInMap = 0; indexInMap < numSitesOnChrom; indexInMap++) {
			int mapPosition = snvMap.getPosition(chrom, indexInMap);
			
			if (currentRegion.beforeRange(chrom, mapPosition)) {
				CompareUtils.ensureTrue(false, "ERROR: Map position must have a corresponding region, not only regions following it!");
			} else if (currentRegion.inRange(chrom, mapPosition)) {
				// Do nothing. This means the site exists within a region (i.e. the current region)
				// System.out.println("Pass: " + chrom + "\t" + mapPosition + "\t" + currentRegion.toString());
			} else if (currentRegion.afterRange(chrom, mapPosition)) {
				currentRegion = regions.get(++regionIndex);
				if (!currentRegion.inRange(chrom, mapPosition)) {
					CompareUtils.ensureTrue(false, "ERROR: At least, the next region must contain the current map position!");
				}
			} else {
				CompareUtils.ensureTrue(false, "ERROR: Impossible state!");
			}
		}
	}
	
	// ========================================================================
	public void clear() {
		mRangeInTargetSet.clear();
		mDummyRange.set(Chrom.c0, 0, 0, false, 1);		
	}

	// ========================================================================
	public RegionAndSiteWalker(NullaryClassFactory<T> theClassFactory) {
		mClassFactory = theClassFactory;
		createAndSetNewBufferRegion(Chrom.c0, 0);
		mRangeInTargetSet = new ArrayList<Boolean>();
	}
	
	// ========================================================================
	private void createAndSetNewBufferRegion(Chrom chrom, int position) {
		mDummyRange = mClassFactory.newInstance();		
		mDummyRange.set(chrom, position);				
	}

	// ========================================================================
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub

	}

}
