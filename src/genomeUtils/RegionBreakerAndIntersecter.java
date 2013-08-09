package genomeUtils;

import static genomeUtils.RegionRange.RegionRangeOverlap.AfterWithOverlap;
import static genomeUtils.RegionRange.RegionRangeOverlap.BeforeWithOverlap;
import static genomeUtils.RegionRange.RegionRangeOverlap.ConsumedByAlignedLeft;
import static genomeUtils.RegionRange.RegionRangeOverlap.ConsumedByTotal;
import static genomeUtils.RegionRange.RegionRangeOverlap.Equals;
import static genomeUtils.RegionRange.RegionRangeOverlap.SubsumesAlignedLeft;
import static genomeUtils.RegionRange.RegionRangeOverlap.SubsumesTotal;
import genomeUtils.RegionRange.RegionRangeOverlap;

import java.util.ArrayList;

import lohcate.CopyNumberRegionRange;
import lohcateEnums.EventType;

import nutils.Cast;
import nutils.CompareUtils;

public class RegionBreakerAndIntersecter {

	//=========================================================================
	public static interface RegionIntersectTester<E> {		
		public boolean isValidRegion(E region);
		public boolean isInvalidRegion(E region);
		public void takeActionOnEqualRegions(E region);
	}

	//=========================================================================
	/** Given an "original" list of regions and an "added" list of regions, this
	 *  takes the union of the regions.  If there are overlaps, the regions
	 *  are broken into component pieces in which the overlapping sections 
	 *  are their own newfound regions, while the non-overlapping subsections
	 *  of the regions are created as newfound regions on their own.  For example:
	 *  
	 *  original:       xxxxxxxx         or             xxx
	 *  new:               xxx                       xxxxxxxx
	 *  resulting:      aaabbbcc                     aaabbbcc 
	 *  
	 *  The block of contiguous a's represent a region, the b's a region, and the
	 *  c's a region.  This is just one example.  There are more:
	 *  
	 *  original:          xxxxxxx        or         xxxxx   xxxxx   
	 *  added:          xxxxx   xxxxx                   xxxxxxx
	 *  resulting:      aaabbcccddeee                aaabbcccddeee
	 *  
	 *  original:         xxx xx xxx      or         xxxxxxxxxxxxxx
	 *  added:          xxxxxxxxxxxxxx                 xxx xx xxx
	 *  resulting:      aabbbcddefffgg               aabbbcddefffgg
	 *  
	 *  @param regionsTarget The list of original regions.  This parameter will
	 *  be modified by the function, so better to pass in a copy if the caller wants
	 *  to keep the original regions.
	 *  @param regionsAdded The list of regions to be added.
	 *  @return The new list of regions.  This, however, is just the original regions modified.
	 *  The function was designed this way (i.e. to modify the original regions) to increase
	 *  runtime efficiency -- namely, to avoid the memory calls to create a new list.  Updating
	 *  an already existing list is more efficient, and if the user wants the original to not
	 *  be modified, then he/she should pass in a copy.
	 */
	public static<E extends RegionRange> 
		ArrayList<E> takeUnionAndBreakDownIntersectingRegions(ArrayList<E> regionsTarget, ArrayList<E> regionsSource, RegionIntersectTester<E> regionTester, Class<E> theClass) {
		
		// Make a pseudo-shallow array copy so we don't alter the caller's master copy 
		regionsSource = new ArrayList<E>(regionsSource);  
		
		int indexTarget = 0;
		int indexSource = 0;
		RegionRange.RegionRangeOverlap overlapTypePredicted = null;
		
		while ((indexTarget < regionsTarget.size()) && (indexSource < regionsSource.size())) {
			//System.out.println(indexTarget + "\t" + regionsTarget.size() + "\t" + indexSource + "\t" + regionsSource.size());
			E regionTarget = regionsTarget.get(indexTarget);
			E regionSource = regionsSource.get(indexSource);

			// Check and make sure that the current target and source regions are valid
			if (regionTester.isInvalidRegion(regionTarget)) {					
				regionsTarget.remove(indexTarget);
				continue;
			} else if (regionTester.isInvalidRegion(regionSource)) {
				regionsSource.remove(indexSource);
				continue;
			}
			
			// Test for what type of overlap we have
			RegionRange.RegionRangeOverlap overlapType = regionTarget.testAndCharacterizeOverlap(regionSource);			
			
			// Some of the following cases reduce down to other cases.  The following comparison
			// checks to make sure that the reduction producres for cases were done properly.  Some
			// cases are endpoints and do not result in reductions, however.
			if ((overlapTypePredicted != null) && (overlapTypePredicted != overlapType)) {
				CompareUtils.throwErrorAndExit("ERROR: Predicted and determined region overlap types don't match!" 
						+ overlapTypePredicted + "\t" + overlapType);
			}
			
			//System.out.println(overlapType);
			switch(overlapType) {
			case Equals:
				regionTester.takeActionOnEqualRegions(regionTarget);				
				indexTarget++;
				indexSource++;
				overlapTypePredicted = null;
				break;
			case BeforeWithoutOverlap: case AdjacentBefore:
				// We leave the block in the target list as is, and we simply increment our target index
				indexTarget++;
				overlapTypePredicted = null;
				break;
			case AfterWithoutOverlap: case AdjacentAfter: {  // use braces to prevent declarations from spilling into other cases					
				regionsTarget.add(indexTarget, Cast.downcast(regionSource.getCopy(), theClass));
				indexTarget++;  // To keep with the correct target region
				indexSource++;  // move to the next region to add
				overlapTypePredicted = null;
				break;
			}
			case BeforeViaDiffChromosome: case AfterViaDiffChromosome: 
				overlapTypePredicted = null;
				CompareUtils.throwErrorAndExit("ERROR: Should be on different chromosomes!");					
				break;
			case SubsumesTotal: case SubsumesAlignedRight: case BeforeWithOverlap: {					
				breakdownThreeCasesHelper(regionsTarget, regionsSource, indexTarget, indexSource, false, true, theClass);
				indexTarget++;  // increment the target array index				

				// The following procedure reduces these cases to: 
				//    SubsumesAlignedLeft, if SubsumesTotal was this case
				//    Equals, if SubsumesAlignedRight was this case
				overlapTypePredicted = (overlapType == SubsumesTotal) ? SubsumesAlignedLeft : 
					((overlapType == BeforeWithOverlap) ? ConsumedByAlignedLeft : Equals);
				break;
			}
			case SubsumesAlignedLeft: {
				// The following procedure reduces the case to Equals, and whatever follows
				breakdownSubsumesConsumedAlignedLeftHelper(regionsTarget, regionsSource, indexTarget, indexSource, false, theClass);
				overlapTypePredicted = RegionRangeOverlap.Equals;
				break;
			}
			case ConsumedByTotal: case ConsumedByAlignedRight: case AfterWithOverlap: {
				// The following procedure reduces these cases to:
				//     ConsumedByAlignedLeft, if ConsumedByTotal was this case
				//     Equals, if ConsumedByAlignedRight was this case
				breakdownThreeCasesHelper(regionsSource, regionsTarget, indexSource, indexTarget, true, false, theClass);
				indexTarget++;  // increment the target array index
				
				overlapTypePredicted = (overlapType == ConsumedByTotal) ? ConsumedByAlignedLeft :
					((overlapType == AfterWithOverlap) ? SubsumesAlignedLeft : Equals);
				break;
			}
			case ConsumedByAlignedLeft: {
				// The following procedure reduces the case to Equals, and whatever follows
				breakdownSubsumesConsumedAlignedLeftHelper(regionsSource, regionsTarget, indexSource, indexTarget, true, theClass);
				overlapTypePredicted = RegionRangeOverlap.Equals;
				break;					
			}		
			default:
				CompareUtils.throwErrorAndExit("ERROR: Invalid option!");
			}				
		}
				

		// Add any remaining regions.  We know that either of the following two loops will execute, but not both
		// since the previous loop was broken by failure of one of the loop conditions.
		// 
		// We do not need to add any regions to the target array, as either the array elements
		// were already traversed, or they already exist in the original array.  However, we need to remove
		// elements that may not match the cluster type desired
		while (indexTarget < regionsTarget.size()) {
			if (regionTester.isInvalidRegion(regionsTarget.get(indexTarget))) {
				regionsTarget.remove(indexTarget);
			} else {
				indexTarget++;
			}
		}
		
		// We only need to add elements (actually, their copies) if more still exist in the source array. 
		for (; indexSource < regionsSource.size(); indexSource++) {
			E regionSource = regionsSource.get(indexSource);
			if (regionTester.isValidRegion(regionSource)) {
				regionsTarget.add(Cast.downcast(regionSource.getCopy(), theClass));
			}
		}
		
		return regionsTarget;
	}

	//=========================================================================
	private static<E extends RegionRange> void breakdownThreeCasesHelper(
			ArrayList<E> regionsList1, 
			ArrayList<E> regionsList2,
			int indexList1, int indexList2, 
			boolean makeCopyOfRegionFromList1,
			boolean insertNewRegionIntoList1,
			Class<E> theClass) {
		
		E region1 = regionsList1.get(indexList1);
		E region2 = regionsList2.get(indexList2);
		
		E cnrrBefore = Cast.downcast(region1.getCopy(), theClass);  // make a copy to keep the stats
		cnrrBefore.setRangeEnd(region2.getRangeStart() - 1);  // make left-adjacent to the toAdd region
		if (insertNewRegionIntoList1) {
			regionsList1.add(indexList1, cnrrBefore);  // add the preceding block to the array at this point
			indexList1++;
		} else {
			regionsList2.add(indexList2, cnrrBefore);  // add the preceding block to the array at this point
			indexList2++;
		}
		
		E cnrrCurrent = (makeCopyOfRegionFromList1) ? 
				Cast.downcast(region1.getCopy(), theClass) : 	// make a copy to keep the stats, and to modify the region without modifying the paramter version 
				region1;                                        // don't make a copy, but allow the region itself to be modified later
	
		cnrrCurrent.setRangeStart(region2.getRangeStart());  // make left-aligned to the toAdd region
		regionsList1.set(indexList1, cnrrCurrent);
	}


	//=========================================================================
	/** A helper function used by the SubsumesAlignedLeft and ConsumedByAlignedLeft cases. */
	private static<E extends RegionRange> void breakdownSubsumesConsumedAlignedLeftHelper(
			ArrayList<E> regionsList1, 
			ArrayList<E> regionsList2,
			int indexList1, int indexList2, 
			boolean makeCopyOfRegionFromList1,
			Class<E> theClass) {
		
		E region1 = regionsList1.get(indexList1);
		E region2 = regionsList2.get(indexList2);
		
		E cnrrAfter = Cast.downcast(region1.getCopy(), theClass);  // make a copy to keep the stats
		cnrrAfter.setRangeStart(region2.getRangeEnd() + 1);   // make right-adjacent to the toAdd region					
		regionsList1.add(indexList1 + 1, cnrrAfter);          // Add the next region
		
		E cnrrCurrent = (makeCopyOfRegionFromList1) ? 
				Cast.downcast(region1.getCopy(), theClass) : 		// make a copy to keep the stats, and to modify the region without modifying the paramter version 
				region1;                    // don't make a copy, but allow the region itself to be modified later
						
		cnrrCurrent.setRangeEnd(region2.getRangeEnd());  // make right-aligned to the toAdd region
		regionsList1.set(indexList1, cnrrCurrent);
	}


	//=========================================================================
	public static ArrayList<CopyNumberRegionRange>
		takeUnionAndBreakDownIntersectingRegionsOld(ArrayList<CopyNumberRegionRange> regionsTarget, ArrayList<CopyNumberRegionRange> regionsSource, final EventType clusterType) {		
		
		// Make a pseudo-shallow array copy so we don't alter the caller's master copy 
		regionsSource = new ArrayList<CopyNumberRegionRange>(regionsSource);  
		
		int indexTarget = 0;
		int indexSource = 0;
		RegionRange.RegionRangeOverlap overlapTypePredicted = null;
		
		while ((indexTarget < regionsTarget.size()) && (indexSource < regionsSource.size())) {
			//System.out.println(indexTarget + "\t" + regionsTarget.size() + "\t" + indexSource + "\t" + regionsSource.size());
			CopyNumberRegionRange regionTarget = regionsTarget.get(indexTarget);
			CopyNumberRegionRange regionSource = regionsSource.get(indexSource);
	
			// Check and make sure that we are using the correct cluster type
			if (regionTarget.mCopyNumberEventType != clusterType) {
				regionsTarget.remove(indexTarget);
				//indexTarget++;
				continue;
			} else if (regionSource.mCopyNumberEventType != clusterType) {
				regionsSource.remove(indexSource);
				//indexSource++;
				continue;
			}
			
			// Test for what type of overlap we have
			RegionRange.RegionRangeOverlap overlapType = regionTarget.testAndCharacterizeOverlap(regionSource);			
			
			// Some of the following cases reduce down to other cases.  The following comparison
			// checks to make sure that the reduction producres for cases were done properly.  Some
			// cases are endpoints and do not result in reductions, however.
			if ((overlapTypePredicted != null) && (overlapTypePredicted != overlapType)) {
				CompareUtils.throwErrorAndExit("ERROR: Predicted and determined region overlap types don't match!" 
						+ overlapTypePredicted + "\t" + overlapType);
			}
			
			//System.out.println(overlapType);
			switch(overlapType) {
			case Equals: 
				regionTarget.mRecurrenceScore += 1.0; //regionSource.mRecurrenceScore; 
				indexTarget++;
				indexSource++;
				overlapTypePredicted = null;
				break;
			case BeforeWithoutOverlap: case AdjacentBefore:
				// We leave the block in the target list as is, and we simply increment our target index
				indexTarget++;
				overlapTypePredicted = null;
				break;
			case AfterWithoutOverlap: case AdjacentAfter: {  // use braces to prevent declarations from spilling into other cases					
				regionsTarget.add(indexTarget, regionSource.getCopy());
				indexTarget++;  // To keep with the correct target region
				indexSource++;  // move to the next region to add
				overlapTypePredicted = null;
				break;
			}
			case BeforeViaDiffChromosome: case AfterViaDiffChromosome: 
				overlapTypePredicted = null;
				CompareUtils.throwErrorAndExit("ERROR: Should be on different chromosomes!");					
				break;
			case SubsumesTotal: case SubsumesAlignedRight: case BeforeWithOverlap: {					
				breakdownThreeCasesHelper(regionsTarget, regionsSource, indexTarget, indexSource, false, true, CopyNumberRegionRange.class);
				indexTarget++;  // increment the target array index				
	
				// The following procedure reduces these cases to: 
				//    SubsumesAlignedLeft, if SubsumesTotal was this case
				//    Equals, if SubsumesAlignedRight was this case
				overlapTypePredicted = (overlapType == SubsumesTotal) ? SubsumesAlignedLeft : 
					((overlapType == BeforeWithOverlap) ? ConsumedByAlignedLeft : Equals);
				break;
			}
			case SubsumesAlignedLeft: {
				// The following procedure reduces the case to Equals, and whatever follows
				breakdownSubsumesConsumedAlignedLeftHelper(regionsTarget, regionsSource, indexTarget, indexSource, false, CopyNumberRegionRange.class);
				overlapTypePredicted = RegionRangeOverlap.Equals;
				break;
			}
			case ConsumedByTotal: case ConsumedByAlignedRight: case AfterWithOverlap: {
				// The following procedure reduces these cases to:
				//     ConsumedByAlignedLeft, if ConsumedByTotal was this case
				//     Equals, if ConsumedByAlignedRight was this case
				breakdownThreeCasesHelper(regionsSource, regionsTarget, indexSource, indexTarget, true, false, CopyNumberRegionRange.class);
				indexTarget++;  // increment the target array index
				
				overlapTypePredicted = (overlapType == ConsumedByTotal) ? ConsumedByAlignedLeft :
					((overlapType == AfterWithOverlap) ? SubsumesAlignedLeft : Equals);
				break;
			}
			case ConsumedByAlignedLeft: {
				// The following procedure reduces the case to Equals, and whatever follows
				breakdownSubsumesConsumedAlignedLeftHelper(regionsSource, regionsTarget, indexSource, indexTarget, true, CopyNumberRegionRange.class);
				overlapTypePredicted = RegionRangeOverlap.Equals;
				break;					
			}		
			default:
				CompareUtils.throwErrorAndExit("ERROR: Invalid option!");
			}				
		}
				
	
		// Add any remaining regions.  We know that either of the following two loops will execute, but not both
		// since the previous loop was broken by failure of one of the loop conditions.
		// 
		// We do not need to add any regions to the target array, as either the array elements
		// were already traversed, or they already exist in the original array.  However, we need to remove
		// elements that may not match the cluster type desired
		while (indexTarget < regionsTarget.size()) {
			if (regionsTarget.get(indexTarget).mCopyNumberEventType != clusterType) {
				regionsTarget.remove(indexTarget);
			} else {
				indexTarget++;
			}
		}
		
		// We only need to add elements (actually, their copies) if more still exist in the source array. 
		for (; indexSource < regionsSource.size(); indexSource++) {
			CopyNumberRegionRange regionSource = regionsSource.get(indexSource);
			if (regionSource.mCopyNumberEventType == clusterType) {
				regionsTarget.add(regionSource.getCopy());
			}
		}
		
		return regionsTarget;
	}

	//=========================================================================
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub

	}
}
