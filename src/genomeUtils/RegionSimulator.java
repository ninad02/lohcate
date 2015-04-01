package genomeUtils;

import genomeEnums.Chrom;

import java.util.ArrayList;

import lohcate.LOHcate;

import nutils.ArrayUtils;
import nutils.CompareUtils;
import nutils.NumberUtils;

public class RegionSimulator {
	
	// ========================================================================
	public static interface SampleInformation<E extends SiteInformation> {
		public int getIndexChromStart(Chrom chrom);
		public int getIndexChromEnd(Chrom chrom);
		public int getIndex(Chrom chrom, int position);
		public E getSiteAtIndex(Chrom chrom, int index);
	}

	// ========================================================================
	// MEMBER VARIABLES
	// ========================================================================
	
	// ========================================================================
	/** Generates a known region given the parameter values and adds it to the list of existing regions. */
	public static<E extends RegionRange<E>, T extends SiteInformation> 
		E generateRegion(int approxLength, E newRegion, ArrayList<E> existingRegions, SampleInformation<T> oneSampleInfo) {
		//ArrayList<CNRegion> regions = new ArrayList<CNRegion>();
		
		int numChromRetries = -1;
		while (true) {
			++numChromRetries;
			
			// First randomly select a chromosome
			Chrom chrom = GenotypeUtils.getRandomAutosomalChromosome();		
			LOHcate.LogOutput.println("Chrom Selected: " + chrom + "\tTry: " + numChromRetries);
			int indexChromFirstPosition = oneSampleInfo.getIndexChromStart(chrom);
			if (indexChromFirstPosition < 0) {
				// Data for this particular chromosome does not exist in this sample, so 
				// repeat the outer loop and pick a new chromosome.
				continue;
			}		

			// Now, we are guaranteed to have a chromosome that exists in the sample
			int indexChromLastPosition = oneSampleInfo.getIndexChromEnd(chrom);
			LOHcate.LogOutput.println("Index Bounds:\t" + chrom + "\t" + indexChromFirstPosition + "\t" + indexChromLastPosition);
			CompareUtils.ensureTrue(indexChromLastPosition >= indexChromFirstPosition, "ERROR: Index last position < Index first position in simulated region!");		
			int chromPositionLast = oneSampleInfo.getSiteAtIndex(chrom, indexChromLastPosition).getPosition(); 

			int numSameChromAttempts = 5;
			for (int i = 0; i < numSameChromAttempts; i++) {
				// Now get a random position	
				int randomIndex = NumberUtils.getRandomInteger(indexChromFirstPosition, indexChromLastPosition);
				int randomPosition = oneSampleInfo.getSiteAtIndex(chrom, randomIndex).getPosition();
				int lengthFromEnd = chromPositionLast - randomPosition;
				int compareFactor = (lengthFromEnd > approxLength) ? 1 : -1;  
				//SequenceLogger.outputPrintln(positionAtIndex + "\tChromEnd:\t" + chromPositionLast + "\tDiff:\t" + lengthFromEnd);

				// Now search for the position that will be the closest index
				int targetPosition = randomPosition + (compareFactor * approxLength);
				//targetPosition = Math.min(Chrom.MaxPositionOnLongestChrom, Math.max(0, targetPosition));
				
				int resultIndex = oneSampleInfo.getIndex(chrom, targetPosition); 				
				resultIndex = ArrayUtils.getInsertPoint(resultIndex);
				
				// We then make sure that resultIndex remains within the bounds of the chromosome
				resultIndex = Math.max(indexChromFirstPosition, Math.min(indexChromLastPosition, resultIndex));
				int resultPosition = oneSampleInfo.getSiteAtIndex(chrom, resultIndex).getPosition();

				//SequenceLogger.outputPrintln("End Pos: " + resultPosition); 			
				//SequenceLogger.outputPrintln("Length: " + (resultPosition - )); 			
				int numSitesInterrogated = Math.abs(resultIndex - randomIndex) + 1;
				if (compareFactor > 0) {
					newRegion.set(chrom, randomPosition, resultPosition, true, numSitesInterrogated);
				} else {
					newRegion.set(chrom, resultPosition, randomPosition, true, numSitesInterrogated);
				}
				//SequenceLogger.outputPrintln("Dummy: " + mDummyCNRegion);
				//SequenceLogger.outputPrintln("");

				boolean overlap = false;
				for (int j = 0; !overlap && j < existingRegions.size(); j++) {
					overlap = newRegion.overlapRange(existingRegions.get(j));
				}

				// If there's no overlap, then we add the region and exit
				if (!overlap) {									
					existingRegions.add(newRegion);
					return newRegion;
				}
			}
		}	
	}
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub

	}

}
