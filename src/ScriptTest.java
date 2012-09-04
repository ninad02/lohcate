import java.util.ArrayList;
import java.util.Collections;
import java.util.Random;

import shared.Utils;

import lohcateEnums.Chrom;
import lohcateEnums.ClusterType;

/**
 * LOHcate --- A software tool for LOH calling and visualization in cancer genomes
 * D Wheeler, N Dewal, & SG Reddy
 * Human Genome Sequencing Center, Baylor College of Medicine (Houston, TX)
 * 
 * This is a class used to develop test functions for methods in its counterpart java file
 * 
 * @author Ninad Dewal
 *
 */
public class ScriptTest extends Script {
	
	public static void TestRegionOverlaps() {
		ArrayList<CopyNumberRegionRange> regionList = new ArrayList<CopyNumberRegionRange>();
		
		regionList.add(new CopyNumberRegionRange(ClusterType.LOH, Chrom.c1, 1, 1000));  //1
		regionList.add(new CopyNumberRegionRange(ClusterType.LOH, Chrom.c2, 1, 1));     //2
		regionList.add(new CopyNumberRegionRange(ClusterType.LOH, Chrom.c2, 1, 1));
		regionList.add(new CopyNumberRegionRange(ClusterType.LOH, Chrom.c2, 1, 10));
		regionList.add(new CopyNumberRegionRange(ClusterType.LOH, Chrom.c2, 4, 5));
		regionList.add(new CopyNumberRegionRange(ClusterType.LOH, Chrom.c2, 9, 10));
		regionList.add(new CopyNumberRegionRange(ClusterType.LOH, Chrom.c2, 5, 9));
		regionList.add(new CopyNumberRegionRange(ClusterType.LOH, Chrom.c2, 6, 8));
		
		for (int i = 0; i < regionList.size(); i++) {
			System.out.print(i);
			for (int j = 0; j < regionList.size(); j++) {
				System.out.print("\t");
				
				CopyNumberRegionRange cnrr1 = regionList.get(i);
				CopyNumberRegionRange cnrr2 = regionList.get(j);
				System.out.print(cnrr1.testAndCharacterizeOverlap(cnrr2));								
			}
			System.out.println("");
		}
	}
	
	public static void TestRegionCombining() {
		ArrayList< ArrayList<CopyNumberRegionRange> > regionsFromSamples = new ArrayList< ArrayList<CopyNumberRegionRange> >();
		int numSamples = 5;
		for (int i = 0; i < numSamples; i++) { regionsFromSamples.add(new ArrayList<CopyNumberRegionRange>()); }

		StringBuilder sb = new StringBuilder(4096);
		String resultStringPrev = null;
		ArrayList<String> outList = new ArrayList<String>();
		
		int numTrials = 500;
		for (int trial = 0; trial < numTrials; trial++) {
			System.out.println("Trial: " + trial);
			// Clear the contained objects
			for (ArrayList<CopyNumberRegionRange> regionsFromOneSample : regionsFromSamples) {
				regionsFromOneSample.clear();				
			}
			
			// Now fill
			regionsFromSamples.get(0).add(new CopyNumberRegionRange(ClusterType.LOH, Chrom.c2, 1, 10));

			regionsFromSamples.get(1).add(new CopyNumberRegionRange(ClusterType.LOH, Chrom.c2, 1, 1));

			regionsFromSamples.get(2).add(new CopyNumberRegionRange(ClusterType.LOH, Chrom.c2, 4, 5));
			regionsFromSamples.get(2).add(new CopyNumberRegionRange(ClusterType.LOH, Chrom.c2, 9, 10));

			regionsFromSamples.get(3).add(new CopyNumberRegionRange(ClusterType.LOH, Chrom.c2, 5, 9));

			regionsFromSamples.get(4).add(new CopyNumberRegionRange(ClusterType.Dup, Chrom.c2, 6, 8));

			
			// Do random shuffling
			Random randomGen = new Random();
			boolean shouldShuffle = true;
			if (shouldShuffle) {
				for (int j = 0; j < 2; j++) {
					int ind1 = randomGen.nextInt(numSamples);
					int ind2;
					do {
						ind2 = randomGen.nextInt(numSamples);
					} while (ind2 == ind1);

					// Now swap
					ArrayList<CopyNumberRegionRange> regionsTemp = regionsFromSamples.get(ind1);
					regionsFromSamples.set(ind1, regionsFromSamples.get(ind2));
					regionsFromSamples.set(ind2, regionsTemp);

				}
			}

			boolean goForward = false;
			int startIndex = goForward ? 0 : numSamples - 1;
			int increment =  goForward ? 1 : -1;

			for (int i = startIndex + increment; (i >= 0) && (i < numSamples); i += increment) {
				//System.out.println("");
				//System.out.println(i);
				Script.takeUnionAndBreakDownIntersectingRegions(regionsFromSamples.get(startIndex), regionsFromSamples.get(i), ClusterType.LOH);
			}

			// Now we print the regions
			outList.clear();			
			for (CopyNumberRegionRange cnrr : regionsFromSamples.get(startIndex)) {
				outList.add(cnrr.toString() + "\t" + cnrr.mCopyNumberClusterType + "\t" + cnrr.mRecurrenceScore);				
			}
			Collections.sort(outList);
			
			// Now construct the string
			sb.delete(0, sb.length());	
			for (String s : outList) {
				sb.append(s).append("\n");
			}
			String resultString = sb.toString();
			
			if (resultStringPrev == null) {
				resultStringPrev = resultString;
				System.out.println(resultStringPrev);
			} else {
				boolean match = resultString.equals(resultStringPrev);
				if (!match) {
					Utils.throwErrorAndExit("ERROR: Doesn't match!\n" + resultString + "\n" + resultStringPrev);
				}
			}
		}
	}

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		//TestRegionOverlaps();
		TestRegionCombining();

	}

}
