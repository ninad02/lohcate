package shared;

import java.io.BufferedWriter;
import java.util.ArrayList;
import java.util.BitSet;

import com.carrotsearch.hppc.FloatArrayList;

import cern.colt.list.DoubleArrayList;

import nutils.IOUtils;
import nutils.NumberUtils;
import nutils.PrimitiveWrapper;
import nutils.StringUtils;
import nutils.counter.DynamicRoundedDoubleCounter;

public class GeneMutationPermutation {

	public static final boolean doByLength = true;
	public static final boolean CollapseGenesIntoOne = true;
	
	public static void doWork(String inFileTotal, String inFileTarget) {
		ArrayList<GeneMutationLengthAndCount> totalList = readGeneMutationLengthAndCountsFromFile(inFileTotal);
		ArrayList<GeneMutationLengthAndCount> targetList = readGeneMutationLengthAndCountsFromFile(inFileTarget);
		
		PrimitiveWrapper.WInteger totalLengthTarget   = new PrimitiveWrapper.WInteger();
		PrimitiveWrapper.WInteger totalMutCountTarget = new PrimitiveWrapper.WInteger();
		double mutFreqTarget = calcMutFrequency(targetList, totalLengthTarget, totalMutCountTarget);
		System.out.printf("Mut Freq Target: %g\n", mutFreqTarget);
		System.out.printf("Mut Count Target: %d\n", totalMutCountTarget.mInt);
		System.out.printf("Mut Length Target: %d\n", totalLengthTarget.mInt);
		System.out.printf("Gene Count: %d\n", targetList.size());
		
		//double mutFreqTotal  = calcMutFrequency(totalList);		
		
		int numIter = 1_000_000;
		DoubleArrayList mutFreqList = new DoubleArrayList(numIter);
		ArrayList<GeneMutationLengthAndCount> listForTrial = new ArrayList<GeneMutationLengthAndCount>();
		PrimitiveWrapper.WInteger totalLengthForTrial   = new PrimitiveWrapper.WInteger();
		PrimitiveWrapper.WInteger totalMutCountForTrial = new PrimitiveWrapper.WInteger();
		
		for (int iter = 0; iter < numIter; iter++) {			
			int criterion = doByLength ? totalLengthTarget.mInt : targetList.size();
			if (iter % 100000 == 0) { System.out.printf("Iter: %d\t%d\n", iter, criterion); }
			performRandomization(totalList, criterion, listForTrial);
			double mutFreq = calcMutFrequency(listForTrial, totalLengthForTrial, totalMutCountForTrial);
			mutFreqList.add(mutFreq);
		}
		
		// Now sort the mutation frequencies
		mutFreqList.sort();
		
		//
		BufferedWriter out = IOUtils.getBufferedWriter(inFileTarget + ".mutFreqs.txt");
		for (int i = 0; i < mutFreqList.size(); i++) {
			IOUtils.writeToBufferedWriter(out, "" + mutFreqList.get(i), true);			
		}
		IOUtils.closeBufferedWriter(out);
		
		// Now go from end to the beginning and find the index whose element is less than the target mutation frequency
		int theIndex = mutFreqList.size() - 1; 
		for (; theIndex >= 0; theIndex--) {
			if (mutFreqList.get(theIndex) < mutFreqTarget) {
				break;
			}
		}
				
		double pVal = 1;
		if (theIndex >= 0) {
			int numElem = mutFreqList.size() - theIndex;
			pVal = (double) numElem / (double) mutFreqList.size();
		}
		System.out.printf("P-Value: %g\n", pVal);
		
		DynamicRoundedDoubleCounter roundedCounter = new DynamicRoundedDoubleCounter(2);
		for (int i = 0; i < mutFreqList.size(); i++) {
			roundedCounter.register(mutFreqList.get(i));
		}
		
		com.carrotsearch.hppc.DoubleArrayList mutFreqUnique = new com.carrotsearch.hppc.DoubleArrayList();
		FloatArrayList mutFreqUniqueProportions = new FloatArrayList();
		roundedCounter.getValuesWithProportions(mutFreqUnique, mutFreqUniqueProportions);
		for (int i = 0; i < mutFreqUnique.size(); i++) {
			System.out.printf("%g\t%g\n", mutFreqUnique.get(i), mutFreqUniqueProportions.get(i));
		}
	}

	// ========================================================================
	private static BitSet ElementDone = new BitSet(65536); 
	public static synchronized void performRandomization(ArrayList<GeneMutationLengthAndCount> theList, final int totalRegionLength, ArrayList<GeneMutationLengthAndCount> returnList) {
		ElementDone.clear();
		returnList.clear();
		
		int lastIndex = theList.size() - 1;
		for (int lengthSoFar = 0; lengthSoFar < totalRegionLength; ) {
			int randomIndex = NumberUtils.getRandomInteger(0, lastIndex);
			if (!ElementDone.get(randomIndex)) {				
				ElementDone.set(randomIndex);
				GeneMutationLengthAndCount glmc = theList.get(randomIndex);
				returnList.add(glmc);
				lengthSoFar += (doByLength ? glmc.mLength : 1);
			}			
		}		
	}
	
	// ========================================================================
	public static double calcMutFrequency(ArrayList<GeneMutationLengthAndCount> gmlcList, PrimitiveWrapper.WInteger totalLength, PrimitiveWrapper.WInteger totalMutCount) {
		totalMutCount.mInt = 0;
		totalLength.mInt = 0;
		
		for (GeneMutationLengthAndCount gmlc : gmlcList) {
			totalMutCount.mInt += gmlc.mCount;
			totalLength.mInt   += gmlc.mLength;
		}
		
		return (double) totalMutCount.mInt / (double) totalLength.mInt;
	}
	
	// ========================================================================
	public static ArrayList<GeneMutationLengthAndCount> readGeneMutationLengthAndCountsFromFile(String inFilename) {
		ArrayList<String> allLines = IOUtils.readAllLinesFromFile(inFilename);
		
		ArrayList<GeneMutationLengthAndCount> gmlcList = new ArrayList<GeneMutationLengthAndCount>(allLines.size() * 5);
		for (String line : allLines) {
			String[] cols = line.split(StringUtils.TabPatternStr);
			//System.out.println(line);
			//if (cols.length < 4) continue;
			
			int theCount = Integer.parseInt(cols[0]);
			if (CollapseGenesIntoOne) {
				GeneMutationLengthAndCount gmlc = new GeneMutationLengthAndCount(cols[1], cols[2], Integer.parseInt(cols[3]), theCount);
				gmlcList.add(gmlc);
			} else {
				for (int i = 0; i < theCount; i++) {
					GeneMutationLengthAndCount gmlc = new GeneMutationLengthAndCount(cols[1], cols[2], Integer.parseInt(cols[3]), 1);
					gmlcList.add(gmlc);
				}		
			}
		}
		return gmlcList;
	}
	
	// ========================================================================
	private static class GeneMutationLengthAndCount {
		public String mGeneName;
		public String mTranscriptID;
		public int mCount;
		public int mLength;
		
		public GeneMutationLengthAndCount(String geneName, String transcriptID, int length, int count) {
			mCount = count;
			mLength = length;
			mGeneName = geneName;
			mTranscriptID = transcriptID;
		}
	}

	// ========================================================================
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		doWork(args[0], args[1]);
	}

}
