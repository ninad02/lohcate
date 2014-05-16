package shared;

import java.io.BufferedWriter;
import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.Hashtable;

import lohcateEnums.EventType;

import nutils.IOUtils;
import nutils.StringUtils;
import nutils.counter.BucketCounterEnum;

public class GeneEventSampleCompare {

	public static enum MatchType {
		MatchCall,
		MatchNoCall,
		MismatchCallDiffCall,
		MismatchCallNoCall,
		MismatchNoCallCall;
	}
	
	
	public static void doWork(String inFilename1, String inFilename2) {
		
		ArrayList<GeneEventPatient> gepList1 = getList(inFilename1);
		ArrayList<GeneEventPatient> gepList2 = getList(inFilename2);
		
		Hashtable<String, BucketCounterEnum<MatchType>> matchCounts = new Hashtable<>();
				
		for (GeneEventPatient gep2 : gepList2) {
			//System.out.println(gep2.mGeneName + "\t" + gep2.mPatient + "\t" + gep2.mEventType);
		}
		
		BufferedWriter out = IOUtils.getBufferedWriter(inFilename1 + ".vs." + (new File(inFilename2)).getName());
		BufferedWriter out2 = IOUtils.getBufferedWriter(inFilename1 + ".vs." + (new File(inFilename2)).getName() + ".GeneTallies.txt");
		StringBuilder sb = new StringBuilder(2048);
		
		int numMatch = 0;
		int numNonExist = 0;
		int numMismatch = 0;
		int numMatchNoCall = 0;
		int numMismatchOneNoCallOtherCall = 0;
		
		String delim = StringUtils.TabStr;

		// Write header
		sb.append("Gene (Nexus)")
		  .append(delim).append("Sample (Nexus)")
		  .append(delim).append("Event (Nexus)")
		  .append(delim).append("Gene (LOHcate)")
		  .append(delim).append("Sample (LOHcate)")
		  .append(delim).append("Event (LOHcate)")
		  .append(delim).append("1=Call exists in LOHcate, 0 otherwise")
		  .append(delim).append("1=Call exists in LOHcate and Mismatch")
		  .append(delim).append("1=Call exists in LOHcate and Match");
		IOUtils.writeToBufferedWriter(out, sb.toString(), true);

		
		for (GeneEventPatient gep1 : gepList1) {
			sb.setLength(0);
			sb.append(gep1.mGeneName).append(delim).append(gep1.mPatient).append(delim).append(gep1.mEventType);
			
			BucketCounterEnum<MatchType> counterForGene = matchCounts.get(gep1.mGeneName);
			if (counterForGene == null) {
				counterForGene = new BucketCounterEnum<>(MatchType.class);
				matchCounts.put(gep1.mGeneName, counterForGene);
			}
			
			int resultIndex = Collections.binarySearch(gepList2, gep1, TwoListComparator);
//			int resultIndex = 0;
//			for (resultIndex = 0; resultIndex < gepList2.size(); resultIndex++) {
//				if (TwoListComparator.compare(gep1, gepList2.get(resultIndex)) == 0) {
//					break;
//				}
//			}
//			
//			if (resultIndex >= gepList2.size()) {
			if (resultIndex < 0) {
				sb.append(delim).append(delim).append(delim);
				sb.append(delim).append("0").append(delim).append("0").append(delim).append("0");
				
				 if (gep1.mEventType != EventType.Ignored) {
					 numNonExist++;
					 counterForGene.increment(MatchType.MismatchCallNoCall);
				 } else {
					 numMatchNoCall++;
					 counterForGene.increment(MatchType.MatchNoCall);
				 }
				
			} else {
				GeneEventPatient gep2 = gepList2.get(resultIndex);
				sb.append(delim).append(gep2.mGeneName).append(delim).append(gep2.mPatient).append(delim).append(gep2.mEventType);
				
				if (gep1.mEventType == gep2.mEventType) {
					numMatch++;
					sb.append(delim).append("1").append(delim).append("0").append(delim).append("1");
					counterForGene.increment(MatchType.MatchCall);
					
				} else if (gep1.mEventType == EventType.Ignored) {
					// Do nothing
					sb.append(delim).append("1").append(delim).append("1").append(delim).append("0");
					numMismatchOneNoCallOtherCall++;
					counterForGene.increment(MatchType.MismatchNoCallCall);
					
				} else {
					numMismatch++;
					sb.append(delim).append("1").append(delim).append("1").append(delim).append("0");
					counterForGene.increment(MatchType.MismatchCallDiffCall);
				}
			}
			
			IOUtils.writeToBufferedWriter(out, sb.toString(), true);
		}
		IOUtils.closeBufferedWriter(out);
		
		// Print header
		sb.setLength(0);
		sb.append("Gene");
		for (MatchType mt : MatchType.values()) {
			sb.append(delim).append(mt);
		}
		IOUtils.writeToBufferedWriter(out2, sb.toString(), true);
		
		for (String geneName : matchCounts.keySet()) {
			BucketCounterEnum<MatchType> counts = matchCounts.get(geneName);
			
			sb.setLength(0);
			sb.append(geneName);
			for (MatchType mt : MatchType.values()) {
				sb.append(delim).append(counts.getCount(mt));
			}

			int matchTotal = counts.getCount(MatchType.MatchCall) + counts.getCount(MatchType.MatchNoCall);
			int mismatchTotal = counts.getCount(MatchType.MismatchCallDiffCall) + counts.getCount(MatchType.MismatchCallNoCall) + counts.getCount(MatchType.MismatchNoCallCall); 
			sb.append(delim).append(matchTotal);
			sb.append(delim).append(mismatchTotal);
			sb.append(delim).append(counts.getCount(MatchType.MatchCall) / (double) (counts.getCount(MatchType.MatchCall) + counts.getCount(MatchType.MismatchCallDiffCall) + counts.getCount(MatchType.MismatchCallNoCall)));
			
			IOUtils.writeToBufferedWriter(out2, sb.toString(), true);
		}
		IOUtils.closeBufferedWriter(out2);
		
		System.out.println("Num Matches:\t" + numMatch);
		System.out.println("Num MisMatches:\t" + numMismatch);
		System.out.println("Num NonExist:\t" + numNonExist);
		System.out.println("Fraction Matches (with Mismatches):\t" + (numMatch) / (double) (numMatch + numMismatch));
		System.out.println("Fraction Matches (with Mismatches and Nonexistants):\t" + (numMatch) / (double) (numMatch + numMismatch + numNonExist));
		
		System.out.println("Num Matches No Call:\t" + numMatchNoCall);
		System.out.println("Num Mismatches No Call:\t" + numMismatchOneNoCallOtherCall);
		System.out.println("Fraction Matches No Call:\t" + (numMatchNoCall / (double) (numMatchNoCall + numMismatchOneNoCallOtherCall)));
		
		int numPerfectMatches = numMatch + numMatchNoCall;
		System.out.println("Number Perfect Matches:\t" + numPerfectMatches);
		
		int numMismatchesTotal = numMismatch + numNonExist + numMismatchOneNoCallOtherCall;
		System.out.println("Number Mismatches:\t" + numMismatchesTotal);
		
		System.out.println("Fraction Match Total:\t" + (numPerfectMatches / (double) (numMismatchesTotal + numPerfectMatches)));
		
	
	}

	// ========================================================================
	private static Comparator<GeneEventPatient> TwoListComparator = new Comparator<GeneEventPatient>() {

		@Override
		public int compare(GeneEventPatient gep1, GeneEventPatient gep2) {
			int result = gep1.mGeneName.compareToIgnoreCase(gep2.mGeneName);
			if (result == 0) {
				result = gep1.mPatient.compareToIgnoreCase(gep2.mPatient);				
			}
			return result;
		}
		
	};
	
	// ========================================================================
	private static ArrayList<GeneEventPatient> getList(String inFilename) {
		//System.out.println("---------------------------------");
		ArrayList<String> allLines = IOUtils.readAllLinesFromFile(inFilename);		
		EventType theEvent;
		ArrayList<GeneEventPatient> list = new ArrayList<>(allLines.size());
		for (String line : allLines) {
			//System.out.println(line);
			String[] cols = line.split(StringUtils.TabPatternStr);
			
			
			theEvent = (cols[1].equalsIgnoreCase("NA") /*|| cols[1].equalsIgnoreCase("cnLOH")*/) ? EventType.Ignored : EventType.valueOf(cols[1]);
			GeneEventPatient gep = new GeneEventPatient(cols[0].trim(), theEvent, cols[2].trim());
			list.add(gep);			
		}
		Collections.sort(list);
		return list;
	}
	
	// ========================================================================	
	private static class GeneEventPatient implements Comparable<GeneEventPatient> {
		String mGeneName;
		EventType mEventType;
		String mPatient;
		
		private GeneEventPatient(String geneName, EventType et, String patient) {
			mGeneName = geneName;
			mEventType = et;
			mPatient = patient;
		}

		@Override
		public int compareTo(GeneEventPatient o) {
			int result = mGeneName.compareToIgnoreCase(o.mGeneName);
			if (result == 0) {
				result = mPatient.compareToIgnoreCase(o.mPatient);
				if (result == 0) {				
					result = mEventType.compareTo(o.mEventType);
				}
			}
			return result;
		}
		
		
	}
	
	// ========================================================================
	public static void main(String[] args) {
		doWork(args[0], args[1]);
	}
	
}
