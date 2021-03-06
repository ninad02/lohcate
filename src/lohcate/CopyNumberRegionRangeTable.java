package lohcate;

import genomeEnums.Chrom;
import genomeUtils.RegionRange;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;

import lohcateEnums.EventType;

import nutils.IOUtils;
import nutils.MapUtils;
import nutils.StringUtils;

public class CopyNumberRegionRangeTable<E extends CopyNumberRegionRangeLOHcate> {
	
	HashMap<String, ArrayList<E>> mTable; 
	
	// ========================================================================
	public CopyNumberRegionRangeTable() {
		mTable = new HashMap<>();
	}

	// ========================================================================
	public static ArrayList<RegionRange<?>> readListOfRegions(String inFilename) {
		ArrayList<String> allLines = IOUtils.readAllLinesFromFile(inFilename);
		ArrayList<RegionRange<?>> allRegions = new ArrayList<>(allLines.size());
		
		for (String line : allLines) {
			String[] cols = line.split(StringUtils.TabPatternStr);
			Chrom chrom = Chrom.getChrom(cols[0]);
			int rangeStart = Integer.parseInt(cols[1]);
			int rangeEnd   = Integer.parseInt(cols[2]);	 
			
			RegionRange.Default rr = new RegionRange.Default(chrom, rangeStart, rangeEnd);
			allRegions.add(rr);
		}
		return allRegions;
	}
	
	// ========================================================================
	public static CopyNumberRegionRangeTable<CopyNumberRegionRangeLOHcate> readTableFromFile(String inFilename, boolean hasHeader) {
		CopyNumberRegionRangeTable<CopyNumberRegionRangeLOHcate> table = new CopyNumberRegionRangeTable<CopyNumberRegionRangeLOHcate>();
		
		BufferedReader in = IOUtils.getBufferedReader(inFilename);
		
		int lineCounter = 0;
		for (String line = null; ((line = IOUtils.getNextLineInBufferedReader(in)) != null); ) {
			String[] cols = line.split(StringUtils.TabPatternStr);
			if (++lineCounter == 1 && hasHeader) continue;  // skip the header line
			
			String patientName = cols[0];
			Chrom chrom = Chrom.getChrom(cols[1]);
			int regionStart = Integer.parseInt(cols[2]);
			int regionEnd   = Integer.parseInt(cols[3]);			
			int numMarkers     = (cols.length > 4) ?   Integer.parseInt(cols[4]) : 0;
			double log2CopyNum = (cols.length > 5) ? Double.parseDouble(cols[5]) : 0;
			EventType event    = (cols.length > 6) ?  EventType.valueOf(cols[6]) : EventType.Ignored;
			
			CopyNumberRegionRangeLOHcate cnrr = new CopyNumberRegionRangeLOHcate(event, chrom, regionStart, regionEnd);		
			cnrr.mCopyNumber = log2CopyNum;
			cnrr.mRecurrenceScore = numMarkers;
			cnrr.set(chrom, regionStart, regionEnd, true, numMarkers);
			
			ArrayList<CopyNumberRegionRangeLOHcate> regionsForPatient = MapUtils.getOrCreateArrayList(patientName, table.mTable);
			regionsForPatient.add(cnrr);			
		}
		
		IOUtils.closeBufferedReader(in);
		return table;
	}

	// ========================================================================
	public static void filterTable(String inFilename1, ArrayList<RegionRange<?>> allRegions, boolean hasHeader1) {
		CopyNumberRegionRangeTable<CopyNumberRegionRangeLOHcate> table1 = readTableFromFile(inFilename1, hasHeader1);		

		String outFilename = inFilename1 + ".filtered.txt";
		BufferedWriter out = IOUtils.getBufferedWriter(outFilename);
		
		String delim = StringUtils.FileExtensionTSV.mDelimiter;
		StringBuilder sb = new StringBuilder(1024);
		for (Map.Entry<String, ArrayList<CopyNumberRegionRangeLOHcate>> keyVal : table1.mTable.entrySet()) {
			
			String patientName = keyVal.getKey();	
			ArrayList<CopyNumberRegionRangeLOHcate> regions = keyVal.getValue();
			
 			for (CopyNumberRegionRangeLOHcate cnrr : regions) {
				
				for (RegionRange<?> regionOther : allRegions) {
					
					if (cnrr.overlapRange(regionOther)) {
						
						sb.setLength(0);
						sb.append(patientName);
						addRegionInfoToStringBuilder(delim, cnrr, sb);						  
						IOUtils.writeToBufferedWriter(out, sb.toString(), true);
						IOUtils.flushBufferedWriter(out);
						break;
					}
				}
			}
		}
		
		
		IOUtils.closeBufferedWriter(out);
	}
	

	// ========================================================================
	private static void addRegionInfoToStringBuilder(String delim, CopyNumberRegionRangeLOHcate cnrr, StringBuilder sb) {
		sb.append(delim)
		  .append(cnrr.getChromosome().getCode())
		  .append(delim)
		  .append(cnrr.getRangeStart())
		  .append(delim)
		  .append(cnrr.getRangeEnd())
		  .append(delim)
		  .append(cnrr.getNumSitesInterrogated())
		  .append(delim)
		  .append(cnrr.mCopyNumber)
		  .append(delim)
		  .append(cnrr.mCopyNumberEventType);
	}
	
	// ========================================================================
	public static void compareTables(String inFilename1, String inFilename2, boolean hasHeader1, boolean hasHeader2) {
		CopyNumberRegionRangeTable<CopyNumberRegionRangeLOHcate> table1 = readTableFromFile(inFilename1, hasHeader1);
		CopyNumberRegionRangeTable<CopyNumberRegionRangeLOHcate> table2 = readTableFromFile(inFilename2, hasHeader2);
		
		String outFilename = inFilename1 + ".vs." + (new File(inFilename2)).getName() + ".txt";
		BufferedWriter out = IOUtils.getBufferedWriter(outFilename);
		
		String delim = StringUtils.FileExtensionTSV.mDelimiter;
		StringBuilder sb = new StringBuilder(1024);
		for (Map.Entry<String, ArrayList<CopyNumberRegionRangeLOHcate>> keyVal : table1.mTable.entrySet()) {
			String patientName = keyVal.getKey();
			ArrayList<CopyNumberRegionRangeLOHcate> regions = keyVal.getValue();
			
			ArrayList<CopyNumberRegionRangeLOHcate> resultRegions = table2.mTable.get(patientName);
			if (resultRegions == null) {
				for (CopyNumberRegionRangeLOHcate cnrr : regions) {
					sb.setLength(0);
					sb.append(patientName);
					addRegionInfoToStringBuilder(delim, cnrr, sb);
					sb.append(delim)
					  .append("0")
					  .append(delim)
					  .append("0")
					  .append(delim)
					  .append("0")
					  .append(delim)
					  .append("0")
					  ;
					IOUtils.writeToBufferedWriter(out, sb.toString(), true);
				}
			} else {				
				for (CopyNumberRegionRangeLOHcate cnrr : regions) {
					boolean oneOverlapFound = false;
					boolean hadGermline = false;
					
					for (CopyNumberRegionRangeLOHcate cnrrOther : resultRegions) {						
						if (cnrr.overlapRange(cnrrOther)) {
							if (cnrrOther.mCopyNumberEventType == EventType.Ignored) continue;
							
							if (cnrrOther.mCopyNumberEventType == EventType.HETGermline) {
								hadGermline = true;
								continue;
							}
							
							oneOverlapFound = true;
							EventType et1 = cnrr.mCopyNumberEventType;
							EventType et2 = cnrrOther.mCopyNumberEventType;
							
							boolean eventsMatch =
									((et1 == EventType.LOH) && (et2 == EventType.LOH || et2 == EventType.DELHom))
								 || (et1 == EventType.GainSomatic && et2 == EventType.GainSomatic); 
									
														
							sb.setLength(0);
							sb.append(patientName);
							addRegionInfoToStringBuilder(delim, cnrr, sb);
							sb.append(delim)
							  .append("1")
							  .append(delim)
							  .append("0")
							  .append(delim)
							  .append("1")
							  .append(delim)
							  .append(eventsMatch ? "1" : "0");
							addRegionInfoToStringBuilder(delim, cnrrOther, sb);
							IOUtils.writeToBufferedWriter(out, sb.toString(), true);
						}
					}
					
					if (!oneOverlapFound) {
						sb.setLength(0);
						sb.append(patientName);
						addRegionInfoToStringBuilder(delim, cnrr, sb);
						sb.append(delim)
						  .append("1")
						  .append(delim)
						  .append(hadGermline ? "1" : "0")
						  .append(delim)
						  .append("0")
						  .append(delim)
						  .append("0")
						  ;
						IOUtils.writeToBufferedWriter(out, sb.toString(), true);
					}
				}
			}
		}
		
		IOUtils.closeBufferedWriter(out);
	}
	
	// ========================================================================
	public static void main(String[] args) {
		boolean filterStep = false;
		if (filterStep) {
			ArrayList<RegionRange<?>> allRegions = readListOfRegions(args[1]);		
			System.out.println("Read regions...\t" + allRegions.size());
			filterTable(args[0], allRegions, false);
		} else {
			compareTables(args[0], args[1], true, true);
		}
	}
	
}
