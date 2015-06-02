package shared;

import genomeEnums.Chrom;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.ArrayList;

import nutils.IOUtils;
import nutils.StringUtils;

public class CompareMAFs {

	public static final int MAF_Col_Chromosome =    4;
	public static final int MAF_Col_Position   =    5;
	public static final int MAF_Col_TumorBarcode   =  15;
	public static final int MAF_Col_NormalBarcode   = 16;
	public static final int MAF_Col_ValidationStatus = 24;
	public static final int MAF_Col_MutationStatus   = 25;
	public static final int MAF_Col_ValidationMethod = 28;
	
	
	public CompareMAFs() {
		// TODO Auto-generated constructor stub
	}
	
	public static void compareMAF(String testMAFFilename, String validationMAFFilename) {
		
		
		BufferedReader inTestMAF = IOUtils.getBufferedReader(testMAFFilename);
		BufferedWriter outTestMAF = IOUtils.getBufferedWriter(testMAFFilename + (new File(validationMAFFilename)).getName());
		
		int numMatch = 0;
		ArrayList<String> validationMAFLines = IOUtils.readAllLinesFromFile(validationMAFFilename);
		
		String line = null;		
		for (int lineCounter = 0; (line = IOUtils.getNextLineInBufferedReader(inTestMAF)) != null; lineCounter++) {
			String[] cols = line.split(StringUtils.TabStr);
			
			if (lineCounter == 0) {
				IOUtils.writeToBufferedWriter(outTestMAF, cols, StringUtils.TabStr, false, true);
				
			} else {
				int resultIndex = compareOneMAFLineWithMAFLines_Helper(line, validationMAFLines);
				if (resultIndex >= 0) {
					String lineValidation = validationMAFLines.get(resultIndex);
					
					for (int i = MAF_Col_ValidationStatus; i <= MAF_Col_ValidationMethod; i++) {
						cols[i] = StringUtils.extractNthColumnValue(lineValidation, i, StringUtils.TabStr);
					}
					
					numMatch++;
				} else {
					cols[MAF_Col_ValidationStatus] = "DeNovo";
				}
				
				IOUtils.writeToBufferedWriter(outTestMAF, cols, StringUtils.TabStr, false, true);
			}		  
		}
		
		IOUtils.closeBufferedWriter(outTestMAF);
		IOUtils.closeBufferedReader(inTestMAF);
		System.out.println(numMatch);
	}
	
	private static int compareOneMAFLineWithMAFLines_Helper(String lineMAF, ArrayList<String> linesMAF) {
		for (int i = 0; i < linesMAF.size(); i++) {
			if (compareTwoMAFLines_Helper(lineMAF, linesMAF.get(i))) {
				return i;
			}
		}
		return -1;
	}
	
	private static boolean compareTwoMAFLines_Helper(String line1, String line2) {
		String delim = StringUtils.TabStr;
		Chrom chrom1 = Chrom.getChrom(StringUtils.extractNthColumnValue(line1, MAF_Col_Chromosome, delim));
		Chrom chrom2 = Chrom.getChrom(StringUtils.extractNthColumnValue(line2, MAF_Col_Chromosome, delim));		
		if (chrom1 != chrom2) return false;
		
		// Next, compare the positions
		int pos1 = Integer.parseInt(StringUtils.extractNthColumnValue(line1, MAF_Col_Position, delim));
		int pos2 = Integer.parseInt(StringUtils.extractNthColumnValue(line2, MAF_Col_Position, delim));
		if (pos1 != pos2) return false;
		
		// Next, compare the Tumor sample IDs
		String tumorName1 = StringUtils.extractNthColumnValue(line1, MAF_Col_TumorBarcode, delim);
		String tumorName2 = StringUtils.extractNthColumnValue(line2, MAF_Col_TumorBarcode, delim);
		if ( ! tumorName1.equalsIgnoreCase(tumorName2)) return false;
		
		// Next, compare Normal sample IDs 
		String normalName1 = StringUtils.extractNthColumnValue(line1, MAF_Col_NormalBarcode, delim);
		String normalName2 = StringUtils.extractNthColumnValue(line2, MAF_Col_NormalBarcode, delim);
		if ( ! normalName1.equalsIgnoreCase(normalName2)) return false;
		
		return true;
	}

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		CompareMAFs.compareMAF(args[0], args[1]);

	}

}
