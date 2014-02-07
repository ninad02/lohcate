package nutils;

import java.nio.file.Files;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;

public class SortUtils {

	public SortUtils() {
		// TODO Auto-generated constructor stub
				
	}

	// ========================================================================
	public static void sortFile(String inFilename, final int[] columnIndices, boolean[] columnsAreNumeric, boolean[] reverse, final String delim, int numHeaders) {
		ArrayList<String> allLines = sortLines(inFilename, columnIndices, columnsAreNumeric, reverse, delim, numHeaders);
		IOUtils.writeOutputFile(inFilename + ".sorted.txt", allLines);
	}
	
	// ========================================================================
	public static ArrayList<String> sortLines(String inFilename, final int[] columnIndices, boolean[] columnsAreNumeric, boolean[] reverse, final String delim, int numHeaders) {
		ArrayList<String> allLines = IOUtils.readAllLinesFromFile(inFilename, false, true);
		return sortLines(allLines, columnIndices, columnsAreNumeric, reverse, delim, numHeaders);
	}
	
	// ========================================================================
	/** Given a list of tab-delimited strings, this sorts the strings based on the order of given columns. 
	 *  Column indices begin at 0.
	 */ 
	public static ArrayList<String> sortLines(ArrayList<String> lines, final int[] columnIndices, boolean[] columnsAreNumeric, final boolean[] reverse, final String delim, int numHeaders) {
	
		// Ensure no negative numbers for headers
		numHeaders = Math.max(numHeaders, 0);  
		
		// Now check whether we have at least two non-header lines 
		if (lines.size() <= (numHeaders + 1)) return lines;
		
		// We extract the headers from the lines
		ArrayList<String> headers = new ArrayList<String>();
		for (int i = 0; i < numHeaders; i++) {			
			headers.add(lines.get(i));		
			lines.set(i, null);
		}
		ArrayUtils.removeNullElements(lines);
		
		// Next, check that we have columnIndices.  If it is null or of zero length, we just execute the normal sort
		if ((columnIndices == null) || (columnIndices.length == 0)) {
			Collections.sort(lines);
			lines.addAll(0, headers);
			return lines;
		}
		
		// Now ensure that we have a value for the numeric column indicator for every column index
		boolean columnsNumericNotMatching = false;
		final boolean[] columnsAreNumericFinal;  // need this for comparator below
		if ((columnsAreNumeric == null) || (columnsNumericNotMatching = (columnsAreNumeric.length < columnIndices.length))) {						
			boolean[] columnsAreNumericTemp = new boolean[columnIndices.length];
			Arrays.fill(columnsAreNumericTemp, false);
			if (columnsNumericNotMatching) {
				ArrayUtils.arrayCopy(columnsAreNumericTemp, columnsAreNumeric, columnsAreNumeric.length);
			}
			columnsAreNumericFinal = columnsAreNumericTemp;
		} else {
			columnsAreNumericFinal = columnsAreNumeric;
		}
				
		// Now define a custom comparator to perform the sorting  		
		Comparator<String> columnComparator = new Comparator<String>() {
			@Override
			public int compare(String s1, String s2) {
				for (int i = 0; i < columnIndices.length; i++) {
					int col = columnIndices[i];
					String s1ColStr = StringUtils.extractNthColumnValue(s1, col, delim);
					String s2ColStr = StringUtils.extractNthColumnValue(s2, col, delim);

					// Set the multiplier for reversing
					final int multiplier = reverse[i] ? -1 : 1; 
					
					if (columnsAreNumericFinal[i]) {
						double s1Val = Double.parseDouble(s1ColStr);
						double s2Val = Double.parseDouble(s2ColStr);
						int result = Double.compare(s1Val, s2Val);
						if (result != 0) return (multiplier * result);
					} else {
						int result = s1ColStr.compareTo(s2ColStr);
						if (result != 0) return (multiplier * result);
					}
				}
				return 0;
			}
			 
		};
		
		// Sort, and then add the header lines back
		Collections.sort(lines, columnComparator);
		lines.addAll(0, headers);
		
		return lines;
	}

	
	// ========================================================================
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		String inFilename = IOUtils.pathConcat(new String[] {"E:", "Research", "Data", "tcga-acc2", "regions", "gistic", "Regions.Matrix.JISTIC.txt"}); 				
		ArrayList<String> sortedLines = SortUtils.sortLines(inFilename, new int[] { 1, 2 }, new boolean[] { true, true }, new boolean[] { false, false }, StringUtils.TabStr, 1);
		StringUtils.removeDuplicates(sortedLines, null);
		IOUtils.writeOutputFile(inFilename + ".sorted.txt", sortedLines);		
		
	}

}
