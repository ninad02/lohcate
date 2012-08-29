package shared;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Random;

import com.sun.xml.internal.ws.wsdl.parser.MexEntityResolver;

import sun.security.jgss.spi.MechanismFactory;

/**
 * LOHcate --- A software tool for LOH calling and visualization in cancer genomes
 * D Wheeler & SG Reddy
 * Human Genome Sequencing Center, Baylor College of Medicine (Houston, TX)
 * 
 * @author Siddharth G. Reddy, Ninad Dewal
 *
 */
public class Utils {
	
	public static final String CommaStr = ",";
	public static final String TabStr        = "\t";
	public static final String TabPatternStr = "\\t";
	public static final String rsPrefix = "rs";
	public static final String rsNegative = "rs-1";
	public static final String SemicolonStr = ";";
	public static final String NAStr  = "N/A";
	
	public static final int NumAutosomes = 22;
	
	public static final FileExtensionAndDelimiter FileExtensionCSV = new FileExtensionAndDelimiter(".csv", CommaStr);
	public static final FileExtensionAndDelimiter FileExtensionTSV = new FileExtensionAndDelimiter(".tsv", TabStr);

	
	public static double log(float param) {
		if (param==0)
			return param;
		else
			return Math.log(param);
	}
	
	public static float normalize(float val, float min, float max, float new_min, float new_max) {
		return new_max - ( (new_max - new_min) * ((max - val) / (max - min)) );
	}
	
	public static String rmLeadingSpaces(String param) {
		while(param.charAt(0)==' ')
			param = param.substring(1, param.length());
		return param;
	}
	
	public static float min(float one, float two) {
		if (one < two)
			return one;
		return two;
	}
	public static float max(float one, float two) {
		if (one > two)
			return one;
		return two;
	}
	
	public static String gClean(String param) {
		String target = "";
		for (int i = 0; i<param.length(); i++) {
			if ((int)param.charAt(i)>=97 && (int)param.charAt(i)<=122)
				target += (char)((int)param.charAt(i)-32);
			else
				target += param.charAt(i);
		}
		return target;
	}
	
	public static int indexOf(ArrayList<Integer> arr, int val) {
		for (int i = 0; i<arr.size(); i++)
			if (arr.get(i)==val)
				return i;
		return -1;
	}
	
	public static int indexOf(ArrayList<String> arr, String elem) {
		for (int i = 0; i<arr.size(); i++)
			if (arr.get(i).equals(elem))
				return i;
		return -1;
	}
	
	// removes last n elements from the list
	public static<T> ArrayList<T> removeLastNElements(ArrayList<T> theList, int numLastElementsToRemove) {
		numLastElementsToRemove = Math.min(numLastElementsToRemove, theList.size());
		int listSize = theList.size();
		for (int i = 0; i < numLastElementsToRemove; i++) {
			theList.remove(listSize - 1);
			listSize--;
		}
		return theList;
	}
	
	/** Returns whether the value is in the range, with an exclusive lower bound and inclusive upper bound. */
	public static boolean inRangeLowerExclusive(float value, float boundLower, float boundUpper) {
		return ((boundLower < value) && (value <= boundUpper));
	}
	
	
	
	// removes any elements set to null in the list
	public static<T> ArrayList<T> removeNullElements(ArrayList<T> theList) {
		
		int indexNow = 0;
		int indexAhead = indexNow + 1;
		int numNullValues = 0;
		
		boolean remainingElementsNonNull = true;
		int size = theList.size();
		for ( ; (indexNow < size) && remainingElementsNonNull; indexNow++) {
			if (theList.get(indexNow) == null) {
				// Currently we have a null element.
				numNullValues++;
				
				// Now we move indexAhead forward until we 
				// reach a non-null value or the end of the list				
				for (indexAhead = Math.max(indexNow + 1, indexAhead); 
						(indexAhead < size) && (theList.get(indexAhead) == null); 
						indexAhead++) {}				
				
				// At this point, indexAhead has either crossed its bounds
				// or has found an element
				if (indexAhead >= size) {    
					// crossed bounds.  Means that there were no non-null elements past 
					// the current index.  All we do now is remove the current and 
					// remaining (null) elements from the list.
					int numLastElementsToRemove = size - indexNow;  					
					removeLastNElements(theList, numLastElementsToRemove);
					numNullValues += numLastElementsToRemove - 1;
					remainingElementsNonNull = false;
				} else {  
					// found a non-null element, so swap them
					theList.set(indexNow, theList.get(indexAhead));
					theList.set(indexAhead, null);
				}
				
			} 
		}
				
		return theList;		
	}
	
	
	private static void Test_removeNullElements() {
		
		ArrayList<Integer> list1 = new ArrayList<Integer>();
		ArrayList<Integer> list2 = new ArrayList<Integer>();
		ArrayList<Integer> listM = new ArrayList<Integer>();  // a master list
		
		ArrayList<Integer> list2IndicesToRemove = new ArrayList<Integer>();
		
		int numTrials = 10000;
		int numElementsInArray = 100;
		int numElementsToDelete = 30;
		Random randomGen = new Random();
		
		for (int i = 0; i < numElementsInArray; i++) {			
			listM.add(new Integer(i));
		}
		
		for (int trial = 0; trial < numTrials; trial++) {
			list1.clear();
			list2.clear();
			list1.addAll(listM);
			list2.addAll(listM);
			list2IndicesToRemove.clear();
			
			// Now select random elements to delete from list1 
			for (int j = 0; j < numElementsToDelete; ) {
				int indexToDelete = randomGen.nextInt(numElementsInArray);
				if (list1.get(indexToDelete) != null) {
					list1.set(indexToDelete, null);
					list2IndicesToRemove.add(indexToDelete);
					j++;						
				}
			}
			
			// Now we must delete from list2, though in reverse order
			Collections.sort(list2IndicesToRemove);
			Collections.reverse(list2IndicesToRemove);
			for (Integer theIndex : list2IndicesToRemove) {
				list2.remove(theIndex);
			}
			
			// Now perform the deletion on the list 1
			removeNullElements(list1);

			// Now compare
			boolean areListsEqual = list1.equals(list2);
			if (areListsEqual) {
				System.out.println(list1);
				System.out.println(list2);
				System.out.println("");
				//break;
			}
		}
	}
	
	/** Given a string that contains an rs#, this returns an rs#.  
	 *  If the string is "rs2342,abcd", this will return: "rs2342"
	 *  @return the number after rs and before a non-numeric character, if all valid
	 *  @return "rs-1", if the rs# is -1 (rs-1), as listed in some files
	 *  @return "rs", if there is no valid rs# after rs
	 *  @return null, if there is no "rs" substring in this line.  
	 */
	public static String extractRsNumberFromLine(String line) {
		int rsIndex = line.indexOf(rsPrefix);
		if (rsIndex >= 0) {
			int startIndex = rsIndex + rsPrefix.length(); 
			for (int i = startIndex; i < line.length(); i++) {
				char ch = line.charAt(i);
				if (!Character.isDigit(ch)) {
					// Two things can happen.  We are either at the end of the rs#, 
					// or we have rs-1, which we must be able to parse					
					if ((ch == '-') 
						&& (i == startIndex) 
						&& (i < line.length() - 1) 
						&& (line.charAt(i + 1) == '1')) {						
						return rsNegative;					
					} else {
						return line.substring(rsIndex, i);
					}
				}
			}
			return line.substring(rsIndex);
		} else {
			return null;
		}
	}
	
	public static void Test_extractRsNumberFromLine() {
		String[] toTest = new String[] {
				"rs-1",
				",rs-1",
				"rs-1,",
				";rs-1;",

				"rs-2",
				",rs-0",
				"rs-0,",
				";rs-0;",
				";rs;-1",
				";rs'123",

				"rs12345",
				"rs12345a",
				",rs12345",
				",rs12345,",
				
				",rs12$345,"
		};
		
		for (int i = 0; i < toTest.length; i++) {
			System.out.println(extractRsNumberFromLine(toTest[i]));
		}
	}
	
	public static void main(String[] args) {
		Test_removeNullElements();
		Test_extractRsNumberFromLine();
	}

	// ========================================================================
	// INNER CLASS
	// ========================================================================
	public static class FileExtensionAndDelimiter {
		public String mExtension;
		public String mDelimiter;
		
		public FileExtensionAndDelimiter(String extension, String delimiter) {
			mExtension = extension;
			mDelimiter = delimiter;
		}
	}
}

