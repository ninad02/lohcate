package shared;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.Random;

import lohcateEnums.ClusterType;
import lohcateEnums.Nuc;

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
	public static final String RangeDash = "-";
	public static final int NumAutosomes = 22;
	
	public static final FileExtensionAndDelimiter FileExtensionCSV = new FileExtensionAndDelimiter(".csv", CommaStr);
	public static final FileExtensionAndDelimiter FileExtensionTSV = new FileExtensionAndDelimiter(".tsv", TabStr);

	
	public static double log(float param) {
		return (param == 0) ? param : Math.log(param);
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
	
	// ========================================================================
	/** Given a target collection, this adds elements from the source collection at the end
	 *  of the target collection.  We create this function because it is ironically more 
	 *  efficient than the addAll() method in some collections, which stupidly allocates
	 *  memory to create an extra and needless temporary array in its implementation.
	 * @param listToWhichToAdd
	 * @param elementsToAdd
	 */
	// We create this function because it is ironically more efficient than the ArrayList.addAll() method,
	// which stupidly allocates memory to create an extra and needless temporary array in its implementaion.
	public static<T> void addAll(ArrayList<T> listToWhichToAdd, ArrayList<T> elementsToAdd) {
		listToWhichToAdd.ensureCapacity( listToWhichToAdd.size() + elementsToAdd.size() );  // Ensure capacity so no memory reallocation later		
		for (T element : elementsToAdd) {
			listToWhichToAdd.add(element);
		}
	}
	
	/** Given an array of ClusterTypes, this returns the counts of the ClusterTypes, 
	 *  according to ordering of the ClusterType enums. 
	 */
	public static int[] getClusterTypeCounts(ClusterType[] clusterTypeArray) {
		int[] counts = new int[ClusterType.values().length];
		Arrays.fill(counts, 0);
		for (ClusterType ct : clusterTypeArray) {
			counts[ct.ordinal()]++;
		}
		
		return counts;
	}
	
	/** Given a string of nucleotides, this returns the fraction of characters that are G or C. */	
	private static final char[] charBuffer = new char[65536];
	
	/** @return -1, if no valid nucleotides are present in the string; or if valid nucleotides are present,
	 *  the fraction of G/C nucleotides
	 */
	public static synchronized double calcFractionGC(String nucleotideString) {		
		int strLen = nucleotideString.length();
		int validLength = 0;
		int numGC = 0;
		nucleotideString.getChars(0, strLen, charBuffer, 0);  // more efficient than calling charAt() each character
		
		for (int i = 0; i < strLen; i++) {
			Nuc theNuc = Nuc.getNucUnsafe(charBuffer[i]);
			if (theNuc != null) {
				validLength++;
				if (theNuc.isGC()) {
					numGC++;
				}
			}
		}
		
		// Now divide
		if (validLength == 0) return -1;
		return (double) numGC / (double) validLength;
	}
	
	/** Given an arraylist of arraylists of a certain type, this adds numNewArraysToAdd new ArrayLists of
	 *  that type to the target list.
	 * @param targetList
	 * @param numNewArraysToAdd
	 */
	public static<T> void addNewEmptyArrayLists(ArrayList<ArrayList<T>> targetList, int numNewArraysToAdd) {
		for (int i = 0; i < numNewArraysToAdd; i++) {
			targetList.add(new ArrayList<T>());
		}
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
	
	/** Given a tab-delimited string, this returns the value (as a String of the nth column, where 0
	 *  represents the first column, 1 represents the second column, and so on.  
	 *  @return the value in the column, or null if the column does not exist.  A blank
	 *  string is returned if the column exists but contains no value (i.e. two consecutive
	 *  tab characters).
	 */
	public static String extractNthColumnValue(String line, int nthColumn, String delim) {		
				
		// Test whether we have a blank line or not
		if (line.isEmpty()) return "";
		
		// The location of the first delimiter
		int indexFirstDelimiter = line.indexOf(delim);
		
		// Test this special case.  If the 0th column is desired, the line may only contain
		// one column with no delimiter.  In addition, if the first column is blank (evidenced
		// by a delimiter character at the beginning of the line), we just return a blank string.
		if (nthColumn == 0) {
			
			// If we have no delimiter existing in the line, we set its index to the string length 
			indexFirstDelimiter = (indexFirstDelimiter < 0) ? line.length() : indexFirstDelimiter;
			
			// Now we return the first token, whether it exists or is blank 
			return ((indexFirstDelimiter == 0) ? "" : line.substring(0, indexFirstDelimiter));
			
		} else {
			
			// Test whether our delimiter exists
			if (indexFirstDelimiter < 0) {
				Utils.throwErrorAndExit("ERROR: The string does not contain at least two columns.  String:\n" + line);
			}
			
			// Now move the first delimiter over to the appropriate column
			for (int i = 1; i < nthColumn; i++) {
				indexFirstDelimiter = line.indexOf(delim, indexFirstDelimiter + 1);
				if (indexFirstDelimiter < 0) {
					Utils.throwErrorAndExit("ERROR: The string does not have " + (nthColumn + 1) + " colum(s).  String:\n" + line);
				}
			}
			
			// Now we have our delimiter placed correctly.  We now check for the second delimiter
			int indexSecondDelimiter = line.indexOf(delim, indexFirstDelimiter + 1);
			
			// If the second delimiter index doesn't exist, then it means we are at the end of the line
			indexSecondDelimiter = (indexSecondDelimiter < 0) ? line.length() : indexSecondDelimiter;
			
			// Now return the substring
			if (1 == indexSecondDelimiter - indexFirstDelimiter) {
				return "";
			} else {
				return line.substring(indexFirstDelimiter + 1, indexSecondDelimiter);
			}
		}
	}
	
	/** Returns whether the value is in the range, with an exclusive lower bound and inclusive upper bound. */
	public static boolean inRangeLowerExclusive(float value, float boundLower, float boundUpper) {
		return ((boundLower < value) && (value <= boundUpper));
	}
	
	/** Compares two long values.  -1 if first argument is smaller, 0 if equal, 1 if first argument is greater. */
	public static int compareLong(final long num1, final long num2) {
		return (num1 == num2) ? 0 : ((num1 < num2) ? -1 : 1);
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
	
	public static void TestExtractNthColumn() {
		String testString = "a b cd efg h i jkl m";		
		testString = testString.replaceAll(" ", "\t");
		System.out.println(testString);
		
		for (int i = 0; i < 8; i++) {
			System.out.println(extractNthColumnValue(testString, i, "\t"));
		}
		
		System.out.println(Character.toChars('a' + 2));
	}
	
	public static void TestExtractNthColumnRobust() {
		int numTrials = 10000000;
		StringBuilder sb = new StringBuilder(4096);
		Random randomGen = new Random();
		String delim = Utils.TabStr;
		
		for (int trial = 0; trial < numTrials; trial++) {
			if (trial % 1000 == 0) System.out.println("Trial: " + trial);
			sb.setLength(0);
			int numColumns = randomGen.nextInt(10) + 1;
			for (int col = 0; col < numColumns; col++) {
				int colLength = randomGen.nextInt(4);
				for (int i = 0; i < colLength; i++) {
					int charCode = 'a' + randomGen.nextInt(22);
					sb.append(Character.toChars(charCode));
				}
				if (col < numColumns - 1) { sb.append(delim); }
			}
			String line = sb.toString();
			//System.out.println(trial + "\t[" + line + "]");
			
			String[] truthCols = line.split(TabPatternStr);
			for (int i = 0; i < truthCols.length; i++) {
				String colExtracted = extractNthColumnValue(line, i, Utils.TabStr);
				if (!truthCols[i].equals(colExtracted)) {
					Utils.throwErrorAndExit("ERROR: Values don't match!\n" + line + "\n" + i + "\t[" + truthCols[i] + "]\t[" + colExtracted + "]");
				}
			}
		}
	}
	
	public static void main(String[] args) {
		//Test_removeNullElements();
		//Test_extractRsNumberFromLine();
		//TestExtractNthColumn();
		TestExtractNthColumnRobust();
	}
	
	// ========================================================================
	public static StringBuilder constructColumnDelimitedString(String[] values, String delimiter, StringBuilder sb, boolean clearStringBuilder) {
		if (clearStringBuilder) {
			sb.setLength(0);
		}
		
		for (int i = 0; i < values.length; i++) {
			if (i > 0) {
				sb.append(delimiter);
			}
			sb.append(values[i]);
		}
		return sb;
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

	public static void throwErrorAndExit(String errorString) {
		Exception e = new Exception(errorString);
		e.printStackTrace();
		System.exit(-1);			
	}

	/** This assumes that the list is sorted.  It tests whether the element already exists in
	 *  the list and inserts the element only if it doesn't already exist in the list. 
	 * @param targetList The list that is tested for the presence of the element in the list  
	 * @param element The element to be added if it doesn't already exist in the list
	 * @return true if the element already existed in the list, false otherwise
	 */
	public static<T> boolean checkInListAndInsertIfMissing(List <? extends Comparable<? super T>> list, T key) {
		int resultIndex = Collections.binarySearch(list, key);
		if (resultIndex < 0) {
			int insertionIndex = -(resultIndex + 1);
			List<T> listCasted = (List<T>) list;
			listCasted.add(insertionIndex, key);			
			return false;
		} else {
			return true;
		}		
	}
}

