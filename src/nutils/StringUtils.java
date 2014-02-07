package nutils;

import java.util.ArrayList;
import java.util.BitSet;
import java.util.Comparator;
import java.util.ListIterator;
import java.util.Random;

import com.sun.org.apache.xml.internal.utils.StringComparable;


public class StringUtils {

	public static final FileExtensionAndDelimiter FileExtensionCSV = new FileExtensionAndDelimiter(".csv", StringUtils.CommaStr);
	public static final FileExtensionAndDelimiter FileExtensionTSV = new FileExtensionAndDelimiter(".tsv", StringUtils.TabStr);
	public static final FileExtensionAndDelimiter FileTextTabDelim = new FileExtensionAndDelimiter(".txt", StringUtils.TabStr);
	
	public static final String TabStr        = "\t";
	public static final String CommaStr = ",";
	public static final String TabPatternStr = "\\t";
	public static final String SemicolonStr = ";";
	public static final String SpaceString = " ";
	public static final String ColonString = ":";
	public static final String DoubleQuoteStr = "\"";
	public static final String DotStr = ".";
	
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

	// ========================================================================
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub

	}
	
	// ========================================================================
	/** This removes any leading whitespace from the start of the string
	 * @return A new string with all leading whitespace removed
	 */
	public static String removeLeadingWhitespace(String s) {
		int indexSubstringStart = 0;
		for (; indexSubstringStart < s.length() && Character.isWhitespace(s.charAt(indexSubstringStart)); indexSubstringStart++) {}
		return s.substring(indexSubstringStart, s.length());
	}
	

	// ========================================================================
	/** Given a string builder, this will convert the characters in the string builder to an integer.  
	 *  The argument errorCode will be set to -1 if there is an error.  Else it will be set to 0. */
	public static int atoi(StringBuilder sb, int startIndex, PrimitiveWrapper.WInteger errorCode) {
		errorCode.mInt = 0;
		int sum = 0;
		int multiplyingFactor = 1;
		
		// Make sure that the buffer has at least one char in it
		if ((sb.length() == 0) || (startIndex >= sb.length())) {
			errorCode.mInt = -2;
			return 0;
		}
		
		for (int i = sb.length() - 1; i >= startIndex; i--) {
			char ch = sb.charAt(i);
		
			if ((i == 0) && (ch == '-')) {
				sum *= -1;
			} else if (Character.isDigit(ch)) {
				int digitValue = ch - '0';
				sum += (multiplyingFactor * digitValue);
			} else {
				errorCode.mInt = -1;
				return 0;
			}
			
			multiplyingFactor *= 10;
		}
		return sum;
	}

	// ========================================================================
	/** Given a string, this returns a string containing the unique characters in
	 *  that string, maintained in the same order as the characters appeared. */	
	private static final BitSet BitSetUnique = new BitSet(65536);  // number of unicode characters 
	public synchronized static String extractUniqueCharacters(String s) {
		StringBuilder sb = new StringBuilder(s.length());		
		
		BitSetUnique.clear();
		for (int i = 0; i < s.length(); i++) {
			char ch = s.charAt(i);
			int bitIndex = Character.getNumericValue(ch);			
			if (!BitSetUnique.get(bitIndex)) {				
				BitSetUnique.set(bitIndex);
				sb.append(ch);
			}
		}
		return sb.toString();		
	}

	// ========================================================================
	/** Returns the index of the first whitespace char in a string.  If no
	 *  whitespace char exists, -1 is returned.
	 */
	public static int getIndexOfFirstWhitespaceChar(String s) {
		for (int i = 0; i < s.length(); i++) {
			if (Character.isWhitespace(s.charAt(i))) {
				return i;
			}
		}
		return -1;
	}

	// ========================================================================
	/** Returns the index of the nth occurrence of the indicated target char in the string. 
	 *  If no such occurrence exists, -1 is returned.
	 */
	public static int getIndexOfNthChar(String s, int startIndex, char targetChar, int nthOccurrence) {
		
		// First, see if the nthOccurrence is 0.  If it is, return 0
		if (nthOccurrence <= 0) return 0;
		
		int sLength = s.length();
		startIndex = Math.min(Math.max(0, startIndex), sLength - 1);
		nthOccurrence = Math.min(nthOccurrence, sLength - startIndex);
		
		int numOccurrences = 0;
		for (int i = startIndex; i < sLength; i++) {
			numOccurrences += (s.charAt(i) == targetChar) ? 1 : 0;
			
			// Test if this is the desired occurrence.
			if (numOccurrences == nthOccurrence) {
				return i;
			}
		}
		
		return -1;
	}
	
	// ========================================================================
	/** This identifies and returns the first whitespace character in a string.  If
	 *  no whitespace characters are found, 0 is returned. */
	public static char identifyFirstWhitespaceChar(String s) {
		for (int i = 0; i < s.length(); i++) {
			char ch = s.charAt(i);
			if (Character.isWhitespace(ch)) {
				return ch;
			}
		}
		return 0;
	}

	// ========================================================================
	/** Returns the number of non-blank elements in a string array. */
	public static int numNonBlankElementsIn(String[] strArray) {
		int count = 0;
		for (int i = 0; i < strArray.length; i++) {
			count = (strArray[i].equals("")) ? count : count + 1;  
		}
		return count;
	}

	// ========================================================================
	/** Given an array of Strings, this tallies the total length of all the strings combined. */
	public static int tallyTotalLength(String[] sArray) {
		int sum = 0;
		for (int i = 0; i < sArray.length; i++) {
			sum += sArray[i].length();
		}
		return sum;
	}

	// ========================================================================
	private static final StringBuilder ConstructStringSB = new StringBuilder(2048);
	/** Given a string and the number of repetitions, this constructs a new string
	 *  consisting of the arg string repeated the specified number of times. 
	 */
	public synchronized static String constructRepetitiveString(String strToRepeat, int numRepetitions) {
		ConstructStringSB.setLength(0);
		for (int i = 0; i < numRepetitions; i++) {
			ConstructStringSB.append(strToRepeat);
		}
		return ConstructStringSB.toString();
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
				CompareUtils.throwErrorAndExit("ERROR: The string does not contain at least two columns.  String:\n" + line);
			}
			
			// Now move the first delimiter over to the appropriate column
			for (int i = 1; i < nthColumn; i++) {
				indexFirstDelimiter = line.indexOf(delim, indexFirstDelimiter + 1);
				if (indexFirstDelimiter < 0) {
					CompareUtils.throwErrorAndExit("ERROR: The string does not have " + (nthColumn + 1) + " colum(s).  String:\n" + line);
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

	private static void TestExtractNthColumn() {
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
		String delim = StringUtils.TabStr;
		
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
			
			String[] truthCols = line.split(StringUtils.TabPatternStr);
			for (int i = 0; i < truthCols.length; i++) {
				String colExtracted = extractNthColumnValue(line, i, StringUtils.TabStr);
				if (!truthCols[i].equals(colExtracted)) {
					CompareUtils.throwErrorAndExit("ERROR: Values don't match!\n" + line + "\n" + i + "\t[" + truthCols[i] + "]\t[" + colExtracted + "]");
				}
			}
		}
	}

	// ========================================================================
	/** Given a list of strings, this removes any duplicate entries.  Assumes that the input is sorted. 
	 * 
	 * @param lines The list of lines to be filtered for duplicates.  This list of lines is changed by the method.
	 * @param comparator Compares two elements of the list.  If set to null, just tests via .equals()
	 * @return This same list of lines, except with duplicates removed.    
	 */
	public static<T> ArrayList<T> removeDuplicates(ArrayList<T> lines, Comparator<? super T> comparator) {
		if (lines.size() <= 1) return lines;

		T elementPrev = lines.get(0);  // Initialize the previous element
		for (int i = 1; i < lines.size(); i++) {
			T elementCurr = lines.get(i);
			if ( (UtilsBasic.isNull(comparator)    && elementPrev.equals(elementCurr)) ||
				 (UtilsBasic.isNotNull(comparator) && (comparator.compare(elementPrev, elementCurr) == 0)) ) {
				lines.set(i, null);
			} else {
				elementPrev = elementCurr;
			}
		}
		
		ArrayUtils.removeNullElements(lines);
		return lines;
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
}
