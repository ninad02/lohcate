package shared;

import java.io.*;
import java.util.*;

/**
 * LOHcate --- A software tool for LOH calling and visualization in cancer genomes
 * D Wheeler & N Dewal & SG Reddy
 * Human Genome Sequencing Center, Baylor College of Medicine (Houston, TX)
 * 
 * @author Ninad Dewal
 *
 */

public class IOUtils {

	private static final Random randomGen = new Random();
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub

	}

	// ========================================================================
	/** Given a directory full path, this ensures that the directory path is created. */
	public static void createDirectoryPath(String directoryFullPath, boolean lastElementIsFilename) {
		File dirPath = new File(directoryFullPath);
		if (lastElementIsFilename) {
			dirPath = dirPath.getParentFile();			
		}
		
		try {
			if (!dirPath.isDirectory()) {
				boolean result = dirPath.mkdirs();
				if (!result) {
					Utils.throwErrorAndExit("ERROR: Could not create directory: " + dirPath.getAbsolutePath());
				}
			}
		} catch (SecurityException e) {
			System.err.println("ERROR: Doesn't have necessary permissions to create directory: " + dirPath.getAbsolutePath());
			System.err.println("       User must create directory him/herself first!");
			e.printStackTrace();			
		}
	}
	
	
	// ========================================================================
	/** Given a filename string and a subdirectory to be created, this creates the subdirectory and returns the filename. */
	public static String createSubdirectory(String inFilename, String newDirectoryName, boolean createDirectoryOnDisk) {
		String parentPath = (new File(inFilename)).getParent();
		parentPath = (parentPath == null) ? "" : parentPath + File.separatorChar;
		newDirectoryName = (newDirectoryName.charAt(newDirectoryName.length() - 1) == File.separatorChar) ? newDirectoryName : newDirectoryName + File.separatorChar;
		String newDirPath = parentPath + newDirectoryName + File.separatorChar;
		if (createDirectoryOnDisk) {
			File newDir = new File(newDirPath);
			if (!newDir.exists()) {			
				try {
					newDir.mkdirs();
				} catch (SecurityException e) {
					e.printStackTrace();
					System.exit(-1);
				}
			}
		}
		return newDirPath;
	}

	// ========================================================================
	/** Given a parent pathname and a filename, this returns the full path. */
	public static String pathConcat(String parentPath, String filename) {
		return parentPath + File.separatorChar + filename;
	}

	// ========================================================================
	/** Given an array of names in a path, this concatenates all strings using
	 *  the appropriate path separator character. */
	public static String pathConcat(String[] pathNameArray) {
		int totalLen = 0;
		for (int i = 0; i < pathNameArray.length; i++) {
			totalLen += pathNameArray[i].length();
		}
		
		StringBuilder sb = new StringBuilder(totalLen + pathNameArray.length + 32);  // 32 is extra buffer just in case		
		for (int i = 0; i < pathNameArray.length; i++) {
			sb.append(pathNameArray[i]);
			if (i < pathNameArray.length - 1) {
				sb.append(File.separatorChar);
			}
		}
		
		return sb.toString();
	}

	// ========================================================================
	/** Deletes a file and returns whether the operation was a successs. */
	public static boolean deleteFile(String filename) {
		return (new File(filename)).delete();
	}
	
	// ========================================================================
	/** Deletes all the files in a list.  Returns whether deleting all the files
	 *  was successful.  Returns false even if one file failed to delete. */
	public static boolean deleteAllFiles(Collection<String> filenames) {
		boolean successfullyDeleted = true;
		for (String tempFile : filenames) {
			successfullyDeleted = deleteFile(tempFile) && successfullyDeleted;
		}	
		return successfullyDeleted;
	}
	
	// ========================================================================
	/** Returns a PrintStream for the filename.  */
	public static PrintStream getPrintStream(String outFilename) {
		try { 
			return new PrintStream(outFilename);
		} catch (FileNotFoundException e) {
			e.printStackTrace();
			System.exit(-1);
		}
		return null;
	}
	
	// ========================================================================
	/** Closes the PrintStream. */
	public static void closePrintStream(PrintStream out) {
		out.close();
	}
	
	// ========================================================================
	/** Returns a BufferedReader object for the filename.  This function exists
	 *  to avoid the use of try-catch blocks, which may speed up performance. 
	 */
	public static BufferedReader getBufferedReader(String inFilename) {
		try { 
			return new BufferedReader(new FileReader(inFilename));
		} catch (IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}
		return null;
	}

	// ========================================================================
	/** Returns a BufferedWriter object for the filename.  This function exists
	 *  to avoid the use of try-catch blocks in the caller code, which may 
	 *  speed up performance.
	 */
	public static BufferedWriter getBufferedWriter(String outFilename) {
		return IOUtils.getBufferedWriter(outFilename, false);
	}

	// ========================================================================
	/** Returns a BufferedWriter object for the filename.  This function exists
	 *  to avoid the use of try-catch blocks in the caller code, which may 
	 *  speed up performance.
	 */
	public static BufferedWriter getBufferedWriter(String outFilename, boolean append) {
		try { 
			return new BufferedWriter(new FileWriter(outFilename, append));
		} catch (IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}
		return null;
	}

	// ========================================================================
	/** Given a BufferedWriter, this writes the next column.  It also writes the delimiter/separator
	 *  if desired, as well as a new line if desired.
	 */
	public static void writeToBufferedWriter(BufferedWriter out, String columnString, String separator, boolean shouldWriteSeparator, boolean writeNewLine) {
		try {
			if (shouldWriteSeparator) {
				out.write(separator);
			}
			out.write(columnString);
			if (writeNewLine) {
				out.newLine();
			}
		} catch (IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}
	}
	
	// ========================================================================
	/** Given a BufferedWriter, this writes the next line.  This function 
	 *  exists to avoid the use of try-catch blocks in the caller code.
	 */
	public static void writeToBufferedWriter(BufferedWriter out, String stringToWrite, boolean writeNewLine) {
		try { 
			out.write(stringToWrite);
			if (writeNewLine) {
				out.newLine();
			}
		} catch (IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}
	}
	
	// ========================================================================
	/** Writes a String array to the output. */
	public static void writeToBufferedWriter(BufferedWriter out, String[] arrayToWrite, String columnSeparator, boolean prefixWithColumnSeparator, boolean writeNewLine) {
		try { 
			for (int i = 0; i < arrayToWrite.length; i++) {
				if ((i > 0) || prefixWithColumnSeparator) {
					out.write(columnSeparator);
				}
				out.write(arrayToWrite[i]);
			}
			
			if (writeNewLine) {
				out.newLine();
			}
		} catch (IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}
		
	}
	
	// ========================================================================
	/** Writes a String array to the output, with each String element written on its unique line. */
	public static void writeToBufferedWriter(BufferedWriter out, String[] arrayToWrite) {
		try {
			for (int i = 0; i < arrayToWrite.length; i++) {
				out.write(arrayToWrite[i]);
				out.newLine();
			}
		} catch (IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}
	}
	

	// ========================================================================
	/** Given a BufferedWriter, this flushes the stream. */
	public static void flushBufferedWriter(BufferedWriter out) {
		try { 
			out.flush();
		} catch (IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}
	}

	// ========================================================================
	/** Given a BufferedWriter, this writes the next line.  This function 
	 *  exists to avoid the use of try-catch blocks in the caller code.
	 */
	public static void appendToBufferedWriter(BufferedWriter out, String stringToWrite, boolean writeNewLine) {
		try { 
			out.append(stringToWrite);
			if (writeNewLine) {
				out.newLine();
			}
		} catch (IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}
	}

	// ========================================================================
	/** Given a BufferedWriter, this closes the stream. */
	public static void closeBufferedWriter(BufferedWriter out) {
		try { 
			out.close();
		} catch (IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}
	}

	// ========================================================================
	/** Given a BufferedReader, this returns the next line.  This function
	 *  exists to avoid the use of try-catch blocks and thus may speed up
	 *  performance.
	 */
	public static String getNextLineInBufferedReader(BufferedReader in) {
		try { 
			return in.readLine();
		} catch (IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}
		return null;
	}

	// ========================================================================
	/** Given a BufferedReader, this closes the stream. */
	public static void closeBufferedReader(BufferedReader in) {
		try { 
			in.close();
		} catch (IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}
	}

	// ========================================================================
	/** Creates a new file and writes to it given the arguments. */
	public static void writeOutputFile(String outFilename, Collection<String> stringsToWrite) {
		try {
			//StringBuilder sb = new StringBuilder(1048576 / 4);			
			
			BufferedWriter out = new BufferedWriter(new FileWriter(outFilename), 1048576);
			//BufferedWriter out = new BufferedWriter(new OutputStreamWriter(new BufferedOutputStream(new FileOutputStream(outFilename), 1048576)), 1048576);
			
			for (Iterator<String> iter = stringsToWrite.iterator(); iter.hasNext(); ) {
				//sb.append(iter.next());
				//sb.append("\n");
				String lineToWrite = iter.next();
				out.write(lineToWrite, 0, lineToWrite.length());
				//out.write(iter.next());
				out.newLine();
			}
			//char[] chars = new char[sb.length()];
			//sb.getChars(0, sb.length(), chars, 0);
			//FileWriter out = new FileWriter(outFilename);
			//out.write(chars);
			out.flush();
			out.close();
		} catch (IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}
	}

	// ========================================================================
	/** Creates a new file and writes to it given the argumetns. */
	public static void writeOutputFile(String outFilename, String[] stringsToWrite) {
		try {
			BufferedWriter out = new BufferedWriter(new FileWriter(outFilename));
			
			for (int i = 0; i < stringsToWrite.length; i++) {
				out.write(stringsToWrite[i]);
				out.newLine();
			}
			out.close();
		} catch (IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}
	}

	// ========================================================================
	/** Given an open BufferedWriter, this writes the lines passed in into the bufferedwriter. */
	public static void writeOutputFile(BufferedWriter out, Collection<String> stringsToWrite) {
		try {
			for (Iterator<String> iter = stringsToWrite.iterator(); iter.hasNext(); ) {
				out.write(iter.next());
				out.newLine();
			}
		} catch (IOException e) {
			e.printStackTrace();
			System.exit(-1);
		} 
	}

	// ========================================================================
	public static void writeOutputFile(BufferedWriter out, String[][] strMatrix, char separatorChar) {
		try {
			for (int i = 0; i < strMatrix.length; i++) {
				String[] rowStrArray = strMatrix[i];
				for (int j = 0; j < rowStrArray.length; j++) {
					out.write(rowStrArray[j]);
					if (j < rowStrArray.length - 1) { out.write(separatorChar); }										
				}
				out.newLine();
			}
		} catch (IOException e) {
			e.printStackTrace();
			System.exit(-1);
		} 
	}

	// ========================================================================
	public static void writeOutputFile(BufferedWriter out, int[][] intMatrix, String[] rowHeader, String[] colHeader, char separatorChar) {
		try {
			// First do the header for the columns
			if (colHeader != null) {
				for (int i = 0; i < colHeader.length; i++) {
					if (i > 0 || rowHeader != null) {  							
						out.write(separatorChar);
					}
					out.write(colHeader[i]);
				}
			}
			out.newLine();
			
			// Now write the matrix
			for (int i = 0; i < intMatrix.length; i++) {								
				// Write out the row header if it isn't null
				if (rowHeader != null) {
					out.write(rowHeader[i]);
					out.write(separatorChar);
				}
				
				int[] rowIntArray = intMatrix[i];
				for (int j = 0; j < rowIntArray.length; j++) {					
					out.write("" + rowIntArray[j]);
					if (j < rowIntArray.length - 1) { out.write(separatorChar); }										
				}
				out.newLine();
			}
		} catch (IOException e) {
			e.printStackTrace();
			System.exit(-1);
		} 
	}

	// ========================================================================
	/** Returns all the lines from the file in the array, including the header string. */
	public static ArrayList<String> readAllLinesFromFile(String inFilename, boolean trimLines, boolean ignoreWhitespaceLines) {
		return readAllLinesFromFile(inFilename, trimLines, ignoreWhitespaceLines, null);
	}
	
	// ========================================================================
	/** Given an infilename, this reads all the lines in the infilename and returns a list of strings. */
	public static ArrayList<String> readAllLinesFromFile(String inFilename, boolean trimLines, boolean ignoreWhitespaceLines, StringBuilder headerString) {
		LinkedList<String> theList = new LinkedList<String>();
		BufferedReader in = getBufferedReader(inFilename);
		
		int lineCounter = -1;
		String line;		
		while ((line = getNextLineInBufferedReader(in)) != null) {			
			String lineTrimmed = trimLines ? line.trim() : line;			
			boolean isBlankLine = lineTrimmed.equals("");
			if (!isBlankLine || (isBlankLine && !ignoreWhitespaceLines)) {
				++lineCounter;
				if ((lineCounter == 0) && (headerString != null)) {
					headerString.setLength(0);
					headerString.append(lineTrimmed);				
				} else {
					theList.add(lineTrimmed);					
				}
			}
		}
		
		closeBufferedReader(in);
		return new ArrayList<String>(theList);
	}

	// ========================================================================
	public static ArrayList<String> readAllLinesFromFile(String inFilename) {
		return readAllLinesFromFile(inFilename, false, false);
	}

	// ========================================================================
	/** Given an infilename with data in columnar format, this reads all the lines
	 *  and returns the data as a matrix. */
	public static String[][] readAllLinesFromFileAsMatrix(String inFilename) {
		return IOUtils.readAllLinesFromFileAsMatrix(inFilename, "\\s");
	}

	// ========================================================================
	public static String[][] readAllLinesFromFileAsMatrix(String inFilename, String regex) {
		int numLinesInFile = IOUtils.countNumberLinesInFile(inFilename);
		String[][] rV = new String[numLinesInFile][];
		
		BufferedReader in = getBufferedReader(inFilename);
		for (int i = 0; i < numLinesInFile; i++) {
			String line = getNextLineInBufferedReader(in);
			rV[i] = line.split(regex);
		}
		closeBufferedReader(in); 
		return rV;
	}

	// ========================================================================
	public static double[][] readAllLinesFromFileAsMatrixDouble(String inFilename, int startingRow) {
		return IOUtils.readAllLinesFromFileAsMatrixDouble(inFilename, "\\s", startingRow);
	}

	// ========================================================================
	public static double[][] readAllLinesFromFileAsMatrixDouble(String inFilename, String regex, int startingRow) {
		int numLinesInFile = IOUtils.countNumberLinesInFile(inFilename);
		double[][] rV = new double[numLinesInFile][];
		
		BufferedReader in = getBufferedReader(inFilename);
		for (int i = 0; i < startingRow; i++) { getNextLineInBufferedReader(in); }  // pass over unwanted lines 
		
		for (int i = 0; i < numLinesInFile - startingRow; i++) {
			String line = getNextLineInBufferedReader(in);
			String[] components = line.split(regex);
			rV[i] = new double[components.length];
			
			for (int j = 0; j < components.length; j++) {
				rV[i][j] = Double.parseDouble(components[j]);
			}
		}
		closeBufferedReader(in); 
		return rV;
	}

	// ========================================================================
	/** Given an input filename with two columns, this creates a hashtable from the 
	 *  file, with the key as the 0th column and the value as the 1st column.
	 */
	public static Hashtable<String, String> readHashtableFromFile(String inFilename, boolean shouldTrim) {
		Hashtable<String, String> newTable = new Hashtable<String, String>(); 
		BufferedReader in = getBufferedReader(inFilename);
		String line;
		while ((line = getNextLineInBufferedReader(in)) != null) {
			String[] components = line.split("\\s");
			if (shouldTrim) {
				newTable.put(components[0].trim(), components[1].trim());
			} else {
				newTable.put(components[0], components[1]);
			}
		}
		closeBufferedReader(in);
		return newTable;
	}

	// ========================================================================
	/** Returns the header row of a file. */
	public static String readHeaderRow(String inFilename) {
		BufferedReader in = getBufferedReader(inFilename);
		String headerLine = getNextLineInBufferedReader(in);
		closeBufferedReader(in);
		return headerLine;
	}

	// ========================================================================
	/** Prints a matrix to screen -- used as a debug routine. */
	public static void printStringMatrix(String[][] theMatrix, String outFilename) {
		BufferedWriter out = getBufferedWriter(outFilename);
		
		for (int row = 0; row < theMatrix.length; row++) {
			String[] line = theMatrix[row];
			for (int col = 0; col < line.length; col++) {
				writeToBufferedWriter(out, line[col], false);
				if (col < line.length - 1) { writeToBufferedWriter(out, "\t", false); }
			}
			writeToBufferedWriter(out, "", true);
		}
		closeBufferedWriter(out);
	}

	// ========================================================================
	/** Counts the number of lines in a file. */
	public static int countNumberLinesInFile(String inFilename) {		
		BufferedReader in = getBufferedReader(inFilename);
		int lineCount = 0;
		while (getNextLineInBufferedReader(in) != null) {
			lineCount++;
		}		
		closeBufferedReader(in);
		return lineCount;
	}

	// ========================================================================
	/** Creates a junk file of the specified length. */
	public static void createJunkFile(int numBytes, String outFilename) {		
		try {
			FileWriter out = new FileWriter(outFilename);
			for (int i = 0; i < numBytes; i++) {
				out.write((byte) randomGen.nextInt(128));
			}
			out.close();
		} catch (IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}
	}

	// ========================================================================
	/** Given a series of files, this takes in subsequent files and merges.  */
	public static void horizontallyConcatenateFiles(String pathName, char columnSeparatorChar) {
		String dirName = ((new File(pathName)).isDirectory()) ? pathName : (new File(pathName)).getParent();
		String[] fileListing = (new File(dirName)).list();
		if (fileListing.length == 0) return;
		String outFilename = dirName + File.separatorChar + "AllFilesAppended.txt";
		BufferedReader[] allIn = new BufferedReader[fileListing.length];
		BufferedWriter out = getBufferedWriter(outFilename);
			
		for (int i = 0; i < allIn.length; i++) {
			allIn[i] = getBufferedReader(dirName + File.separatorChar + fileListing[i]);			
		}
		
		String columnSeparatorStr = "" + columnSeparatorChar;
		boolean loopFlag = true;
		while (loopFlag) {
			for (int i = 0; i < allIn.length; i++) {
				String line = getNextLineInBufferedReader(allIn[i]);
				if (line == null) {
					loopFlag = false;
					break;  // don't let anything else be added
				} else {
					if (i > 0) {
						writeToBufferedWriter(out, columnSeparatorStr, false);
					}
					writeToBufferedWriter(out, line, (i == allIn.length - 1));					
				}
			}			
		}
	
		for (int i = 0; i < allIn.length; i++) {
			closeBufferedReader(allIn[i]);
		}
		closeBufferedWriter(out);
	}

	// ========================================================================
	/** Given a series of 2-column files, this takes in subsequent files and merges.  This assumes
	 *  that the files have the same number of rows. 
	 *  @param pathName The name of the directory containing all the files to concatenate
	 */
	public static void horizontallyConcatenateFiles(String pathName, ArrayList<Integer> columnsToAppend, char columnSeparatorChar) {
		
		String dirName = ((new File(pathName)).isDirectory()) ? pathName : (new File(pathName)).getParent();
		String[] fileListing = (new File(dirName)).list();
		if (fileListing.length == 0) return;
		String outFilename = dirName + File.separatorChar + "AllFilesAppended.txt";
		StringBuilder[] sbLines = null;	  // the output data structure
		columnsToAppend = (columnsToAppend ==  null) ? new ArrayList<Integer>() : columnsToAppend;
		
		String line;
		for (int i = 0; i < fileListing.length; i++) {
			String inFilename = dirName + File.separatorChar + fileListing[i];
			System.out.println(inFilename);
			BufferedReader in = getBufferedReader(inFilename);
			
			// Initialize the output data structures
			if (i == 0) {
				ArrayList<String> allLinesFirstFile = readAllLinesFromFile(inFilename);
				sbLines = new StringBuilder[allLinesFirstFile.size()];
				for (int j = 0; j < sbLines.length; j++) {					
					int lineLength = allLinesFirstFile.get(j).length();
					int predictedLength = lineLength * fileListing.length;
					sbLines[j] = new StringBuilder(predictedLength);
					sbLines[j].append(allLinesFirstFile.get(j));
				}
			} else {
				int sbIndex = 0;
				while ((line = getNextLineInBufferedReader(in)) != null) {
					if (columnsToAppend.size() > 0) {
						String[] components = line.split("\t");
						for (int j = 0; j < columnsToAppend.size(); j++) {
							int theCol = columnsToAppend.get(j).intValue();
							String stringToAdd = theCol < components.length ? components[theCol] : "NA";
							sbLines[sbIndex].append(columnSeparatorChar);
							sbLines[sbIndex].append(stringToAdd);
						}						
					} else {
						sbLines[sbIndex].append(columnSeparatorChar);
						sbLines[sbIndex].append(line);
					}
					sbIndex++;
				}
			}			
			
			closeBufferedReader(in);
		}
		
		// Now, write to output
		BufferedWriter out = getBufferedWriter(outFilename);
		for (int i = 0; i < sbLines.length; i++) {
			writeToBufferedWriter(out, sbLines[i].toString(), true);
		}
		closeBufferedWriter(out);
	}

	// ========================================================================
	/** Runs an external command.  It returns two BufferedReaders.  The first
	 *  (at index 0) is the standard input stream, while the second (at index 1)
	 *  is the standard error.
	 */
	public static BufferedReader[] runExternalCommand(String command) {            
		try {	
			Process p = Runtime.getRuntime().exec(command);
			BufferedReader stdInput = new BufferedReader(new InputStreamReader(p.getInputStream()));
			BufferedReader stdError = new BufferedReader(new InputStreamReader(p.getErrorStream()));
			
			BufferedReader[] rV = new BufferedReader[2];
			rV[0] = stdInput;
			rV[1] = stdError;
			return rV;
			
		} catch (IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}
		return null;
	}

}

