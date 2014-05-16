package shared;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import java.io.PrintStream;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Collections;

public class nexusDxTreeParser {

	private static final Character Dot        = '.';
	private static final Character Underscore = '_';
	private static final Character Hyphen     = '-';
	private static final Character FileSeparatorUnixChr = '/';
	//private static final String    FileSeparatorUnixStr = "/";
	//private static final Character PathSeparatorUnix = ':';
	private static final char      DoubleQuote = '"';
	
	private static final int NumSpacesInIndent = 4;
	
	// ========================================================================
	public static ListingsPackage parse(String fileToParse) {
		ArrayList<String> allLines = parseFile(fileToParse);

		// Init the return value
		ListingsPackage listingsPackage = new ListingsPackage();
		
		// Initialize the parent path string with the separator
		ArrayList<String> fullPath = new ArrayList<String>(50);
		int levelPrev = 0;
		StringBuilder sb = new StringBuilder(1024);
		DxTreeFile fileCurrent = listingsPackage.mRoot;
		
		for (int row = 0; row < allLines.size(); row++) {
			String line = allLines.get(row);
			
			//System.out.println(line);			
			if ((line.length() == 1) && (line.charAt(0) == Dot)) {
				// Do nothing -- use this as a placeholder for debugging
				
			} else {				
				// We determine what the subdir level is heuristically
				int levelTest = 0;
				int strIndex = 0;
				
				// One-line for-loop
				for ( ; (strIndex < line.length() && !isValidFirstChar(line.charAt(strIndex))); strIndex += NumSpacesInIndent, levelTest++);
				
				// If we don't find a valid first filename character on the line, we error and exit
				if (strIndex >= line.length()) {
					System.err.printf("ERROR: Could not find valid entry on row (%d) in file (%s): %s\n", row, fileToParse, line);
					System.exit(-1);
				}
				
				// Pull the file or directory name and check the level difference of the previous
				String fileOrDirName = line.substring(strIndex);
				boolean hasWhitespace = containsWhitespace(fileOrDirName); 
				
				int levelDiff = levelTest - levelPrev;
				levelPrev = levelTest;
				int fullPathLastIndex = fullPath.size() - 1;
			
				//System.out.println("level: " + levelDiff);
				if (levelDiff > 1) {
					System.err.println("Badly constructed directory tree.  Subdirectory cannot be more than 1 level deeper than parent directory!");
					System.exit(-1);
					
				} else if (levelDiff == 1) {
					fullPath.add(fileOrDirName);
					
					DxTreeFile child = new DxTreeFile(fileOrDirName, fileCurrent);
					fileCurrent.addChildSorted(child);
					fileCurrent = child;
					
				} else if (levelDiff == 0) {
					// Replace the last element of the full path with this one
					fullPath.set(fullPathLastIndex, fileOrDirName);
										
					DxTreeFile child = new DxTreeFile(fileOrDirName, fileCurrent.getParent());
					fileCurrent.getParent().addChildSorted(child);
					fileCurrent = child;
					
				} else {
					// Remove the revelent components of the parent paths and add the new subdir/file
					for (int i = fullPathLastIndex; i >= fullPathLastIndex + levelDiff; i--) {
						fullPath.remove(i);
						fileCurrent = fileCurrent.getParent();
					}
					
					fullPath.add(fileOrDirName);
					DxTreeFile child = new DxTreeFile(fileOrDirName, fileCurrent);
					fileCurrent.addChildSorted(child);
					fileCurrent = child;
				}
				
				String fullPathStr = constructFullPath(fullPath, FileSeparatorUnixChr, hasWhitespace, sb);
				listingsPackage.mFullPathsAll.add(fullPathStr);
			}
		}
		
		for (String s: listingsPackage.mFullPathsAll) {
			System.out.println("PATH: " + s);
		}
		
		listingsPackage.mRoot.print(System.out, "");
		
		return listingsPackage;
	}

	// ========================================================================
	public static boolean containsWhitespace(String s) {
		int length = s.length();
		for (int i = 0; i < length; i++) {
			if (Character.isWhitespace(s.charAt(i))) return true;
		}
		return false;
	}

	// ========================================================================
	public static String constructFullPath(ArrayList<String> fullPath, Character fileSeparator, boolean placeDoubleQuotes) {		
		return constructFullPath(fullPath, fileSeparator, placeDoubleQuotes, new StringBuilder(1024));
	}
	
	// ========================================================================
	public static String constructFullPath(ArrayList<String> fullPath, Character fileSeparator, boolean placeDoubleQuotes, StringBuilder sb) {
		sb.setLength(0);
		if (placeDoubleQuotes) sb.append(DoubleQuote);
		
		for (String fileOrDirName : fullPath) {
			sb.append(fileSeparator).append(fileOrDirName);			
		}		
		
		if (placeDoubleQuotes) sb.append(DoubleQuote);
		return sb.toString();
	}
	
	// ========================================================================
	private static boolean isValidFirstChar(char ch) {
		return
			   Character.isAlphabetic(ch) 
			|| Character.isDigit(ch)
			|| (ch == Dot)
			|| (ch == Underscore)
			|| (ch == Hyphen)
			;	
	}
	

	// ========================================================================
	private static ArrayList<String> parseFile(String dxTreeFileToParse) {
		BufferedReader in = null;
		BufferedWriter out = null;
		
		Charset utf8 = Charset.forName("UTF-8");
		
		try {
			in = Files.newBufferedReader(Paths.get(dxTreeFileToParse), utf8);
			out = Files.newBufferedWriter(Paths.get(dxTreeFileToParse + ".out.txt"), utf8);
		} catch (IOException e) {
			printStackTraceAndExit(e);
		}
		
		// Open file
		/*
		try {
			in = new BufferedReader(new InputStreamReader(new FileInputStream(dxTreeFileToParse), "UTF-8"));
		} catch (FileNotFoundException e) {
			printStackTraceAndExit(e);
		}*/
		
		ArrayList<String> lines = new ArrayList<String>(1000);
		String line;
		
		// Ideally, shouldn't do much processing within try/catch, but using
		// bad programming style as a shortcut here
		try {
			while ((line = in.readLine()) != null) {
				lines.add(line);
				out.write(line);
				out.newLine();
			}
				
			in.close();
			out.close();
		} catch (IOException e) {
			printStackTraceAndExit(e);	
		}
		
		return lines;
	}
	

	// ========================================================================
	private static void printStackTraceAndExit(Exception e) {		
		e.printStackTrace();
		System.exit(-1);
	}
	
	// ========================================================================
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		parse(args[0]);

	}
	
	// ========================================================================
	public static enum Traversal {
		WentBelow, WentAbove, StayedSame;
	}
	
	// ========================================================================
	public static class DxTreeFile implements Comparable<DxTreeFile> {
		String mName;
		DxTreeFile mParent;
		ArrayList<DxTreeFile> mChildren;
		
		public DxTreeFile(String name, DxTreeFile parent) {
			mName = name;
			mParent = parent;
			mChildren = new ArrayList<DxTreeFile>();
		}

		public DxTreeFile getParent() { return mParent; }
		public boolean hasChildren() { return !mChildren.isEmpty(); }

		public DxTreeFile find(String nameOfChild) {
			for (DxTreeFile dtf : mChildren) {
				if (dtf.mName.equals(nameOfChild)) return dtf;
			}
			return null;
		}
		
		public boolean addChildSorted(DxTreeFile child) {
			int resultIndex = Collections.binarySearch(mChildren, child);
			if (resultIndex < 0) {
				int insertIndex = -(resultIndex + 1);
				mChildren.add(insertIndex, child);
				return true;
			} else {
				return false;
			}			
		}
		
		public void print(PrintStream out, String prefix) {
			out.println(prefix + mName);
			for (DxTreeFile dtf : mChildren) {
				dtf.print(out, prefix + "\t");
			}
		}
		
		@Override
		public int compareTo(DxTreeFile rhs) {
			return mName.compareTo(rhs.mName);
		}
		
		/* (non-Javadoc)
		 * @see java.lang.Object#equals(java.lang.Object)
		 */
		@Override
		public boolean equals(Object obj) {
			if (this == obj)
				return true;
			if (obj == null)
				return false;
			if (getClass() != obj.getClass())
				return false;
			DxTreeFile other = (DxTreeFile) obj;
			if (mName == null) {
				if (other.mName != null)
					return false;
			} else if (!mName.equals(other.mName))
				return false;
			return true;
		}
	}
	
	// ========================================================================
	public static class ListingsPackage {
		ArrayList<String> mFullPathsAll;
		DxTreeFile mRoot;
		
		public ListingsPackage() {
			mFullPathsAll = new ArrayList<String>(1024);
			mRoot = new DxTreeFile("", null);
		}
	}

}


