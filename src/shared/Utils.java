package shared;

import genomeUtils.GenotypeUtils;

import java.util.Arrays;
import java.util.Collection;

import nutils.NumberUtils;


import lohcateEnums.ClusterType;

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
	public static final String SpaceString = " ";
	public static final String ColonString = ":";
	public static final String DoubleQuoteStr = "\"";
	public static final String DotStr = ".";
	public static final int NumAutosomes = 22;
	
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
	
	public static void TestModularity() {
		
		long numIter = Integer.MAX_VALUE - 1;

		for (int i = 0; i < 5000000; i++);

		int x = 0;
		long time3 = System.currentTimeMillis();
		for (int i = 0; i < numIter; i++) {
			x = NumberUtils.incrementNum(x);
		}
		long time4 = System.currentTimeMillis();

		
		x = 0;
		long time1 = System.currentTimeMillis();		
		for (int i = 0; i < numIter; i++) {
			x = x + 1;
		}
		long time2 = System.currentTimeMillis();
		
		
		System.out.println("Time for no function:\t" + (time2 - time1));
		System.out.println("Time for function:   \t" + (time4 - time3));
	}
	
	
	public static void main(String[] args) {
		//Test_removeNullElements();
		//Test_extractRsNumberFromLine();
		//TestExtractNthColumn();
		//TestExtractNthColumnRobust();
		//TestModularity();
		System.out.println(GenotypeUtils.calcFractionGC("CCTTGCGCAGGTG"));
	}
}

