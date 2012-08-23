package shared;

import java.util.ArrayList;

/**
 * LOHcate --- A software tool for LOH calling and visualization in cancer genomes
 * D Wheeler & SG Reddy
 * Human Genome Sequencing Center, Baylor College of Medicine (Houston, TX)
 * 
 * @author Siddharth G. Reddy
 *
 */
public class Utils {
	
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
	
}

