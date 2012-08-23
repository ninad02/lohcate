import java.util.ArrayList;

import shared.Utils;
/**
 * LOHcate --- A software tool for LOH calling and visualization in cancer genomes
 * D Wheeler & SG Reddy
 * Human Genome Sequencing Center, Baylor College of Medicine (Houston, TX)
 * 
 * @author Siddharth G. Reddy
 *
 */
public class Gene {
	public int[] arr, counts;
	public int min, max;
	public String lbl, chr;
	public ArrayList<ArrayList<String>> patients;
	public Gene(String param, String chr_param) {
		patients = new ArrayList<ArrayList<String>>();
		lbl = param;
		chr = chr_param;
		arr = new int[4]; //stores hit counts for nonsynonymous, synonymous, germline, somatic variants
		counts = new int[Script.cluster_names.length]; //stores hit counts for LOH, dup, &c. variants
		for (int i = 0; i<counts.length; i++)
			counts[i] = 0;
		for (int i = 0; i<arr.length; i++)
			arr[i] = 0;
		min = Integer.MAX_VALUE;
		max = Integer.MIN_VALUE;
	}
	public float getRecurrence(int ind, int total) { return (float)patients.get(ind).size() / (float)total; }
	public float getArrElem(int ind) { 
		float rtn = (float)arr[ind] / (float)(max - min); 
		if (Float.toString(rtn).equals("Infinity")) //div by 0
			rtn = arr[ind]; //pretend that the denominator is 1 (I think this is OK to do)
		if (Float.toString(rtn).equals("NaN")) //div results in an indeterminate form
			rtn = 0; //let's just call it 0, since this only happens when numerator (hit counts) == 0
		return rtn;
	}
	public String toString() {
		String rtn = "";
		for (int i = 0; i<arr.length; i++)
			rtn += Utils.log(arr[i]) + ","; //log normalize to make histogram readable
		for (int i = 0; i<counts.length; i++)
			rtn += Utils.log(counts[i]) + ","; //...
		for (int i = 0; i<patients.size(); i++)
			rtn += Integer.toString(patients.get(i).size()) + ","; //# patients ==> recurrence figure, straight-up (not a fraction)
		return rtn.substring(0, rtn.length() - 1);
	}
}