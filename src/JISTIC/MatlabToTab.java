package JISTIC;
import java.io.*;
import java.util.*;

class MatlabToTab {

    public static void main(String[] args) {
    try {
	if(args.length != 1 && args.length != 2){
	    System.err.println("usage: MatlabToTab <Matlab output name> [<locations file name>]");
	    return;
	}

	Map<String, String> locations = new HashMap<String, String>();
	if(args.length == 2){
	    BufferedReader locfile = new BufferedReader(new FileReader(args[1]));
	    locfile.readLine();
	    for(String next = locfile.readLine(); next != null; next = locfile.readLine()){
		String[] tokens = next.split("\\t");
		locations.put(tokens[0], tokens[1] + '\t' + tokens[2]);
	    }
	}
	BufferedReader datafile = new BufferedReader(new FileReader(args[0] + "dat.txt"));
	BufferedReader samplesfile = new BufferedReader(new FileReader(args[0] + "samples.txt"));
	BufferedReader markerfile = new BufferedReader(new FileReader(args[0] + "markers.txt"));
	BufferedReader chrfile = new BufferedReader(new FileReader(args[0] + "markerchr.txt"));
	BufferedReader posfile = new BufferedReader(new FileReader(args[0] + "markerpos.txt"));
	String topline = "SNP marker\tChromosome\tPosition";
	for(String next = samplesfile.readLine(); next != null; next = samplesfile.readLine())
	    topline += "\t" + next;
	System.out.println(topline);
	for(String markername = markerfile.readLine(); markername != null; markername = markerfile.readLine()){
	    String chr = chrfile.readLine();
	    int pos = (new Double(posfile.readLine())).intValue();
	    String loc = locations.get(markername);
	    if(loc == null)
		loc = chr + '\t' + pos;
	    System.out.println(markername + '\t' + loc + '\t' + datafile.readLine());
	}
    } catch (Exception e) {
	e.printStackTrace();
	return;
    }
    }
}
