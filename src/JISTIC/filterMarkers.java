package JISTIC;
import java.util.*;
import java.io.*;

public class filterMarkers {

    static class CopyNumberData extends GenomeLocation {
	CopyNumberData(String data) {
	    super(Gene.translateChr(data.split("\\t")[1]),
		(int)Math.round(Double.parseDouble(data.split("\\t")[2])),
		(int)Math.round(Double.parseDouble(data.split("\\t")[2])));
	    String[] tokens = data.split("\\t");	
	    name = tokens[0];
	    this.data = data;
        }

	String name;
	String data;
	boolean bad = false;
    }

    public static void main(String[] args)
    {
    try {
      String[] bad_region_file_names=null;
      BufferedReader bad_region_file = null;
      String next_bad_region_line = null;
      String input_matrix = null;
      String output_file=null;
      for(int i=0; i < args.length; i++)
	if(args[i].toLowerCase().startsWith("badregions="))
		bad_region_file_names=args[i].substring(11).split(","); 	
	else if(args[i].toLowerCase().startsWith("excludedregions="))
		bad_region_file_names=args[i].substring(16).split(","); 	
	else if(args[i].toLowerCase().startsWith("outputfile="))
		output_file=args[i].substring(11);
	else if(input_matrix == null)
		input_matrix = args[i];	
	else {
	    System.err.println("Invalid command-line argument: " + args[i]);
	    System.exit(1);
	}

      
      if(bad_region_file_names == null)
	    System.err.println("no excluded regions specified");

      PrintStream resultStream=null;
		if(output_file!=null)
		{
			resultStream=new PrintStream(new FileOutputStream(output_file));
		}else
		{
			resultStream=System.out;
		}
		
		
      BufferedReader markerfile = new BufferedReader(
			input_matrix == null ?
				new InputStreamReader(System.in) :
				new FileReader(input_matrix));
      resultStream.println(markerfile.readLine());
      CopyNumberData next_marker = null;
      for(int chromosome = 1; chromosome <= Marker.numchrom; chromosome++){
	List<CopyNumberData> markers = new Vector<CopyNumberData>();
	if(next_marker != null)
	    markers.add(next_marker);
	for(String next = markerfile.readLine(); next != null; next = markerfile.readLine()){
	    CopyNumberData marker = new CopyNumberData(next);
	    if(marker.chromosome == chromosome)
		markers.add(marker);
	    else {
		next_marker = marker;
		break;
	    }
	}

	for (int j = 0; j < bad_region_file_names.length; j++) {
	bad_region_file = new BufferedReader(new FileReader(bad_region_file_names[j]));	
	bad_region_file.readLine(); 		
	for(String next = bad_region_file.readLine(); next != null; next = bad_region_file.readLine())
	{
	    String[] tokens = next.split("\\t");
	    CopyNumberData bad_region_start = new CopyNumberData("dummy\t" + tokens[1] + '\t' + tokens[2]);
	    if(bad_region_start.chromosome != chromosome)
		continue;
	    int start = Collections.binarySearch(markers, bad_region_start);
	    if(start < 0)
		start = - start - 1;
	    
	    int end = Collections.binarySearch(markers,
				new CopyNumberData("dummy\t" + tokens[1] + '\t' + tokens[3]));
	    if(end < 0)
		end = - end - 2;
	    for(int i = start; i <= end; i++)
		markers.get(i).bad = true;
	}
	bad_region_file.close();
	}

	for(Iterator<CopyNumberData> iter = markers.iterator(); iter.hasNext(); ){
	    CopyNumberData marker = iter.next();
	    if(marker.bad)
		continue;
	    resultStream.println(marker.data);
	}
      }
    } catch (Exception e) {
	e.printStackTrace();
    }
    }
}
