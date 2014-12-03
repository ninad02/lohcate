package JISTIC;
import java.util.*;
import java.io.*;

public class convertSEG {

	static boolean defaultMarkerStatus;

	static class CopyNumberData extends GenomeLocation {
		CopyNumberData(String data) {
			super(Gene.translateChr(data.split("\\t")[1]),
					(int)Math.round(Double.parseDouble(data.split("\\t")[2])),
					(int)Math.round(Double.parseDouble(data.split("\\t")[2])));
			String[] tokens = data.split("\\t");	
			name = tokens[0];
		}

		float[] copy_number;
		String name;
		boolean bad = convertSEG.defaultMarkerStatus;
	}

	public static void main(String[] args)
	{
		try {					
			String logfile="convertSEGCall.log";
			String[] bad_region_file_names=null;
			String[] good_region_file_names=null;
			String output_file=null;
			boolean remove_zero = true;
			boolean verbose = false;
			int maxsegsize=-1;
			int minsegsize=-1;
			String sample_list_file = null;
			for(int i=2; i < args.length; i++) {
				if(args[i].equalsIgnoreCase("excluderemovedmarkers"))
					remove_zero = true;
				else if(args[i].equalsIgnoreCase("includeremovedmarkers"))
					remove_zero = false;
				else if(args[i].equalsIgnoreCase("verbose"))
					verbose = true;
				else if(args[i].toLowerCase().startsWith("badregions="))
					bad_region_file_names=args[i].substring(11).split(","); 	
				else if(args[i].toLowerCase().startsWith("excludedregions="))
					bad_region_file_names=args[i].substring(16).split(",");
				else if(args[i].toLowerCase().startsWith("maxsegsize="))
					maxsegsize=Integer.parseInt(args[i].substring(11));
				else if(args[i].toLowerCase().startsWith("minsegsize="))
					minsegsize=Integer.parseInt(args[i].substring(11));
				else if(args[i].toLowerCase().startsWith("includedregions="))
					good_region_file_names=args[i].substring(16).split(","); 	
				else if (args[i].toLowerCase().startsWith("numchrom=")) {
					Marker.numchrom = Integer.parseInt(args[i].substring(9));
				}else if(args[i].toLowerCase().startsWith("samplestouse="))
					sample_list_file=args[i].substring(13);
				else if(args[i].toLowerCase().startsWith("outputfile="))
				{
					output_file=args[i].substring(11);
					logfile=output_file+".log";
				}
				else {
					System.err.println("Invalid command-line argument: " + args[i]);
					System.exit(1);
				}
			}

			
			PrintStream fout = new PrintStream(new FileOutputStream(logfile));
			for (int i = 0; i < args.length; i++) {
				fout.print(" "+args[i]);
			}
			fout.close();
			
			PrintStream resultStream=null;
			if(output_file!=null)
			{
				resultStream=new PrintStream(new FileOutputStream(output_file));
			}else
			{
				resultStream=System.out;
			}
			
			if(bad_region_file_names != null && good_region_file_names != null){
				System.err.println("Can specify either excluded region files  or included region files, not both");
				System.exit(1);
			}

			defaultMarkerStatus = good_region_file_names != null;

			if(verbose){
				System.err.println("SEG file - " + args[0]);
				System.err.println("probe locations - " + args[1]);
				if(bad_region_file_names != null || good_region_file_names != null){
					System.err.print(defaultMarkerStatus ? "included regions - " : "excluded regions - ");
					String[] file_names = defaultMarkerStatus ? good_region_file_names : bad_region_file_names;
					for(int i=0; i < file_names.length; i++){
						if(i > 0)
							System.err.print(",");
						System.err.print(file_names[i]);
					}
					System.err.println();
					if(!remove_zero)
						System.err.println((defaultMarkerStatus ? "all others " : "") + "zeroed and included in output");
				}
				System.err.println("samples to use - " + (sample_list_file == null ? "all" : sample_list_file));
			}
			if(sample_list_file == null)
				sample_list_file = args[0];
			List<String> samples = new Vector<String>();
			BufferedReader infile = new BufferedReader(new FileReader(sample_list_file));
			infile.readLine();
			for(String next = infile.readLine(); next != null; next = infile.readLine())
				if(!samples.contains(next.split("\\t")[0]))
					samples.add(next.split("\\t")[0]);
			infile.close();
			String header_line = "Marker\tChromosome\tPosition";
			for(Iterator<String> iter = samples.iterator(); iter.hasNext(); )
				header_line += "\t" + iter.next();
			resultStream.println(header_line);			
			BufferedReader markerfile = new BufferedReader(new FileReader(args[1]));
			// skip header line
			markerfile.readLine();
			CopyNumberData next_marker = null;
			int total_markers = 0, total_bad = 0;
			for(int chromosome = 1; chromosome <= Marker.numchrom; chromosome++){
				infile = new BufferedReader(new FileReader(args[0]));
				infile.readLine();
				List<CopyNumberData> markers = new Vector<CopyNumberData>();
				if(next_marker != null)
					markers.add(next_marker);
				for(String next = markerfile.readLine(); next != null; next = markerfile.readLine()){
					CopyNumberData marker = new CopyNumberData(next);
					marker.copy_number = new float[samples.size()];
					Arrays.fill(marker.copy_number, 0);
					if(marker.chromosome == chromosome)
						markers.add(marker);
					else {
						next_marker = marker;
						break;
					}
				}

				if(bad_region_file_names != null || good_region_file_names != null){
					String[] file_names = defaultMarkerStatus ? good_region_file_names : bad_region_file_names;
					for (int j = 0; j < file_names.length; j++) {


						BufferedReader region_file = new BufferedReader(new FileReader(file_names[j]));	
						region_file.readLine(); 		
						for(String next = region_file.readLine(); next != null; next = region_file.readLine())
						{
							String[] tokens = next.split("\\t");
							if(tokens.length<3)
								continue;
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
								markers.get(i).bad = !defaultMarkerStatus;
						}
						region_file.close();
					}
				}
				
				int current_pos = -1, current_sample = -1;
				for(String line1 = infile.readLine(), line2 = infile.readLine();
						line1 != null; 
						line1 = line2, line2 = infile.readLine())
				{
					String[] tokens = line1.split("\\t");

					while(line2!=null && (line2.split("\\t")[5].equalsIgnoreCase("NA")|| (minsegsize!=-1 && Integer.parseInt(line2.split("\\t")[3])-Integer.parseInt(line2.split("\\t")[2])<minsegsize)))
					{
						line2=infile.readLine();	    	
					}

					int sample = samples.indexOf(tokens[0]);
					if(Gene.translateChr(tokens[1]) != chromosome || sample < 0){
						current_sample = -1;
						continue;
					}
					float val = Float.parseFloat(tokens[5]);
					
					
					if(maxsegsize>-1 && Integer.parseInt(tokens[3])-Integer.parseInt(tokens[2])>maxsegsize)
					{
						val=0;
					}
					
					int start = Collections.binarySearch(markers,
							new CopyNumberData("dummy\t" + tokens[1] + '\t' + tokens[2]));
					if(start < 0)
						start = - start - 1;
					if(sample == current_sample && start>=markers.size())
					{
						System.err.println("Warning, segment out of markers boundaries:\n\t"+line1);
						current_sample = -1;
						continue;
					}
						
					if(sample == current_sample && !markers.get(start).bad)
						start = current_pos+1;
					current_sample = sample;
					int end_pos = Integer.parseInt(tokens[3]);
					int end = Collections.binarySearch(markers,
							new CopyNumberData("dummy\t" + tokens[1] + '\t' + end_pos));
					if(end < 0) {
						System.out.println("XX: End: " + end_pos + "\t" + markers.get(0).toString() + "\t" + markers.get(markers.size() - 1).toString());
						end = - end - 2;
					}

					if(line2 != null && !markers.get(end).bad){
						tokens = line2.split("\\t");
						if(Gene.translateChr(tokens[1]) == chromosome &&
								tokens[0].equals(samples.get(sample)))
						{
							int next_start = Collections.binarySearch(markers,
									new CopyNumberData("dummy\t" + tokens[1] + '\t' + tokens[2]));
							if(next_start < 0)
								next_start = - next_start - 1;	
							
							if(next_start>=markers.size())
								end= markers.size()-1;
							else if(markers.get(next_start).bad)
								end = next_start - 1;
							else {
								end_pos = (end_pos + Integer.parseInt(tokens[2])) / 2;
								end = Collections.binarySearch(markers,
										new CopyNumberData("dummy\t" + tokens[1] + '\t' + end_pos));
								if(end < 0)
									end = - end - 2;
							}
							
						}
					}
					current_pos = end;
					for(int i = start; i <= end; i++)
						if(!markers.get(i).bad)
							markers.get(i).copy_number[sample] = val;
				}
				int chr_bad = 0;
				total_markers += markers.size();
				for(Iterator<CopyNumberData> iter = markers.iterator(); iter.hasNext(); ){
					CopyNumberData marker = iter.next();
					if(marker.bad){
						chr_bad++;
						total_bad++;
					}
					if(remove_zero && marker.bad)
						continue;
					String marker_line = marker.name + '\t' + marker.chromosome + '\t' + marker.start;
					for(int i=0; i < marker.copy_number.length; i++){
						float val = marker.copy_number[i];
						marker_line += "\t" + val;
					}
					resultStream.println(marker_line);
				}
				infile.close();
				if(verbose)
					System.err.println("chromosome " + chromosome + ' ' + markers.size() + " markers, " + chr_bad + " excluded");
			}
			if(verbose)
				System.err.println("total " + total_markers + " markers, " + total_bad + " excluded");
		} catch (Exception e) {
			e.printStackTrace(System.err);
			System.exit(1);
		}
	}
}
