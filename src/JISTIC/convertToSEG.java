package JISTIC;

import java.util.*;
import java.io.*;

public class convertToSEG {

	public static void main(String[] args) {
		String output_file = null;
		try {
			String[] bad_region_file_names = null;
			for (int i = 1; i < args.length; i++)
				if (args[i].toLowerCase().startsWith("outputfile=")) {
					output_file = args[i].substring(11);
				}else if (args[i].toLowerCase().startsWith("numchrom=")) {
					Marker.numchrom = Integer.parseInt(args[i].substring(9));
				}else if(args[i].toLowerCase().startsWith("excludedregions="))
				{
					bad_region_file_names=args[i].substring(16).split(","); 	
				}else
				{
					System.err.println("Invalid command-line argument: "
							+ args[i]);
					System.exit(1);
				}

			PrintStream resultStream = null;
			if (output_file != null) {
				resultStream = new PrintStream(
						new FileOutputStream(output_file));
			} else {
				resultStream = System.out;
			}

			resultStream
					.println("barcode\tchromosome\tstart\tstop\tnum.mark\tseg.mean");

			Distribution.by_chromosome = true;

			SNPdata markers = new SNPdata(args[0], null, null,null);
			int[] ones=new int[markers.copy_number_samples.length];
			Arrays.fill(ones, 1);
			for (int chromosome = 1; chromosome <= Marker.numchrom; chromosome++) {
				markers.InMemory(chromosome);
				if (bad_region_file_names != null) {
					for (int j = 0; j < bad_region_file_names.length; j++) {

						BufferedReader region_file = new BufferedReader(
								new FileReader(bad_region_file_names[j]));
						region_file.readLine();
						for (String next = region_file.readLine(); next != null; next = region_file
								.readLine()) {
							Marker bad_region = new Marker(next,false,ones);
							if (bad_region.chromosome != chromosome)
								continue;
							markers.add(-Collections.binarySearch(markers,
									bad_region) - 1, bad_region);
						}
						region_file.close();
					}
				}

				for (int sample = 0; sample < markers.copy_number_samples.length; sample++)
					for (int region_start = -1;;) {
						while (++region_start < markers.size()
								&& markers.get(region_start).copy_number == null)
							;
						if (region_start >= markers.size())
							break;
						int next_region = region_start;
						while (++next_region < markers.size()
								&& markers.get(next_region).copy_number != null
								&& markers.get(next_region).copy_number[sample] == markers
										.get(region_start).copy_number[sample])
							;
						resultStream
								.println(markers.copy_number_samples[sample]
										+ '\t'
										+ chromosome
										+ '\t'
										+ markers.get(region_start).start
										+ '\t'
										+ markers.get(next_region - 1).start
										+ '\t'
										+ (next_region - region_start)
										+ '\t'
										+ markers.get(region_start).copy_number[sample]);
						region_start = next_region - 1;
					}
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
