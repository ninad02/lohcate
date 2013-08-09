package JISTIC;

import java.util.*;
import java.io.*;

public class convertBirdsuite {
	public static void main(String[] args) {
		try {
			List<String> samples = new Vector<String>();
			BufferedReader infile = new BufferedReader(new FileReader(args[0]));
			for (String next = infile.readLine(); next.split("\\t")[3]
					.equals("1"); next = infile.readLine())
				if (!samples.contains(next.split("\\t")[0]))
					samples.add(next.split("\\t")[0]);
			infile.close();
			String header_line = "Marker\tChromosome\tPosition";
			for (Iterator<String> iter = samples.iterator(); iter.hasNext();)
				header_line += "\t" + iter.next();
			System.out.println(header_line);
			int[] ones=new int[samples.size()];
			Arrays.fill(ones, 1);
			BufferedReader markerfile = new BufferedReader(new FileReader(
					args[1]));
			markerfile.readLine();
			infile = new BufferedReader(new FileReader(args[0]));
			Marker next_marker = null;			
			String next_line = infile.readLine();
			for (int chromosome = 1; chromosome <= Marker.numchrom; chromosome++) {
				List<Marker> markers = new Vector<Marker>();
				if (next_marker != null)
					markers.add(next_marker);
				for (String next = markerfile.readLine(); next != null; next = markerfile
						.readLine()) {
					Marker marker = new Marker(next, false,ones);
					marker.copy_number = new float[samples.size()];
					Arrays.fill(marker.copy_number, 0);
					if (marker.chromosome == chromosome)
						markers.add(marker);
					else {
						next_marker = marker;
						break;
					}
				}
				for (String next = next_line; next != null; next = infile
						.readLine()) {
					String[] tokens = next.split("\\t");
					if (Gene.translateChr(tokens[3]) != chromosome) {
						next_line = next;
						break;
					}
					int sample = samples.indexOf(tokens[0]);
					int val = Integer.parseInt(tokens[2]) - 2;
					int start = Collections.binarySearch(markers, new Marker(
							"dummy\t" + tokens[3] + '\t' + tokens[4], false,ones));
					if (start < 0)
						start = -start - 1;
					int end = Collections.binarySearch(markers, new Marker(
							"dummy\t" + tokens[3] + '\t' + tokens[5], false,ones));
					if (end < 0)
						end = -end - 2;
					for (int i = start; i <= end; i++)
						markers.get(i).copy_number[sample] = val;
				}
				for (Iterator<Marker> iter = markers.iterator(); iter.hasNext();) {
					Marker marker = iter.next();
					String marker_line = marker.name + '\t' + marker.chromosome
							+ '\t' + marker.start;
					boolean found_aberrations = false;
					for (int i = 0; i < marker.copy_number.length; i++) {
						int val = Math.round(marker.copy_number[i]);
						marker_line += "\t" + val;
						if (val != 0)
							found_aberrations = true;
					}
					if (found_aberrations)
						System.out.println(marker_line);
				}
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
