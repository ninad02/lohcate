package JISTIC;

import java.util.*;
import java.io.*;

import JISTIC.Distribution.sparse_cell;
import JISTIC.Genome.Band;

public class SNPdata extends Vector<Marker> {

	public SNPdata(String cpfile, String lohfile, List<Genome.Band> bands, String differentialfile)
			throws IOException {
		// start differential groups
		int[] diffgroupsones=null;
		TreeSet<Marker> markers = new TreeSet<Marker>();
		byScore = new Vector<List<Integer>>();
		sample_distributions = new Vector<List<TreeSet<Distribution.sparse_cell>>>();
		for (int i = 0; i < Marker.aberrationtypes.length; i++) {
			Marker.aberration ab = Marker.aberrationtypes[i];
			byScore.add(new Vector<Integer>());
			sample_distributions
					.add(new Vector<TreeSet<Distribution.sparse_cell>>());
		}

		Marker.initial_columns = 3;
		// Sets of markers discovered that are either distal or proximal to
		// genome bands, arranged by chromosome
		Map<Integer, SortedSet<Marker>> proximal_markers = new HashMap<Integer, SortedSet<Marker>>(), distal_markers = new HashMap<Integer, SortedSet<Marker>>();
		if (cpfile != null) {
			BufferedReader infile = new BufferedReader(new FileReader(cpfile));
			String[] tokens = infile.readLine().trim().split("\\t");
			while (tokens[Marker.initial_columns].startsWith("Genetic")
					|| tokens[Marker.initial_columns].startsWith("Score"))
				Marker.initial_columns++;
			copy_number_samples = new String[tokens.length
					- Marker.initial_columns];
			System.arraycopy(tokens, Marker.initial_columns,
					copy_number_samples, 0, copy_number_samples.length);
			for (int i = 0; i < copy_number_samples.length; i++) {
				sample_distributions.get(0).add(
						new TreeSet<Distribution.sparse_cell>());
				sample_distributions.get(1).add(
						new TreeSet<Distribution.sparse_cell>());
			}
			diffgroupsones=new int[copy_number_samples.length];
			Arrays.fill(diffgroupsones, 1);
			diffgroups=diffgroupsones.clone();		
			
			// load differential groups if they exist
			if(differentialfile!=null)
			{
				int numgroup2=0;
				differentialGscore=true;
				BufferedReader ifile = new BufferedReader(new FileReader(differentialfile));
				String line=ifile.readLine();
				line=ifile.readLine();
				while(line!=null)
				{					
					if(!line.equals(""))
					{
						String[] toks = line.split("\\t");
						int index=Arrays.asList(copy_number_samples).indexOf(toks[0]);
						if(index==-1)
						{
							System.out.println("Error, differential sample not found: "+toks[0]);
							System.exit(1);
						}
						if(diffgroups[index]!=-1)
						{
							diffgroups[index]=-1;
							numgroup2++;
						}
					}
					line=ifile.readLine();
				}
				int numgroup1=diffgroups.length-numgroup2;
				// now weight by group size
				for (int i = 0; i < diffgroups.length; i++) {
					if(diffgroups[i]>0)
					{
						diffgroups[i]*=numgroup2;
					}else
					{
						diffgroups[i]*=numgroup1;
					}
				}
				
			}
			
			Marker.per_sample_thresholds[0] = new float[copy_number_samples.length];
			Marker.per_sample_thresholds[1] = new float[copy_number_samples.length];
			if (focal_thresholds) {
				Arrays.fill(Marker.per_sample_thresholds[0],
						Float.NEGATIVE_INFINITY);
				Arrays.fill(Marker.per_sample_thresholds[1],
						Float.POSITIVE_INFINITY);
				Iterator<Genome.Band> band_iter = bands.iterator();

				PrintStream broad_aberrations_out = new PrintStream(
						new FileOutputStream(Distribution.dir + "broad.matrix"));
				String arm = "";
				ArrayList<String> arms = new ArrayList<String>();
				for (Iterator iterator = bands.iterator(); iterator.hasNext();) {
					Band band = (Band) iterator.next();
					String newarm = band.chromosome + "" + band.ID.charAt(0);
					if (!newarm.equals(arm)) {
						arms.add(newarm);
						arm = newarm;
					}
				}

				float[][] per_sample_broad = new float[arms.size()][copy_number_samples.length];
				Genome.Band current_band = band_iter.next();
				Genome.Band previous_band = current_band;
				float[][] per_sample_values = new float[copy_number_samples.length][];
				per_sample_values[0] = new float[max_arm_markers];
				Arrays.fill(per_sample_values[0], 0);
				for (int i = 1; i < copy_number_samples.length; i++)
					per_sample_values[i] = per_sample_values[0].clone();
				
				for (int arm_num_markers = 0;; arm_num_markers++) {
					String next = infile.readLine();
					//System.out.println(next);					
					Marker new_marker = next == null ? null : new Marker(next
							.trim(), false,diffgroupsones);
					if (new_marker != null
							&& new_marker.compareTo(current_band) < 0)
						throw new Error(
								"Marker in input not arranged in order; "
										+ new_marker.toString()
										+ " should come before "
										+ current_band.toString() + "\t" + current_band.start + "\t" + current_band.end);
					while (new_marker != null
							&& new_marker
									.compareTo(current_band.end_location()) > 0)
						current_band = band_iter.next();

					// if finished with arm calculate the threshold
					if (new_marker == null
							|| !current_band.SameArm(previous_band)) 
					{
						if (arm_num_markers > 0) {
							String narm = previous_band.chromosome + ""
									+ previous_band.ID.charAt(0);
							int armindex = arms.indexOf(narm);
							for (int i = 0; i < copy_number_samples.length; i++) {
								per_sample_broad[armindex][i] = per_sample_values[i][arm_num_markers / 2];
								Arrays.sort(per_sample_values[i], 0,
										arm_num_markers);
								Marker.per_sample_thresholds[0][i] = Math
										.max(
												Marker.per_sample_thresholds[0][i],
												per_sample_values[i][arm_num_markers / 2]);
								Marker.per_sample_thresholds[1][i] = Math
										.min(
												Marker.per_sample_thresholds[1][i],
												per_sample_values[i][arm_num_markers / 2]);
								if (Region.debug)
									System.err
											.println("for arm "
													+ previous_band.toString()
													+ " sample "
													+ i
													+ " median "
													+ per_sample_values[i][arm_num_markers / 2]);
							}
						}
						previous_band = current_band;
						// if last marker leave loop
						if (new_marker == null)
							break;
						arm_num_markers = 0;
					}
					// if only one chrom and it is not this one skip marker
					if(Distribution.chrom2run!=-1 && Distribution.chrom2run!=current_band.chromosome)
					{
						continue;
					}
					// save the copy number for the marker
					for (int i = 0; i < copy_number_samples.length; i++)
						per_sample_values[i][arm_num_markers] = new_marker.copy_number[i];
				}
				infile.close();
				broad_aberrations_out.print("ARM");
				for (int i = 0; i < copy_number_samples.length; i++)
					broad_aberrations_out.print("\t" + copy_number_samples[i]);

				for (int i = 0; i < per_sample_broad.length; i++) {
					broad_aberrations_out.print("\n" + arms.get(i).toString());
					for (int j = 0; j < per_sample_broad[i].length; j++) {
						broad_aberrations_out.print("\t"
								+ per_sample_broad[i][j]);
					}
				}
				broad_aberrations_out.close();
				if (Region.debug)
					System.err.println("per sample thresholds AMP "
							+ Arrays.toString(Marker.per_sample_thresholds[0])
							+ " DEL "
							+ Arrays.toString(Marker.per_sample_thresholds[1]));
				infile = new BufferedReader(new FileReader(cpfile));
				infile.readLine();
			} else {
				Arrays.fill(Marker.per_sample_thresholds[0], 0);
				Arrays.fill(Marker.per_sample_thresholds[1], 0);
			}
			Marker previous_marker = null;
			for (String next = infile.readLine(); next != null; next = infile
					.readLine(), cp_num_markers++) {
				Marker new_marker = new Marker(next.trim(), false,diffgroups);
				if(Distribution.chrom2run!=-1 && Distribution.chrom2run!=new_marker.chromosome)
				{
					cp_num_markers--;
					continue;
				}
				if (Distribution.by_chromosome && previous_marker != null
						&& new_marker.chromosome < previous_marker.chromosome)
					throw new Error("Marker in input not arranged in order; "
							+ new_marker.toString() + " should come before "
							+ previous_marker.toString());
				if (bands != null) {
					// keep track of markers that are in a fake band,
					// i.e. outside the bands /in the band file.
					int marker_band = Collections.binarySearch(bands,
							new_marker);
					if (marker_band < 0)
						marker_band = -marker_band - 2;
					if (bands.get(marker_band).fake) {
						Map<Integer, SortedSet<Marker>> map = bands
								.get(marker_band).start == 0 ? proximal_markers
								: distal_markers;
						SortedSet<Marker> outside_markers = map
								.get(new_marker.chromosome);
						if (outside_markers == null)
							map.put(new_marker.chromosome,
									outside_markers = new TreeSet<Marker>());
						outside_markers.add(new_marker);
					}
				}
				previous_marker = new_marker;
				if (!Distribution.by_chromosome)
					markers.add(new_marker);
							
				byScore.get(0).add(new_marker.Gscore[0]);
				byScore.get(1).add(new_marker.Gscore[1]);
				for (int i = 0; i < copy_number_samples.length; i++)
					for (int j = 0; j < 2; j++) {
						float val = new_marker.copy_number[i];
						if (j == 0
								&& val <= Marker.per_sample_thresholds[0][i]
										+ Marker.thresholds[0]
								|| j == 1
								&& val >= Marker.per_sample_thresholds[1][i]
										+ Marker.thresholds[1])
							val = 0;
						Distribution.sparse_cell new_cell = sparse_cell.getSingletonObject(diffgroups[i]*Math.round(Math.abs(val)/ Distribution.binsize));
						Distribution.sparse_cell old_cell = sample_distributions
								.get(j).get(i).floor(new_cell);
						if (new_cell.equals(old_cell))
							old_cell.add(new_cell);
						else
							sample_distributions.get(j).get(i).add(new sparse_cell(new_cell.gscore));
					}
			}
			infile.close();
			if (Distribution.by_chromosome) {
				this.cpfile = new BufferedReader(new FileReader(cpfile));
				this.cpfile.readLine();
			}
		}
		if (lohfile != null) {
			int num_loh_markers = 0;
			Marker.initial_columns = 3;
			BufferedReader infile = new BufferedReader(new FileReader(lohfile));
			String[] tokens = infile.readLine().trim().split("\\t");
			while (tokens[Marker.initial_columns].startsWith("Genetic")
					|| tokens[Marker.initial_columns].startsWith("Score"))
				Marker.initial_columns++;
			loh_samples = new String[tokens.length - Marker.initial_columns];
			System.arraycopy(tokens, Marker.initial_columns, loh_samples, 0,
					loh_samples.length);
			for (int i = 0; i < loh_samples.length; i++)
				sample_distributions.get(2).add(
						new TreeSet<Distribution.sparse_cell>());
			
			diffgroupsones=new int[loh_samples.length];
			Arrays.fill(diffgroupsones, 1);
			
			Marker previous_marker = null;
			for (String next = infile.readLine(); next != null; next = infile
					.readLine(), loh_num_markers++) {
				Marker new_marker = new Marker(next.trim(), true,diffgroupsones);
				if(Distribution.chrom2run!=-1 && Distribution.chrom2run!=new_marker.chromosome)
				{
					loh_num_markers--;
					continue;
				}
				if (Distribution.by_chromosome && previous_marker != null
						&& new_marker.chromosome < previous_marker.chromosome)
					throw new Error("Marker in input not arranged in order; "
							+ new_marker.toString() + " should come before "
							+ previous_marker.toString());
				if (bands != null) {
					// keep track of markers that are in a fake band,
					// i.e. outside the bands /in the band file.
					int marker_band = Collections.binarySearch(bands,
							new_marker);
					if (marker_band < 0)
						marker_band = -marker_band - 2;
					if (bands.get(marker_band).fake) {
						Map<Integer, SortedSet<Marker>> map = bands
								.get(marker_band).start == 0 ? proximal_markers
								: distal_markers;
						SortedSet<Marker> outside_markers = map
								.get(new_marker.chromosome);
						if (outside_markers == null)
							map.put(new_marker.chromosome,
									outside_markers = new TreeSet<Marker>());
						outside_markers.add(new_marker);
					}
				}
				previous_marker = new_marker;
				if (!Distribution.by_chromosome) {
					Marker old_marker = markers.floor(new_marker);
					if (new_marker.equals(old_marker))
						old_marker.add(new_marker);
					else {
						if (cpfile != null) {
							throw new Error(
									"LOH and copy-number with different marker sets: currently does not work correctly");
						}
						if (!Distribution.by_chromosome)
							markers.add(new_marker);
					}
				}
				num_loh_markers++;
				byScore.get(2).add(new_marker.Gscore[2]);
				for (int i = 0; i < copy_number_samples.length; i++) {
					float val = new_marker.isLOH[i];
					Distribution.sparse_cell new_cell = sparse_cell.getSingletonObject(
							Math.round(Math.abs(val) / Distribution.binsize));
					Distribution.sparse_cell old_cell = sample_distributions
							.get(2).get(i).floor(new_cell);
					if (new_cell.equals(old_cell))
						old_cell.add(new_cell);
					else
						sample_distributions.get(2).get(i).add(new sparse_cell(new_cell.gscore));
				}
			}
			infile.close();
			// if (num_loh_markers < markers.size())
			if (num_loh_markers != cp_num_markers)
				throw new Error(
						"LOH and copy-number with different marker sets: currently does not work correctly");
			if (Distribution.by_chromosome) {
				this.lohfile = new BufferedReader(new FileReader(lohfile));
				this.lohfile.readLine();
			}
		}
		for (int i = 0; i < Marker.aberrationtypes.length; i++) {
			Collections.sort(byScore.get(i));
			Collections.reverse(byScore.get(i));
		}
		addAll(markers);

		// Issue warnings for markers that are outside the bands
		for (int side = 1; side >= -1; side -= 2)
			for (Iterator<SortedSet<Marker>> iter = (side == 1 ? proximal_markers
					: distal_markers).values().iterator(); iter.hasNext();) {
				SortedSet<Marker> outside_markers = iter.next();
				int band_index = Collections.binarySearch(bands,
						outside_markers.first());
				if (band_index < 0)
					band_index = -band_index - 2;
				Genome.Band band = bands.get(band_index + side);
				System.err.println(outside_marker_warning(outside_markers
						.first(), band, side));
				if (outside_markers.size() > 2)
					System.err.println("(" + (outside_markers.size() - 2)
							+ " messages suppressed)");
				if (outside_markers.size() > 1)
					System.err.println(outside_marker_warning(outside_markers
							.last(), band, side));
			}
	}

	static String outside_marker_warning(Marker marker, Genome.Band band,
			int side) {
		return "Warning: marker "
				+ marker
				+ " is "
				+ (side == 1 ? "proximal to the first" : "distal to the last")
				+ " band ("
				+ band.ID
				+ ") on chromosome "
				+ band.chromosome
				+ ", which "
				+ (side == 1 ? "starts at " + band.start : "ends at "
						+ band.end)
				+ ".  If this is not what you expect, check you are using the correct files and/or modify the band file located at "
				+ Distribution.bands_file + '.';
	}

	public Collection<Marker> closest_markers(GenomeLocation loc) {
		int start_index = Collections.binarySearch(this, loc);
		if (start_index < 0)
			start_index = -start_index - 1;
		while (start_index > 0 && get(start_index - 1).overlaps(loc))
			--start_index;
		int end_index = Collections.binarySearch(this, loc.end_location());
		if (end_index < 0)
			end_index = -end_index - 1;
		while (end_index < size() && get(end_index).overlaps(loc))
			++end_index;
		Collection<Marker> result = new Vector<Marker>();
		for (int i = start_index; i < end_index; i++)
			result.add(get(i));
		if (result.isEmpty()) {
			start_index = Collections.binarySearch(this, loc);
			if (start_index < 0)
				start_index = -start_index - 1;
			end_index = start_index;
			if (end_index < size())
				++end_index;
			// in case we are after the last marker
			if (start_index >= this.elementCount && this.elementCount > 0)
				start_index = this.elementCount - 1;
			for (int i = start_index; i < end_index; i++)
				result.add(get(i));
		}
		return result;
	}

	public Collection<Marker> closest_markers(Marker.aberration type,
			GenomeLocation loc) {
		int start_index = Collections.binarySearch(this, loc);
		if (start_index < 0)
			start_index = -start_index - 1;
		while (start_index > 0 && get(start_index - 1).overlaps(loc))
			--start_index;
		int end_index = Collections.binarySearch(this, loc.end_location());
		if (end_index < 0)
			end_index = -end_index - 1;
		while (end_index < size() && get(end_index).overlaps(loc))
			++end_index;
		Collection<Marker> result = new Vector<Marker>();
		for (int i = start_index; i < end_index; i++)
			if (get(i).Gscore[type.ordinal()] >= 0)
				result.add(get(i));
		if (result.isEmpty()) {
			while (end_index < size()
					&& get(end_index).Gscore[type.ordinal()] < 0)
				++end_index;
			int forward_distance = end_index == size() ? Integer.MAX_VALUE
					: get(end_index).distance(loc);
			while (start_index > 0
					&& get(start_index - 1).Gscore[type.ordinal()] < 0)
				--start_index;
			int backward_distance = start_index == 0 ? Integer.MAX_VALUE : get(
					start_index - 1).distance(loc);
			if (forward_distance != Integer.MAX_VALUE
					|| backward_distance != Integer.MAX_VALUE)
				if (forward_distance < backward_distance)
					result.add(get(end_index));
				else
					result.add(get(start_index - 1));
		}
		return result;
	}

	List<SignificantRegion> significantRegions(Distribution d, Genome genome,
			Genome miRNA_genome) {
		List<SignificantRegion> result = new Vector<SignificantRegion>();
		int end_chromosome = 0;
		while (end_chromosome < size()) {
			int start_chromosome = end_chromosome;
			end_chromosome = start_chromosome;
			while (end_chromosome < size()
					&& get(end_chromosome).chromosome == get(start_chromosome).chromosome)
				end_chromosome++;
			result.addAll(SignificantRegion.findRegions(d, subList(
					start_chromosome, end_chromosome), genome, miRNA_genome,
					(d.type == Marker.aberration.LOH ? loh_samples
							: copy_number_samples).length));
		}
		return result;
	}

	public void InMemory(int chromosome) throws IOException {
		sample_distributions = null;
		byScore = null;
		if (!Distribution.by_chromosome)
			return;
		clear();
		TreeSet<Marker> chr_markers = new TreeSet<Marker>();
		if (cpmarker != null && cpmarker.equals(lohmarker)) {
			cpmarker.add(lohmarker);
			lohmarker = null;
		}
		if (cpmarker != null && cpmarker.chromosome == chromosome)
			chr_markers.add(cpmarker);
		if (lohmarker != null && lohmarker.chromosome == chromosome)
			chr_markers.add(lohmarker);
		try {
			if (cpfile != null)
			{
				if(Distribution.chrom2run==-1)
				{
					for (;;) {
						String next = cpfile.readLine();
						if (next == null)
							break;
						cpmarker = new Marker(next.trim(), false,diffgroups);
						if (cpmarker.chromosome != chromosome)
							break;
						chr_markers.add(cpmarker);
					}					
				}else
				{
					for (;;) {
						String next = cpfile.readLine();
						if (next == null)
							break;
						cpmarker = new Marker(next.trim(), false,diffgroups);
						if (cpmarker.chromosome == chromosome)
							chr_markers.add(cpmarker);
					}	
				}				
			}
			if (lohfile != null)
				if(Distribution.chrom2run==-1)
				{
					for (;;) {
						String next = lohfile.readLine();
						if (next == null)
							break;
						lohmarker = new Marker(next.trim(), true,diffgroups);
						if (lohmarker.chromosome != chromosome)
							break;
						Marker old_marker = chr_markers.floor(lohmarker);
						if (lohmarker.equals(old_marker))
							old_marker.add(lohmarker);
						else
							chr_markers.add(lohmarker);
					}
				}else
				{
					for (;;) {
						String next = lohfile.readLine();
						if (next == null)
							break;
						lohmarker = new Marker(next.trim(), true, diffgroups);
						if (lohmarker.chromosome != chromosome)
							continue;
						Marker old_marker = chr_markers.floor(lohmarker);
						if (lohmarker.equals(old_marker))
							old_marker.add(lohmarker);
						else
							chr_markers.add(lohmarker);
					}
				}
		} catch (OutOfMemoryError e) {
			System.err.println("While reading " + chromosome
					+ " into memory, after reading " + chr_markers.size());
			throw e;
		}
		addAll(chr_markers);
		/*
		 * if(lohfile != null && cpfile != null) for(Marker.aberration type :
		 * Marker.aberration.values()) for(ListIterator<Marker> marker_iter =
		 * listIterator(); marker_iter.hasNext(); ){ Marker marker =
		 * marker_iter.next(); if(marker.Gscore[type.ordinal()] >= 0) continue;
		 * int marker_index = marker_iter.previousIndex(); int floor =
		 * marker_index, ceiling = marker_index; while(floor >= 0 &&
		 * get(floor).Gscore[type.ordinal()] < 0) floor--; while(ceiling <
		 * size() && get(ceiling).Gscore[type.ordinal()] < 0) ceiling++;
		 * marker.imputeNeighborScores(type, floor < 0 ? null : get(floor),
		 * ceiling >= size() ? null : get(ceiling)); }
		 */
	}

	public int[] diffgroups=null;
	public boolean differentialGscore=false;
	public String[] copy_number_samples = null;
	public String[] loh_samples = null;

	public List<List<Integer>> byScore;

	public List<List<TreeSet<Distribution.sparse_cell>>> sample_distributions;

	static public boolean focal_thresholds = false;

	static final int max_arm_markers = 150000;

	BufferedReader cpfile = null, lohfile = null;
	Marker cpmarker = null, lohmarker = null;
	public int cp_num_markers = 0, loh_num_markers = 0;
}
