package JISTIC;

import java.util.*;

public class Region extends GenomeLocation {

	Marker.aberration type;
	List<Marker> markers;
	List<Gene> genes;
	List<Gene> miRNA;
	double max_qvalue;
	double min_qvalue;
	boolean loh_only;
	int arm = 0;
	Genome.Band[] bands = null;
	static boolean debug = false;

	public static int gene_distance_max = 100000;

	Region(Marker.aberration type, List<Marker> markers, Genome genome,
			Genome miRNA_genome) {
		super(markers.get(0).chromosome, markers.get(0).start, markers
				.get(markers.size() - 1).end);
		this.type = type;
		this.markers = markers;
		max_qvalue = 0;
		min_qvalue = 1;
		loh_only = type == Marker.aberration.LOH;
		for (Iterator<Marker> iter = markers.iterator(); iter.hasNext();) {
			Marker marker = iter.next();
			max_qvalue = Math.max(max_qvalue, marker.qvalue[type.ordinal()]);
			min_qvalue = Math.min(min_qvalue, marker.qvalue[type.ordinal()]);
			if (marker.qvalue[Marker.aberration.AMP.ordinal()] < Distribution.qvalthreshold
					|| marker.qvalue[Marker.aberration.DEL.ordinal()] < Distribution.qvalthreshold)
				loh_only = false;
		}

		genes = get_associated_genes(genome);
		miRNA = get_associated_genes(miRNA_genome);

		if (genome.bands != null) {
			int starting_band = Collections.binarySearch(genome.bands, markers
					.get(0));
			if (starting_band < 0)
				starting_band = -starting_band - 1;
			while (starting_band > 0
					&& genome.bands.get(starting_band - 1).overlaps(
							markers.get(0)))
				starting_band--;
			int ending_band = Collections.binarySearch(genome.bands, markers
					.get(markers.size() - 1));
			if (ending_band < 0)
				ending_band = -ending_band - 1;
			bands = new Genome.Band[2];
			bands[0] = genome.bands.get(starting_band);
			bands[1] = ending_band == starting_band + 1 ? null : genome.bands
					.get(ending_band - 1);
			arm = bands[0].Arm();
			if (bands[1] != null && bands[1].Arm() != arm)
				arm = -1;
		}
	}

	private List<Gene> get_associated_genes(Genome genome) {
		if (genome == null)
			return new Vector<Gene>();

		int genome_start_index = Collections.binarySearch(genome, markers
				.get(0));
		if (genome_start_index < 0)
			genome_start_index = -genome_start_index - 1;
		if (debug)
			System.out.println(toString() + " start gene "
					+ genome.get(genome_start_index).toString() + " after "
					+ genome.get(genome_start_index - 1).toString());
		while (genome_start_index > 0
				&& genome.get(genome_start_index - 1).overlaps(markers.get(0))) {
			genome_start_index--;
			if (debug)
				System.out.println("retreat to "
						+ genome.get(genome_start_index).toString() + " after "
						+ genome.get(genome_start_index - 1).toString());
		}
		int genome_end_index = Collections.binarySearch(genome, markers
				.get(markers.size() - 1));
		if (genome_end_index < 0)
			genome_end_index = -genome_end_index - 1;
		if (debug)
			System.out.println(toString() + " end gene "
					+ genome.get(genome_end_index).toString() + " after "
					+ genome.get(genome_end_index - 1).toString());
		while (genome_end_index < genome.size()
				&& genome.get(genome_end_index).overlaps(
						markers.get(markers.size() - 1))) {
			genome_end_index++;
			if (debug)
				System.out.println("forward to "
						+ genome.get(genome_end_index).toString() + " after "
						+ genome.get(genome_end_index - 1).toString());
		}

		int forward_distance = genome_end_index == genome.size() ? Integer.MAX_VALUE
				: genome.get(genome_end_index).distance(this);
		int backward_distance = genome_start_index == 0 ? Integer.MAX_VALUE
				: genome.get(genome_start_index - 1).distance(this);
		if (forward_distance < gene_distance_max)
			do
				genome_end_index++;
			while (genome_end_index < genome.size()
					&& genome.get(genome_end_index).distance(this) <= forward_distance);

		if (backward_distance < gene_distance_max)
			do
				genome_start_index--;
			while (genome_start_index >= 0
					&& genome.get(genome_start_index).distance(this) <= backward_distance);

		List<Gene> result = genome.subList(genome_start_index, genome_end_index);	

		while (genome_start_index > 0) {
			int dist = genome.get(--genome_start_index).distance(this);
			if (dist == Integer.MAX_VALUE)
				break;
			else if (dist == 0) {
				result = new Vector<Gene>(result);
				result.add(genome.get(genome_start_index));
			}
		}

		while (genome_end_index < genome.size()) {
			int dist = genome.get(genome_end_index++).distance(this);
			if (dist == Integer.MAX_VALUE)
				break;
			else if (dist == 0) {
				result = new Vector<Gene>(result);
				result.add(genome.get(genome_end_index - 1));
			}
		}
		return result;
	}

	public String toString() {
		String result = (loh_only ? "LOH" : Distribution.name(type)) + ':'
				+ super.toString();
		if (bands != null) {
			result += "(" + bands[0].toString();
			if (bands[1] != null)
				result += "-" + bands[1].toString();
			result += ')';
		}
		return result;
	}
}
