package JISTIC;

import java.util.*;

public class Gene extends GenomeLocation implements Comparator<Region> {
	public int[] exonstart = null;
	public int[] exonend = null;
	public String ID;
	public String symbol;
	double[] qvalue = null;
	int[] Gscore = null;

	public static int translateChr(String s) {
		if (s.charAt(0) == 'X')
			return 23;
		else if (s.charAt(0) == 'Y')
			return 24;
		else
			return Math.round(Float.parseFloat(s));
	}

	public Gene(String geneline) {
		super(translateChr(geneline.split("\\t")[3]), Integer.parseInt(geneline
				.split("\\t")[5]), Integer.parseInt(geneline.split("\\t")[6]));
		ID = geneline.split("\\t")[0];
		symbol = geneline.split("\\t")[1];
	}

	Collection<Marker> findAberrations(SNPdata data) {
		if (qvalue != null)
			return null;
		Collection<Marker> closest_markers = data.closest_markers(this);		
		if(closest_markers.size()!=0)
		{
			qvalue = new double[Marker.aberrationtypes.length];
			Gscore = new int[Marker.aberrationtypes.length];
			for (int i = 0; i < Marker.aberrationtypes.length; i++) {
				Marker.aberration ab = Marker.aberrationtypes[i];
				// closest_markers.add(data.closest_markers(ab, this));
				// qvalue[ab.ordinal()] =
				// closest_markers.get(ab.ordinal()).isEmpty() ? 1 : 0;
				for (Iterator<Marker> iter = closest_markers.iterator(); iter
						.hasNext();) {
					Marker m = iter.next();
					Gscore[ab.ordinal()] += m.Gscore[ab.ordinal()];
					// qvalue[ab.ordinal()] += m.qvalue[ab.ordinal()] /closest_markers.get(ab.ordinal()).size();
				}
				Gscore[ab.ordinal()] = Math.round(((float)Gscore[ab.ordinal()])/ closest_markers.size());
				Map.Entry<Integer, Double> entry = Distribution.distributions[ab
						.ordinal()].gscore_to_qvalue
						.floorEntry(Gscore[ab.ordinal()]);
				qvalue[ab.ordinal()] = entry == null ? 1. : entry.getValue()
						.doubleValue();
			}
		}
		return closest_markers;
	}

	boolean Aberration(Marker.aberration type) {
		return qvalue[type.ordinal()] < Distribution.qvalthreshold;
	}

	boolean hasAberrations() {
		for (int i = 0; i < Marker.aberrationtypes.length; i++) {
			Marker.aberration ab = Marker.aberrationtypes[i];
			if (Aberration(ab))
				return true;
		}
		return false;
	}

	boolean isInRegion(Region r) {
		return r.genes.contains(this) || r.miRNA.contains(this);
	}

	public int compare(Region r1, Region r2) {
		if (overlaps(r1) && !overlaps(r2))
			return -1;
		if (overlaps(r2) && !overlaps(r1))
			return 1;
		if (r1.type != Marker.aberration.LOH
				&& r2.type == Marker.aberration.LOH)
			return -1;
		if (r1.type == Marker.aberration.LOH
				&& r2.type != Marker.aberration.LOH)
			return 1;
		if ((r1 instanceof PeakRegion) && !(r2 instanceof PeakRegion))
			return -1;
		if (!(r1 instanceof PeakRegion) && (r2 instanceof PeakRegion))
			return 1;
		int result = overlaps(r1) ? inCommon(r2) - inCommon(r1) : distance(r1)
				- distance(r2);
		if (result == 0 && r1.type != r2.type)
			result = r1.type == Marker.aberration.AMP ? -1 : 1;
		return result == 0 ? r1.start = r2.start : result;
	}

	public boolean equals(Object o) {
		return (o instanceof Gene) && same_location((Gene) o)
				&& ID.equals(((Gene) o).ID);
	}

	// Override comparisons so two genes with the same location but different
	// names are not considered equal
	public int compareTo(GenomeLocation o) {
		int result = super.compareTo(o);
		return result == 0 && (o instanceof Gene) ? ID.compareTo(((Gene) o).ID)
				: result;
	}

	public String toString() {
		return ID + ':' + super.toString();
	}
}
