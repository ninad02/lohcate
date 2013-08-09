package JISTIC;

import java.util.Arrays;
import java.util.Comparator;
import java.util.Map;

public class Marker extends GenomeLocation {
	public String name;
	public float[] copy_number = null;
	public float[] copy_number_uncapped = null;
	public float[] isLOH = null;
	public int Gscore[];
	public double qvalue[];

	public static enum aberration {
		AMP, DEL, LOH
	};

	public static aberration[] aberrationtypes = null;

	public static int initial_columns = 3;

	public static class GscoreOrdering implements Comparator<Marker> {
		Marker.aberration type;

		GscoreOrdering(Marker.aberration t) {
			type = t;
		}

		public int compare(Marker m1, Marker m2) {
			return m2.Gscore[type.ordinal()] - m1.Gscore[type.ordinal()];
		}
	}

	public static float[] thresholds = { (float) .3, (float) -.3 };
	public static float[] highlevelthresholds = { (float) .9, (float) -1.3 };
	public static float[] caps = { Float.POSITIVE_INFINITY,
			Float.NEGATIVE_INFINITY };
	public static float[][] per_sample_thresholds = new float[2][];
	
	public static int numchrom=24;

	public static GscoreOrdering[] order = {
			new GscoreOrdering(aberration.AMP),
			new GscoreOrdering(aberration.DEL),
			new GscoreOrdering(aberration.LOH) };

	public Marker(String data, boolean islohdata, int[] weights) {
		super((int) Math.round(Double.parseDouble(data.split("\\t")[1])),
				(int) Math.round(Double.parseDouble(data.split("\\t")[2])),
				(int) Math.round(Double.parseDouble(data.split("\\t")[2])));
		String[] tokens = data.split("\\t");
		name = tokens[0];
		Gscore = new int[aberrationtypes.length];
		qvalue = new double[aberrationtypes.length];
		Arrays.fill(Gscore, -1);
		Arrays.fill(qvalue, Double.NaN);
		int amp_sum = 0, del_sum = 0;
		if (tokens.length > initial_columns)
			// if(Character.isLetter(tokens[initial_columns].charAt(0))){
			if (islohdata) {
				isLOH = new float[tokens.length - initial_columns];
				int numLOH = 0;
				Arrays.fill(isLOH, 0);
				for (int i = 0; i < isLOH.length; i++) {
					isLOH[i] = Float.parseFloat(tokens[i + initial_columns]);
					numLOH += isLOH[i];
				}
				/*
				 * for(int i=0; i < isLOH.length; i++)
				 * if(tokens[i+initial_columns].charAt(0) == 'L'){ isLOH[i] = 1;
				 * numLOH++; }
				 */
				Gscore[aberration.LOH.ordinal()] = Math.round(numLOH
						/ Distribution.binsize);
			} else {
				copy_number = new float[tokens.length - initial_columns];
				copy_number_uncapped = new float[tokens.length
						- initial_columns];
				for (int i = 0; i < copy_number.length; i++) {
					float val = Float.parseFloat(tokens[i + initial_columns]);
					if (Distribution.mixed_aberrations && val < 0)
						val *= -1;
					copy_number_uncapped[i] = val;
					copy_number[i] = Math.max(caps[aberration.DEL.ordinal()],
							Math.min(caps[aberration.AMP.ordinal()], val));
					if (copy_number[i] > per_sample_thresholds[0][i]
							+ thresholds[aberration.AMP.ordinal()])
					{
						amp_sum += weights[i]*Math.round(copy_number[i]/ Distribution.binsize);						
					}
					else if (copy_number[i] < per_sample_thresholds[1][i]
							+ thresholds[aberration.DEL.ordinal()])
						del_sum -= weights[i]*Math.round(copy_number[i]/ Distribution.binsize);
				}
				Gscore[aberration.AMP.ordinal()] = amp_sum;
				Gscore[aberration.DEL.ordinal()] = del_sum;
			}
		for (int i = 0; i < Marker.aberrationtypes.length; i++) {
			Marker.aberration type = Marker.aberrationtypes[i];
			if ((Gscore[type.ordinal()] >= 0 || Distribution.differentialfile!=null)
					&& Distribution.distributions != null) {
				Map.Entry<Integer, Double> entry = Distribution.distributions[type
						.ordinal()] == null ? null
						: Distribution.distributions[type.ordinal()].gscore_to_qvalue
								.floorEntry(Gscore[type.ordinal()]);
				qvalue[type.ordinal()] = entry == null ? 1. : entry.getValue()
						.doubleValue();
			}
		}
		if (Distribution.tumorfile == null
				&& Distribution.distributions != null && !Distribution.keep_all
				&& Gscore[0] == 0 && Gscore[1] == 0 && Gscore[2] == 0) {
			copy_number = null;
			isLOH = null;
		}
		/*
		 * if(qvalue[2] > Distribution.qvalthreshold) isLOH = null;
		 */
	}

	public boolean equals(Object o) {
		return (o instanceof Marker) && ((Marker) o).name.equals(name);
	}

	public int compareTo(GenomeLocation o) {
		int result = super.compareTo(o);
		return result == 0 && (o instanceof Marker) ? name
				.compareTo(((Marker) o).name) : result;
	}

	public void add(Marker o) {
		if (o.copy_number != null) {
			copy_number = o.copy_number;
			Gscore[aberration.AMP.ordinal()] = o.Gscore[aberration.AMP
					.ordinal()];
			Gscore[aberration.DEL.ordinal()] = o.Gscore[aberration.DEL
					.ordinal()];
		}
		if (o.isLOH != null) {
			isLOH = o.isLOH;
			Gscore[aberration.LOH.ordinal()] = o.Gscore[aberration.LOH
					.ordinal()];
			qvalue[aberration.LOH.ordinal()] = o.qvalue[aberration.LOH
					.ordinal()];
		}
	}

	public String toString() {
		String result = name + '\t' + super.toString();
		for (int i = 0; i < qvalue.length; i++) {
			result += '\t';
			if (!Double.isNaN(qvalue[i]))
				result += qvalue[i];

			result += '\t';
			if (!Double.isNaN(qvalue[i]))
				result += Gscore[i];

		}
		return result;
	}

	public void imputeNeighborScores(aberration type, Marker n1, Marker n2) {
		if (distance(n1) == Integer.MAX_VALUE) {
			Gscore[type.ordinal()] = n2.Gscore[type.ordinal()];
			qvalue[type.ordinal()] = n2.qvalue[type.ordinal()];
			if (type == aberration.LOH)
				isLOH = n2.isLOH;
			else
				copy_number = n2.copy_number;
		} else if (distance(n2) == Integer.MAX_VALUE) {
			Gscore[type.ordinal()] = n1.Gscore[type.ordinal()];
			qvalue[type.ordinal()] = n1.qvalue[type.ordinal()];
			if (type == aberration.LOH)
				isLOH = n1.isLOH;
			else
				copy_number = n1.copy_number;
		} else {
			int dist1 = distance(n1);
			int dist2 = distance(n2);
			if (dist1 == 0 && dist2 == 0)
				dist1 = dist2 = 1;
			Gscore[type.ordinal()] = (dist2 * n1.Gscore[type.ordinal()] + dist1
					* n2.Gscore[type.ordinal()])
					/ (dist1 + dist2);
			Map.Entry<Integer, Double> entry = Distribution.distributions[type
					.ordinal()].gscore_to_qvalue.floorEntry(Gscore[type
					.ordinal()]);
			qvalue[type.ordinal()] = entry == null ? 1. : entry.getValue()
					.doubleValue();
			if (type == aberration.LOH) {
				isLOH = new float[n1.isLOH.length];
				for (int i = 0; i < isLOH.length; i++)
					isLOH[i] = Math.round((dist2 * n1.isLOH[i] + dist1
							* n2.isLOH[i])
							/ (dist1 + dist2));
			} else {
				copy_number = new float[n1.copy_number.length];
				for (int i = 0; i < copy_number.length; i++)
					copy_number[i] = (dist2 * n1.copy_number[i] + dist1
							* n2.copy_number[i])
							/ (dist1 + dist2);
			}
		}
	}

	int aberrationStrength(aberration type, int sample) {
		switch (type) {
		case AMP:
			if (copy_number != null)
				if (copy_number[sample] > per_sample_thresholds[0][sample]
						+ highlevelthresholds[type.ordinal()])
					return 2;
				else if (copy_number[sample] > per_sample_thresholds[0][sample]
						+ thresholds[type.ordinal()])
					return 1;
			break;

		case DEL:
			if (copy_number != null)
				if (copy_number[sample] < per_sample_thresholds[1][sample]
						+ highlevelthresholds[type.ordinal()])
					return 2;
				else if (copy_number[sample] < per_sample_thresholds[1][sample]
						+ thresholds[type.ordinal()])
					return 1;
			break;

		case LOH:
			if (isLOH != null && isLOH[sample] > 0)
				return 1;
			break;
		}
		return 0;
	}

	int GscoreContribution(aberration type, int sample) {
		return aberrationStrength(type, sample) == 0 ? 0 : Math.abs(Math
				.round((type == aberration.LOH ? isLOH : copy_number)[sample]
						/ Distribution.binsize));
	}

	
	/*
	int GscoreContributionCeil(aberration type, int sample) {
		return aberrationStrength(type, sample) == 0 ? 0 : (int) Math.ceil(Math
				.abs((type == aberration.LOH ? isLOH : copy_number)[sample]
						/ Distribution.binsize));
	}*/

	float copy_number_val(int sample) {
		return copy_number == null ? 0 : copy_number[sample];
	}

	float copy_number_val_uncapped(int sample) {
		return copy_number == null ? 0 : copy_number_uncapped[sample];
	}
}
