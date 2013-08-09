package JISTIC;

import java.util.*;

public class SignificantRegion extends Region {

	static boolean process_one_aberrant_sample = false;

	List<Region> subregions;
	List<PeakRegion> peaks;

	SignificantRegion(Distribution distribution, List<Marker> markers,
			Genome genome, Genome miRNA_genome) {
		super(distribution.type, markers, genome, miRNA_genome);
		subregions = new Vector<Region>();
		int region_start = 0;
		while (region_start < markers.size()) {
			int region_end = region_start;
			while (region_end < markers.size()
					&& markers.get(region_end).Gscore[type.ordinal()] == markers
							.get(region_start).Gscore[type.ordinal()])
				region_end++;
			subregions.add(new Region(type, markers.subList(region_start,
					region_end), genome, miRNA_genome));
			region_start = region_end;
		}
		peaks = new Vector<PeakRegion>();
	}

	static List<SignificantRegion> findRegions(Distribution distribution,
			List<Marker> markers, Genome genome, Genome miRNA_genome,
			int num_samples) {
		Marker.aberration type = distribution.type;
		List<SignificantRegion> regions = new Vector<SignificantRegion>();
		List<Set<Integer>> peeled_off = new Vector<Set<Integer>>();
		for (int i = 0; i < markers.size(); i++)
			peeled_off.add(new HashSet<Integer>());
		int[] Gscores = new int[markers.size()];
		double[] qvalues = new double[markers.size()];
		int[] peel_off_Gscores = new int[markers.size()];
		double[] peel_off_qvalues = new double[markers.size()];
		double peak_qvalue = 1;
		for (Iterator<Marker> iter = markers.iterator(); iter.hasNext();)
			peak_qvalue = Math.min(peak_qvalue, iter.next().qvalue[type
					.ordinal()]);

		if (debug)
			System.err.println("finding regions for " + type.toString()
					+ ":chr" + markers.get(0).chromosome + " peak qvalue "
					+ peak_qvalue);

		for (ListIterator<Marker> iter = markers.listIterator(); iter.hasNext();) {
			Marker marker = iter.next();
			Gscores[iter.previousIndex()] = peel_off_Gscores[iter
					.previousIndex()] = marker.Gscore[type.ordinal()];
			qvalues[iter.previousIndex()] = peel_off_qvalues[iter
					.previousIndex()] = marker.qvalue[type.ordinal()];
		}

		peak_loop: while (peak_qvalue < Distribution.qvalthreshold) {
			int raw_peak_start = 0;
			while (peel_off_qvalues[raw_peak_start] > peak_qvalue)
				raw_peak_start++;
			int raw_peak_end = raw_peak_start;
			while (raw_peak_end < markers.size()
					&& peel_off_qvalues[raw_peak_end] <= peak_qvalue)
				raw_peak_end++;
			int peak_start = raw_peak_start, peak_end = raw_peak_end;
			for (int sample = 0; sample < num_samples; sample++)
				if (!peeled_off.get(peak_start).contains(sample)) {
					int one_out_peak_gscore = 0;
					for (int marker = peak_start; marker < peak_end; marker++)
						one_out_peak_gscore = Math.max(one_out_peak_gscore,
								peel_off_Gscores[marker]
										- markers.get(marker).GscoreContribution(type,sample));

					if (one_out_peak_gscore <= 0)
						if (process_one_aberrant_sample) {
							for (int marker = peak_start; marker < peak_end; marker++) {
								peel_off_Gscores[marker] -= markers.get(marker).GscoreContribution(type, sample);
								Map.Entry<Integer, Double> entry = distribution.gscore_to_qvalue
										.floorEntry(peel_off_Gscores[marker]);
								peel_off_qvalues[marker] = entry == null ? 1.
										: entry.getValue().doubleValue();
							}
							peak_qvalue = 1;
							for (int i = 0; i < markers.size(); i++)
								peak_qvalue = Math.min(peak_qvalue,
										peel_off_qvalues[i]);
							continue peak_loop;
						} else {
							Region peak = new Region(type, markers.subList(
									peak_start, peak_end), genome, miRNA_genome);
							System.err
									.println(peak.toString()
											+ " has only one aberrant sample: "
											+ distribution.data.copy_number_samples[sample]);
							System.err
									.println("This means this data is probably unsuitable for the JISTIC algorithm,");
							System.err
									.println("either because there are too few samples or because the samples have too few aberrations");
							System.err
									.println("If you are convinced you want to run JISTIC on this data,");
							System.err.println("please add the line");
							System.err.println();
							System.err.println("FixForSingleAberrantSample");
							System.err.println();
							System.err.println("to the spec file and re-run");
							System.exit(1);
						}

					while (peak_start > 0
							&& peel_off_Gscores[peak_start - 1]
									- markers.get(peak_start - 1)
											.GscoreContribution(type, sample) >= one_out_peak_gscore)
						peak_start--;
					while (peak_end < markers.size()
							&& peel_off_Gscores[peak_end]
									- markers.get(peak_end).GscoreContribution(
											type, sample) >= one_out_peak_gscore)
						peak_end++;
				}
			PeakRegion peak = new PeakRegion(type, markers.subList(peak_start,
					peak_end), genome, miRNA_genome, peak_qvalue);
			SignificantRegion region = null;
			int region_index = Collections.binarySearch(regions, peak);
			if (region_index < 0)
				region_index = -region_index - 2;
			if (region_index >= 0 && regions.get(region_index).overlaps(peak))
				region = regions.get(region_index);
			else {
				int start = peak_start, end = peak_end;
				while (start > 0
						&& !(qvalues[start - 1] >= Distribution.qvalthreshold))
					start--;
				while (Double.isNaN(qvalues[start]))
					start++;
				while (end < markers.size()
						&& !(qvalues[end] >= Distribution.qvalthreshold))
					end++;
				while (Double.isNaN(qvalues[end - 1]))
					end--;
				regions.add(region_index + 1, region = new SignificantRegion(
						distribution, markers.subList(start, end), genome,
						miRNA_genome));
			}
			if (debug) {
				System.err.println("found peak " + peak.toString()
						+ " in region " + region.toString()
						+ ", previously had " + region.peaks.size() + " peaks");
				System.err.print("before peel_off, peak Gscores are");
				for (int i = peak_start; i < peak_end; i++)
					System.err.print(" " + peel_off_Gscores[i]);
				System.err.println();
				System.err.print("before peel_off, peak qvalues are");
				for (int i = peak_start; i < peak_end; i++)
					System.err.print(" " + peel_off_qvalues[i]);
				System.err.println();
			}

			region.peaks.add(-Collections.binarySearch(region.peaks, peak) - 1,
					peak);
			if (region.isBroad()) {
				int peak_gscore = 0;
				for (int i = peak_start; i < peak_end; i++)
					peak_gscore = Math.max(peak_gscore, Gscores[i]);
				if (peak_gscore > distribution.Gsig * 2) {
					int broad_region_start = peak_start;
					while (broad_region_start > 0
							&& Gscores[broad_region_start - 1] > peak_gscore
									- distribution.Gsig)
						broad_region_start--;
					int broad_region_end = peak_end;
					while (broad_region_end < markers.size()
							&& Gscores[broad_region_end] > peak_gscore
									- distribution.Gsig)
						broad_region_end++;
					peak.has_focal = (markers.get(broad_region_end - 1).end - markers
							.get(broad_region_start).start) < broad_threshold[region.arm + 1][region.chromosome - 1];
				}
			}

			// limits for peeloff
			int max_peeloff[]=new int[2];
			if(limited_peeloff)
			{
				// aberration remaining from peak for each sample 
				int[] remain=new int[num_samples];
				// aberration added on top of remaining for each sample
				int[] added=new int[num_samples];
				int totalremain=0;
				int totaladded=0;
				
				// obtaining forward limit 
				ListIterator<Marker> iter = markers.listIterator(peak_end-1);
				max_peeloff[1]=iter.nextIndex();
				
				// initialize with peak end marker
				Marker marker=iter.next();
				for (int sample = 0; sample < num_samples; sample++)
				{
					remain[sample]=0;
					if (!peeled_off.get(max_peeloff[1]).contains(sample))
						remain[sample]=marker.GscoreContribution(type, sample);
					totalremain+=remain[sample];	
				}
				
				// iterate forward while under threshold
				while(iter.hasNext())
				{
					totalremain=0;
					totaladded=0;					
					max_peeloff[1]=iter.nextIndex();
					marker=iter.next();
					for (int sample = 0; sample < num_samples; sample++)
					{
						int markerCont=0;
						if (!peeled_off.get(max_peeloff[1]).contains(sample))
						{
							markerCont=marker.GscoreContribution(type, sample);												
							if(markerCont<remain[sample])
							{
								int end=Math.min(max_peeloff[1]+sparam, markers.size());
								for (int lfmid = max_peeloff[1]; lfmid < end; lfmid++) {
									markerCont=Math.min(Math.max(markerCont, markers.get(lfmid).GscoreContribution(type, sample)),remain[sample]);
								}
							}
						}
						
						added[sample]=Math.max(markerCont-remain[sample],0);
						remain[sample]=markerCont-added[sample];
						totalremain+=remain[sample];
						totaladded+=Math.max(markerCont-remain[sample], 0);						
					}
					
					// finished if cutoff or no peak remaining
					if(totalremain==0 || totaladded>=peeloff_gscore_thres[type.ordinal()])
					{
						break;
					}
					
				}
				
				
				// obtaining backward limit 
				iter = markers.listIterator(peak_start+1);
				
				// initialize with peak start marker
				max_peeloff[0]=iter.previousIndex();
				marker=iter.previous();
				totalremain=0;
				for (int sample = 0; sample < num_samples; sample++)
				{
					remain[sample]=0;
					if (!peeled_off.get(max_peeloff[0]).contains(sample))
						remain[sample]=marker.GscoreContribution(type, sample);
					totalremain+=remain[sample];
				}
				
				// iterate backwards while under threshold
				while(iter.hasPrevious())
				{
					totalremain=0;
					totaladded=0;		
					max_peeloff[0]=iter.previousIndex();
					marker=iter.previous();
					for (int sample = 0; sample < num_samples; sample++)
					{
						
						int markerCont=0;					
						if (!peeled_off.get(max_peeloff[0]).contains(sample))
						{
							markerCont=marker.GscoreContribution(type, sample);												
							if(markerCont<remain[sample])
							{
								int end=Math.max(max_peeloff[0]-sparam, 0);
								for (int lfmid = max_peeloff[0]; lfmid >= end; lfmid--) {
									markerCont=Math.min(Math.max(markerCont, markers.get(lfmid).GscoreContribution(type, sample)),remain[sample]);
								}
								/*
								int end=Math.min(max_peeloff[0]+sparam, markers.size());
								for (int lfmid = max_peeloff[0]; lfmid < end; lfmid++) {
									markerCont=Math.min(Math.max(markerCont, markers.get(lfmid).GscoreContribution(type, sample)),remain[sample]);
								}
								*/
							}
						}
						added[sample]=Math.max(markerCont-remain[sample],0);
						remain[sample]=markerCont-added[sample];
						totalremain+=remain[sample];
						totaladded+=Math.max(markerCont-remain[sample], 0);						
					}														
					
					// finished if cutoff or no peak remaining
					if(totalremain==0 || totaladded>=peeloff_gscore_thres[type.ordinal()])
					{
						break;
					}
					
				}																
			}
			
			peel_off_loop: for (int sample = 0; sample < num_samples; sample++)
				for (ListIterator<Marker> iter = markers
						.listIterator(peak_start); iter.nextIndex() < peak_end;)
					if (iter.next().aberrationStrength(type, sample) > 0
							&& !peeled_off.get(iter.previousIndex()).contains(
									sample)) {
						int peeloff_limits[];
						
						if (arm_peeloff)
							peeloff_limits = genome.ArmLimits(iter);
						else {
							peeloff_limits = new int[2];
							while (iter.hasPrevious()
									&& iter.previous().aberrationStrength(type,sample) > 0)											
								;
							if (iter.hasNext()
									&& iter.next().aberrationStrength(type,
											sample) > 0)
								iter.previous();
							peeloff_limits[0] = iter.nextIndex();
							
							while (iter.hasNext()
									&& iter.next().aberrationStrength(type,
											sample) > 0)
								;
							if (iter.hasPrevious()
									&& iter.previous().aberrationStrength(type,
											sample) > 0)
								iter.next();
							peeloff_limits[1] = iter.previousIndex();
							
							if(limited_peeloff)
							{
								peeloff_limits[0] = Math.max(peeloff_limits[0],max_peeloff[0]);
								peeloff_limits[1] = Math.min(peeloff_limits[1],max_peeloff[1]);								
							}
							
							
						}
						if (debug)
							System.err
									.println("for sample "
											+ sample
											+ " removed "
											+ Arrays.toString(peeloff_limits)
											+ " : "
											+ markers.get(peeloff_limits[0])
													.toString()
											+ " thru "
											+ markers.get(peeloff_limits[1])
													.toString());
						for (int marker = peeloff_limits[0]; marker <= peeloff_limits[1]; marker++)
							if (peeled_off.get(marker).add(sample))
								peel_off_Gscores[marker] -= markers.get(marker)
										.GscoreContribution(type, sample);
						continue peel_off_loop;
					}

			peak_qvalue = 1;
			for (int i = 0; i < markers.size(); i++) {
				Map.Entry<Integer, Double> entry = distribution.gscore_to_qvalue
						.floorEntry(peel_off_Gscores[i]);
				peel_off_qvalues[i] = entry == null ? 1. : entry.getValue()
						.doubleValue();
				peak_qvalue = Math.min(peak_qvalue, peel_off_qvalues[i]);
			}

			if (debug) {
				System.err.print("after peel_off, peak Gscores are");
				for (int i = peak_start; i < peak_end; i++)
					System.err.print(" " + peel_off_Gscores[i]);
				System.err.println();
				System.err.print("after peel_off, peak qvalues are");
				for (int i = peak_start; i < peak_end; i++)
					System.err.print(" " + peel_off_qvalues[i]);
				System.err.print(" new peak " + peak_qvalue);
				System.err.println();
			}
			
		}
		return regions;
	} 

	static int sparam = 11;
	static boolean arm_peeloff = false;	
	static boolean limited_peeloff = false;
	static int[] peeloff_gscore_thres={0,0,0};
	static float peeloff_qval_thres = Float.NaN;


	static int[][] broad_threshold = null;
	/*{
			{ 50000000, 50000000, 40000000, 40000000, 35000000, 35000000, 30000000, 30000000, 30000000, 30000000, 30000000, 25000000, 25000000, 20000000, 20000000, 20000000,
					15000000, 15000000, 12000000, 12000000, 10000000, 10000000, 30000000, 12000000 },
			{ 50000000, 50000000, 40000000, 40000000, 35000000, 35000000, 30000000, 30000000, 30000000, 30000000, 30000000, 25000000, 25000000, 20000000, 20000000, 20000000,
					15000000, 15000000, 12000000, 12000000, 10000000, 10000000, 30000000, 12000000 },
			{ 50000000, 50000000, 40000000, 40000000, 35000000, 35000000, 30000000, 30000000, 30000000, 30000000, 30000000, 25000000, 25000000, 20000000, 20000000, 20000000,
					15000000, 15000000, 12000000, 12000000, 10000000, 10000000, 30000000, 12000000 } };*/

	boolean isBroad() {
		return (end - start)>= broad_threshold[arm + 1][chromosome - 1];
	}

	public String toString() {
		return (isBroad() ? "BROAD:" : "") + super.toString();
	}
}
