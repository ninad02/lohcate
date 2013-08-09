package JISTIC;
import java.util.*;
import java.io.*;

public class Distribution {
	static class sparse_cell implements Comparable<sparse_cell> {
		int gscore;
		int count;
		double logcount;

				

	    public static sparse_cell getSingletonObject(int gscore)
	    {
	      if (singletonref == null)
	          singletonref = new sparse_cell(gscore);
	      singletonref.gscore=gscore;
	      singletonref.count=1;
	      singletonref.logcount=0;
	      return singletonref;
	    }
	    
	    public static sparse_cell getSingletonObject(sparse_cell c1, sparse_cell c2)
	    {
	      if (singletonref == null)
	          // it's ok, we can call this constructor
	    	  singletonref = new sparse_cell(c1.gscore + c2.gscore);
	      singletonref.gscore=c1.gscore + c2.gscore;
	      singletonref.count=-1;
	      singletonref.logcount=c1.logcount + c2.logcount;
	      return singletonref;
	    }

	    private static sparse_cell singletonref;

	    
		sparse_cell(int i)
		{ gscore = i; count = 1; logcount = Double.NaN; }

		sparse_cell(sparse_cell c1, sparse_cell c2)
		{ gscore = c1.gscore + c2.gscore; count = -1; logcount = c1.logcount + c2.logcount; }

		sparse_cell(sparse_cell c1)
		{ gscore = c1.gscore; count = c1.count; logcount = c1.logcount; }
		
		
		void add(sparse_cell o)
		{
			if(count > 0)
				count += o.count;
			else
			{
				logcount = addLogs(logcount, o.logcount);
			}
		}

		static double addLogs(double a, double b)
		{
			double larger = Math.max(a, b);
			double smaller = Math.min(a, b);
			return Math.log1p(Math.exp(smaller - larger)) + larger;
		}

		public boolean equals(Object o)
		{
			return (o instanceof sparse_cell) && ((sparse_cell)o).gscore == gscore;
		}

		public int compareTo(sparse_cell o)
		{
			return o.gscore - gscore;
		}

		public String toString()
		{ return String.valueOf(gscore) + ':' + Math.round(logcount); }

		public void computeLog()
		{ logcount = Math.log(count); }
	}

	List< TreeSet<sparse_cell> > sample_distributions;

	TreeSet<sparse_cell> gscore_distribution;

	TreeMap<Integer, Double> gscore_to_qvalue = null;

	public static float binsize = (float).001;

	public static float qvalthreshold = (float).25;

	public static boolean by_chromosome = false;

	public static Distribution[] distributions = null;

	public static boolean skip_all_normal = false;

	public static boolean mixed_aberrations = false;
	
	public static int chrom2run=-1; 
		

	double log_total_perms;

	int Gsig = 0;

	int num_markers;

	Marker.aberration type;

	SNPdata data;

	static TreeSet<sparse_cell> convolution(TreeSet<sparse_cell> d1, TreeSet<sparse_cell> d2)
	{
		TreeSet<sparse_cell> result = new TreeSet<sparse_cell>();
		for(Iterator<sparse_cell> iter1 = d1.iterator(); iter1.hasNext(); )
		{
			sparse_cell c1 = iter1.next();
			for(Iterator<sparse_cell> iter2 = d2.iterator(); iter2.hasNext(); )
			{
				// obtain cell with sum of gscores and product of counts
				sparse_cell new_cell = sparse_cell.getSingletonObject(c1, iter2.next());
				// if cell exists add counts, if not add new cell
				sparse_cell old_cell = result.floor(new_cell);
				if(new_cell.equals(old_cell))
					old_cell.add(new_cell);
				else
					result.add(new sparse_cell(new_cell));
			}
		}
		return result;
	}

	public TreeSet<sparse_cell> convolution(List<TreeSet<sparse_cell>> dataconv)
	{
		TreeSet<sparse_cell> result=dataconv.get(0);
		for(ListIterator< TreeSet<sparse_cell> > iter = dataconv.listIterator(1); iter.hasNext(); )
		{
			TreeSet<sparse_cell> nextitem=iter.next();
			result = convolution(result, nextitem);
		}
		return result;
	}
	
	public Distribution(SNPdata data, Marker.aberration type)
	{
		this.type = type;
		this.data = data;
		sample_distributions = data.sample_distributions.get(type.ordinal());
		int num_samples = sample_distributions.size();
		num_markers = type == Marker.aberration.LOH ? data.loh_num_markers : data.cp_num_markers;
		// compute log for each cell
		for(Iterator<TreeSet<sparse_cell> > distr_iter = sample_distributions.iterator(); distr_iter.hasNext(); )
			for(Iterator<sparse_cell> iter = distr_iter.next().iterator(); iter.hasNext(); )
				iter.next().computeLog();
		// convolution across samples
		if(data.differentialGscore)
		{
			ArrayList<TreeSet<sparse_cell>> dists1=new ArrayList<TreeSet<sparse_cell>>();
			ArrayList<TreeSet<sparse_cell>> dists2=new ArrayList<TreeSet<sparse_cell>>();
			for (int i = 0; i < data.diffgroups.length; i++) {
				if(data.diffgroups[i]>0)
				{
					dists1.add(sample_distributions.get(i));
				}else
				{
					dists2.add(sample_distributions.get(i));
				}
			}
			TreeSet<sparse_cell> gscore_dist1=convolution(dists1);
			TreeSet<sparse_cell> gscore_dist2=convolution(dists2);
			gscore_distribution=convolution(gscore_dist1,gscore_dist2);			
		}else
		{
			gscore_distribution=convolution(sample_distributions);
		}
		Iterator<sparse_cell> cell_iter = gscore_distribution.iterator();
		sparse_cell previous_cell = cell_iter.next();
		while(cell_iter.hasNext()){
			sparse_cell next_cell = cell_iter.next();
			next_cell.add(previous_cell);
			previous_cell = next_cell;
		}
		log_total_perms = Math.log(num_markers)*num_samples;	
	}

	void calculateSignificance(SNPdata data) {
		gscore_to_qvalue = new TreeMap<Integer, Double>();
		List<Integer> byScore = data.byScore.get(type.ordinal());
		ListIterator<Integer> iter = byScore.listIterator();

		int significant_markers = 0;
		/* FELIX DANGER, I JUST CHANGED THIS !!!
		while(significant_markers < byScore.size() &&
				byScore.get(significant_markers).intValue() >= 0)
				*/
		while(significant_markers < byScore.size())
		{
			sparse_cell gscore_cell = gscore_distribution.floor(sparse_cell.getSingletonObject(iter.next()));			
			double pval = Math.exp(gscore_cell.logcount - log_total_perms);
			while(significant_markers < byScore.size() &&
					iter.previous().equals(byScore.get(significant_markers)))
			{
				significant_markers++;
				iter.next();
			}
			if(significant_markers == byScore.size())
				iter.previous();
			double qval = (pval * num_markers) / significant_markers;
			Integer m = null;
			while(iter.nextIndex() < significant_markers)
				m = iter.next();
			gscore_to_qvalue.put(m.intValue(), qval);
		}

		double previous_qval = 1;
		Map.Entry<Integer, Double> entry;
		for(Iterator<Map.Entry<Integer, Double> > map_iter = gscore_to_qvalue.entrySet().iterator(); map_iter.hasNext(); )
			if((entry = map_iter.next()).getValue().doubleValue() > previous_qval)
				entry.setValue(previous_qval);
			else if((previous_qval = entry.getValue().doubleValue()) < qvalthreshold && Gsig == 0)
				Gsig = entry.getKey().intValue();

		if(!by_chromosome)
			for(ListIterator<Marker> marker_iter = data.listIterator(); marker_iter.hasNext(); ){
				Marker marker = marker_iter.next();
				if(marker.Gscore[type.ordinal()] >= 0 || data.differentialGscore){
					entry = Distribution.distributions[type.ordinal()].gscore_to_qvalue.floorEntry(marker.Gscore[type.ordinal()]);
					marker.qvalue[type.ordinal()] = entry==null ? 1. : entry.getValue().doubleValue();
				}
				/*
	    if(marker.Gscore[type.ordinal()] >= 0)
		continue;
	    int marker_index = marker_iter.previousIndex();
	    int floor = marker_index, ceiling = marker_index;
	    while(floor >= 0 && data.get(floor).Gscore[type.ordinal()] < 0)
		floor--;
	    while(ceiling < data.size() && data.get(ceiling).Gscore[type.ordinal()] < 0)
		ceiling++;
	    marker.imputeNeighborScores(type, floor < 0 ? null : data.get(floor),
	    				ceiling >= data.size() ? null : data.get(ceiling));
				 */
			}
	}

	static String cpfile = null, lohfile = null, genomefile = null, bands_file = null, tumorfile = null,differentialfile=null, spec_file = null, dir = "", miRNAfile = null,
	distribution_output = "GScoreQValueDistribution";


	static boolean keep_all = true;

	public static String name(Marker.aberration ab) {
		return mixed_aberrations && ab.equals(Marker.aberration.AMP)
		? "ABER" : ab.name();
	}

	public static void main(String[] args)
	{
		try {
			
			
			
			for(int i=0; i < args.length; i++)
				if(args[i].toLowerCase().startsWith("copynumber="))
					cpfile = args[i].substring(11);
				else if(args[i].toLowerCase().startsWith("loh="))
					lohfile = args[i].substring(4);
				else if(args[i].toLowerCase().startsWith("locations="))
					genomefile = args[i].substring(10);
				else if(args[i].toLowerCase().startsWith("genes="))
					genomefile = args[i].substring(6);
				else if(args[i].toLowerCase().startsWith("mirna="))
					miRNAfile = args[i].substring(6);
				else if(args[i].toLowerCase().startsWith("bands="))
					bands_file = args[i].substring(6);
				else if(args[i].toLowerCase().startsWith("tumors="))
					tumorfile = args[i].substring(7);
				else if(args[i].toLowerCase().startsWith("differential="))
					differentialfile = args[i].substring(13);			
				else if (args[i].toLowerCase().startsWith("numchrom=")) 
					Marker.numchrom = Integer.parseInt(args[i].substring(9));				
				else if(args[i].toLowerCase().startsWith("singlechrom="))
				{
					chrom2run = Integer.parseInt(args[i].substring(12));
					by_chromosome = true;
				}else if(args[i].toLowerCase().startsWith("spec="))
					spec_file = args[i].substring(5);
				else{
					System.err.println("unidentified argument " + args[i]);
					System.exit(1);
				}
			if(cpfile == null || genomefile == null){
				System.err.println("usage: JISTIC.Distribution [copynumber=<filename>] [loh=<filename>] locations=<filename>");
				System.exit(1);
			}
			if(spec_file != null){
				File sp = new File(spec_file);
				BufferedReader infile = new BufferedReader(new FileReader(sp));
				for(String next = infile.readLine(); next != null; next = infile.readLine()){
					String[] tokens = next.split("\\s+");
					if(tokens[0].equalsIgnoreCase("histogramBinsize"))
						binsize = Float.parseFloat(tokens[1]);
					else if(tokens[0].equalsIgnoreCase("qvalthreshold"))
						qvalthreshold = Float.parseFloat(tokens[1]);
					else if(tokens[0].equalsIgnoreCase("genedistancemax"))
						Region.gene_distance_max = Integer.parseInt(tokens[1]);
					else if(tokens[0].equalsIgnoreCase("ByChromosome"))
						by_chromosome = true;
					else if(tokens[0].equalsIgnoreCase("KeepAll"))
						keep_all = true;
					else if(tokens[0].equalsIgnoreCase("NotKeepAll"))
						keep_all = false;
					else if(tokens[0].equalsIgnoreCase("FixForSingleAberrantSample"))
						SignificantRegion.process_one_aberrant_sample = true;
					else if(tokens[0].equalsIgnoreCase("SkipAllNormal"))
						skip_all_normal = true;
					else if(tokens[0].equalsIgnoreCase("MixedAberrations"))
						mixed_aberrations = true;
					else if(tokens[0].equalsIgnoreCase("debug"))
						Region.debug = true;
					else if(tokens[0].equalsIgnoreCase("LimitedPeeloff"))
					{	
						if(tokens.length==4)
						{
							SignificantRegion.sparam=Integer.parseInt(tokens[3].split("=")[1]);
						}
						try {
							SignificantRegion.limited_peeloff = true;
							SignificantRegion.peeloff_gscore_thres[0] = Integer.parseInt(tokens[1]);
							SignificantRegion.peeloff_gscore_thres[1] = Integer.parseInt(tokens[2]);
						} catch (NumberFormatException e) {
							String distr_file_name = tokens[1];
							SignificantRegion.peeloff_qval_thres = Float.parseFloat(tokens[2]);
							if(tokens.length > 2)
								for(Marker.aberration ab : Marker.aberration.values())
								{
									String distrabbfilename=distr_file_name + '.' + ab.name();
									try {
										BufferedReader distr_file = new BufferedReader(new FileReader(distrabbfilename));
										String[] tokens2=null;
										do tokens2 = distr_file.readLine().split("=");
										while(Double.parseDouble(tokens2[1]) > SignificantRegion.peeloff_qval_thres);
										SignificantRegion.peeloff_gscore_thres[ab.ordinal()] = Integer.parseInt(tokens2[0]);
									} catch (FileNotFoundException e1) {
										System.err.println("Could not open limited peel-off q-value file: "+distrabbfilename);														
										System.err.println("This problem needs to be fixed only if you are working with aberrations of type " + ab.name());
									}
								}
						}
					}else if(tokens[0].equalsIgnoreCase("ArmPeeloff")){
						if(bands_file == null){
							System.err.println("arm peeloff option required bands file");
							System.exit(1);
						}
						SignificantRegion.arm_peeloff = true;
					} else if(tokens[0].equalsIgnoreCase("aberrationThresholds")){
						int AMP_index = 1;
						if(tokens[1].equalsIgnoreCase("focal")){
							if(bands_file == null){
								System.err.println("focus thresholds require bands file");
								System.exit(1);
							}
							AMP_index++;
							SNPdata.focal_thresholds = true;
						}
						if(tokens.length > AMP_index){
							Marker.thresholds[0] = Float.parseFloat(tokens[AMP_index]);
							Marker.thresholds[1] = Float.parseFloat(tokens[AMP_index+1]);
						}
					} else if(tokens[0].equalsIgnoreCase("highLevelThresholds")){
						Marker.highlevelthresholds[0] = Float.parseFloat(tokens[1]);
						Marker.highlevelthresholds[1] = Float.parseFloat(tokens[2]);
					} else if(tokens[0].equalsIgnoreCase("caps")){
						Marker.caps[0] = Float.parseFloat(tokens[1]);
						Marker.caps[1] = Float.parseFloat(tokens[2]);
					} else if(tokens[0].equalsIgnoreCase("DistributionOutput"))
						distribution_output = tokens[1];
				}
				infile.close();
				if(sp.getParent() != null)
					dir = sp.getParent() + '/';
			}

			PrintStream fout = new PrintStream(new FileOutputStream(dir+"DistributionCall.log"));
			for (int i = 0; i < args.length; i++) {
				fout.print(" "+args[i]);
			}
			fout.close();
			
			Marker.aberration[] types = Marker.aberration.values();
			if(cpfile == null){
				types = new Marker.aberration[1];
				types[0] = Marker.aberration.LOH;
			} else {
				if(lohfile == null){
					types = new Marker.aberration[2];
					System.arraycopy(Marker.aberration.values(), 0, types, 0, 2);
				}
				if(mixed_aberrations){
					types = new Marker.aberration[types.length-1];
					types[0] = Marker.aberration.AMP;
					if(types.length > 1)
						types[1] = Marker.aberration.LOH;
				}
			}
			Marker.aberrationtypes=types;
			
			System.out.println("Loading genome locations file");
			Genome gene_genome = new Genome(genomefile, bands_file);
			System.out.println("Loading miRNA locations file");			
			Genome miRNA_genome = miRNAfile == null ? null : new Genome(miRNAfile, bands_file);			
			System.out.println("Loading copy number data");
			SNPdata data = new SNPdata(cpfile, lohfile, gene_genome.bands,differentialfile);
			System.out.println("Calculating null distributions");
			distributions = new Distribution[types.length];
			Arrays.fill(distributions, null);
			for(Marker.aberration ab : types){
				distributions[ab.ordinal()] = new Distribution(data, ab);
				distributions[ab.ordinal()].calculateSignificance(data);
				if(distribution_output != null)
					(new PrintStream(dir + distribution_output + '.' + ab.name())).println(distributions[ab.ordinal()].gscore_to_qvalue.toString().replaceAll("[{}]","").replaceAll(", ","\n"));
				if(SignificantRegion.peeloff_gscore_thres[ab.ordinal()] == 0 &&
						!Double.isNaN(SignificantRegion.peeloff_qval_thres))
					if(Math.abs(SignificantRegion.peeloff_qval_thres - qvalthreshold) < .0001)
						SignificantRegion.peeloff_gscore_thres[ab.ordinal()] = distributions[ab.ordinal()].Gsig;
					else {
						Map.Entry<Integer,Double> entry;
						Iterator< Map.Entry<Integer,Double> > iter = distributions[ab.ordinal()].gscore_to_qvalue.entrySet().iterator();
						do
							SignificantRegion.peeloff_gscore_thres[ab.ordinal()] = (entry = iter.next()).getKey().intValue();
						while(entry.getValue().doubleValue() > SignificantRegion.peeloff_qval_thres);
					}
			}
			if(SNPdata.focal_thresholds && data.copy_number_samples != null){
				PrintStream thresholdsfile = new PrintStream(new FileOutputStream(dir + "thresholds.txt"));
				thresholdsfile.println("SampleName\tAMP\tDEL");
				for(int i=0; i < data.copy_number_samples.length; i++)
					thresholdsfile.println(data.copy_number_samples[i] + '\t'
							+ (Marker.per_sample_thresholds[0][i] + Marker.thresholds[0]) + '\t'
							+ (Marker.per_sample_thresholds[1][i] + Marker.thresholds[1]));
				thresholdsfile.close();
			}

						
			PrintStream markerfile = new PrintStream(new FileOutputStream(dir + "markers.txt"));
			String marker_header = "marker\tlocation";
			if(cpfile != null)
				marker_header += "\tAMP qvalue\tAMP Gscore\tDEL qvalue\tDEL Gscore";
			else
				marker_header += "\t\t";
			if(lohfile != null)
				marker_header += "\tLOH qvalue\tLOH Gscore";
			markerfile.println(marker_header);

			/*
			PrintStream genesfile = new PrintStream(new FileOutputStream(dir + "significantGenes.txt"));
			genesfile.println("gene ID\taberration\tqvalue");
			PrintStream miRNA_output_file = null;
			if(miRNA_genome != null){
				miRNA_output_file = new PrintStream(new FileOutputStream(dir + "significantmiRNA.txt"));
				miRNA_output_file.println("gene ID\taberration\tqvalue");
			}
			*/

			PrintStream regionsfile = new PrintStream(new FileOutputStream(dir + "regions.txt"));
			regionsfile.println("region\tqvalue");

			PrintStream peaksfile = new PrintStream(new FileOutputStream(dir + "peaks.txt"));
			peaksfile.println("peak\tqvalue\tbroad\tfocal");
			PrintStream gene_file = null;
			PrintStream miRNA_file = null;
			PrintStream biolearn_gene_file_genatomy = null;
			PrintStream biolearn_miRNA_file_genatomy = null;
			PrintStream biolearn_gene_file = null;
			PrintStream biolearn_miRNA_file = null;
			PrintStream biolearn_file = null;
			PrintStream biolearn_continuous_uncapped_file = null;	
			PrintStream biolearn_continuous_uncapped_gene_file = null;
			//PrintStream biolearn_continuous_uncapped_mean_gene_file = null;
			PrintStream biolearn_continuous_uncapped_miRNA_file = null;
			//PrintStream continuous_uncapped_mean_gene_file = null;
			PrintStream continuous_uncapped_gene_file = null;
			PrintStream continuous_uncapped_miRNA_file = null;
			List< List<Integer> > tumor_mapping = null;
			List<Integer> copy_mapping = new Vector<Integer>();
			List<Integer> loh_mapping = new Vector<Integer>();

			String samples_heading = "Gene";
			for(int i=0; i < data.copy_number_samples.length; i++)
			    samples_heading += "\t" + data.copy_number_samples[i];

				gene_file = new PrintStream(new FileOutputStream(dir + "gene.matrix"));
				gene_file.println(samples_heading);
				continuous_uncapped_gene_file = new PrintStream(new FileOutputStream(dir + "gene.continuous.matrix"));
				continuous_uncapped_gene_file.println(samples_heading);
				//continuous_uncapped_mean_gene_file = new PrintStream(new FileOutputStream(dir + "gene.continuous.mean.matrix"));
				//continuous_uncapped_mean_gene_file.println(samples_heading);
				if(miRNA_genome != null){
					miRNA_file = new PrintStream(new FileOutputStream(dir + "miRNA.matrix"));
					miRNA_file.println(samples_heading);
					continuous_uncapped_miRNA_file = new PrintStream(new FileOutputStream(dir + "miRNA.continuous.matrix"));
					continuous_uncapped_miRNA_file.println(samples_heading);
				}

			
				
			if(tumorfile != null){
				String biolearn_heading = "Region";
				BufferedReader infile = new BufferedReader(new FileReader(tumorfile));
				infile.readLine();
				for(String next = infile.readLine(); next != null; next = infile.readLine())
				{
					String[] tokens = next.split("\\t");
					if(tokens.length == 1 || tokens[1].isEmpty())
						continue;
					biolearn_heading += "\t" + tokens[tokens.length-1];
					if(data.copy_number_samples != null)
						copy_mapping.add(Arrays.asList(data.copy_number_samples).indexOf(tokens[0]));
					if(data.loh_samples != null)
						loh_mapping.add(Arrays.asList(data.loh_samples).indexOf(tokens[0]));
				}
				infile.close();
				tumor_mapping = new Vector< List<Integer> >();
				tumor_mapping.add(copy_mapping);
				tumor_mapping.add(copy_mapping);
				tumor_mapping.add(loh_mapping);
				biolearn_file = new PrintStream(new FileOutputStream(dir + "biolearn.matrix"));
				biolearn_file.println(biolearn_heading);
				biolearn_continuous_uncapped_file = new PrintStream(new FileOutputStream(dir + "biolearn.continuous.matrix"));
				biolearn_continuous_uncapped_file.println(biolearn_heading);
				biolearn_gene_file = new PrintStream(new FileOutputStream(dir + "biolearn.gene.matrix"));
				biolearn_gene_file.println(biolearn_heading);
				
				biolearn_gene_file_genatomy = new PrintStream(new FileOutputStream(dir + "genatomy.gene.matrix"));
				biolearn_gene_file_genatomy.println(biolearn_heading);

				biolearn_continuous_uncapped_gene_file = new PrintStream(new FileOutputStream(dir + "biolearn.gene.continuous.matrix"));
				biolearn_continuous_uncapped_gene_file.println(biolearn_heading);
				//biolearn_continuous_uncapped_mean_gene_file = new PrintStream(new FileOutputStream(dir + "biolearn.gene.continuous.mean.matrix"));
				//biolearn_continuous_uncapped_mean_gene_file.println(biolearn_heading);
				if(miRNA_genome != null){
					biolearn_miRNA_file = new PrintStream(new FileOutputStream(dir + "biolearn.miRNA.matrix"));
					biolearn_miRNA_file.println(biolearn_heading);
					biolearn_miRNA_file_genatomy = new PrintStream(new FileOutputStream(dir + "genatomy.miRNA.matrix"));
					biolearn_miRNA_file_genatomy.println(biolearn_heading);					
					biolearn_continuous_uncapped_miRNA_file = new PrintStream(new FileOutputStream(dir + "biolearn.miRNA.continuous.matrix"));
					biolearn_continuous_uncapped_miRNA_file.println(biolearn_heading);
				}
			}

			PrintStream[] tumor_to_region = new PrintStream[Marker.aberrationtypes.length];
			PrintStream[] continuous_tumors_uncapped = new PrintStream[Marker.aberrationtypes.length];
			PrintStream[] geneoutput = new PrintStream[Marker.aberrationtypes.length];
			PrintStream[] mirnaoutput = new PrintStream[Marker.aberrationtypes.length];			
			
			for(Marker.aberration ab : types){
				tumor_to_region[ab.ordinal()] = new PrintStream(new FileOutputStream(dir + name(ab) + ".tumors.matrix"));
				continuous_tumors_uncapped[ab.ordinal()] = new PrintStream(new FileOutputStream(dir + name(ab) + ".tumors.continuous.matrix"));
				geneoutput[ab.ordinal()]  = new PrintStream(new FileOutputStream(dir + "gene."+name(ab)));
				geneoutput[ab.ordinal()].println("GENEID\tGSCORE\tQVAL");
				String[] tumors = ab == Marker.aberration.LOH
				? data.loh_samples : data.copy_number_samples;
				String tumor_line = "Region";
				for(int tumor=0; tumor < tumors.length; tumor++)
					tumor_line += "\t" + tumors[tumor];
				tumor_to_region[ab.ordinal()].println(tumor_line);
				continuous_tumors_uncapped[ab.ordinal()].println(tumor_line);
				if(miRNA_genome != null){
					mirnaoutput[ab.ordinal()]  = new PrintStream(new FileOutputStream(dir + "mirna."+name(ab)));
					mirnaoutput[ab.ordinal()].println("MIRNAID\tGSCORE\tQVAL");
				}
			}

			
			System.out.println("Obtaining significant regions");
			Map<Gene, List<Region> > gene_regions = new HashMap<Gene, List<Region> >();
			Map<Gene, List<Region> > miRNA_regions = new HashMap<Gene, List<Region> >();
			for(int chromosome = 1; chromosome <= (by_chromosome ? Marker.numchrom : 1); chromosome++){
				if(Distribution.chrom2run!=-1 && Distribution.chrom2run!=chromosome)
				{
					continue;
				}
				data.InMemory(chromosome);
				for(Iterator<Marker> marker_iter = data.iterator(); marker_iter.hasNext(); )
					markerfile.println(marker_iter.next().toString());
				for (int ia = 0; ia < Marker.aberrationtypes.length; ia++) {
					Marker.aberration ab = Marker.aberrationtypes[ia];	
					List<SignificantRegion> regions=new ArrayList<SignificantRegion>();
					if(differentialfile==null)
						regions = data.significantRegions(distributions[ab.ordinal()], gene_genome, miRNA_genome);
					List<Region> important_regions = new Vector<Region>();
					for(Iterator<SignificantRegion> iter = regions.iterator(); iter.hasNext(); )
					{
						SignificantRegion region = iter.next();
						if(!region.equals(region.peaks.get(0)))
							important_regions.add(region);
						for(Iterator<Region> subiter = region.subregions.iterator(); subiter.hasNext(); ){
							Region subregion = subiter.next();
							regionsfile.println(subregion.toString() + '\t' + subregion.max_qvalue);
						}
						regionsfile.println("");
						for(Iterator<PeakRegion> peakiter = region.peaks.iterator(); peakiter.hasNext(); ){
							PeakRegion peak = peakiter.next();
							peaksfile.println(peak.toString() + '\t' + peak.peak_qvalue
									+ '\t' + (region.isBroad() ? 'Y' : 'N') + '\t'
									+ (peak.has_focal ? 'Y' : 'N'));
							important_regions.add(peak);
						}
					}

					String region_line = "";
					for(Iterator<Region> iter = important_regions.iterator(); iter.hasNext(); )
						region_line += "\t" + iter.next().toString();

					String matrix_file = dir + name(ab) + ".genes";
					if(by_chromosome)
						matrix_file += ".chr" + chromosome;
					matrix_file += ".matrix";
					PrintStream genes_to_region = new PrintStream(new FileOutputStream(matrix_file));
					genes_to_region.println(region_line);
					PrintStream miRNA_to_region = null;
					if(miRNA_genome != null){
						miRNA_to_region = new PrintStream(new FileOutputStream(matrix_file.replaceFirst("genes","miRNA")));
						miRNA_to_region.println(region_line);
					}

					for(int genome_index = 0; genome_index < 2; genome_index++){
						Genome genome = genome_index == 0 ? gene_genome : miRNA_genome;
						if(genome == null)
							continue;
						Map<Gene, List<Region> > regions_map = genome_index == 0 ? gene_regions : miRNA_regions;
						for(int gene_index = 0; gene_index < genome.size() && (!by_chromosome || genome.get(gene_index).chromosome <= chromosome); gene_index++)
							if(!by_chromosome || genome.get(gene_index).chromosome == chromosome){
								Gene gene = genome.get(gene_index);								
								if(ab.ordinal() == 0){
									Collection<Marker> closest_markers = gene.findAberrations(data);
									if(closest_markers.size()==0)
									{
										continue;
									}																		
									String line = gene.ID;
									String linegenatomy = gene.ID+",CNV";									
									String continuous_line_uncapped = gene.ID;
									//String continuous_line_uncapped_mean = gene.ID;
									String all_samples_line = gene.ID;
									String all_samples_continuous_line_uncapped = gene.ID;
									//String all_samples_continuous_line_uncapped_mean = gene.ID;
									boolean had_abnormalities = false;
									for(Iterator<Integer> tumor_iter = copy_mapping.iterator(); tumor_iter.hasNext(); ){
										int tumor = tumor_iter.next().intValue();
										float avg_copy_number = Float.POSITIVE_INFINITY;
										float avg_copy_number_uncapped = Float.POSITIVE_INFINITY;		
										float mean_avg_copy_number_uncapped = 0;
										
										for(Iterator<Marker> marker_iter = closest_markers.iterator(); marker_iter.hasNext(); )
										{
											Marker m=marker_iter.next();												
											mean_avg_copy_number_uncapped  += m.copy_number_val_uncapped(tumor) / closest_markers.size();											
											avg_copy_number = Math.min(m.copy_number_val(tumor), avg_copy_number);
											avg_copy_number_uncapped = Math.min(m.copy_number_val_uncapped(tumor), avg_copy_number_uncapped);																							
										}
										
										if(avg_copy_number == Float.POSITIVE_INFINITY)
											avg_copy_number=0;
										if(avg_copy_number_uncapped == Float.POSITIVE_INFINITY)
											avg_copy_number_uncapped=0;
										
										//continuous_line_uncapped_mean += "\t" + mean_avg_copy_number_uncapped;
										continuous_line_uncapped += "\t" + avg_copy_number_uncapped;
										line += '\t';
										linegenatomy += '\t';
										if(avg_copy_number <= Marker.per_sample_thresholds[0][tumor] + Marker.thresholds[0] &&
												avg_copy_number >= Marker.per_sample_thresholds[1][tumor] + Marker.thresholds[1])
										{
											line += 'N';
											linegenatomy+='1';
										}else
										{
											line += avg_copy_number > Marker.thresholds[0]? 'A' : 'D';
											linegenatomy += avg_copy_number > Marker.thresholds[0]? '2' : '0';
										}
									}
									for(int tumor = 0; tumor < data.copy_number_samples.length; tumor++){
										float avg_copy_number = Float.POSITIVE_INFINITY;
										float avg_copy_number_uncapped = Float.POSITIVE_INFINITY;	
										float mean_avg_copy_number_uncapped = 0;
										if(ab==Marker.aberration.DEL)
										{
											avg_copy_number = Float.POSITIVE_INFINITY;
											avg_copy_number_uncapped = Float.POSITIVE_INFINITY;											
										}
										
										for(Iterator<Marker> marker_iter = closest_markers.iterator(); marker_iter.hasNext(); )
										{
											Marker m=marker_iter.next();
											mean_avg_copy_number_uncapped += m.copy_number_val_uncapped(tumor) / closest_markers.size();
											avg_copy_number = Math.min(m.copy_number_val(tumor), avg_copy_number);
											avg_copy_number_uncapped = Math.min(m.copy_number_val_uncapped(tumor), avg_copy_number_uncapped);																						
										}
										
										if(avg_copy_number == Float.POSITIVE_INFINITY)
											avg_copy_number=0;
										if(avg_copy_number_uncapped == Float.POSITIVE_INFINITY)
											avg_copy_number_uncapped=0;
										
										//all_samples_continuous_line_uncapped_mean += "\t" + mean_avg_copy_number_uncapped;
										all_samples_continuous_line_uncapped += "\t" + avg_copy_number_uncapped;
										all_samples_line += '\t';
										if(avg_copy_number <= Marker.per_sample_thresholds[0][tumor] + Marker.thresholds[0] &&
												avg_copy_number >= Marker.per_sample_thresholds[1][tumor] + Marker.thresholds[1])
											all_samples_line += 'N';
										else {
											had_abnormalities = true;
											all_samples_line += avg_copy_number > Marker.thresholds[0]
											                                            ? 'A' : 'D';
										}
									}
									if(!skip_all_normal || had_abnormalities){			
									    if(tumorfile != null){
											(genome_index == 0 ? biolearn_gene_file : biolearn_miRNA_file).println(line);
											(genome_index == 0 ? biolearn_gene_file_genatomy : biolearn_miRNA_file_genatomy).println(linegenatomy);
											(genome_index == 0 ? biolearn_continuous_uncapped_gene_file : biolearn_continuous_uncapped_miRNA_file).println(continuous_line_uncapped);
											/*
											if(genome_index == 0)
												biolearn_continuous_uncapped_mean_gene_file.println(continuous_line_uncapped_mean);
												*/
									    }
									    (genome_index == 0 ? gene_file : miRNA_file).println(all_samples_line);
									    (genome_index == 0 ? continuous_uncapped_gene_file : continuous_uncapped_miRNA_file).println(all_samples_continuous_line_uncapped);
									    /*
									    if(genome_index == 0)
									    	continuous_uncapped_mean_gene_file.println(all_samples_continuous_line_uncapped_mean);
									    	*/
									}
								}
								if(gene.Gscore==null)
									continue;
								String linegene=gene.ID+"\t"+gene.Gscore[ab.ordinal()]+"\t"+gene.qvalue[ab.ordinal()];
								(genome_index == 0 ? geneoutput[ab.ordinal()] : mirnaoutput[ab.ordinal()]).println(linegene);								

								String gene_line = gene.ID;
								List<Region> regionsForGene = new Vector<Region>();
								for(Iterator<Region> region_iter = important_regions.iterator(); region_iter.hasNext(); )
								{
									Region region = region_iter.next();
									if(gene.isInRegion(region))
									{
										gene_line += "\t1";
										regionsForGene.add(region);
									} else
										gene_line += "\t0";
								}

								if(!regionsForGene.isEmpty()){
									(genome_index == 0 ? genes_to_region : miRNA_to_region).println(gene_line);
									List<Region> old_genes = regions_map.get(gene);
									if(old_genes == null)
										regions_map.put(gene, regionsForGene);
									else
										old_genes.addAll(regionsForGene);
								}
							}
					}
					genes_to_region.close();
					if(miRNA_to_region != null)
						miRNA_to_region.close();

					for(Iterator<Region> iter = important_regions.iterator(); iter.hasNext(); ){
						Region region = iter.next();
						String continuous_line = region_line = region.toString();
						String continuous_line_uncapped = region_line = region.toString();
						String[] tumors = ab == Marker.aberration.LOH
						? data.loh_samples : data.copy_number_samples;
						for(int tumor=0; tumor < tumors.length; tumor++){
							int[] strength_count = new int[3];
							Arrays.fill(strength_count, 0);
							float total_strength = 0;
							float total_strength_uncapped = 0;
							for(Iterator<Marker> marker_iter = region.markers.iterator(); marker_iter.hasNext(); )
							{
								Marker marker = marker_iter.next();
								strength_count[marker.aberrationStrength(ab, tumor)]++;
								total_strength += (ab == Marker.aberration.LOH
										? marker.isLOH : marker.copy_number)
										[tumor];
								total_strength_uncapped += (ab == Marker.aberration.LOH
										? marker.isLOH : marker.copy_number_uncapped)
										[tumor];
							}
							region_line += "\t" +
							(strength_count[2] > region.markers.size()/2 ? '2' :
								(strength_count[0] >= (region.markers.size()+1)/2 ? '0' : '1'));
							continuous_line += "\t" + (total_strength / region.markers.size());
							continuous_line_uncapped += "\t" + (total_strength_uncapped / region.markers.size());
						}
						tumor_to_region[ab.ordinal()].println(region_line);
						continuous_tumors_uncapped[ab.ordinal()].println(continuous_line_uncapped);
						if(biolearn_file != null){
							String biolearn_line = region.toString();
							String biolearn_continuous_line = region.toString();
							String biolearn_continuous_uncapped_line = region.toString();
							for(Iterator<Integer> tumor_iter = tumor_mapping.get(ab.ordinal()).iterator(); tumor_iter.hasNext(); ){
								int tumor = tumor_iter.next().intValue();
								int[] strength_count = new int[3];
								Arrays.fill(strength_count, 0);
								float total_strength = 0;
								float total_strength_uncapped = 0;		    
								if(tumor < 0)
									strength_count[0] = region.markers.size();
								else
									for(Iterator<Marker> marker_iter = region.markers.iterator(); marker_iter.hasNext(); )
									{
										Marker marker = marker_iter.next();
										strength_count[marker.aberrationStrength(ab, tumor)]++;
										total_strength += (ab == Marker.aberration.LOH
												? marker.isLOH : marker.copy_number)
												[tumor];
										total_strength_uncapped += (ab == Marker.aberration.LOH
												? marker.isLOH : marker.copy_number_uncapped)
												[tumor];
									}
								biolearn_line += "\t" +
								(strength_count[2] > region.markers.size()/2 ? '2' :
									(strength_count[0] >= region.markers.size()/2 ? '0' : '1'));
								biolearn_continuous_line += "\t" + (total_strength / region.markers.size());
								biolearn_continuous_uncapped_line += "\t" + (total_strength_uncapped / region.markers.size());
							}
							biolearn_file.println(biolearn_line);
							biolearn_continuous_uncapped_file.println(biolearn_continuous_uncapped_line);
						}
						region.markers = null;
					}
				}
			}
			
			// closing all files						
			markerfile.close();
			/*
			genesfile.close();
			if(miRNA_genome != null){
				miRNA_output_file.close();
			}*/
			regionsfile.close();
			peaksfile.close();
			gene_file.close();
			continuous_uncapped_gene_file.close();
			//continuous_uncapped_mean_gene_file.close();
			if(miRNA_genome != null){
				miRNA_file.close();
				continuous_uncapped_miRNA_file.close();
			}
			
			if(tumorfile != null){
				biolearn_file.close();
				biolearn_continuous_uncapped_file.close();
				biolearn_gene_file.close();
				biolearn_gene_file_genatomy.close();
				biolearn_continuous_uncapped_gene_file.close();
				//biolearn_continuous_uncapped_mean_gene_file.close();
				if(miRNA_genome != null){
					biolearn_miRNA_file.close();
					biolearn_miRNA_file_genatomy.close();
					biolearn_continuous_uncapped_miRNA_file.close();
				}
			}
			

			for (int ia = 0; ia < Marker.aberrationtypes.length; ia++) {
				Marker.aberration ab = Marker.aberrationtypes[ia];
				tumor_to_region[ab.ordinal()].close();
				geneoutput[ab.ordinal()].close();
				if(miRNA_genome != null){
					mirnaoutput[ab.ordinal()].close();
				}
				continuous_tumors_uncapped[ab.ordinal()].close();				
			}			
			
			System.out.println("Significant regions obtained");
			printGeneRegion(gene_regions, "geneToRegion.txt");
			printGeneRegion(miRNA_regions, "miRNAToRegion.txt");																
		} catch (Exception e) {
			e.printStackTrace(System.err);
		}
	}

	static void printGeneRegion(Map<Gene, List<Region> > map, String file)
	throws Exception
	{
		if(map.isEmpty())
			return;
		PrintStream geneToRegion = new PrintStream(new FileOutputStream(dir + file));
		for(Iterator<Map.Entry<Gene, List<Region> > > iter = map.entrySet().iterator(); iter.hasNext(); )
		{
			Map.Entry<Gene, List<Region> > entry = iter.next();
			Collections.sort(entry.getValue(), entry.getKey());
			if(entry.getValue().size() > 1 &&
					entry.getKey().compare(entry.getValue().get(0), entry.getValue().get(1)) == 0)
				throw new Exception(entry.getKey().ID + " is in two different regions, " + entry.getValue().get(0).toString() + " and " + entry.getValue().get(1).toString());
			geneToRegion.println(entry.getKey().ID + '\t' + 
					entry.getValue().get(0).toString());
		}
		geneToRegion.close();
	}
}
