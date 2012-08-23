import java.io.File;
import java.util.ArrayList;

import shared.FileOps;
import shared.Utils;

/**
 * LOHcate --- A software tool for LOH calling and visualization in cancer genomes
 * D Wheeler & SG Reddy
 * Human Genome Sequencing Center, Baylor College of Medicine (Houston, TX)
 * 
 * This runner class is the big enchilada. See main method for 'switchboard' / pipeline overview.
 * 
 * @author Siddharth G. Reddy
 *
 */
public class Script {
	
	public static String[] cluster_names = {"dup", "loh", "roc-loh", "het"}; //het needs to be the last element of cluster_names, but more elem.s can be added to the 'front' (as long as you handle them in the getClusters() method )
	
	private static final float NAF_STRIP_EXPANDER = 1.25f; //# of std. deviations to move away from the mean (when defining the thickness of the horizontal band containing HET ball, LOH sidelobes, &c.)
	private static final float HET_BALL_EPS = 0.035f, DUP_WEDGE_LASSO = 0.015f; //DBSCAN parameters for HET ball / DUP wedge detection
	private static final int HET_BALL_MINPTS = 100, DUP_WEDGE_MINPTS = 100; //DBSCAN parameters for HET ball / DUP wedge detection
	
	private static final int REGION_SEGMENTATION_DIST_THRESHOLD = 2000000; //greatest possible distance between 2 'adjacent' points of LOH in a region of 'contiguous' LOH
	
	//Vv.vV is useful for DBSCAN parameter tuning. once you've selected a couple of mouse ear plots that you'd like to be part of your test set, 'whitelist' them by uncommenting the extension of the patient-loop conditional in the curateSNPCalls() method (and, of course, uncommenting the wList() method)
	/*public static boolean wList(String fname) {
		//target-aml --> String[] wlist = {"PAKKBK-04", "PAMYAS-04", "PANLIZ-03", "PANLRE-04"};////target-aml --> 
		//target-all --> String[] wlist = {"PAPJIB-03", "PAPLDL-03", "PARAKF-03", "PASLZM-03"};
		//pvera --> String[] wlist = {"-WS-G", "-LPC-G", "-TC-G"};
		//hepato --> String[] wlist = {"021-T", "044-T", "024-T", "023-T"};
		//renal-we --> String[] wlist = {"CW-5590", "CZ-5463", "CZ-5453", "CJ-5681"};
		for (String elem : wlist)
			if (fname.indexOf(elem)!=-1)
				return true;
		return false;
	}*/
	
	/**
	 * Curate SNP calls by clustering data points into HET ball, DUP wedge, LOH sidelobes, &c., and by grabbing dbsnp population allele frequencies when possible
	 * @param inDir naf-taf-inputs
	 * @param opt 0::Illumina, 1::SOLiD
	 */
	public static void curateSNPCalls(String inDir, String outDir, int opt) {
		File[] files = (new File(inDir)).listFiles();
		String[] load, split;
		String som_load, germ_load, germ_som, toWrite, mut_type, gene, temp;
		int cluster, len_gload, som_start;
		int[] cluster_get;
		for (File file : files) {
			if (file.getName().indexOf("germline")!=-1) {// && wList(file.getName())) {
				System.out.println(file.getName().replace(".germline.txt", ""));
				som_load = FileOps.loadFromFile(file.getAbsolutePath().replace("germline", "somatic"));
				germ_load = FileOps.loadFromFile(file.getAbsolutePath()); len_gload = germ_load.split("\n").length;
				try {
					load = (germ_load + som_load.replace(som_load.split("\n")[0] + "\n", "")).split("\n");
					som_start = germ_load.split("\n").length;
				} catch (Exception e) { 
					try {
						load = germ_load.split("\n");
						som_start = load.length + 1;
					} catch (Exception e2) {
						load = som_load.split("\n");
						som_start = 0;
					}
				}
				som_load = ""; germ_load = "";
				FileOps.writeToFile(outDir + "/" + file.getName().replace(".germline.txt", ".csv"), "chr,pos,n_vaf,t_vaf,allele_freq,gene,mutation_type,germ_som,cluster\n");
				cluster_get = getClusters(load, som_start, opt); //get cluster assignments (HET ball, LOH sidelobes, DUP wedge, &c.)
				for (int i = 1; i<load.length; i++) {
					if (load[i].indexOf("refName")==-1 && load[i].indexOf("coord")==-1) { //upstream pipelines will randomly spit header lines into the middle of a naf-taf-input file. we're just avoiding those
						try {
							cluster = cluster_get[i-1];//assignCluster(Float.parseFloat(load[i].split("\t")[7]), Float.parseFloat(load[i].split("\t")[8]));
							
							if (opt==1) //SOLiD
								toWrite = "chr" + load[i].split("\t")[0].replace("chr", "") + "," + load[i].split("\t")[1] + "," + (Float.parseFloat(load[i].split("\t")[6]) / Float.parseFloat(load[i].split("\t")[5])) + "," + (Float.parseFloat(load[i].split("\t")[4]) / Float.parseFloat(load[i].split("\t")[3])) + ","; //chr,pos,n_vaf,t_vaf
							else //Illumina
								toWrite = load[i].split("\t")[0] + "," + load[i].split("\t")[1] + "," + load[i].split("\t")[7] + "," + load[i].split("\t")[8] + ","; //chr,pos,n_vaf,t_vaf
							
							if (opt==1) { //SOLiD
								if (load[i].split("\t")[7].indexOf("novel")!=-1)
									temp = "novel";
								else //strangely, the variant base is n/a in the SOLiD naf-taf-inputs (and there's not much point in looking up the reference base's allele frequency)
									temp = "N/A";
							}
							else { //Illumina
								if (load[i].split("\t")[9].indexOf("novel")!=-1)
									temp = "novel";
								else {
									try {
										temp = load[i].split("\t")[9].split(load[i].split("\t")[2] + "\\|")[1];
										try {
											temp = temp.split(",")[0];
										} catch (Exception e) { temp = temp.split(";")[0]; }
										try {
											temp = temp.split(";")[0];
										} catch (Exception e) { temp = temp.split(",")[0]; }
									} catch (Exception e) { temp = "N/A"; }
								}
							}
							toWrite += temp + ","; //allele_freqs
							germ_som = "germline";
							if (i > len_gload)
								germ_som = "somatic";
							
							if (opt==1) { //SOLiD
								try {
									mut_type = load[i].split("\t")[9].replace("syn", "synonymous").replace("nsynonymous", "nonsynonymous");
								} catch (Exception e) { mut_type = "N/A"; }
								try {
									gene = load[i].split("\t")[8];
								} catch (Exception e) { gene = "N/A"; }
								FileOps.appendToFile(outDir + "/" + file.getName().replace(".germline.txt", ".csv"), toWrite + gene + "," + mut_type + "," + germ_som + "," + Integer.toString(cluster) + "\n"); //gene, mutation_type, cluster
								if (cluster==2) //include roc-loh points in general loh pool
									FileOps.appendToFile(outDir + "/" + file.getName().replace(".germline.txt", ".csv"), toWrite + gene + "," + mut_type + "," + germ_som + ",1\n"); //gene, mutation_type, cluster
							}
							else { //Illumina
								FileOps.appendToFile(outDir + "/" + file.getName().replace(".germline.txt", ".csv"), toWrite + load[i].split("\t")[11] + "," + load[i].split("\t")[10].split("_SNV")[0] + "," + germ_som + "," + Integer.toString(cluster) + "\n"); //gene, mutation_type, cluster
								if (cluster==2) //include roc-loh points in general loh pool
									FileOps.appendToFile(outDir + "/" + file.getName().replace(".germline.txt", ".csv"), toWrite + load[i].split("\t")[11] + "," + load[i].split("\t")[10].split("_SNV")[0] + "," + germ_som + ",1\n"); //gene, mutation_type, cluster
							}
						} catch (Exception e) { }//e.printStackTrace(); }
					}
				}
			}
		}
	}
	
	//Vv.vV implements our old window-based clustering method. ugly, isn't it?
	/*public static int assignCluster(float N, float T) {
		if ((.35 < N && N < .65) && ((.2 < T && T < .375))) //left-of-center duplication
			return 0;
		else if ((.35 < N && N < .65) && ((.625 < T && T < .75))) //right-of-center duplication
			return 0;
		else if ((.35 < N && N < .65) && T < .2) //left-of-center loss of heterozygosity
			return 1;
		else if ((.35 < N && N < .65) && T > .75) //right-of-center loss of heterozygosity
			return 2;
		else if ((.35 < N && N < .65) && (.2 < T && T < .75)) //central, spherical cluster (heterozygotes)
			return 3;
		else
			return cluster_names.length;
	}*/
	
	/**
	 * A helper method for curateSNPCalls()
	 * @param load line-split FileOps.loadFromFile of naf-taf-input
	 * @param som_start start index of .somatic.txt naf-taf-input data in 'load' parameter
	 */
	public static int[] getClusters(String[] load, int som_start, int opt) {
		int[] rtn = new int[load.length - 1];
		
		//NAF {het, loh, dup} FRAME definition via peak detection and parameter-tuned standard deviation expansion
		float NAF_frame_lbound = 0.2f, NAF_frame_ubound = 0.8f; //we have to avoid the often hugely dense peak of homozygous mutations (AF > 0.8) and the occasionally hugely dense peak of neg. tail noise / somatics / &c. (AF < 0.2)
		float bin_size = 0.025f; //smoothing parameter
		float std_dev = 0f; //standard deviation in NAF-coord across {0.2 < NAF <= 0.8} variants in VAF plot
		float count_mean = 0f;
		int[] bin_count = new int[(int) ((NAF_frame_ubound - NAF_frame_lbound) / bin_size)];
		int count = 0;
		for (float i = NAF_frame_lbound + bin_size; i<= NAF_frame_ubound; i+= bin_size) { //iterate through bins (<=> horizontal strips of thickness 'bin_size' on VAF plot)
			bin_count[count] = 0;
			for (int row = 1; row<load.length; row++) { //iterate through points
				try {
					if (opt==1) {//SOLiD
						if (i - bin_size < (Float.parseFloat(load[row].split("\t")[6]) / Float.parseFloat(load[row].split("\t")[5])) && (Float.parseFloat(load[row].split("\t")[6]) / Float.parseFloat(load[row].split("\t")[5])) <= i)
							bin_count[count]++;
					}
					else { //Illumina
						if (i - bin_size < Float.parseFloat(load[row].split("\t")[7]) && Float.parseFloat(load[row].split("\t")[7]) <= i) //counting points that fall in current bin
							bin_count[count]++;
					}
				} catch (Exception e) { }
			}
			count_mean += i * bin_count[count];
			count++;
		}
		int tot_count = 0;
		for (int elem : bin_count)
			tot_count += elem;
		count_mean /= (float)tot_count;
		for (int i = 0; i<bin_count.length; i++) //calculate std. deviation of # points in each bin
			std_dev += Math.pow(count_mean - ((i+1) * bin_size + NAF_frame_lbound) * bin_count[i], 2);
		std_dev = (float)Math.sqrt(std_dev) / (float)tot_count;
		NAF_frame_lbound = count_mean - (NAF_STRIP_EXPANDER * std_dev);
		NAF_frame_ubound = count_mean + (NAF_STRIP_EXPANDER * std_dev);
		//System.out.println("MEAN = " + count_mean + " | STD_DEV = " + std_dev);
		//System.out.println("NAF FRAME = " + (count_mean - std_dev) + " - " + (count_mean + std_dev));
		
		//apply DBScan to points within NAF frame
		ArrayList<Floint> points = new ArrayList<Floint>();
		for (int i = 1; i<load.length; i++) {
			try {
				if (opt==1) { //SOLiD
					if ((Float.parseFloat(load[i].split("\t")[6]) / Float.parseFloat(load[i].split("\t")[5])) > NAF_frame_lbound && (Float.parseFloat(load[i].split("\t")[6]) / Float.parseFloat(load[i].split("\t")[5])) <= NAF_frame_ubound)
						points.add(new Floint((Float.parseFloat(load[i].split("\t")[4]) / Float.parseFloat(load[i].split("\t")[3])), (Float.parseFloat(load[i].split("\t")[6]) / Float.parseFloat(load[i].split("\t")[5])), i-1));
				}
				else { //Illumina
					if (Float.parseFloat(load[i].split("\t")[7]) > NAF_frame_lbound && Float.parseFloat(load[i].split("\t")[7]) <= NAF_frame_ubound)
						points.add(new Floint(Float.parseFloat(load[i].split("\t")[8]), Float.parseFloat(load[i].split("\t")[7]), i-1));
				}
			} catch (Exception e) { }
		}
		Floint[] param = new Floint[points.size()];
		for (int i = 0; i<param.length; i++)
			param[i] = points.get(i);
		
		//Vv.vV well/poorly tuned parameters for different data sets
		//data set --> (NAF_STRIP_EXPANDER, (HET_BALL_EPS, HET_BALL_MINPTS), (DUP_WEDGE_LASSO, DUP_WEDGE_MINPTS))
		//target-aml --> (1, (0.035, 100), (0.015, 100))
		//target-all --> (1, (0.035, 100), (0.01, 100))
		//pvera --> (1.25, (0.035, 100), (0.015, 100))
		//hepato --> (1.25, (0.035, 100), (0.015, 100))
		//renal-we --> (1, (0.05, 500), (0.02, 350))
		
		DBSCAN obj = new DBSCAN(param, HET_BALL_EPS, HET_BALL_MINPTS); //parameters for capturing HET ball//KDBSCAN(param, 10, 0.0325f);
		obj.cluster();
		int het_ball = obj.getLargestCluster();
		int[] dbscan = obj.getClustAssignments();
		obj = new DBSCAN(param, HET_BALL_EPS + DUP_WEDGE_LASSO, DUP_WEDGE_MINPTS); //loosened parameter for capturing DUP lasso (around HET ball)//KDBSCAN(param, 10, 0.027f);
		obj.cluster();
		int[] dup_inc = obj.getClustAssignments();
		for (int i = 0; i<dbscan.length; i++)
			if (dbscan[i]!=dup_inc[i])
				dbscan[i] = -1;

		count = 0;
		for (int i = 1; i<load.length; i++) {
			try {
				if (i > som_start) //if we are now iterating through points coming from the *.somatic.txt input file
					rtn[i-1] = -1; //somatic
				else if ((opt==1 && (NAF_frame_lbound < (Float.parseFloat(load[i].split("\t")[6]) / Float.parseFloat(load[i].split("\t")[5])) && (Float.parseFloat(load[i].split("\t")[6]) / Float.parseFloat(load[i].split("\t")[5])) <= NAF_frame_ubound)) //SOLiD
						|| (opt!=1 && (NAF_frame_lbound < Float.parseFloat(load[i].split("\t")[7]) && Float.parseFloat(load[i].split("\t")[7]) <= NAF_frame_ubound))) { //Illumina
					if (dbscan[count]==het_ball)
						rtn[i-1] = cluster_names.length - 1; //HET
					else if (dbscan[count]==-1)
						rtn[i-1] = 0; //DUP
					else //anything not in the HET ball / DUP wedge is considered part of a LOH sidelobe
						rtn[i-1] = 1; //LOH
					if (opt==1) { //SOLiD
						if (rtn[i-1]==1 && (Float.parseFloat(load[i].split("\t")[4]) / Float.parseFloat(load[i].split("\t")[3])) > 0.5)
							rtn[i-1] = 2; //right-of-center LOH
					}
					else { //Illumina
						if (rtn[i-1]==1 && Float.parseFloat(load[i].split("\t")[8]) > 0.5)
							rtn[i-1] = 2; //right-of-center LOH
					}
					count++;
				}
				else
					rtn[i-1] = cluster_names.length; //outside NAF frame (<=> 'other')
			} catch (Exception e) { }
		}
		return rtn;
	}
	
	/**
	 * Define 'contiguous' regions of LOH, given our curated SNP calls.
	 * @param inDir curated SNP calls
	 */
	public static void segmentRegions(String inDir, String outDir) {
		File[] files = (new File(inDir)).listFiles();
		String[] split;
		String[][] toWrite = new String[22][cluster_names.length - 1]; //[chromosome][cluster]
		String load;
		for (int i = 0 ; i<toWrite.length; i++)
			for (int j = 0; j<toWrite[i].length; j++)
				toWrite[i][j] = "";
		String chr;
		int start, prev = -1, min, min_ind, prev_min;
		//boolean het_break;
		ArrayList<Integer> sorted;
		for (File file : files) { //iterate through patients
			if (file.getName().indexOf(".csv")!=-1) {
				System.out.println(file.getName());
				for (int i = 0; i<toWrite.length; i++)
					for (int j = 0; j<toWrite[i].length; j++)
						toWrite[i][j] += file.getName().replace(".csv", "") + ",";
				load = FileOps.loadFromFile(file.getAbsolutePath());
				for (int i = 1; i<=22; i++) { //iterate through chromosomes
					try {
						if (i<=21)
							chr = Integer.toString(i);
						else
							chr = "X";
						System.out.print(chr + " ");
						split = load.split("\nchr" + chr + ",");
						for (int k = 0; k<cluster_names.length - 1; k++) { //iterate through {dup, loh, roc-loh}
							start = -1;
							prev = -1;
							//het_break = false;
							sorted = new ArrayList<Integer>();
							prev_min = Integer.MIN_VALUE;
							while (sorted.size() < split.length - 1) { //sort SNPs by 'geographic' order on chromosome
								//System.out.println(sorted.size() + " of " + split.length);
								min_ind = -1;
								min = Integer.MAX_VALUE;
								for (int j = 1; j<split.length; j++) {
									if (prev_min < Integer.parseInt(split[j].split(",")[0]) && Integer.parseInt(split[j].split(",")[0]) < min) {
										min = Integer.parseInt(split[j].split(",")[0]);
										min_ind = j;
									}
								}
								sorted.add(min_ind);
								prev_min = min;
								
								//instead of waiting for the points to get sorted then iterating through the sorted list, we're just going to handle points as they're grabbed
																								
								//if (Integer.parseInt(split[min_ind].split("\n")[0].split(",")[7])==cluster_names.length - 1)
									//het_break = true;
								if (Integer.parseInt(split[min_ind].split("\n")[0].split(",")[7])==k) { //if point is part of cluster[k]
									if (start==-1) //init
										start = Integer.parseInt(split[min_ind].split(",")[0]);
									else {
										if ((Integer.parseInt(split[min_ind].split(",")[0]) - prev > REGION_SEGMENTATION_DIST_THRESHOLD) && Integer.parseInt(split[min_ind].split(",")[0]) > start) {//2nd condition prevents nonsensical data (i.e. two points, same location) from throwing exceptions) // || het_break) && Integer.parseInt(split[min_ind].split(",")[0]) > start) {
											if (prev > start) //we need this condition for patients that may only have 1 mutation (unlikely, but should be caught before throwing an exception)
												toWrite[i-1][k] += start + "-" + prev + ","; //add region to row
											start = Integer.parseInt(split[min_ind].split(",")[0]); //begin defining a new region
											//het_break = false;
										}
										prev = Integer.parseInt(split[min_ind].split(",")[0]); //keep track of previous endpoint of region
									}
								}
							}	
							if (start!=-1 && prev!=-1) //if there's only 1 region, then the while loop will only run once => there won't be any point beyond the 1 region we're looking at => toWrite will never be accessed
								toWrite[i-1][k] += start + "-" + prev + ",";
						}
					} catch (Exception e) { e.printStackTrace(); }
				}
				System.out.println();
				for (int i = 0; i<toWrite.length; i++) { //iterate through chromosomes
					if (i<21)
						chr = Integer.toString(i+1);
					else
						chr = "X";
					for (int j = 0; j<toWrite[i].length; j++) { //iterate through clusters
						if (toWrite[i][j].split(",").length > 1)
							FileOps.appendToFile(outDir + "/chr" + chr + "/" + cluster_names[j] + ".csv", toWrite[i][j].substring(0, toWrite[i][j].length() - 1) + "\n");
						toWrite[i][j] = "";
					}
				}
			}
		}
	}
	
	/**
	 * Now that we have our regions of LOH, let's 'pick up' any somatic mutations that fall within these regions (we'll call them LOH too). somatics don't separate into nice sidelobes in allele-fraction plots.
	 * @param snps_inDir curated SNP calls
	 * @param reg_inDir region segmentation data
	 */
	public static void pickupSomatics(String snps_inDir, String reg_inDir, String outDir) {
		File[] files = (new File(snps_inDir)).listFiles();
		String[] load, reg_load;
		boolean toBreak, toAppend;
		for (File file : files) { //iterate through patients
			if (file.getName().indexOf("csv")!=-1) {
				System.out.println(file.getName().replace(".csv", ""));
				load = FileOps.loadFromFile(file.getAbsolutePath()).split("\n");
				FileOps.writeToFile(outDir + "/" + file.getName(), load[0] + "\n");
				for (int i = 1; i<load.length; i++) { //iterate through variants
					toAppend = true;
					if (load[i].split(",")[8].equals("-1")) { //if somatic
						System.out.println(i + " of " + load.length);
						//find out if variant lies within region of continuous DUP/LOH...
						for (int k = 0; k<cluster_names.length - 1; k++) { //iterate through {dup, loh, roc-loh}
							reg_load = FileOps.loadFromFile(reg_inDir + "/chr" + load[i].split(",")[0].replace("chr", "") + "/" + cluster_names[k] + ".csv").split("\n");
							toBreak = false;
							for (int row = 0; row<reg_load.length; row++) { //iterate through patients
								for (int col = 1; col<reg_load[row].split(",").length; col++) { //iterate through regions
									if (Integer.parseInt(reg_load[row].split(",")[col].split("-")[0]) < Integer.parseInt(load[i].split(",")[1])
											&& Integer.parseInt(load[i].split(",")[1]) <= Integer.parseInt(reg_load[row].split(",")[col].split("-")[1])) { //if point falls within region
										FileOps.appendToFile(outDir + "/" + file.getName(), load[i].replace(",-1\n", ",1\n"));
										toAppend = false;
										toBreak = true; //since the point has already been positively identified, we don't need to keep iterating through patients, regions
										break; //ditto
									}
								}
								if (toBreak)
									break;
							}
						}
					}
					if (toAppend)
						FileOps.appendToFile(outDir + "/" + file.getName(), load[i] + "\n");
				}
			}
		}
	}
	
	/**
	 * Score regions of 'contiguous' LOH based on recurrence across patients, variant density, &c. Key feature of this method is that it separates regions into their 'least common denominators', which allows for detection of regions that, while smaller than the ones defined by segmentRegions(), may be more recurrent across patients.
	 * @param regions_inDir region segmentation data
	 * @param snps_inDir curated SNP calls
	 */
	public static void scoreRegions(String regions_inDir, String snps_inDir, String outDir) {
		File[] dirs = (new File(regions_inDir)).listFiles();
		String[] load, pat_load;
		String toWrite;
		Point xth, yth;
		ArrayList<Point> regions;
		ArrayList<ArrayList<String>> region_pats; //stores patients in which each region is found
		ArrayList<Double> recur; //stores region recurrence scores
		float event_density, het_density;
		int ind, temp, big_ind, big_count, big_het_count;
		for (File dir : dirs) { //iterate through chromosomes
			if (dir.getName().indexOf("chr")!=-1) {				
				System.out.println(dir.getName());
				for (int i = 0; i<cluster_names.length - 1; i++) { //iterate through {dup, loh, roc-loh}
					regions = new ArrayList<Point>(); //Point.x --> region start, Point.y --> region end
					recur = new ArrayList<Double>();
					region_pats = new ArrayList<ArrayList<String>>();
					try {
						load = FileOps.loadFromFile(dir.getAbsolutePath() + "/" + cluster_names[i] + ".csv").split("\n");
						//calculate recurrence score
						for (int x = 0; x<load.length; x++) { //iterate through patients
							System.out.println(x + " of " + load.length);
							for (int l = 1; l<load[x].split(",").length; l++) { //iterate through regions
								xth = new Point(Integer.parseInt(load[x].split(",")[l].split("-")[0]), Integer.parseInt(load[x].split(",")[l].split("-")[1]));
								big_ind = indexOf(xth, regions);
								if (big_ind==-1) { //if region hasn't been seen before
									regions.add(xth);
									recur.add(0.0);
									region_pats.add(new ArrayList<String>());
									region_pats.get(region_pats.size() - 1).add(load[x].split(",")[0]);
									if (regions.size()==0)
										big_ind = 0;
									else
										big_ind = regions.size() - 1;
								}
								for (int y = 0; y<load.length; y++) { //iterate through patients, excluding x-th
									if (y!=x) {
										for (int m = 1; m<load[y].split(",").length; m++) { //iterate through y-th patient's regions
											yth = new Point(Integer.parseInt(load[y].split(",")[m].split("-")[0]), Integer.parseInt(load[y].split(",")[m].split("-")[1]));
											if (yth.x < xth.x && yth.y > xth.y) { //if region Y surrounds region X
												recur.set(big_ind, recur.get(big_ind) + 1.0); //increase region X's recurrence score
												region_pats.get(big_ind).add(load[y].split(",")[0]);
												break;
											}
											else if (yth.x < xth.x && (xth.x < yth.y && yth.y < xth.y)) { //if region X slides to the right of region Y (they are still intersecting)
												ind = indexOf((new Point(xth.x, yth.y)), regions);
												if (ind!=-1) { //if we have seen region overlap(X, Y) before
													recur.set(ind, recur.get(ind) + 1.0);
													if (Utils.indexOf(region_pats.get(ind), load[x].split(",")[0])==-1)
														region_pats.get(ind).add(load[x].split(",")[0]);
													if (Utils.indexOf(region_pats.get(ind), load[y].split(",")[0])==-1)
														region_pats.get(ind).add(load[y].split(",")[0]);
												}
												else {
													regions.add(new Point(xth.x, yth.y));
													recur.add(1.0);
													region_pats.add(new ArrayList<String>());
													region_pats.get(region_pats.size() - 1).add(load[x].split(",")[0]);
													region_pats.get(region_pats.size() - 1).add(load[y].split(",")[0]);
												}
												break;
											}
											else if (xth.x == yth.x && xth.y == yth.y) { //if regions X and Y have the same bounds
												recur.set(big_ind, recur.get(big_ind) + 1.0);
												region_pats.get(big_ind).add(load[y].split(",")[0]);
												break;
											}
										}
									}
								}	
							}
						}
						//calculate variant/het density scores...
						FileOps.writeToFile(outDir + "/" + dir.getName() + "/" + cluster_names[i] + ".csv", "region,size,recurrence,num_events,num_hets,event_density,het_density\n");//,exclusivity\n");
						for (int r = 0; r<regions.size(); r++) {
							System.out.println(r + " of " + regions.size());
							big_count = 0;
							big_het_count = 0;
							for (String pat : region_pats.get(r)) { //iterate through patients in which region is found
								pat_load = FileOps.loadFromFile(snps_inDir + "/" + pat + ".csv").split("\n" + dir.getName() + ",");
								for (int j = 1; j<pat_load.length; j++) { //iterate through positions
									if (Integer.parseInt(pat_load[j].split("\n")[0].split(",")[7])==i) { //if point belongs to cluster i
										if (Integer.parseInt(pat_load[j].split(",")[0]) >= regions.get(r).x && Integer.parseInt(pat_load[j].split(",")[0]) <= regions.get(r).y)
											big_count++;
									}
									else if (Integer.parseInt(pat_load[j].split("\n")[0].split(",")[7])==cluster_names.length - 1) { //if point is part of HET ball
										if (Integer.parseInt(pat_load[j].split(",")[0]) >= regions.get(r).x && Integer.parseInt(pat_load[j].split(",")[0]) <= regions.get(r).y)
											big_het_count++;
									}
								}
							}
							event_density = ((big_count / region_pats.get(r).size()) / (float)(regions.get(r).y - regions.get(r).x));
							if (Float.toString(event_density).equals("Infinity")) //div by 0 (usually occurs when # events == 1, which implies event range == 0)
								event_density = -1;
							else if (Float.toString(event_density).equals("NaN")) //indeterminate form (usually occurs when numerator, # events == 0)
								event_density = 0;
							het_density = ((big_het_count / region_pats.get(r).size()) / (float)(regions.get(r).y - regions.get(r).x));
							if (Float.toString(het_density).equals("Infinity")) //ditto
								het_density = -1;
							else if (Float.toString(het_density).equals("NaN")) //ditto
								het_density = 0;
							FileOps.appendToFile(outDir + "/" + dir.getName() + "/" + cluster_names[i] + ".csv", regions.get(r).x + "-" + regions.get(r).y + "," + (regions.get(r).y - regions.get(r).x) + "," + (recur.get(r) / (float)(load.length - 1)) + "," + big_count + "," + big_het_count + "," + event_density + "," + het_density + "\n");
						}
					} catch (Exception e) { e.printStackTrace(); }
				}				
			}
		}
	}
	
	public static int indexOf(Point param, ArrayList<Point> arr) {
		for (int i = 0; i<arr.size(); i++)
			if (arr.get(i).x==param.x && arr.get(i).y==param.y)
				return i;
		return -1;
	}
	
	public static int indexOf(Gene gene, ArrayList<Gene> arr) {
		for (int i = 0; i<arr.size(); i++)
			if (arr.get(i).lbl.equals(gene.lbl))
				return i;
		return -1;
	}
	
	/**
	 * Parse curated SNP calls and region segmentation data into BED files that can be uploaded to the UCSC genome browser.
	 * @param inDir curated SNP calls
	 * @param regions_inDir region segmentation data
	 */
	public static void genBrowserTracks(String inDir, String regions_inDir, String outDir) {
		File[] files = (new File(inDir)).listFiles();
		ArrayList<String> pats = new ArrayList<String>();
		for (File file : files)
			if (file.getName().indexOf("csv")!=-1)
				pats.add(file.getName());
		String load, chr;
		String[] split, r_split = null, germ_som = {"germline", "somatic"};
		String[][][][] toWrite = new String[22][cluster_names.length - 1][pats.size()][2]; //...[germline/somatic]
		int[][] big_min = new int[22][cluster_names.length - 1], big_max = new int[22][cluster_names.length - 1]; //[chromosome][cluster]
		int[][][] event_count = new int[22][cluster_names.length - 1][pats.size()];
		for (int x = 0; x<22; x++) { //chromosomes
			for (int y = 0; y<cluster_names.length - 1; y++) { //clusters
				big_min[x][y] = Integer.MAX_VALUE;
				big_max[x][y] = Integer.MIN_VALUE;
				for (int z = 0; z<pats.size(); z++) { //patients
					for (int w = 0; w<2; w++) //germline, somatic
						toWrite[x][y][z][w] = "";
					event_count[x][y][z] = 0;
				}
			}
		}
		int min_val, score;
		int[] rgb = new int[3];
		String blockStarts, blockSizes;
		boolean bool, inRegion = false;
		ArrayList<Integer> bstarrlist; int[] bstarr;
		for (int p = 0; p<pats.size(); p++) { //iterate through patients
			System.out.println(p + " of " + pats.size());
			load = FileOps.loadFromFile(inDir + "/" + pats.get(p));
			for (int k = 0; k<cluster_names.length - 1; k++) { //iterate through {dup, loh, roc-loh}
				for (int i = 1; i<=22; i++) { //iterate through chromosomes
					if (i<=21)
						chr = Integer.toString(i);
					else
						chr = "X";
					System.out.print(chr + " ");
						
					for (int w = 0; w<2; w++) { //iterate through germline, somatic
						bool = true;
						try {
							r_split = FileOps.loadFromFile(regions_inDir + "/chr" + chr + "/" + cluster_names[k] + ".csv").split(pats.get(p).replace(".csv", "") + ",")[1].split("\n")[0].split(",");
						} catch (Exception e) { bool = false; }
						if (bool) {
							for (int j = 0; j<r_split.length; j++) //iterate through regions
								toWrite[i-1][k][p][w] += "chr" + chr + " " + r_split[j].split("-")[0] + " " + r_split[j].split("-")[1] + " row" + toWrite[i-1][k][p][w].split("\n").length + " 0 + " + r_split[j].split("-")[0] + " " + r_split[j].split("-")[1] + " 205,201,201\n";
							split = load.split("\nchr" + chr + ",");
							for (int j = 1; j<split.length; j++) { //iterate through positions
								if (split[j].split(",")[6].equals(germ_som[w])) {
									if (Integer.parseInt(split[j].split(",")[7].split("\n")[0])==k) {
										if (Integer.parseInt(split[j].split(",")[0]) < big_min[i-1][k]) //look for overall chromosomal window start point
											big_min[i-1][k] = Integer.parseInt(split[j].split(",")[0]);
										if (Integer.parseInt(split[j].split(",")[0]) > big_max[i-1][k]) //look for overall chromosomal window end point
											big_max[i-1][k] = Integer.parseInt(split[j].split(",")[0]);
									}
									score = 0;
									if (Integer.parseInt(split[j].split(",")[7].split("\n")[0])==k) { //if point is part of cluster k
										score = 900; //used with grayscale BEDs (historical artifact)
										event_count[i-1][k][p]++;
										if (split[j].split("\n")[0].indexOf("nonsynonymous")!=-1) { //if nonsynonymous
											rgb[0] = 255; rgb[1] = 0; rgb[2] = 0; //red
										}
										else { //if synonymous
											rgb[0] = 0; rgb[1] = 100; rgb[2] = 0; //dark green
										}	
										
										if (split[j].split("\n")[0].indexOf("somatic")!=-1) {
											rgb[0] = 0; rgb[1] = 26; rgb[2] = 255; //blue
										}
									}
									else if (Integer.parseInt(split[j].split(",")[7].split("\n")[0])==cluster_names.length - 1) { //if point is part of HET ball
										score = 300; //used with grayscale BEDs (historical artifact)
										if (split[j].split("\n")[0].indexOf("nonsynonymous")!=-1) { //if nonsynonymous
											rgb[0] = 255; rgb[1] = 165; rgb[2] = 0; //orange
										}
										else { //if synonymous
											rgb[0] = 50; rgb[1] = 205; rgb[2] = 50; //lime green
										}
									}
									if (score!=0) //if point is part of cluster k or HET ball
										toWrite[i-1][k][p][w] += "chr" + chr + " " + split[j].split(",")[0] + " " + (Integer.parseInt(split[j].split(",")[0]) + 1) + " row" + toWrite[i-1][k][p][w].split("\n").length + " 0 + " + split[j].split(",")[0] + " " + (Integer.parseInt(split[j].split(",")[0]) + 1) + " " + rgb[0] + "," + rgb[1] + "," + rgb[2] + "\n";
								}
							}
						}
					}
				}
				System.out.println();
			}
		}
		System.out.println("writing to file...");
		String manifest;
		ArrayList<Integer> patsort, blacklist;
		int temp_max, temp_max_ind;
		for (int i = 1; i<=22; i++) { //chromosomes
			if (i<=21)
				chr = Integer.toString(i);
			else
				chr = "X";
			for (int k = 0; k<cluster_names.length - 1; k++) { //clusters
				manifest = "";
				patsort = new ArrayList<Integer>(); 
				blacklist = new ArrayList<Integer>();
				while (patsort.size() < pats.size()) { //sort patients by event count
					temp_max = Integer.MIN_VALUE;
					temp_max_ind = -1;
					for (int t = 0; t<pats.size(); t++) {
						if (Utils.indexOf(blacklist, t)==-1) {
							if (event_count[i-1][k][t] > temp_max) {
								temp_max = event_count[i-1][k][t];
								temp_max_ind = t;
							}
						}
					}
					patsort.add(temp_max_ind);
					blacklist.add(temp_max_ind);
				}
				for (int p : patsort) { //patients
					toWrite[i-1][k][p][0] += toWrite[i-1][k][p][1]; //combine germline and somatic files
					toWrite[i-1][k][p][0] = "browser hide all\n"
						+ "track name=\"" + pats.get(p).replace(".germline.csv", "") + "\" description=\" \" visibility=dense itemRgb=\"on\"\n" + toWrite[i-1][k][p][0];
					toWrite[i-1][k][p][0] = "browser position chr" + chr + ":" + big_min[i-1][k] + "-" + big_max[i-1][k] + "\n" + toWrite[i-1][k][p][0];
					FileOps.writeToFile(outDir + "/chr" + chr + "_" + cluster_names[k] + "/" + pats.get(p).replace(".germline", ""), toWrite[i-1][k][p][0]);
					manifest += "https://raw.github.com/sidreddy96/lohcate/master/" + inDir.split("/")[inDir.split("/").length - 2] + "/browser_tracks/chr" + chr + "_" + cluster_names[k] + "/" + pats.get(p).replace(".germline", "") + "\n"; //you'll have to deal with this line as you please
				}
				FileOps.writeToFile(outDir + "/chr" + chr + "_" + cluster_names[k] + ".txt", manifest);
			}
		}
	}
	
	/**
	 * Parse region score data into BED files that can be uploaded to the UCSC genome browser.
	 * @param inDir region score data
	 */
	public static void addScoreTracks(String inDir, String outDir) {
		String chr;
		String[] load;
		String[] col_names = {"recurrence", "het_density", "event_density"};
		int[] cols = {2, 5, 6}; //recurrence, exclusivity, variant density
		String[][][] toWrite = new String[22][cluster_names.length - 1][cols.length];
		int[][] big_min = new int[22][cluster_names.length - 1], big_max = new int[22][cluster_names.length - 1];
		for (int x = 0; x<toWrite.length; x++) {
			for (int y = 0; y<toWrite[x].length; y++) {
				big_min[x][y] = Integer.MAX_VALUE;
				big_max[x][y] = Integer.MIN_VALUE;
				for (int z = 0; z<toWrite[x][y].length; z++)
					toWrite[x][y][z] = "";
			}
		}
		float[] min = new float[cols.length], max = new float[cols.length];
		int r_ind;
		float temp_score;
		ArrayList<Point> regions;
		for (int i = 1; i<=22; i++) { //iterate through chromosomes
			if (i<=21)
				chr = Integer.toString(i);
			else
				chr = "X";
			System.out.println("chr" + chr);
			for (int k = 0; k<cluster_names.length - 1; k++) { //iterate through clusters
				System.out.print(k + " ");
				load = FileOps.loadFromFile(inDir + "/chr" + chr + "/" + cluster_names[k] + ".csv").split("\n");
				for (int c = 0; c<cols.length; c++) {
					min[c] = Float.MAX_VALUE;
					max[c] = Float.MIN_VALUE;
				}
				for (int j = 1; j<load.length; j++) { //iterate through regions
					for (int c = 0; c<cols.length; c++) { //iterate through cols:{recurrence, exclusivity, variant density}
						if (Float.parseFloat(load[j].split(",")[cols[c]]) >= 0) {
							if (Float.parseFloat(load[j].split(",")[cols[c]]) < min[c]) //grab min. col. value
								min[c] = Float.parseFloat(load[j].split(",")[cols[c]]);
							if (Float.parseFloat(load[j].split(",")[cols[c]]) > max[c]) //grab max. col. value
								max[c] = Float.parseFloat(load[j].split(",")[cols[c]]);
						}
					}
					if (Integer.parseInt(load[j].split(",")[0].split("-")[0]) < big_min[i-1][k]) //grab overall chromosomal start point
						big_min[i-1][k] = Integer.parseInt(load[j].split(",")[0].split("-")[0]);
					if (Integer.parseInt(load[j].split(",")[0].split("-")[0]) > big_max[i-1][k]) //grab overall chromosomal end point
						big_max[i-1][k] = Integer.parseInt(load[j].split(",")[0].split("-")[0]);
				}
				
				for (int c = 0; c<cols.length; c++) { //iterate through cols:{recurrence, exclusivity, variant density}
					regions = new ArrayList<Point>();
					for (int j = 1; j<load.length; j++) { //iterate through regions
						//System.out.println(j + " of " + load.length);
						if (regions.size()==0)
							regions.add(new Point(Integer.parseInt(load[j].split(",")[0].split("-")[0]), Integer.parseInt(load[j].split(",")[0].split("-")[1]), Float.parseFloat(load[j].split(",")[cols[c]])));
						else
							regions = splice(new Point(Integer.parseInt(load[j].split(",")[0].split("-")[0]), Integer.parseInt(load[j].split(",")[0].split("-")[1]), Float.parseFloat(load[j].split(",")[cols[c]])), regions);
						//|^| we need to 'splice' each region into the list 'regions', since our scoring system allows for overlaps between regions (and because at a given chromosomal point, we are only interested in considering the most recurrent region -- giving the score track the highest possible 'resolution' in the genome browser)
					}
					//System.out.println();
					
					for (Point region : regions) { //int j = 1; j<load.length; j++) //iterate through regions
						temp_score = region.score;
						if (temp_score==-1) //if region has density == Infinity (which happens when the region only contains 1 variant)
							temp_score = max[c]; //we don't want to modify region.score, so we have a temp var. take its place
						toWrite[i-1][k][c] += "chr" + chr + " " + region.x + " " + region.y + " row" + toWrite[i-1][k][c].split("\n").length + " " + Utils.normalize(temp_score, min[c], max[c], 200, 900) + " + " + region.x + " " + region.y + "\n";
					}
				}
			}
			System.out.println();
		}
		System.out.println("writing to file...");
		String manifest;
		for (int i = 1; i<=22; i++) { //chromosomes
			if (i<=21)
				chr = Integer.toString(i);
			else
				chr = "X";
			for (int k = 0; k<cluster_names.length - 1; k++) { //clusters
				manifest = "";
				for (int c = 0; c<cols.length; c++) { //cols:{recurrence, het density, variant density}
					if (toWrite[i-1][k][c].length() > 0) {
						toWrite[i-1][k][c] = "browser hide all\n"
							+ "track name=\"" + col_names[c] + "\" description=\" \" visibility=dense useScore=1\n" + toWrite[i-1][k][c];
						toWrite[i-1][k][c] = "browser position chr" + chr + ":" + big_min[i-1][k] + "-" + big_max[i-1][k] + "\n" + toWrite[i-1][k][c];
						FileOps.writeToFile(outDir + "/chr" + chr + "_" + cluster_names[k] + "/" + col_names[c] + ".csv", toWrite[i-1][k][c]);
						manifest += "https://raw.github.com/sidreddy96/lohcate/master/" + inDir.split("/")[inDir.split("/").length - 2] + "/score_tracks/chr" + chr + "_" + cluster_names[k] + "/" + col_names[c] + ".csv" + "\n";
					}
				}
				FileOps.writeToFile(outDir + "/chr" + chr + "_" + cluster_names[k] + ".txt", manifest);
			}
		}
	}
	
	/**
	 * helper method for addScoreTracks()
	 */
	public static ArrayList<Point> splice(Point r, ArrayList<Point> arr) {
		ArrayList<Point> rtn = new ArrayList<Point>(), toAdd = new ArrayList<Point>();
		toAdd.add(r);
		for (Point elem : arr) { //iterate through 'template' arr
			if (intersecting(elem, r)) { 
				if (r.score >= elem.score) { //if takeover
					if (r.x - elem.x > 0)
						rtn.add(new Point(elem.x, r.x, elem.score)); //add non-zero left region
					if (elem.y - r.y > 0)
						rtn.add(new Point(r.y, elem.y, elem.score)); //add non-zero right region
				}
				else //ninad, what if a more highly recurrent sub-region (within the bounds of Point r) is already in List<Point> arr ??!! we wouldn't want to splice our less-recurrent, perhaps somewhat larger Point r into a slot that is already occupied by a more highly recurrent subregion
					toAdd = split(elem, toAdd); //this takes care of the above case
			}
			else
				rtn.add(elem);
		}
		for (Point elem : toAdd)
			rtn.add(elem);
		return rtn;
	}
	
	/**
	 * helper method for splice()
	 */
	public static ArrayList<Point> split(Point r, ArrayList<Point> arr) { //splice, w/out insertion at end (because the 'dominant' block fragment will already be present in our origional region list)
		ArrayList<Point> rtn = new ArrayList<Point>();
		for (Point elem : arr) {
			if (!intersecting(r, elem))
				rtn.add(elem);
			else {
				if (r.x - elem.x > 0)
					rtn.add(new Point(elem.x, r.x, elem.score)); //add non-zero left region
				if (elem.y - r.y > 0)
					rtn.add(new Point(r.y, elem.y, elem.score)); //add non-zero right region
			}
		}
		return rtn;
	}
	
	/**
	 * Returns true if the bounds of the two regions overlap.
	 */
	public static boolean intersecting(Point a, Point b) {
		Point smaller = a, bigger = b;
		if (bigger.y - bigger.x < smaller.y - smaller.x) {
			smaller = b;
			bigger = a;
		}
		return ((bigger.x <= smaller.x && smaller.x <= bigger.y) || (bigger.x <= smaller.y && smaller.y <= bigger.y));
	}
	
	/**
	 * Generate 'master' gene enrichment table (used to generate histograms).
	 * @param inDir curated SNP calls
	 */
	public static void getGeneEnrichment(String inDir, String outDir) {
		File[] files = (new File(inDir)).listFiles();
		String[] load;
		int ind;
		String name;
		ArrayList<Gene> genes = new ArrayList<Gene>(); 
		for (File file : files) { //iterate through curated SNP calls
			if (file.getName().indexOf("csv")!=-1) {
				System.out.println(file.getName());
				load = FileOps.loadFromFile(file.getAbsolutePath()).split("\n");
				for (int i = 1; i<load.length; i++) { //iterate through variants
					if (i%10000==0)
						System.out.println(i + " of " + load.length);
					name = load[i].split(",")[5];
					if (name.length() > 0) { //avoid ".", which crops up a lot
						if (name.indexOf("dist")!=-1) //gene names can sometimes come with an uninteresting/irrelevant prefix
							name = name.split("\\(")[0];
						ind = indexOf(new Gene(name, load[i].split(",")[0]), genes);
						if (ind==-1) { //if we haven't seen this gene before
							genes.add(new Gene(name, load[i].split(",")[0]));
							ind = genes.size() - 1;
							for (int k = 0; k<cluster_names.length; k++)
								genes.get(ind).patients.add(new ArrayList<String>());
						}
						if (Integer.parseInt(load[i].split(",")[1]) > genes.get(ind).max) //get right-bound of gene's range of variants
							genes.get(ind).max = Integer.parseInt(load[i].split(",")[1]); //...left-bound...
						if (Integer.parseInt(load[i].split(",")[1]) < genes.get(ind).min)
							genes.get(ind).min = Integer.parseInt(load[i].split(",")[1]);
						if (load[i].split(",")[6].indexOf("nonsynonymous")!=-1) //increment nonsynonymous variant count
							genes.get(ind).arr[0]++;
						else if (load[i].split(",")[6].indexOf("synonymous")!=-1) //...synonymous...
							genes.get(ind).arr[1]++;
						if (load[i].split(",")[7].equals("germline")) //...germline...
							genes.get(ind).arr[2]++;
						else if (load[i].split(",")[7].equals("somatic")) //...somatic...
							genes.get(ind).arr[3]++;
						if (Integer.parseInt(load[i].split(",")[8]) < cluster_names.length && Integer.parseInt(load[i].split(",")[8]) >= 0)
							genes.get(ind).counts[Integer.parseInt(load[i].split(",")[8])]++; //increment LOH/DUP/&c. count
						if (Integer.parseInt(load[i].split(",")[8])==2) //include right-of-center LOH count in general LOH count
							genes.get(ind).counts[1]++;
						for (int k = 0; k<cluster_names.length; k++)
							if (Utils.indexOf(genes.get(ind).patients.get(k), file.getName())==-1 && Integer.parseInt(load[i].split(",")[8])==k)
								genes.get(ind).patients.get(k).add(file.getName());
					}
				}
			}
		}
		String header =  "chr,range,gene,nonsynonymous,synonymous,germline,somatic,";
		for (int i = 0; i<cluster_names.length; i++)
			header += cluster_names[i] + ",";
		for (int i = 0; i<cluster_names.length; i++)
			header += cluster_names[i] + "_recurrence,";
		FileOps.writeToFile(outDir, header.substring(0, header.length() - 1) + "\n");
		for (Gene gene : genes)
			FileOps.appendToFile(outDir, gene.chr + "," + gene.min + "-" + gene.max + "," + gene.lbl + "," + gene.toString() + "\n");
	}
	
	public static void main(String[] args) {
		long sys_time_init = System.currentTimeMillis();
		
		String root = args[0]; //project directory
		switch (Integer.parseInt(args[1])) { //args[1] --> 'switchboard' parameter
			case 0:
				curateSNPCalls(root + "/naf-taf-inputs", root + "/snps", Integer.parseInt(args[2])); //args[2] --> 0::Illumina, 1::SOLiD
				break;
			case 1:
				segmentRegions(root + "/snps", root + "/regions");
				pickupSomatics(root + "/snps", root + "/regions", root + "/curated_snps");
				break;
			case 2:
				scoreRegions(root + "/regions", root + "/curated_snps", root + "/scored_regions");
				break;
			case 3:
				genBrowserTracks(root + "/curated_snps", root + "/regions", root + "/browser_tracks");
				addScoreTracks(root + "/scored_regions", root + "/score_tracks"); //come back to this
				break;
			case 4:
				getGeneEnrichment(root + "/curated_snps", root + "/gene_enrichment.csv");
				Enrichment.genTopGeneLists(root + "/gene_enrichment.csv", root + "/gene_enrichment");
				break;
			case 5:
				for (int i = 0; i<cluster_names.length - 1; i++)
					Enrichment.getPathwayEnrichment(root + "/gene_enrichment.csv", root + "/kegg_pathways_roster.tsv", root + "/kegg/pathway_enrichment/" + cluster_names[i] + ".csv", i);
				break;
			case 6:
				Enrichment.annotatePathways(root + "/gene_enrichment.csv", root + "/kegg_pathways_roster.tsv", root + "/KEGG");
				break;
			case 7:
				Enrichment.getGOTermCounts(root + "/gene_enrichment.csv", root + "/gene_association.goa_human", root + "/GO/counts/go_term_counts.csv");
				for (int i = 0; i<cluster_names.length - 1; i++) {
					Enrichment.getGOTermCounts(root + "/gene_enrichment/gene_enrichment_top_" + cluster_names[i] + ".csv", root + "/gene_association.goa_human", root + "/GO/counts/go_term_counts_top_" + cluster_names[i] + ".csv");
					Enrichment.getGOTermEnrichment(root + "/GO/counts/go_term_counts.csv", root + "/GO/counts/go_term_counts_top_" + Script.cluster_names[i] + ".csv", root + "/GO/enrichment/go_term_enrichment_top_" + Script.cluster_names[i] + ".csv");
					
					Enrichment.getGOTermCounts(root + "/gene_enrichment/g1/gene_enrichment_top_" + cluster_names[i] + ".csv", root + "/gene_association.goa_human", root + "/GO/g1/counts/go_term_counts_top_" + cluster_names[i] + ".csv");
					Enrichment.getGOTermEnrichment(root + "/GO/g1/counts/go_term_counts.csv", root + "/GO/g1/counts/go_term_counts_top_" + Script.cluster_names[i] + ".csv", root + "/GO/g1/enrichment/go_term_enrichment_top_" + Script.cluster_names[i] + ".csv");
				}
				break;
		}
		
		System.out.println("Time elapsed: " + (System.currentTimeMillis()-sys_time_init)/1000 + " seconds");
	}

}
