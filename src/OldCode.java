import java.util.ArrayList;

import shared.FileOps;
import shared.Utils;


public class OldCode {

	/**
	 * Define 'contiguous' regions of LOH, given our curated SNP calls.
	 * @param inDir curated SNP calls
	 */
	/*
	public static void segmentRegionsOld(String inDir, String outDir) {
		File[] files = (new File(inDir)).listFiles();
		String[] split;
		String[][] toWrite = new String[22][Script.cluster_names.length - 1]; //[chromosome][cluster]
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
						chr = Integer.toString(i);
						System.out.print(chr + " ");
						split = load.split("\nchr" + chr + ",");
						for (int k = 0; k<Script.cluster_names.length - 1; k++) { //iterate through {dup, loh, roc-loh}
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
																								
								//if (Integer.parseInt(split[min_ind].split("\n")[0].split(",")[7])==Script.cluster_names.length - 1)
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
							FileOps.appendToFile(outDir + "/chr" + chr + "/" + Script.cluster_names[j] + ".csv", toWrite[i][j].substring(0, toWrite[i][j].length() - 1) + "\n");
						toWrite[i][j] = "";
					}
				}
			}
		}
	}
	*/
	
	
	/**
	 * Now that we have our regions of LOH, let's 'pick up' any somatic mutations that fall within these regions (we'll call them LOH too). somatics don't separate into nice sidelobes in allele-fraction plots.
	 * @param snps_inDir curated SNP calls
	 * @param reg_inDir region segmentation data
	 */
	/* DEPRECATED
	public static void pickupSomatics(String snps_inDir, String reg_inDir, String outDir) {
		File[] files = (new File(snps_inDir)).listFiles();				
		
		for (File file : files) { //iterate through patients
			if (file.getName().indexOf("csv")!=-1) {
				System.out.println(file.getName().replace(".csv", ""));
				
				String[] load = FileOps.loadFromFile(file.getAbsolutePath()).split("\n");
				
				FileOps.writeToFile(outDir + "/" + file.getName(), load[0] + "\n");
				
				for (int i = 1; i < load.length; i++) { //iterate through variants
					boolean useUnmodifiedRow = true;
					String[] columns = load[i].split(",");
					
					if (columns[8].equals("-1")) { //if somatic
						System.out.println(i + " of " + load.length);
						
						//find out if variant lies within region of continuous DUP/LOH...
						for (int k = 0; k<Script.cluster_names.length - 1; k++) { //iterate through {dup, loh, roc-loh}
							String[] reg_load = FileOps.loadFromFile(reg_inDir + "/chr" + columns[0].replace("chr", "") + "/" + Script.cluster_names[k] + ".csv").split("\n");
							
							boolean toBreak = false;
							for (int row = 0; row < reg_load.length; row++) { //iterate through patients
								String[] regions = reg_load[row].split(",");
								
								for (int col = 1; col < regions.length; col++) { //iterate through regions
									if (Integer.parseInt(regions[col].split("-")[0]) < Integer.parseInt(columns[1])
											&& Integer.parseInt(columns[1]) <= Integer.parseInt(regions[col].split("-")[1])) { //if point falls within region
										
										FileOps.appendToFile(outDir + "/" + file.getName(), load[i].replace(",-1\n", ",1\n"));
										useUnmodifiedRow = false;
										toBreak = true; //since the point has already been positively identified, we don't need to keep iterating through patients, regions
										break; //ditto
									}
								}
								if (toBreak)
									break;
							}
						}
					}
					if (useUnmodifiedRow) FileOps.appendToFile(outDir + "/" + file.getName(), load[i] + "\n");
				}
			}
		}
	}
	*/
	
	/**
	 * Score regions of 'contiguous' LOH based on recurrence across patients, variant density, &c. Key feature of this method is that it separates regions into their 'least common denominators', which allows for detection of regions that, while smaller than the ones defined by segmentRegions(), may be more recurrent across patients.
	 * @param regions_inDir region segmentation data
	 * @param snps_inDir curated SNP calls
	 */
	/*
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
				for (int i = 0; i<Script.cluster_names.length - 1; i++) { //iterate through {dup, loh, roc-loh}
					regions = new ArrayList<Point>(); //Point.x --> region start, Point.y --> region end
					recur = new ArrayList<Double>();
					region_pats = new ArrayList<ArrayList<String>>();
					try {
						load = FileOps.loadFromFile(dir.getAbsolutePath() + "/" + Script.cluster_names[i] + ".csv").split("\n");
						//calculate recurrence score
						for (int x = 0; x<load.length; x++) { //iterate through patients
							
							System.out.println(x + " of " + load.length);
							String columns[] = load[x].split(",");
							
							for (int l = 1; l < columns.length; l++) { // iterate
																		// through
																		// regions
								xth = new Point(
										Integer.parseInt(columns[l].split("-")[0]),
										Integer.parseInt(columns[l].split("-")[1]));
								big_ind = indexOf(xth, regions);
								if (big_ind == -1) { // if region hasn't been
														// seen before
									regions.add(xth);
									recur.add(0.0);
									region_pats.add(new ArrayList<String>());
									region_pats.get(region_pats.size() - 1).add(columns[0]);
									big_ind = regions.size() - 1;
								}
								
								for (int y = 0; y<load.length; y++) { //iterate through patients, excluding x-th
									if (y != x) {
										String columnsOther[] = load[y].split(",");
										
										for (int m = 1; m<columnsOther.length; m++) { //iterate through y-th patient's regions
											yth = new Point(Integer.parseInt(columnsOther[m].split("-")[0]), Integer.parseInt(columnsOther[m].split("-")[1]));
											if (yth.x < xth.x && yth.y > xth.y) { //if region Y surrounds region X
												recur.set(big_ind, recur.get(big_ind) + 1.0); //increase region X's recurrence score
												region_pats.get(big_ind).add(columnsOther[0]);
												break;
											} else if (yth.x < xth.x && (xth.x < yth.y && yth.y < xth.y)) { //if region X slides to the right of region Y (they are still intersecting)
												ind = indexOf((new Point(xth.x, yth.y)), regions);
												if (ind!=-1) { //if we have seen region overlap(X, Y) before
													recur.set(ind, recur.get(ind) + 1.0);
													if (Utils.indexOf(region_pats.get(ind), columns[0])==-1)
														region_pats.get(ind).add(columns[0]);
													if (Utils.indexOf(region_pats.get(ind), columnsOther[0])==-1)
														region_pats.get(ind).add(columnsOther[0]);
												} else {
													regions.add(new Point(xth.x, yth.y));
													recur.add(1.0);
													ArrayList<String> newPats = new ArrayList<String>();
													region_pats.add(newPats);
													newPats.add(columns[0]);
													newPats.add(columnsOther[0]);
												}
												break;
											} else if (xth.x == yth.x && xth.y == yth.y) { //if regions X and Y have the same bounds
												recur.set(big_ind, recur.get(big_ind) + 1.0);
												region_pats.get(big_ind).add(columnsOther[0]);
												break;
											}
										}
									}
								}	
							}
						}
						//calculate variant/het density scores...
						FileOps.writeToFile(outDir + "/" + dir.getName() + "/" + Script.cluster_names[i] + ".csv", "region,size,recurrence,num_events,num_hets,event_density,het_density\n");//,exclusivity\n");
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
									else if (Integer.parseInt(pat_load[j].split("\n")[0].split(",")[7])==Script.cluster_names.length - 1) { //if point is part of HET ball
										if (Integer.parseInt(pat_load[j].split(",")[0]) >= regions.get(r).x && Integer.parseInt(pat_load[j].split(",")[0]) <= regions.get(r).y)
											big_het_count++;
									}
								}
							}
							int regionLength = (regions.get(r).y - regions.get(r).x) + 1;
							event_density = ((big_count / region_pats.get(r).size()) / (float) regionLength);
							if (Float.isInfinite(event_density)) {
								Utils.throwErrorAndExit("ERROR: Cannot have infinite density!");
							} else if (Float.isNaN(event_density)) {
								Utils.throwErrorAndExit("ERROR: Cannot have NaN density!");
							}

							het_density = ((big_het_count / region_pats.get(r).size()) / (float) regionLength);
							if (Float.isInfinite(het_density)) {
								Utils.throwErrorAndExit("ERROR: Cannot have infinite density!");
							} else if (Float.isNaN(het_density)) {
								Utils.throwErrorAndExit("ERROR: Cannot have NaN density!");
							}
							
							FileOps.appendToFile(outDir + "/" + dir.getName() + "/" + Script.cluster_names[i] + ".csv", regions.get(r).x + "-" + regions.get(r).y + "," + (regions.get(r).y - regions.get(r).x) + "," + (recur.get(r) / (float)(load.length - 1)) + "," + big_count + "," + big_het_count + "," + event_density + "," + het_density + "\n");
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
			if (arr.get(i).mLabel.equals(gene.mLabel))
				return i;
		return -1;
	}*/
	
	/**
	 * Parse curated SNP calls and region segmentation data into BED files that can be uploaded to the UCSC genome browser.
	 * @param inDir curated SNP calls
	 * @param regions_inDir region segmentation data
	 */
	/*
	public static void genBrowserTracks(String inDir, String regions_inDir, String outDir) {
		File[] files = (new File(inDir)).listFiles();
		ArrayList<String> pats = new ArrayList<String>();
		for (File file : files)
			if (file.getName().indexOf("csv")!=-1)
				pats.add(file.getName());
		String load, chr;
		String[] split, r_split = null, germ_som = {"germline", "somatic"};
		String[][][][] toWrite = new String[22][Script.cluster_names.length - 1][pats.size()][2]; //...[germline/somatic]
		int[][] big_min = new int[22][Script.cluster_names.length - 1], big_max = new int[22][Script.cluster_names.length - 1]; //[chromosome][cluster]
		int[][][] event_count = new int[22][Script.cluster_names.length - 1][pats.size()];
		for (int x = 0; x<22; x++) { //chromosomes
			for (int y = 0; y<Script.cluster_names.length - 1; y++) { //clusters
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
			for (int k = 0; k<Script.cluster_names.length - 1; k++) { //iterate through {dup, loh, roc-loh}
				for (int i = 1; i<=22; i++) { //iterate through chromosomes
					chr = Integer.toString(i);
					System.out.print(chr + " ");
						
					for (int w = 0; w<2; w++) { //iterate through germline, somatic
						bool = true;
						try {
							r_split = FileOps.loadFromFile(regions_inDir + "/chr" + chr + "/" + Script.cluster_names[k] + ".csv").split(pats.get(p).replace(".csv", "") + ",")[1].split("\n")[0].split(",");
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
										} else { //if synonymous
											rgb[0] = 0; rgb[1] = 100; rgb[2] = 0; //dark green
										}	
										
										if (split[j].split("\n")[0].indexOf("somatic")!=-1) {
											rgb[0] = 0; rgb[1] = 26; rgb[2] = 255; //blue
										}
									}
									else if (Integer.parseInt(split[j].split(",")[7].split("\n")[0])==Script.cluster_names.length - 1) { //if point is part of HET ball
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
			chr = Integer.toString(i);
			for (int k = 0; k<Script.cluster_names.length - 1; k++) { //clusters
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
					FileOps.writeToFile(outDir + "/chr" + chr + "_" + Script.cluster_names[k] + "/" + pats.get(p).replace(".germline", ""), toWrite[i-1][k][p][0]);
					manifest += "https://raw.github.com/sidreddy96/lohcate/master/" + inDir.split("/")[inDir.split("/").length - 2] + "/browser_tracks/chr" + chr + "_" + Script.cluster_names[k] + "/" + pats.get(p).replace(".germline", "") + "\n"; //you'll have to deal with this line as you please
				}
				FileOps.writeToFile(outDir + "/chr" + chr + "_" + Script.cluster_names[k] + ".txt", manifest);
			}
		}
	}*/
	
	
	
	/**
	 * Parse region score data into BED files that can be uploaded to the UCSC genome browser.
	 * @param inDir region score data
	 */
	public static void addScoreTracks(String inDir, String outDir) {
			
		String[] col_names = {"recurrence", "het_density", "event_density"};
		int[] cols = {2, 5, 6}; //recurrence, exclusivity, variant density
		String[][][] toWrite = new String[22][Script.cluster_names.length - 1][cols.length];
		
		int[][] big_min = new int[22][Script.cluster_names.length - 1],
				big_max = new int[22][Script.cluster_names.length - 1];
		
		for (int x = 0; x<toWrite.length; x++) {
			for (int y = 0; y<toWrite[x].length; y++) {
				big_min[x][y] = Integer.MAX_VALUE;
				big_max[x][y] = Integer.MIN_VALUE;
				for (int z = 0; z<toWrite[x][y].length; z++)
					toWrite[x][y][z] = "";
			}
		}
		
		float[] min = new float[cols.length], 
				max = new float[cols.length];
		
		int r_ind;
		float temp_score;
		ArrayList<Point> regions;
		for (int i = 1; i<=22; i++) { //iterate through chromosomes
			String chr = Integer.toString(i);
			
			System.out.println("chr" + chr);
			
			for (int k = 0; k<Script.cluster_names.length - 1; k++) { //iterate through clusters
				System.out.print(k + " ");
				String[] load = FileOps.loadFromFile(inDir + "/chr" + chr + "/" + Script.cluster_names[k] + ".csv").split("\n");
				for (int c = 0; c<cols.length; c++) {
					min[c] = Float.MAX_VALUE;
					max[c] = Float.MIN_VALUE;
				}
				
				for (int j = 1; j<load.length; j++) { //iterate through regions
					String[] regionCols = load[j].split(",");
					for (int c = 0; c<cols.length; c++) { //iterate through cols:{recurrence, exclusivity, variant density}
						
						if (Float.parseFloat(regionCols[cols[c]]) >= 0) {
							if (Float.parseFloat(regionCols[cols[c]]) < min[c]) //grab min. col. value
								min[c] = Float.parseFloat(regionCols[cols[c]]);
							
							if (Float.parseFloat(regionCols[cols[c]]) > max[c]) //grab max. col. value
								max[c] = Float.parseFloat(regionCols[cols[c]]);
						}
					}
					if (Integer.parseInt(regionCols[0].split("-")[0]) < big_min[i-1][k]) //grab overall chromosomal start point
						big_min[i-1][k] = Integer.parseInt(regionCols[0].split("-")[0]);
					if (Integer.parseInt(regionCols[0].split("-")[0]) > big_max[i-1][k]) //grab overall chromosomal end point
						big_max[i-1][k] = Integer.parseInt(regionCols[0].split("-")[0]);
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
			String chr = Integer.toString(i);
			
			for (int k = 0; k<Script.cluster_names.length - 1; k++) { //clusters
				manifest = "";
				for (int c = 0; c<cols.length; c++) { //cols:{recurrence, het density, variant density}
					if (toWrite[i-1][k][c].length() > 0) {
						toWrite[i-1][k][c] = "browser hide all\n"
							+ "track name=\"" + col_names[c] + "\" description=\" \" visibility=dense useScore=1\n" + toWrite[i-1][k][c];
						toWrite[i-1][k][c] = "browser position chr" + chr + ":" + big_min[i-1][k] + "-" + big_max[i-1][k] + "\n" + toWrite[i-1][k][c];
						FileOps.writeToFile(outDir + "/chr" + chr + "_" + Script.cluster_names[k] + "/" + col_names[c] + ".csv", toWrite[i-1][k][c]);
						manifest += "https://raw.github.com/sidreddy96/lohcate/master/" + inDir.split("/")[inDir.split("/").length - 2] + "/score_tracks/chr" + chr + "_" + Script.cluster_names[k] + "/" + col_names[c] + ".csv" + "\n";
					}
				}
				FileOps.writeToFile(outDir + "/chr" + chr + "_" + Script.cluster_names[k] + ".txt", manifest);
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
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub

	}

}
