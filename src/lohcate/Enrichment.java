package lohcate;
import java.awt.Color;
import java.awt.image.BufferedImage;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;

import nutils.ArrayUtils;
import nutils.StringUtils;


import shared.*;
import enrichmentanalysis.Gene;
import enrichmentanalysis.HeatMap;
import enrichmentanalysis.Pathway;
import enrichmentanalysis.PathwayElement;
import enrichmentanalysis.Rectangle;

/**
 * LOHcate --- A software tool for LOH calling and visualization in cancer genomes
 * D Wheeler & SG Reddy
 * Human Genome Sequencing Center, Baylor College of Medicine (Houston, TX)
 * 
 * Runner class for enrichment analysis. Separate from main runner, Script.java, for historical reasons.
 * 
 * @author Siddharth G. Reddy
 *
 */
public class Enrichment {
	
	/**
	 * Generates table that lists enrichment of pathways in KEGG Pathways (all that have at least 1 gene in our list) with arbitrary figure
	 * @param inDir gene enrichment table
	 * @param kegg_pathway_gene_list_inDir KEGG Pathways gene::pathway roster table
	 */
	public static void getPathwayEnrichment(String inDir, String kegg_pathway_gene_list_inDir, String outDir, int clust) {
		
		String toWrite, roster_load = FileOps.loadFromFile(kegg_pathway_gene_list_inDir);
		String[] load = FileOps.loadFromFile(inDir).split("\n");
		String split = FileOps.loadFromFile(inDir);
		ArrayList<String> pathways = new ArrayList<String>(), pathways_grab = new ArrayList<String>();
		ArrayList<ArrayList<String>> path_genes = new ArrayList<ArrayList<String>>(); //pathways:genes
		int ind;
		//fill path_genes
		//String wlist = ",05200,04510,04310,04010,04340,04630,04060,04350,04370,04110,04210,04115,04330,";
		for (int i = 1; i<load.length; i++) { //iterate through genes
			if (i%1000==0)
				System.out.println(i + " of " + load.length);
			if (load[i].split(",")[2].length() > 1) { //avoid the "." gene name
				pathways_grab = getPathways(load[i].split(",")[2], roster_load); //grab pathways in which gene is implicated
				for (String pathway : pathways_grab) {
					//if (wlist.indexOf("," + pathway + ",")!=-1) {
						ind = ArrayUtils.linearSearch(pathways, pathway);
						if (ind==-1) {
							pathways.add(pathway);
							path_genes.add(new ArrayList<String>());
							path_genes.get(path_genes.size() - 1).add(load[i].split(",")[2]);
						}
						else
							path_genes.get(ind).add(load[i].split(",")[2]);
					//}
				}
			}
		}
		
		float tau, big_mean = 0f, mean;
		int g_p, g_x = 0;
		for (int i = 1; i<load.length; i++) //grab mean variant count
			big_mean += Float.valueOf(load[i].split(",")[7 + clust]); //there are 7 columns in the gene enrichment table before 'cluster' variant counts start getting listed
		big_mean /= (float)(load.length - 1);
		for (int i = 1; i<load.length; i++) //count # of genes with val. above mean
			if (Float.valueOf(load[i].split(",")[7 + clust]) >= big_mean)
				g_x++;
		FileOps.writeToFile(outDir, "pathway,kegg_id,tau,|g_p|,|P|,|g_x|,|X|\n");
		for (int i = 0; i<pathways.size(); i++) { //iterate through pathways
			System.out.println(i + " of " + pathways.size());
			g_p = 0;
			mean = 0f;
			for (int j = 1; j<path_genes.get(i).size(); j++) //grab mean variant count for genes in pathway
				mean += Float.valueOf(split.split("," + path_genes.get(i).get(j) + ",")[1].split(",")[4 + clust]); //we're splitting by gene name, so we don't have to account for the first three columsn in the table (hence the 7 --> 4)
			mean /= (float)path_genes.get(i).size();
			for (int j = 1; j<path_genes.get(i).size(); j++) //count # of genes in pathway with variant counts above mean (as calculated from general pool)
				if (Float.valueOf(split.split("," + path_genes.get(i).get(j) + ",")[1].split(",")[4 + clust]) >= big_mean)
					g_p++;
			tau = ((float)g_p / (float)g_x) * ((float)(load.length - 1) / (float)path_genes.get(i).size()); //fisher's exact test statistic
			try { //we can try to web-scrape the full name of the pathway, given that we know the 5-digit KEGG Pathways ID
				toWrite = getPathName(pathways.get(i)).replace(",", ";") + "," + pathways.get(i) + "," + tau + "," + g_p + "," + path_genes.get(i).size() + "," + g_x + "," + (load.length - 1) + "\n";
			} catch (Exception e) { 
				toWrite = "," + pathways.get(i) + "," + tau + "," + g_p + "," + path_genes.get(i).size() + "," + g_x + "," + (load.length - 1) + "\n";
			}
			FileOps.appendToFile(outDir, toWrite);
		}		
	}

	/**
	 * Web-scrapes full pathway name from KEGG Pathways, given a 5-digit ID (i.e. path == "05200" will return "Pathways in Cancer (Homo Sapiens)", or something like that)
	 */
	public static String getPathName(String path) throws IOException {
		String load = FileOps.getHTML("http://www.genome.jp/dbget-bin/www_bget?pathway+hsa" + path);
		return load.split("<nobr>Name</nobr>")[1].split("<br>")[0].split(">")[load.split("<nobr>Name</nobr>")[1].split("<br>")[0].split(">").length - 1];
	}

	/**
	 * Generates a 'preliminary' table of pathway/gene scores that are used to annotate KEGG Pathways gene-wiring diagrams. Exists for historical reasons.
	 * @param inDir gene enrichment table
	 * @param kegg_paths_inDir KEGG Pathways gene::pathway roster table
	 */
	public static void genPrelimTable(String inDir, String kegg_paths_inDir, String outDir, int clust) {
		String[] split = FileOps.loadFromFile(inDir).split("\n");
		String roster_load = FileOps.loadFromFile(kegg_paths_inDir), load = FileOps.loadFromFile(inDir);		
		ArrayList<String> pathways = new ArrayList<String>(), pathways_grab = new ArrayList<String>();
		ArrayList<ArrayList<String>> path_genes = new ArrayList<ArrayList<String>>(); //pathways:genes
		int ind;
		//String wlist = ",05200,";//",05200,04510,04310,04010,04340,04630,04060,04350,04370,04110,04210,04115,04330,";
		//fill path_genes
		for (int i = 1; i<split.length; i++) { //iterate through genes
			if (i%1000==0)
				System.out.println(i + " of " + split.length);
			pathways_grab = getPathways(split[i].split(",")[2], roster_load); //grab pathways in which gene is implicated
			for (String pathway : pathways_grab) {
				//if (wlist.indexOf("," + pathway + ",")!=-1) {
					ind = ArrayUtils.linearSearch(pathways, pathway);
					if (ind==-1) {
						pathways.add(pathway);
						path_genes.add(new ArrayList<String>());
						if (split[i].split(",")[2].length() > 1)
							path_genes.get(path_genes.size() - 1).add( (split[i].split(",")[2]).toUpperCase() );
					}
					else
						if (split[i].split(",")[2].length() > 1)
							path_genes.get(ind).add(Utils.gClean(split[i].split(",")[2]));
				//}
			}
		}
		
		//write to file
		String toWrite = "";
		FileOps.writeToFile(outDir, "");
		for (int i = 0; i<pathways.size(); i++) {
			System.out.println(i + " of " + pathways.size());
			toWrite = pathways.get(i) + ",";
			for (String elem : path_genes.get(i)) { //iterate through genes in pathway
				try {
					toWrite += elem + ";" + Utils.log(Float.parseFloat(load.split("," + elem + ",")[1].split("\n")[0].split(",")[4 + clust])) + "|" + load.split("," + elem + ",")[1].split("\n")[0].split(",")[4 + Enrichment.cluster_names.length - 1] + ",";
					//|^| I've tentatively added a feature to the gene-wiring diagram annotator that colors each gene with two colors, so the preliminary tables will have to include two 'scores' for each gene. But for now, I've commented out the downstream two-color option, so you will only see the first score manifest itself. Even so, I don't want to have to come back and change genPrelimTable(), so let's leave this the way it is.
				} catch (Exception e) { System.out.println("throwing " + elem); }
			}
			FileOps.appendToFile(outDir, toWrite.substring(0, toWrite.length() - 1) + "\n");
		}
	}
	
	/**
	 * Using a .tsv file downloaded from KEGG Pathways as a 'roster', grab all pathways that include the parameter-specified gene.
	 */
	public static ArrayList<String> getPathways(String gene, String kegg_pathway_gene_list) {
		ArrayList<String> pre = new ArrayList<String>();
		String[] split = kegg_pathway_gene_list.split("\t" + gene + "\t");
		for (int i = 0; i<split.length - 1; i++) { //iterate through instances of 'gene' in table
			pre.add(split[i].split("hsa")[split[i].split("hsa").length - 1]);
			try { //formatting
				pre.set(i, pre.get(i).split("\t")[0]);
			} catch (Exception e) { }
		}
		return pre;
	}
	
	/**
	 * Annotates gene-wiring diagrams from KEGG Pathways, given a prelim. table
	 * @param inDir prelim. table
	 */
	public static void annotateDiagrams(String inDir, String outDir) {
		String[] lines = FileOps.loadFromFile(inDir).split("\n");
		HashMap<String, ArrayList<Gene>> path_genes = new HashMap<String, ArrayList<Gene>>(); //<pathway, genes>
		Gene tmpgene;
		for (String line : lines) { //iterate through pathways
			path_genes.put(line.split(",")[0], new ArrayList<Gene>()); //add pathway
			for (int col = 1; col<line.split(",").length; col++) { //iterate through genes
				if (col%1000==0)
					System.out.println(col + " of " + line.split(",").length);
				tmpgene = new Gene(line.split(",")[col].split(";")[0]); //load gene name
				tmpgene.rank = Float.valueOf(line.split(",")[col].split(";")[1].split("\\|")[0]); //load score
				tmpgene.rank_2 = Float.valueOf(line.split(",")[col].split(";")[1].split("\\|")[1]); //load score 2
				if (tmpgene.rank!=0) //filter for blanks
					path_genes.get(line.split(",")[0]).add(tmpgene); //add gene to pathway
			}
		}
		genModDiagrams(path_genes, outDir); //now that we've compiled a list of pathways::genes, let's start writing to images!
	}
	
	/**
	 * Colours gene-wiring diagrams, given pathways::genes.
	 * @param path_genes pathway_KEGG_ID::genes
	 */
	public static void genModDiagrams(HashMap<String, ArrayList<Gene>> path_genes, String outDir) {
		HashMap<Rectangle, ArrayList<Gene>> mod_hash = new HashMap<Rectangle, ArrayList<Gene>>(); //<rectangle, genes>
		ArrayList<Gene> genes;
		String pathway;
		Iterator big_it = path_genes.entrySet().iterator();
		int big_counter = 0, small_count;
		float[] tmp_arr = new float[1];
		
		while (big_it.hasNext()) { //iterate through pathways
			
			System.out.println("pathway " + big_counter + " of " + path_genes.size()); big_counter++;
			
			Map.Entry big_pairs = (Map.Entry)big_it.next();
			genes = (ArrayList)big_pairs.getValue();
			pathway = (String)big_pairs.getKey();
			try {
				mod_hash = new HashMap<Rectangle, ArrayList<Gene>>(); //every bounded box in the KEGG Pathway diagram is a 'Rectangle', and may have multiple genes associated with it
				
				Pathway tmpath = new Pathway(pathway);
				tmpath.setDiagram(outDir + "/" + pathway + ".png", outDir + "/backend/" + pathway + ".html"); //grab pathway diagram+HTML image-map annotations from KEGG
				small_count = 0;
				for (Gene gene : genes) { //iterate through genes
					gene.setPathElem(tmpath, outDir + "/backend/" + pathway + ".html"); //set pathway diagram loc.
					small_count++;
				}
				
				HashMap<Rectangle, ArrayList<Gene>> temp_mod_hash = new HashMap<Rectangle, ArrayList<Gene>>();
				Iterator it;
				PathwayElement temp_2;
				Rectangle temp;
				boolean flag;
				for (Gene gene : genes) { //iterate through genes and place rectangles in hash
					for (Rectangle box : gene.getPathElem().getBoxes()) {
						flag = false;
						it = temp_mod_hash.entrySet().iterator();
						while (it.hasNext()) {
					        Map.Entry pairs = (Map.Entry)it.next();
					        temp = (Rectangle)pairs.getKey();
						    try {
						    	if (box.equals(temp)) {
						    		flag = true;
							        break;
							    }
						    } catch (Exception e) { }
						}
						if (!flag) //flagging prevents duplicates from entering the hash
							temp_mod_hash.put(box, new ArrayList<Gene>());
					}
				}
				
				it = temp_mod_hash.entrySet().iterator();
				int counter = 0;
				float score,max_pos = 0,max_neg = 0,min_pos = 0,min_neg = 0;
				float[] median, median_2;
				ArrayList<Float> pre_arr = new ArrayList<Float>(), pre_arr_2 = new ArrayList<Float>();
				float score_2 = 0,max_pos_2 = 0,max_neg_2 = 0,min_pos_2 = 0,min_neg_2 = 0;
				ArrayList<Gene> val, temp_arr;
				Color quartile_color = null;
				
				it = temp_mod_hash.entrySet().iterator();
			    while (it.hasNext()) { //iterate through rectangles
			        Map.Entry pairs = (Map.Entry)it.next();
			        temp = (Rectangle)pairs.getKey();
			        for (Gene gene : genes) { //iterate through genes
				        for (Rectangle box : gene.getPathElem().getBoxes()) {
					       	try {
						       	if (box.equals(temp)) { //if gene's pathway location is the same as the current iteration's box
						       		if ((gene.rank>max_pos && gene.rank>0) || (gene.rank!=0 && max_pos==0))
										max_pos = gene.rank;
									if ((Math.abs(gene.rank)>Math.abs(max_neg) && gene.rank<0) || (gene.rank!=0 && max_neg==0))
										max_neg = gene.rank;
									if ((gene.rank<min_pos && gene.rank>0) || (gene.rank!=0 && min_pos==0))
										min_pos = gene.rank;
									if ((Math.abs(gene.rank)<Math.abs(min_neg) && gene.rank<0) || (gene.rank!=0 && min_neg==0))
										min_neg = gene.rank;
									
									if ((gene.rank_2>max_pos_2 && gene.rank_2>0) || (gene.rank_2!=0 && max_pos_2==0))
										max_pos_2 = gene.rank_2;
									if ((Math.abs(gene.rank_2)>Math.abs(max_neg_2) && gene.rank_2<0) || (gene.rank_2!=0 && max_neg_2==0))
										max_neg_2 = gene.rank_2;
									if ((gene.rank_2<min_pos_2 && gene.rank_2>0) || (gene.rank_2!=0 && min_pos_2==0))
										min_pos_2 = gene.rank_2;
									if ((Math.abs(gene.rank_2)<Math.abs(min_neg_2) && gene.rank_2<0) || (gene.rank_2!=0 && min_neg_2==0))
										min_neg_2 = gene.rank_2;

						       		pre_arr.add(gene.rank);
						       		pre_arr_2.add(gene.rank_2);
						       	}
					       	} catch (Exception e) { }
				        }
			        }
			    }
			    
			    median = new float[pre_arr.size()];
			    median_2 = new float[median.length];
			    for (int i = 0; i<pre_arr.size(); i++) {
			    	median[i] = pre_arr.get(i) + Math.abs(max_neg);
			    	median_2[i] = pre_arr_2.get(i) + Math.abs(max_neg_2);
			    }
			    Arrays.sort(median);
			    Arrays.sort(median_2);
				
				//System.out.println("upper-bound: " + (Utils.max(max_pos, min_neg) + Math.abs(max_neg)) + " || median: " + median[median.length / 2] + " || lower-bound: 0");// + Utils.min(max_neg, min_pos));
				//System.out.println("2::upper-bound: " + (Utils.max(max_pos_2, min_neg_2) + Math.abs(max_neg_2)) + " || median: " + median_2[median_2.length / 2] + " || lower-bound: 0");// + Utils.min(max_neg, min_pos));
				
				it = temp_mod_hash.entrySet().iterator();
			    while (it.hasNext()) { //iterate through rectangles
			        Map.Entry pairs = (Map.Entry)it.next();
			        temp = (Rectangle)pairs.getKey();
			        score = -1 * Float.MAX_VALUE;
			        score_2 = -1 * Float.MAX_VALUE;
			        for (Gene gene : genes) { //iterate through genes
				        for (Rectangle box : gene.getPathElem().getBoxes()) {
					       	try {
						       	if (box.equals(temp)) { //if gene's pathway location is the same as the current iteration's box
						       		temp_mod_hash.get(temp).add(gene);
						       		if (gene.rank>score) //grab max score
						       			score = gene.rank;
						       		if (gene.rank_2>score_2) //grab max score
						       			score_2 = gene.rank_2;
						       	}
					       	} catch (Exception e) { }
				        }
			        }
			        
			        Color high = Color.RED, low = Color.GREEN, mid = Color.GRAY;
			        
			        //we're shifting the values + value gradient to min=0
			        
			        temp_arr = temp_mod_hash.get(temp);
			        temp.setRank(score); //set rank
			        temp.setRank2(score_2);
		        	HeatMap hmap = new HeatMap();
		        	tmp_arr[0] = Math.max(max_pos, min_neg) + Math.abs(max_neg); //upper bound
		        	hmap.setGradientPos(high, tmp_arr);
		        	tmp_arr[0] = Math.min(max_neg, min_pos) + Math.abs(max_neg);//lower bound
		    		hmap.setGradientNeg(low, tmp_arr);
		    		hmap.setMidVal(median[median.length / 2]);
		    		hmap.setZeroColor(mid);
		    		
		    		HeatMap hmap_2 = new HeatMap();
		        	tmp_arr[0] = Math.max(max_pos_2, min_neg_2) + Math.abs(max_neg_2); //upper bound
		        	hmap_2.setGradientPos(low, tmp_arr);
		        	tmp_arr[0] = Math.min(max_neg_2, min_pos_2) + Math.abs(max_neg_2); //lower bound
		    		hmap_2.setGradientNeg(high, tmp_arr);
		    		hmap.setMidVal(median_2[median_2.length / 2]);
		    		hmap_2.setZeroColor(mid);
		        	if (score!=0) {
			        	try {
			        		System.out.println((score + Math.abs(max_neg)) + " | " + (score_2 + Math.abs(max_neg_2)));
			        		temp.setRankColor(hmap.getColor(score + Math.abs(max_neg), 0));
			        		temp.setRank2Color(hmap_2.getColor(score_2 + Math.abs(max_neg_2), 0));
			        	} catch (Exception e) { 
			        		e.printStackTrace();
			        		System.out.println(score);
			        		temp.setRankColor(Color.WHITE);
			        		temp.setRank2Color(Color.WHITE); 
			        	}
		        	}
			        else {
			        	temp.setRankColor(Color.WHITE);
			        	temp.setRank2Color(Color.WHITE);
			        }
		        	mod_hash.put(temp, temp_arr);
			        
			        counter++;
			    }
			    
			    Rectangle key;
			    it = mod_hash.entrySet().iterator();
			    int count = 0;
			    BufferedImage img = FileOps.loadImage(outDir + "/" + pathway + ".png");
			    while (it.hasNext()) { //iterate through <rectangle, genes>
			    	count++;
			        Map.Entry pairs = (Map.Entry)it.next();
			        key = (Rectangle)pairs.getKey();
			        val = (ArrayList)pairs.getValue();
			        if (mod_hash.get(key).size()>0) {
				        try { //export diagram to file
				        	img = key.writeToDiagram(img);
				        } catch (Exception e) { e.printStackTrace(); }
			    	}
			    }
			    FileOps.writeImageToFile(img, outDir + "/" + pathway + ".png");
			} catch (Exception e) { e.printStackTrace(); }
		}
		
		System.out.println("done.");
	}
	
	/**
	 * Runner method for pathway annotation.
	 * @param en_inDir gene enrichment table
	 * @param roster_inDir KEGG Pathways gene::pathway roster table
	 */
	public static void annotatePathways(String en_inDir, String roster_inDir, String outDir) {
		for (int i = 0; i<Enrichment.cluster_names.length - 1; i++) {
			genPrelimTable(en_inDir, roster_inDir, outDir + "/" + Enrichment.cluster_names[i] + ".csv", i);
			annotateDiagrams(outDir + "/" + Enrichment.cluster_names[i] + ".csv", outDir);
		}
	}
	
	/**
	 * Count # times each Gene Ontology (GO) term occurs in list of genes
	 * @param inDir gene enrichment table
	 * @param go_inDir gene::GO term association table
	 */
	public static void getGOTermCounts(String inDir, String go_inDir, String outDir) {
		String go_load = FileOps.loadFromFile(go_inDir), toWrite = "";
		String[] load = FileOps.loadFromFile(inDir).split("\n");
		ArrayList<String> lbls = new ArrayList<String>();
		ArrayList<Integer> counts = new ArrayList<Integer>();
		int ind;
		ArrayList<String> go_grab;
		for (int i = 1; i<load.length; i++) { //iterate through genes
			System.out.println(i + " of " + load.length);
			try {
				go_grab = getGO(load[i].split(",")[2], go_load);
				for (String grab : go_grab) {
					ind = ArrayUtils.linearSearch(lbls, grab);
					if (ind==-1) {
						lbls.add(grab);
						counts.add(1);
					}
					else
						counts.set(ind, counts.get(ind) + 1);
				}
			} catch (Exception e) { }
		}
		
		for (int i = 0; i<lbls.size(); i++) {
			System.out.println(i + " of " + lbls.size());
			try {
				toWrite += lbls.get(i) + "," + counts.get(i) + "\n";
			} catch (Exception e) { e.printStackTrace(); System.out.println("throwing " + lbls.get(i)); }
		}
		FileOps.writeToFile(outDir, toWrite);
	}
	
	/**
	 * Grabs GO terms associated with parameter-specified gene ('query')
	 */
	public static ArrayList<String> getGO(String query, String go_load) {
		ArrayList<String> rtn = new ArrayList<String>();
		String[] split = go_load.split("\t" + query + "\tGO:"), line, prev_line = null;
		String grab, pre_rtn;
		for (int i = 1; i<split.length; i++) { //iterate through GO terms (associated with 'query')
			line = split[i].split("\n")[0].split("\t");
			if (i>1)
				if (!line[1].equals(prev_line[1])) //avoid dupes
					rtn.add("GO:" + line[1]);
			else	
				rtn.add("GO:" + line[1]);
			prev_line = line;
		}
		return rtn;
	}
	
	/**
	 * Grab top 20 genes by different variant counts. g1 lists require that genes have a region of variance of size > 0 (since variant density is biased towards those div-by-zero => Infinity values)
	 * @param inDir gene enrichment table
	 */
	public static void genTopGeneLists(String inDir, String outDir) {
		String[] load = FileOps.loadFromFile(inDir).split("\n");
		ArrayList<Integer> top, g1top;
		float temp_max, temp_max_g1;
		int temp_max_ind, temp_max_ind_g1;
		for (int i = 0; i<Enrichment.cluster_names.length - 1; i++) { //iterate through {dup, loh}
			System.out.println(Enrichment.cluster_names[i]);
			//FileOps.writeToFile(outDir + "/gene_enrichment_top_" + Script.cluster_names[i] + ".csv", trim(i, load[0]) + "\n");
			//FileOps.writeToFile(outDir + "/g1/gene_enrichment_top_" + Script.cluster_names[i] + ".csv", trim(i, load[0]) + "\n");
			ArrayList<Integer> blacklist = new ArrayList<Integer>();
			top = new ArrayList<Integer>();
			g1top = new ArrayList<Integer>();
			ArrayList<Integer> g1blacklist = new ArrayList<Integer>();
			
			int baseColOfClusterTypes = 9;
			while (top.size() < 20) { //grab top 20 genes by cluster[i] count
				temp_max = Float.MIN_VALUE; 
				temp_max_g1 = Float.MIN_VALUE;
				temp_max_ind = -1; 
				temp_max_ind_g1 = -1;
				
				for (int row = 1; row < load.length; row++) { //iterate through genes
					String[] columns = load[row].split(StringUtils.TabStr);
					
					if (columns[4].equals(Script.MissingGeneNameValue)) { //avoid "." gene name
						if (ArrayUtils.linearSearch(blacklist, row)==-1 && Float.parseFloat(columns[baseColOfClusterTypes]) > temp_max) {
							temp_max = Float.parseFloat(columns[baseColOfClusterTypes]);
							temp_max_ind = row;
						}
						if (ArrayUtils.linearSearch(g1blacklist, row)==-1 && Float.parseFloat(columns[baseColOfClusterTypes]) > temp_max_g1
								&& Float.parseFloat(columns[baseColOfClusterTypes]) <= 1) {//!columns[1].split("-")[0].equals(columns[1].split("-")[1])) { //handle genes with variant range > 0
							temp_max_g1 = Float.parseFloat(columns[baseColOfClusterTypes]);
							temp_max_ind_g1 = row;
						}
					}
				}
				g1blacklist.add(temp_max_ind_g1);
				blacklist.add(temp_max_ind);
				top.add(temp_max_ind);
				g1top.add(temp_max_ind_g1);
			}			
			
			for (int g = 0; g<top.size(); g++) { //iterate through top genes
				try {
					FileOps.appendToFile(outDir + "/gene_enrichment_top_" + Enrichment.cluster_names[i] + ".csv", trim(i, load[top.get(g)]) + "\n");
					FileOps.appendToFile(outDir + "/g1/gene_enrichment_top_" + Enrichment.cluster_names[i] + ".csv", trim(i, load[g1top.get(g)]) + "\n");
				} catch (Exception e) { }//e.printStackTrace(); }
			}
		}
	}
	
	public static String trim(int cluster, String param) {
		String rtn = "";
		for (int i = 0; i<param.split(",").length - (2 * Enrichment.cluster_names.length); i++) //avoid recurrence/variant count columns
			rtn += param.split(",")[i] + ",";
		if (param.indexOf("chr,range,")!=-1)
			rtn += "event_count,het_count,recurrence\n";
		else
			rtn += param.split(",")[param.split(",").length - Enrichment.cluster_names.length - (Enrichment.cluster_names.length - cluster)] + "," //avoid the recurrence columns (there's 1 for each cluster) & move left [cluster_names.length - cluster] columns
					+ param.split(",")[param.split(",").length - Enrichment.cluster_names.length - 1] + "," //this method is a little more robust than just using [baseColOfClusterTypes], but it doesn't matter that much which indexing method you choose
					+ param.split(",")[param.split(",").length - (Enrichment.cluster_names.length - cluster)] + "\n";
		return rtn;
	}
	
	/**
	 * Compare GO term counts in general pool of genes and top 20 lists, and calculate fold change for each comparison.
	 * @param inDir term counts for general pool
	 * @param spec_inDir term counts for targeted pool
	 */
	public static void getGOTermEnrichment(String inDir, String spec_inDir, String outDir) {
		String[] spec_load = FileOps.loadFromFile(spec_inDir).split("\n");
		String toWrite = "term,enrichment\n", load = FileOps.loadFromFile(inDir);
		int count = 0;
		for (String row : spec_load) {
			System.out.println(count + " of " + spec_load.length); count++;
			try {
				toWrite += getGOText(row.split(",")[0]).replace(",", ";") + "," + (Float.parseFloat(row.split(",")[1]) / Float.parseFloat(load.split(row.split(",")[0] + ",")[1].split("\n")[0]));
			} catch (Exception e) { e.printStackTrace(); System.out.println("throwing " + row.split(",")[0]); }
		}
		FileOps.writeToFile(outDir, toWrite);
	}
	
	/**
	 * Web-scrape full name of GO term, given "GO:***...*" ('query')
	 */
	public static String getGOText(String query) throws Exception {
		String grab = FileOps.getHTML("http://www.ebi.ac.uk/QuickGO/GSearch?q=" + query);
		return StringUtils.removeLeadingWhitespace(grab.split("<span class=\"obsolete_false\">")[1].split("</span>")[0].replace("\n", ""));
	}

	public static String[] cluster_names = {"dup", "loh", "roc-loh", "het"}; //het needs to be the last element of cluster_names, but more elem.s can be added to the 'front' (as long as you handle them in the getClusters() method )
	
}
