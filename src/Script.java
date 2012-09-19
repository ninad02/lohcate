import genomeUtils.RegionRange;
import genomeUtils.RegionRange.RegionRangeOverlap;
import static genomeUtils.RegionRange.RegionRangeOverlap.*;

import java.awt.Color;
import java.awt.Font;
import java.awt.geom.Ellipse2D;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.Date;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYDotRenderer;
import org.jfree.chart.renderer.xy.XYItemRenderer;
import org.jfree.chart.renderer.xy.XYShapeRenderer;
import org.jfree.chart.title.LegendTitle;
import org.jfree.data.xy.DefaultXYDataset;
import org.jfree.data.xy.XYDataset;

import lohcateEnums.Chrom;
import lohcateEnums.ClusterType;
import lohcateEnums.ColorPastel;
import lohcateEnums.SNVType;
import lohcateEnums.SeqPlatform;
import lohcateEnums.VariantLocation;
import shared.BucketCounter;
import shared.FileOps;
import shared.GraphUtils;
import shared.IOUtils;
import shared.PrimitiveWrapper;
import shared.Utils;
import shared.Utils.FileExtensionAndDelimiter;

/**
 * LOHcate --- A software tool for LOH calling and visualization in cancer genomes
 * D Wheeler & SG Reddy
 * Human Genome Sequencing Center, Baylor College of Medicine (Houston, TX)
 * 
 * This runner class is the big enchilada. See main method for 'switchboard' / pipeline overview.
 * 
 * @author Siddharth G. Reddy, Ninad Dewal
 *
 */
public class Script {
	
	public static String[] cluster_names = {"dup", "loh", "roc-loh", "het"}; //het needs to be the last element of cluster_names, but more elem.s can be added to the 'front' (as long as you handle them in the getClusters() method )
	
	private static final float NAF_STRIP_EXPANDER = 2; //1.25f; //# of std. deviations to move away from the mean (when defining the thickness of the horizontal band containing HET ball, LOH sidelobes, &c.)
	private static final float HET_BALL_EPS = 0.035f, DUP_WEDGE_LASSO = 0.015f; //DBSCAN parameters for HET ball / DUP wedge detection
	private static final float NON_HET_BALL_EPS = 0.065f;
	private static final int NON_HET_BALL_MINPTS = 30;
	private static final int HET_BALL_MINPTS = 100, DUP_WEDGE_MINPTS = 100; //DBSCAN parameters for HET ball / DUP wedge detection
	private static final float ClusterDiagonalLeeway = (float) 0.2;
	public static final int DefaultDiploidCopyNumber = 2;
	
	private static final int REGION_SEGMENTATION_DIST_THRESHOLD = 2000000; //greatest possible distance between 2 'adjacent' points of LOH in a region of 'contiguous' LOH
	
	private static final String GermlineSuffix = "." + VariantLocation.Germline.toLowerCase();	
	private static final String NovelStr  = "novel";	
	private static final String ChromPrefix = "chr";
	public static final String MissingGeneNameValue = ".";
	private static final float MaxVariantAlleleFrequency = 1.0f;
	public static final String GenBrowserTrack = ".genBrowserTrack";

	// Column constants for the curated TSV files (files that have a cluster column)
	private static final int Col_NAFTAFInput_Chrom = 0;
	private static final int Col_NAFTAFInput_Position = 1;
	private static final int Col_NAFTAFInput_FlankingStringNormal = 5;
	private static final int Col_NAFTAFInput_FlankingStringTumor  = 6;
	private static final int Col_NAFTAFInput_TotalCoverageNormal  = 7;
	private static final int Col_NAFTAFInput_TotalCoverageTumor   = 8;
	private static final int Col_NAFTAFInput_VariantRatioNormal   = 11;
	private static final int Col_NAFTAFInput_VariantRatioTumor    = 12;
	private static final int Col_NAFTAFInput_DbSNPString    = 13;
	private static final int Col_NAFTAFInput_MutationType   = 14;
	private static final int Col_NAFTAFInput_HugoSymbol     = 15;
	
	private static final double GCContentThresholdLow  = 0.05;
	private static final double GCContentThresholdHigh = 0.80;
	
	// Column constants for the curated TSV files (files that have a cluster column)
	private static final int ColCuratedTSV_Chrom = 0;
	private static final int ColCuratedTSV_Position = 1;
	//private static final int ColCuratedTSV_VariantBaseTumor = 2;
	private static final int ColCuratedTSV_VafTumor = 3;
	private static final int ColCuratedTSV_Gene = 5;
	private static final int ColCuratedTSV_MutationType = 6;
	private static final int ColCuratedTSV_VariantLocation = 7;
	private static final int ColCuratedTSV_Cluster = 8;
	
	public static final GenomicCoordinateComparatorInTextFileLine LineComparatorTab = new GenomicCoordinateComparatorInTextFileLine();
	
	private static ArrayList<String> curateSNPCalls_removeHeaderLinesFromRows(ArrayList<String> rows) {		
		for (int rowIndex = 0; rowIndex < rows.size(); rowIndex++) {
			String row = rows.get(rowIndex);
			if (row.indexOf("refName") >= 0 && row.indexOf("coord") >= 0) {
				rows.set(rowIndex, null);
			} else {
				String gcString = Utils.extractNthColumnValue(row, Col_NAFTAFInput_FlankingStringTumor, Utils.FileExtensionTSV.mDelimiter); 
				double fractionGCNormal = Utils.calcFractionGC(gcString);
				if ((fractionGCNormal < GCContentThresholdLow) || (fractionGCNormal >= GCContentThresholdHigh)) {
					rows.set(rowIndex, null);
				}
			}
		}
		return Utils.removeNullElements(rows);		
	}
	
	// ========================================================================
	// INNER COMPARATOR CLASS
	// ========================================================================
	/** We define a comparator for sorting rows.  This assumes that the columns are tab
	 *  delimited and that the chromosome resides in the first column (in the form of an 
	 *  integer or chrN, where N represents an integer), while an integer genomic coordinate 
	 *  resides in the second column. 
	 */
	public static class GenomicCoordinateComparatorInTextFileLine implements Comparator<String> {
		public int compare(String line1, String line2) {
			long line1Coordinates = readChromNumAndBasePairPosition(line1);
			long line2Coordinates = readChromNumAndBasePairPosition(line2);
			return Utils.compareLong(line1Coordinates, line2Coordinates);			
		}
			
		/** This function reads the chromosome number and base pair position and packs
		 *  them into a long variable, with the chromosome number taking up the first 32
		 *  most significant bits and the base pair position taking up the last 32 bits.
		 *  This allows for returning multiple integer values and for easy sorting therafter
		 */
		private long readChromNumAndBasePairPosition(String line) {
			int firstTabIndexLine = line.indexOf(Utils.TabStr);			
			if (firstTabIndexLine < 0) Utils.throwErrorAndExit("ERROR: Columns are not tab delimited in input line 1!");

			int secondTabIndexLine = line.indexOf(Utils.TabStr, firstTabIndexLine + 1);
			if (secondTabIndexLine < 0) Utils.throwErrorAndExit("ERROR: Columns are not tab delimited in input line 1!");

			String chromString = "";
			int indexChromPrefix = line.indexOf(ChromPrefix);
			if (indexChromPrefix < 0) { // it does not exist
				chromString = line.substring(0, firstTabIndexLine);
			} else if (indexChromPrefix < firstTabIndexLine) {
				chromString = line.substring(indexChromPrefix + ChromPrefix.length(), firstTabIndexLine);				
			} else if (indexChromPrefix < secondTabIndexLine) {
				Utils.throwErrorAndExit("ERROR: " + ChromPrefix + " string in an incorrect place on line!");
			}  // We don't care if the string exists elsewhere
			
			int chromNum = Chrom.getChrom(chromString).getCode();					
			int basePairPosition = Integer.parseInt(line.substring(firstTabIndexLine + 1, secondTabIndexLine));
			//System.out.println(chromNum + "\t" + basePairPosition);
			return ( 0L | (((long) chromNum) << 32) | ((long) basePairPosition) );   // pack two values into one long	
		}
	}
	
	// ========================================================================
	// STAGE 1: Curate SNP Calls and Cluster with DBSCAN
	// ========================================================================
	
	/**
	 * Curate SNP calls by clustering data points into HET ball, DUP wedge, LOH sidelobes, &c., and by grabbing dbsnp population allele frequencies when possible
	 * @param inDir naf-taf-inputs
	 * @param opt 0::Illumina, 1::SOLiD
	 */
	public static void curateSNPCalls(String inDir, String outDir, String vafComparisonPlotDir, SeqPlatform platform) {
		File[] files = (new File(inDir)).listFiles();
		StringBuilder sb = new StringBuilder(8192);
		Utils.FileExtensionAndDelimiter fileExtDelim = Utils.FileExtensionTSV;		
		
		// Create output directory
		IOUtils.createDirectoryPath(outDir, false);
		IOUtils.createDirectoryPath(vafComparisonPlotDir, false);
		
		String[] columnHeaders = new String[] { "chr", "pos", "n_vaf", "t_vaf", "allele_freq", "gene", "mutation_type", "germ_som", "cluster" };
		String headerStr = Utils.constructColumnDelimitedString(columnHeaders, fileExtDelim.mDelimiter, sb, true).toString();		
				
		int fileIndex = 0;
		for (File file : files) {			
			int indexOfSubstring = file.getName().indexOf(GermlineSuffix);
			if (indexOfSubstring >= 0) {
															// && wList(file.getName())) {
				String samplenameRoot = file.getName().substring(0, indexOfSubstring);  				
				System.out.println("Processing (" + ++fileIndex + "): " + file.getName());
				
				String somaticFilename = file.getAbsolutePath().replace(VariantLocation.Germline.toLowerCase(), VariantLocation.Somatic.toLowerCase());
				ArrayList<String> somaticSpecificVariantRows  = FileOps.readAllLinesFromFile(somaticFilename);
				ArrayList<String> germlineSpecificVariantRows = FileOps.readAllLinesFromFile(file.getAbsolutePath());
				
				String headerStringGermline = germlineSpecificVariantRows.get(0);
				String headerStringSomatic  = somaticSpecificVariantRows.get(0);
				
				//upstream pipelines will randomly spit header lines into the 
				// middle of a naf-taf-input file. we're just avoiding those
				germlineSpecificVariantRows = curateSNPCalls_removeHeaderLinesFromRows(germlineSpecificVariantRows);
				somaticSpecificVariantRows  = curateSNPCalls_removeHeaderLinesFromRows(somaticSpecificVariantRows);
				
				// Create a combined set of rows
				ArrayList<String> allVariantRows = new ArrayList<String>(germlineSpecificVariantRows.size() + somaticSpecificVariantRows.size());
				allVariantRows.addAll(germlineSpecificVariantRows);
				int indexFirstSomaticRowInAllVariants = allVariantRows.size();  // save the size with just the germline variants added
				allVariantRows.addAll(somaticSpecificVariantRows);
				
				//Collections.sort(germlineSpecificVariantRows, LineComparatorTab);
				//Collections.sort(somaticSpecificVariantRows,  LineComparatorTab);
				Collections.sort(allVariantRows,              LineComparatorTab);
								
				String outFilename = samplenameRoot + fileExtDelim.mExtension;
				String outFilenameFullPath = outDir + File.separator + outFilename; 
				
				BufferedWriter out = IOUtils.getBufferedWriter(outFilenameFullPath);
				IOUtils.writeToBufferedWriter(out, headerStr, true);

				// Now get the clusters
				float[] copyNumRatios = getRoughCopyNumberRatioPerSite(allVariantRows, getAvgCoverageRatioPerChrom(allVariantRows));				
				
				ClusterType[] clustersGermline = getClusters(allVariantRows, copyNumRatios, outFilenameFullPath, indexFirstSomaticRowInAllVariants, platform); //get cluster assignments (HET ball, LOH sidelobes, DUP wedge, &c.)
				
//				PrintStream outStream = IOUtils.getPrintStream(outFilenameFullPath + ".withCopyNum.txt");
//				for (int r = 0; r < allVariantRows.size(); r++) {
				//	outStream.print(allVariantRows.get(r));
//					outStream.println("\t" + copyNumRatios[r] + "\t" + (copyNumRatios[r] * 2));
//				}
//				IOUtils.closePrintStream(outStream);

				
				System.out.println("Got clusters");
				
				// Now initialize the data structure needed to plot
				int[] clusterTypeCounts = Utils.getClusterTypeCounts(clustersGermline);
				double[][][] clusterCoordinates = new double[clusterTypeCounts.length][2][];
				for (int clusterTypeIndex = 0; clusterTypeIndex < clusterCoordinates.length; clusterTypeIndex++) {
					clusterCoordinates[clusterTypeIndex][0] = new double[ clusterTypeCounts[clusterTypeIndex] ];
					clusterCoordinates[clusterTypeIndex][1] = new double[ clusterTypeCounts[clusterTypeIndex] ];
				}
				int[] clusterTypeIndex = new int[clusterTypeCounts.length];
				Arrays.fill(clusterTypeIndex, -1);  
				
				// Do the post-processing
				int startingRowGermlineOrSomaticOrAll = 0;  // 0 because header line has been stripped away
				for (int i = startingRowGermlineOrSomaticOrAll; i < allVariantRows.size(); i++) {
					String strRow = allVariantRows.get(i);
					ClusterType clusterType = clustersGermline[i];
					int indexForClusterType = ++(clusterTypeIndex[clusterType.ordinal()]);  // increase the index of the data structure of the cluster type
										
					String[] germCols = strRow.split(Utils.TabStr);
										
					float vafNormal = extractVAFNormal(germCols, platform);
					float vafTumor  = extractVAFTumor (germCols, platform);
					clusterCoordinates[ clusterType.ordinal() ][0][indexForClusterType] = vafTumor;
					clusterCoordinates[ clusterType.ordinal() ][1][indexForClusterType] = vafNormal;
					
					String chromStr = (platform == SeqPlatform.Illumina) ? germCols[Col_NAFTAFInput_Chrom] : ChromPrefix + germCols[Col_NAFTAFInput_Chrom].replace(ChromPrefix, "");
					
					sb.setLength(0); // Clear the string builder
					sb.append(chromStr)
					  .append(fileExtDelim.mDelimiter).append(germCols[Col_NAFTAFInput_Position])
					  .append(fileExtDelim.mDelimiter).append(vafNormal)
					  .append(fileExtDelim.mDelimiter).append(vafTumor)
					  .append(fileExtDelim.mDelimiter); //chr,pos,n_vaf,t_vaf
					
					String variantAnnotation = "";
					if (platform == SeqPlatform.SOLiD) {
						variantAnnotation = (germCols[7].indexOf(NovelStr) >= 0) ? NovelStr : Utils.NAStr;
						// Sidd: For the NAStr case, strangely, the variant base is n/a in the SOLiD naf-taf-inputs 
						// (and there's not much point in looking up the reference base's allele frequency)
					} else if (platform == SeqPlatform.Illumina) {
						String dbsnpStr = germCols[Col_NAFTAFInput_DbSNPString]; 
						if (dbsnpStr.indexOf(NovelStr) >= 0) {
							variantAnnotation = NovelStr;
						} else {
							variantAnnotation = Utils.extractRsNumberFromLine(dbsnpStr);
							if ((variantAnnotation == null) || 
								 variantAnnotation.equals(Utils.rsPrefix) || 
								 variantAnnotation.equals(Utils.rsNegative)) {
								
								 variantAnnotation = Utils.NAStr;
							}
						}
					} else {
						variantAnnotation = Utils.NAStr;  // Could change if another platform added later
					}
					sb.append(variantAnnotation);
					
					//System.out.println("Got variant annotations");
					
					// Determine whether we are at a germline or somatic row
					String targetTissue = (i >= indexFirstSomaticRowInAllVariants) ? 
							VariantLocation.Somatic.toLowerCase() : VariantLocation.Germline.toLowerCase(); 
										
					String gene = "";
					String mutationType = "";
					
					if (platform == SeqPlatform.SOLiD) {
						mutationType = germCols[9].replace("syn", "synonymous").replace("nsynonymous", "nonsynonymous");
						mutationType = (mutationType == "") ? Utils.NAStr : mutationType;
						gene = germCols[8];
						gene = (gene == "") ? Utils.NAStr : gene;
					} else if (platform == SeqPlatform.Illumina) {
						gene         = germCols[Col_NAFTAFInput_HugoSymbol];
						mutationType = germCols[Col_NAFTAFInput_MutationType];
					}
					
					sb.append(fileExtDelim.mDelimiter).append(gene)
					  .append(fileExtDelim.mDelimiter).append(mutationType)
					  .append(fileExtDelim.mDelimiter).append(targetTissue)
					  .append(fileExtDelim.mDelimiter).append(clusterType)
					  .append(fileExtDelim.mDelimiter).append(copyNumRatios[i] * DefaultDiploidCopyNumber);
					
					
					IOUtils.writeToBufferedWriter(out, sb.toString(), true);
					IOUtils.flushBufferedWriter(out);
				}				
				IOUtils.closeBufferedWriter(out);
				
				// Now let's create the datasets needed to 
				
				DefaultXYDataset xyDataset = new DefaultXYDataset();
				for (ClusterType ct : ClusterType.values()) {
					double[][] xyDataSeries = clusterCoordinates[ct.ordinal()];
					xyDataset.addSeries(ct.name(), xyDataSeries);					
				}
				plotVAFComparison(xyDataset, vafComparisonPlotDir + File.separator + samplenameRoot, samplenameRoot);
				
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
	
	/** Plots the VAF (variant allele frequency) of the normal tissue comapred to the tumor tissue. */
	public static void plotVAFComparison(XYDataset xyDataset, String outFilenameRoot, String sampleName) {
		String title = "VAF Comparison: " + sampleName;
		String xAxisLabel = "VAF Tumor";
		String yAxisLabel = "VAF Normal";
		JFreeChart theChart = 
				ChartFactory.createScatterPlot(title, xAxisLabel, yAxisLabel, xyDataset, PlotOrientation.VERTICAL, true, false, false);
		XYPlot xyPlot = (XYPlot) theChart.getPlot();
		
		XYItemRenderer itemRenderer = getXYItemRendererHelper(5);		
		itemRenderer.setSeriesPaint(ClusterType.Dup.ordinal(), ColorPastel.Dark_Red.getColor());
		itemRenderer.setSeriesPaint(ClusterType.LOH.ordinal(), ColorPastel.Yellow_Green.getColor());
		itemRenderer.setSeriesPaint(ClusterType.HETGermline.ordinal(), ColorPastel.Dark_Yellow.getColor());
		itemRenderer.setSeriesPaint(ClusterType.HETSomatic.ordinal(), ColorPastel.Gray_80.getColor());
		itemRenderer.setSeriesPaint(ClusterType.Noise.ordinal(), ColorPastel.CMYK_Cyan.getColor());
		itemRenderer.setSeriesPaint(ClusterType.Null.ordinal(), ColorPastel.Light_Magenta_Red.getColor());		
		xyPlot.setRenderer(itemRenderer);		
		
		xyPlot.setBackgroundPaint(ColorPastel.Gray_15.getColor());
		xyPlot.setDomainGridlinePaint(Color.white);
		xyPlot.setRangeGridlinePaint(Color.white);		
		
		LegendTitle legendTitle = theChart.getLegend();
		legendTitle.setID("Clusters");
		legendTitle.setItemFont(new Font("Arial", Font.BOLD, 20));
		
		Font rangeAxisLabelFont = new Font("Arial", Font.BOLD, 20);
		Font rangeAxisTickFont = new Font("Arial", Font.BOLD, 20);
		xyPlot.getRangeAxis().setRange(0, 1.02);
		xyPlot.getDomainAxis().setRange(0, 1.02);
		xyPlot.getRangeAxis().setLabelFont(rangeAxisLabelFont);		
		xyPlot.getRangeAxis().setTickLabelFont(rangeAxisTickFont);
		xyPlot.getDomainAxis().setLabelFont(rangeAxisLabelFont);
		xyPlot.getDomainAxis().setTickLabelFont(rangeAxisTickFont);		
	
		// Now write the plot		
		int width = 800;
		int height = 800;		
		GraphUtils.saveChartAsPNG(outFilenameRoot + ".VAFComparison", theChart, width, height);
	}
	
	public static XYItemRenderer getXYItemRendererHelper(int size) {
		XYShapeRenderer xyShapeRend = new XYShapeRenderer();				
		
		//XYDotRenderer xyDotRend = new XYDotRenderer();
		Ellipse2D.Double theEllipse = new Ellipse2D.Double(0, 0, size, size);
		//xyDotRend.setBaseShape(theEllipse, true);
		xyShapeRend.setBaseShape(theEllipse);
		//xyDotRend.setDotWidth(size);
		//xyDotRend.setDotHeight(size);
		//return xyDotRend;
		return xyShapeRend;
	}

	// ========================================================================
	public static final float DefaultTumorNormalRatio = 1.0f;
	public static final float TumorNormalRatioOfSomaticSite = -1.0f;
	public static float[] getRoughCopyNumberRatioPerSite(ArrayList<String> rows, float[] ratioPerChrom) {
		float[] copyNumRatio    =   new float[rows.size()];
		boolean[] isSomaticSite = new boolean[rows.size()];
		
		String prevGene = "";
		int rowOfFirstInstanceOfGene = -1;
		
		float ratioSum = 0;
		int numRowsWithSameGene = 0;
		
		for (int row = 0; row < rows.size(); row++) {
			String line = rows.get(row);
			
			Chrom chrom =      Chrom.getChrom( Utils.extractNthColumnValue(line, Col_NAFTAFInput_Chrom,      Utils.FileExtensionTSV.mDelimiter) );
			String currentGene =               Utils.extractNthColumnValue(line, Col_NAFTAFInput_HugoSymbol, Utils.FileExtensionTSV.mDelimiter);
			
			int covgNormal = Integer.parseInt( Utils.extractNthColumnValue(line, Col_NAFTAFInput_TotalCoverageNormal, Utils.FileExtensionTSV.mDelimiter) );
			int covgTumor  = Integer.parseInt( Utils.extractNthColumnValue(line, Col_NAFTAFInput_TotalCoverageTumor,  Utils.FileExtensionTSV.mDelimiter) );
			
			isSomaticSite[row] = (covgNormal <= 0);
			float ratio = isSomaticSite[row] ? TumorNormalRatioOfSomaticSite : ((float) covgTumor / (float) covgNormal);
			
			if (currentGene.equalsIgnoreCase(prevGene)) {
				if (!isSomaticSite[row]) {
					numRowsWithSameGene++;
					ratioSum += ratio;
				} else {
					// We have a somatic site, we take no action in this block
				}
			} else {
				// New gene listed.  We must now write the ratios for the previous gene
				if (row > 0) {
					boolean onlySomaticRowsForGene = (numRowsWithSameGene <= 0);
					float averageRatio = onlySomaticRowsForGene ? TumorNormalRatioOfSomaticSite : 
						((numRowsWithSameGene == 1) ? DefaultTumorNormalRatio : 
													  ((ratioSum / (float) numRowsWithSameGene) / ratioPerChrom[chrom.ordinal()]));

					for (int i = rowOfFirstInstanceOfGene; i <= row; i++) {					
						copyNumRatio[i] = isSomaticSite[i] ? TumorNormalRatioOfSomaticSite : averageRatio;
					}
				}
				
				// New gene listed.  Thus, set a new row of first gene			
				rowOfFirstInstanceOfGene = row;
				prevGene = currentGene;
				if (isSomaticSite[row]) {
					ratioSum = numRowsWithSameGene = 0;
				} else {
					numRowsWithSameGene = 1;
					ratioSum = ratio;					 					
				}
			}			
		}
		
		return copyNumRatio;
	}
	
	// ========================================================================
	/** Tumor : Normal avg coverage ratio per chromosome. */
	public static float[] getAvgCoverageRatioPerChrom(ArrayList<String> rows) {			
		float[] ratios = new float[Chrom.values().length];
		Arrays.fill(ratios, DefaultTumorNormalRatio);  // Fill with 1.0 since that's an equal ratio
		
		int totalCoverageTumor = 0;
		int totalCoverageNormal = 0;
		Chrom prevChrom = null;
		
		for (int row = 0; row < rows.size(); row++) {
			String line = rows.get(row);
			
			Chrom chrom = Chrom.getChrom( Utils.extractNthColumnValue(line, Col_NAFTAFInput_Chrom, Utils.FileExtensionTSV.mDelimiter) );
			
			if ((prevChrom == null) || (chrom != prevChrom)) {
				if (prevChrom != null) {
					if (totalCoverageNormal > 0) {  // avoid div/0
						ratios[prevChrom.ordinal()] = (float) totalCoverageTumor / (float) totalCoverageNormal;
					}
				}
				
				totalCoverageTumor = totalCoverageNormal = 0;
				prevChrom = chrom;
			}
			
			int covgNormal = Integer.parseInt( Utils.extractNthColumnValue(line, Col_NAFTAFInput_TotalCoverageNormal, Utils.FileExtensionTSV.mDelimiter) );
			int covgTumor  = Integer.parseInt( Utils.extractNthColumnValue(line, Col_NAFTAFInput_TotalCoverageTumor,  Utils.FileExtensionTSV.mDelimiter) );
			
			totalCoverageNormal += covgNormal; 
			totalCoverageTumor  += covgTumor; 			
		}
		
		// We finished the loop without writing the ratio of the last chromosome
		if (totalCoverageNormal > 0) {
			ratios[prevChrom.ordinal()] = (float) totalCoverageTumor / (float) totalCoverageNormal;
		}
		
		return ratios;
	}
	
	
	
	
	// ========================================================================
	/*
	public static void graphDistProbabilities(double[][] distProb, String chartTitleStr, String outFilenamePrefix) {
		DefaultXYDataset xyDataset = new DefaultXYDataset();
		Dimension chartDimension = new Dimension(800, 800);
		xyDataset.addSeries("Graph",  distProb);		
			
		boolean showLegend = true;
		
		JFreeChart theChart = 
			ChartFactory.createScatterPlot(chartTitleStr, "X-coordinate", "Probability", 
					xyDataset, PlotOrientation.VERTICAL, showLegend, false, false);
		XYPlot xyPlot = (XYPlot) theChart.getPlot();
		//xyPlot.getDomainAxis().setTickLabelsVisible(false);
			
		XYDotRenderer xyDotRend = getXYDotRendererHelper(4);
		xyDotRend.setSeriesPaint(0, Color.blue);
		xyDotRend.setSeriesPaint(1, Color.green);
		xyPlot.setRenderer(xyDotRend);
		xyPlot.setBackgroundPaint(Color.white);
		xyPlot.setDomainGridlinePaint(Color.gray);
		xyPlot.setRangeGridlinePaint(Color.gray);
		
		// Set the chart title font
		TextTitle title = theChart.getTitle();
		Font titleFont = new Font("Arial", Font.BOLD, 44);
		title.setFont(titleFont);
		
		// Set the chart legend font
		LegendTitle legendTitle = theChart.getLegend();
		Font legendFont = new Font("Arial", Font.BOLD, 24);
		if (showLegend) { legendTitle.setItemFont(legendFont); }
		
		// Set the chart range axis font
		Font rangeAxisLabelFont = new Font("Arial", Font.BOLD, 32);
		Font rangeAxisTickFont = new Font("Arial", Font.BOLD, 22);
		xyPlot.getRangeAxis().setLabelFont(rangeAxisLabelFont);		
		xyPlot.getRangeAxis().setTickLabelFont(rangeAxisTickFont);
		xyPlot.getDomainAxis().setLabelFont(rangeAxisLabelFont);
		xyPlot.getDomainAxis().setTickLabelFont(rangeAxisTickFont);
	
		// 
		String outFilename = outFilenamePrefix + ".ChartDist.jpg"; 
		saveChartAsJPEG(outFilename, theChart, chartDimension.width, chartDimension.height);
		*/
	
	// ========================================================================
	// INNER CLASS
	// ========================================================================
	/** An inner class that calculates statistics for a set of sites (with variant 
	 *  allele frequencies information for tumor and matched normal). */
	public static class AlleleFrequencyStatsForSample {
		
		//NAF {het, loh, dup} FRAME definition via peak detection and parameter-tuned standard deviation expansion
		//we have to avoid the often hugely dense peak of homozygous mutations (AF > 0.8) and the occasionally hugely dense peak of neg. tail noise / somatics / &c. (AF < 0.2)
		float mVAFNormalFrameLower = 0.2f;
		float mVAFNormalFrameUpper = 0.8f; 
		float mBinSize             = 0.025f; //smoothing parameter
		int   mNumBins;
		int[]   mBinCount;    // The bins in which counts are binned and stored
		float[] mBinValue;    // The value each bin represents
		
		// Statistics values
		float mCountMean = -1;   
		float mVariance = -1;
		float mStdDev = -1;
		
		public AlleleFrequencyStatsForSample() {
			mNumBins = (int) ((mVAFNormalFrameUpper - mVAFNormalFrameLower) / mBinSize) + 1;
			mBinCount = new   int[mNumBins];
			mBinValue = new float[mNumBins];
			deduceBinValues();
			Arrays.fill(mBinCount, 0);			
		}
		
		public void tabulateAndPerformStatistics(ArrayList<String> rows, SeqPlatform platform) {			
			tallyVariantAlleleFrequenciesIntoBins(rows, platform, true);
			calculateSummaryStatistics();
		}
		
		private void deduceBinValues() {
			for (int i = 0; i < mBinValue.length; i++) {
				mBinValue[i] = mVAFNormalFrameLower + ((i + 1) * mBinSize);
			}
		}
		
		private void tallyVariantAlleleFrequenciesIntoBins(ArrayList<String> rows, SeqPlatform platform, boolean clearBinCountBins) {
			if (clearBinCountBins) {
				Arrays.fill(mBinCount, 0);
			}
			
			// First, tally the variant allele frequencies into bins
			for (int row = 0; row < rows.size(); row++) {
				float vafNormal = extractVAFNormal(rows.get(row), platform);
				if (Utils.inRangeLowerExclusive(vafNormal, mVAFNormalFrameLower, mVAFNormalFrameUpper)) {
					int binNumber = (int) ((vafNormal - mVAFNormalFrameLower) / mBinSize);
					mBinCount[binNumber]++;
				}
			}
		}
				
		private void calculateSummaryStatistics() {
			// Next, calculate summary statistics (mean, variance, standard deviation) 
			float total = 0;
			int numElementsTotal = 0;
						
			for (int i = 0; i < mBinCount.length; i++) {
				total            += (mBinCount[i] * mBinValue[i]);
				numElementsTotal +=  mBinCount[i];
			}	
					
			mCountMean = (float) total / (float) numElementsTotal;
			
			float varianceTotal = 0f;
			for (int i = 0; i < mBinValue.length; i++) { //calculate std. deviation of # points in each bin
				float diff = mCountMean - mBinValue[i];
				float diffSquared = diff * diff;
				float totalForBin = diffSquared * mBinCount[i];
				varianceTotal += totalForBin;			
			}
			
			mVariance = varianceTotal / (float) (numElementsTotal - 1);  // -1 for degrees of freedom		
			mStdDev = (float) Math.sqrt(mVariance);  //standard deviation in NAF-coord across {0.2 < NAF <= 0.8} variants in VAF plot
		}
	}
	
	
	/**
	 * A helper method for curateSNPCalls()
	 * @param load line-split FileOps.loadFromFile of naf-taf-input
	 * @param som_start start index of .somatic.txt naf-taf-input data in 'load' parameter
	 * This assumes header has been stripped out
	 * 
	 *  //Vv.vV well/poorly tuned parameters for different data sets
		//data set --> (NAF_STRIP_EXPANDER, (HET_BALL_EPS, HET_BALL_MINPTS), (DUP_WEDGE_LASSO, DUP_WEDGE_MINPTS))
		//target-aml --> (1, (0.035, 100), (0.015, 100))
		//target-all --> (1, (0.035, 100), (0.01, 100))
		//pvera --> (1.25, (0.035, 100), (0.015, 100))
		//hepato --> (1.25, (0.035, 100), (0.015, 100))
		//renal-we --> (1, (0.05, 500), (0.02, 350))
	 */
	public static ClusterType[] getClusters(ArrayList<String> rows, float[] copyNumRatios, String outFilenameFullPath, int startingRowGermlineOrSomaticOrAll, SeqPlatform platform) {
		
		// Get the allele frequency statistics, and adjust the frames based on the resulting standard deviation
		AlleleFrequencyStatsForSample afStatsSample = new AlleleFrequencyStatsForSample();
		afStatsSample.tabulateAndPerformStatistics(rows, platform);
		float vafNormalFrameAdjustedLower = afStatsSample.mCountMean - (NAF_STRIP_EXPANDER * afStatsSample.mStdDev);
		float vafNormalFrameAdjustedUpper = afStatsSample.mCountMean + (NAF_STRIP_EXPANDER * afStatsSample.mStdDev);
		
		//apply DBScan to points within NAF frame
		double scalingFactor = DefaultDiploidCopyNumber;
		//double scalingFactor = DefaultDiploidCopyNumber;
		ArrayList<Floint> points = getValidPointsListForClustering(rows, copyNumRatios, scalingFactor, platform, vafNormalFrameAdjustedLower, vafNormalFrameAdjustedUpper);
	
		System.out.println("Begin clustering algorithm: " + (new Date()).toString());
		
		//DBScanFaster dbscanner = new DBScanFaster(points, HET_BALL_EPS + DUP_WEDGE_LASSO, DUP_WEDGE_MINPTS, 0, 0, 1, 1); //parameters for capturing HET ball//KDBSCAN(param, 10, 0.0325f);
		DBScanFaster dbscanner = new DBScanFaster(points, HET_BALL_EPS, DUP_WEDGE_MINPTS * 2, 0, 0, 1, 1); //parameters for capturing HET ball//KDBSCAN(param, 10, 0.0325f);
		
		dbscanner.cluster();
		int clusterIDofHetBall = dbscanner.getCentralClusterID();
		Boolean[] pointsWithinRadius = dbscanner.getPointsWithinMinRadiusOfCluster(clusterIDofHetBall);
		int[] clusterAssignments = dbscanner.getClustAssignments();  // save and cache	
		
		int[] clusterTypeIDsFromAlgorithm = new int[ClusterType.values().length];		
		clusterTypeIDsFromAlgorithm[ClusterType.Dup.ordinal()]  = -1;
		clusterTypeIDsFromAlgorithm[ClusterType.Null.ordinal()] = -2;
		clusterTypeIDsFromAlgorithm[ClusterType.HETGermline.ordinal()] = clusterIDofHetBall;
		clusterTypeIDsFromAlgorithm[ClusterType.Noise.ordinal()] = dbscanner.getClusterIDOfNoise();
		
		
//		PrintStream outStream = IOUtils.getPrintStream(outFilenameFullPath + ".withCopyNum.txt");
//		for (int r = 0; r < rows.size(); r++) {
//			outStream.print(rows.get(r));
//			outStream.println("\t" + copyNumRatios[r] + "\t" + (copyNumRatios[r] * scalingFactor + "\t" + clusterAssign));
//		}
//		IOUtils.closePrintStream(outStream);
		
		// Now, re-run DBScan, but with changed parameters		
		//dbscanner.changeParams(HET_BALL_EPS + DUP_WEDGE_LASSO, DUP_WEDGE_MINPTS);		
		//dbscanner.cluster();
		//int clusterIDofHetBallWithWedge = dbscanner.getCentralClusterID();
		//int[] clusterAssignmentsWithWedge = dbscanner.getClustAssignments();

		// We now do a third pass, this time allowing only those points that are not part of the het cluster
		// We then perform a lower neighbor threshold on these points so that they can be declared as LOH
		ArrayList<Floint> nonHetPoints = new ArrayList<Floint>(points.size());
		boolean[] isNonHetPoint = new boolean[points.size()];		
		
		double threshold = 0.05;
		for (int i = 0; i < clusterAssignments.length; i++) {
			isNonHetPoint[i] = (clusterAssignments[i] != clusterIDofHetBall); 
			if (isNonHetPoint[i]) {
				nonHetPoints.add(points.get(i));
				if (points.get(i).mZ < scalingFactor - threshold) {
					clusterAssignments[i] = clusterTypeIDsFromAlgorithm[ClusterType.Dup.ordinal()];
				} else if (points.get(i).mZ > scalingFactor + threshold) {
					clusterAssignments[i] = 5;
				} else {
					clusterAssignments[i] = clusterIDofHetBall;
				}
			}
		}
		DBScanFaster dbscannerNonHet = new DBScanFaster(nonHetPoints, NON_HET_BALL_EPS, NON_HET_BALL_MINPTS, 0, 0, 1, 1);
		dbscannerNonHet.cluster();
		int[] clusterAssignmentsNonHet = dbscannerNonHet.getClustAssignments();
				
		System.out.println("End clustering algorithm: " + (new Date()).toString());
		
		int nonHetIndex = -1;
		for (int i = 0; i < clusterAssignments.length; i++) {
			
			if (isNonHetPoint[i]) {
				++nonHetIndex;
				
				if (points.get(i).mZ > scalingFactor) {
					//clusterAssignments[i] = clusterTypeIDsFromAlgorithm[ClusterType.Dup.ordinal()];
				} else {
					//clusterAssignments[i] = 5;
				}
				/*
				if (clusterAssignmentsNonHet[nonHetIndex] == clusterIDofHetBall) {
					// There is a collision of numeric cluster assignments.  This is due to the
					// fact that we run DBScan on different input lists, thus leading to respective
					// clusterIDs being assigned to each cluster for each run.  We need to ensure
					// that we are performing no collisions.  If the hetball cluster is the first
					// cluster, then we take the next cluster ID, else, we subract 1 from the hetball 
					// cluster ID.  Yes, this is kind of a hack.						
					clusterAssignmentsNonHet[nonHetIndex] = 
						(clusterIDofHetBall == DBSCAN2.ClusterIDOfNoise + 1) ? clusterIDofHetBall + 1 : clusterIDofHetBall - 1; 
				}
				//clusterAssignments[i] = clusterAssignmentsNonHet[nonHetIndex];  // assign for lower thresholds
				clusterAssignments[i] = clusterTypeIDsFromAlgorithm[ClusterType.Dup.ordinal()];  // assign for lower thresholds
				*/
			} else {	
				/*
				if (pointsWithinRadius[i] == null) {
					Utils.throwErrorAndExit("ERROR: Cannot possibly be a non-het point!");
				} else if (pointsWithinRadius[i] == Boolean.FALSE) {					
					if (pointOnDiagonal(points.get(i), ClusterDiagonalLeeway)) {
						clusterAssignments[i] = clusterTypeIDsFromAlgorithm[ClusterType.Null.ordinal()];
					} else {
						clusterAssignments[i] = clusterTypeIDsFromAlgorithm[ClusterType.Dup.ordinal()];   // only change if we're in a het ball region
					}
				}*/
			}
			
			/*
			if (isNonHetPoint[i]) {
				++nonHetIndex;
				if (clusterAssignmentsNonHet[nonHetIndex] == clusterIDofHetBall) {
					// There is a collision of numeric cluster assignments.  This is due to the
					// fact that we run DBScan on different input lists, thus leading to respective
					// clusterIDs being assigned to each cluster for each run.  We need to ensure
					// that we are performing no collisions.  If the hetball cluster is the first
					// cluster, then we take the next cluster ID, else, we subract 1 from the hetball 
					// cluster ID.  Yes, this is kind of a hack.						
					clusterAssignmentsNonHet[nonHetIndex] = 
						(clusterIDofHetBall == DBSCAN2.ClusterIDOfNoise + 1) ? clusterIDofHetBall + 1 : clusterIDofHetBall - 1; 
				}
				clusterAssignments[i] = clusterAssignmentsNonHet[nonHetIndex];  // assign for lower thresholds
			} else {								
				if (pointsWithinRadius[i] == null) {
					Utils.throwErrorAndExit("ERROR: Cannot possibly be a non-het point!");
				} else if (pointsWithinRadius[i] == Boolean.FALSE) {					
					if (pointOnDiagonal(points.get(i), ClusterDiagonalLeeway)) {
						clusterAssignments[i] = clusterTypeIDsFromAlgorithm[ClusterType.Null.ordinal()];
					} else {
						clusterAssignments[i] = clusterTypeIDsFromAlgorithm[ClusterType.Dup.ordinal()];   // only change if we're in a het ball region
					}
				}
			}
			*/
			
			/*
			if (isNonHetPoint[i]) { ++nonHetIndex; }  // make sure we're always at the correct index
						
			if (clusterAssignments[i] == DBSCAN2.ClusterIDOfNoise) {
				if (isNonHetPoint[i]) {					
					if (clusterAssignmentsNonHet[nonHetIndex] == clusterIDofHetBall) {
						// There is a collision of numeric cluster assignments.  This is due to the
						// fact that we run DBScan on different input lists, thus leading to respective
						// clusterIDs being assigned to each cluster for each run.  We need to ensure
						// that we are performing no collisions.  If the hetball cluster is the first
						// cluster, then we take the next cluster ID, else, we subract 1 from the hetball 
						// cluster ID.  Yes, this is kind of a hack.						
						clusterAssignmentsNonHet[nonHetIndex] = 
							(clusterIDofHetBall == DBSCAN2.ClusterIDOfNoise + 1) ? clusterIDofHetBall + 1 : clusterIDofHetBall - 1; 
					}
					clusterAssignments[i] = clusterAssignmentsNonHet[nonHetIndex];  // assign for lower thresholds
				} else {
					clusterAssignments[i] = clusterIDofDup;   // only change if we're in a het ball region
				}
			}
			*/		
		}

		// Now assign the cluster types
		
		return assignClusters(rows, platform, vafNormalFrameAdjustedLower, vafNormalFrameAdjustedUpper, clusterAssignments, clusterTypeIDsFromAlgorithm);
	}
	
	// ========================================================================
	// Convenience private function to break up caller function
	private static ClusterType[] assignClusters(ArrayList<String> rows, 
												SeqPlatform platform,
												//int startingRowSomatic,												
												float vafBoundLower, 
												float vafBoundUpper,
												int[] clusterAssignments,
												int[] clusterTypeIDsFromAlgorithm) {
		
		int indexInClusterAssignments = -1;  // we need to keep a special index, since not all rows are used.
		ClusterType[] returnClusters = new ClusterType[rows.size()];
		
		for (int row = 0; row < rows.size(); row++) {
			
			String line = rows.get(row);
			float vafNormal = extractVAFNormal(line, platform);
			boolean vafInRangeNormal = Utils.inRangeLowerExclusive(vafNormal, vafBoundLower, vafBoundUpper);
			
			if (/*(row >= startingRowSomatic) &&*/ (!vafInRangeNormal)) {
				// The vafNormal is either very low (homozygous reference) or very high (homozygous common variant).
				// We do some very simple decision making now (which should be replaced by formal clustering later)
				// to partition the calls.
				float justBelowZero = -0.0001f;				
				float hetBoundaryLower = 0.3333f;
				float hetBoundaryUpper = 0.6667f;
				float vafTumor = extractVAFTumor(line, platform);
				if (Utils.inRangeLowerExclusive(vafNormal, justBelowZero, vafBoundLower)) {
					// We are equal to or below the lower frame boundary
					if (vafTumor <= hetBoundaryLower) {
						// Normal: AA, Tumor: AA [Thus homozygous reference in both, no events]
						returnClusters[row] = ClusterType.Null;
					} else if (vafTumor > hetBoundaryUpper) {
						// Normal: AA, Tumor: BB or CC [made by: AA -> AB or AC (somatic het mutation) -> B or C (LOH, loss of A)]
						returnClusters[row] = ClusterType.LOH;
					} else {
						// Normal: AA, Tumor: AB [made by: AA -> AB (somatic het mutation)
						returnClusters[row] = ClusterType.HETSomatic;
					}					
				} else if (Utils.inRangeLowerExclusive(vafNormal, vafBoundUpper, MaxVariantAlleleFrequency)) {
					// We are above the upper frame boundary
					if (vafTumor <= hetBoundaryLower) {
						// Normal: BB, Tumor: AA [made by: BB -> AB (reverse somatic het mutation) -> A (LOH, loss of B)]
						returnClusters[row] = ClusterType.LOH;
					} else if (vafTumor > hetBoundaryUpper) {
						// Normal: BB, Tumor: BB or CC (ambiguous until we know exact variant for tumor)
						// TODO - Leave as Null for now, but will need to change later to resolve the
						// ambiguity mentioned above
						returnClusters[row] = ClusterType.Null;
					} else {
						// Normal: BB, Tumor: AB or CB [made by: BB -> AB (reverse somatic het mutation) or BB -> CB (somatic het mutation)
						returnClusters[row] = ClusterType.HETSomatic;
					}
				} else {
					Utils.throwErrorAndExit("ERROR: Contradiction - variant allele frequency cannot be in and out of bounds simultanteously!" + vafNormal);
				}
				
			} else {	// Our vaf-normal is in range in somatic sites, or we're at a germline site regardless of vaf-normal value	
				if (vafInRangeNormal) {
					++indexInClusterAssignments;
					final int assignedClusterID = clusterAssignments[indexInClusterAssignments]; // we create a local variable for fast test modifications 
					
					if (assignedClusterID        == clusterTypeIDsFromAlgorithm[ClusterType.HETGermline.ordinal()]) {
						returnClusters[row] = ClusterType.HETGermline; //HET
						
					} else if (assignedClusterID == clusterTypeIDsFromAlgorithm[ClusterType.Dup.ordinal()]) {
						returnClusters[row] = ClusterType.Dup; //DUP
						
					} else if (assignedClusterID == clusterTypeIDsFromAlgorithm[ClusterType.Null.ordinal()]) {
						returnClusters[row] = ClusterType.Null;  // we're on a diagonal
						
					} else if (assignedClusterID == clusterTypeIDsFromAlgorithm[ClusterType.Noise.ordinal()]) { 
						returnClusters[row] = ClusterType.Noise;
						
					} else { //anything not in the HET ball / DUP wedge is considered part of a LOH sidelobe
						//float vafTumor = extractVAFTumor(components, platform);						
						//returnClusters[row] = (vafTumor >= 0.5) ? ClusterType.LOHref : ClusterType.LOH; //LOH
						returnClusters[row] = ClusterType.LOH;
					}
					
				} else {
					returnClusters[row] = ClusterType.Null;   //outside NAF frame (<=> 'other')
				}
			}
		}
		
		return returnClusters;
	}
	
	// ========================================================================
	// Convenience function 
	private static ArrayList<Floint> 
		getValidPointsListForClustering(ArrayList<String> rows, float[] copyNumRatios, double scalingFactor, SeqPlatform platform, float vafBoundLower, float vafBoundUpper) {
		ArrayList<Floint> points = new ArrayList<Floint>(rows.size());
		
		for (int row = 0; row < rows.size(); row++) {
			String line = rows.get(row);			
			float vafNormal = extractVAFNormal(line, platform);
			if (Utils.inRangeLowerExclusive(vafNormal, vafBoundLower, vafBoundUpper)) {	
				Floint thePoint = new Floint(extractVAFTumor(line, platform), vafNormal, copyNumRatios[row] * (float) scalingFactor);
				points.add(thePoint);			
			}
		}
		
		return points;
	}
	
	// ========================================================================	
	// Tests whether the point lies roughly on the x = y diagonal, given a leeway from a ratio of 1.0
	private static boolean pointOnDiagonal(Floint point, double leeway) {
		double ratio = point.mY / point.mX;
		double equalRatio = 1.0;
		leeway = Math.min(leeway, equalRatio);
		return ( ((equalRatio - leeway) <= ratio) && (ratio <= (equalRatio + leeway)) );
	}
	
	// ========================================================================
	// Extracts and calculates the variant allele frequency depending on the platform 
	private static float extractVAFNormal(String[] components, SeqPlatform platform) {
		switch(platform) {
		case Illumina: return  Float.parseFloat(components[Col_NAFTAFInput_VariantRatioNormal]);
		case SOLiD:    return (Float.parseFloat(components[6]) / Float.parseFloat(components[5]));
		}
		return Float.NaN;
	}
	
	// ========================================================================
	// Extracts and calculates the variant allele frequency in the normal depending on the platform tab-delimited file
	private static float extractVAFNormal(String line, SeqPlatform platform) {
		switch(platform) {
		case Illumina: return  Float.parseFloat(Utils.extractNthColumnValue(line, Col_NAFTAFInput_VariantRatioNormal, Utils.TabStr));
		case SOLiD:    return (Float.parseFloat(Utils.extractNthColumnValue(line, 6, Utils.TabStr)) / 
				               Float.parseFloat(Utils.extractNthColumnValue(line, 5, Utils.TabStr)));
		}
		return Float.NaN;		
	}
	
	// ========================================================================
	private static float extractVAFTumor(String[] components, SeqPlatform platform) {
		switch(platform) {
		case Illumina: return  Float.parseFloat(components[Col_NAFTAFInput_VariantRatioTumor]);
		case SOLiD:    return (Float.parseFloat(components[4]) / Float.parseFloat(components[3]));
		}
		return Float.NaN;
	}
	
	// ========================================================================
	// Extracts and calculates the variant allele frequency in the tumor depending on the platform tab-delimited file
	private static float extractVAFTumor(String line, SeqPlatform platform) {
		switch(platform) {
		case Illumina: return  Float.parseFloat(Utils.extractNthColumnValue(line, Col_NAFTAFInput_VariantRatioTumor, Utils.TabStr));
		case SOLiD:    return (Float.parseFloat(Utils.extractNthColumnValue(line, 4, Utils.TabStr)) / 
				               Float.parseFloat(Utils.extractNthColumnValue(line, 3, Utils.TabStr)));
		}
		return Float.NaN;		
	}
	
	// ========================================================================
	// Stage 2: Segmentation
	// ========================================================================
	
	// ========================================================================
	/** 
	 * Define contiguous regions of copy number abberration (including LOH),
	 * given the results of our curated SNPs and DBScan clusters.
	 * @author Ninad Dewal (re-implementation or original by Siddharth Reddy)
	 * @param inDir curated SNP calls
	 * @param outDir segmentation results
	 */
	public static void segmentRegionsAllFiles(String inDir, String outDir, String outDirBrowserTracks) {
		Utils.FileExtensionAndDelimiter fileExtDelim = Utils.FileExtensionTSV;
		File[] files = (new File(inDir)).listFiles();
		ArrayList<CopyNumberRegionsByChromosome> regionsInSamples = new ArrayList<CopyNumberRegionsByChromosome>();
		ArrayList<File> validFilesList = new ArrayList<File>(files.length);
		
		// Create our output directories
		IOUtils.createDirectoryPath(outDir, false);
		IOUtils.createDirectoryPath(outDirBrowserTracks, false);
		
		// First we determine the regions from each sample
		for (File file : files) {
			CopyNumberRegionsByChromosome regionsInOneSample = segmentRegionsOneFile(file, outDir);
			if (regionsInOneSample != null) {
				regionsInSamples.add(regionsInOneSample);
				validFilesList.add(file);
			}
		}
		
		// Now we have all the contiguous regions from all the samples.  Find the regions of the cluster type
		
		// Declare the maximum stretch of a region for a particular cluster type 
		int[] maxBasePairsContiguousRegion = new int[] {
				REGION_SEGMENTATION_DIST_THRESHOLD,
				REGION_SEGMENTATION_DIST_THRESHOLD,
				Integer.MAX_VALUE
		};  
		
		ArrayList<ArrayList<CopyNumberRegionsByChromosome>> regionsInSamplesPerClusterType = new ArrayList<ArrayList<CopyNumberRegionsByChromosome>>();
		for (int clusterIndex = 0; clusterIndex < ClusterType.DupLOHHetG.length; clusterIndex++) {
			ClusterType clusterType = ClusterType.DupLOHHetG[clusterIndex];
			ArrayList<CopyNumberRegionsByChromosome> regionsInSamplesForOneClusterType = new ArrayList<CopyNumberRegionsByChromosome>();
			regionsInSamplesPerClusterType.add(regionsInSamplesForOneClusterType);			
			int maxBasePairsContiguousRegionForCluster = maxBasePairsContiguousRegion[clusterIndex];
			
			for (CopyNumberRegionsByChromosome regionsInOneSample : regionsInSamples) {
				CopyNumberRegionsByChromosome regionsInOneSampleMerged = mergeRegionsWithConstraints(regionsInOneSample, clusterType, maxBasePairsContiguousRegionForCluster);
				regionsInSamplesForOneClusterType.add(regionsInOneSampleMerged);
			}
		}

		PrintStream out = IOUtils.getPrintStream(outDir + File.separator + "testOut.txt");
		
		// Now we want to determine the recurrent regions
		ArrayList<CopyNumberRegionsByChromosome> regionsRecurrentPerClusterType = new ArrayList<CopyNumberRegionsByChromosome>();
		for (int clusterIndex = 0; clusterIndex < ClusterType.DupLOHHetG.length; clusterIndex++) {
			ClusterType clusterType = ClusterType.DupLOHHetG[clusterIndex];
			
			CopyNumberRegionsByChromosome recurrentRegionsForOneClusterType = 
					determineRecurrentRegions(regionsInSamplesPerClusterType.get(clusterIndex), clusterType);
			regionsRecurrentPerClusterType.add(recurrentRegionsForOneClusterType);
			recurrentRegionsForOneClusterType.print(out, fileExtDelim.mDelimiter);
		}
		
		
		System.out.println("Scoring regions...");
		// Now, we need to go through the recurrent regions and score them, based on the 
		// contents of the individual samples that make up the recurrent regions.
		for (int clusterIndex = 0; clusterIndex < ClusterType.DupLOHHetG.length; clusterIndex++) {
			ClusterType clusterType = ClusterType.DupLOHHetG[clusterIndex];

			CopyNumberRegionsByChromosome recurrentRegionsForOneClusterType = regionsRecurrentPerClusterType.get(clusterIndex);
			ArrayList<CopyNumberRegionsByChromosome> regionsInSamplesForOneClusterType = regionsInSamplesPerClusterType.get(clusterIndex);
			
			countClusterTypesInRegions(recurrentRegionsForOneClusterType, regionsInSamplesForOneClusterType, validFilesList);
			recurrentRegionsForOneClusterType.print(out, fileExtDelim.mDelimiter);
		}
				
		IOUtils.closePrintStream(out);

		System.out.println("Generating Browser tracks...");
		// Generatebrowser tracks
		for (int clusterIndex = 0; clusterIndex < ClusterType.DupLOH.length; clusterIndex++) {
			ClusterType clusterType = ClusterType.DupLOH[clusterIndex];
			ArrayList<CopyNumberRegionsByChromosome> regionsInSamplesForOneClusterType = regionsInSamplesPerClusterType.get(clusterIndex);
			CopyNumberRegionsByChromosome recurrentRegionsForOneClusterType = regionsRecurrentPerClusterType.get(clusterIndex);
			
			String outDirBrowserTracksForCluster = outDirBrowserTracks + File.separator + clusterType.name();
			IOUtils.createDirectoryPath(outDirBrowserTracksForCluster, false);
			genBrowserTracks(regionsInSamplesForOneClusterType, validFilesList, recurrentRegionsForOneClusterType, clusterType, outDirBrowserTracksForCluster);
		}
		
		
	}
	
	// ========================================================================
	/** Given a list of regions (call it list1) and a list of samplename files, this parses each sample.
	 *  For each site in the sample, it examines the cluster type assigned to the site.  If the site falls
	 *  within a region from list1, this increments list1's region's count of the site's cluster type. 
	 */
	public static void countClusterTypesInRegions(CopyNumberRegionsByChromosome regionsGenomeWide, 
												  ArrayList<CopyNumberRegionsByChromosome> regionsGenomeWideSampleSpecific, 
												  ArrayList<File> sampleFiles) {
		regionsGenomeWide.clearClusterCounts();
		
		StringBuilder sb = new StringBuilder(4096);

		int sampleIndex = -1;
		for (File inFile : sampleFiles) {
			
			// Load all lines into memory, extract the header row, and then sort by chrom/position
			ArrayList<String> allLines = IOUtils.readAllLinesFromFile(inFile.getAbsolutePath(), false, true, sb);
			String headerString = sb.toString();
			Collections.sort(allLines, LineComparatorTab);
						
			CopyNumberRegionsByChromosome regionsOneSample = regionsGenomeWideSampleSpecific.get(++sampleIndex);
			if (inFile.getName().indexOf(regionsOneSample.mSampleName) < 0) {
				Utils.throwErrorAndExit("ERROR: Samples don't match up: " + regionsOneSample.mSampleName + "\t" + inFile.getName());
			}
			
			// We initialize some indices for efficiency purposes
			Chrom prevChrom = null;
			ArrayList<CopyNumberRegionRange> regionsInChr = null;
			ArrayList<CopyNumberRegionRange> regionsInChrSampleSpecific = null;
			PrimitiveWrapper.WInteger regionIndexInChr               = new PrimitiveWrapper.WInteger(-1);
			PrimitiveWrapper.WInteger regionIndexInChrSampleSpecific = new PrimitiveWrapper.WInteger(-1);
			
			// Go through all the lines
			for (int lineIndex = 0; lineIndex < allLines.size(); lineIndex++) {
				String line = allLines.get(lineIndex);
				
				final Chrom chrom  = Chrom.getChrom                       (Utils.extractNthColumnValue(line, ColCuratedTSV_Chrom,    Utils.TabStr));					
				final int position = Integer.parseInt                     (Utils.extractNthColumnValue(line, ColCuratedTSV_Position, Utils.TabStr));
				final ClusterType clusterType = ClusterType.getClusterType(Utils.extractNthColumnValue(line, ColCuratedTSV_Cluster,  Utils.TabStr));
				
				if ((prevChrom == null) || (chrom != prevChrom)) {
					regionsInChr =                regionsGenomeWide.mRegionsByChrom.get(chrom.ordinal());
					regionsInChrSampleSpecific =   regionsOneSample.mRegionsByChrom.get(chrom.ordinal());					
					regionIndexInChr.mInt = 0;
					regionIndexInChrSampleSpecific.mInt = 0;

					prevChrom = chrom;  //  Must set
				}
				
				if (regionsInChr.isEmpty() || regionsInChrSampleSpecific.isEmpty()) {
					// no regions stored for this chromosome, so just move to the next line.  If
					// the next line for the sample is also the same chromosome, then we'll just
					// continue again, and so on... until we reach a line whose chromosome contains
					// a non-empty list in the passed in regions
					continue;  
				}
				
				// We have a list of regions in place, with an index set to >= 0.  We test 
				// whether the site position resides with the range, else we move on
				Boolean inRange = scanRegionsForPoint(regionsInChr, regionIndexInChr, chrom, position);		
				
				// Test our results
				if (inRange == null) {
					// This means that we could not find a region on this chromosome
					// that either: contains the site, or exists after the site.  The
					// site (and all other sites following the site), exist after the 
					// last region stored on this chromosome.  We just go through the
					// sites until we reach a site on the next chromosome
					continue;
				} else if (inRange == Boolean.FALSE) {
					// This means that the site lies prior to the next block in the
					// list.  We just need to go through the sites until we reach a
					// site that exists within the region.
					continue;
				} else {
					Boolean inRangeSampleSpecific = scanRegionsForPoint(regionsInChrSampleSpecific, regionIndexInChrSampleSpecific, chrom, position);
					if (inRangeSampleSpecific == null) {
						// We could find a region on the sample's chromosome that either:
						// contains the site, or exists downstream of the site.  The site
						// (and all other sites following the site), exist after the last 
						// region stored on this chromosome for the particular sample.  We
						// therefore take no action on this particular site.
					} else if (inRangeSampleSpecific == Boolean.FALSE) {
						// The site lies upstream of the next region in the sample's list
						// of regions.  We therefore take no action on this particular site.
					} else {
						CopyNumberRegionRange region = regionsInChr.get(regionIndexInChr.mInt);
						region.mClusterTypeCounts.increment(clusterType.ordinal());
					}
				}
			}
		}
	}

	// ========================================================================
	/** Given a sorted list of regions and an index at which to start, this scans across the list of 
	 *  regions to find the region that contains the given (chrom, position).  If we reach a region
	 *  that contains the given (chrom, position), we stop the scanning and have two return values:
	 *  the index of the region in the list, and Boolean.TRUE.  If we exhaust all regions and still 
	 *  have not reached the (chrom, position), it means that all regions exist upstream of the 
	 *  (chrom, position); we return two values: the index of the region in the list (which should 
	 *  be equal to the list size, as we have exhausted all blocks), and null.  Lastly, if we cannot 
	 *  find any region that contains the (chrom, position) but reach a region that lies downstream
	 *  of (chrom, position), we stop and return two values: Boolean.FALSE and the index of the 
	 *  downstream region.
	 * @param regionsInChromosome list of regions
	 * @param regionIndex index in which to start scanning.  This value is changed during list execution
	 * @param chrom 
	 * @param position
	 * @return Described above
	 */
	public static Boolean scanRegionsForPoint(ArrayList<? extends CopyNumberRegionRange> regionsInChromosome, 
											  PrimitiveWrapper.WInteger regionIndex, 
											  Chrom chrom, int position) {
		
		for (; regionIndex.mInt < regionsInChromosome.size(); regionIndex.mInt++) {			
			CopyNumberRegionRange region = regionsInChromosome.get(regionIndex.mInt);
			
			if (region.inRange(chrom, position)) {
				return Boolean.TRUE;				
			} else if (region.getRangeStart() > position) {
				return Boolean.FALSE;				
			}
		}
		return null;
	}
	
	// ========================================================================
	/** This takes in a list of regions per sample and intersects all the regions.  It reports
	 *  for each intersecting region the recurrence score across samples.
	 *   
	 */
	public static CopyNumberRegionsByChromosome determineRecurrentRegions(ArrayList<CopyNumberRegionsByChromosome> regionsInSamples, ClusterType clusterType) {
		
		// Test for a null parameter
		if (regionsInSamples.isEmpty()) return null;
		
		// Create a copy that serves as the target (or "sink") for all the intersecting regions
		CopyNumberRegionsByChromosome target = regionsInSamples.get(0).getDeepCopy();
		
		// Now traverse the rest of the samples
		for (int i = 1; i < regionsInSamples.size(); i++) {
			takeUnionAndBreakDownIntersectingRegions(target, regionsInSamples.get(i), clusterType);
		}
		
		return target;
	}
	
	// ========================================================================
	/** Given the segmented regions in one sample, this method merges the regions of one type (LOH, Amp, etc) 
	 *  in the following manner.  LOH/Amp regions cannot be longer than maxLengthContiguousRegion long. 
	 */
	public static CopyNumberRegionsByChromosome mergeRegionsWithConstraints(CopyNumberRegionsByChromosome regionsInSample, ClusterType clusterType, int maxBasePairsContiguousRegion) {
		
		// Create an empty return object
		CopyNumberRegionsByChromosome regionsInSampleMerged = new CopyNumberRegionsByChromosome(regionsInSample.mSampleName); 				
		
		// Go chromosome by chromosome
		for (Chrom chrom : Chrom.Autosomes) {				
			ArrayList<CopyNumberRegionRange> regionsInChromOriginal =       regionsInSample.mRegionsByChrom.get(chrom.getCode());
			ArrayList<CopyNumberRegionRange> regionsInChromMerged   = regionsInSampleMerged.mRegionsByChrom.get(chrom.getCode());
						
			// We declare a stored region that can be extended.  Initialize to null for now
			CopyNumberRegionRange regionToExtend = null;
			
			// Iterate through the regions for this chromosome
			for (CopyNumberRegionRange currentRegion : regionsInChromOriginal) {

				if (currentRegion.mCopyNumberClusterType == clusterType) {
					// Check if there's a region already waiting for extension.  
					// If not, create a new one (and a copy at that), and add to array
					if (regionToExtend == null) {
						regionToExtend = currentRegion.getCopy();
						regionsInChromMerged.add(regionToExtend);  // add this to the new array
						
					} else {
						if (currentRegion.mCopyNumberClusterType != regionToExtend.mCopyNumberClusterType) {
							Utils.throwErrorAndExit("ERROR: Must have same cluster type!\t" + currentRegion.mCopyNumberClusterType + "\t" + regionToExtend.mCopyNumberClusterType);
						}
						
						int maxEndIndexInclusive = regionToExtend.getRangeStart() + maxBasePairsContiguousRegion - 1;						
						if (currentRegion.getRangeStart() <= maxEndIndexInclusive) {
							regionToExtend.setRangeEnd(currentRegion.getRangeEnd());
						} else {
							// The current region is out of bounds.  We simply set
							// the current region as the new region to extend.
							regionToExtend = currentRegion.getCopy();
							regionsInChromMerged.add(regionToExtend);  // add this to the new array
						}
					}	
				}
			}		
		}
		
		return regionsInSampleMerged;
	}
	
	// ========================================================================
	/** The continuation of the @method segmentRegionsAllFiles method, but by individual. */ 
	public static CopyNumberRegionsByChromosome segmentRegionsOneFile(File inFile, String outDir) {
		FileExtensionAndDelimiter fileExtAndDelim = Utils.FileExtensionTSV;		
		
		// First check that the file is a file of the desired extension		
		int indexOfDelimiter = inFile.getName().indexOf(fileExtAndDelim.mExtension); 
		if (indexOfDelimiter < 0) return null;
		String samplenameRoot = inFile.getName().substring(0, indexOfDelimiter);
		
		StringBuilder sb = new StringBuilder(4096);
		
		// Load all lines into memory, extract the header row, and then sort by chrom/position
		ArrayList<String> allLines = IOUtils.readAllLinesFromFile(inFile.getAbsolutePath(), false, true, sb);
		String headerString = sb.toString();
		Collections.sort(allLines, LineComparatorTab);	
		
		// Have an array of regions for amplifications and LOH
		CopyNumberRegionsByChromosome regionsByChrom = new CopyNumberRegionsByChromosome(samplenameRoot);		 	
		CopyNumberRegionRange currentRegion = null;
		Chrom chromPreviousRow = null; // used to determine whether chrom has changed
		
		// We start at index 0 assuming no header and that the rows are sorted by chrom/position
		for (int row = 0; row < allLines.size(); row++) {			
			String line = allLines.get(row);			
			//String[] columns = line.split(Utils.TabPatternStr);
			
//			final Chrom chrom  = Chrom.getChrom  (columns[ColCuratedTSV_Chrom]);
//			final int position = Integer.parseInt(columns[ColCuratedTSV_Position]);
//			final ClusterType clusterType = ClusterType.getClusterType(columns[ColCuratedTSV_Cluster]);
			
			final Chrom chrom  = Chrom.getChrom                       (Utils.extractNthColumnValue(line, ColCuratedTSV_Chrom,    Utils.TabStr));					
			final int position = Integer.parseInt                     (Utils.extractNthColumnValue(line, ColCuratedTSV_Position, Utils.TabStr));
			final ClusterType clusterType = ClusterType.getClusterType(Utils.extractNthColumnValue(line, ColCuratedTSV_Cluster,  Utils.TabStr));
			
			// If the chromosome has changed, we set that we have no current region
			if ((chrom != chromPreviousRow) && (currentRegion != null)) {
				currentRegion.makeFinalized();
				currentRegion = null;
			}
			
			// Now determine whether we create or extend regions
			if (currentRegion == null) {
				if (segmentRegionsOneFile_isValidCluster(clusterType)) {
					currentRegion = new CopyNumberRegionRange(clusterType, chrom, position);
					regionsByChrom.addRegion(chrom, currentRegion);
				}
			} else {
				// Compare the position, make sure it's after the current range end
				if (position < currentRegion.getRangeEnd()) {
					Utils.throwErrorAndExit("ERROR: segmentRegionsOneFile(): Rows not sorted!");					
				} else if (position == currentRegion.getRangeEnd()) {
					System.out.println("WARNING: Duplicate coordinates: " + chrom + "\t" + position);
					continue;  // We ignore duplicate positions
				}
				
				// Now we know the position is after the current range end
				boolean isLOHCluster = (clusterType == ClusterType.LOH);
				if ((currentRegion.mCopyNumberClusterType == clusterType) && 
					(currentRegion.getChromosome() == chrom) &&
					(!isLOHCluster || (isLOHCluster && (position - currentRegion.getRangeEnd() < REGION_SEGMENTATION_DIST_THRESHOLD))) ) {
					boolean result = currentRegion.extendRange(chrom, position);
					if (!result) {
						Utils.throwErrorAndExit("ERROR: Could not extend range! " + currentRegion.toString() + "\t" + chrom + "\t" + position + "\t" + inFile.getName());
					}
				} else {
					// The cluster types are different, or the chrom didn't match, 
					// or the position was too far away.  Create a new region and add it
					// if it is not null or it is not noise.					
					if  (segmentRegionsOneFile_isValidCluster(clusterType)) {
						currentRegion.makeFinalized();
						currentRegion = new CopyNumberRegionRange(clusterType, chrom, position);
						regionsByChrom.addRegion(chrom, currentRegion);
					}
				}
			}
			
			chromPreviousRow = chrom;  // Save the chromosome of this row
		}
		
		if (currentRegion != null) currentRegion.makeFinalized();
		
		System.out.println(inFile.getName());				
		return regionsByChrom;
	}
	
	// ========================================================================
	private static boolean segmentRegionsOneFile_isValidCluster(ClusterType clusterType) {
		return (clusterType != ClusterType.Noise && clusterType != ClusterType.Null && clusterType != ClusterType.HETSomatic);
	}
	
	// ========================================================================
	/** Given an individual's copy number regions by chromosome, this partitions the
	 *  regions into three such objects based on the clustering type (Dup, LOH, HET)
	 */
	//public static ArrayList<CopyNumberRegionsByChromosome>
	
	// ========================================================================
	// INNER CLASS
	// ========================================================================
	public static class CopyNumberRegionRange extends RegionRange {
		public ClusterType mCopyNumberClusterType;
		public float mRecurrenceScore;
		protected BucketCounter mClusterTypeCounts;
		
		public CopyNumberRegionRange(ClusterType copyNumberClusterType, Chrom chrom, int regionStart) {
			super(chrom, regionStart);
			mCopyNumberClusterType = copyNumberClusterType;
			mRecurrenceScore = 1.0f;
			mClusterTypeCounts = new BucketCounter(ClusterType.getNumClusterTypes(), 0);
		}		
		
		public CopyNumberRegionRange(ClusterType copyNumberClusterType, Chrom chrom, int regionStart, int regionEnd) {
			super(chrom, regionStart, regionEnd);
			mCopyNumberClusterType = copyNumberClusterType;
			mRecurrenceScore = 1.0f;
			mClusterTypeCounts = new BucketCounter(ClusterType.getNumClusterTypes(), 0);
		}
		
		public CopyNumberRegionRange(CopyNumberRegionRange rhs) {
			super(rhs);
			mCopyNumberClusterType = rhs.mCopyNumberClusterType;
			mRecurrenceScore       = rhs.mRecurrenceScore;
			mClusterTypeCounts     = rhs.mClusterTypeCounts.getCopy();
		}
		
		public CopyNumberRegionRange getCopy() { return new CopyNumberRegionRange(this); }
	}
	
	// ========================================================================
	// INNER CLASS
	// ========================================================================
	/** Stores the regions for a given sample, split by chromosome.  
	 *  Chromosomes are indexed starting at 1 */
	public static class CopyNumberRegionsByChromosome {
		private static final ArrayList<CopyNumberRegionRange> dummyListForChrom0 = new ArrayList<CopyNumberRegionRange>();
		
		String mSampleName;
		ArrayList< ArrayList<CopyNumberRegionRange> > mRegionsByChrom;
		
		public CopyNumberRegionsByChromosome(String sampleName) {
			mRegionsByChrom = createRegionsByChromList();
			mSampleName = sampleName;
		}
		
		public static ArrayList<ArrayList<CopyNumberRegionRange>> createRegionsByChromList() {		
			ArrayList<ArrayList<CopyNumberRegionRange>> regionsByChrom = new ArrayList<ArrayList<CopyNumberRegionRange>>();
			
			// Create a dummy entry for chromosome 0
			for (int i = 0; i <= Chrom.values().length; i++) {
				regionsByChrom.add(new ArrayList<CopyNumberRegionRange>());
			}
			return regionsByChrom;		
		}
		
		public void addRegion(Chrom chrom, CopyNumberRegionRange region) {
			mRegionsByChrom.get(chrom.getCode()).add(region);
		}
		
		public void clearClusterCounts() {
			for (ArrayList<CopyNumberRegionRange> regionsInChrom : mRegionsByChrom) {
				for (CopyNumberRegionRange oneRegionInChrom : regionsInChrom) {
					oneRegionInChrom.mClusterTypeCounts.clear();
				}
			}
		}
		
		/** Returns an exact replica of this entire object, including copies of any contained regions. */ 
		public CopyNumberRegionsByChromosome getDeepCopy() {
			CopyNumberRegionsByChromosome newCopy = new CopyNumberRegionsByChromosome(mSampleName);
			
			for (int chromIndex = 0; chromIndex < mRegionsByChrom.size(); chromIndex++) {
				ArrayList<CopyNumberRegionRange> cnrrList    =         mRegionsByChrom.get(chromIndex);
				ArrayList<CopyNumberRegionRange> cnrrListNew = newCopy.mRegionsByChrom.get(chromIndex); 
				
				for (CopyNumberRegionRange cnrr : cnrrList) {
					cnrrListNew.add(cnrr.getCopy());
				}				
			}
			
			return newCopy;
		}
		
		/** Prints the contents to the PrintStream. */
		public void print(PrintStream out, String delimiter) {
			StringBuilder sb = new StringBuilder(65536);

			for (ArrayList<CopyNumberRegionRange> regionsOnChr : mRegionsByChrom) {
				for (CopyNumberRegionRange cnrr : regionsOnChr) {
					sb.setLength(0);
					sb.append(cnrr.mCopyNumberClusterType);
					sb.append(delimiter).append(cnrr.mRecurrenceScore);
					sb.append(delimiter).append(cnrr.toString());
					
					double densityClusterType =        ((double) cnrr.mClusterTypeCounts.getCount(cnrr.mCopyNumberClusterType.ordinal()) / (double) cnrr.getRangeLength());
					double densityClusterHetGermline = ((double) cnrr.mClusterTypeCounts.getCount(ClusterType.HETGermline.ordinal())     / (double) cnrr.getRangeLength());
					sb.append(delimiter).append(densityClusterType);
					sb.append(delimiter).append(densityClusterHetGermline);
					
					cnrr.mClusterTypeCounts.constructString(sb, false, delimiter);		
					out.println(sb.toString());
				}
			}
		}		
	}
	// ========================================================================
	// ========================================================================
	
	
	/** Given an "original" list of regions and an "added" list of regions, this
	 *  takes the union of the regions.  If there are overlaps, the regions
	 *  are broken into component pieces in which the overlapping sections 
	 *  are their own newfound regions, while the non-overlapping subsections
	 *  of the regions are created as newfound regions on their own.  For example:
	 *  
	 *  original:       xxxxxxxx         or             xxx
	 *  new:               xxx                       xxxxxxxx
	 *  resulting:      aaabbbcc                     aaabbbcc 
	 *  
	 *  The block of contiguous a's represent a region, the b's a region, and the
	 *  c's a region.  This is just one example.  There are more:
	 *  
	 *  original:          xxxxxxx        or         xxxxx   xxxxx   
	 *  added:          xxxxx   xxxxx                   xxxxxxx
	 *  resulting:      aaabbcccddeee                aaabbcccddeee
	 *  
	 *  original:         xxx xx xxx      or         xxxxxxxxxxxxxx
	 *  added:          xxxxxxxxxxxxxx                 xxx xx xxx
	 *  resulting:      aabbbcddefffgg               aabbbcddefffgg
	 *  
	 *  @param regionsTarget The list of original regions.  This parameter will
	 *  be modified by the function, so better to pass in a copy if the caller wants
	 *  to keep the original regions.
	 *  @param regionsAdded The list of regions to be added.
	 *  @return The new list of regions.  This, however, is just the original regions modified.
	 *  The function was designed this way (i.e. to modify the original regions) to increase
	 *  runtime efficiency -- namely, to avoid the memory calls to create a new list.  Updating
	 *  an already existing list is more efficient, and if the user wants the original to not
	 *  be modified, then he/she should pass in a copy.
	 */
	public static CopyNumberRegionsByChromosome
		takeUnionAndBreakDownIntersectingRegions(CopyNumberRegionsByChromosome regionsTarget, CopyNumberRegionsByChromosome regionsSource, final ClusterType clusterType) {
		
		// First iterate over the chromosomes
		// for (int chromIndex = 1; chromIndex < regionsTarget.mRegionsByChrom.size(); chromIndex++) {
		for (Chrom chrom : Chrom.Autosomes) {
			ArrayList<CopyNumberRegionRange> regionsChrTarget = regionsTarget.mRegionsByChrom.get(chrom.getCode());
			ArrayList<CopyNumberRegionRange> regionsChrSource = regionsSource.mRegionsByChrom.get(chrom.getCode());
			takeUnionAndBreakDownIntersectingRegions(regionsChrTarget, regionsChrSource, clusterType);	
		}

		return regionsTarget;
	}
	
	public static ArrayList<CopyNumberRegionRange>
		takeUnionAndBreakDownIntersectingRegions(ArrayList<CopyNumberRegionRange> regionsTarget, ArrayList<CopyNumberRegionRange> regionsSource, final ClusterType clusterType) {		
		
		// Make a pseudo-shallow array copy so we don't alter the caller's master copy 
		regionsSource = new ArrayList<CopyNumberRegionRange>(regionsSource);  
		
		int indexTarget = 0;
		int indexSource = 0;
		RegionRange.RegionRangeOverlap overlapTypePredicted = null;
		
		while ((indexTarget < regionsTarget.size()) && (indexSource < regionsSource.size())) {
			//System.out.println(indexTarget + "\t" + regionsTarget.size() + "\t" + indexSource + "\t" + regionsSource.size());
			CopyNumberRegionRange regionTarget = regionsTarget.get(indexTarget);
			CopyNumberRegionRange regionSource = regionsSource.get(indexSource);

			// Check and make sure that we are using the correct cluster type
			if (regionTarget.mCopyNumberClusterType != clusterType) {
				regionsTarget.remove(indexTarget);
				//indexTarget++;
				continue;
			} else if (regionSource.mCopyNumberClusterType != clusterType) {
				regionsSource.remove(indexSource);
				//indexSource++;
				continue;
			}
			
			// Test for what type of overlap we have
			RegionRange.RegionRangeOverlap overlapType = regionTarget.testAndCharacterizeOverlap(regionSource);			
			
			// Some of the following cases reduce down to other cases.  The following comparison
			// checks to make sure that the reduction producres for cases were done properly.  Some
			// cases are endpoints and do not result in reductions, however.
			if ((overlapTypePredicted != null) && (overlapTypePredicted != overlapType)) {
				Utils.throwErrorAndExit("ERROR: Predicted and determined region overlap types don't match!" 
						+ overlapTypePredicted + "\t" + overlapType);
			}
			
			//System.out.println(overlapType);
			switch(overlapType) {
			case Equals: 
				regionTarget.mRecurrenceScore += 1.0; //regionSource.mRecurrenceScore; 
				indexTarget++;
				indexSource++;
				overlapTypePredicted = null;
				break;
			case BeforeWithoutOverlap: case AdjacentBefore:
				// We leave the block in the target list as is, and we simply increment our target index
				indexTarget++;
				overlapTypePredicted = null;
				break;
			case AfterWithoutOverlap: case AdjacentAfter: {  // use braces to prevent declarations from spilling into other cases					
				regionsTarget.add(indexTarget, regionSource.getCopy());
				indexTarget++;  // To keep with the correct target region
				indexSource++;  // move to the next region to add
				overlapTypePredicted = null;
				break;
			}
			case BeforeViaDiffChromosome: case AfterViaDiffChromosome: 
				overlapTypePredicted = null;
				Utils.throwErrorAndExit("ERROR: Should be on different chromosomes!");					
				break;
			case SubsumesTotal: case SubsumesAlignedRight: case BeforeWithOverlap: {					
				breakdownThreeCasesHelper(regionsTarget, regionsSource, indexTarget, indexSource, false, true);
				indexTarget++;  // increment the target array index				

				// The following procedure reduces these cases to: 
				//    SubsumesAlignedLeft, if SubsumesTotal was this case
				//    Equals, if SubsumesAlignedRight was this case
				overlapTypePredicted = (overlapType == SubsumesTotal) ? SubsumesAlignedLeft : 
					((overlapType == BeforeWithOverlap) ? ConsumedByAlignedLeft : Equals);
				break;
			}
			case SubsumesAlignedLeft: {
				// The following procedure reduces the case to Equals, and whatever follows
				breakdownSubsumesConsumedAlignedLeftHelper(regionsTarget, regionsSource, indexTarget, indexSource, false);
				overlapTypePredicted = RegionRangeOverlap.Equals;
				break;
			}
			case ConsumedByTotal: case ConsumedByAlignedRight: case AfterWithOverlap: {
				// The following procedure reduces these cases to:
				//     ConsumedByAlignedLeft, if ConsumedByTotal was this case
				//     Equals, if ConsumedByAlignedRight was this case
				breakdownThreeCasesHelper(regionsSource, regionsTarget, indexSource, indexTarget, true, false);
				indexTarget++;  // increment the target array index
				
				overlapTypePredicted = (overlapType == ConsumedByTotal) ? ConsumedByAlignedLeft :
					((overlapType == AfterWithOverlap) ? SubsumesAlignedLeft : Equals);
				break;
			}
			case ConsumedByAlignedLeft: {
				// The following procedure reduces the case to Equals, and whatever follows
				breakdownSubsumesConsumedAlignedLeftHelper(regionsSource, regionsTarget, indexSource, indexTarget, true);
				overlapTypePredicted = RegionRangeOverlap.Equals;
				break;					
			}		
			default:
				Utils.throwErrorAndExit("ERROR: Invalid option!");
			}				
		}
				

		// Add any remaining regions.  We know that either of the following two loops will execute, but not both
		// since the previous loop was broken by failure of one of the loop conditions.
		// 
		// We do not need to add any regions to the target array, as either the array elements
		// were already traversed, or they already exist in the original array.  However, we need to remove
		// elements that may not match the cluster type desired
		while (indexTarget < regionsTarget.size()) {
			if (regionsTarget.get(indexTarget).mCopyNumberClusterType != clusterType) {
				regionsTarget.remove(indexTarget);
			} else {
				indexTarget++;
			}
		}
		
		// We only need to add elements (actually, their copies) if more still exist in the source array. 
		for (; indexSource < regionsSource.size(); indexSource++) {
			CopyNumberRegionRange regionSource = regionsSource.get(indexSource);
			if (regionSource.mCopyNumberClusterType == clusterType) {
				regionsTarget.add(regionSource.getCopy());
			}
		}
		
		return regionsTarget;
	}
	
	/** A helper function used by the SubsumesAlignedLeft and ConsumedByAlignedLeft cases. */
	private static void breakdownSubsumesConsumedAlignedLeftHelper(
			ArrayList<CopyNumberRegionRange> regionsList1, 
			ArrayList<CopyNumberRegionRange> regionsList2,
			int indexList1, int indexList2, 
			boolean makeCopyOfRegionFromList1) {
		
		CopyNumberRegionRange region1 = regionsList1.get(indexList1);
		CopyNumberRegionRange region2 = regionsList2.get(indexList2);
		
		CopyNumberRegionRange cnrrAfter = region1.getCopy();  // make a copy to keep the stats
		cnrrAfter.setRangeStart(region2.getRangeEnd() + 1);   // make right-adjacent to the toAdd region					
		regionsList1.add(indexList1 + 1, cnrrAfter);          // Add the next region
		
		CopyNumberRegionRange cnrrCurrent = (makeCopyOfRegionFromList1) ? 
				region1.getCopy() : 		// make a copy to keep the stats, and to modify the region without modifying the paramter version 
				region1;                    // don't make a copy, but allow the region itself to be modified later
						
		cnrrCurrent.setRangeEnd(region2.getRangeEnd());  // make right-aligned to the toAdd region
		regionsList1.set(indexList1, cnrrCurrent);
	}
	
	private static void breakdownThreeCasesHelper(
			ArrayList<CopyNumberRegionRange> regionsList1, 
			ArrayList<CopyNumberRegionRange> regionsList2,
			int indexList1, int indexList2, 
			boolean makeCopyOfRegionFromList1,
			boolean insertNewRegionIntoList1) {
		
		CopyNumberRegionRange region1 = regionsList1.get(indexList1);
		CopyNumberRegionRange region2 = regionsList2.get(indexList2);
		
		CopyNumberRegionRange cnrrBefore = region1.getCopy();  // make a copy to keep the stats
		cnrrBefore.setRangeEnd(region2.getRangeStart() - 1);  // make left-adjacent to the toAdd region
		if (insertNewRegionIntoList1) {
			regionsList1.add(indexList1, cnrrBefore);  // add the preceding block to the array at this point
			indexList1++;
		} else {
			regionsList2.add(indexList2, cnrrBefore);  // add the preceding block to the array at this point
			indexList2++;
		}
		
		CopyNumberRegionRange cnrrCurrent = (makeCopyOfRegionFromList1) ? 
				region1.getCopy() : 		// make a copy to keep the stats, and to modify the region without modifying the paramter version 
				region1;                    // don't make a copy, but allow the region itself to be modified later

		cnrrCurrent.setRangeStart(region2.getRangeStart());  // make left-aligned to the toAdd region
		regionsList1.set(indexList1, cnrrCurrent);
	}

	// ========================================================================
	private static void getBrowserTracksHelper_Helper_findMinMaxPositionsAcrossSamplesPerChromosome(int[] minPerChrom, int[] maxPerChrom, ArrayList<File> sampleFiles, ClusterType clusterType) {
		Arrays.fill(minPerChrom, Integer.MAX_VALUE);
		Arrays.fill(maxPerChrom, Integer.MIN_VALUE);
		
		for (File inFile : sampleFiles) {
			int lineCounter = -1;
			String line;
			BufferedReader in = IOUtils.getBufferedReader(inFile.getAbsolutePath());
			
			while ((line = IOUtils.getNextLineInBufferedReader(in)) != null) {
				if (++lineCounter == 0) { continue; }
				
				final Chrom chrom  = Chrom.getChrom  (Utils.extractNthColumnValue(line, ColCuratedTSV_Chrom,    Utils.TabStr));					
				final int position = Integer.parseInt(Utils.extractNthColumnValue(line, ColCuratedTSV_Position, Utils.TabStr));
				final ClusterType clusterTypeForSite = ClusterType.getClusterType(Utils.extractNthColumnValue(line, ColCuratedTSV_Cluster,      Utils.TabStr));

				if (clusterTypeForSite == clusterType) {
					// Set min, max indices across samples for this chromosome
					minPerChrom[chrom.ordinal()] = Math.min(minPerChrom[chrom.ordinal()], position);
					maxPerChrom[chrom.ordinal()] = Math.max(maxPerChrom[chrom.ordinal()], position);
				}				
			}
			IOUtils.closeBufferedReader(in);
		}		
	}
	
	// ========================================================================
	public static void genBrowserTracks(ArrayList<CopyNumberRegionsByChromosome> samplesWithRegions, 
										ArrayList<File> sampleFiles, 
										CopyNumberRegionsByChromosome recurrentRegionsForOneClusterType,
										ClusterType clusterType, String outDir) {		

		// First, we want to find the minimum and maximum positions of sites across samples for each chromosome
		int numChromArrayElements = Chrom.values().length;		
		int[] minStartPositionOnChromosomeAcrossSamples = Utils.newIntArray(numChromArrayElements, Integer.MAX_VALUE); 				
		int[] maxStartPositionOnChromosomeAcrossSamples = Utils.newIntArray(numChromArrayElements, Integer.MIN_VALUE);
		getBrowserTracksHelper_Helper_findMinMaxPositionsAcrossSamplesPerChromosome
			(minStartPositionOnChromosomeAcrossSamples, maxStartPositionOnChromosomeAcrossSamples, sampleFiles, clusterType);

		// Initialize color array and event count 2D array
		int[] rgb = Utils.newIntArray(3, 0);  				
		int[][] eventCount = Utils.newIntArray2D(numChromArrayElements, sampleFiles.size(), 0); 
		
		int sampleIndex = -1;
		int greyscaleBEDScore = 900;
		StringBuilder sb = new StringBuilder(4096);
		
		// Now go through the files
		for (File inFile : sampleFiles) {

			// Load all lines into memory, extract the header row, and then sort by chrom/position
			ArrayList<String> allLines = IOUtils.readAllLinesFromFile(inFile.getAbsolutePath(), false, true, sb);
			String headerString = sb.toString();
			Collections.sort(allLines, LineComparatorTab);			

			CopyNumberRegionsByChromosome regionsOneSample = samplesWithRegions.get(++sampleIndex);
			if (inFile.getName().indexOf(regionsOneSample.mSampleName) < 0) {
				Utils.throwErrorAndExit("ERROR: Samples don't match up: " + regionsOneSample.mSampleName + "\t" + inFile.getName());
			}
			String sampleNameRoot = regionsOneSample.mSampleName;
			
			System.out.println("Processing sample (" + sampleIndex + " / " + sampleFiles.size() + "):\t" + sampleNameRoot);

			// This is the list of regions
			int[] rowNum = Utils.newIntArray(numChromArrayElements, 0);			
			
			String snowString = "205,201,201";  // assign color to each row of region in cluster type
			ArrayList<ArrayList<String>> outList = new ArrayList<ArrayList<String>>();						 
			
			// We generate strings for the list of regions that we will print.  Just before that, however,
			// we will print some header information.
			for (Chrom chrom : Chrom.values()) {								
				ArrayList<String> outListForChrom = new ArrayList<String>(allLines.size());
				outList.add(outListForChrom);
				if (chrom.isInvalid()) continue;  // skip rest of code if invalid chrom (after of course, having added a dummy list)
					
				// Store header information
				constructHeaderString_BED(sb, chrom, minStartPositionOnChromosomeAcrossSamples[chrom.ordinal()], 
						                             maxStartPositionOnChromosomeAcrossSamples[chrom.ordinal()]);
				outListForChrom.add(sb.toString());  // replace the dummy placeholder string
				outListForChrom.add("browser hide all");								
				outListForChrom.add( constructTrackNameString_BED(sampleNameRoot, sb, true).toString() );
	
				for (CopyNumberRegionRange cnrr : regionsOneSample.mRegionsByChrom.get(chrom.ordinal())) {
					constructRowString_BEDFormat(chrom, sb, cnrr.getRangeStart(), cnrr.getRangeEnd() + 1, ++rowNum[chrom.ordinal()], snowString);
					outListForChrom.add(sb.toString());
				}
			}
						
			// Go through all the lines
			for (int lineIndex = 0; lineIndex < allLines.size(); lineIndex++) {
				String line = allLines.get(lineIndex);

				final Chrom chrom  = Chrom.getChrom                       (Utils.extractNthColumnValue(line, ColCuratedTSV_Chrom,    Utils.TabStr));					
				final int position = Integer.parseInt                     (Utils.extractNthColumnValue(line, ColCuratedTSV_Position, Utils.TabStr));
				final ClusterType clusterTypeForSite = ClusterType.getClusterType(Utils.extractNthColumnValue(line, ColCuratedTSV_Cluster,      Utils.TabStr));
				final SNVType mutationType = SNVType.getSNVType                  (Utils.extractNthColumnValue(line, ColCuratedTSV_MutationType, Utils.TabStr));
				final VariantLocation varLoc = VariantLocation.getVariantLocation(Utils.extractNthColumnValue(line, ColCuratedTSV_VariantLocation, Utils.TabStr));

				int score = 0;
				ColorPastel theColorMutationType = ColorPastel.Pastel_Cyan;  // set as default
				
				if (clusterTypeForSite == clusterType) {					
					// Increment the event count for this chromosome and sample
					eventCount[chrom.ordinal()][sampleIndex]++;					
					score = greyscaleBEDScore; //used with grayscale BEDs (historical artifact)
										
					if (mutationType != null) {
						switch(mutationType) {
						case NonSynonymous_SNV: theColorMutationType = ColorPastel.RGB_Red;    break;
						case Synonymous_SNV:    theColorMutationType = ColorPastel.Dark_Green; break;					
						}
					}
					
					if (varLoc == VariantLocation.Somatic) {
						theColorMutationType = ColorPastel.Blue;
					}
					
				} else if  (clusterTypeForSite == ClusterType.HETGermline) {					
					score = 300; //used with grayscale BEDs (historical artifact)
					
					if (mutationType != null) {
						switch(mutationType) {
						case NonSynonymous_SNV: theColorMutationType = ColorPastel.Yellow_Orange; break;
						case Synonymous_SNV:    theColorMutationType = ColorPastel.Light_Green;   break;						
						}
					}
				}

				if (score == greyscaleBEDScore) {
					theColorMutationType.getRGBValues(rgb);					
					sb.setLength(0);
					sb.append(rgb[0]).append(Utils.CommaStr).append(rgb[1]).append(Utils.CommaStr).append(rgb[2]);
					constructRowString_BEDFormat(chrom, sb, position, position + 1, ++rowNum[chrom.ordinal()], sb.toString());
					outList.get(chrom.ordinal()).add(sb.toString());
				}								
			}
			
			// Now write output to disk
			for (Chrom chrom : Chrom.values()) {
				if (chrom.isInvalid()) continue;
				
				String outFilename = 
						outDir + File.separator + constructGenBrowserTrackFilenamePerSampleAndCluster(sampleNameRoot, clusterType, chrom);
				IOUtils.writeOutputFile(outFilename, outList.get(chrom.ordinal()));	
			}
		}
		
		// Now we create the output files for the recurrence track
		String recurrenceNameRoot = clusterType.name() + "_Recurrence";
		float bedScoreMin   = 200f;
		float bedScoreRange = 900f - bedScoreMin;		
		
		for (Chrom chrom : Chrom.values()) {				
			if (chrom.isInvalid()) continue;
			
			ArrayList<String> outListForChrom = new ArrayList<String>(10000);				
			
			// Store header information
			constructHeaderString_BED(sb, chrom, minStartPositionOnChromosomeAcrossSamples[chrom.ordinal()], 
								                 maxStartPositionOnChromosomeAcrossSamples[chrom.ordinal()]);
			outListForChrom.add(sb.toString());  // replace the dummy placeholder string
			outListForChrom.add("browser hide all");										
			outListForChrom.add( constructTrackNameString_BED(recurrenceNameRoot, sb, false).toString() );
			
			// Now write the regions	
			int rowNum = 0;
			float recurrenceMin = Float.MAX_VALUE;
			float recurrenceMax = Float.MIN_VALUE;
			
			// First get the min and max values to scale
			for (CopyNumberRegionRange cnrr : recurrentRegionsForOneClusterType.mRegionsByChrom.get(chrom.ordinal())) {
				recurrenceMin = Math.min(recurrenceMin, cnrr.mRecurrenceScore);
				recurrenceMax = Math.max(recurrenceMax, cnrr.mRecurrenceScore);
			}
			float recurrenceMinMaxRange = recurrenceMax - recurrenceMin + 1;
			
			// Scale accordingly -- this idea was inspired from Siddharth Reddy
			for (CopyNumberRegionRange cnrr : recurrentRegionsForOneClusterType.mRegionsByChrom.get(chrom.ordinal())) {
				float score = cnrr.mRecurrenceScore;
				float fraction = (score - recurrenceMin) / recurrenceMinMaxRange;
				float newScoreRaw = fraction * bedScoreRange;
				int newScoreAdj = (int) (newScoreRaw + bedScoreMin);
				
				constructRowString_BEDFormat(chrom, sb, cnrr.getRangeStart(), cnrr.getRangeEnd() + 1, ++rowNum, newScoreAdj, "");
				outListForChrom.add(sb.toString());
			}
			
			// Now write output to disk
			String outFilename = 
					outDir + File.separator + constructGenBrowserTrackFilenamePerSampleAndCluster(recurrenceNameRoot, clusterType, chrom);
			IOUtils.writeOutputFile(outFilename, outListForChrom);	
		}
		
		// Now, we want to write the manifest file
		writeManifestFile(outDir, clusterType, eventCount, samplesWithRegions, recurrenceNameRoot);		
	}
	
	// ========================================================================
	public static StringBuilder constructHeaderString_BED(StringBuilder sb, Chrom chrom, int rangeStart, int rangeEnd) {
		sb.setLength(0);
		sb.append("browser position chr").append(chrom.ordinal()).append(Utils.ColonString)
		  .append(rangeStart).append("-").append(rangeEnd);
		return sb;
	}
	
	// ========================================================================
	public static StringBuilder constructTrackNameString_BED(String sampleNameRoot, StringBuilder sb, boolean useRGB) {
		sb.setLength(0);				
		sb.append("track name=").append(Utils.DoubleQuoteStr).append(sampleNameRoot).append(Utils.DoubleQuoteStr).append(Utils.SpaceString)
		  .append("description=").append(Utils.DoubleQuoteStr).append(Utils.SpaceString).append(Utils.DoubleQuoteStr).append(Utils.SpaceString)
		  .append("visibility=dense").append(Utils.SpaceString);
		
		if (useRGB) {
			sb.append("itemRgb=").append(Utils.DoubleQuoteStr).append("on").append(Utils.DoubleQuoteStr);
		} else {
			sb.append("useScore=1");
		}
		return sb;
	}

	// ========================================================================
	public static StringBuilder constructRowString_BEDFormat(Chrom chrom, StringBuilder sb, int rangeStart, int rangeEnd, int rowNum, String colorString) {
		return constructRowString_BEDFormat(chrom, sb, rangeStart, rangeEnd, rowNum, 0, colorString);
	}
	
	public static StringBuilder constructRowString_BEDFormat(Chrom chrom, StringBuilder sb, int rangeStart, int rangeEnd, int rowNum, int score, String colorString) {
		sb.setLength(0);
		sb.append(ChromPrefix).append(chrom.ordinal()).
		   append(Utils.SpaceString).append(rangeStart).
		   append(Utils.SpaceString).append(rangeEnd).
		   append(Utils.SpaceString).append("row").append(rowNum).
		   append(Utils.SpaceString).append(score).
		   append(Utils.SpaceString).append("+").				   
		   append(Utils.SpaceString).append(rangeStart).
		   append(Utils.SpaceString).append(rangeEnd);
		
		if (!colorString.isEmpty()) {
			sb.append(Utils.SpaceString).append(colorString);
		}
		
		return sb;
	}
	
	
	// ========================================================================
	// Convenience function to write the manifest file
	private static void writeManifestFile(String outDir, ClusterType clusterType, int[][] eventCount, 
										  ArrayList<CopyNumberRegionsByChromosome> samplesWithRegions,
										  String recurrenceNameRoot) {
		
		String filePrefix = outDir + File.separator;
		//String filePrefix = "h"
		
		for (Chrom chrom : Chrom.values()) {
			if (chrom.isInvalid()) continue;

			String outFilename = 
					outDir + File.separator + clusterType.name() + Utils.DotStr + ChromPrefix + chrom.ordinal() + Utils.DotStr + "Manifest.txt";
			
			// Now we need to add to the manifest file (specific to a chrom and cluster type)
			// We need to sort the samples based on the event counts.  We use a trick using
			// bit-shifting to sort the coutn and the sample names together
			int[] eventCountForChrom = eventCount[chrom.ordinal()];
			long[] eventCountForChromWithSampleID = new long[eventCountForChrom.length];

			for (int indexSample = 0; indexSample < eventCountForChrom.length; indexSample++) {
				long combinedValue = 0L | (((long) eventCountForChrom[indexSample]) << Integer.SIZE) | ((long) indexSample);
				eventCountForChromWithSampleID[indexSample] = combinedValue;
			}

			// Now sort on the combined value.  Note that the sorting will happen effectively 
			// by the event count, as it resides in the MSBs.
			Arrays.sort(eventCountForChromWithSampleID);			

			// Create the buffered writer and insert the filenames of the recurrence tracks
			BufferedWriter outManifest = IOUtils.getBufferedWriter(outFilename);
			String recurrenceFilename = 
					filePrefix + constructGenBrowserTrackFilenamePerSampleAndCluster(recurrenceNameRoot, clusterType, chrom);
			IOUtils.writeToBufferedWriter(outManifest, recurrenceFilename, true);
						
			// Now go through and figure out the sampleIndices in sorted order, from highest to lowest			
			for (int i = eventCountForChromWithSampleID.length - 1; i >= 0; i--) {
				int indexSample = (int) (eventCountForChromWithSampleID[i] & 0xFFFFFFFF);
				String theSampleName = samplesWithRegions.get(indexSample).mSampleName;
				String samplePath = filePrefix + constructGenBrowserTrackFilenamePerSampleAndCluster(theSampleName, clusterType, chrom);
				IOUtils.writeToBufferedWriter(outManifest, samplePath, true);
			}
			IOUtils.closeBufferedWriter(outManifest);
		}
	}

	// ========================================================================
	// Convenience function for constructing an output filename string
	private static String constructGenBrowserTrackFilenamePerSampleAndCluster(String sampleName, ClusterType clusterType, Chrom chrom) {
		return clusterType.name() + Utils.DotStr + ChromPrefix + chrom.ordinal() + Utils.DotStr 
				                  + sampleName + GenBrowserTrack + Utils.FileExtensionTSV.mExtension;
	}
	
	// ========================================================================

	
	/**
	 * Generate 'master' gene enrichment table (used to generate histograms).
	 * @param inDir curated SNP calls
	 */
	public static void getGeneEnrichment(String inDir, String outFilename) {
		File[] files = (new File(inDir)).listFiles();
		FileExtensionAndDelimiter fileExtDelim = Utils.FileExtensionTSV; 
				
		StringBuilder sb = new StringBuilder(4096);
		
		// A list of genes that we will keep binary sorted for efficieny purposes
		ArrayList<Gene> genes = new ArrayList<Gene>(); 
		Gene dummyGene = new Gene("", Chrom.c0);  // We'll use this gene for binary searching
		
		for (File file : files) { //iterate through curated SNP calls
			if (file.getName().indexOf(fileExtDelim.mExtension) >= 0) {
				System.out.println(file.getName());
				
				// Load all lines into memory, extract the header row, and then sort by chrom/position
				ArrayList<String> allLines = IOUtils.readAllLinesFromFile(file.getAbsolutePath(), false, true, sb);
				String headerString = sb.toString();
				Collections.sort(allLines, LineComparatorTab);					
				
				int startingRow = 0;  // because the header was stripped away
				for (int row = startingRow; row < allLines.size(); row++) {
				
					//if (row % 10000 == 0) System.out.println(row + " of " + allLines.size());
					
					String components[] = allLines.get(row).split(Utils.TabStr); 
					String geneName =               components[ColCuratedTSV_Gene];
					Chrom chrom  = Chrom.getChrom(  components[ColCuratedTSV_Chrom]);
					int position = Integer.parseInt(components[ColCuratedTSV_Position]);

					
					if (!geneName.equals(MissingGeneNameValue)) { //avoid ".", which crops up a lot
						if (geneName.indexOf("dist") >= 0) { //gene names can sometimes come with an uninteresting/irrelevant prefix
							geneName = geneName.split("\\(")[0];
						}
					
						// Set the dummy gene for binary search
						dummyGene.mLabel = geneName;
						dummyGene.mChrom = chrom;
						
						int resultIndex = Collections.binarySearch(genes, dummyGene);
						if (resultIndex < 0) {  
							// we haven't seen this gene before
							resultIndex = -(resultIndex + 1);  // calculate the proper insertion point 							 
							genes.add(resultIndex, new Gene(geneName, chrom));							
						}
						Gene currentGene = genes.get(resultIndex);
						
						if (position > currentGene.mMaxBasePairPosition) { //get right-bound of gene's range of variants
							currentGene.mMaxBasePairPosition = position; 
						}
						if (position < currentGene.mMinBasePairPosition) { //...left-bound...
							currentGene.mMinBasePairPosition = position;
						}
												
						SNVType mutationType = SNVType.getSNVType(components[ColCuratedTSV_MutationType]);
						if (mutationType != null) {
							currentGene.incrementCountForMutationType(mutationType);  // increment synonymous or nonsynonymous variant count
						}
						
						VariantLocation variantLocation = VariantLocation.getVariantLocation(components[ColCuratedTSV_VariantLocation]);
						if (variantLocation == null) { Utils.throwErrorAndExit("ERROR: Invalid variantLocation: " + components[ColCuratedTSV_VariantLocation]); }
						currentGene.incrementCountForVariantLocation(variantLocation);
																
						ClusterType clusterType = ClusterType.getClusterType(components[ColCuratedTSV_Cluster]);
						if (clusterType == null) { Utils.throwErrorAndExit("ERROR: Invalid cluster type: " + components[ColCuratedTSV_Cluster]); }
						currentGene.incrementCountForClusterType(clusterType);  // increment LOH/DUP/&c. count
						double vafTumor = Double.parseDouble(components[ColCuratedTSV_VafTumor]);
						if ((clusterType == ClusterType.LOH) && (vafTumor > 0.5)) {
							currentGene.mCountLOHreferenceLost++;
						}										
						currentGene.addPatientIfNotAlreadyAdded(file.getName(), clusterType);
					}
				}
			}
		}
		
		String logStr = "_log";
		String recurrenceStr = "_recurrence";
		String densityStr = "_density";
		String[] columnHeaders = new String[] { 
				"chr", "bp_start", "bp_end", "length", "gene", 
				SNVType.NonSynonymous_SNV.toLowerCase(), 
				SNVType.Synonymous_SNV.toLowerCase(), 
				VariantLocation.Germline.toLowerCase(), 
				VariantLocation.Somatic.toLowerCase(),
				
				ClusterType.Dup.name(),
				ClusterType.LOH.name(), 
				ClusterType.LOH.name() + "_refLost", 
				ClusterType.HETGermline.name(),
				ClusterType.HETSomatic.name(),
				
				ClusterType.Dup.name() + logStr,
				ClusterType.LOH.name() + logStr, 
				ClusterType.LOH.name() + "_refLost" + logStr, 
				ClusterType.HETGermline.name() + logStr,
				ClusterType.HETSomatic.name() + logStr,
				
				ClusterType.Dup.name()              + densityStr,
				ClusterType.LOH.name()              + densityStr, 
				ClusterType.LOH.name() + "_refLost" + densityStr, 
				ClusterType.HETGermline.name()      + densityStr,
				ClusterType.HETSomatic.name()       + densityStr,
				
				ClusterType.Dup.name()         + recurrenceStr, 
				ClusterType.LOH.name()         + recurrenceStr, 
				ClusterType.HETGermline.name() + recurrenceStr,
				ClusterType.HETSomatic.name()  + recurrenceStr
		};
		String headerStr = Utils.constructColumnDelimitedString(columnHeaders, fileExtDelim.mDelimiter, sb, true).toString();

		BufferedWriter out = IOUtils.getBufferedWriter(outFilename);
		IOUtils.writeToBufferedWriter(out, headerStr, true);
		for (Gene gene : genes) {			
			sb.setLength(0);
			sb.append(gene.mChrom.getCode())
				.append(fileExtDelim.mDelimiter).append(gene.mMinBasePairPosition)
				.append(fileExtDelim.mDelimiter).append(gene.mMaxBasePairPosition)
				.append(fileExtDelim.mDelimiter).append(gene.getRangeLength())
				.append(fileExtDelim.mDelimiter).append(gene.mLabel)
				.append(fileExtDelim.mDelimiter).append(gene.countsToString(fileExtDelim.mDelimiter));
			IOUtils.writeToBufferedWriter(out, sb.toString(), true);
			IOUtils.flushBufferedWriter(out);
		}
		IOUtils.closeBufferedWriter(out);
	}
	
	
	// ========================================================================
	// ENTRY POINT
	// ========================================================================
	public static void main(String[] args) {
		long sys_time_init = System.currentTimeMillis();
		
		String nafTafInptus     = "naf_taf_inputs";
		String classifiedSites  = "sites_classified";
		String vafPlots         = "vafPlots";
		String regions          = "regions";
		String browserTracks    = "browser_tracks";
		
		
		String root = args[0]; //project directory
		switch (Integer.parseInt(args[1])) { //args[1] --> 'switchboard' parameter
			case 0:
				SeqPlatform platform = SeqPlatform.getPlatform(Integer.parseInt(args[2]));
				curateSNPCalls(root + File.separator + nafTafInptus, 
						       root + File.separator + classifiedSites, 
						       root + File.separator + vafPlots, 
						       platform); //args[2] --> 0::Illumina, 1::SOLiD
				break;
			case 1:
				segmentRegionsAllFiles(root + File.separator + classifiedSites, 
						               root + File.separator + regions,
						               root + File.separator + browserTracks);
				break;
			case 2:
				getGeneEnrichment(root + File.separator + classifiedSites, 
								  root + File.separator + "gene_enrichment.csv");
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
