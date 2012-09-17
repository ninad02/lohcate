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
				ClusterType[] clustersGermline = getClusters(allVariantRows, indexFirstSomaticRowInAllVariants, platform); //get cluster assignments (HET ball, LOH sidelobes, DUP wedge, &c.)
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
						gene = germCols[Col_NAFTAFInput_HugoSymbol];
						mutationType = germCols[Col_NAFTAFInput_MutationType].split("_SNV")[0];
					}
					
					sb.append(fileExtDelim.mDelimiter).append(gene)
					  .append(fileExtDelim.mDelimiter).append(mutationType)
					  .append(fileExtDelim.mDelimiter).append(targetTissue)
					  .append(fileExtDelim.mDelimiter).append(clusterType);
					
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
				String[] components = rows.get(row).split("\t");
				float vafNormal = extractVAFNormal(components, platform);
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
	 */
	public static ClusterType[] getClusters(ArrayList<String> rows, int startingRowGermlineOrSomaticOrAll, SeqPlatform platform) {
		
		AlleleFrequencyStatsForSample afStatsSample = new AlleleFrequencyStatsForSample();
		afStatsSample.tabulateAndPerformStatistics(rows, platform);
		
		// Now we adjust the frames based on the standard deviation
		float vafNormalFrameAdjustedLower = afStatsSample.mCountMean - (NAF_STRIP_EXPANDER * afStatsSample.mStdDev);
		float vafNormalFrameAdjustedUpper = afStatsSample.mCountMean + (NAF_STRIP_EXPANDER * afStatsSample.mStdDev);
		//System.out.println("MEAN = " + count_mean + " | STD_DEV = " + std_dev);
		//System.out.println("NAF FRAME = " + (count_mean - std_dev) + " - " + (count_mean + std_dev));
		
		//apply DBScan to points within NAF frame
		ArrayList<Floint> points = new ArrayList<Floint>(rows.size());
		for (int row = 0; row < rows.size(); row++) {
			String line = rows.get(row);			
			float vafNormal = extractVAFNormal(line, platform);
			if (Utils.inRangeLowerExclusive(vafNormal, vafNormalFrameAdjustedLower, vafNormalFrameAdjustedUpper)) {	
				points.add(  new Floint(extractVAFTumor(line, platform), vafNormal)  );			
			}
		}
	
		//Vv.vV well/poorly tuned parameters for different data sets
		//data set --> (NAF_STRIP_EXPANDER, (HET_BALL_EPS, HET_BALL_MINPTS), (DUP_WEDGE_LASSO, DUP_WEDGE_MINPTS))
		//target-aml --> (1, (0.035, 100), (0.015, 100))
		//target-all --> (1, (0.035, 100), (0.01, 100))
		//pvera --> (1.25, (0.035, 100), (0.015, 100))
		//hepato --> (1.25, (0.035, 100), (0.015, 100))
		//renal-we --> (1, (0.05, 500), (0.02, 350))
		
		System.out.println("Begin clustering algorithm: " + (new Date()).toString());
		
		DBScanFaster dbscanner = new DBScanFaster(points, /* HET_BALL_EPS, HET_BALL_MINPTS*/ HET_BALL_EPS + DUP_WEDGE_LASSO, DUP_WEDGE_MINPTS, 0, 0, 1, 1); //parameters for capturing HET ball//KDBSCAN(param, 10, 0.0325f);
		//DBScanFast dbscanner = new DBScanFast(points, HET_BALL_EPS, HET_BALL_MINPTS); //parameters for capturing HET ball//KDBSCAN(param, 10, 0.0325f);
		dbscanner.cluster();
		int clusterIDofHetBall = dbscanner.getCentralClusterID();
		Boolean[] pointsWithinRadius = dbscanner.getPointsWithinMinRadiusOfCluster(clusterIDofHetBall);
		int[] clusterAssignments = dbscanner.getClustAssignments();  // save and cache	
		
		
		// Now, re-run DBScan, but with changed parameters		
		//dbscanner.changeParams(HET_BALL_EPS + DUP_WEDGE_LASSO, DUP_WEDGE_MINPTS);		
		//dbscanner.cluster();
		//int clusterIDofHetBallWithWedge = dbscanner.getCentralClusterID();
		//int[] clusterAssignmentsWithWedge = dbscanner.getClustAssignments();

		// We now do a third pass, this time allowing only those points that are not part of the het cluster
		// We then perform a lower neighbor threshold on these points so that they can be declared as LOH
		ArrayList<Floint> nonHetPoints = new ArrayList<Floint>(points.size());
		boolean[] isNonHetPoint = new boolean[points.size()];		
		
		for (int i = 0; i < clusterAssignments.length; i++) {
			isNonHetPoint[i] = (clusterAssignments[i] != clusterIDofHetBall); 
			if (isNonHetPoint[i]) {
				nonHetPoints.add(points.get(i));
			}
		}
		DBScanFaster dbscannerNonHet = new DBScanFaster(nonHetPoints, NON_HET_BALL_EPS, NON_HET_BALL_MINPTS, 0, 0, 1, 1);
		dbscannerNonHet.cluster();
		int[] clusterAssignmentsNonHet = dbscannerNonHet.getClustAssignments();
				
		System.out.println("End clustering algorithm: " + (new Date()).toString());

		int clusterIDofDup = -1;
		int clusterIDofNull = -2;
		int nonHetIndex = -1;
		for (int i = 0; i < clusterAssignments.length; i++) {
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
						clusterAssignments[i] = clusterIDofNull;
					} else {
						clusterAssignments[i] = clusterIDofDup;   // only change if we're in a het ball region
					}
				}
			}
			
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
		int indexInClusterAssignments = -1;  // we need to keep a special index, since not all rows are used.
		ClusterType[] returnClusters = new ClusterType[rows.size()];
		
		for (int row = 0; row < rows.size(); row++) {
			
			//String[] components = rows.get(row).split("\t");
			String line = rows.get(row);
			float vafNormal = extractVAFNormal(line, platform);
			boolean vafInRangeNormal = Utils.inRangeLowerExclusive(vafNormal, vafNormalFrameAdjustedLower, vafNormalFrameAdjustedUpper);
			
			if ((row >= startingRowGermlineOrSomaticOrAll) && (!vafInRangeNormal)) {
				// The vafNormal is either very low (homozygous reference) or very high (homozygous common variant).
				// We do some very simple decision making now (which should be replaced by formal clustering later)
				// to partition the calls.
				float justBelowZero = -0.0001f;				
				float hetBoundaryLower = 0.3333f;
				float hetBoundaryUpper = 0.6667f;
				float vafTumor = extractVAFTumor(line, platform);
				if (Utils.inRangeLowerExclusive(vafNormal, justBelowZero, vafNormalFrameAdjustedLower)) {
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
				} else if (Utils.inRangeLowerExclusive(vafNormal, vafNormalFrameAdjustedUpper, MaxVariantAlleleFrequency)) {
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
					int assignedClusterID = clusterAssignments[indexInClusterAssignments]; // we create a local variable for fast test modifications 
					
					if (assignedClusterID == clusterIDofHetBall) {
						returnClusters[row] = ClusterType.HETGermline; //HET
						
					} else if (assignedClusterID == clusterIDofDup) {
						returnClusters[row] = ClusterType.Dup; //DUP
						
					} else if (assignedClusterID == clusterIDofNull) {
						returnClusters[row] = ClusterType.Null;  // we're on a diagonal
						
					} else if (assignedClusterID == dbscanner.getClusterIDOfNoise()) { 
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
	
	// Tests whether the point lies roughly on the x = y diagonal, given a leeway from a ratio of 1.0
	private static boolean pointOnDiagonal(Floint point, double leeway) {
		double ratio = point.mY / point.mX;
		double equalRatio = 1.0;
		leeway = Math.min(leeway, equalRatio);
		return ( ((equalRatio - leeway) <= ratio) && (ratio <= (equalRatio + leeway)) );
	}
	
	// Extracts and calculates the variant allele frequency depending on the platform 
	private static float extractVAFNormal(String[] components, SeqPlatform platform) {
		switch(platform) {
		case Illumina: return  Float.parseFloat(components[Col_NAFTAFInput_VariantRatioNormal]);
		case SOLiD:    return (Float.parseFloat(components[6]) / Float.parseFloat(components[5]));
		}
		return Float.NaN;
	}
	
	// Extracts and calculates the variant allele frequency in the normal depending on the platform tab-delimited file
	private static float extractVAFNormal(String line, SeqPlatform platform) {
		switch(platform) {
		case Illumina: return  Float.parseFloat(Utils.extractNthColumnValue(line, Col_NAFTAFInput_VariantRatioNormal, Utils.TabStr));
		case SOLiD:    return (Float.parseFloat(Utils.extractNthColumnValue(line, 6, Utils.TabStr)) / 
				               Float.parseFloat(Utils.extractNthColumnValue(line, 5, Utils.TabStr)));
		}
		return Float.NaN;		
	}
	
	private static float extractVAFTumor(String[] components, SeqPlatform platform) {
		switch(platform) {
		case Illumina: return  Float.parseFloat(components[Col_NAFTAFInput_VariantRatioTumor]);
		case SOLiD:    return (Float.parseFloat(components[4]) / Float.parseFloat(components[3]));
		}
		return Float.NaN;
	}
	
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
		for (int clusterIndex = 0; clusterIndex < ClusterType.DupLOHHetG.length; clusterIndex++) {
			ClusterType clusterType = ClusterType.DupLOHHetG[clusterIndex];
			ArrayList<CopyNumberRegionsByChromosome> regionsInSamplesForOneClusterType = regionsInSamplesPerClusterType.get(clusterIndex);
			
			genBrowserTracks(regionsInSamplesForOneClusterType, validFilesList, clusterType, outDirBrowserTracks);
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
	public static void genBrowserTracks(ArrayList<CopyNumberRegionsByChromosome> regionsPerSample, 
										ArrayList<File> sampleFiles, 
										ClusterType clusterType, String outDir) {		

		// First, we want to find the minimum and maximum positions of sites across samples for each chromosome
		int numChromArrayElements = Chrom.values().length;
		StringBuilder sb = new StringBuilder(4096);
		int[] minStartPositionOnChromosomeAcrossSamples = new int[numChromArrayElements];
		int[] maxStartPositionOnChromosomeAcrossSamples = new int[numChromArrayElements];
		getBrowserTracksHelper_Helper_findMinMaxPositionsAcrossSamplesPerChromosome
			(minStartPositionOnChromosomeAcrossSamples, maxStartPositionOnChromosomeAcrossSamples, sampleFiles, clusterType);

		// Initialize color array
		int[] rgb = new int[3];
				
		int[][] eventCount = new int[numChromArrayElements][sampleFiles.size()];
		for (int[] eventsPerChrom : eventCount) {
			Arrays.fill(eventsPerChrom, 0);
		}
		
		int sampleIndex = -1;
		int greyscaleBEDScore = 900;
		
		// Now go through the files
		for (File inFile : sampleFiles) {

			// Load all lines into memory, extract the header row, and then sort by chrom/position
			ArrayList<String> allLines = IOUtils.readAllLinesFromFile(inFile.getAbsolutePath(), false, true, sb);
			String headerString = sb.toString();
			Collections.sort(allLines, LineComparatorTab);			

			CopyNumberRegionsByChromosome regionsOneSample = regionsPerSample.get(++sampleIndex);
			if (inFile.getName().indexOf(regionsOneSample.mSampleName) < 0) {
				Utils.throwErrorAndExit("ERROR: Samples don't match up: " + regionsOneSample.mSampleName + "\t" + inFile.getName());
			}
			String sampleNameRoot = regionsOneSample.mSampleName;

			// This is the list of regions
			int rowNum = 0;			
			ArrayList<ArrayList<String>> outList = new ArrayList<ArrayList<String>>();
			outList.add(new ArrayList<String>());  // add dummy list for chr 0
			String snowString = "205,201,201";  // assign color to each row of region in cluster type 
			
			// We generate strings for the list of regions that we will print.  Just before that, however,
			// we will print some header information.
			for (Chrom chrom : Chrom.Autosomes) {
				ArrayList<String> outListForChrom = new ArrayList<String>();
				outList.add(outListForChrom);
					
				// Store header information
				sb.setLength(0);
				sb.append("browser position chr").append(chrom.ordinal()).append(Utils.ColonString)
				  .append(minStartPositionOnChromosomeAcrossSamples[chrom.ordinal()]).append("-")
				  .append(maxStartPositionOnChromosomeAcrossSamples[chrom.ordinal()]);
				outListForChrom.add(sb.toString());				
				outListForChrom.add("browser hide all");
								
				sb.setLength(0);				
				sb.append("track name=").append(Utils.DoubleQuoteStr).append(sampleNameRoot).append(Utils.DoubleQuoteStr).append(Utils.SpaceString)
				  .append("description=").append(Utils.DoubleQuoteStr).append(Utils.SpaceString).append(Utils.DoubleQuoteStr).append(Utils.SpaceString)
				  .append("visibility=dense").append(Utils.SpaceString).append("itemRgb=").append(Utils.DoubleQuoteStr).append("on").append(Utils.DoubleQuoteStr);
				outListForChrom.add(sb.toString());
	
				for (CopyNumberRegionRange cnrr : regionsOneSample.mRegionsByChrom.get(chrom.ordinal())) {
					sb.setLength(0);
					sb.append(ChromPrefix).append(chrom.ordinal()).
					   append(Utils.SpaceString).append(cnrr.getRangeStart()).
					   append(Utils.SpaceString).append(cnrr.getRangeEnd()).
					   append(Utils.SpaceString).append("row").append(++rowNum).
					   append(Utils.SpaceString).append("0 +").
					   append(Utils.SpaceString).append(cnrr.getRangeStart()).
					   append(Utils.SpaceString).append(cnrr.getRangeEnd()).
					   append(Utils.SpaceString).append(snowString);
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
				ColorPastel theColor = null;
				if (clusterTypeForSite == clusterType) {
					// Increment the event count for this chromosome and sample
					eventCount[chrom.ordinal()][sampleIndex]++;					
					score = greyscaleBEDScore; //used with grayscale BEDs (historical artifact)
					
					switch(mutationType) {
					case NonSynonymous: theColor = ColorPastel.RGB_Red;    break;
					case Synonymous:    theColor = ColorPastel.Dark_Green; break;
					default:            theColor = ColorPastel.White;      break; 
					}
					
					if (varLoc == VariantLocation.Somatic) {
						theColor = ColorPastel.Blue;
					}
					
				} else if  (clusterTypeForSite == ClusterType.HETGermline) {					
					score = 300; //used with grayscale BEDs (historical artifact)
					
					switch(mutationType) {
					case NonSynonymous: theColor = ColorPastel.Yellow_Orange; break;
					case Synonymous:    theColor = ColorPastel.Light_Green;   break;
					default:            theColor = ColorPastel.White;         break; 
					}
				}

				if (score > 0) {
					theColor.getRGBValues(rgb);
					sb.setLength(0);
					sb.append(ChromPrefix).append(chrom.ordinal()).
					   append(Utils.SpaceString).append(position).
					   append(Utils.SpaceString).append(position).
					   append(Utils.SpaceString).append("row").append(++rowNum).
					   append(Utils.SpaceString).append("0 +").
					   append(Utils.SpaceString).append(position).
					   append(Utils.SpaceString).append(position).
					   append(Utils.SpaceString).append(rgb[0]).append(Utils.CommaStr).append(rgb[1]).append(Utils.CommaStr).append(rgb[2]);
					outList.get(chrom.ordinal()).add(sb.toString());
				}								
			}
			
			// Now write output to disk
			for (Chrom chrom : Chrom.Autosomes) {
				String outFilename = outDir + File.separator + clusterType.name() + Utils.DotStr + ChromPrefix + Utils.DotStr + chrom.ordinal() +  
						Utils.DotStr + sampleNameRoot + GenBrowserTrack + Utils.FileExtensionTSV.mExtension;
				IOUtils.writeOutputFile(outFilename, outList.get(chrom.ordinal()));
			}
		}
	}

	
	/*
	

		}
		
		
		System.out.println("writing to file...");
		String manifest;
		ArrayList<Integer> patsort, blacklist;
		int temp_max, temp_max_ind;
		for (int i = 1; i<=22; i++) { //chromosomes
			chr = Integer.toString(i);
			for (int k = 0; k<cluster_names.length - 1; k++) { //clusters
				manifest = "";
				patsort = new ArrayList<Integer>(); 
				blacklist = new ArrayList<Integer>();
				
				/* The following block jsut sorts patients by event count very inefficiently. *
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
				
				/* NPD: Won't need to combine germline and somatic in new code since it's already combined. *
				for (int p : patsort) { //patients
					toWrite[i-1][k][p][0] += toWrite[i-1][k][p][1]; //combine germline and somatic files
					toWrite[i-1][k][p][0] = 
							"browser hide all\n"
								+ "track name=\"" + pats.get(p).replace(".germline.csv", "") 
								+ "\" description=\" \" visibility=dense itemRgb=\"on\"\n" 
								+ toWrite[i-1][k][p][0];
					
					toWrite[i-1][k][p][0] = 
							"browser position chr" 
							+ chr + ":" 
							+ big_min[i-1][k] + "-" + big_max[i-1][k] + "\n" 
							+ toWrite[i-1][k][p][0];
					
					FileOps.writeToFile(outDir + "/chr" + chr + "_" + cluster_names[k] + "/" + pats.get(p).replace(".germline", ""), toWrite[i-1][k][p][0]);
					manifest += "https://raw.github.com/sidreddy96/lohcate/master/" + inDir.split("/")[inDir.split("/").length - 2] + "/browser_tracks/chr" + chr + "_" + cluster_names[k] + "/" + pats.get(p).replace(".germline", "") + "\n"; //you'll have to deal with this line as you please
				}
				FileOps.writeToFile(outDir + "/chr" + chr + "_" + cluster_names[k] + ".txt", manifest);
			}
		}
	}
	*/
	
	
	/**
	 * Define 'contiguous' regions of LOH, given our curated SNP calls.
	 * @param inDir curated SNP calls
	 */
	/*
	public static void segmentRegionsOld(String inDir, String outDir) {
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
						chr = Integer.toString(i);
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
						for (int k = 0; k<cluster_names.length - 1; k++) { //iterate through {dup, loh, roc-loh}
							String[] reg_load = FileOps.loadFromFile(reg_inDir + "/chr" + columns[0].replace("chr", "") + "/" + cluster_names[k] + ".csv").split("\n");
							
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
				for (int i = 0; i<cluster_names.length - 1; i++) { //iterate through {dup, loh, roc-loh}
					regions = new ArrayList<Point>(); //Point.x --> region start, Point.y --> region end
					recur = new ArrayList<Double>();
					region_pats = new ArrayList<ArrayList<String>>();
					try {
						load = FileOps.loadFromFile(dir.getAbsolutePath() + "/" + cluster_names[i] + ".csv").split("\n");
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
					chr = Integer.toString(i);
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
										} else { //if synonymous
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
			chr = Integer.toString(i);
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
	}*/
	
	
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
			chr = Integer.toString(i);
			
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
			chr = Integer.toString(i);
			
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
				SNVType.NonSynonymous.toLowerCase(), 
				SNVType.Synonymous.toLowerCase(), 
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
		
		String root = args[0]; //project directory
		switch (Integer.parseInt(args[1])) { //args[1] --> 'switchboard' parameter
			case 0:
				SeqPlatform platform = SeqPlatform.getPlatform(Integer.parseInt(args[2]));
				curateSNPCalls(root + "/naf-taf-inputs", root + "/snps", root + "/mouse_ears", platform); //args[2] --> 0::Illumina, 1::SOLiD
				break;
			case 1:
				segmentRegionsAllFiles(root + "/snps", root + "/regions", root + "/browser_tracks");
				
				//No need for this function now: 
				//pickupSomatics(root + "/snps", root + "/regions", root + "/curated_snps");
				break;
			case 2:
				//scoreRegions(root + "/regions", root + "/curated_snps", root + "/scored_regions");
				break;
			case 3:
				//genBrowserTracks(root + "/curated_snps", root + "/regions", root + "/browser_tracks");
				addScoreTracks(root + "/scored_regions", root + "/score_tracks"); //come back to this
				break;
			case 4:
				//getGeneEnrichment(root + "/curated_snps", root + "/gene_enrichment.csv");
				getGeneEnrichment(root + "/snps", root + "/gene_enrichment.csv");
				//Enrichment.genTopGeneLists(root + "/gene_enrichment.csv", root + "/gene_enrichment");
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
