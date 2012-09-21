import java.awt.Color;
import java.awt.Font;
import java.awt.geom.Ellipse2D;
import java.awt.geom.Ellipse2D.Double;
import java.io.BufferedWriter;
import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Date;

import lohcateEnums.Chrom;
import lohcateEnums.ClusterType;
import lohcateEnums.ColorPastel;
import lohcateEnums.SeqPlatform;
import lohcateEnums.VariantLocation;

import org.apache.commons.math3.distribution.BinomialDistribution;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYItemRenderer;
import org.jfree.chart.renderer.xy.XYShapeRenderer;
import org.jfree.chart.title.LegendTitle;
import org.jfree.data.xy.DefaultXYDataset;
import org.jfree.data.xy.XYDataset;

import com.carrotsearch.hppc.IntArrayList;

import shared.ArrayUtils;
import shared.FileOps;
import shared.GraphUtils;
import shared.IOUtils;
import shared.NumberUtils;
import shared.Utils;
import shared.Utils.FileExtensionAndDelimiter;


public class Clustering {
	
	// ========================================================================
	// STAGE 1: Curate SNP Calls and Cluster with DBSCAN
	// ========================================================================
	
	private static final float NAF_STRIP_EXPANDER = 2; //1.25f; //# of std. deviations to move away from the mean (when defining the thickness of the horizontal band containing HET ball, LOH sidelobes, &c.)
	private static final float HET_BALL_EPS = 0.035f; //DBSCAN parameters for HET ball / DUP wedge detection	
	private static final float DUP_WEDGE_LASSO = 0.015f;
	
	private static final float NON_HET_BALL_EPS = 0.065f;
	private static final int NON_HET_BALL_MINPTS = 30;
	
	private static final int HET_BALL_MINPTS = 100; //DBSCAN parameters for HET ball / DUP wedge detection
	private static final int DUP_WEDGE_MINPTS = 100;
	
	private static final float ClusterDiagonalLeeway = (float) 0.2;	
	
	public static final boolean Doing3D = false;
	public static final boolean UsePValuePlane = true;	
	public static final boolean MultipleTestingCorrect = false;
	private static final boolean AssignAmplifications = true;
	private static final boolean ForcePointsOnDiagonalAsNull = false;
	public static final float PValueBinDistAlpha_UpperPlaneThresh = 0.025f;
	public static float ScalingFactor = 1.0f;
	public static final float AmplificationThreshold = 2.5f;
	
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
		
		public float getValueNStandardDeviationsAway(float n) { return (mCountMean + (n * mStdDev)); }
		
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
				float vafNormal = Clustering.extractVAFNormal(rows.get(row), platform);
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
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub

	}

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

		String[] columnHeaders = new String[] { "chr", "pos", "n_vaf", "t_vaf", "dbsnp", "gene", "mutation_type", "germ_som", "cluster" };
		String headerStr = Utils.constructColumnDelimitedString(columnHeaders, fileExtDelim.mDelimiter, sb, true).toString();		

		int fileIndex = 0;
		for (File file : files) {			
			int indexOfSubstring = file.getName().indexOf(Script.GermlineSuffix);
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
				germlineSpecificVariantRows = Script.curateSNPCalls_removeHeaderLinesFromRows(germlineSpecificVariantRows);
				somaticSpecificVariantRows  = Script.curateSNPCalls_removeHeaderLinesFromRows(somaticSpecificVariantRows);

				// Create a combined set of rows
				ArrayList<String> allVariantRows = new ArrayList<String>(germlineSpecificVariantRows.size() + somaticSpecificVariantRows.size());
				allVariantRows.addAll(germlineSpecificVariantRows);
				int indexFirstSomaticRowInAllVariants = allVariantRows.size();  // save the size with just the germline variants added
				allVariantRows.addAll(somaticSpecificVariantRows);

				//Collections.sort(germlineSpecificVariantRows, LineComparatorTab);
				//Collections.sort(somaticSpecificVariantRows,  LineComparatorTab);
				Collections.sort(allVariantRows,              Script.LineComparatorTab);

				String outFilename = samplenameRoot + fileExtDelim.mExtension;
				String outFilenameFullPath = outDir + File.separator + outFilename; 

				BufferedWriter out = IOUtils.getBufferedWriter(outFilenameFullPath);
				IOUtils.writeToBufferedWriter(out, headerStr, true);

				// Now get the clusters
				float[] copyNumRatios = getRoughCopyNumberRatioPerSite(allVariantRows, getAvgCoverageRatioPerChrom(allVariantRows));
				double[] imbalancePValues = Clustering.getPValuesTumorImbalance(allVariantRows);
				ClusterType[] clusters = null;
				
				if (UsePValuePlane) {
					clusters = Clustering.getClusters_withPlane(allVariantRows, imbalancePValues, copyNumRatios, outFilenameFullPath, indexFirstSomaticRowInAllVariants, platform); //get cluster assignments (HET ball, LOH sidelobes, DUP wedge, &c.)
				} else {
					clusters = Clustering.getClusters_Old(allVariantRows, imbalancePValues, copyNumRatios, outFilenameFullPath, indexFirstSomaticRowInAllVariants, platform); //get cluster assignments (HET ball, LOH sidelobes, DUP wedge, &c.)
				}

				boolean toPrintWithCopyNum = false;
				if (toPrintWithCopyNum) {
					PrintStream outStream = IOUtils.getPrintStream(outFilenameFullPath + ".withCopyNum.txt");
					for (int r = 0; r < allVariantRows.size(); r++) {
						outStream.print(allVariantRows.get(r));					
						outStream.println("\t" + copyNumRatios[r] + "\t" + (copyNumRatios[r] * 2));
					}
					IOUtils.closePrintStream(outStream);
				}


				System.out.println("Got clusters");

				// Now initialize the data structure needed to plot
				int[] clusterTypeCounts = Utils.getClusterTypeCounts(clusters);
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
					ClusterType clusterType = clusters[i];
					int indexForClusterType = ++(clusterTypeIndex[clusterType.ordinal()]);  // increase the index of the data structure of the cluster type

					String[] germCols = strRow.split(Utils.TabStr);

					float vafNormal = Clustering.extractVAFNormal(germCols, platform);
					float vafTumor  = Clustering.extractVAFTumor (germCols, platform);
					clusterCoordinates[ clusterType.ordinal() ][0][indexForClusterType] = vafTumor;
					clusterCoordinates[ clusterType.ordinal() ][1][indexForClusterType] = vafNormal;

					String chromStr = (platform == SeqPlatform.Illumina) ? germCols[Script.Col_NAFTAFInput_Chrom] : Script.ChromPrefix + germCols[Script.Col_NAFTAFInput_Chrom].replace(Script.ChromPrefix, "");

					sb.setLength(0); // Clear the string builder
					sb.append(chromStr)
					.append(fileExtDelim.mDelimiter).append(germCols[Script.Col_NAFTAFInput_Position])
					.append(fileExtDelim.mDelimiter).append(vafNormal)
					.append(fileExtDelim.mDelimiter).append(vafTumor)
					.append(fileExtDelim.mDelimiter); //chr,pos,n_vaf,t_vaf

					String variantAnnotation = "";
					if (platform == SeqPlatform.SOLiD) {
						variantAnnotation = (germCols[7].indexOf(Script.NovelStr) >= 0) ? Script.NovelStr : Utils.NAStr;
						// Sidd: For the NAStr case, strangely, the variant base is n/a in the SOLiD naf-taf-inputs 
						// (and there's not much point in looking up the reference base's allele frequency)
					} else if (platform == SeqPlatform.Illumina) {
						String dbsnpStr = germCols[Script.Col_NAFTAFInput_DbSNPString]; 
						if (dbsnpStr.indexOf(Script.NovelStr) >= 0) {
							variantAnnotation = Script.NovelStr;
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
								gene         = germCols[Script.Col_NAFTAFInput_HugoSymbol];
								mutationType = germCols[Script.Col_NAFTAFInput_MutationType];
							}

							sb.append(fileExtDelim.mDelimiter).append(gene)
							.append(fileExtDelim.mDelimiter).append(mutationType)
							.append(fileExtDelim.mDelimiter).append(targetTissue)
							.append(fileExtDelim.mDelimiter).append(clusterType)
							//.append(fileExtDelim.mDelimiter).append(copyNumRatios[i] * DefaultDiploidCopyNumber);
							//.append(fileExtDelim.mDelimiter).append(imbalancePValues[i]);
							;


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
				Clustering.plotVAFComparison(xyDataset, vafComparisonPlotDir + File.separator + samplenameRoot, samplenameRoot);

			}
		}		
	}

	/** Plots the VAF (variant allele frequency) of the normal tissue comapred to the tumor tissue. */
	public static void plotVAFComparison(XYDataset xyDataset, String outFilenameRoot, String sampleName) {
		String title = "VAF Comparison: " + sampleName;
		String xAxisLabel = "VAF Tumor";
		String yAxisLabel = "VAF Normal";
		JFreeChart theChart = 
				ChartFactory.createScatterPlot(title, xAxisLabel, yAxisLabel, xyDataset, PlotOrientation.VERTICAL, true, false, false);
		XYPlot xyPlot = (XYPlot) theChart.getPlot();
		
		XYItemRenderer itemRenderer = Clustering.getXYItemRendererHelper(5);		
		itemRenderer.setSeriesPaint(ClusterType.Amp.ordinal(), ColorPastel.Dark_Red.getColor());
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
	public static double[] getPValuesTumorImbalance(ArrayList<String> rows) {
		double[] pValues = new double[rows.size()];
		
		for (int row = 0; row < rows.size(); row++) {
			String line = rows.get(row);
			
			Chrom chrom =      Chrom.getChrom( Utils.extractNthColumnValue(line, Script.Col_NAFTAFInput_Chrom,      Utils.FileExtensionTSV.mDelimiter) );
			int varCovgTumor = Integer.parseInt(Utils.extractNthColumnValue(line, Script.Col_NAFTAFInput_VariantCoverageTumor, Utils.FileExtensionTSV.mDelimiter));
			int totCovgTumor = Integer.parseInt(Utils.extractNthColumnValue(line, Script.Col_NAFTAFInput_TotalCoverageTumor, Utils.FileExtensionTSV.mDelimiter));
			
			int maxRefOrVarCovg = Math.max(varCovgTumor, totCovgTumor - varCovgTumor);
			pValues[row] = NumberUtils.cumulativeProbabilitySuccess(totCovgTumor, maxRefOrVarCovg, 0.5);
		}
		
		return pValues;
	}

	// ========================================================================
	public static float[] getRoughCopyNumberRatioPerSite(ArrayList<String> rows, float[] ratioPerChrom) {
		float[] copyNumRatio    =   new float[rows.size()];
		boolean[] isSomaticSite = new boolean[rows.size()];
		
		String prevGene = "";
		int rowOfFirstInstanceOfGene = -1;
		
		float ratioSum = 0;
		int numRowsWithSameGene = 0;
		
		for (int row = 0; row < rows.size(); row++) {
			String line = rows.get(row);
			
			Chrom chrom =      Chrom.getChrom( Utils.extractNthColumnValue(line, Script.Col_NAFTAFInput_Chrom,      Utils.FileExtensionTSV.mDelimiter) );
			String currentGene =               Utils.extractNthColumnValue(line, Script.Col_NAFTAFInput_HugoSymbol, Utils.FileExtensionTSV.mDelimiter);
			
			int covgNormal = Integer.parseInt( Utils.extractNthColumnValue(line, Script.Col_NAFTAFInput_TotalCoverageNormal, Utils.FileExtensionTSV.mDelimiter) );
			int covgTumor  = Integer.parseInt( Utils.extractNthColumnValue(line, Script.Col_NAFTAFInput_TotalCoverageTumor,  Utils.FileExtensionTSV.mDelimiter) );
			
			isSomaticSite[row] = (covgNormal <= 0);
			float ratio = isSomaticSite[row] ? Script.TumorNormalRatioOfSomaticSite : ((float) covgTumor / (float) covgNormal);
			
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
					float averageRatio = onlySomaticRowsForGene ? Script.TumorNormalRatioOfSomaticSite : 
						((numRowsWithSameGene == 1) ? Script.DefaultTumorNormalRatio : 
													  ((ratioSum / (float) numRowsWithSameGene) / ratioPerChrom[chrom.ordinal()]));
	
					for (int i = rowOfFirstInstanceOfGene; i <= row; i++) {					
						copyNumRatio[i] = isSomaticSite[i] ? Script.TumorNormalRatioOfSomaticSite : averageRatio;
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
		Arrays.fill(ratios, Script.DefaultTumorNormalRatio);  // Fill with 1.0 since that's an equal ratio
		
		int totalCoverageTumor = 0;
		int totalCoverageNormal = 0;
		Chrom prevChrom = null;
		
		for (int row = 0; row < rows.size(); row++) {
			String line = rows.get(row);
			
			Chrom chrom = Chrom.getChrom( Utils.extractNthColumnValue(line, Script.Col_NAFTAFInput_Chrom, Utils.FileExtensionTSV.mDelimiter) );
			
			if ((prevChrom == null) || (chrom != prevChrom)) {
				if (prevChrom != null) {
					if (totalCoverageNormal > 0) {  // avoid div/0
						ratios[prevChrom.ordinal()] = (float) totalCoverageTumor / (float) totalCoverageNormal;
					}
				}
				
				totalCoverageTumor = totalCoverageNormal = 0;
				prevChrom = chrom;
			}
			
			int covgNormal = Integer.parseInt( Utils.extractNthColumnValue(line, Script.Col_NAFTAFInput_TotalCoverageNormal, Utils.FileExtensionTSV.mDelimiter) );
			int covgTumor  = Integer.parseInt( Utils.extractNthColumnValue(line, Script.Col_NAFTAFInput_TotalCoverageTumor,  Utils.FileExtensionTSV.mDelimiter) );
			
			totalCoverageNormal += covgNormal; 
			totalCoverageTumor  += covgTumor; 			
		}
		
		// We finished the loop without writing the ratio of the last chromosome
		if (totalCoverageNormal > 0) {
			ratios[prevChrom.ordinal()] = (float) totalCoverageTumor / (float) totalCoverageNormal;
		}
		
		return ratios;
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
	public static ClusterType[] getClusters_Old(ArrayList<String> rows, double[] imbalancePValues, float[] copyNumRatios, String outFilenameFullPath, int startingRowGermlineOrSomaticOrAll, SeqPlatform platform) {

		// Get the allele frequency statistics, and adjust the frames based on the resulting standard deviation
		AlleleFrequencyStatsForSample afStatsSample = new AlleleFrequencyStatsForSample();
		afStatsSample.tabulateAndPerformStatistics(rows, platform);
		float vafNormalFrameAdjustedLower = afStatsSample.mCountMean - (Clustering.NAF_STRIP_EXPANDER * afStatsSample.mStdDev);
		float vafNormalFrameAdjustedUpper = afStatsSample.mCountMean + (Clustering.NAF_STRIP_EXPANDER * afStatsSample.mStdDev);

		//apply DBScan to points within NAF frame		
		//double scalingFactor = DefaultDiploidCopyNumber;		
		ArrayList<Floint> points = Clustering.getValidPointsListForClustering(rows, imbalancePValues, ScalingFactor, platform, vafNormalFrameAdjustedLower, vafNormalFrameAdjustedUpper);

		System.out.println("Begin clustering algorithm: " + (new Date()).toString());

		//DBScanFaster dbscanner = new DBScanFaster(points, HET_BALL_EPS + DUP_WEDGE_LASSO, DUP_WEDGE_MINPTS, 0, 0, 1, 1); //parameters for capturing HET ball//KDBSCAN(param, 10, 0.0325f);
		DBScanFaster dbscanner = new DBScanFaster(points, Clustering.HET_BALL_EPS, Clustering.HET_BALL_MINPTS, 0, 0, 1, 1); //parameters for capturing HET ball//KDBSCAN(param, 10, 0.0325f);

		dbscanner.cluster();
		int clusterIDofHetBall = dbscanner.getCentralClusterID();
		Boolean[] pointsWithinRadius = dbscanner.getPointsWithinMinRadiusOfCluster(clusterIDofHetBall);
		int[] clusterAssignments = dbscanner.getClustAssignments();  // save and cache	

		int[] clusterTypeIDsFromAlgorithm = new int[ClusterType.values().length];		
		clusterTypeIDsFromAlgorithm[ClusterType.Amp.ordinal()]  = -1;
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
			}
		}
		DBScanFaster dbscannerNonHet = new DBScanFaster(nonHetPoints, Clustering.NON_HET_BALL_EPS, Clustering.NON_HET_BALL_MINPTS, 0, 0, 1, 1);
		dbscannerNonHet.cluster();
		int[] clusterAssignmentsNonHet = dbscannerNonHet.getClustAssignments();

		System.out.println("End clustering algorithm: " + (new Date()).toString());

		int nonHetIndex = -1;
		for (int i = 0; i < clusterAssignments.length; i++) {

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
					//clusterAssignments[i] = clusterTypeIDsFromAlgorithm[ClusterType.Dup.ordinal()];  // assign for lower thresholds

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
					}
				}
			 */

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
					clusterAssignments[i] = clusterTypeIDsFromAlgorithm[ClusterType.Dup.ordinal()];   // only change if we're in a het ball region
				}
			}
			*/
			
			
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
					if (pointOnDiagonal(points.get(i), ClusterDiagonalLeeway)) {
						clusterAssignments[i] = clusterTypeIDsFromAlgorithm[ClusterType.Null.ordinal()];
					} else {
						clusterAssignments[i] = clusterAssignmentsNonHet[nonHetIndex];   
					}
				} else {
					clusterAssignments[i] = clusterTypeIDsFromAlgorithm[ClusterType.Amp.ordinal()];   // only change if we're in a het ball region
				}
			}

		}

		// Now assign the cluster types

		return Clustering.assignClusters(rows, platform, vafNormalFrameAdjustedLower, vafNormalFrameAdjustedUpper, clusterAssignments, clusterTypeIDsFromAlgorithm);
	}
	
	
	// ========================================================================
	public static ClusterType[] getClusters_withPlane(ArrayList<String> rows, double[] imbalancePValues, float[] copyNumRatios, String outFilenameFullPath, int startingRowGermlineOrSomaticOrAll, SeqPlatform platform) {

		// Get the allele frequency statistics, and adjust the frames based on the resulting standard deviation
		AlleleFrequencyStatsForSample afStatsSample = new AlleleFrequencyStatsForSample();
		afStatsSample.tabulateAndPerformStatistics(rows, platform);
		float vafNormalFrameAdjustedLower = afStatsSample.getValueNStandardDeviationsAway(-NAF_STRIP_EXPANDER); 								
		float vafNormalFrameAdjustedUpper = afStatsSample.getValueNStandardDeviationsAway( NAF_STRIP_EXPANDER);

		//apply DBScan to points within NAF frame
		ArrayList<Floint> pointsUpperPlane = new ArrayList<Floint>();
		ArrayList<Floint> pointsLowerPlane = new ArrayList<Floint>();
		int[] indexMapToUpperLowerPlane    = new int[rows.size()];
		IntArrayList mapToRowsFromUpper = new IntArrayList(rows.size());
		IntArrayList mapToRowsFromLower = new IntArrayList(rows.size());

		getValidPointsForClustering(rows, pointsUpperPlane, pointsLowerPlane, indexMapToUpperLowerPlane, 
				                    mapToRowsFromUpper, mapToRowsFromLower, 
				                    imbalancePValues, ScalingFactor, copyNumRatios, 
				                    platform, 
				                    vafNormalFrameAdjustedLower, vafNormalFrameAdjustedUpper);		

		System.out.println("Begin clustering algorithm: " + (new Date()).toString());

		//DBScanFaster dbscanner = new DBScanFaster(points, HET_BALL_EPS + DUP_WEDGE_LASSO, DUP_WEDGE_MINPTS, 0, 0, 1, 1); //parameters for capturing HET ball//KDBSCAN(param, 10, 0.0325f);
		
		// Scan the lower plane
		DBScanFaster dbscannerLowerPlane = new DBScanFaster(pointsLowerPlane, Clustering.HET_BALL_EPS, Clustering.HET_BALL_MINPTS, 0, 0, 1, 1); //parameters for capturing HET ball//KDBSCAN(param, 10, 0.0325f);

		dbscannerLowerPlane.cluster();
		int clusterIDofHetBall       = dbscannerLowerPlane.getCentralClusterID();
		Boolean[] pointsWithinRadius = dbscannerLowerPlane.getPointsWithinMinRadiusOfCluster(clusterIDofHetBall);
		int[] clusterAssignments     = dbscannerLowerPlane.getClustAssignments();  // save and cache	

		int[] clusterTypeIDsFromAlgorithm = new int[ClusterType.values().length];		
		clusterTypeIDsFromAlgorithm[ClusterType.Amp.ordinal()]  = -1;
		clusterTypeIDsFromAlgorithm[ClusterType.Null.ordinal()] = -2;
		clusterTypeIDsFromAlgorithm[ClusterType.HETGermline.ordinal()] = clusterIDofHetBall;
		clusterTypeIDsFromAlgorithm[ClusterType.Noise.ordinal()] = dbscannerLowerPlane.getClusterIDOfNoise();


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
		ArrayList<Floint> nonHetPoints = new ArrayList<Floint>(pointsLowerPlane.size());
		boolean[] isNonHetPoint                  = new boolean[pointsLowerPlane.size()];		
		
		for (int i = 0; i < clusterAssignments.length; i++) {
			isNonHetPoint[i] = (clusterAssignments[i] != clusterIDofHetBall); 
			if (isNonHetPoint[i]) {
				nonHetPoints.add(pointsLowerPlane.get(i));
			}
		}
		
		DBScanFaster dbscannerLowerPlaneNonHet = new DBScanFaster(nonHetPoints, Clustering.NON_HET_BALL_EPS, Clustering.NON_HET_BALL_MINPTS, 0, 0, 1, 1);
		dbscannerLowerPlaneNonHet.cluster();
		int[] clusterAssignmentsNonHet = dbscannerLowerPlaneNonHet.getClustAssignments();

		// UPPER PLANE
		DBScanFaster dbscannerUpperPlane = new DBScanFaster(pointsUpperPlane, Clustering.NON_HET_BALL_EPS, Clustering.NON_HET_BALL_MINPTS, 0, 0, 1, 1);
		dbscannerUpperPlane.cluster();
		int[] clusterAssignmentsUpperPlane = dbscannerUpperPlane.getClustAssignments();
		for (int ind = 0; ind < clusterAssignmentsUpperPlane.length; ind++) {
			if (clusterAssignmentsUpperPlane[ind] == clusterIDofHetBall) {
				clusterAssignmentsUpperPlane[ind] = 
					(clusterIDofHetBall == DBSCAN2.ClusterIDOfNoise + 1) ? clusterIDofHetBall + 1 : clusterIDofHetBall - 1;
			}			
			int indexIntoOriginalRows = mapToRowsFromUpper.get(ind);
			float copyNum = copyNumRatios[indexIntoOriginalRows] * Script.DefaultDiploidCopyNumber;
			if (AssignAmplifications && (copyNum > AmplificationThreshold)) {
				clusterAssignmentsUpperPlane[ind] = clusterTypeIDsFromAlgorithm[ClusterType.Amp.ordinal()];
			}
			
			if (ForcePointsOnDiagonalAsNull) {
				if (pointOnDiagonal(pointsUpperPlane.get(ind), ClusterDiagonalLeeway)) {
					clusterAssignmentsUpperPlane[ind] = clusterTypeIDsFromAlgorithm[ClusterType.Null.ordinal()];
				}
			}
		}
		
		System.out.println("End clustering algorithm: " + (new Date()).toString());

		int nonHetIndex = -1;
		for (int i = 0; i < clusterAssignments.length; i++) {

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
					//clusterAssignments[i] = clusterTypeIDsFromAlgorithm[ClusterType.Dup.ordinal()];  // assign for lower thresholds

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
					}
				}
			 */

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
					clusterAssignments[i] = clusterTypeIDsFromAlgorithm[ClusterType.Dup.ordinal()];   // only change if we're in a het ball region
				}
			}
			*/
			
			
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
					if (pointOnDiagonal(pointsLowerPlane.get(i), ClusterDiagonalLeeway)) {
						clusterAssignments[i] = clusterTypeIDsFromAlgorithm[ClusterType.Null.ordinal()];
					} else {
						clusterAssignments[i] = clusterAssignmentsNonHet[nonHetIndex];   
					}
				} else {
					clusterAssignments[i] = clusterTypeIDsFromAlgorithm[ClusterType.Amp.ordinal()];   // only change if we're in a het ball region
				}
			}

		}

		
		
		// Now merge the two planes into one array
		ArrayUtils.IntArray clusterAssignmentsFinal = new ArrayUtils.IntArray(pointsLowerPlane.size() + pointsUpperPlane.size());  		
		for (int row = 0; row < rows.size(); row++) {
			if (indexMapToUpperLowerPlane[row] != Integer.MAX_VALUE) {				
				if (indexMapToUpperLowerPlane[row] >= 0) {
					int trueIndex = indexMapToUpperLowerPlane[row];
					clusterAssignmentsFinal.add(clusterAssignmentsUpperPlane[trueIndex]);
				} else {
					int trueIndex = -(indexMapToUpperLowerPlane[row] + 1);
					clusterAssignmentsFinal.add(clusterAssignments[trueIndex]);
				}
			}
		}
		

		return Clustering.assignClusters(rows, platform, vafNormalFrameAdjustedLower, vafNormalFrameAdjustedUpper, clusterAssignmentsFinal.mArray, clusterTypeIDsFromAlgorithm);
	}

	// ========================================================================
	// Convenience private function to break up caller function
	static ClusterType[] assignClusters(ArrayList<String> rows, 
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
			float vafNormal = Clustering.extractVAFNormal(line, platform);
			boolean vafInRangeNormal = Utils.inRangeLowerExclusive(vafNormal, vafBoundLower, vafBoundUpper);

			if (/*(row >= startingRowSomatic) &&*/ (!vafInRangeNormal)) {
				// The vafNormal is either very low (homozygous reference) or very high (homozygous common variant).
				// We do some very simple decision making now (which should be replaced by formal clustering later)
				// to partition the calls.
				float justBelowZero = -0.0001f;				
				float hetBoundaryLower = 0.3333f;
				float hetBoundaryUpper = 0.6667f;
				float vafTumor = Clustering.extractVAFTumor(line, platform);
				if (Utils.inRangeLowerExclusive(vafNormal, justBelowZero, vafBoundLower)) {
					// We are equal to or below the lower frame boundary
					if (vafNormal <= 0.05 && vafTumor >= 0.10) {
						returnClusters[row] = ClusterType.HETSomatic;
					} else if (vafTumor <= hetBoundaryLower) {
						// Normal: AA, Tumor: AA [Thus homozygous reference in both, no events]
						returnClusters[row] = ClusterType.Null;
					} else if (vafTumor > hetBoundaryUpper) {
						// Normal: AA, Tumor: BB or CC [made by: AA -> AB or AC (somatic het mutation) -> B or C (LOH, loss of A)]
						returnClusters[row] = ClusterType.LOH;
					} else {
						// Normal: AA, Tumor: AB [made by: AA -> AB (somatic het mutation)
						returnClusters[row] = ClusterType.Null;
					}					
				} else if (Utils.inRangeLowerExclusive(vafNormal, vafBoundUpper, Script.MaxVariantAlleleFrequency)) {
					// We are above the upper frame boundary
					if (vafNormal >= (1.0 - 0.05) && vafTumor >= 0.10 && vafTumor <= 0.90) {
						returnClusters[row] = ClusterType.HETSomatic;
					} else if (vafTumor <= hetBoundaryLower) {
						// Normal: BB, Tumor: AA [made by: BB -> AB (reverse somatic het mutation) -> A (LOH, loss of B)]
						returnClusters[row] = ClusterType.LOH;
					} else if (vafTumor > hetBoundaryUpper) {
						// Normal: BB, Tumor: BB or CC (ambiguous until we know exact variant for tumor)
						// TODO - Leave as Null for now, but will need to change later to resolve the
						// ambiguity mentioned above
						returnClusters[row] = ClusterType.Null;
					} else {
						// Normal: BB, Tumor: AB or CB [made by: BB -> AB (reverse somatic het mutation) or BB -> CB (somatic het mutation)
						returnClusters[row] = ClusterType.Null;
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

					} else if (assignedClusterID == clusterTypeIDsFromAlgorithm[ClusterType.Amp.ordinal()]) {
						returnClusters[row] = ClusterType.Amp; //DUP

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
	private static void getValidPointsForClustering(ArrayList<String> rows,                // input
													ArrayList<Floint> pointsUpperPlane,    // output
													ArrayList<Floint> pointsLowerPlane,    // output
													int[] indexMapToUpperLowerPlane,       // output - if positive, to upper plane, if negative, to lower plane, if does not map, contains Integer.MAX_VALUE													
													IntArrayList mapToRowsFromUpper, 
													IntArrayList mapToRowsFromLower,
													double[] verticalFactor, double scalingFactor,
													float[] copyNumRatios,
													SeqPlatform platform, 
													float vafBoundLower, float vafBoundUpper) {
		
		pointsUpperPlane.clear();
		pointsLowerPlane.clear();
		double alphaAdjusted = MultipleTestingCorrect ? 
				(PValueBinDistAlpha_UpperPlaneThresh * 0.008 /* Light FDR */) : PValueBinDistAlpha_UpperPlaneThresh;
		for (int row = 0; row < rows.size(); row++) {
			String line = rows.get(row);			
			float vafNormal = Clustering.extractVAFNormal(line, platform);
			if (Utils.inRangeLowerExclusive(vafNormal, vafBoundLower, vafBoundUpper)) {
				Floint thePoint = new Floint(Clustering.extractVAFTumor(line, platform), vafNormal, (float) (verticalFactor[row] * scalingFactor));
				
				// Assume p-value as vertical row factor
				if (verticalFactor[row] <= alphaAdjusted && !isCopyNumInDiploidRange(copyNumRatios[row])) {
					indexMapToUpperLowerPlane[row] =  pointsUpperPlane.size();
					pointsUpperPlane.add(thePoint);
					mapToRowsFromUpper.add(row);
					//System.out.printf("(%g, %g, %g)\n", thePoint.mX, thePoint.mY, thePoint.mZ);										
				} else {
					indexMapToUpperLowerPlane[row] = -pointsLowerPlane.size() - 1;
					pointsLowerPlane.add(thePoint);
					mapToRowsFromLower.add(row);
					
				}
				
			} else {				
				indexMapToUpperLowerPlane[row] = Integer.MAX_VALUE;
			}
		}
	}
	
	// ========================================================================
	private static boolean isCopyNumInDiploidRange(float copyNumRatio) {
		float copyNum = copyNumRatio * Script.DefaultDiploidCopyNumber;
		float threshold = 0.2f;
		return ((Script.DefaultDiploidCopyNumber - threshold) < copyNum && copyNum < (Script.DefaultDiploidCopyNumber + threshold));
	}
	
	// ========================================================================
	// Convenience function 
	static ArrayList<Floint> 
	getValidPointsListForClustering(ArrayList<String> rows, double[] verticalFactor, double scalingFactor, SeqPlatform platform, float vafBoundLower, float vafBoundUpper) {
		ArrayList<Floint> points = new ArrayList<Floint>(rows.size());

		for (int row = 0; row < rows.size(); row++) {
			String line = rows.get(row);			
			float vafNormal = Clustering.extractVAFNormal(line, platform);
			if (Utils.inRangeLowerExclusive(vafNormal, vafBoundLower, vafBoundUpper)) {	
				Floint thePoint = new Floint(Clustering.extractVAFTumor(line, platform), vafNormal, (float) (verticalFactor[row] * scalingFactor));
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
	static float extractVAFNormal(String[] components, SeqPlatform platform) {
		switch(platform) {
		case Illumina: return  Float.parseFloat(components[Script.Col_NAFTAFInput_VariantRatioNormal]);
		case SOLiD:    return (Float.parseFloat(components[6]) / Float.parseFloat(components[5]));
		}
		return Float.NaN;
	}

	// ========================================================================
	// Extracts and calculates the variant allele frequency in the normal depending on the platform tab-delimited file
	static float extractVAFNormal(String line, SeqPlatform platform) {
		switch(platform) {
		case Illumina: return  Float.parseFloat(Utils.extractNthColumnValue(line, Script.Col_NAFTAFInput_VariantRatioNormal, Utils.TabStr));
		case SOLiD:    return (Float.parseFloat(Utils.extractNthColumnValue(line, 6, Utils.TabStr)) / 
				               Float.parseFloat(Utils.extractNthColumnValue(line, 5, Utils.TabStr)));
		}
		return Float.NaN;		
	}

	// ========================================================================
	static float extractVAFTumor(String[] components, SeqPlatform platform) {
		switch(platform) {
		case Illumina: return  Float.parseFloat(components[Script.Col_NAFTAFInput_VariantRatioTumor]);
		case SOLiD:    return (Float.parseFloat(components[4]) / Float.parseFloat(components[3]));
		}
		return Float.NaN;
	}

	// ========================================================================
	// Extracts and calculates the variant allele frequency in the tumor depending on the platform tab-delimited file
	static float extractVAFTumor(String line, SeqPlatform platform) {
		switch(platform) {
		case Illumina: return  Float.parseFloat(Utils.extractNthColumnValue(line, Script.Col_NAFTAFInput_VariantRatioTumor, Utils.TabStr));
		case SOLiD:    return (Float.parseFloat(Utils.extractNthColumnValue(line, 4, Utils.TabStr)) / 
				               Float.parseFloat(Utils.extractNthColumnValue(line, 3, Utils.TabStr)));
		}
		return Float.NaN;		
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

}
