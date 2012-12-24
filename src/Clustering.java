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

import javax.swing.text.html.HTMLDocument.HTMLReader.IsindexAction;

import lohcateEnums.Chrom;
import lohcateEnums.ClusterType;
import lohcateEnums.ColorPastel;
import lohcateEnums.SeqPlatform;
import lohcateEnums.VariantLocation;

import org.apache.commons.math3.distribution.BinomialDistribution;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYItemRenderer;
import org.jfree.chart.renderer.xy.XYShapeRenderer;
import org.jfree.chart.title.LegendTitle;
import org.jfree.data.xy.DefaultXYDataset;
import org.jfree.data.xy.XYDataset;

import com.carrotsearch.hppc.DoubleArrayList;
import com.carrotsearch.hppc.IntArrayList;
import com.martiansoftware.jsap.FlaggedOption;
import com.martiansoftware.jsap.JSAP;
import com.martiansoftware.jsap.JSAPResult;

import shared.ArgumentParserUtils;
import shared.ArrayUtils;
import shared.ArrayUtils.DynamicBucketCounter;
import shared.FileOps;
import shared.GraphUtils;
import shared.IOUtils;
import shared.NumberUtils;
import shared.PrimitiveWrapper;
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
	
	public static final boolean Doing3D = false;	
	public static final boolean UsePValuePlane = true;	
		
	public static final boolean UseBidrectionalAdditiveOffset = true;
	public static final boolean AllowGeneWithOneSite = true;
	
	private static final boolean ForcePointsOnDiagonalAsNull = true;
	private static final float ClusterDiagonalLeeway = (float) 0.2;	
	
	public static final float PValueBinDistAlpha_UpperPlaneThresh = 0.025f;
	public static float ScalingFactor = 1.0f;			
	public static final float ExpectedVAFNormalTrisomy = 2.0f / 3.0f;
	
	public static boolean correctAllelicBias() { return !ParamsBool.IgnoreAllelicBias.getValue(); } 
	
	// ========================================================================
	// PARAMETER Registration
	// ========================================================================
	
	public static enum ParamsNum {
		AmplificationThreshold(2.3f, 'a', "ampThresh"),
		DeletionThreshold     (1.7f, 'd', "delThresh"),
		GermlineTrisomyThreshold(1.48f, JSAP.NO_SHORTFLAG, "germlineAneuploidyGainThreshold")
		;
		
		private Number mDefaultValue;
		private Number mValue;
		private char mFlagShort;
		private String mFlagLong;
		
		private ParamsNum(Number defaultValue, char flagShort, String flagLong) {
			mValue = mDefaultValue = defaultValue;
			mFlagShort = flagShort;
			mFlagLong  = flagLong;
		}
		
		public Number getValue() { return mValue; } 
		public Number getDefault() { return mDefaultValue; }
		public char getShortFlag() { return mFlagShort; }
		public String getLongFlag() { return mFlagLong; }
		public void setValue(Number theValue) { mValue = theValue; }
		
	}
	
	// ========================================================================
	public static enum ParamsBool {
		IgnoreAllelicBias(false, JSAP.NO_SHORTFLAG, "ignoreAllelicBias"),
		IgnoreMultipleTesting(false, JSAP.NO_SHORTFLAG, "ignoreMultipleTesting")
		;
		
		private boolean mValue;
		private boolean mDefaultValue;
		private char mFlagShort;
		private String mFlagLong;

		private ParamsBool(boolean defaultValue, char flagShort, String flagLong) {
			mValue = mDefaultValue = defaultValue;
			mFlagShort = flagShort;
			mFlagLong  = flagLong;
		}
		
		public boolean getValue() { return mValue; } 
		public boolean getDefault() { return mDefaultValue; }
		public char getShortFlag() { return mFlagShort; }
		public String getLongFlag() { return mFlagLong; }
		public void setValue(boolean value) { mValue = value; }
	}
	
	
	// ========================================================================
	public static JSAP registerClusteringParameters(JSAP jsap) {
		jsap = (jsap == null) ? new JSAP() : jsap;
		
		//ArgumentParserUtils.createSwitch(ParamsBool.IgnoreAllelicBias.name(), ParamsBool.IgnoreAllelicBias.getShortFlag(), ParamsBool.IgnoreAllelicBias.getLongFlag(), "Specifies whether allelic bias is to be ignored and not corrected", jsap);
		ArgumentParserUtils.createSwitch(ParamsBool.IgnoreMultipleTesting.name(), ParamsBool.IgnoreMultipleTesting.getShortFlag(), ParamsBool.IgnoreMultipleTesting.getLongFlag(), "Specifies whether multiple testing correction is not to be done", jsap);
		
		FlaggedOption ampThreshold = new FlaggedOption(
				ParamsNum.AmplificationThreshold.name()).setStringParser(JSAP.DOUBLE_PARSER)
				.setDefault(ParamsNum.AmplificationThreshold.getDefault().toString())
				.setShortFlag(ParamsNum.AmplificationThreshold.getShortFlag())
				.setLongFlag(ParamsNum.AmplificationThreshold.getLongFlag())
				.setUsageName("Copy_Number_Threshold_for_Amplification");
		ArgumentParserUtils.registerJSAPParameter(jsap, ampThreshold);
		
		FlaggedOption delThreshold = new FlaggedOption(
				ParamsNum.DeletionThreshold.name()).setStringParser(JSAP.DOUBLE_PARSER)
				.setDefault(ParamsNum.DeletionThreshold.getDefault().toString())
				.setShortFlag(ParamsNum.DeletionThreshold.getShortFlag())
				.setLongFlag(ParamsNum.DeletionThreshold.getLongFlag())
				.setUsageName("Copy_Number_Threshold_for_Deletion");
		ArgumentParserUtils.registerJSAPParameter(jsap, delThreshold);


		FlaggedOption germlineChromGainThreshold = new FlaggedOption(
				ParamsNum.GermlineTrisomyThreshold.name()).setStringParser(JSAP.DOUBLE_PARSER)
				.setDefault(ParamsNum.GermlineTrisomyThreshold.getDefault().toString())
				.setShortFlag(ParamsNum.GermlineTrisomyThreshold.getShortFlag())
				.setLongFlag(ParamsNum.GermlineTrisomyThreshold.getLongFlag())
				.setUsageName("Copy_Number_Threshold_for_Germline_Chromosomal_Gain");
		ArgumentParserUtils.registerJSAPParameter(jsap, germlineChromGainThreshold);		
		
		return jsap;
	}
	
	// ========================================================================
	public static void configureParameters(JSAPResult config) {
		ParamsNum.AmplificationThreshold.setValue(  config.getDouble(ParamsNum.AmplificationThreshold.name()));
		ParamsNum.GermlineTrisomyThreshold.setValue(config.getDouble(ParamsNum.GermlineTrisomyThreshold.name()));
		//ParamsBool.IgnoreAllelicBias.setValue(     config.getBoolean(ParamsBool.IgnoreAllelicBias.name()));
	}
	
	// ========================================================================
	// INNER CLASS
	// ========================================================================
	/** An inner class that calculates statistics for a set of sites (with variant 
	 *  allele frequencies information for tumor and matched normal). */
	public static class AlleleFrequencyStatsForSample {
		
		//NAF {het, loh, dup} FRAME definition via peak detection and parameter-tuned standard deviation expansion
		//we have to avoid the often hugely dense peak of homozygous mutations (AF > 0.8) and the occasionally hugely dense peak of neg. tail noise / somatics / &c. (AF < 0.2)
		public static float VAFNormalFrameLower = correctAllelicBias() ? 0.1f : 0.2f;
		public static float VAFNormalFrameUpper = correctAllelicBias() ? 0.9f : 0.8f; 
		public static float BinSize             = 0.025f; //smoothing parameter
		
		int   mNumBins;
		int[]   mBinCount;    // The bins in which counts are binned and stored
		float[] mBinValue;    // The value each bin represents
		
		// Statistics values
		float mCountMean = -1;   
		float mVariance = -1;
		float mStdDev = -1;
		
		public AlleleFrequencyStatsForSample() {
			mNumBins = (int) ((VAFNormalFrameUpper - VAFNormalFrameLower) / BinSize) + 1;
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
				mBinValue[i] = VAFNormalFrameLower + ((i + 1) * BinSize);
			}
		}
		
		private void tallyVariantAlleleFrequenciesIntoBins(ArrayList<String> rows, SeqPlatform platform, boolean clearBinCountBins) {
			if (clearBinCountBins) {
				Arrays.fill(mBinCount, 0);
			}
			
			// First, tally the variant allele frequencies into bins
			for (int row = 0; row < rows.size(); row++) {
				float vafNormal = Clustering.extractVAFNormal(rows.get(row), platform);
				if (Utils.inRangeLowerExclusive(vafNormal, VAFNormalFrameLower, VAFNormalFrameUpper)) {
					int binNumber = (int) ((vafNormal - VAFNormalFrameLower) / BinSize);
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
	public static void classifySites(String inDir, String allelicBiasInFile, String outDir, String vafComparisonPlotDir, String vafWaterfallPlotDir, String copyNumberPlotDir, SeqPlatform platform) {
		File[] files = (new File(inDir)).listFiles();
		StringBuilder sb = new StringBuilder(8192);
		Utils.FileExtensionAndDelimiter fileExtDelim = Utils.FileExtensionTSV;		

		// TODO -- make column names static constants
		System.out.println("Reading Allelic Bias file...");
		AllelicBiasTable allelicBiasTable = correctAllelicBias() ? AllelicBiasTable.readFileAndConstructTable(allelicBiasInFile, 3, 4) : null;		
		System.out.println("Finished Reading Allelic Bias file...");
		
		// Create output directory
		IOUtils.createDirectoryPath(outDir, false);
		IOUtils.createDirectoryPath(vafComparisonPlotDir, false);
		IOUtils.createDirectoryPath(vafWaterfallPlotDir, false);
		IOUtils.createDirectoryPath(copyNumberPlotDir, false);
		

		String[] columnHeaders = new String[] { "chr", "pos", "VAF_Normal", "VAF_Tumor", "dbsnp", "gene", "mutation_type", "germ_som", "cluster" };
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
				System.out.println("\tRead All Lines...");
				System.out.printf("\tNum Sites Total: %d\n", germlineSpecificVariantRows.size() + somaticSpecificVariantRows.size());				
				
				//String headerStringGermline = germlineSpecificVariantRows.get(0);
				//String headerStringSomatic  = somaticSpecificVariantRows.get(0);

				//upstream pipelines will randomly spit header lines into the 
				// middle of a naf-taf-input file. we're just avoiding those
				germlineSpecificVariantRows = Script.curateSNPCalls_removeHeaderLinesFromRows(germlineSpecificVariantRows);
				somaticSpecificVariantRows  = Script.curateSNPCalls_removeHeaderLinesFromRows(somaticSpecificVariantRows);
				System.out.printf("\tNum Sites Retained after Header and GC Removal: %d\n", germlineSpecificVariantRows.size() + somaticSpecificVariantRows.size());
				
				// Create a combined set of rows
				ArrayList<String> allVariantRows = new ArrayList<String>(germlineSpecificVariantRows.size() + somaticSpecificVariantRows.size());
				allVariantRows.addAll(germlineSpecificVariantRows);
				int indexFirstSomaticRowInAllVariants = allVariantRows.size();  // save the size with just the germline variants added
				allVariantRows.addAll(somaticSpecificVariantRows);

				//Collections.sort(germlineSpecificVariantRows, LineComparatorTab);
				//Collections.sort(somaticSpecificVariantRows,  LineComparatorTab);
				Collections.sort(allVariantRows,              Script.LineComparatorTab);
				System.out.println("\tAll Sites Finished Sorting...");
				

				String outFilename = samplenameRoot + fileExtDelim.mExtension;
				String outFilenameFullPath = outDir + File.separator + outFilename; 

				BufferedWriter out = IOUtils.getBufferedWriter(outFilenameFullPath);
				IOUtils.writeToBufferedWriter(out, headerStr, true);

				// First, get the copy number ratios at a per-chromosome level				
				int[]   numSitesPerChrom           = new   int[Chrom.values().length];
				float[] avgReadCountPerChromNormal = new float[Chrom.values().length];
				float[] copyNumRatioPerChromNormal = new float[Chrom.values().length];
				PrimitiveWrapper.WFloat coverageRatioTumorToNormal = new PrimitiveWrapper.WFloat(0);
				float[] tumorNormalCopyNumRatiosPerChrom = calcAvgCoverageRatioPerChrom(allVariantRows, copyNumRatioPerChromNormal, numSitesPerChrom, avgReadCountPerChromNormal, coverageRatioTumorToNormal);
				
				// Now get the copy number ratios at a segmented sub-chromosomal level (by gene)
				float[] tumorCopyNumRatiosPerGene = calcRoughCopyNumberRatioPerSite(allVariantRows, coverageRatioTumorToNormal.mFloat);				
				
				// Now get the clusters
				float[] adjustedVAFNormal = new float[allVariantRows.size()];
				float[] adjustedVAFTumor  = new float[allVariantRows.size()];
				calculateAdjustedVAFs(allVariantRows, allelicBiasTable, adjustedVAFNormal, adjustedVAFTumor, copyNumRatioPerChromNormal, platform);

				double[] imbalancePValuesTumor  = new double[allVariantRows.size()];
				double[] imbalancePValuesNormal = new double[allVariantRows.size()];				
				getPValuesImbalance(allVariantRows, imbalancePValuesTumor, imbalancePValuesNormal, adjustedVAFNormal, adjustedVAFTumor);
				
				ClusterType[] clusters = null;
				System.out.println("\tInferred Copy Number and Allelic Imbalance Per Site...");
				
				if (UsePValuePlane) {
					clusters = Clustering.getClusters_withPlane(allVariantRows, imbalancePValuesTumor, imbalancePValuesNormal, tumorCopyNumRatiosPerGene, copyNumRatioPerChromNormal, adjustedVAFNormal, adjustedVAFTumor, outFilenameFullPath, indexFirstSomaticRowInAllVariants, platform); //get cluster assignments (HET ball, LOH sidelobes, DUP wedge, &c.)
				} else {
					clusters = Clustering.getClusters_Old(allVariantRows, imbalancePValuesTumor, tumorCopyNumRatiosPerGene, outFilenameFullPath, indexFirstSomaticRowInAllVariants, platform); //get cluster assignments (HET ball, LOH sidelobes, DUP wedge, &c.)
				}

				boolean toPrintWithCopyNum = true;
				if (toPrintWithCopyNum) {
					PrintStream outStream = IOUtils.getPrintStream(outFilenameFullPath + ".withCopyNum.txt");
					for (int row = 0; row < allVariantRows.size(); row++) {
						outStream.print(allVariantRows.get(row));					
						outStream.printf("\t%s", clusters[row].name());
						outStream.printf("\t%g\t%g", tumorCopyNumRatiosPerGene[row], (tumorCopyNumRatiosPerGene[row] * Script.DefaultDiploidCopyNumber));  
						outStream.printf("\t%g\t%g", adjustedVAFTumor[row], adjustedVAFNormal[row]);
						outStream.printf("\t%g\t%g", imbalancePValuesTumor[row], imbalancePValuesNormal[row]);
						outStream.println("");
					}
					IOUtils.closePrintStream(outStream);
				}


				System.out.println("\tGot clusters");

				// Now initialize the data structure needed to plot
				int[] clusterTypeCounts = Utils.getClusterTypeCounts(clusters);
				double[][][] clusterCoordinates  = getCoordinateArraysPerClusterType(clusterTypeCounts);
				double[][][] waterfallPlotTumor  = getCoordinateArraysPerClusterType(clusterTypeCounts);
				double[][][] waterfallPlotNormal = getCoordinateArraysPerClusterType(clusterTypeCounts);
				double[][][] copyNumPlot         = getCoordinateArraysPerClusterType(clusterTypeCounts);
				
				int[] clusterTypeIndex = new int[clusterTypeCounts.length];
				Arrays.fill(clusterTypeIndex, -1);  

				// Do the post-processing
				Chrom prevChrom = Chrom.c0;
				int prevPosition = 0;
				long positionGenomeWide = 0;
				
				
				DoubleArrayList chromBoundaryXValue = new DoubleArrayList();
				DoubleArrayList chromBoundaryYValue = new DoubleArrayList();
				
				int startingRowGermlineOrSomaticOrAll = 0;  // 0 because header line has been stripped away
				for (int i = startingRowGermlineOrSomaticOrAll; i < allVariantRows.size(); i++) {
					String strRow = allVariantRows.get(i);
					ClusterType clusterType = clusters[i];
					int indexForClusterType = ++(clusterTypeIndex[clusterType.ordinal()]);  // increase the index of the data structure of the cluster type

					String[] germCols = strRow.split(Utils.TabStr);

					//=IF(A3=A2,B3-B2+L2,B3+L2)
							
					String chromStr = (platform == SeqPlatform.Illumina) ? germCols[Script.Col_NAFTAFInput_Chrom] : Script.ChromPrefix + germCols[Script.Col_NAFTAFInput_Chrom].replace(Script.ChromPrefix, "");
					Chrom chrom = Chrom.getChrom(chromStr);
					int position = Integer.parseInt(germCols[Script.Col_NAFTAFInput_Position]);
		
					if (chrom == prevChrom) {
						positionGenomeWide += (position - prevPosition);
					} else {
						positionGenomeWide += position;
						for (double d = 0; d <= 5.0; d += 0.02) {
							chromBoundaryXValue.add(positionGenomeWide);
							chromBoundaryYValue.add(d);
						}
					}
					prevPosition = position;
					prevChrom = chrom;
					
					float vafNormal = adjustedVAFNormal[i]; //Clustering.extractVAFNormal(germCols, platform);
					float vafTumor  = adjustedVAFTumor[i]; //Clustering.extractVAFTumor (germCols, platform);
					clusterCoordinates[ clusterType.ordinal() ][0][indexForClusterType] = vafTumor;
					clusterCoordinates[ clusterType.ordinal() ][1][indexForClusterType] = vafNormal;
					
					waterfallPlotTumor[ clusterType.ordinal() ][0][indexForClusterType] = positionGenomeWide;
					waterfallPlotTumor[ clusterType.ordinal() ][1][indexForClusterType] = vafTumor;
					
					waterfallPlotNormal[ clusterType.ordinal() ][0][indexForClusterType] = positionGenomeWide;
					waterfallPlotNormal[ clusterType.ordinal() ][1][indexForClusterType] = vafNormal;	
					
					copyNumPlot[ clusterType.ordinal() ][0][indexForClusterType] = positionGenomeWide;
					copyNumPlot[ clusterType.ordinal() ][1][indexForClusterType] = (tumorCopyNumRatiosPerGene[i] * Script.DefaultDiploidCopyNumber);	
					

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

				DefaultXYDataset xyDatasetVAFPlot = new DefaultXYDataset();
				DefaultXYDataset xyDatasetWaterfallPlotTumor  = new DefaultXYDataset();
				DefaultXYDataset xyDatasetWaterfallPlotNormal = new DefaultXYDataset();							
				DefaultXYDataset xyDatasetCopyNumber = new DefaultXYDataset();
				
				for (ClusterType ct : ClusterType.values()) {					
					xyDatasetVAFPlot.addSeries(ct.name(), clusterCoordinates[ct.ordinal()]);
					xyDatasetWaterfallPlotTumor.addSeries( ct.name(), waterfallPlotTumor[ct.ordinal()]);
					xyDatasetWaterfallPlotNormal.addSeries(ct.name(), waterfallPlotNormal[ct.ordinal()]);
					xyDatasetCopyNumber.addSeries(ct.name(), copyNumPlot[ct.ordinal()]);
				}
				double[][] boundaryArrays = combineTwoDynamicArraysIntoOneStatic(chromBoundaryXValue, chromBoundaryYValue);
				xyDatasetWaterfallPlotTumor.addSeries("Boundary", boundaryArrays);
				xyDatasetWaterfallPlotNormal.addSeries("Boundary", boundaryArrays);
				xyDatasetCopyNumber.addSeries("Boundary", boundaryArrays);
				
				Clustering.plotVAFComparison(xyDatasetVAFPlot,            vafComparisonPlotDir + File.separator + samplenameRoot + ".VAFComparison", samplenameRoot);
				Clustering.plotVAFGenomeWide(xyDatasetWaterfallPlotTumor,  vafWaterfallPlotDir + File.separator + samplenameRoot + ".VAF_GenomeWide_Tumor",  samplenameRoot, true);				
				Clustering.plotVAFGenomeWide(xyDatasetWaterfallPlotNormal, vafWaterfallPlotDir + File.separator + samplenameRoot + ".VAF_GenomeWide_Normal", samplenameRoot, false);
				Clustering.plotCopyNumGenomeWide(xyDatasetCopyNumber,        copyNumberPlotDir + File.separator + samplenameRoot + ".CopyNumber_GenomeWide", samplenameRoot);
			}
		}		
	}
	
	private static double[][] combineTwoDynamicArraysIntoOneStatic(DoubleArrayList list1, DoubleArrayList list2) {
		return new double[][] { list1.toArray(), list2.toArray() };
	}
	
	private static double[][][] getCoordinateArraysPerClusterType(int[] clusterTypeCounts) {
		double[][][] clusterCoordinates = new double[clusterTypeCounts.length][2][];
		
		for (int clusterTypeIndex = 0; clusterTypeIndex < clusterCoordinates.length; clusterTypeIndex++) {
			clusterCoordinates[clusterTypeIndex][0] = new double[ clusterTypeCounts[clusterTypeIndex] ];
			clusterCoordinates[clusterTypeIndex][1] = new double[ clusterTypeCounts[clusterTypeIndex] ];
		}
		
		return clusterCoordinates;
	}

	/** Plots the VAF (variant allele frequency) of a tissue vs the genomic position. */
	public static void plotCopyNumGenomeWide(XYDataset xyDataset, String outFilenameRoot, String sampleName) {
		String xAxisLabel = "Position";
		String yAxisLabel = "Copy Number";
		String title = "Copy Number: " + sampleName;				
		
		JFreeChart theChart = 
				ChartFactory.createScatterPlot(title, xAxisLabel, yAxisLabel, xyDataset, PlotOrientation.VERTICAL, true, false, false);
		XYPlot xyPlot = (XYPlot) theChart.getPlot();
		
		XYItemRenderer itemRenderer = Clustering.getXYItemRendererHelper(5);		
		setSeriesPaintPerCluster(itemRenderer);
		xyPlot.setRenderer(itemRenderer);		
		
		xyPlot.setBackgroundPaint(ColorPastel.White.getColor());
		xyPlot.setDomainGridlinePaint(Color.white);
		xyPlot.setRangeGridlinePaint(Color.white);		
		
		LegendTitle legendTitle = theChart.getLegend();
		legendTitle.setID("Clusters");
		legendTitle.setItemFont(new Font("Arial", Font.BOLD, 20));

		//xyPlot.getRangeAxis().setRange(0, 1.02);		
		
		Font rangeAxisLabelFont = new Font("Arial", Font.BOLD, 20);
		Font rangeAxisTickFont = new Font("Arial", Font.BOLD, 20);
		xyPlot.getRangeAxis().setLabelFont(rangeAxisLabelFont);		
		xyPlot.getRangeAxis().setTickLabelFont(rangeAxisTickFont);
		xyPlot.getDomainAxis().setLabelFont(rangeAxisLabelFont);
		xyPlot.getDomainAxis().setTickLabelFont(rangeAxisTickFont);	
		GraphUtils.saveChartAsPNG(outFilenameRoot, theChart, 24000, 800);
	}
	
	/** Plots the VAF (variant allele frequency) of a tissue vs the genomic position. */
	public static void plotVAFGenomeWide(XYDataset xyDataset, String outFilenameRoot, String sampleName, boolean isTumor) {
		String xAxisLabel = "Position";
		String yAxisLabel = "VAF";
		String title = "VAF GenomeWide: " + sampleName + (isTumor ? " [Tumor]" : " [Normal]");
		
		JFreeChart theChart = 
				ChartFactory.createScatterPlot(title, xAxisLabel, yAxisLabel, xyDataset, PlotOrientation.VERTICAL, true, false, false);
		XYPlot xyPlot = (XYPlot) theChart.getPlot();
		
		XYItemRenderer itemRenderer = Clustering.getXYItemRendererHelper(5);		
		setSeriesPaintPerCluster(itemRenderer);
		xyPlot.setRenderer(itemRenderer);		
		
		xyPlot.setBackgroundPaint(ColorPastel.White.getColor());
		xyPlot.setDomainGridlinePaint(Color.white);
		xyPlot.setRangeGridlinePaint(Color.white);		
		
		LegendTitle legendTitle = theChart.getLegend();
		legendTitle.setID("Clusters");
		legendTitle.setItemFont(new Font("Arial", Font.BOLD, 20));

		xyPlot.getRangeAxis().setRange(0, 1.02);		
		
		Font rangeAxisLabelFont = new Font("Arial", Font.BOLD, 20);
		Font rangeAxisTickFont = new Font("Arial", Font.BOLD, 20);
		xyPlot.getRangeAxis().setLabelFont(rangeAxisLabelFont);		
		xyPlot.getRangeAxis().setTickLabelFont(rangeAxisTickFont);
		xyPlot.getDomainAxis().setLabelFont(rangeAxisLabelFont);
		xyPlot.getDomainAxis().setTickLabelFont(rangeAxisTickFont);	
		GraphUtils.saveChartAsPNG(outFilenameRoot, theChart, 2400, 800);
	}
	
	/** Plots the VAF (variant allele frequency) of the normal tissue comapred to the tumor tissue. */
	public static void plotVAFComparison(XYDataset xyDataset, String outFilenameRoot, String sampleName) {
		String xAxisLabel = "VAF Tumor";
		String yAxisLabel = "VAF Normal";
		String title = "VAF Comparison: " + sampleName;
		
		JFreeChart theChart = 
				ChartFactory.createScatterPlot(title, xAxisLabel, yAxisLabel, xyDataset, PlotOrientation.VERTICAL, true, false, false);
		XYPlot xyPlot = (XYPlot) theChart.getPlot();
		
		XYItemRenderer itemRenderer = Clustering.getXYItemRendererHelper(5);		
		setSeriesPaintPerCluster(itemRenderer);
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
		int width  = 900;
		int height = 900;
		GraphUtils.saveChartAsPNG(outFilenameRoot, theChart, width, height);
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
	
	private static void setSeriesPaintPerCluster(XYItemRenderer itemRenderer) {
		itemRenderer.setSeriesPaint(ClusterType.GainGermline.ordinal(), ColorPastel.Violet.getColor());
		itemRenderer.setSeriesPaint(ClusterType.GainSomatic.ordinal(), ColorPastel.Dark_Red.getColor());
		itemRenderer.setSeriesPaint(ClusterType.LOH.ordinal(), ColorPastel.Yellow_Green.getColor());
		itemRenderer.setSeriesPaint(ClusterType.cnLOH.ordinal(), ColorPastel.CMYK_Cyan.getColor());
		itemRenderer.setSeriesPaint(ClusterType.HETGermline.ordinal(), ColorPastel.Gray_50.getColor());
		itemRenderer.setSeriesPaint(ClusterType.HETSomatic.ordinal(), ColorPastel.Red_Orange.getColor());
		itemRenderer.setSeriesPaint(ClusterType.Noise.ordinal(), ColorPastel.Gray_60.getColor());
		itemRenderer.setSeriesPaint(ClusterType.Null.ordinal(), ColorPastel.Gray_30.getColor());	
		itemRenderer.setSeriesPaint(ClusterType.Null.ordinal() + 1, ColorPastel.Black.getColor());
	}

	// ========================================================================
	private static void calculateAdjustedVAFs(ArrayList<String> rows, AllelicBiasTable allelicBiasTable, 
											  float[] adjustedVAFNormal, float[] adjustedVAFTumor,
											  float[] copyNumRatioPerChromNormal,
											  SeqPlatform platform) {
		
		final float defaultVAFNormal = 0.50f;
		String delim = Utils.FileExtensionTSV.mDelimiter;
		
		for (int row = 0; row < rows.size(); row++) {
			String line = rows.get(row);
			
			Chrom chrom =      Chrom.getChrom( Utils.extractNthColumnValue(line, Script.Col_NAFTAFInput_Chrom,    delim) );
			int position =   Integer.parseInt( Utils.extractNthColumnValue(line, Script.Col_NAFTAFInput_Position, delim) );
			float vafNormal = extractVAFNormal(line, platform);
			float vafTumor  = extractVAFTumor(line, platform);
			float vafNormalExpected = defaultVAFNormal;
			
			// Set Default values
			float offset = 0.0f;  // Default
			float adjustmentFactor = 1.0f;
			
			// Perform adjustment calculations
			if (correctAllelicBias()) {
				if (Utils.inRangeLowerExclusive(vafNormal, AlleleFrequencyStatsForSample.VAFNormalFrameLower, AlleleFrequencyStatsForSample.VAFNormalFrameUpper)) {
					boolean isGermlineChromGain = copyNumRatioPerChromNormal[chrom.ordinal()] > ParamsNum.GermlineTrisomyThreshold.getValue().doubleValue();
					float avgVAFNormal = allelicBiasTable.getAvgVAF(chrom, position);
					if (avgVAFNormal > 0) {
						// Site exists in table
						if (isGermlineChromGain) {
							if (Math.abs(vafNormal - ExpectedVAFNormalTrisomy) < Math.abs(vafNormal - (1.0 - ExpectedVAFNormalTrisomy))) {
								vafNormalExpected = ExpectedVAFNormalTrisomy;
							} else {
								vafNormalExpected = 1.0f - ExpectedVAFNormalTrisomy;
							}
						} 

						adjustmentFactor = vafNormalExpected / avgVAFNormal;
						float absDiff = Math.abs(vafNormalExpected - avgVAFNormal);
						offset = (vafNormal > vafNormalExpected) ? -absDiff : ((vafNormal == vafNormalExpected) ? 0 : absDiff);												
					}
				}
			}
						
			// Implement adjustments
			if (UseBidrectionalAdditiveOffset) {
				adjustedVAFNormal[row] = vafNormal + offset;							
				adjustedVAFTumor[row]  = vafTumor  + offset;				
			} else {
				adjustedVAFNormal[row] = vafNormal * adjustmentFactor;			
				adjustedVAFTumor[row]  = vafTumor  * adjustmentFactor;				
			}
			
			// Ensure within bounds
			adjustedVAFNormal[row] = Math.min(1.0f, Math.max(0, adjustedVAFNormal[row]));							
			adjustedVAFTumor[row]  = Math.min(1.0f, Math.max(0, adjustedVAFTumor[row]));
		}		
	}
	// ========================================================================
	public static void getPValuesImbalance(ArrayList<String> rows, 
			double[] imbalancePValuesTumor, double[] imbalancePValuesNormal,
			float[] adjustedVAFNormal, float[] adjustedVAFTumor) {
		
		for (int row = 0; row < rows.size(); row++) {
			String line = rows.get(row);
									
			int totCovgTumor = Integer.parseInt(Utils.extractNthColumnValue(line, Script.Col_NAFTAFInput_TotalCoverageTumor,   Utils.FileExtensionTSV.mDelimiter));
			int varCovgTumor = (int) Math.round(totCovgTumor * adjustedVAFTumor[row]);

			//int varCovgNormal = Integer.parseInt(Utils.extractNthColumnValue(line, Script.Col_NAFTAFInput_VariantCoverageNormal, Utils.FileExtensionTSV.mDelimiter));
			int totCovgNormal = Integer.parseInt(Utils.extractNthColumnValue(line, Script.Col_NAFTAFInput_TotalCoverageNormal,   Utils.FileExtensionTSV.mDelimiter));
			int varCovgNormal = (int) Math.round(totCovgNormal * adjustedVAFNormal[row]);
			
			imbalancePValuesTumor[row]  = getPValuesImbalanceTissue(totCovgTumor, varCovgTumor);
			imbalancePValuesNormal[row] = getPValuesImbalanceTissue(totCovgNormal, varCovgNormal);		
		}
	}

	// ------------------------------------------------------------------------
	private static double getPValuesImbalanceTissue(int coverageTotal, int coverageVariant) {
		int maxRefOrVarCovg = Math.max(coverageVariant, coverageTotal - coverageVariant);
		return NumberUtils.cumulativeProbabilitySuccess(coverageTotal, maxRefOrVarCovg, 0.5);
	}

	// ========================================================================
	public static float[] calcRoughCopyNumberRatioPerSite(ArrayList<String> rows, float coverageRatioTumorNormalGenomeWide) {
		float[] copyNumRatios    =   new float[rows.size()];
		boolean[] isSomaticSite = new boolean[rows.size()];
		
		String prevGene = "";
		int rowOfFirstInstanceOfGene = -1;
		
		float ratioSum = 0;
		int readCountSumTumor = 0;
		int readCountSumNormal = 0;
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
					readCountSumTumor  += covgTumor;
					readCountSumNormal += covgNormal;
				} else {
					// We have a somatic site, we take no action in this block
				}
			} else {
				// New gene listed.  We must now write the ratios for the previous gene
				if (row > 0) {					
					float averageRatio = Script.TumorNormalRatioOfSomaticSite;
					
					if (numRowsWithSameGene == 0) {
						Utils.ensureTrue(false, "ERROR: Gene must > 0 sites representative of that gene!");
					} else if (numRowsWithSameGene == 1 && !AllowGeneWithOneSite) {					
						averageRatio = Script.DefaultTumorNormalRatio;
					} else if (numRowsWithSameGene > 1 || AllowGeneWithOneSite) {
						//averageRatio = ((ratioSum / (float) numRowsWithSameGene) / coverageRatioTumorNormalGenomeWide);
						averageRatio = ((float) readCountSumTumor / (float) readCountSumNormal) / coverageRatioTumorNormalGenomeWide;
					}
	
					for (int i = rowOfFirstInstanceOfGene; i <= row; i++) {					
						copyNumRatios[i] = isSomaticSite[i] ? Script.TumorNormalRatioOfSomaticSite : averageRatio;
					}
				}
				
				// New gene listed.  Thus, set a new row of first gene			
				rowOfFirstInstanceOfGene = row;
				prevGene = currentGene;
				if (isSomaticSite[row]) {
					ratioSum = numRowsWithSameGene = readCountSumTumor = readCountSumNormal = 0;
				} else {
					numRowsWithSameGene = 1;
					ratioSum = ratio;					 					
					readCountSumTumor  = covgTumor;
					readCountSumNormal = covgNormal;
				}
			}			
		}
		
		// Normalize to diploid
		DescriptiveStatistics ds = new DescriptiveStatistics();
		for (int row = 0; row < rows.size(); row++) {
			ds.addValue(copyNumRatios[row]);
		}
		double median = ds.getPercentile(50);
		System.out.println("Median: " + median);
		double diff = 1.0 - median;
		
		for (int row = 0; row < rows.size(); row++) {
			copyNumRatios[row] += (float) diff;
		}
		
		return copyNumRatios;
	}

	// ========================================================================
	private static final int roundToNearest5(int num) {
		return (int) (Math.round(num / 5.0) * 5);
	}
	
	// ========================================================================
	/** Tumor : Normal avg coverage ratio per chromosome. */
	public static float[] calcAvgCoverageRatioPerChrom(ArrayList<String> rows, float[] copyNumRatioPerChromNormal, int[] numSitesPerChrom, float[] avgReadCountPerChromNormal, PrimitiveWrapper.WFloat coverageRatioTumorToNormal) {			
		float[] tumorNormalRatioPerChrom = new float[Chrom.values().length];
		int[] totalReadCountPerChromNormal = new int[Chrom.values().length];
		int[] totalReadCountPerChromTumor  = new int[Chrom.values().length];
		
		Arrays.fill(tumorNormalRatioPerChrom, Script.DefaultTumorNormalRatio);  // Fill with 1.0 since that's an equal ratio
		Arrays.fill(copyNumRatioPerChromNormal, Script.DefaultNormalRatio);
		
		Arrays.fill(numSitesPerChrom, 0);  // initialize counts
		Arrays.fill(avgReadCountPerChromNormal, 0); // initialize
		Arrays.fill(totalReadCountPerChromNormal, 0);
		Arrays.fill(totalReadCountPerChromTumor,  0);
		
		DynamicBucketCounter readCountTalliesTumor  = new DynamicBucketCounter();
		DynamicBucketCounter readCountTalliesNormal = new DynamicBucketCounter();
		
		for (int row = 0; row < rows.size(); row++) {
			String line = rows.get(row);
			
			Chrom chrom    =   Chrom.getChrom( Utils.extractNthColumnValue(line, Script.Col_NAFTAFInput_Chrom,               Utils.FileExtensionTSV.mDelimiter) );
			int covgNormal = Integer.parseInt( Utils.extractNthColumnValue(line, Script.Col_NAFTAFInput_TotalCoverageNormal, Utils.FileExtensionTSV.mDelimiter) );
			int covgTumor  = Integer.parseInt( Utils.extractNthColumnValue(line, Script.Col_NAFTAFInput_TotalCoverageTumor,  Utils.FileExtensionTSV.mDelimiter) );

			int chromIndex = chrom.ordinal();
			totalReadCountPerChromNormal[chromIndex] += covgNormal;
			totalReadCountPerChromTumor[chromIndex]  += covgTumor;
			++numSitesPerChrom[chromIndex];
			
			//=ROUND(A1*2/10,0)*10/2
					
			readCountTalliesTumor.incrementCount(roundToNearest5(covgTumor));
			readCountTalliesNormal.incrementCount(roundToNearest5(covgNormal));			
		}
		
		// Calculate the genome-wide average read count in the normal
		int totalReadCountNormal = ArrayUtils.arraySum(totalReadCountPerChromNormal);
		int totalReadCountTumor  = ArrayUtils.arraySum(totalReadCountPerChromTumor);
		int totalNumSites        = ArrayUtils.arraySum(numSitesPerChrom);
		float avgCoverageNormal = (float) totalReadCountNormal / (float) totalNumSites;
		float avgCoverageTumor  = (float) totalReadCountTumor  / (float) totalNumSites;
		//coverageRatioTumorToNormal.mFloat = avgCoverageTumor / avgCoverageNormal;
		coverageRatioTumorToNormal.mFloat = (float) readCountTalliesTumor.getKeyWithMaxCount() / (float) readCountTalliesNormal.getKeyWithMaxCount(); 
		System.out.println("Tumor Normal Genome-wide ratio: " + coverageRatioTumorToNormal.mFloat);
				
		System.out.printf("Average Read Count Normal: %g\n", avgCoverageNormal);
		
		for (Chrom chrom : Chrom.values()) {
			int chromIndex = chrom.ordinal();
			
			if (totalReadCountPerChromNormal[chromIndex] > 0) {
				tumorNormalRatioPerChrom[chromIndex] = (float) totalReadCountPerChromTumor[chromIndex] / (float) totalReadCountPerChromNormal[chromIndex];
			}
			
			if (numSitesPerChrom[chromIndex] > 0) {
				avgReadCountPerChromNormal[chromIndex] = totalReadCountPerChromNormal[chromIndex] / numSitesPerChrom[chromIndex];
				copyNumRatioPerChromNormal[chromIndex] = avgReadCountPerChromNormal[chromIndex] / avgCoverageNormal;
			}	
			
			System.out.printf("Chrom: %d\tNormal Ratio:%g\tTumor-Normal Ratio %g\n", chromIndex, copyNumRatioPerChromNormal[chromIndex], tumorNormalRatioPerChrom[chromIndex]);
		}
					
		return tumorNormalRatioPerChrom;
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
		clusterTypeIDsFromAlgorithm[ClusterType.GainSomatic.ordinal()]  = -1;
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
					clusterAssignments[i] = clusterTypeIDsFromAlgorithm[ClusterType.GainSomatic.ordinal()];   // only change if we're in a het ball region
				}
			}

		}

		// Now assign the cluster types

		return Clustering.assignClusters(rows, platform, vafNormalFrameAdjustedLower, vafNormalFrameAdjustedUpper, clusterAssignments, clusterTypeIDsFromAlgorithm);
	}
	
	
	// ========================================================================
	public static ClusterType[] getClusters_withPlane(ArrayList<String> rows, 
			double[] imbalancePValuesTumor, double[] imbalancePValuesNormal, 
			float[] copyNumRatios,  float[] copyNumRatioPerChromNormal,
			float[] adjustedVAFNormal, float[] adjustedVAFTumor, 
			String outFilenameFullPath, int startingRowGermlineOrSomaticOrAll, SeqPlatform platform) {

		// Get the allele frequency statistics, and adjust the frames based on the resulting standard deviation
		AlleleFrequencyStatsForSample afStatsSample = new AlleleFrequencyStatsForSample();
		afStatsSample.tabulateAndPerformStatistics(rows, platform);
		float vafNormalFrameAdjustedLower = afStatsSample.getValueNStandardDeviationsAway(-NAF_STRIP_EXPANDER); 								
		float vafNormalFrameAdjustedUpper = afStatsSample.getValueNStandardDeviationsAway( NAF_STRIP_EXPANDER);
		
		if (correctAllelicBias()) {
			vafNormalFrameAdjustedLower = AlleleFrequencyStatsForSample.VAFNormalFrameLower; 								
			vafNormalFrameAdjustedUpper = AlleleFrequencyStatsForSample.VAFNormalFrameUpper;
		}

		//apply DBScan to points within NAF frame
		ArrayList<Floint> pointsUpperPlane = new ArrayList<Floint>();
		ArrayList<Floint> pointsLowerPlane = new ArrayList<Floint>();
		int[] indexMapToUpperLowerPlane    = new int[rows.size()];
		IntArrayList mapToRowsFromUpper = new IntArrayList(rows.size());
		IntArrayList mapToRowsFromLower = new IntArrayList(rows.size());

		getValidPointsForClustering(rows, pointsUpperPlane, pointsLowerPlane, indexMapToUpperLowerPlane, 
				                    mapToRowsFromUpper, mapToRowsFromLower, 
				                    imbalancePValuesTumor, imbalancePValuesNormal, 
				                    ScalingFactor, 
				                    copyNumRatios, copyNumRatioPerChromNormal,
				                    adjustedVAFNormal, adjustedVAFTumor, 
				                    platform, 
				                    vafNormalFrameAdjustedLower, vafNormalFrameAdjustedUpper);		

		System.out.printf("\tPoints Lower Plane: %d\n", pointsLowerPlane.size());
		System.out.printf("\tPoints Upper Plane: %d\n", pointsUpperPlane.size());
		System.out.println("\tBegin clustering algorithm: " + (new Date()).toString());

		//DBScanFaster dbscanner = new DBScanFaster(points, HET_BALL_EPS + DUP_WEDGE_LASSO, DUP_WEDGE_MINPTS, 0, 0, 1, 1); //parameters for capturing HET ball//KDBSCAN(param, 10, 0.0325f);
			
		int[] clusterTypeIDsFromAlgorithm = new int[ClusterType.values().length];
		clusterTypeIDsFromAlgorithm[ClusterType.GainGermline.ordinal()]   = -ClusterType.GainGermline.ordinal() - 1;
		clusterTypeIDsFromAlgorithm[ClusterType.GainSomatic.ordinal()]    = -ClusterType.GainSomatic.ordinal()  - 1;
		clusterTypeIDsFromAlgorithm[ClusterType.cnLOH.ordinal()]          = -ClusterType.cnLOH.ordinal()        - 1;
		clusterTypeIDsFromAlgorithm[ClusterType.Null.ordinal()]           = -ClusterType.Null.ordinal()         - 1;		
		clusterTypeIDsFromAlgorithm[ClusterType.HETGermline.ordinal()]    = -ClusterType.HETGermline.ordinal()  - 1;
		clusterTypeIDsFromAlgorithm[ClusterType.Noise.ordinal()] = DBSCAN2.getClusterIDOfNoise();
		
		// ----------------------- Scan the lower plane
		DBScanFaster dbscannerLowerPlane = new DBScanFaster(pointsLowerPlane, Clustering.HET_BALL_EPS, Clustering.HET_BALL_MINPTS, 0, 0, 1, 1); //parameters for capturing HET ball//KDBSCAN(param, 10, 0.0325f);
		dbscannerLowerPlane.cluster();

		// Now, we find the central heterozygous cluster ID and replace it with a fixed, known number.   
		// In parallel, we track which points in the lower plane were not part of the central cluster.
		// These points are clustered in another pass with lower neighbor density/size parameters.
		ArrayList<Floint> nonHetPoints = new ArrayList<Floint>(pointsLowerPlane.size());
		boolean[] isNonHetPoint                  = new boolean[pointsLowerPlane.size()];		
		int clusterIDofHetBall = dbscannerLowerPlane.getCentralClusterID();
		int[] clusterAssignmentsLowerPlane = dbscannerLowerPlane.getClustAssignments();  // save and cache	
		
		for (int i = 0; i < clusterAssignmentsLowerPlane.length; i++) {
			isNonHetPoint[i] = (clusterAssignmentsLowerPlane[i] != clusterIDofHetBall); 
			if (isNonHetPoint[i]) {
				//nonHetPoints.add(pointsLowerPlane.get(i));
				//clusterAssignmentsLowerPlane[i] = clusterTypeIDsFromAlgorithm[ClusterType.Null.ordinal()];
			} else {
				clusterAssignmentsLowerPlane[i] = clusterTypeIDsFromAlgorithm[ClusterType.HETGermline.ordinal()];
			}
		}
		
		System.out.println("\tPoints Lower Plane Non-Het: " + nonHetPoints.size());
		
		//DBScanFaster dbscannerLowerPlaneNonHet = new DBScanFaster(nonHetPoints, Clustering.NON_HET_BALL_EPS, Clustering.NON_HET_BALL_MINPTS, 0, 0, 1, 1);
		//dbscannerLowerPlaneNonHet.cluster();
		//int[] clusterAssignmentsNonHet = dbscannerLowerPlaneNonHet.getClustAssignments();

		// -------------------- UPPER PLANE
		DBScanFaster dbscannerUpperPlane = new DBScanFaster(pointsUpperPlane, Clustering.NON_HET_BALL_EPS, Clustering.NON_HET_BALL_MINPTS, 0, 0, 1, 1);
		dbscannerUpperPlane.cluster();
		int[] clusterAssignmentsUpperPlane = dbscannerUpperPlane.getClustAssignments();
		
		for (int ind = 0; ind < clusterAssignmentsUpperPlane.length; ind++) {	
			int indexIntoOriginalRows = mapToRowsFromUpper.get(ind);
			float copyNum = copyNumRatios[indexIntoOriginalRows] * Script.DefaultDiploidCopyNumber;
			
			if (isCopyNumAmplified(copyNum) && clusterAssignmentsUpperPlane[ind] != DBSCAN2.ClusterIDOfNoise) {
				clusterAssignmentsUpperPlane[ind] = clusterTypeIDsFromAlgorithm[ClusterType.GainSomatic.ordinal()];
			} else if (isCopyNumInDiploidRange(copyNum) && clusterAssignmentsUpperPlane[ind] != DBSCAN2.ClusterIDOfNoise) {
				clusterAssignmentsUpperPlane[ind] = clusterTypeIDsFromAlgorithm[ClusterType.cnLOH.ordinal()];
			}
			
			String line = rows.get(indexIntoOriginalRows);
			Chrom chrom = Chrom.getChrom( Utils.extractNthColumnValue(line, Script.Col_NAFTAFInput_Chrom,    Utils.FileExtensionTSV.mDelimiter) );
			if (copyNumRatioPerChromNormal[chrom.ordinal()] > ParamsNum.GermlineTrisomyThreshold.getValue().doubleValue()) {
				clusterAssignmentsUpperPlane[ind] = clusterTypeIDsFromAlgorithm[ClusterType.GainGermline.ordinal()];
			} else {
				if (ForcePointsOnDiagonalAsNull) {
					if (pointOnDiagonal(pointsUpperPlane.get(ind), ClusterDiagonalLeeway) && isCopyNumInDiploidRange(copyNum)) {
						clusterAssignmentsUpperPlane[ind] = clusterTypeIDsFromAlgorithm[ClusterType.Null.ordinal()];
					}
				}
				
			}
		}
		
		System.out.println("\tEnd clustering algorithm: " + (new Date()).toString());

		// --- LOWER PLANE ADJUSTMENT
		/*
		int nonHetIndex = -1;
		for (int i = 0; i < clusterAssignmentsLowerPlane.length; i++) {
			if (isNonHetPoint[i]) {
				++nonHetIndex;
								
				float copyNum = copyNumRatios[mapToRowsFromLower.get(i)] * Script.DefaultDiploidCopyNumber;

				// Now adjust
				if (ForcePointsOnDiagonalAsNull && pointOnDiagonal(pointsLowerPlane.get(i), ClusterDiagonalLeeway)) {
					clusterAssignmentsLowerPlane[i] = clusterTypeIDsFromAlgorithm[ClusterType.Null.ordinal()];
				} else if (isCopyNumAmplified(copyNum) && clusterAssignmentsLowerPlane[i] != DBSCAN2.ClusterIDOfNoise) {
					clusterAssignmentsLowerPlane[i] = clusterTypeIDsFromAlgorithm[ClusterType.GainSomatic.ordinal()];
				} else if (isCopyNumInDiploidRange(copyNum) && clusterAssignmentsLowerPlane[i] != DBSCAN2.ClusterIDOfNoise) {
					clusterAssignmentsLowerPlane[i] = clusterTypeIDsFromAlgorithm[ClusterType.cnLOH.ordinal()];
				} else {
					clusterAssignmentsLowerPlane[i] = clusterAssignmentsNonHet[nonHetIndex];   
				}
			}
		}*/
		
		
		// Now merge the two planes into one array
		ArrayUtils.IntArray clusterAssignmentsFinal = new ArrayUtils.IntArray(pointsLowerPlane.size() + pointsUpperPlane.size());  		
		for (int row = 0; row < rows.size(); row++) {
			if (indexMapToUpperLowerPlane[row] != Integer.MAX_VALUE) {				
				if (indexMapToUpperLowerPlane[row] >= 0) {
					int trueIndex = indexMapToUpperLowerPlane[row];
					clusterAssignmentsFinal.add(clusterAssignmentsUpperPlane[trueIndex]);
				} else {
					int trueIndex = -(indexMapToUpperLowerPlane[row] + 1);
					clusterAssignmentsFinal.add(clusterAssignmentsLowerPlane[trueIndex]);
				}
			}
		}
		

		return Clustering.assignClusters(rows, platform, vafNormalFrameAdjustedLower, vafNormalFrameAdjustedUpper, clusterAssignmentsFinal.mArray, clusterTypeIDsFromAlgorithm);
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

					} else if (assignedClusterID == clusterTypeIDsFromAlgorithm[ClusterType.GainSomatic.ordinal()]) {
						returnClusters[row] = ClusterType.GainSomatic; //DUP

					} else if (assignedClusterID == clusterTypeIDsFromAlgorithm[ClusterType.GainGermline.ordinal()]) {
						returnClusters[row] = ClusterType.GainGermline; //DUP
						
					} else if (assignedClusterID == clusterTypeIDsFromAlgorithm[ClusterType.cnLOH.ordinal()]) {
						returnClusters[row] = ClusterType.cnLOH; //DUP
						
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
													double[] imbalancePValuesTumor,
													double[] imbalancePValuesNormal,
													double scalingFactor,
													float[] copyNumRatios, float[] copyNumRatioPerChromNormal,
													float[] adjustedVAFNormal, float[] adjustedVAFTumor, 
													SeqPlatform platform, 
													float vafBoundLower, float vafBoundUpper) {
		
		pointsUpperPlane.clear();
		pointsLowerPlane.clear();
		double alphaAdjusted = ParamsBool.IgnoreMultipleTesting.getValue() ? 
				 PValueBinDistAlpha_UpperPlaneThresh :
				(PValueBinDistAlpha_UpperPlaneThresh * 0.008 /* Light FDR */);
				
		for (int row = 0; row < rows.size(); row++) {
			String line = rows.get(row);			
			
			float vafNormal = Clustering.extractVAFNormal(line, platform);
			if (Utils.inRangeLowerExclusive(vafNormal, vafBoundLower, vafBoundUpper)) {				
				Floint thePoint = new Floint(adjustedVAFTumor[row], adjustedVAFNormal[row], (float) (imbalancePValuesTumor[row] * scalingFactor));
				
				// Assume p-value as vertical row factor
				boolean tumorSigImbalanced = (imbalancePValuesTumor[row] <= alphaAdjusted); // && !isCopyNumRatioInDiploidRange(copyNumRatios[row]);
				Chrom chrom = Chrom.getChrom( Utils.extractNthColumnValue(line, Script.Col_NAFTAFInput_Chrom,    Utils.FileExtensionTSV.mDelimiter) );
				boolean normalSigImbalanced = (imbalancePValuesNormal[row] <= alphaAdjusted) && (copyNumRatioPerChromNormal[chrom.ordinal()] > ParamsNum.GermlineTrisomyThreshold.getValue().doubleValue());
				
				if (tumorSigImbalanced || normalSigImbalanced) {
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
	private static boolean isCopyNumAmplified(float copyNum) {
		return (copyNum >= ParamsNum.AmplificationThreshold.getValue().doubleValue());
	}
	
	// ========================================================================
	private static boolean isCopyNumInDiploidRange(float copyNum) {
		return (ParamsNum.DeletionThreshold.getValue().doubleValue() < copyNum) && (copyNum < ParamsNum.AmplificationThreshold.getValue().doubleValue());		
	}
	
	// ========================================================================
	private static boolean isCopyNumRatioInDiploidRange(float copyNumRatio) {
		float copyNum = copyNumRatio * Script.DefaultDiploidCopyNumber;
		return isCopyNumInDiploidRange(copyNum);
	}
	
	// ========================================================================
	// Convenience function 
	private static ArrayList<Floint> 
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
