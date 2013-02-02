import genomeEnums.Chrom;
import genomeEnums.Nuc;
import genomeEnums.VariantLocation;
import genomeUtils.GenotypeUtils;
import genomeUtils.RegionSimulator;

import java.awt.Color;
import java.awt.Font;
import java.awt.geom.Ellipse2D;
//import java.awt.geom.Ellipse2D.Double;
import java.io.BufferedWriter;
import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.Date;

import javax.swing.text.html.HTMLDocument.HTMLReader.IsindexAction;

import lohcateEnums.ClusterType;
import lohcateEnums.ColorPastel;
import lohcateEnums.MutationType;
import lohcateEnums.SeqPlatform;

import nutils.ArgumentParserUtils;
import nutils.ArrayUtils;
import nutils.ArrayUtils.DoubleParallelArray;
import nutils.BitSetUtils;
import nutils.CompareUtils;
import nutils.ContingencyTable;
import nutils.GraphUtils;
import nutils.IOUtils;
import nutils.NullaryClassFactory;
import nutils.PrimitiveWrapper;
import nutils.BitSetUtils.BitShiftAndMask;
import nutils.ContingencyTable.ContingencyTableValue;
import nutils.StringUtils;
import nutils.StringUtils.FileExtensionAndDelimiter;
import nutils.counter.DynamicBucketCounter;

import org.apache.commons.math3.distribution.BinomialDistribution;
import org.apache.commons.math3.distribution.PoissonDistribution;
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

import nutils.NumberUtils;
import shared.Utils;
import sun.text.CompactByteArray;


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
		AmplificationThreshold(2.3f, 'a', "ampThresh", "Copy_Number_Threshold_for_Amplification"),
		DeletionThreshold     (1.7f, 'd', "delThresh", "Copy_Number_Threshold_for_Deletion"),
		GermlineTrisomyThreshold(1.48f, JSAP.NO_SHORTFLAG, "germlineAneuploidyGainThreshold", "Copy_Number_Threshold_for_Germline_Chromosomal_Gain"),
		FDRAlpha              (0.01f, JSAP.NO_SHORTFLAG, "FDR_Alpha", "FDR_Alpha_Value_for_Allelic_Fraction_Imbalance")
		;
		
		private Number mDefaultValue;
		private Number mValue;
		private char mFlagShort;
		private String mFlagLong;
		private String mUsageName;
		
		private ParamsNum(Number defaultValue, char flagShort, String flagLong, String usageName) {
			mValue = mDefaultValue = defaultValue;
			mFlagShort = flagShort;
			mFlagLong  = flagLong;
			mUsageName = usageName;
		}
		
		public Number getValue() { return mValue; } 
		public Number getDefault() { return mDefaultValue; }
		public char getShortFlag() { return mFlagShort; }
		public String getLongFlag() { return mFlagLong; }
		public void setValue(Number theValue) { mValue = theValue; }
		public String getUsageName() { return mUsageName; }		
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
		
		for (ParamsNum pn : ParamsNum.values()) {
			FlaggedOption fo = new FlaggedOption(pn.name())
					.setStringParser(JSAP.DOUBLE_PARSER)
					.setDefault(pn.getDefault().toString())
					.setShortFlag(pn.getShortFlag())
					.setLongFlag(pn.getLongFlag())
					.setUsageName(pn.getUsageName());
			ArgumentParserUtils.registerJSAPParameter(jsap, fo);
		}
				
		return jsap;
	}
	
	// ========================================================================
	public static void configureParameters(JSAPResult config) {
		for (ParamsNum pn : ParamsNum.values()) {
			pn.setValue(config.getDouble(pn.name()));
		}		
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
		
		public void tabulateAndPerformStatistics(ArrayList<ClusteringInputOneSite> rows, SeqPlatform platform) {			
			tallyVariantAlleleFrequenciesIntoBins(rows, platform, true);
			calculateSummaryStatistics();
		}
		
		public float getValueNStandardDeviationsAway(float n) { return (mCountMean + (n * mStdDev)); }
		
		private void deduceBinValues() {
			for (int i = 0; i < mBinValue.length; i++) {
				mBinValue[i] = VAFNormalFrameLower + ((i + 1) * BinSize);
			}
		}
		
		private void tallyVariantAlleleFrequenciesIntoBins(ArrayList<ClusteringInputOneSite> rows, SeqPlatform platform, boolean clearBinCountBins) {
			if (clearBinCountBins) {
				Arrays.fill(mBinCount, 0);
			}
			
			// First, tally the variant allele frequencies into bins
			for (int row = 0; row < rows.size(); row++) {
				float vafNormal = rows.get(row).calcVAFNormal(); 
				if (NumberUtils.inRangeLowerExclusive(vafNormal, VAFNormalFrameLower, VAFNormalFrameUpper)) {
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
		ClusteringInputOneSite.TestClusteringInputOneSite_Robust();
	}
	
	// ========================================================================
	private static void classifySitesHelper_MakeSubDirs(String outDir, String vafComparisonPlotDir, String vafWaterfallPlotDir, String copyNumberPlotDir) {
		IOUtils.createDirectoryPath(outDir, false);
		IOUtils.createDirectoryPath(vafComparisonPlotDir, false);
		IOUtils.createDirectoryPath(vafWaterfallPlotDir, false);
		IOUtils.createDirectoryPath(copyNumberPlotDir, false);
	}
	
	// ========================================================================
	private static ClusteringInputOneSample readLinesFromFiles(File file) {
		String somaticFilename = file.getAbsolutePath().replace(VariantLocation.Germline.toLowerCase(), VariantLocation.Somatic.toLowerCase());
		ArrayList<String> somaticSpecificVariantRows  = IOUtils.readAllLinesFromFile(somaticFilename);
		ArrayList<String> germlineSpecificVariantRows = IOUtils.readAllLinesFromFile(file.getAbsolutePath());
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
		ArrayList<String> allVariantRowsStr = new ArrayList<String>(germlineSpecificVariantRows.size() + somaticSpecificVariantRows.size());
		allVariantRowsStr.addAll(germlineSpecificVariantRows);
		int indexFirstSomaticRowInAllVariants = allVariantRowsStr.size();  // save the size with just the germline variants added
		allVariantRowsStr.addAll(somaticSpecificVariantRows);

		//Collections.sort(germlineSpecificVariantRows, LineComparatorTab);
		//Collections.sort(somaticSpecificVariantRows,  LineComparatorTab);
		Collections.sort(allVariantRowsStr,              Script.LineComparatorTab);
		System.out.println("\tAll Sites Finished Sorting...");
		
		// Now parse into numeric values
		ClusteringInputOneSample oneSampleData = new ClusteringInputOneSample(allVariantRowsStr);				
		allVariantRowsStr.clear();
		
		return oneSampleData;
	}

	// ========================================================================
	// ========================================================================
	private static class ClusteringInputOneSampleMetaData {

		float[]  mAdjustedVAFNormal;
		float[]  mAdjustedVAFTumor;
		double[] mImbalancePValuesTumor;
		double[] mImbalancePValuesNormal;		
		
		PrimitiveWrapper.WFloat mCoverageRatioTumorToNormal;
		DynamicBucketCounter mReadCountTalliesTumor;
		DynamicBucketCounter mReadCountTalliesNormal;
		
		int[]   mNumSitesPerChrom;
		float[] mAvgReadCountPerChromNormal;
		float[] mCopyNumRatioPerChromNormal;
		float[] mTumorCopyNumRatiosPerGene;
		
		double mFDRNormal = 0;
		double mFDRTumor = 0;
		
		public ClusteringInputOneSampleMetaData(int numSites) {
			mAdjustedVAFNormal         = new float[numSites];
			mAdjustedVAFTumor          = new float[numSites];
			mImbalancePValuesTumor    = new double[numSites];
			mImbalancePValuesNormal   = new double[numSites];
			mTumorCopyNumRatiosPerGene = new float[numSites];
			
			mCoverageRatioTumorToNormal = new PrimitiveWrapper.WFloat(0);
			mReadCountTalliesTumor  = new DynamicBucketCounter();
			mReadCountTalliesNormal = new DynamicBucketCounter();
			
			mNumSitesPerChrom           = new   int[Chrom.values().length];
			mAvgReadCountPerChromNormal = new float[Chrom.values().length];
			mCopyNumRatioPerChromNormal = new float[Chrom.values().length];
		}		
	}
	
	// ========================================================================
	// ========================================================================
	/** Performs pre-processing (copy number calling, allelic imbalace p-values) on the input data
	 *  and stores the results in the metadata
	 */
	private static void classifySitesHelper_preprocessMetaData(ClusteringInputOneSample oneSampleData, ClusteringInputOneSampleMetaData metaData, AllelicBiasTable allelicBiasTable, SeqPlatform platform) {
		
		// First, get the copy number ratios at a per-chromosome level				
		float[] tumorNormalCopyNumRatiosPerChrom = calcAvgCoverageRatioPerChrom(oneSampleData.mInfoSites, metaData.mCopyNumRatioPerChromNormal, metaData.mNumSitesPerChrom, metaData.mAvgReadCountPerChromNormal, metaData.mCoverageRatioTumorToNormal, metaData.mReadCountTalliesTumor, metaData.mReadCountTalliesNormal);
						
		// Now get the copy number ratios at a segmented sub-chromosomal level (by gene)
		calcRoughCopyNumberRatioPerSite(oneSampleData.mInfoSites, metaData.mTumorCopyNumRatiosPerGene, metaData.mCoverageRatioTumorToNormal.mFloat, metaData.mReadCountTalliesTumor, metaData.mReadCountTalliesNormal);				
		
		// Now adjust the VAF values basead on biases
		calculateAdjustedVAFs(oneSampleData.mInfoSites, allelicBiasTable, metaData.mAdjustedVAFNormal, metaData.mAdjustedVAFTumor, metaData.mCopyNumRatioPerChromNormal, platform);

		getPValuesImbalance(oneSampleData.mInfoSites, metaData.mImbalancePValuesTumor, metaData.mImbalancePValuesNormal, metaData.mAdjustedVAFNormal, metaData.mAdjustedVAFTumor);
						
		// Calculate the FDR for tumor and normal
		boolean pointsIndependent = true;
		metaData.mFDRTumor  = getFDR_BenjaminiHochberg(metaData.mImbalancePValuesTumor,  ParamsNum.FDRAlpha.getValue().doubleValue(), pointsIndependent);
		metaData.mFDRNormal = getFDR_BenjaminiHochberg(metaData.mImbalancePValuesNormal, ParamsNum.FDRAlpha.getValue().doubleValue(), pointsIndependent);
	}

	// ========================================================================
	private static Boolean classifySites_addToContingencyTables(ClusterType eventTest, ClusterType eventTruth, ArrayList<ContingencyTable> eventTables) {		
		
		if (CompareUtils.isNull(eventTest) || CompareUtils.isNull(eventTruth)) {
			return null;
		}
		
		if (eventTest == ClusterType.Null) {
			return Boolean.FALSE;
		}
		
		for (ClusterType ct : ClusterType.values()) {
			ContingencyTable table = eventTables.get(ct.ordinal());
			
			boolean matchTruthAndEvent = (eventTruth == ct);
			boolean matchTestAndEvent  = (eventTest  == ct);
			
			if (matchTruthAndEvent && matchTestAndEvent) {
				table.incrCount(ContingencyTableValue.TruePositive);
			} else if (matchTruthAndEvent && !matchTestAndEvent) {
				table.incrCount(ContingencyTableValue.FalseNegative);
			} else if (!matchTruthAndEvent && matchTestAndEvent) {
				table.incrCount(ContingencyTableValue.FalsePositive);
			} else {
				table.incrCount(ContingencyTableValue.TrueNegative);
			}			
		}
		
		return Boolean.TRUE;
	}
	
	// ========================================================================
	/**
	 * Curate SNP calls by clustering data points into HET ball, DUP wedge, LOH sidelobes, &c., and by grabbing dbsnp population allele frequencies when possible
	 * @param inDir naf-taf-inputs
	 * @param opt 0::Illumina, 1::SOLiD
	 */
	public static void classifySites(String inDir, String allelicBiasInFile, String outDir, String vafComparisonPlotDir, String vafWaterfallPlotDir, String copyNumberPlotDir, SeqPlatform platform) {
		File[] files = (new File(inDir)).listFiles();
		StringBuilder sb = new StringBuilder(8192);
		StringUtils.FileExtensionAndDelimiter fileExtDelim = StringUtils.FileExtensionTSV;		

		// TODO -- make column names static constants
		System.out.println("Reading Allelic Bias file...");
		AllelicBiasTable allelicBiasTable = correctAllelicBias() ? AllelicBiasTable.readFileAndConstructTable(allelicBiasInFile, 3, 4) : null;		
		System.out.println("Finished Reading Allelic Bias file...");
		
		// Create output directory
		classifySitesHelper_MakeSubDirs(outDir, vafComparisonPlotDir, vafWaterfallPlotDir, copyNumberPlotDir);

		String[] columnHeaders = new String[] { "chr", "pos", "VAF_Normal", "VAF_Tumor", "dbsnp", "gene", "mutation_type", "germ_som", "cluster" };
		String headerStr = StringUtils.constructColumnDelimitedString(columnHeaders, fileExtDelim.mDelimiter, sb, true).toString();		

		ArrayList<DynamicBucketCounter> dbcByChromLOH = NullaryClassFactory.newList(ArrayList.class, DynamicBucketCounter.ClassFactory, Chrom.values().length); 
		
		int countAtPositionMax = -1; 
		int fileIndex = 0;
		for (File file : files) {			
			int indexOfSubstring = file.getName().indexOf(Script.GermlineSuffix);
			if (indexOfSubstring >= 0) {
				String samplenameRoot = file.getName().substring(0, indexOfSubstring);  				
				System.out.println("Processing (" + ++fileIndex + "): " + file.getName());
				ClusteringInputOneSample oneSampleData = readLinesFromFiles(file);
				
				String outFilename = samplenameRoot + fileExtDelim.mExtension;
				String outFilenameFullPath = outDir + File.separator + outFilename; 
				BufferedWriter out = IOUtils.getBufferedWriter(outFilenameFullPath);
				IOUtils.writeToBufferedWriter(out, headerStr, true);
				
				ArrayList<ContingencyTable> eventTables = NullaryClassFactory.newList(ArrayList.class, ContingencyTable.ClassFactory, ClusterType.values().length);
				System.out.println("Num Event Tables: " + eventTables.size());
								
				int numSimulationIter = 10;
				for (int simulationIter = 0; simulationIter < numSimulationIter; simulationIter++) {
					LOHcateSimulator simulator = new LOHcateSimulator();
					LOHcateSimulator.LOHcateSimulatorGoldStandard goldStandard = new LOHcateSimulator.LOHcateSimulatorGoldStandard(oneSampleData.getNumSites());
					simulator.generateSimulatedDataForSample(new LOHcateSimulator.LOHcateSimulatorParams(), oneSampleData, goldStandard);

					// Declare meta-data and fill it
					ClusteringInputOneSampleMetaData metaData = new ClusteringInputOneSampleMetaData(oneSampleData.mInfoSites.size());
					classifySitesHelper_preprocessMetaData(oneSampleData, metaData, allelicBiasTable, platform);				

					System.out.printf("FDR Alpha:\t%g\nFDR Tumor p-value: \t%g\nFDR Normal p-value:\t%g\n", ParamsNum.FDRAlpha.getValue().doubleValue(), metaData.mFDRTumor, metaData.mFDRNormal);
					System.out.println("\tInferred Copy Number and Allelic Imbalance Per Site...");					

					ClusterType[] clusters = null;
					int indexFirstSomaticRowInAllVariants = 0;  // Set as invalid value for now
					if (UsePValuePlane) {
						clusters = Clustering.getClusters_withPlane(oneSampleData.mInfoSites, metaData.mImbalancePValuesTumor, metaData.mImbalancePValuesNormal, metaData.mFDRTumor, metaData.mFDRNormal, metaData.mTumorCopyNumRatiosPerGene, metaData.mCopyNumRatioPerChromNormal, metaData.mAdjustedVAFNormal, metaData.mAdjustedVAFTumor, outFilenameFullPath, indexFirstSomaticRowInAllVariants, platform); //get cluster assignments (HET ball, LOH sidelobes, DUP wedge, &c.)
					} else {
						clusters = Clustering.getClusters_Old(oneSampleData.mInfoSites, metaData.mImbalancePValuesTumor, metaData.mTumorCopyNumRatiosPerGene, outFilenameFullPath, indexFirstSomaticRowInAllVariants, platform); //get cluster assignments (HET ball, LOH sidelobes, DUP wedge, &c.)
					}

					boolean toPrintWithCopyNum = true;
					if (toPrintWithCopyNum) {
						PrintStream outStream = IOUtils.getPrintStream(outFilenameFullPath + ".withCopyNum.txt");
						for (int row = 0; row < oneSampleData.mInfoSites.size(); row++) {
							ClusteringInputOneSite oneSiteInfo = oneSampleData.mInfoSites.get(row);
							outStream.print(oneSiteInfo.printToString(sb, true, StringUtils.FileExtensionTSV.mDelimiter).toString());										
							outStream.printf("\t%s", clusters[row].name());
							outStream.printf("\t%g\t%g", metaData.mTumorCopyNumRatiosPerGene[row], (metaData.mTumorCopyNumRatiosPerGene[row] * Script.DefaultDiploidCopyNumber));  
							outStream.printf("\t%g\t%g", metaData.mAdjustedVAFTumor[row], metaData.mAdjustedVAFNormal[row]);
							outStream.printf("\t%g\t%g", metaData.mImbalancePValuesTumor[row], metaData.mImbalancePValuesNormal[row]);
							outStream.println("");
						}
						IOUtils.closePrintStream(outStream);
					}

					System.out.println("\tGot clusters");

					// Now initialize the data structure needed to plot
					int[] clusterTypeCounts = ArrayUtils.getEnumTypeCounts(clusters, ClusterType.values().length);
					DoubleParallelArray[] clusterCoordinates  = getCoordinateArraysPerClusterType(clusterTypeCounts);
					DoubleParallelArray[] waterfallPlotTumor  = getCoordinateArraysPerClusterType(clusterTypeCounts);
					DoubleParallelArray[] waterfallPlotNormal = getCoordinateArraysPerClusterType(clusterTypeCounts);
					DoubleParallelArray[] copyNumPlot         = getCoordinateArraysPerClusterType(clusterTypeCounts);

					// Do the post-processing
					Chrom prevChrom = Chrom.c0;
					int prevPosition = 0;
					long positionGenomeWide = 0;

					DoubleArrayList chromBoundaryXValue = new DoubleArrayList();
					DoubleArrayList chromBoundaryYValue = new DoubleArrayList();

					int startingRowGermlineOrSomaticOrAll = 0;  // 0 because header line has been stripped away
					for (int row = startingRowGermlineOrSomaticOrAll; row < oneSampleData.mInfoSites.size(); row++) {
						ClusteringInputOneSite oneSiteInfo = oneSampleData.mInfoSites.get(row);
						ClusterType clusterType = clusters[row];

						classifySites_addToContingencyTables(clusterType, goldStandard.mSomaticEvents.get(row), eventTables);

						//=IF(A3=A2,B3-B2+L2,B3+L2)
						Chrom chrom = oneSiteInfo.getChrom(); 							
						int position = oneSiteInfo.getPosition(); 							

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

						// Add recurrence count
						if (clusterType.isLOH()) {
							int countAtPosition = dbcByChromLOH.get(chrom.ordinal()).incrementCount(position);
							countAtPositionMax = Math.max(countAtPositionMax, countAtPosition);
						}

						float vafNormal = metaData.mAdjustedVAFNormal[row];
						float vafTumor  = metaData.mAdjustedVAFTumor[row]; 
						clusterCoordinates [clusterType.ordinal()].add( vafTumor, vafNormal );
						waterfallPlotTumor [clusterType.ordinal()].add( positionGenomeWide, vafTumor);
						waterfallPlotNormal[clusterType.ordinal()].add( positionGenomeWide, vafNormal);
						copyNumPlot        [clusterType.ordinal()].add( positionGenomeWide, (metaData.mTumorCopyNumRatiosPerGene[row] * Script.DefaultDiploidCopyNumber));

						//IOUtils.writeToBufferedWriter(out, sb.toString(), true);
						IOUtils.flushBufferedWriter(out);
					}									

					// Now let's create the datasets needed to
					double[][] boundaryArrays = combineTwoDynamicArraysIntoOneStatic(chromBoundaryXValue, chromBoundaryYValue);
					DefaultXYDataset xyDatasetVAFPlot             = classifySitesHelper_createAndFillXYData(clusterCoordinates, null);
					DefaultXYDataset xyDatasetWaterfallPlotTumor  = classifySitesHelper_createAndFillXYData(waterfallPlotTumor,  boundaryArrays);
					DefaultXYDataset xyDatasetWaterfallPlotNormal = classifySitesHelper_createAndFillXYData(waterfallPlotNormal, boundaryArrays);
					DefaultXYDataset xyDatasetCopyNumber          = classifySitesHelper_createAndFillXYData(copyNumPlot, boundaryArrays);

					Clustering.plotVAFComparison(xyDatasetVAFPlot,            vafComparisonPlotDir + File.separator + samplenameRoot + ".VAFComparison", samplenameRoot);
					Clustering.plotVAFGenomeWide(xyDatasetWaterfallPlotTumor,  vafWaterfallPlotDir + File.separator + samplenameRoot + ".VAF_GenomeWide_Tumor",  samplenameRoot, true);				
					Clustering.plotVAFGenomeWide(xyDatasetWaterfallPlotNormal, vafWaterfallPlotDir + File.separator + samplenameRoot + ".VAF_GenomeWide_Normal", samplenameRoot, false);
					Clustering.plotCopyNumGenomeWide(xyDatasetCopyNumber,        copyNumberPlotDir + File.separator + samplenameRoot + ".CopyNumber_GenomeWide", samplenameRoot);				
				}
				IOUtils.closeBufferedWriter(out);
				
				for (ClusterType eventType : ClusterType.values()) {
					ContingencyTable table = eventTables.get(eventType.ordinal());
					System.out.printf("Event:\t%10s\tSensitivity:\t%g\tSpecificity:\t%g\n", eventType, table.getSensitivity(), table.getSpecificity());						
				}	
			}
		}
		
		// Now plot the LOH recurrence across samples
		plotRecurrenceGenomeWide(dbcByChromLOH, copyNumberPlotDir, countAtPositionMax, ClusterType.LOH);

	}

	// ========================================================================
	private static DefaultXYDataset classifySitesHelper_createAndFillXYData(DoubleParallelArray[] coordinatesByEvent, double[][] boundaryArrays) {
		DefaultXYDataset xyDataset = new DefaultXYDataset();
		
		for (ClusterType eventType : ClusterType.values()) {
			xyDataset.addSeries(eventType.name(), coordinatesByEvent[eventType.ordinal()].mArray);			
		}
		
		if (boundaryArrays != null) {
			xyDataset.addSeries("Boundary", boundaryArrays);
		}
		
		return xyDataset;
	}
	
	// ========================================================================
	/** Calculates the FDR via the Benjamini-Hochberg correction and returns the p-value
	 *  threshold.
	 *  
	 *  BH correction works as follows.  Sort the p-values from least to greatest.  Then,
	 *  starting from greatest p-value (end of array), move to least p-value (beginning 
	 *  of array) such that you reach a p-value that meets the condition:
	 *  
	 *  p-value[k] <= (k / array-length) * alpha
	 *  
	 *  NOTE: Since this is mathematics, k starts at 1, not 0
	 *  
	 */
	public static double getFDR_BenjaminiHochberg(double[] pValues, double fdrAlpha, boolean pointsIndependent) {
		double[] pValuesClone = pValues.clone();
		Arrays.sort(pValuesClone);
		
		double c_m = 1.0;  // Default under condition of independence
		
		if (!pointsIndependent) {
			double eulerMascheroniConstant = 0.577215664901532;
			double e_m = 1.0 / (2 * pValuesClone.length);
			c_m = Math.log(pValuesClone.length) + eulerMascheroniConstant + e_m;
		}
		
		for (int k = pValuesClone.length - 1; k >= 0; k--) {
			double threshold = ((double) (k + 1) / ((double) pValuesClone.length * c_m)) * fdrAlpha;
			if (pValuesClone[k] <= threshold) {
				return pValuesClone[k];
			}
		}
		
		return 0;
	}
	
	// ========================================================================
	/** Plots the event recurrence genome-wide. */
	public static void plotRecurrenceGenomeWide(ArrayList<DynamicBucketCounter> dbcCountByEvent, String outDir, int maxValue, ClusterType clusterType) {
		DoubleArrayList positionsGenomeWideAllSamples = new DoubleArrayList();
		DoubleArrayList lohCountAllSamples            = new DoubleArrayList();
		double lastPositionOnPrevChrom = 0;
		
		DoubleArrayList chromBoundaryXValue = new DoubleArrayList();
		DoubleArrayList chromBoundaryYValue = new DoubleArrayList();
		PrimitiveWrapper.WInteger lastKey = new PrimitiveWrapper.WInteger(0);
				
		for (Chrom chrom : Chrom.values()) {
			if (chrom.isAutosomal()) {
				DoubleArrayList[] positionAndCount = dbcCountByEvent.get(chrom.ordinal()).toArrayListDouble();
				for (int i = 0; i < positionAndCount[0].size(); i++) {
					positionAndCount[0].set(i, positionAndCount[0].get(i) + lastPositionOnPrevChrom);
				}
				
				positionsGenomeWideAllSamples.addAll(positionAndCount[0]);
				lohCountAllSamples.addAll           (positionAndCount[1]);
				
				boolean lastKeyExists = dbcCountByEvent.get(chrom.ordinal()).getKeyLast(lastKey);
				lastPositionOnPrevChrom += (lastKeyExists) ? lastKey.mInt : 30000000;
				
				double increment = (maxValue) / 50.0;
				for (double d = 0; d <= maxValue; d += increment) {
					chromBoundaryXValue.add(lastPositionOnPrevChrom);
					chromBoundaryYValue.add(d);
				}
			}
		}
		DefaultXYDataset posAndLOHCountDataset = new DefaultXYDataset();
		posAndLOHCountDataset.addSeries(ClusterType.LOH.name(), combineTwoDynamicArraysIntoOneStatic(positionsGenomeWideAllSamples, lohCountAllSamples));
		posAndLOHCountDataset.addSeries("Boundary",             combineTwoDynamicArraysIntoOneStatic(chromBoundaryXValue, chromBoundaryYValue));
		Clustering.plotEventSampleRecurrence(posAndLOHCountDataset, outDir + File.separator + "All_Samples." + clusterType.name());
	}
	
	// ========================================================================
	private static double[][] combineTwoDynamicArraysIntoOneStatic(DoubleArrayList list1, DoubleArrayList list2) {
		return new double[][] { list1.toArray(), list2.toArray() };
	}
	
	// ========================================================================
	private static DoubleParallelArray[] getCoordinateArraysPerClusterType(int[] clusterTypeCounts) {
		
		DoubleParallelArray[] clusterCoordinates = new DoubleParallelArray[clusterTypeCounts.length];
		
		for (int clusterTypeIndex = 0; clusterTypeIndex < clusterCoordinates.length; clusterTypeIndex++) {
			clusterCoordinates[clusterTypeIndex] = new DoubleParallelArray(clusterTypeCounts[clusterTypeIndex]);
		}
		
		return clusterCoordinates;
	}

	// ========================================================================
	public static void plotEventSampleRecurrence(XYDataset xyDataset, String outFilenameRoot) {
		String xAxisLabel = "Position";
		String yAxisLabel = "# Samples with Event";
		String title = "Recurrence over Samples for Event: " + xyDataset.getSeriesKey(0);				
		
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
		legendTitle.setID("Event");
		legendTitle.setItemFont(new Font("Arial", Font.BOLD, 20));

		//xyPlot.getRangeAxis().setRange(0, 1.02);		
		
		Font rangeAxisLabelFont = new Font("Arial", Font.BOLD, 20);
		Font rangeAxisTickFont = new Font("Arial", Font.BOLD, 20);
		xyPlot.getRangeAxis().setLabelFont(rangeAxisLabelFont);		
		xyPlot.getRangeAxis().setTickLabelFont(rangeAxisTickFont);
		xyPlot.getDomainAxis().setLabelFont(rangeAxisLabelFont);
		xyPlot.getDomainAxis().setTickLabelFont(rangeAxisTickFont);	
		GraphUtils.saveChartAsPNG(outFilenameRoot, theChart, 2400, 800);
		
	}
	
	// ========================================================================
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
		GraphUtils.saveChartAsPNG(outFilenameRoot, theChart, 2400, 800);
	}
	
	// ========================================================================
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
	
	// ========================================================================
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

	// ========================================================================
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
	private static void setSeriesPaintPerCluster(XYItemRenderer itemRenderer) {
		boolean allGray = false;
		
		if (allGray) {
			for (ClusterType eventType : ClusterType.values()) {
				itemRenderer.setSeriesPaint(eventType.ordinal(), ColorPastel.Gray_50.getColor());
			}
		} else {
			itemRenderer.setSeriesPaint(ClusterType.GainGermline.ordinal(), ColorPastel.Violet.getColor());
			itemRenderer.setSeriesPaint(ClusterType.GainSomatic.ordinal(), ColorPastel.Dark_Red.getColor());
			itemRenderer.setSeriesPaint(ClusterType.LOH.ordinal(), ColorPastel.RGB_Blue.getColor());
			itemRenderer.setSeriesPaint(ClusterType.cnLOH.ordinal(), ColorPastel.CMYK_Yellow.getColor());
			itemRenderer.setSeriesPaint(ClusterType.HETGermline.ordinal(), ColorPastel.Gray_50.getColor());
			itemRenderer.setSeriesPaint(ClusterType.HETSomatic.ordinal(), ColorPastel.Red_Orange.getColor());
			itemRenderer.setSeriesPaint(ClusterType.Noise.ordinal(), ColorPastel.Gray_60.getColor());
			itemRenderer.setSeriesPaint(ClusterType.Null.ordinal(), ColorPastel.Gray_30.getColor());				
		}
		itemRenderer.setSeriesPaint(ClusterType.Null.ordinal() + 1, ColorPastel.Black.getColor());
	}

	// ========================================================================
	private static void calculateAdjustedVAFs(ArrayList<ClusteringInputOneSite> rows, AllelicBiasTable allelicBiasTable, 
											  float[] adjustedVAFNormal, float[] adjustedVAFTumor,
											  float[] copyNumRatioPerChromNormal,
											  SeqPlatform platform) {
		
		final float defaultVAFNormal = 0.50f;
		String delim = StringUtils.FileExtensionTSV.mDelimiter;
		
		for (int row = 0; row < rows.size(); row++) {
			ClusteringInputOneSite oneSiteInfo = rows.get(row);
			
			Chrom chrom  = oneSiteInfo.getChrom();
			int position = oneSiteInfo.getPosition();
			float vafNormal = oneSiteInfo.calcVAFNormal();
			float vafTumor  = oneSiteInfo.calcVAFTumor();
			float vafNormalExpected = defaultVAFNormal;
			
			// Set Default values
			float offset = 0.0f;  // Default
			float adjustmentFactor = 1.0f;
			
			// Perform adjustment calculations
			if (correctAllelicBias()) {
				if (NumberUtils.inRangeLowerExclusive(vafNormal, AlleleFrequencyStatsForSample.VAFNormalFrameLower, AlleleFrequencyStatsForSample.VAFNormalFrameUpper)) {
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
						} else {
							adjustmentFactor = vafNormalExpected / avgVAFNormal;
							float absDiff = Math.abs(vafNormalExpected - avgVAFNormal);
							offset = (vafNormal > vafNormalExpected) ? -absDiff : ((vafNormal == vafNormalExpected) ? 0 : absDiff);
						}
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
	public static void getPValuesImbalance(ArrayList<ClusteringInputOneSite> rows, 
			double[] imbalancePValuesTumor, double[] imbalancePValuesNormal,
			float[] adjustedVAFNormal, float[] adjustedVAFTumor) {
		
		for (int row = 0; row < rows.size(); row++) {
			ClusteringInputOneSite oneSiteInfo = rows.get(row);
												
			int varCovgTumor =  (int) Math.round(oneSiteInfo.mCovgTotalTumor  * adjustedVAFTumor[row]);
			int varCovgNormal = (int) Math.round(oneSiteInfo.mCovgTotalNormal * adjustedVAFNormal[row]);
			
			imbalancePValuesTumor[row]  = getPValuesImbalanceTissue(oneSiteInfo.mCovgTotalTumor,  varCovgTumor);
			imbalancePValuesNormal[row] = getPValuesImbalanceTissue(oneSiteInfo.mCovgTotalNormal, varCovgNormal);		
		}
	}

	// ------------------------------------------------------------------------
	private static double getPValuesImbalanceTissue(int coverageTotal, int coverageVariant) {
		int maxRefOrVarCovg = Math.max(coverageVariant, coverageTotal - coverageVariant);
		return nutils.NumberUtils.cumulativeProbabilitySuccess(coverageTotal, maxRefOrVarCovg, 0.5);
	}

	// ========================================================================
	public static float[] calcRoughCopyNumberRatioPerSite(ArrayList<ClusteringInputOneSite> rows, float[] copyNumRatios, float coverageRatioTumorNormalGenomeWide, DynamicBucketCounter readCountTalliesTumor, DynamicBucketCounter readCountTalliesNormal) {		
		boolean[] isSomaticSite = new boolean[rows.size()];
		
		String prevGene = "";
		int rowOfFirstInstanceOfGene = -1;
		
		float ratioSum = 0;
		int readCountSumTumor = 0;
		int readCountSumNormal = 0;
		int numRowsWithSameGene = 0;
		
		int avgCoverageNormal = readCountTalliesNormal.getKeyWithMaxCount();
		int avgCoverageTumor  = readCountTalliesTumor.getKeyWithMaxCount();
		double avgCoverageHaploidNormal = avgCoverageNormal / 2.0;
		double avgCoverageHaploidTumor  = avgCoverageTumor  / 2.0;
		
		int numPoissons = 10;
		PoissonDistribution[] pdNormal = new PoissonDistribution[numPoissons];
		PoissonDistribution[] pdTumor  = new PoissonDistribution[numPoissons];
		for (int i = 0; i < numPoissons; i++) {
			pdNormal[i] = new PoissonDistribution(avgCoverageHaploidNormal * (i + 1));
			pdTumor[i]  = new PoissonDistribution(avgCoverageHaploidTumor  * (i + 1));
		}		
		
		for (int row = 0; row < rows.size(); row++) {
			ClusteringInputOneSite oneSiteInfo = rows.get(row);

			String currentGene = oneSiteInfo.mHugoSymbol;					

			isSomaticSite[row] = (oneSiteInfo.mCovgTotalNormal <= 0);
			float ratio = isSomaticSite[row] ? Script.TumorNormalRatioOfSomaticSite : ((float) oneSiteInfo.mCovgTotalTumor / (float) oneSiteInfo.mCovgTotalNormal);			
			
			if (currentGene.equalsIgnoreCase(prevGene)) {
				if (!isSomaticSite[row]) {
					numRowsWithSameGene++;
					ratioSum += ratio;
					readCountSumTumor  += oneSiteInfo.mCovgTotalTumor;
					readCountSumNormal += oneSiteInfo.mCovgTotalNormal;
				} else {
					// We have a somatic site, we take no action in this block
				}
			} else {
				// New gene listed.  We must now write the ratios for the previous gene
				if (row > 0) {					
					float averageRatio = Script.TumorNormalRatioOfSomaticSite;
					
					if (numRowsWithSameGene == 0) {
						CompareUtils.ensureTrue(false, "ERROR: Gene must > 0 sites representative of that gene!");
					} else if (numRowsWithSameGene == 1 && !AllowGeneWithOneSite) {					
						averageRatio = Script.DefaultTumorNormalRatio;
					} else if (numRowsWithSameGene > 1 || AllowGeneWithOneSite) {
						double readCountAverageTumor  = (double) readCountSumTumor  / (double) numRowsWithSameGene;
						double readCountAverageNormal = (double) readCountSumNormal / (double) numRowsWithSameGene;
						int indexDistTumor  = getIndexOfPoissonWithHighestProbability((int) Math.round(readCountAverageTumor),  pdTumor);
						int indexDistNormal = getIndexOfPoissonWithHighestProbability((int) Math.round(readCountAverageNormal), pdNormal);
						averageRatio = (float) (indexDistTumor + 1) / (float) (indexDistNormal + 1);
						
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
					readCountSumTumor  = oneSiteInfo.mCovgTotalTumor;
					readCountSumNormal = oneSiteInfo.mCovgTotalNormal;
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
	private static final int getIndexOfPoissonWithHighestProbability(int k, PoissonDistribution[] dists) {
		double maxProb = -Double.MAX_VALUE;
		int    maxProbIndex = -1;
				
		for (int i = 0; i < dists.length; i++) {
			double prob = dists[i].probability(k);
			if (prob > maxProb) {
				maxProb = prob;
				maxProbIndex = i;
			}
		}
		
		return maxProbIndex;
	}
	
	// ========================================================================
	private static final int roundToNearest5(int num) {
		return (int) (Math.round(num / 5.0) * 5);
	}
	
	// ========================================================================
	/** Tumor : Normal avg coverage ratio per chromosome. */
	public static float[] calcAvgCoverageRatioPerChrom(ArrayList<ClusteringInputOneSite> rows, float[] copyNumRatioPerChromNormal, int[] numSitesPerChrom, float[] avgReadCountPerChromNormal, PrimitiveWrapper.WFloat coverageRatioTumorToNormal, DynamicBucketCounter readCountTalliesTumor, DynamicBucketCounter readCountTalliesNormal) {			
		float[] tumorNormalRatioPerChrom = new float[Chrom.values().length];
		int[] totalReadCountPerChromNormal = new int[Chrom.values().length];
		int[] totalReadCountPerChromTumor  = new int[Chrom.values().length];
		
		Arrays.fill(tumorNormalRatioPerChrom, Script.DefaultTumorNormalRatio);  // Fill with 1.0 since that's an equal ratio
		Arrays.fill(copyNumRatioPerChromNormal, Script.DefaultNormalRatio);
		
		Arrays.fill(numSitesPerChrom, 0);  // initialize counts
		Arrays.fill(avgReadCountPerChromNormal, 0); // initialize
		Arrays.fill(totalReadCountPerChromNormal, 0);
		Arrays.fill(totalReadCountPerChromTumor,  0);
		
		readCountTalliesTumor.clear();
		readCountTalliesNormal.clear();
		
		for (int row = 0; row < rows.size(); row++) {
			ClusteringInputOneSite oneSiteInfo = rows.get(row);
			
			int chromIndex = oneSiteInfo.getChrom().ordinal();
			totalReadCountPerChromNormal[chromIndex] += oneSiteInfo.mCovgTotalNormal;
			totalReadCountPerChromTumor[chromIndex]  += oneSiteInfo.mCovgTotalTumor;
			++numSitesPerChrom[chromIndex];
			
			//=ROUND(A1*2/10,0)*10/2
			
			readCountTalliesNormal.incrementCount(roundToNearest5(oneSiteInfo.mCovgTotalNormal));			
			readCountTalliesTumor.incrementCount(roundToNearest5(oneSiteInfo.mCovgTotalTumor));
					
		}
		
		// Calculate the genome-wide average read count in the normal
		int totalReadCountNormal = ArrayUtils.arraySum(totalReadCountPerChromNormal);
		int totalReadCountTumor  = ArrayUtils.arraySum(totalReadCountPerChromTumor);
		int totalNumSites        = ArrayUtils.arraySum(numSitesPerChrom);
		float avgCoverageNormal = (float) totalReadCountNormal / (float) totalNumSites;
		float avgCoverageTumor  = (float) totalReadCountTumor  / (float) totalNumSites;
		
		// Re-assign to highest count
		avgCoverageNormal = (float) readCountTalliesNormal.getKeyWithMaxCount();
		avgCoverageTumor  = (float) readCountTalliesTumor.getKeyWithMaxCount();
		
		coverageRatioTumorToNormal.mFloat = avgCoverageTumor / avgCoverageNormal;		 
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

	// ========================================================================
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
	public static ClusterType[] getClusters_Old(ArrayList<ClusteringInputOneSite> rows, double[] imbalancePValues, float[] copyNumRatios, String outFilenameFullPath, int startingRowGermlineOrSomaticOrAll, SeqPlatform platform) {

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
	public static ClusterType[] getClusters_withPlane(ArrayList<ClusteringInputOneSite> rows, 
			double[] imbalancePValuesTumor, double[] imbalancePValuesNormal, 
			double fdrTumor, double fdrNormal,
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
				                    fdrTumor, fdrNormal,
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
			
			ClusteringInputOneSite oneSiteInfo = rows.get(indexIntoOriginalRows);
			Chrom chrom = oneSiteInfo.getChrom();
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
	private static ClusterType[] assignClusters(ArrayList<ClusteringInputOneSite> rows, 
			SeqPlatform platform,
			//int startingRowSomatic,												
			float vafBoundLower, 
			float vafBoundUpper,
			int[] clusterAssignments,
			int[] clusterTypeIDsFromAlgorithm) {

		int indexInClusterAssignments = -1;  // we need to keep a special index, since not all rows are used.
		ClusterType[] returnClusters = new ClusterType[rows.size()];

		for (int row = 0; row < rows.size(); row++) {

			ClusteringInputOneSite oneSiteInfo = rows.get(row);
			float vafNormal = oneSiteInfo.calcVAFNormal();
			boolean vafInRangeNormal = NumberUtils.inRangeLowerExclusive(vafNormal, vafBoundLower, vafBoundUpper);

			if (/*(row >= startingRowSomatic) &&*/ (!vafInRangeNormal)) {
				// The vafNormal is either very low (homozygous reference) or very high (homozygous common variant).
				// We do some very simple decision making now (which should be replaced by formal clustering later)
				// to partition the calls.
				float justBelowZero = -0.0001f;				
				float hetBoundaryLower = 0.3333f;
				float hetBoundaryUpper = 0.6667f;
				float vafTumor = oneSiteInfo.calcVAFTumor();
				if (NumberUtils.inRangeLowerExclusive(vafNormal, justBelowZero, vafBoundLower)) {
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
				} else if (NumberUtils.inRangeLowerExclusive(vafNormal, vafBoundUpper, Script.MaxVariantAlleleFrequency)) {
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
					CompareUtils.throwErrorAndExit("ERROR: Contradiction - variant allele frequency cannot be in and out of bounds simultanteously!" + vafNormal);
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
	private static void getValidPointsForClustering(ArrayList<ClusteringInputOneSite> rows, // input
													ArrayList<Floint> pointsUpperPlane,    // output
													ArrayList<Floint> pointsLowerPlane,    // output
													int[] indexMapToUpperLowerPlane,       // output - if positive, to upper plane, if negative, to lower plane, if does not map, contains Integer.MAX_VALUE													
													IntArrayList mapToRowsFromUpper, 
													IntArrayList mapToRowsFromLower,
													double[] imbalancePValuesTumor,
													double[] imbalancePValuesNormal,
													double fdrTumor, 
													double fdrNormal,
													double scalingFactor,
													float[] copyNumRatios, float[] copyNumRatioPerChromNormal,
													float[] adjustedVAFNormal, float[] adjustedVAFTumor, 
													SeqPlatform platform, 
													float vafBoundLower, float vafBoundUpper) {
		
		pointsUpperPlane.clear();
		pointsLowerPlane.clear();
//		double alphaAdjusted = ParamsBool.IgnoreMultipleTesting.getValue() ? 
//				 PValueBinDistAlpha_UpperPlaneThresh :
//				(PValueBinDistAlpha_UpperPlaneThresh * 0.008 /* Light FDR */);
				
		for (int row = 0; row < rows.size(); row++) {
			ClusteringInputOneSite oneSiteInfo = rows.get(row);			
			
			float vafNormal = oneSiteInfo.calcVAFNormal();
			if (NumberUtils.inRangeLowerExclusive(vafNormal, vafBoundLower, vafBoundUpper)) {				
				Floint thePoint = new Floint(adjustedVAFTumor[row], adjustedVAFNormal[row], (float) (imbalancePValuesTumor[row] * scalingFactor));
				
				// Assume p-value as vertical row factor
				boolean tumorSigImbalanced = (imbalancePValuesTumor[row] <= fdrTumor); // && !isCopyNumRatioInDiploidRange(copyNumRatios[row]);
				Chrom chrom = oneSiteInfo.getChrom();
				boolean normalSigImbalanced = (imbalancePValuesNormal[row] <= fdrNormal) && (copyNumRatioPerChromNormal[chrom.ordinal()] > ParamsNum.GermlineTrisomyThreshold.getValue().doubleValue());
				
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
		getValidPointsListForClustering(ArrayList<ClusteringInputOneSite> rows, double[] verticalFactor, double scalingFactor, SeqPlatform platform, float vafBoundLower, float vafBoundUpper) {
		
		ArrayList<Floint> points = new ArrayList<Floint>(rows.size());

		for (int row = 0; row < rows.size(); row++) {
			ClusteringInputOneSite oneSiteInfo = rows.get(row);			
			float vafNormal = oneSiteInfo.calcVAFNormal();
			if (NumberUtils.inRangeLowerExclusive(vafNormal, vafBoundLower, vafBoundUpper)) {	
				Floint thePoint = new Floint(oneSiteInfo.calcVAFTumor(), vafNormal, (float) (verticalFactor[row] * scalingFactor));
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
		case Illumina: return  Float.parseFloat(StringUtils.extractNthColumnValue(line, Script.Col_NAFTAFInput_VariantRatioNormal, Utils.TabStr));
		case SOLiD:    return (Float.parseFloat(StringUtils.extractNthColumnValue(line, 6, Utils.TabStr)) / 
				               Float.parseFloat(StringUtils.extractNthColumnValue(line, 5, Utils.TabStr)));
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
		case Illumina: return  Float.parseFloat(StringUtils.extractNthColumnValue(line, Script.Col_NAFTAFInput_VariantRatioTumor, Utils.TabStr));
		case SOLiD:    return (Float.parseFloat(StringUtils.extractNthColumnValue(line, 4, Utils.TabStr)) / 
				               Float.parseFloat(StringUtils.extractNthColumnValue(line, 3, Utils.TabStr)));
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

	// ========================================================================
	// ========================================================================
	public static class ClusteringInputOneSample implements RegionSimulator.SampleInformation<ClusteringInputOneSite> {
		ArrayList<ClusteringInputOneSite> mInfoSites;
		ClusteringInputOneSite mDummySite;
		
		public ClusteringInputOneSample(int numSitesEstimated) {
			mInfoSites = new ArrayList<ClusteringInputOneSite>(numSitesEstimated);
			constructorCommon();
		}
		
		public ClusteringInputOneSample(ArrayList<String> rows) {
			mInfoSites = new ArrayList<ClusteringInputOneSite>(rows.size());
			constructorCommon();
			parseLines(rows);
		}
		
		private void constructorCommon() {
			mDummySite = new ClusteringInputOneSite();
		}
		
		private void parseLines(ArrayList<String> rows) {			
			for (String line : rows) {
				ClusteringInputOneSite cios = new ClusteringInputOneSite(line, SeqPlatform.Illumina);
				mInfoSites.add(cios);
			}
		}
		
		/** Returns the index of the first site (on the chromosome) in this sample.  If
		 *  the chromosome does not exist, a negative index is returned in accordance to
		 *  the definition from Collections.binarySearch()
		 * @param chrom
		 * @return
		 */
		public synchronized int getIndexChromStart(Chrom chrom) {
			int resultIndex = getIndex(chrom, 1);
			if (resultIndex >= 0) {
				return resultIndex;
			} else {				
				int insertPoint = ArrayUtils.getInsertPoint(resultIndex);
				if (insertPoint >= mInfoSites.size()) return resultIndex;
				
				ClusteringInputOneSite oneSiteInfo = mInfoSites.get(insertPoint);
				return (oneSiteInfo.getChrom().equals(chrom) ? insertPoint : resultIndex);
			}
		}
		
		public synchronized int getIndexChromEnd(Chrom chrom) {
			Chrom nextChrom = chrom.getNextChrom();
			if (nextChrom == null) {
				// There is no next chromosome
				return mInfoSites.size() - 1;
			} else {
				// There is a position for a next chromosome
				int resultIndex = getIndexChromStart(nextChrom);
				if (resultIndex < 0) {
					// This means the next chromosome doesn't exist in the dataset.  
					// Try the following chromosome via recursion
					return getIndexChromEnd(nextChrom);  
				} else {
					return resultIndex - 1;
				}
			}
		}
		
		public synchronized int getIndex(Chrom chrom, int position) {
			mDummySite.setChrom(chrom);
			mDummySite.setPosition(position);
			return Collections.binarySearch(mInfoSites, mDummySite, ClusteringInputOneSite.ClusteringInputOneSiteComparator);			
		}
		
		public synchronized ClusteringInputOneSite getSiteAtIndex(Chrom chrom, int index) {
			return getSiteAtIndex(index);
		}
		
		public synchronized ClusteringInputOneSite getSiteAtIndex(int index) {
			return mInfoSites.get(index);
		}
		
		public synchronized int getNumSites() { return mInfoSites.size(); }
	}
	
	// ========================================================================
	// ========================================================================
	
	public static class ClusteringInputOneSite implements Comparable<ClusteringInputOneSite>, RegionSimulator.SiteInformation {
		
		public static Comparator<ClusteringInputOneSite> ClusteringInputOneSiteComparator = new Comparator<ClusteringInputOneSite>() {
			public int compare(ClusteringInputOneSite site1, ClusteringInputOneSite site2) {
				//System.out.println(site1.getChrom() + "\t" + site1.getPosition() + "\t" + site2.getChrom() + "\t" + site2.getPosition());
				int result = site1.getChrom().compareTo(site2.getChrom()); 
				if (result == 0) {
					result = Integer.compare(site1.getPosition(), site2.getPosition());
				}
				return result;
			}
		};
		
		// ========================================================================
		private long mDataUnit_ChromProsRevVarAllelesMutType;
		
		private static final BitSetUtils.BitShiftAndMask bsmChrom = new BitShiftAndMask(5, 58);
		private static final BitSetUtils.BitShiftAndMask bsmPos   = BitShiftAndMask.createBitShiftAndMaskInChain(28, bsmChrom);
		private static final BitSetUtils.BitShiftAndMask bsmAlleleRef  = BitShiftAndMask.createBitShiftAndMaskInChain(3, bsmPos);
		private static final BitSetUtils.BitShiftAndMask bsmAlleleVarN = BitShiftAndMask.createBitShiftAndMaskInChain(3, bsmAlleleRef);
		private static final BitSetUtils.BitShiftAndMask bsmAlleleVarT = BitShiftAndMask.createBitShiftAndMaskInChain(3, bsmAlleleVarN);
		private static final BitSetUtils.BitShiftAndMask bsmAlleleVarPop = BitShiftAndMask.createBitShiftAndMaskInChain(3, bsmAlleleVarT);		
		private static final BitSetUtils.BitShiftAndMask bsmMutType    = BitShiftAndMask.createBitShiftAndMaskInChain(5, bsmAlleleVarPop);
		
		public String mFlankingNormal;
		public String mFlankingTumor;
		public String mHugoSymbol;
		
		public short mCovgTotalNormal;
		public short mCovgTotalTumor;
		public short mCovgVarNormal;
		public short mCovgVarTumor;
		public int mRsID; 
		
		public ClusteringInputOneSite() {
			clear();
		}			
		
		public ClusteringInputOneSite(String line, SeqPlatform platform) {
			mFlankingNormal = StringUtils.extractNthColumnValue(line, Script.Col_NAFTAFInput_FlankingStringNormal,  StringUtils.FileExtensionTSV.mDelimiter);
			mFlankingTumor  = StringUtils.extractNthColumnValue(line, Script.Col_NAFTAFInput_FlankingStringTumor,   StringUtils.FileExtensionTSV.mDelimiter);
			
			mCovgTotalNormal = Short.parseShort(StringUtils.extractNthColumnValue(line, Script.Col_NAFTAFInput_TotalCoverageNormal,   StringUtils.FileExtensionTSV.mDelimiter));
			mCovgTotalTumor  = Short.parseShort(StringUtils.extractNthColumnValue(line, Script.Col_NAFTAFInput_TotalCoverageTumor,    StringUtils.FileExtensionTSV.mDelimiter));
			mCovgVarNormal   = Short.parseShort(StringUtils.extractNthColumnValue(line, Script.Col_NAFTAFInput_VariantCoverageNormal, StringUtils.FileExtensionTSV.mDelimiter));
			mCovgVarTumor    = Short.parseShort(StringUtils.extractNthColumnValue(line, Script.Col_NAFTAFInput_VariantCoverageTumor,  StringUtils.FileExtensionTSV.mDelimiter));

			// Now for rsID
			if (platform == SeqPlatform.SOLiD) {
				//variantAnnotation = (germCols[7].indexOf(Script.NovelStr) >= 0) ? Script.NovelStr : Utils.NAStr;
				// Sidd: For the NAStr case, strangely, the variant base is n/a in the SOLiD naf-taf-inputs 
				// (and there's not much point in looking up the reference base's allele frequency)
			} else if (platform == SeqPlatform.Illumina) {
				String dbsnpStr = StringUtils.extractNthColumnValue(line, Script.Col_NAFTAFInput_DbSNPString, StringUtils.FileExtensionTSV.mDelimiter); 
				if (dbsnpStr.indexOf(Script.NovelStr) >= 0) {
					mRsID = GenotypeUtils.RsID_Novel;
				} else {
					String rsNumTemp = GenotypeUtils.extractRsNumberFromLine(dbsnpStr);
					
					if ((rsNumTemp == null) || 
						 rsNumTemp.equals(Utils.rsPrefix) || 
						 rsNumTemp.equals(Utils.rsNegative)) {
						
						mRsID = GenotypeUtils.RsID_Unknown;
					} else {
						mRsID = GenotypeUtils.getNumberFromRsId(rsNumTemp);
					}
				}
			} else {
				mRsID = GenotypeUtils.RsID_Unknown;  // Could change if another platform added later
			}
		
			// Now for Hugo symbol and mutation type		
			if (platform == SeqPlatform.SOLiD) {
				String mutationTypeStr = StringUtils.extractNthColumnValue(line, 9, StringUtils.FileExtensionTSV.mDelimiter);
				mutationTypeStr.replace("syn", "synonymous").replace("nsynonymous", "nonsynonymous");
				mutationTypeStr = (mutationTypeStr == "") ? Utils.NAStr : mutationTypeStr;			
				mHugoSymbol = StringUtils.extractNthColumnValue(line, 8, StringUtils.FileExtensionTSV.mDelimiter);
				mHugoSymbol = (mHugoSymbol == "") ? Utils.NAStr : mHugoSymbol;

			} else if (platform == SeqPlatform.Illumina) {			
				String mutationTypeStr = StringUtils.extractNthColumnValue(line, Script.Col_NAFTAFInput_MutationType, StringUtils.FileExtensionTSV.mDelimiter);
				int mutTypeCode = MutationType.getSNVType(mutationTypeStr).ordinal();
				mDataUnit_ChromProsRevVarAllelesMutType = bsmMutType.setValueInCompactUnit(mutTypeCode, mDataUnit_ChromProsRevVarAllelesMutType);
				
				mHugoSymbol     = StringUtils.extractNthColumnValue(line, Script.Col_NAFTAFInput_HugoSymbol, StringUtils.FileExtensionTSV.mDelimiter);
				
			}
		
			// Now set the chrom, position, and alleles
			Chrom chrom  = Chrom.getChrom(   StringUtils.extractNthColumnValue(line, Script.Col_NAFTAFInput_Chrom,     StringUtils.FileExtensionTSV.mDelimiter) );
			int position = Integer.parseInt( StringUtils.extractNthColumnValue(line, Script.Col_NAFTAFInput_Position,  StringUtils.FileExtensionTSV.mDelimiter) );
			
			char nucChar = StringUtils.extractNthColumnValue(line, Script.Col_NAFTAFInput_AlleleRef,  StringUtils.FileExtensionTSV.mDelimiter).charAt(0);
			Nuc aRef     = nucChar == Script.MissingAllele ? Nuc.N : Nuc.getNuc(nucChar);
			
			nucChar      = StringUtils.extractNthColumnValue(line, Script.Col_NAFTAFInput_AlleleVarN, StringUtils.FileExtensionTSV.mDelimiter).charAt(0);
			Nuc aVarN    = nucChar == Script.MissingAllele ? Nuc.N : Nuc.getNuc(nucChar);
			
			nucChar      = StringUtils.extractNthColumnValue(line, Script.Col_NAFTAFInput_AlleleVarT, StringUtils.FileExtensionTSV.mDelimiter).charAt(0);
			Nuc aVarT    = nucChar == Script.MissingAllele ? Nuc.N : Nuc.getNuc(nucChar);
			
			nucChar      = StringUtils.extractNthColumnValue(line, Script.Col_NAFTAFInput_AlleleVarPop, StringUtils.FileExtensionTSV.mDelimiter).charAt(0);
			Nuc aVarPop  = nucChar == Script.MissingAllele ? Nuc.N : Nuc.getNuc(nucChar);
			
			setChrom(chrom);
			setPosition(position);			
			mDataUnit_ChromProsRevVarAllelesMutType = bsmAlleleRef.setValueInCompactUnit( aRef.getCode(),  mDataUnit_ChromProsRevVarAllelesMutType);
			mDataUnit_ChromProsRevVarAllelesMutType = bsmAlleleVarN.setValueInCompactUnit(aVarN.getCode(), mDataUnit_ChromProsRevVarAllelesMutType);
			mDataUnit_ChromProsRevVarAllelesMutType = bsmAlleleVarT.setValueInCompactUnit(aVarT.getCode(), mDataUnit_ChromProsRevVarAllelesMutType);
			mDataUnit_ChromProsRevVarAllelesMutType = bsmAlleleVarPop.setValueInCompactUnit(aVarPop.getCode(), mDataUnit_ChromProsRevVarAllelesMutType);
		}

		// ====================================================================
		public void setChrom(Chrom chrom) {
			mDataUnit_ChromProsRevVarAllelesMutType = bsmChrom.setValueInCompactUnit(chrom.ordinal(), mDataUnit_ChromProsRevVarAllelesMutType);
		}
		
		// ====================================================================
		public void setPosition(int position) {			
			mDataUnit_ChromProsRevVarAllelesMutType = bsmPos.setValueInCompactUnit(position,  mDataUnit_ChromProsRevVarAllelesMutType);
		}
		
		// ====================================================================
		public void setCovgTotalNormal(short mCovgTotalNormal) {
			this.mCovgTotalNormal = mCovgTotalNormal;
		}

		// ====================================================================
		public void setCovgTotalTumor(short mCovgTotalTumor) {
			this.mCovgTotalTumor = mCovgTotalTumor;
		}

		// ====================================================================
		public void setCovgVarNormal(short mCovgVarNormal) {
			this.mCovgVarNormal = mCovgVarNormal;
		}

		// ====================================================================
		public void setCovgVarTumor(short mCovgVarTumor) {
			this.mCovgVarTumor = mCovgVarTumor;
		}
		
		// ====================================================================
		public float calcVAFNormal() { return (float) mCovgVarNormal / (float) mCovgTotalNormal; }
		public float calcVAFTumor()  { return (float) mCovgVarTumor  / (float) mCovgTotalTumor;  }
		
		public Chrom getChrom()  { return Chrom.getChrom((byte) bsmChrom.getValueInCompactUnit(mDataUnit_ChromProsRevVarAllelesMutType)); }
		public int getPosition() { return                ((int)   bsmPos.getValueInCompactUnit(mDataUnit_ChromProsRevVarAllelesMutType)); }
		
		public Nuc getReferenceAllele()     { return Nuc.getAllele(bsmAlleleRef.getValueInCompactUnit(mDataUnit_ChromProsRevVarAllelesMutType));  }
		public Nuc getVariantAlleleNormal() { return Nuc.getAllele(bsmAlleleVarN.getValueInCompactUnit(mDataUnit_ChromProsRevVarAllelesMutType)); }
		public Nuc getVariantAlleleTumor()  { return Nuc.getAllele(bsmAlleleVarT.getValueInCompactUnit(mDataUnit_ChromProsRevVarAllelesMutType)); }
		public Nuc getVariantAllelePopulation() { return Nuc.getAllele(bsmAlleleVarPop.getValueInCompactUnit(mDataUnit_ChromProsRevVarAllelesMutType)); }
		
		public MutationType getMutationType() {
			return MutationType.getMutationType((int) bsmMutType.getValueInCompactUnit(mDataUnit_ChromProsRevVarAllelesMutType));
		}

		// ====================================================================
		public void clear() {
			mDataUnit_ChromProsRevVarAllelesMutType = 0;
			mFlankingTumor = mFlankingNormal = mHugoSymbol = "";
			mCovgTotalNormal = mCovgTotalTumor = mCovgVarNormal = mCovgVarTumor = 0;
			mRsID = GenotypeUtils.RsID_Unknown;
		}
		
		// ====================================================================
		public int compareTo(ClusteringInputOneSite rhs) {
			// The compact unit variable has chromosome and position ordered such that it's easy to sort
			return Long.compare(mDataUnit_ChromProsRevVarAllelesMutType, rhs.mDataUnit_ChromProsRevVarAllelesMutType);
		}
		
		// ====================================================================
		
		public StringBuilder printToString(StringBuilder sb, boolean clearStringBuilder, String delimiter) {
			if (clearStringBuilder) {
				sb.setLength(0);
			}
			
			sb.append("chr").append(getChrom().ordinal())
			  .append(delimiter).append(getPosition())
			  .append(delimiter).append(getReferenceAllele())
			  .append(delimiter).append(getVariantAllelePopulation())
			  .append(delimiter).append(getVariantAlleleNormal())
			  .append(delimiter).append(getVariantAlleleTumor())
			  .append(delimiter).append(mFlankingNormal)
			  .append(delimiter).append(mFlankingTumor)
			  .append(delimiter).append(mCovgTotalNormal)
			  .append(delimiter).append(mCovgTotalTumor)
			  .append(delimiter).append(mCovgVarNormal)
			  .append(delimiter).append(mCovgVarTumor)
			  .append(delimiter).append(calcVAFNormal())
			  .append(delimiter).append(calcVAFTumor());
			
			if (mRsID == GenotypeUtils.RsID_Novel) {
				sb.append(delimiter).append(Script.NovelStr);
			} else {
				sb.append(delimiter).append(Utils.rsPrefix).append(mRsID);
			}
			
			sb.append(delimiter).append(getMutationType().getPrintName());
			sb.append(delimiter).append(mHugoSymbol);
			return sb;
		}
		
		// ====================================================================		
		public static void TestClusteringInputOneSite_Robust() {
			int numTrials = 10000000;
			
			ClusteringInputOneSite oneSiteInfo = new ClusteringInputOneSite();
			for (int trial = 0; trial < numTrials; trial++) {
				Chrom chrom = Chrom.getChrom((byte) NumberUtils.getRandomInteger(0, Chrom.values().length - 1));
				int position = NumberUtils.getRandomInteger(0, (int) BitSetUtils.getMask(28));
				oneSiteInfo.setChrom(chrom);
				oneSiteInfo.setPosition(position);
				
				//System.out.println(trial);
				CompareUtils.ensureTrue(oneSiteInfo.getChrom().equals(chrom), "Chrom not equal!");
				CompareUtils.ensureTrue(position == oneSiteInfo.getPosition(), "Position not equal!");
			}
		}
	}
}
