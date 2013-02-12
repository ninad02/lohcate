package lohcate.clustering;
import genomeEnums.Chrom;
import genomeEnums.VariantLocation;
import genomeUtils.ChromPositionTracker;

//import java.awt.geom.Ellipse2D.Double;
import java.io.BufferedWriter;
import java.io.File;
import java.io.PrintStream;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Date;
import java.util.EnumMap;

import javax.swing.text.html.HTMLDocument.HTMLReader.IsindexAction;

import kMeans.DataPoint;
import kMeans.JCA;

import lohcate.AllelicBiasTable;
import lohcate.LOHcateSimulator;
import lohcate.Script;
import lohcate.LOHcateSimulator.LOHcateSimulatorGoldStandard;
import lohcate.LOHcateSimulator.LOHcateSimulatorParams;
import lohcateEnums.ClusterType;
import lohcateEnums.SeqPlatform;

import nutils.ArgumentParserUtils;
import nutils.ArrayUtils;
import nutils.ArrayUtils.ParallelArrayDouble;
import nutils.CompareUtils;
import nutils.ContingencyTable;
import nutils.IOUtils;
import nutils.NullaryClassFactory;
import nutils.PrimitiveWrapper;
import nutils.ContingencyTable.ContingencyTableValue;
import nutils.StringUtils;
import nutils.StringUtils.FileExtensionAndDelimiter;
import nutils.counter.BucketCounterEnum;
import nutils.counter.DynamicBucketCounter;

import org.apache.commons.math3.distribution.BinomialDistribution;
import org.apache.commons.math3.distribution.PoissonDistribution;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.jfree.data.xy.DefaultXYDataset;

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
	
	public static boolean correctAllelicBias() { return !ClusteringParams.GlobalClusteringParams.mIgnoreAllelicBias.getValue(); }			

	// ========================================================================
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
		//int indexFirstSomaticRowInAllVariants = allVariantRowsStr.size();  // save the size with just the germline variants added
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
		boolean pointsIndependent = false;
		metaData.mFDRTumor  = NumberUtils.getFDR_BenjaminiHochberg(metaData.mImbalancePValuesTumor,  ClusteringParams.GlobalClusteringParams.mFDRAlpha.getValue(), pointsIndependent);
		metaData.mFDRNormal = NumberUtils.getFDR_BenjaminiHochberg(metaData.mImbalancePValuesNormal, ClusteringParams.GlobalClusteringParams.mFDRAlpha.getValue(), pointsIndependent);
	}

	// ========================================================================
	private static Boolean classifySites_addToContingencyTables(ClusterType eventTest, ClusterType eventTruth, EnumMap<ClusterType, ContingencyTable> eventTablesMap) {		
		
		if (CompareUtils.isNull(eventTest) || CompareUtils.isNull(eventTruth)) {
			return null;
		}
		
		if (eventTest == ClusterType.Null) {
			return Boolean.FALSE;
		}

		for (ClusterType ct : ClusterType.values()) {
			ContingencyTable table = eventTablesMap.get(ct);
			
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
	private static Boolean classifySites_addToEventCountTables(ClusterType eventTest, BucketCounterEnum<ClusterType> counter) {
		if (CompareUtils.isNull(eventTest)) return null;
		
		if (eventTest == ClusterType.Null) return Boolean.FALSE;
		
		counter.increment(eventTest);
		return Boolean.TRUE;
	}
	
	// ========================================================================
	/**
	 * Curate SNP calls by clustering data points into HET ball, DUP wedge, LOH sidelobes, &c., and by grabbing dbsnp population allele frequencies when possible
	 * @param inDir naf-taf-inputs
	 * @param opt 0::Illumina, 1::SOLiD
	 */
	public static void classifySites(String inDir, String allelicBiasInFile, String sitesClassifiedDir, String vafComparisonPlotDir, String vafWaterfallPlotDir, String copyNumberPlotDir, SeqPlatform platform) {
		File[] files = (new File(inDir)).listFiles();
		
		// TODO -- make column names static constants
		System.out.println("Reading Allelic Bias file...");
		AllelicBiasTable allelicBiasTable = correctAllelicBias() ? AllelicBiasTable.readFileAndConstructTable(allelicBiasInFile, 3, 4) : null;		
		System.out.println("Finished Reading Allelic Bias file...");
		
		// Create output directory
		classifySitesHelper_MakeSubDirs(sitesClassifiedDir, vafComparisonPlotDir, vafWaterfallPlotDir, copyNumberPlotDir);

		//0.486486	0.532609	0.424775	0.301206
	
		EnumMap<Chrom, DynamicBucketCounter> dbcByChromLOH = DynamicBucketCounter.ClassFactory.newEnumMap(Chrom.class);		
		LOHcateSimulator.LOHcateSimulatorParams simulatorParams = null; //new LOHcateSimulator.LOHcateSimulatorParams();
		
		int countAtPositionMax = -1; 
		int fileIndex = 0;
		for (File file : files) {			
			int indexOfSubstring = file.getName().indexOf(Script.GermlineSuffix);
			if (indexOfSubstring >= 0) {
				String filename = file.getName();
				String sampleNameRoot = filename.substring(0, indexOfSubstring);  	
				String extension = filename.substring(filename.lastIndexOf(Utils.DotStr), filename.length());
				System.out.println("Processing (" + ++fileIndex + "): " + file.getName());
				
				classifySitesOneSample(file, sampleNameRoot, extension, allelicBiasTable, dbcByChromLOH, simulatorParams, sitesClassifiedDir, vafComparisonPlotDir, vafWaterfallPlotDir, copyNumberPlotDir, platform);

			}
		}
		
		// Now plot the LOH recurrence across samples
		ClusteringPlotting.plotRecurrenceGenomeWide(dbcByChromLOH, copyNumberPlotDir, ClusterType.LOH);

	}

	// ========================================================================
	public static void classifySitesOneSample(File file, String sampleNameRoot, String sampleFilenameExtension, 
			AllelicBiasTable allelicBiasTable,
			EnumMap<Chrom, DynamicBucketCounter> dbcByChromLOH, 
			LOHcateSimulatorParams simulatorParams,
			String sitesClassifiedDir, String vafComparisonPlotDir, String vafWaterfallPlotDir, String copyNumberPlotDir, 
			SeqPlatform platform) {
		
		StringUtils.FileExtensionAndDelimiter fileExtDelim = StringUtils.FileExtensionTSV;	
		StringBuilder sb = new StringBuilder(8192);
		
		// Get the file and data
		File inFile = file;
		if (!inFile.exists()) {
			CompareUtils.ensureTrue(false, "File does not exist: " + inFile.getAbsolutePath());
		}		
		ClusteringInputOneSample oneSampleData = readLinesFromFiles(inFile);

		String[] columnHeaders = new String[] {
			"chr", "coordinate", "refBase", "varBase", "variantBase-N", "variantBase-T", "refEnv-N", "refEnv-T", "Q20_TotCov_N", "Q20_TotCov_T", "Q20_VarCov_N", "Q20_VarCov_T", "Q20_VariantRatio_N", "Q20_VariantRatio_T", "dbsnp", "MutationType", "Hugo_Symbol", 
			"Event", "TumorNormalCopyNumRatio", "TumorCopyNum", "Adjusted_Q20_VariantRatio_T", "Adjusted_Q20_VariantRatio_N", "pVal_Imbalance_Tumor", "pVal_Imbalance_Normal" 
		};
		String headerStr = StringUtils.constructColumnDelimitedString(columnHeaders, fileExtDelim.mDelimiter, sb, true).toString();
		
		// Initialize output file
		String outFilename = sampleNameRoot + ".withCopyNum" + sampleFilenameExtension;
		String outFilenameFullPath = sitesClassifiedDir + File.separator + outFilename;
		PrintWriter out = new PrintWriter(IOUtils.getBufferedWriter(outFilenameFullPath));
		out.println(headerStr);
		
		// Initialize output file for simulation
		boolean isSimulation = simulatorParams != null;
		String outSimFilename = sampleNameRoot + ".simStats" + sampleFilenameExtension;
		String outSimFilenameFullPath = sitesClassifiedDir + File.separator + outSimFilename;
		PrintWriter outSim = isSimulation ? new PrintWriter(IOUtils.getBufferedWriter(outSimFilenameFullPath)) : null;
		
		// Initialize tables for simulation-evaluation tallying
		EnumMap<ClusterType, ContingencyTable> eventTables = ContingencyTable.ClassFactory.newEnumMap(ClusterType.class);
		EnumMap<ClusterType, BucketCounterEnum<ClusterType>> testCounts = new EnumMap<ClusterType, BucketCounterEnum<ClusterType>>(ClusterType.class);
		for (ClusterType ct : ClusterType.values()) {
			testCounts.put(ct, new BucketCounterEnum<ClusterType>(ClusterType.class));
		}

		//System.out.println("Num Event Tables: " + eventTables.size());
		ChromPositionTracker chromPosTracker = new ChromPositionTracker();
		
		int numLoopIterations = isSimulation ? simulatorParams.getNumIterations() : 1;
		for (int loopIter = 0; loopIter < numLoopIterations; loopIter++) {
			
			LOHcateSimulator simulator = new LOHcateSimulator();
			LOHcateSimulator.LOHcateSimulatorGoldStandard goldStandard = new LOHcateSimulator.LOHcateSimulatorGoldStandard(oneSampleData.getNumSites());
			if (isSimulation) {
				simulator.generateSimulatedDataForSample(simulatorParams, oneSampleData, goldStandard);
			}

			// Declare meta-data and fill it
			ClusteringInputOneSampleMetaData metaData = new ClusteringInputOneSampleMetaData(oneSampleData.mInfoSites.size());
			classifySitesHelper_preprocessMetaData(oneSampleData, metaData, allelicBiasTable, platform);				

			System.out.printf("FDR Alpha:\t%g\nFDR Tumor p-value: \t%g\nFDR Normal p-value:\t%g\n", ClusteringParams.GlobalClusteringParams.mFDRAlpha.getValue(), metaData.mFDRTumor, metaData.mFDRNormal);
			System.out.println("\tInferred Copy Number and Allelic Imbalance Per Site...");					

			int indexFirstSomaticRowInAllVariants = 0;  // Set as invalid value for now
			ClusterType[] clusters = Clustering.getClusters_withPlane(oneSampleData.mInfoSites, metaData.mImbalancePValuesTumor, metaData.mImbalancePValuesNormal, metaData.mFDRTumor, metaData.mFDRNormal, metaData.mTumorCopyNumRatiosPerGene, metaData.mCopyNumRatioPerChromNormal, metaData.mAdjustedVAFNormal, metaData.mAdjustedVAFTumor, outFilenameFullPath, indexFirstSomaticRowInAllVariants, platform); //get cluster assignments (HET ball, LOH sidelobes, DUP wedge, &c.)
			System.out.println("\tGot clusters");

			// Now initialize the data structure needed to plot
			int[] clusterTypeCounts = ArrayUtils.getEnumTypeCounts(clusters, ClusterType.values().length);
			ParallelArrayDouble[] clusterCoordinates  = getCoordinateArraysPerClusterType(clusterTypeCounts);
			ParallelArrayDouble[] waterfallPlotTumor  = getCoordinateArraysPerClusterType(clusterTypeCounts);
			ParallelArrayDouble[] waterfallPlotNormal = getCoordinateArraysPerClusterType(clusterTypeCounts);
			ParallelArrayDouble[] copyNumPlot         = getCoordinateArraysPerClusterType(clusterTypeCounts);

			// Do the post-processing
			chromPosTracker.clear();
			
			DoubleArrayList chromBoundaryXValue = new DoubleArrayList();
			DoubleArrayList chromBoundaryYValue = new DoubleArrayList();

			int startingRowGermlineOrSomaticOrAll = 0;  // 0 because header line has been stripped away
			for (int row = startingRowGermlineOrSomaticOrAll; row < oneSampleData.mInfoSites.size(); row++) {
				ClusteringInputOneSite oneSiteInfo = oneSampleData.mInfoSites.get(row);
				ClusterType clusterType = clusters[row];

				// Print to output
				if (!isSimulation) {
					out.print(oneSiteInfo.printToString(sb, true, StringUtils.FileExtensionTSV.mDelimiter).toString());
					out.printf("\t%s\t", clusterType.name());
					metaData.printSiteInformation(out, row, true);
					out.flush();
				}
				
				ClusterType eventTruth = goldStandard.getEvent(row);
				classifySites_addToContingencyTables(clusterType, eventTruth, eventTables);
				classifySites_addToEventCountTables(clusterType, testCounts.get(eventTruth));

				boolean chromCrossed = chromPosTracker.chromCrossedWithCurrentCoordinates(oneSiteInfo.getChrom(), oneSiteInfo.getPosition());						
				for (double d = 0; chromCrossed && (d <= 5.0); d += 0.02) {
					chromBoundaryXValue.add(chromPosTracker.getPositionGenomeWide());
					chromBoundaryYValue.add(d);
				}

				// Add recurrence count
				if (clusterType.isLOH()) {
					dbcByChromLOH.get(oneSiteInfo.getChrom()).incrementCount(oneSiteInfo.getPosition());					
				}

				clusterCoordinates [clusterType.ordinal()].add( metaData.mAdjustedVAFTumor[row],         metaData.mAdjustedVAFNormal[row] );
				waterfallPlotTumor [clusterType.ordinal()].add( chromPosTracker.getPositionGenomeWide(), metaData.mAdjustedVAFTumor[row]);
				waterfallPlotNormal[clusterType.ordinal()].add( chromPosTracker.getPositionGenomeWide(), metaData.mAdjustedVAFNormal[row]);
				copyNumPlot        [clusterType.ordinal()].add( chromPosTracker.getPositionGenomeWide(), (metaData.mTumorCopyNumRatiosPerGene[row] * Script.DefaultDiploidCopyNumber));
			}									

			// Now let's create the datasets needed to
			double[][] boundaryArrays = ArrayUtils.combineTwoDynamicArraysIntoOneStatic(chromBoundaryXValue, chromBoundaryYValue);
			DefaultXYDataset xyDatasetVAFPlot             = classifySitesHelper_createAndFillXYData(clusterCoordinates, null);
			DefaultXYDataset xyDatasetWaterfallPlotTumor  = classifySitesHelper_createAndFillXYData(waterfallPlotTumor,  boundaryArrays);
			DefaultXYDataset xyDatasetWaterfallPlotNormal = classifySitesHelper_createAndFillXYData(waterfallPlotNormal, boundaryArrays);
			DefaultXYDataset xyDatasetCopyNumber          = classifySitesHelper_createAndFillXYData(copyNumPlot, boundaryArrays);

			ClusteringPlotting.plotVAFComparison(xyDatasetVAFPlot,            vafComparisonPlotDir + File.separator + sampleNameRoot + ".VAFComparison", sampleNameRoot);
			ClusteringPlotting.plotVAFGenomeWide(xyDatasetWaterfallPlotTumor,  vafWaterfallPlotDir + File.separator + sampleNameRoot + ".VAF_GenomeWide_Tumor",  sampleNameRoot, true);				
			ClusteringPlotting.plotVAFGenomeWide(xyDatasetWaterfallPlotNormal, vafWaterfallPlotDir + File.separator + sampleNameRoot + ".VAF_GenomeWide_Normal", sampleNameRoot, false);
			ClusteringPlotting.plotCopyNumGenomeWide(xyDatasetCopyNumber,        copyNumberPlotDir + File.separator + sampleNameRoot + ".CopyNumber_GenomeWide", sampleNameRoot);				
		}
		out.close();
		
		// Print contingency table statistics
		for (ClusterType eventType : ClusterType.values()) {
			ContingencyTable table = eventTables.get(eventType);
			System.out.printf("Event:\t%10s\tSensitivity:\t%g\tSpecificity:\t%g\n", eventType, table.getSensitivity(), table.getSpecificity());						
		}	
		
		// Print matrix breakdown statistics (distribtion of test events for each truth event)
		for (ClusterType eventTruth : ClusterType.values()) {			
			BucketCounterEnum<ClusterType> testCountsForTruth = testCounts.get(eventTruth);
			System.out.printf("Matrix:\t%12s", eventTruth.name());
			for (ClusterType eventTest : ClusterType.values()) {
				System.out.printf("\t%d", testCountsForTruth.getCount(eventTest));				
			}
			System.out.println("");
		}
	}

	
	
	// ========================================================================
	private static DefaultXYDataset classifySitesHelper_createAndFillXYData(ParallelArrayDouble[] coordinatesByEvent, double[][] boundaryArrays) {
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
	private static ParallelArrayDouble[] getCoordinateArraysPerClusterType(int[] clusterTypeCounts) {
		
		ParallelArrayDouble[] clusterCoordinates = new ParallelArrayDouble[clusterTypeCounts.length];
		
		for (int clusterTypeIndex = 0; clusterTypeIndex < clusterCoordinates.length; clusterTypeIndex++) {
			clusterCoordinates[clusterTypeIndex] = new ParallelArrayDouble(clusterTypeCounts[clusterTypeIndex]);
		}
		
		return clusterCoordinates;
	}

	// ========================================================================
	private static void calculateAdjustedVAFs(ArrayList<ClusteringInputOneSite> rows, AllelicBiasTable allelicBiasTable, 
											  float[] adjustedVAFNormal, float[] adjustedVAFTumor,
											  float[] copyNumRatioPerChromNormal,
											  SeqPlatform platform) {
		
		final float defaultVAFNormal = 0.50f;
			
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
					boolean isGermlineChromGain = copyNumRatioPerChromNormal[chrom.ordinal()] > ClusteringParams.GlobalClusteringParams.mGermlineTrisomyThreshold.getValue();
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
						int indexDistTumor  = NumberUtils.getIndexOfPoissonWithHighestProbability((int) Math.round(readCountAverageTumor),  pdTumor);
						int indexDistNormal = NumberUtils.getIndexOfPoissonWithHighestProbability((int) Math.round(readCountAverageNormal), pdNormal);
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
			
			readCountTalliesNormal.incrementCount(NumberUtils.roundToNearest5(oneSiteInfo.mCovgTotalNormal));			
			readCountTalliesTumor.incrementCount(NumberUtils.roundToNearest5(oneSiteInfo.mCovgTotalTumor));
					
		}
		
		// Calculate the genome-wide average read count in the normal
		int totalReadCountNormal = ArrayUtils.arraySum(totalReadCountPerChromNormal);
		int totalReadCountTumor  = ArrayUtils.arraySum(totalReadCountPerChromTumor);
		int totalNumSites        = ArrayUtils.arraySum(numSitesPerChrom);
		float avgCoverageNormal = (float) totalReadCountNormal / (float) totalNumSites;
		float avgCoverageTumor  = (float) totalReadCountTumor  / (float) totalNumSites;
		
		// Re-assign to highest count
		// TODO might change
		//avgCoverageNormal = (float) readCountTalliesNormal.getKeyWithMaxCount();
		//avgCoverageTumor  = (float) readCountTalliesTumor.getKeyWithMaxCount();
		
		coverageRatioTumorToNormal.mFloat = avgCoverageTumor / avgCoverageNormal;		 
		System.out.println("Tumor Normal Genome-wide ratio: " + coverageRatioTumorToNormal.mFloat);				
		System.out.printf("Average Read Count Normal: %g\n", avgCoverageNormal);
		System.out.printf("Average Read Count Tumor:  %g\n", avgCoverageTumor);
		
		for (Chrom chrom : Chrom.values()) {
			if (chrom.isInvalid()) continue; 
			
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
			if (copyNumRatioPerChromNormal[chrom.ordinal()] > ClusteringParams.GlobalClusteringParams.mGermlineTrisomyThreshold.getValue()) {
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
						returnClusters[row] = ClusterType.HETGermline; 

					} else if (assignedClusterID == clusterTypeIDsFromAlgorithm[ClusterType.GainSomatic.ordinal()]) {
						returnClusters[row] = ClusterType.GainSomatic; 

					} else if (assignedClusterID == clusterTypeIDsFromAlgorithm[ClusterType.GainGermline.ordinal()]) {
						returnClusters[row] = ClusterType.GainGermline; 
						
					} else if (assignedClusterID == clusterTypeIDsFromAlgorithm[ClusterType.cnLOH.ordinal()]) {
						returnClusters[row] = ClusterType.cnLOH; 
						
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
				boolean normalSigImbalanced = (imbalancePValuesNormal[row] <= fdrNormal) && (copyNumRatioPerChromNormal[chrom.ordinal()] > ClusteringParams.GlobalClusteringParams.mGermlineTrisomyThreshold.getValue());
				
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
		return (copyNum >= ClusteringParams.GlobalClusteringParams.mAmplificationThreshold.getValue()); 				
	}
	
	// ========================================================================
	private static boolean isCopyNumInDiploidRange(float copyNum) {
		return (ClusteringParams.GlobalClusteringParams.mDeletionThreshold.getValue() < copyNum) && 
			   (copyNum < ClusteringParams.GlobalClusteringParams.mAmplificationThreshold.getValue());		
	}
	
	// ========================================================================
	private static boolean isCopyNumRatioInDiploidRange(float copyNumRatio) {
		float copyNum = copyNumRatio * Script.DefaultDiploidCopyNumber;
		return isCopyNumInDiploidRange(copyNum);
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
	
	// ========================================================================		
	// ========================================================================
	public static ArrayList<CentroidPackage> performKMeansClustering(int numKMeansClusters, int numIter, ArrayList<DataPoint> dataPoints, double multiplier) {
		JCA jca = new JCA(numKMeansClusters, numIter, dataPoints);
	    jca.startAnalysis();
	    jca.sortClusters();
	    
	    ArrayList<CentroidPackage> centroids = new ArrayList<CentroidPackage>(); 
	    for (int k = numKMeansClusters - 1; k >= 0; k--) {
	    	CentroidPackage cp = new CentroidPackage();
	    	cp.mCentroidYValue = jca.getCluster(k).getCentroid().getCy() / multiplier;
	    	cp.mNumDataPoints  = jca.getCluster(k).getNumDataPoints();
	    	centroids.add(cp);
	    }
	    return centroids;
	}
	
	// ========================================================================
	public static class CentroidPackage {
		double mCentroidYValue;
		int mNumDataPoints;
	}
	
	// ========================================================================

}
