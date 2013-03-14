package lohcate.clustering;
import genomeEnums.Chrom;
import genomeEnums.TissueType;
import genomeEnums.VariantLocation;
import genomeUtils.ChromPositionTracker;
import genomeUtils.GenomeConstants;
import genomeUtils.ObjectWalkerTracker;
import genomeUtils.SNVMap;

import java.io.File;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Date;
import java.util.ListIterator;

import kMeans.DataPoint;
import kMeans.JCA;

import lohcate.AllelicBiasTable;
import lohcate.CopyNumberRegionRange;
import lohcate.CopyNumberRegionsByChromosome;
import lohcate.LOHcateSimulator;
import lohcate.Script;
import lohcate.LOHcateSimulator.LOHcateSimulatorParams;
import lohcateEnums.ClusterType;
import lohcateEnums.SeqPlatform;

import nutils.ArrayUtils;
import nutils.ArrayUtils.ParallelArrayDouble;
import nutils.CompareUtils;
import nutils.ContingencyTable;
import nutils.EnumMapSafe;
import nutils.IOUtils;
import nutils.PrimitiveWrapper;
import nutils.ContingencyTable.ContingencyTableValue;
import nutils.RangeDouble;
import nutils.StringUtils;
import nutils.StringUtils.FileExtensionAndDelimiter;
import nutils.counter.BucketCounterEnum;
import nutils.counter.DynamicBucketCounter;
import nutils.math.PoissonDistributionList;

import org.apache.commons.math3.distribution.PoissonDistribution;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.jfree.data.xy.DefaultXYDataset;

import com.carrotsearch.hppc.DoubleArrayList;
import com.carrotsearch.hppc.IntArrayList;
import nutils.NumberUtils;


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
	
	public static final int AllelicBiasTable_Col_NumSamples   = 3;
	public static final int AllelicBiasTable_Col_AvgVAFNormal = 4;
	
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
	/** Returns the number of duplicate sites. */
	private static int removeDuplicateSites(ArrayList<String> rows) {
		int positionPrev = 0;
		Chrom chromPrev = null;
		int numDuplicateSites = 0;
		int lastIndex = rows.size() - 1;
		
		for (int i = lastIndex; i >= 0; i--) {
			String line = rows.get(i);
			
			final Chrom chrom  = Chrom.getChrom(   StringUtils.extractNthColumnValue(line, Script.Col_NAFTAFInput_Chrom,     StringUtils.FileExtensionTSV.mDelimiter) );
			final int position = Integer.parseInt( StringUtils.extractNthColumnValue(line, Script.Col_NAFTAFInput_Position,  StringUtils.FileExtensionTSV.mDelimiter) );

			if (i < lastIndex) {
				if (chrom == chromPrev && position == positionPrev) {
					numDuplicateSites++;
					rows.set(i + 1, null);
				}
			}
			chromPrev = chrom;
			positionPrev = position;				
		}
		
		ArrayUtils.removeNullElements(rows);
		return numDuplicateSites;
	}

	// ========================================================================
	private static ClusteringInputOneSample readLinesFromFiles(File file) {	
		ArrayList<String> allVariantRowsStr = readLinesFromFilesAsStringList(file, true, Script.EliminateHighDensitySNVs, Script.EliminateExtremeGCSites);
		ClusteringInputOneSample oneSampleData = new ClusteringInputOneSample(allVariantRowsStr);				
		allVariantRowsStr.clear();
		return oneSampleData;
	}
	
	// ========================================================================	
	private static ArrayList<String> readLinesFromFilesAsStringList(File file, boolean removeDuplicateRows, boolean removeHighDensitySNVSites, boolean removeExtremeGCSites) {
		String somaticFilename = file.getAbsolutePath().replace(VariantLocation.Germline.toLowerCase(), VariantLocation.Somatic.toLowerCase());
		ArrayList<String> somaticSpecificVariantRows  = IOUtils.readAllLinesFromFile(somaticFilename);
		ArrayList<String> germlineSpecificVariantRows = IOUtils.readAllLinesFromFile(file.getAbsolutePath());
		System.out.println("\tRead All Lines...");
		System.out.printf("\tNum Sites Total: %d\n", germlineSpecificVariantRows.size() + somaticSpecificVariantRows.size());				
		
		//String headerStringGermline = germlineSpecificVariantRows.get(0);
		//String headerStringSomatic  = somaticSpecificVariantRows.get(0);

		//upstream pipelines will randomly spit header lines into the 
		// middle of a naf-taf-input file. we're just avoiding those
		germlineSpecificVariantRows = Script.curateSNPCalls_removeHeaderLinesFromRows(germlineSpecificVariantRows, removeHighDensitySNVSites, removeExtremeGCSites);
		somaticSpecificVariantRows  = Script.curateSNPCalls_removeHeaderLinesFromRows(somaticSpecificVariantRows,  removeHighDensitySNVSites, removeExtremeGCSites);
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
		
		if (removeDuplicateRows) {
			int numDuplicateRows = removeDuplicateSites(allVariantRowsStr);
			System.out.printf("\tRemoved %d duplicate Rows.  Num Rows left: %d\n", numDuplicateRows, allVariantRowsStr.size());
		}
		
		return allVariantRowsStr;
	}

	// ========================================================================
	// ========================================================================
	/** Performs pre-processing (copy number calling, allelic imbalace p-values) on the input data
	 *  and stores the results in the metadata
	 */
	private static void classifySitesHelper_preprocessMetaData(ClusteringInputOneSample oneSampleData, ClusteringInputOneSampleMetaData metaData, AllelicBiasTable allelicBiasTable, SeqPlatform platform) {
		
		// Get the allele frequency statistics, and adjust the frames based on the resulting standard deviation
		if (correctAllelicBias()) {
			metaData.mVAFNormalHetRange.setBoundLower(AlleleFractionStatsForSample.VAFNormalFrameLower);
			metaData.mVAFNormalHetRange.setBoundUpper(AlleleFractionStatsForSample.VAFNormalFrameUpper);
		} else {
			AlleleFractionStatsForSample afStatsSample = new AlleleFractionStatsForSample();
			afStatsSample.tabulateAndPerformStatistics(oneSampleData.mInfoSites, platform);
			metaData.mVAFNormalHetRange.setBoundLower( afStatsSample.getValueNStandardDeviationsAway(-NAF_STRIP_EXPANDER) );
			metaData.mVAFNormalHetRange.setBoundUpper( afStatsSample.getValueNStandardDeviationsAway( NAF_STRIP_EXPANDER) );
		}		
		
		// First, get the copy number ratios at a per-chromosome level				
		float[] tumorNormalCopyNumRatiosPerChrom = calcAvgCoverageRatioPerChrom(oneSampleData, metaData);
						
		// Now get the copy number ratios at a segmented sub-chromosomal level (by gene)
		calcRoughCopyNumberRatioPerSite(oneSampleData.mInfoSites, metaData);	
		
		// Now smooth the copy numbers
		smoothCopyNumbers(oneSampleData.mInfoSites, metaData);
		
		// Now adjust the VAF values basead on biases
		calculateAdjustedVAFs(oneSampleData.mInfoSites, metaData, allelicBiasTable, platform);

		// Calculate the imbalance p-values
		getPValuesImbalance(oneSampleData.mInfoSites, metaData);
		
		// Calculate the FDR for tumor and normal
		boolean pointsIndependent = false;
		double[] tumorPValWithoutGermline = getSitesPvaluesAfterElimintingGermlineGain(oneSampleData, metaData, metaData.mVAFNormalHetRange);
		//metaData.mFDRTumor  = NumberUtils.getFDR_BenjaminiHochberg(metaData.mImbalancePValuesTumor,  ClusteringParams.GlobalClusteringParams.mFDRAlpha.getValue(), pointsIndependent);
		metaData.mFDRTumor  = NumberUtils.getFDR_BenjaminiHochberg(tumorPValWithoutGermline,  ClusteringParams.GlobalClusteringParams.mFDRAlpha.getValue(), pointsIndependent);
		metaData.mFDRNormal = NumberUtils.getFDR_BenjaminiHochberg(metaData.mImbalancePValuesNormal, ClusteringParams.GlobalClusteringParams.mFDRAlpha.getValue(), pointsIndependent);
	}

	// ========================================================================
	// Conve
	private static double[] getSitesPvaluesAfterElimintingGermlineGain(ClusteringInputOneSample oneSampleData, ClusteringInputOneSampleMetaData metaData, RangeDouble vafNormalRange) {		
		int numSites = oneSampleData.getNumSites();
		DoubleArrayList pValues = new DoubleArrayList(numSites);
		
		for (int row = 0; row < numSites; row++) {
			ClusteringInputOneSite oneSiteInfo = oneSampleData.getSiteAtIndex(row); 
			
			if (!metaData.chromHasGermlineGain(oneSiteInfo.getChrom()) && vafNormalRange.inRangeLowerExclusive(oneSiteInfo.calcVAFNormal())) {
				pValues.add(metaData.mImbalancePValuesTumor[row]);
			}
		}
		System.out.printf("Shorted Pvalue List:\t%d\t%d\n", numSites, pValues.size());
		return pValues.toArray();
	}
	
	// ========================================================================
	private static Boolean classifySites_addToContingencyTables(ClusterType eventTest, ClusterType eventTruth, EnumMapSafe<ClusterType, ContingencyTable> eventTablesMap) {		
		
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
	public static AllelicBiasTable constrcutAllelicBiasTable(ArrayList<File> files, AllelicBiasTable biasTable, FileExtensionAndDelimiter fileExtAndDelim, SeqPlatform platform) {
		biasTable = CompareUtils.isNull(biasTable) ? new AllelicBiasTable() : biasTable;
		for (File inFile : files) {
			ArrayList<String> allLines = readLinesFromFilesAsStringList(inFile, false, false, false);
			for (String line : allLines) { 
				final Chrom chrom  = Chrom.getChrom(   StringUtils.extractNthColumnValue(line, Script.Col_NAFTAFInput_Chrom,    fileExtAndDelim.mDelimiter) );
				final int position = Integer.parseInt( StringUtils.extractNthColumnValue(line, Script.Col_NAFTAFInput_Position, fileExtAndDelim.mDelimiter) );
				final float vafNormal = extractVAFNormal(line, platform);
				
				if (AlleleFractionStatsForSample.VAFNormalRange.inRangeBothExclusive(vafNormal)) {
					biasTable.registerSite(chrom, position, vafNormal);
				}
			}
		}
		return biasTable;
	}
	
	// ========================================================================
	/**
	 * Curate SNP calls by clustering data points into HET ball, DUP wedge, LOH sidelobes, &c., and by grabbing dbsnp population allele frequencies when possible
	 * @param inDir naf-taf-inputs
	 * @param opt 0::Illumina, 1::SOLiD
	 */
	public static void classifySites(String inDir, String allelicBiasInFile, String sitesClassifiedDir, String vafComparisonPlotDir, String vafWaterfallPlotDir, String copyNumberPlotDir, SeqPlatform platform) {
		File[] files = (new File(inDir)).listFiles();
		
		// Get the list of valid files
		ArrayList<File> validFiles = new ArrayList<File>(files.length);		
		for (File file : files) {
			int indexOfSubstring = file.getName().indexOf(Script.GermlineSuffix);
			if (indexOfSubstring >= 0) {
				validFiles.add(file);
			}
		}
		
		// We now handle the construction or reading of the allelic bias table.  If we are to correct
		// allelic biases, we either construct it from the input files themselves, or we read pre-calculated
		// allelic biases from an independent file.
		AllelicBiasTable allelicBiasTable = null;
		if (correctAllelicBias()) {
			if (allelicBiasInFile == null) {
				System.out.println((new Date()).toString() + " Calculating Allelic Biases...");
				allelicBiasTable = constrcutAllelicBiasTable(validFiles, null, StringUtils.FileExtensionTSV, platform);
				System.out.println((new Date()).toString() + " Finished Calculating Allelic Biases...");
			} else {
				System.out.println((new Date()).toString() + " Reading Allelic Bias file...");
				allelicBiasTable = AllelicBiasTable.readFileAndConstructTable(allelicBiasInFile, AllelicBiasTable_Col_NumSamples, AllelicBiasTable_Col_AvgVAFNormal);
				System.out.println((new Date()).toString() + " Finished Reading Allelic Bias file...");
			}
		}
		
		// Create output directory
		classifySitesHelper_MakeSubDirs(sitesClassifiedDir, vafComparisonPlotDir, vafWaterfallPlotDir, copyNumberPlotDir);
	
		EnumMapSafe<Chrom, DynamicBucketCounter> dbcByChromLOH = DynamicBucketCounter.ClassFactory.newEnumMap(Chrom.class);		
		//LOHcateSimulator.LOHcateSimulatorParams simulatorParams = new LOHcateSimulator.LOHcateSimulatorParams();
		LOHcateSimulator.LOHcateSimulatorParams simulatorParams = null;
		
		int fileIndex = 0;		
		for (File file : validFiles) {			
			String filename = file.getName();
			String sampleNameRoot = filename.substring(0, file.getName().indexOf(Script.GermlineSuffix));  	
			String extension = filename.substring(filename.lastIndexOf(StringUtils.DotStr), filename.length());
			System.out.println("Processing (" + ++fileIndex + " / " + validFiles.size() + "): " + file.getName());				
			classifySitesOneSample(file, sampleNameRoot, extension, allelicBiasTable, dbcByChromLOH, simulatorParams, sitesClassifiedDir, vafComparisonPlotDir, vafWaterfallPlotDir, copyNumberPlotDir, platform);
		}
		
		// Now plot the LOH recurrence across samples
		ClusteringPlotting.plotRecurrenceGenomeWide(dbcByChromLOH, copyNumberPlotDir, ClusterType.LOH);

	}

	// ========================================================================
	public static void classifySitesOneSample(File file, String sampleNameRoot, String sampleFilenameExtension, 
			AllelicBiasTable allelicBiasTable,
			EnumMapSafe<Chrom, DynamicBucketCounter> dbcByChromLOH, 
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
		EnumMapSafe<ClusterType, ContingencyTable> eventTables = ContingencyTable.ClassFactory.newEnumMap(ClusterType.class);
		EnumMapSafe<ClusterType, BucketCounterEnum<ClusterType>> testCounts = new EnumMapSafe<ClusterType, BucketCounterEnum<ClusterType>>(ClusterType.class);
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
			ClusteringInputOneSampleMetaData metaData = new ClusteringInputOneSampleMetaData(oneSampleData.mInfoSites.size(), sampleNameRoot);
			classifySitesHelper_preprocessMetaData(oneSampleData, metaData, allelicBiasTable, platform);				

			System.out.printf("FDR Alpha:\t%g\nFDR Tumor p-value: \t%g\nFDR Normal p-value:\t%g\n", ClusteringParams.GlobalClusteringParams.mFDRAlpha.getValue(), metaData.mFDRTumor, metaData.mFDRNormal);
			System.out.println("\tInferred Copy Number and Allelic Imbalance Per Site...");					

			// KEY STEP: Obtain the clusters
			int indexFirstSomaticRowInAllVariants = 0;  // Set as invalid value for now
			ArrayList<ClusterType> events = Clustering.getClusters_withPlane(oneSampleData.mInfoSites, metaData, metaData.mVAFNormalHetRange, outFilenameFullPath, indexFirstSomaticRowInAllVariants, platform); //get cluster assignments (HET ball, LOH sidelobes, DUP wedge, &c.)
			System.out.println("\tGot clusters");

			// BEGIN POST-PROCESSING
			boolean doOverride = true;
			if (isSimulation || doOverride) {
				SNVMap snvMap = new SNVMap();
				CopyNumberRegionsByChromosome regionsByChrom = Script.segmentRegionsOneSample(oneSampleData, events, null, snvMap);						
			
				EnumMapSafe<ClusterType, CopyNumberRegionsByChromosome> regionsInSamplePerEventType = 
						new EnumMapSafe<ClusterType, CopyNumberRegionsByChromosome>(ClusterType.class);

				for (ClusterType eventType : ClusterType.AmpLOHcnLOH) {
					CopyNumberRegionsByChromosome regionsInOneSampleMerged = Script.mergeRegionsWithConstraints(regionsByChrom, eventType, Script.REGION_SEGMENTATION_DIST_THRESHOLD);
					//regionsInOneSampleMerged.removeSingletonRegions();
					regionsInSamplePerEventType.put(eventType, regionsInOneSampleMerged);
					//regionsInOneSampleMerged.print(System.out, StringUtils.FileExtensionTSV.mDelimiter);
					
					for (Chrom chrom : Chrom.values()) {
						ArrayList<CopyNumberRegionRange> regionsOnChrom = regionsInOneSampleMerged.getRegions(chrom);
						for (CopyNumberRegionRange region : regionsOnChrom) {
							int mostLikelyList = 0;
							int indexStart = oneSampleData.getIndex(chrom, region.getRangeStart());
							CompareUtils.ensureTrue(indexStart >= 0, "ERROR: Starting index must be > 0");
							int indexEnd   = oneSampleData.getIndex(chrom, region.getRangeEnd());
							CompareUtils.ensureTrue(indexEnd >= indexStart, "ERROR: Ending index must be >= starting index!");
							
							if (!region.spansOneSite()) {
								String listOfListOfMeansStr = "({0.5};{0.3333,0.6667})";
								ArrayList<double[]> listOfListOfMeans = ArrayUtils.getListOfDoubleListsFromStringForm(listOfListOfMeansStr, true);							
								mostLikelyList = determineCopyGainVAFMaxLikelihood(oneSampleData, metaData, region, listOfListOfMeans, TissueType.Tumor);								
							} 
							
							if (mostLikelyList == 0) {
								for (int row = indexStart; row <= indexEnd; row++) {									
									ClusterType event = events.get(row);
									if (event == eventType) {
										events.set(row, ClusterType.HETGermline);
									}
								}
							
							}
						}
					}
					
				}			
				
				
				if (isSimulation) {
				for (int row = 0; row < oneSampleData.mInfoSites.size(); row++) {
					ClusteringInputOneSite oneSiteInfo = oneSampleData.mInfoSites.get(row);		
					final ClusterType event = isSimulation ? goldStandard.getEvent(row) : events.get(row);

					CopyNumberRegionsByChromosome regionsForEvent = regionsInSamplePerEventType.get(event);
					if (CompareUtils.isNotNull(regionsForEvent)) {
						CopyNumberRegionRange region = regionsForEvent.getRegion(oneSiteInfo.getChrom(), oneSiteInfo.getPosition());
						if (CompareUtils.isNotNull(region)) {
							CompareUtils.ensureTrue(region.mCopyNumberClusterType == event, "ERROR: Internal region parsing error!");
							if (region.spansOneSite()) {
								events.set(row, ClusterType.HETGermline);
							} else {
								events.set(row, region.mCopyNumberClusterType);
							}							
						}
					}
				}
				}
			}
			
			// Now initialize the data structure needed to plot
			int[] clusterTypeCounts = ArrayUtils.getEnumTypeCounts(events, ClusterType.class);					
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
				ClusterType eventType = events.get(row);

				// Print to output
				if (numLoopIterations == 1) {
					out.print(oneSiteInfo.printToString(sb, true, StringUtils.FileExtensionTSV.mDelimiter).toString());
					out.printf("\t%s\t", eventType.name());
					metaData.printSiteInformation(out, row, true);
					out.flush();
				}
								
				ClusterType eventTruth = goldStandard.getEvent(row);
				classifySites_addToContingencyTables(eventType, eventTruth, eventTables);
				classifySites_addToEventCountTables(eventType, testCounts.get(eventTruth));

				boolean chromCrossed = chromPosTracker.chromCrossedWithCurrentCoordinates(oneSiteInfo.getChrom(), oneSiteInfo.getPosition());						
				for (double d = 0; chromCrossed && (d <= 5.0); d += 0.02) {
					chromBoundaryXValue.add(chromPosTracker.getPositionGenomeWide());
					chromBoundaryYValue.add(d);
				}

				// Add recurrence count
				if (eventType.isLOH()) {
					dbcByChromLOH.get(oneSiteInfo.getChrom()).incrementCount(oneSiteInfo.getPosition());					
				}

				clusterCoordinates [eventType.ordinal()].add( metaData.mAdjustedVAFTumor[row],         metaData.mAdjustedVAFNormal[row] );
				waterfallPlotTumor [eventType.ordinal()].add( chromPosTracker.getPositionGenomeWide(), metaData.mAdjustedVAFTumor[row]);
				waterfallPlotNormal[eventType.ordinal()].add( chromPosTracker.getPositionGenomeWide(), metaData.mAdjustedVAFNormal[row]);
				copyNumPlot        [eventType.ordinal()].add( chromPosTracker.getPositionGenomeWide(), (metaData.mTumorCopyNumRatiosPerGene[row] * Script.DefaultDiploidCopyNumber));
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
			System.out.printf("Event:\t%12s\tSensitivity:\t%g\tSpecificity:\t%g\n", eventType, table.getSensitivity(), table.getSpecificity());						
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
		oneSampleData.clear();
	}

	// ========================================================================
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
	private static void calculateAdjustedVAFs(ArrayList<ClusteringInputOneSite> rows,
											  ClusteringInputOneSampleMetaData metaData,
											  AllelicBiasTable allelicBiasTable, 											  
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
				if (RangeDouble.inRangeLowerExclusive(vafNormal, AlleleFractionStatsForSample.VAFNormalFrameLower, AlleleFractionStatsForSample.VAFNormalFrameUpper)) {
					boolean isGermlineChromGain = metaData.chromHasGermlineGain(chrom); 							
					float avgVAFNormal = allelicBiasTable.getAvgVAF(chrom, position, ClusteringParams.GlobalClusteringParams.mAllelicBiasMinNumSamples.getValue());
					if (avgVAFNormal > 0) {
						// Site exists in table
						if (isGermlineChromGain) {
							if (Math.abs(vafNormal - ExpectedVAFNormalTrisomy) < Math.abs(vafNormal - (1.0 - ExpectedVAFNormalTrisomy))) {
								vafNormalExpected = ExpectedVAFNormalTrisomy;
							} else {
								vafNormalExpected = 1.0f - ExpectedVAFNormalTrisomy;
							}
						} 
						
						if (!isGermlineChromGain) {
							adjustmentFactor = vafNormalExpected / avgVAFNormal;
							float absDiff = Math.abs(vafNormalExpected - avgVAFNormal);
							offset = (vafNormal > vafNormalExpected) ? -absDiff : ((vafNormal == vafNormalExpected) ? 0 : absDiff);
						}
					}
				}
			}
						
			// Implement adjustments
			if (UseBidrectionalAdditiveOffset) {
				metaData.mAdjustedVAFNormal[row] = vafNormal + offset;							
				metaData.mAdjustedVAFTumor[row]  = vafTumor  + offset;				
			} else {
				metaData.mAdjustedVAFNormal[row] = vafNormal * adjustmentFactor;			
				metaData.mAdjustedVAFTumor[row]  = vafTumor  * adjustmentFactor;				
			}
			
			// Ensure within bounds
			metaData.mAdjustedVAFNormal[row] = Math.min(1.0f, Math.max(0, metaData.mAdjustedVAFNormal[row]));							
			metaData.mAdjustedVAFTumor[row]  = Math.min(1.0f, Math.max(0, metaData.mAdjustedVAFTumor[row]));
		}		
	}
	// ========================================================================
	public static void getPValuesImbalance(ArrayList<ClusteringInputOneSite> rows, ClusteringInputOneSampleMetaData metaData) { 
		
		for (int row = 0; row < rows.size(); row++) {
			ClusteringInputOneSite oneSiteInfo = rows.get(row);
												
			int varCovgTumor =  (int) Math.round(oneSiteInfo.mCovgTotalTumor  * metaData.mAdjustedVAFTumor[row]);
			int varCovgNormal = (int) Math.round(oneSiteInfo.mCovgTotalNormal * metaData.mAdjustedVAFNormal[row]);
			
			metaData.mImbalancePValuesTumor[row]  = getPValuesImbalanceTissue(oneSiteInfo.mCovgTotalTumor,  varCovgTumor);
			metaData.mImbalancePValuesNormal[row] = getPValuesImbalanceTissue(oneSiteInfo.mCovgTotalNormal, varCovgNormal);		
		}
	}

	// ------------------------------------------------------------------------
	private static double getPValuesImbalanceTissue(int coverageTotal, int coverageVariant) {
		int maxRefOrVarCovg = Math.max(coverageVariant, coverageTotal - coverageVariant);
		return nutils.NumberUtils.cumulativeProbabilitySuccess(coverageTotal, maxRefOrVarCovg, 0.5);
	}

	// ========================================================================
	private static void smoothCopyNumbers(ArrayList<ClusteringInputOneSite> rows, ClusteringInputOneSampleMetaData metaData) {
		
		ArrayList<CopyNumberRegionRange> setOfRegions = new ArrayList<CopyNumberRegionRange>();
		int positionDiffThreshold = 5_000_000;
		
		for (Chrom chrom : Chrom.values()) {
			if (chrom.isInvalid()) continue;
		
			setOfRegions.clear();
			ListIterator<CopyNumberRegionRange> regionsOnChrom = metaData.mGeneRegions.getIteratorForChromosome(chrom);
			
			while (regionsOnChrom.hasNext()) {
				CopyNumberRegionRange currentRegion = regionsOnChrom.next();
				
				if (!setOfRegions.isEmpty()) {
					int positionDiff = rows.get(currentRegion.getRangeStart()).getPosition() - rows.get(setOfRegions.get(0).getRangeStart()).getPosition();
					if (positionDiff > positionDiffThreshold) {
						//System.out.println("Num Regions: " + setOfRegions.size());
						double copyNumber = smoothCopyNumbersHelper_findAverage(setOfRegions);
						smoothCopyNumbersHelper_applyCopyNumberToSites(copyNumber, setOfRegions, metaData);					
						setOfRegions.clear();
					}
				}
				
				// Add this to the current set of regions
				setOfRegions.add(currentRegion);
			}
			
			// Now take care of any remaining regions
			double copyNumber = smoothCopyNumbersHelper_findAverage(setOfRegions);
			smoothCopyNumbersHelper_applyCopyNumberToSites(copyNumber, setOfRegions, metaData);					
			
		}		
	}

	// ========================================================================
	private static void smoothCopyNumbersHelper_applyCopyNumberToSites(double copyNumber, ArrayList<CopyNumberRegionRange> setOfRegions, ClusteringInputOneSampleMetaData metaData) {
		for (CopyNumberRegionRange cnrr : setOfRegions) {
			for (int index = cnrr.getRangeStart(); index <= cnrr.getRangeEnd(); index++) {
				metaData.mTumorCopyNumRatiosPerGene[index] = (float) copyNumber / GenomeConstants.DefaultDiploidCopyNumber;
			}
		}
	}
	
	// ========================================================================
	private static double smoothCopyNumbersHelper_findAverage(ArrayList<CopyNumberRegionRange> setOfRegions) {
		double sum = 0;
		for (CopyNumberRegionRange cnrr : setOfRegions) {
			sum += cnrr.mCopyNumber;
		}
		return sum / setOfRegions.size();
	}
	
	// ========================================================================
	public static float[] calcRoughCopyNumberRatioPerSite(ArrayList<ClusteringInputOneSite> rows, ClusteringInputOneSampleMetaData metaData) {		
		
		ObjectWalkerTracker<String> prevGene = new ObjectWalkerTracker<String>("", String.CASE_INSENSITIVE_ORDER);
		int rowOfFirstInstanceOfGene = -1;
		int positionPrev = 0;
		
		float ratioSum = 0;
		int readCountSumTumor = 0;
		int readCountSumNormal = 0;
		int numRowsWithSameGene = 0;
		CopyNumberRegionRange geneRegion = null;
		
		int avgCoverageNormal = metaData.mReadCountTalliesNormal.getKeyWithMaxCount();
		int avgCoverageTumor  = metaData.mReadCountTalliesTumor.getKeyWithMaxCount();
		double avgCoverageHaploidNormal = avgCoverageNormal / 2.0;
		double avgCoverageHaploidTumor  = avgCoverageTumor  / 2.0;
		
		int numPoissons = 10;
		PoissonDistribution[] pdNormal = new PoissonDistribution[numPoissons];
		PoissonDistribution[] pdTumor  = new PoissonDistribution[numPoissons];
		for (int i = 0; i < numPoissons; i++) {
			pdNormal[i] = new PoissonDistribution(avgCoverageHaploidNormal * (i + 1));
			pdTumor[i]  = new PoissonDistribution(avgCoverageHaploidTumor  * (i + 1));
		}		
		
		Arrays.fill(metaData.mTumorCopyNumRatiosPerGene, Script.DefaultTumorNormalRatio);  // Initialize
		
		for (int row = 0; row < rows.size(); row++) {
			ClusteringInputOneSite oneSiteInfo = rows.get(row);
			int positionDifference = oneSiteInfo.getPosition() - positionPrev;
			positionPrev = oneSiteInfo.getPosition();

			metaData.mIsSomaticSite[row] = (oneSiteInfo.mCovgTotalNormal <= 0);
			float ratio = metaData.mIsSomaticSite[row] ? Script.TumorNormalRatioOfSomaticSite : ((float) oneSiteInfo.mCovgTotalTumor / (float) oneSiteInfo.mCovgTotalNormal);						
			boolean validGeneName = !oneSiteInfo.mHugoSymbol.equals(StringUtils.DotStr); 
			
			if (!prevGene.hasChanged(oneSiteInfo.mHugoSymbol, true).hasChanged() 
				&& (geneRegion.getChromosome() == oneSiteInfo.getChrom()) 
				&& (validGeneName || (positionDifference <= 1000))) { 
				
				boolean extendResult = geneRegion.extendRange(oneSiteInfo.getChrom(), row);
				//System.out.println(geneRegion.toString() + "\t" + prevGene.getPrevious() + "\t" + oneSiteInfo.getChrom() + "\t" + row + "\t" + oneSiteInfo.mHugoSymbol);
				CompareUtils.ensureTrue(extendResult, "ERROR: Could not extend result!");	
				
				if (metaData.mIsSomaticSite[row]) {
					// We have a somatic site in the same gene as previous site, we take no action in this block
				} else {
					numRowsWithSameGene++;
					ratioSum += ratio;
					readCountSumTumor  += oneSiteInfo.mCovgTotalTumor;
					readCountSumNormal += oneSiteInfo.mCovgTotalNormal;
				}
				
			} else {
				// New gene listed.  We must now write the ratios for the previous gene
				if (row > 0) {					
					float averageRatio = Script.TumorNormalRatioOfSomaticSite;
					
					if (numRowsWithSameGene == 0) {
						CompareUtils.ensureTrue(metaData.mIsSomaticSite[row - 1], "ERROR: Gene must > 0 sites representative of that gene!");
						
					} else if (numRowsWithSameGene == 1 && !AllowGeneWithOneSite) {					
						averageRatio = Script.DefaultTumorNormalRatio;
						
					} else if (numRowsWithSameGene > 1 || AllowGeneWithOneSite) {
						double readCountAverageTumor  = (double) readCountSumTumor  / (double) numRowsWithSameGene;
						double readCountAverageNormal = (double) readCountSumNormal / (double) numRowsWithSameGene;
						int indexDistTumor  = NumberUtils.getIndexOfPoissonWithHighestProbability((int) Math.round(readCountAverageTumor),  pdTumor);
						int indexDistNormal = NumberUtils.getIndexOfPoissonWithHighestProbability((int) Math.round(readCountAverageNormal), pdNormal);
						averageRatio = (float) (indexDistTumor + 1) / (float) (indexDistNormal + 1);
						
						//averageRatio = ((ratioSum / (float) numRowsWithSameGene) / coverageRatioTumorNormalGenomeWide);
						averageRatio = ((float) readCountSumTumor / (float) readCountSumNormal) / metaData.mCoverageRatioTumorToNormal.mFloat;								
					}
	
					for (int i = rowOfFirstInstanceOfGene; i < row; i++) {					
						//metaData.mTumorCopyNumRatiosPerGene[i] = metaData.mIsSomaticSite[i] ? Script.TumorNormalRatioOfSomaticSite : averageRatio;
						metaData.mTumorCopyNumRatiosPerGene[i] = averageRatio;
						// This allows ratios to be assigned to somatic sites within the gene (that has other variant sites).
						// In other words, the somatic site borrows the copy 
					}
					
					// Now set the copy number for the region
					geneRegion.mCopyNumber = averageRatio * GenomeConstants.DefaultDiploidCopyNumber;
				}
				
				// New gene listed.  Thus, set a new row of first gene			
				rowOfFirstInstanceOfGene = row;				
				geneRegion = new CopyNumberRegionRange(ClusterType.Null, oneSiteInfo.getChrom(), row);
				geneRegion.mCopyNumber = GenomeConstants.DefaultDiploidCopyNumber;
				metaData.mGeneRegions.addRegion(geneRegion);				
				if (metaData.mIsSomaticSite[row]) {
					ratioSum = numRowsWithSameGene = readCountSumTumor = readCountSumNormal = 0;
					// We set the read counts to 0 because of the following case: If the next site
					// resides within the same gene as this site, then the copy number calculations
					// for the next site should not be affected by the read counts at this somatic
					// mutation site (as tumor reads will be non-zero, while normal reads are zero).
				} else {
					ratioSum = ratio;
					numRowsWithSameGene = 1;										 					
					readCountSumTumor  = oneSiteInfo.mCovgTotalTumor;
					readCountSumNormal = oneSiteInfo.mCovgTotalNormal;
				}
			}			
		}
		
		// Normalize to diploid
		DescriptiveStatistics ds = new DescriptiveStatistics();
		for (int row = 0; row < rows.size(); row++) {
			ds.addValue(metaData.mTumorCopyNumRatiosPerGene[row]);
		}
		double median = ds.getPercentile(50);
		System.out.println("Median: " + median);
		double diff = 1.0 - median;
		
		for (int row = 0; row < rows.size(); row++) {
			metaData.mTumorCopyNumRatiosPerGene[row] += (float) diff;
		}
		
		return metaData.mTumorCopyNumRatiosPerGene;
	}

	// ========================================================================
	/** Tumor : Normal avg coverage ratio per chromosome. */
	public static float[] calcAvgCoverageRatioPerChrom(ClusteringInputOneSample oneSampleInfo, ClusteringInputOneSampleMetaData metaData) {		
		
		float[] tumorNormalRatioPerChrom = new float[Chrom.values().length];
		int[] totalReadCountPerChromNormal = new int[Chrom.values().length];
		int[] totalReadCountPerChromTumor  = new int[Chrom.values().length];
		
		Arrays.fill(tumorNormalRatioPerChrom, Script.DefaultTumorNormalRatio);  // Fill with 1.0 since that's an equal ratio
		Arrays.fill(metaData.mCopyNumRatioPerChromNormal, Script.DefaultNormalRatio);
		
		Arrays.fill(metaData.mNumSitesPerChrom, 0);  // initialize counts
		Arrays.fill(metaData.mAvgReadCountPerChromNormal, 0); // initialize
		Arrays.fill(totalReadCountPerChromNormal, 0);
		Arrays.fill(totalReadCountPerChromTumor,  0);
		
		metaData.mReadCountTalliesTumor.clear();
		metaData.mReadCountTalliesNormal.clear();
		
		for (int row = 0; row < oneSampleInfo.mInfoSites.size(); row++) {
			ClusteringInputOneSite oneSiteInfo = oneSampleInfo.mInfoSites.get(row);
			
			int chromIndex = oneSiteInfo.getChrom().ordinal();
			totalReadCountPerChromNormal[chromIndex] += oneSiteInfo.mCovgTotalNormal;
			totalReadCountPerChromTumor[chromIndex]  += oneSiteInfo.mCovgTotalTumor;
			++metaData.mNumSitesPerChrom[chromIndex];
			
			//=ROUND(A1*2/10,0)*10/2
			
			metaData.mReadCountTalliesNormal.incrementCount(NumberUtils.roundToNearest5(oneSiteInfo.mCovgTotalNormal));			
			metaData.mReadCountTalliesTumor.incrementCount(NumberUtils.roundToNearest5(oneSiteInfo.mCovgTotalTumor));
					
		}
		
		// Calculate the genome-wide average read count in the normal
		int totalReadCountNormal = ArrayUtils.arraySum(totalReadCountPerChromNormal);
		int totalReadCountTumor  = ArrayUtils.arraySum(totalReadCountPerChromTumor);
		int totalNumSites        = ArrayUtils.arraySum(metaData.mNumSitesPerChrom);
		float avgCoverageNormal = (float) totalReadCountNormal / (float) totalNumSites;
		float avgCoverageTumor  = (float) totalReadCountTumor  / (float) totalNumSites;
		
		// Re-assign to highest count
		// TODO might change
		//avgCoverageNormal = (float) readCountTalliesNormal.getKeyWithMaxCount();
		//avgCoverageTumor  = (float) readCountTalliesTumor.getKeyWithMaxCount();
		
		metaData.mCoverageRatioTumorToNormal.mFloat = avgCoverageTumor / avgCoverageNormal;		 
		System.out.println("Tumor Normal Genome-wide ratio: " + metaData.mCoverageRatioTumorToNormal.mFloat);				
		System.out.printf("Average Read Count Normal: %g\n", avgCoverageNormal);
		System.out.printf("Average Read Count Tumor:  %g\n", avgCoverageTumor);
		
		for (Chrom chrom : Chrom.values()) {
			if (chrom.isInvalid()) continue; 
			
			int chromIndex = chrom.ordinal();
			
			if (totalReadCountPerChromNormal[chromIndex] > 0) {
				tumorNormalRatioPerChrom[chromIndex] = (float) totalReadCountPerChromTumor[chromIndex] / (float) totalReadCountPerChromNormal[chromIndex];
			}
			
			if (metaData.mNumSitesPerChrom[chromIndex] > 0) {
				metaData.mAvgReadCountPerChromNormal[chromIndex] = totalReadCountPerChromNormal[chromIndex] / metaData.mNumSitesPerChrom[chromIndex];
				metaData.mCopyNumRatioPerChromNormal[chromIndex] = metaData.mAvgReadCountPerChromNormal[chromIndex] / avgCoverageNormal;				
			}	
			
			System.out.printf("Chrom: %d\tNormal Ratio:%g\tTumor-Normal Ratio %g\tNum Sites: %d\n", chromIndex, metaData.mCopyNumRatioPerChromNormal[chromIndex], tumorNormalRatioPerChrom[chromIndex], metaData.mNumSitesPerChrom[chromIndex]);
		}
					
		determineGermlineAneuploidiesByVAF(oneSampleInfo, metaData);
		return tumorNormalRatioPerChrom;
	}	

	// ========================================================================
	public static int determineCopyGainVAFMaxLikelihood(ClusteringInputOneSample oneSampleInfo, ClusteringInputOneSampleMetaData metaData, CopyNumberRegionRange range, ArrayList<double[]> listOfListOfMeans, TissueType targetTissue) {
		Chrom chrom = range.getChromosome();
		if (chrom.isInvalid()) return -1;  // Return if not on a valid chromosome

		int indexChromStart = oneSampleInfo.getIndex(chrom, range.getRangeStart());
		int indexChromEnd   = oneSampleInfo.getIndex(chrom, range.getRangeEnd());		

		indexChromStart = Math.max(0, indexChromStart);
		indexChromEnd = Math.min(indexChromEnd, oneSampleInfo.getNumSites() - 1);
		
		// Construct our lists
		int multiplier = 100;
		ArrayList<PoissonDistributionList> pdListofLists = new ArrayList<PoissonDistributionList>();				
		for (double[] doubleArray : listOfListOfMeans) {
			PoissonDistributionList pdList = new PoissonDistributionList();
			pdListofLists.add(pdList);
			for (double d : doubleArray) {
				pdList.registerMean(multiplier * d);
			}			
		}
		
		double[] prob = new double[pdListofLists.size()];		
		Arrays.fill(prob, 0);
		for (int row = indexChromStart; row <= indexChromEnd; row++) {
			ClusteringInputOneSite oneSiteInfo = oneSampleInfo.mInfoSites.get(row);				

			// Get the relevant variant allele fraction
			double vafTissue = (targetTissue == TissueType.Normal) ? oneSiteInfo.calcVAFNormal() : oneSiteInfo.calcVAFTumor();
			
			// Skip homozygous or near-homozygous sites
			if (oneSiteInfo.refOrVarHasZeroReadCount()) continue;  // skip homozygous sites
			if (!metaData.mVAFNormalHetRange.inRangeLowerExclusive(vafTissue)) continue;  // skip almost homozygous sites

			
			PrimitiveWrapper.WDouble probVar = new PrimitiveWrapper.WDouble(0);
			for (int indexDist = 0; indexDist < pdListofLists.size(); indexDist++) {
				PoissonDistributionList pdList = pdListofLists.get(indexDist);
				pdList.getMeanOfMostLikelyDistribution((int) Math.round(multiplier * vafTissue), probVar);					
				prob[indexDist] += probVar.mDouble;
			}
		}

		return ArrayUtils.getIndexOfMaxElement(prob, 0);
	}
	
	
	// ========================================================================
	public static void determineGermlineAneuploidiesByVAF(ClusteringInputOneSample oneSampleInfo, ClusteringInputOneSampleMetaData metaData) {
		String listOfListOfMeansStr = "({0.5};{0.6667,0.3333})";
		ArrayList<double[]> listOfListOfMeans = ArrayUtils.getListOfDoubleListsFromStringForm(listOfListOfMeansStr, true);
		
		for (Chrom chrom : Chrom.values()) {
			if (chrom.isInvalid()) continue;  // move on if not on a valid chromosome
			
			int indexChromStart = oneSampleInfo.getIndexChromStart(chrom);
			if (indexChromStart < 0) continue;

			int indexChromEnd = oneSampleInfo.getIndexChromEnd(chrom);
			CompareUtils.ensureTrue(indexChromEnd >= 0, "ERROR: Chromosome must have last position!");
						
			CopyNumberRegionRange range = new CopyNumberRegionRange(ClusterType.Null, chrom, oneSampleInfo.getSiteAtIndex(indexChromStart).getPosition(), oneSampleInfo.getSiteAtIndex(indexChromEnd).getPosition());
			int indexOfMostLikelyList = determineCopyGainVAFMaxLikelihood(oneSampleInfo, metaData, range, listOfListOfMeans, TissueType.Normal);
			
			metaData.mChromHasGermlineGain[chrom.ordinal()] = (indexOfMostLikelyList > 0);					
			//System.out.printf("GermProb\tchr%d\t%g\t%g\t%b\n", chrom.ordinal(), prob[0], prob[1], metaData.mChromHasGermlineGain[chrom.ordinal()]);
		}
	}

	
	// ========================================================================
	public static void determineGermlineAneuploidies(ClusteringInputOneSample oneSampleInfo, ClusteringInputOneSampleMetaData metaData) {
		ArrayList<PoissonDistributionList> pdListofLists = new ArrayList<PoissonDistributionList>();
		double[] prob = new double[2];		
		
		for (Chrom chrom : Chrom.values()) {
			if (chrom.isInvalid()) continue;
			
			int indexChromStart = oneSampleInfo.getIndexChromStart(chrom);
			if (indexChromStart < 0) continue;
			
			int indexChromEnd = oneSampleInfo.getIndexChromEnd(chrom);
			CompareUtils.ensureTrue(indexChromEnd >= 0, "ERROR: Chromosome must have last position!");
			
			double avgReadCountOnChrom = metaData.mAvgReadCountPerChromNormal[chrom.ordinal()];
			pdListofLists.clear();

			PoissonDistributionList pdList1 = new PoissonDistributionList();
			pdList1.registerMean(avgReadCountOnChrom / GenomeConstants.DefaultDiploidCopyNumber);
			pdListofLists.add(pdList1);
			
			PoissonDistributionList pdList2 = new PoissonDistributionList();
			pdList2.registerMean(avgReadCountOnChrom * (2.0 / 3.0));
			pdList2.registerMean(avgReadCountOnChrom * (1.0 / 3.0));
			pdListofLists.add(pdList2);
						
			Arrays.fill(prob, 0);
			for (int row = indexChromStart; row <= indexChromEnd; row++) {
				ClusteringInputOneSite oneSiteInfo = oneSampleInfo.mInfoSites.get(row);
				int covgVar = oneSiteInfo.mCovgVarNormal;
				int covgRef = oneSiteInfo.mCovgTotalNormal - covgVar;
				
				if (oneSiteInfo.refOrVarHasZeroReadCount()) continue; // skip homozygous sites
				if (!metaData.mVAFNormalHetRange.inRangeLowerExclusive(oneSiteInfo.calcVAFNormal())) continue;  // skip almost homozygous sites
				
				PrimitiveWrapper.WDouble probRef = new PrimitiveWrapper.WDouble(0);
				PrimitiveWrapper.WDouble probVar = new PrimitiveWrapper.WDouble(0);
				
				for (int indexDist = 0; indexDist < pdListofLists.size(); indexDist++) {
					PoissonDistributionList pdList = pdListofLists.get(indexDist);
					pdList.getMeanOfMostLikelyDistribution(covgRef, probRef);				
					pdList.getMeanOfMostLikelyDistribution(covgVar, probVar);
					//prob[indexDist] += probRef.mDouble + probVar.mDouble;
					prob[indexDist] += probVar.mDouble;
				}

				//int maxIndex = ArrayUtils.getIndexOfMaxElement(prob, 0);				
			}
			
			metaData.mChromHasGermlineGain[chrom.ordinal()] = prob[1] > prob[0];
			System.out.println("Gain\t" + chrom + "\t" + metaData.mChromHasGermlineGain[chrom.ordinal()]);
		}
		
	}
	
	// ========================================================================
	public static ArrayList<ClusterType> getClusters_withPlane(
			ArrayList<ClusteringInputOneSite> rows, 
			ClusteringInputOneSampleMetaData metaData,
			RangeDouble vafNormalRange,
			String outFilenameFullPath, int startingRowGermlineOrSomaticOrAll, SeqPlatform platform) {

		//apply DBScan to points within NAF frame
		ArrayList<Floint> pointsUpperPlane = new ArrayList<Floint>();
		ArrayList<Floint> pointsLowerPlane = new ArrayList<Floint>();
		int[] indexMapToUpperLowerPlane    = new int[rows.size()];
		IntArrayList mapToRowsFromUpper = new IntArrayList(rows.size());
		IntArrayList mapToRowsFromLower = new IntArrayList(rows.size());

		getValidPointsForClustering(rows, pointsUpperPlane, pointsLowerPlane, indexMapToUpperLowerPlane, 
				                    mapToRowsFromUpper, mapToRowsFromLower, 
				                    metaData,				                    
				                    ScalingFactor, 				                     
				                    platform, 
				                    vafNormalRange);		

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
			float copyNum = metaData.mTumorCopyNumRatiosPerGene[indexIntoOriginalRows] * Script.DefaultDiploidCopyNumber;
			
			if (isCopyNumAmplified(copyNum) && clusterAssignmentsUpperPlane[ind] != DBSCAN2.ClusterIDOfNoise) {
				clusterAssignmentsUpperPlane[ind] = clusterTypeIDsFromAlgorithm[ClusterType.GainSomatic.ordinal()];
			} else if (isCopyNumInDiploidRange(copyNum) && clusterAssignmentsUpperPlane[ind] != DBSCAN2.ClusterIDOfNoise) {
				clusterAssignmentsUpperPlane[ind] = clusterTypeIDsFromAlgorithm[ClusterType.cnLOH.ordinal()];
			}
			
			ClusteringInputOneSite oneSiteInfo = rows.get(indexIntoOriginalRows);
			Chrom chrom = oneSiteInfo.getChrom();
			if (metaData.chromHasGermlineGain(chrom)) {
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
					int trueIndex = ArrayUtils.getInsertPoint(indexMapToUpperLowerPlane[row]);							
					clusterAssignmentsFinal.add(clusterAssignmentsLowerPlane[trueIndex]);
				}
			}
		}
		

		return Clustering.assignClusters(rows, platform, vafNormalRange, clusterAssignmentsFinal.mArray, clusterTypeIDsFromAlgorithm);
	}

	// ========================================================================
	// Convenience private function to break up caller function
	private static ArrayList<ClusterType> assignClusters(ArrayList<ClusteringInputOneSite> rows, 
			SeqPlatform platform,
			//int startingRowSomatic,												
			RangeDouble vafNormalRange,
			int[] clusterAssignments,
			int[] clusterTypeIDsFromAlgorithm) {

		int indexInClusterAssignments = -1;  // we need to keep a special index, since not all rows are used.
		ArrayList<ClusterType> returnEvents = ArrayUtils.addToCollection(new ArrayList<ClusterType>(rows.size()), ClusterType.Null, rows.size(), true);		

		for (int row = 0; row < rows.size(); row++) {

			ClusteringInputOneSite oneSiteInfo = rows.get(row);
			float vafNormal = oneSiteInfo.calcVAFNormal();
			boolean vafInRangeNormal = vafNormalRange.inRangeLowerExclusive(vafNormal); 					
			ClusterType eventAtSite = null;

			if (/*(row >= startingRowSomatic) &&*/ (!vafInRangeNormal)) {
				// The vafNormal is either very low (homozygous reference) or very high (homozygous common variant).
				// We do some very simple decision making now (which should be replaced by formal clustering later)
				// to partition the calls.
				float justBelowZero = -0.0001f;				
				float hetBoundaryLower = 0.3333f;
				float hetBoundaryUpper = 0.6667f;
				float vafTumor = oneSiteInfo.calcVAFTumor();
				if (RangeDouble.inRangeLowerExclusive(vafNormal, justBelowZero, vafNormalRange.getBoundLower())) {
					// We are equal to or below the lower frame boundary
					if (vafNormal <= 0.05 && vafTumor >= 0.10) {
						eventAtSite = ClusterType.HETSomatic;
					} else if (vafTumor <= hetBoundaryLower) {
						// Normal: AA, Tumor: AA [Thus homozygous reference in both, no events]
						eventAtSite = ClusterType.Null;
					} else if (vafTumor > hetBoundaryUpper) {
						// Normal: AA, Tumor: BB or CC [made by: AA -> AB or AC (somatic het mutation) -> B or C (LOH, loss of A)]
						eventAtSite = ClusterType.LOH;
					} else {
						// Normal: AA, Tumor: AB [made by: AA -> AB (somatic het mutation)
						eventAtSite = ClusterType.Null;
					}					
				} else if (RangeDouble.inRangeLowerExclusive(vafNormal, vafNormalRange.getBoundUpper(), Script.MaxVariantAlleleFrequency)) {
					// We are above the upper frame boundary
					if (vafNormal >= (1.0 - 0.05) && vafTumor >= 0.10 && vafTumor <= 0.90) {
						eventAtSite = ClusterType.HETSomatic;
					} else if (vafTumor <= hetBoundaryLower) {
						// Normal: BB, Tumor: AA [made by: BB -> AB (reverse somatic het mutation) -> A (LOH, loss of B)]
						eventAtSite = ClusterType.LOH;
					} else if (vafTumor > hetBoundaryUpper) {
						// Normal: BB, Tumor: BB or CC (ambiguous until we know exact variant for tumor)
						// TODO - Leave as Null for now, but will need to change later to resolve the
						// ambiguity mentioned above
						eventAtSite = ClusterType.Null;
					} else {
						// Normal: BB, Tumor: AB or CB [made by: BB -> AB (reverse somatic het mutation) or BB -> CB (somatic het mutation)
						eventAtSite = ClusterType.Null;
					}
				} else {
					CompareUtils.throwErrorAndExit("ERROR: Contradiction - variant allele frequency cannot be in and out of bounds simultanteously!" + vafNormal);
				}

			} else {	// Our vaf-normal is in range in somatic sites, or we're at a germline site regardless of vaf-normal value	
				if (vafInRangeNormal) {
					++indexInClusterAssignments;
					final int assignedClusterID = clusterAssignments[indexInClusterAssignments]; // we create a local variable for fast test modifications 

					if (assignedClusterID        == clusterTypeIDsFromAlgorithm[ClusterType.HETGermline.ordinal()]) {
						eventAtSite = ClusterType.HETGermline; 

					} else if (assignedClusterID == clusterTypeIDsFromAlgorithm[ClusterType.GainSomatic.ordinal()]) {
						eventAtSite = ClusterType.GainSomatic; 

					} else if (assignedClusterID == clusterTypeIDsFromAlgorithm[ClusterType.GainGermline.ordinal()]) {
						eventAtSite = ClusterType.GainGermline; 
						
					} else if (assignedClusterID == clusterTypeIDsFromAlgorithm[ClusterType.cnLOH.ordinal()]) {
						eventAtSite = ClusterType.cnLOH; 
						
					} else if (assignedClusterID == clusterTypeIDsFromAlgorithm[ClusterType.Null.ordinal()]) {
						eventAtSite = ClusterType.Null;  // we're on a diagonal

					} else if (assignedClusterID == clusterTypeIDsFromAlgorithm[ClusterType.Noise.ordinal()]) { 
						eventAtSite = ClusterType.Noise;

					} else { //anything not in the HET ball / DUP wedge is considered part of a LOH sidelobe
						//float vafTumor = extractVAFTumor(components, platform);						
						//returnCluster = (vafTumor >= 0.5) ? ClusterType.LOHref : ClusterType.LOH; //LOH
						eventAtSite = ClusterType.LOH;
					}

				} else {
					eventAtSite = ClusterType.Null;   //outside NAF frame (<=> 'other')
				}
			}
			
			// Now finally assign the cluster.
			CompareUtils.ensureNotNull(eventAtSite, "ERROR: Event type for point cannot be unassigned!");
			returnEvents.set(row, eventAtSite);
		}

		return returnEvents;
	}

	// ========================================================================	
	private static void getValidPointsForClustering(ArrayList<ClusteringInputOneSite> rows, // input
													ArrayList<Floint> pointsUpperPlane,    // output
													ArrayList<Floint> pointsLowerPlane,    // output
													int[] indexMapToUpperLowerPlane,       // output - if positive, to upper plane, if negative, to lower plane, if does not map, contains Integer.MAX_VALUE													
													IntArrayList mapToRowsFromUpper, 
													IntArrayList mapToRowsFromLower,
													ClusteringInputOneSampleMetaData metaData,
													double scalingFactor,
													SeqPlatform platform, 
													RangeDouble vafNormalRange) {
		
		pointsUpperPlane.clear();
		pointsLowerPlane.clear();
//		double alphaAdjusted = ParamsBool.IgnoreMultipleTesting.getValue() ? 
//				 PValueBinDistAlpha_UpperPlaneThresh :
//				(PValueBinDistAlpha_UpperPlaneThresh * 0.008 /* Light FDR */);
				
		for (int row = 0; row < rows.size(); row++) {
			ClusteringInputOneSite oneSiteInfo = rows.get(row);			
			
			float vafNormal = oneSiteInfo.calcVAFNormal();
			if (vafNormalRange.inRangeLowerExclusive(vafNormal)) {				
				Floint thePoint = new Floint(metaData.mAdjustedVAFTumor[row], metaData.mAdjustedVAFNormal[row], (float) (metaData.mImbalancePValuesTumor[row] * scalingFactor));
				
				// Assume p-value as vertical row factor
				boolean tumorSigImbalanced = (metaData.mImbalancePValuesTumor[row] <= metaData.mFDRTumor);
				Chrom chrom = oneSiteInfo.getChrom();
				boolean normalSigImbalanced = /* (imbalancePValuesNormal[row] <= fdrNormal) && */ (metaData.chromHasGermlineGain(chrom));
				
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
		case Illumina: return  Float.parseFloat(StringUtils.extractNthColumnValue(line, Script.Col_NAFTAFInput_VariantRatioNormal, StringUtils.TabStr));
		case SOLiD:    return (Float.parseFloat(StringUtils.extractNthColumnValue(line, 6, StringUtils.TabStr)) / 
				               Float.parseFloat(StringUtils.extractNthColumnValue(line, 5, StringUtils.TabStr)));
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
		case Illumina: return  Float.parseFloat(StringUtils.extractNthColumnValue(line, Script.Col_NAFTAFInput_VariantRatioTumor, StringUtils.TabStr));
		case SOLiD:    return (Float.parseFloat(StringUtils.extractNthColumnValue(line, 4, StringUtils.TabStr)) / 
				               Float.parseFloat(StringUtils.extractNthColumnValue(line, 3, StringUtils.TabStr)));
		}
		return Float.NaN;		
	}
	
	// ========================================================================
	/** Returns whether the variant-allele fraction in the normal indicates a germline or somatic allele. */
	public static boolean isVariantInGermline(double vafNormal) {
		return vafNormal >= 0.05;
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
