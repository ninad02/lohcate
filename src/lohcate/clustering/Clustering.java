package lohcate.clustering;
import genomeEnums.Chrom;
import genomeEnums.Nuc;
import genomeEnums.TissueType;
import genomeEnums.VariantLocation;
import genomeUtils.ChromPositionTracker;
import genomeUtils.ElementPlaneSplit;
import genomeUtils.GenomeConstants;
import genomeUtils.GenotypeUtils;
import genomeUtils.MaxPositionsOnChromosome;
import genomeUtils.ObjectWalkerTracker;
import genomeUtils.RegionRange;
import genomeUtils.SNVMap;
import genomeUtils.SiteList;

import java.io.BufferedWriter;
import java.io.File;
import java.io.PrintStream;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.Collections;
import java.util.Date;
import java.util.ListIterator;

import kMeans.Cluster;
import kMeans.DataPoint;
import kMeans.JCA;

import lohcate.AllelicBiasTable;
import lohcate.CopyNumberRegionRange;
import lohcate.CopyNumberRegionsByChromosome;
import lohcate.LOHcate;
import lohcate.LOHcate.SubdirsDefault;
import lohcate.LOHcateSimulator;
import lohcate.Regions;
import lohcate.LOHcateSimulator.LOHcateSimulatorParams;
import lohcateEnums.EventType;
import lohcateEnums.SeqPlatform;

import nutils.ArgumentParserUtils.InputParameterInteger;
import nutils.ArrayUtils;
import nutils.ArrayUtils.ParallelArrayDouble;
import nutils.ArrayUtils.ParallelArrayDoubleDynamic;
import nutils.Cast;
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
import nutils.counter.BucketCounterEnumMatrix;
import nutils.counter.DynamicBucketCounter;
import nutils.counter.DynamicRoundedDoubleCounter;
import nutils.math.BinomialTestPermutationValues;
import nutils.math.PoissonDistributionList;

import org.apache.commons.math3.distribution.PoissonDistribution;
import org.apache.commons.math3.genetics.PermutationChromosome;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.stat.Frequency;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.apache.commons.math3.stat.regression.SimpleRegression;
import org.jfree.data.xy.DefaultXYDataset;
import org.jfree.data.xy.XYSeries;

import JISTIC.convertSEG;

import com.carrotsearch.hppc.DoubleArrayList;
import com.carrotsearch.hppc.IntArrayList;
import com.carrotsearch.hppc.cursors.DoubleCursor;
import com.martiansoftware.jsap.JSAP;
import com.sun.org.apache.xerces.internal.dom.ElementDefinitionImpl;

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
	
	public static final String GISTIC_Input_HeaderString = "SampleID\tchrom\tposStart\tposEnd\tnumMarkers\tSegmentMean";			
	
	public static boolean correctAllelicBias() { return !ClusteringParams.GlobalClusteringParams.mIgnoreAllelicBias.getValue(); }	
	public static double[] sigPValues = BinomialTestPermutationValues.initialize(new double[65000], -1);
	

	// ========================================================================
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		ClusteringInputOneSite.TestClusteringInputOneSite_Robust();
	}
	
	// ========================================================================
	private static void classifySitesHelper_MakeSubDirs(LOHcate.Subdirs lohcateDirs) {
		lohcateDirs.createDirectoryPathOnDisk(SubdirsDefault.SitesClassified);
		lohcateDirs.createDirectoryPathOnDisk(SubdirsDefault.Plots_VAF_2D);
		lohcateDirs.createDirectoryPathOnDisk(SubdirsDefault.Plots_VAF_GenomeWide);
		lohcateDirs.createDirectoryPathOnDisk(SubdirsDefault.Plots_CopyNumber);
		lohcateDirs.createDirectoryPathOnDisk(SubdirsDefault.Plots_Recurrence);
		lohcateDirs.createDirectoryPathOnDisk(SubdirsDefault.Simulation);
		
		lohcateDirs.createDirectoryPathOnDisk(SubdirsDefault.Regions_GISTIC);
		lohcateDirs.createDirectoryPathOnDisk(SubdirsDefault.Regions_JISTIC);
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
			
			final Chrom chrom  = Chrom.getChrom(   StringUtils.extractNthColumnValue(line, Regions.Col_NAFTAFInput_Chrom,     StringUtils.FileExtensionTSV.mDelimiter) );
			final int position = Integer.parseInt( StringUtils.extractNthColumnValue(line, Regions.Col_NAFTAFInput_Position,  StringUtils.FileExtensionTSV.mDelimiter) );

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
		ArrayList<String> allVariantRowsStr = readLinesFromFilesAsStringList(file, true, Regions.EliminateHighDensitySNVs, Regions.EliminateExtremeGCSites);
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
		System.out.printf("\tNum Lines in File: %d\n", germlineSpecificVariantRows.size() + somaticSpecificVariantRows.size());				
		
		//String headerStringGermline = germlineSpecificVariantRows.get(0);
		//String headerStringSomatic  = somaticSpecificVariantRows.get(0);

		//upstream pipelines will randomly spit header lines into the 
		// middle of a naf-taf-input file. we're just avoiding those
		germlineSpecificVariantRows = Regions.curateSNPCalls_removeHeaderLinesFromRows(germlineSpecificVariantRows, removeHighDensitySNVSites, removeExtremeGCSites);
		somaticSpecificVariantRows  = Regions.curateSNPCalls_removeHeaderLinesFromRows(somaticSpecificVariantRows,  removeHighDensitySNVSites, removeExtremeGCSites);
		System.out.printf("\tNum Sites Retained after Header and GC Removal: %d\n", germlineSpecificVariantRows.size() + somaticSpecificVariantRows.size());
		
		// Create a combined set of rows
		ArrayList<String> allVariantRowsStr = new ArrayList<String>(germlineSpecificVariantRows.size() + somaticSpecificVariantRows.size());
		allVariantRowsStr.addAll(germlineSpecificVariantRows);
		//int indexFirstSomaticRowInAllVariants = allVariantRowsStr.size();  // save the size with just the germline variants added
		allVariantRowsStr.addAll(somaticSpecificVariantRows);

		//Collections.sort(germlineSpecificVariantRows, LineComparatorTab);
		//Collections.sort(somaticSpecificVariantRows,  LineComparatorTab);
		Collections.sort(allVariantRowsStr,              Regions.LineComparatorTab);
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
			afStatsSample.tabulateAndPerformStatistics(oneSampleData, platform);
			metaData.mVAFNormalHetRange.setBoundLower( afStatsSample.getValueNStandardDeviationsAway(-NAF_STRIP_EXPANDER) );
			metaData.mVAFNormalHetRange.setBoundUpper( afStatsSample.getValueNStandardDeviationsAway( NAF_STRIP_EXPANDER) );
		}		
		
		
/*		for (Chrom chrom : Chrom.values()) {
			if (!chrom.isInvalid() && chrom.isAutosomal()) {
				System.out.println("Chrom Index: " + chrom + "\t" + oneSampleData.getIndexChromStart(chrom) + "\t" + oneSampleData.getIndexChromEnd(chrom));
				System.out.println(oneSampleData.getSiteAtIndex(oneSampleData.getIndexChromStart(chrom)).toString());
				System.out.println(oneSampleData.getSiteAtIndex(oneSampleData.getIndexChromEnd(chrom)).toString());
			}
		}*/
		
		// First, get the copy number ratios at a per-chromosome level				
		float[] tumorNormalCopyNumRatiosPerChrom = calcAvgCoverageRatioPerChrom(oneSampleData, metaData);
						
		// Now get the copy number ratios at a segmented sub-chromosomal level (by gene)
		calcRoughCopyNumberRatioPerSite(oneSampleData, metaData);	
		
		// Now smooth the copy numbers
		smoothCopyNumbers(oneSampleData, metaData);
		
		// Now adjust the VAF values basead on biases
		calculateAdjustedVAFs(oneSampleData, metaData, allelicBiasTable, platform);

		// Calculate the imbalance p-values
		getPValuesImbalance(oneSampleData, metaData);
		
		// Calculate the FDR for tumor and normal
		boolean pointsIndependent = false;
		double[] tumorPValWithoutGermline = getSitesPvaluesAfterEliminatingGermlineGain(oneSampleData, metaData, metaData.mVAFNormalHetRange);
		//metaData.mFDRTumor  = NumberUtils.getFDR_BenjaminiHochberg(metaData.mImbalancePValuesTumor,  ClusteringParams.GlobalClusteringParams.mFDRAlpha.getValue(), pointsIndependent);
		metaData.mFDRTumor  = NumberUtils.getFDR_BenjaminiHochberg(tumorPValWithoutGermline,         ClusteringParams.GlobalClusteringParams.mFDRAlpha.getValue(), pointsIndependent);
		metaData.mFDRNormal = NumberUtils.getFDR_BenjaminiHochberg(metaData.mImbalancePValuesNormal, ClusteringParams.GlobalClusteringParams.mFDRAlpha.getValue(), pointsIndependent);
		
		// Calculate permutation-based significance thresholds
		calculatePermutationThresholds(oneSampleData, metaData);
	}

	// ========================================================================
	private static void calculatePermutationThresholds(ClusteringInputOneSample oneSampleData, ClusteringInputOneSampleMetaData metaData) {
		int numIterations = 10000;
		double[] pValueBuffer = new double[numIterations];
		
		System.out.println("Permuting for Significant P-Value Calculation...");
		int numSites = oneSampleData.getNumSites();
		for (int row = 0; row < numSites; row++) {
			ClusteringInputOneSite oneSiteInfo = oneSampleData.getSiteAtIndex(row); 
						
			if (!metaData.chromHasGermlineGain(oneSiteInfo.getChrom()) && metaData.mVAFNormalHetRange.inRangeLowerExclusive(oneSiteInfo.calcVAFNormal())) {
				if (sigPValues[oneSiteInfo.mCovgTotalTumor] < 0) {
					sigPValues[oneSiteInfo.mCovgTotalTumor] = GenotypeUtils.calcSigPValueByPermutation(oneSiteInfo.mCovgTotalTumor, 0.05, pValueBuffer); 
				}				
				metaData.mSigPValuesPerSite[row] = sigPValues[oneSiteInfo.mCovgTotalTumor];
			} else {
				metaData.mSigPValuesPerSite[row] = -1;
			}
			//System.out.println(row + "\t" + oneSiteInfo.mCovgTotalTumor + "\t" + metaData.mSigPValuesPerSite[row]);
		}
		System.out.println("Finished Permuting for Significant P-Value Calculation...");
	}
	
	// ========================================================================
	// Conve
	private static double[] getSitesPvaluesAfterEliminatingGermlineGain(ClusteringInputOneSample oneSampleData, ClusteringInputOneSampleMetaData metaData, RangeDouble vafNormalRange) {		
		int numSites = oneSampleData.getNumSites();
		DoubleArrayList pValues = new DoubleArrayList(numSites);
		
		for (int row = 0; row < numSites; row++) {
			ClusteringInputOneSite oneSiteInfo = oneSampleData.getSiteAtIndex(row); 
			
			if (!metaData.chromHasGermlineGain(oneSiteInfo.getChrom()) /* && vafNormalRange.inRangeLowerExclusive(oneSiteInfo.calcVAFNormal()) */) {
				pValues.add(metaData.mImbalancePValuesTumor[row]);
			}
		}
		System.out.printf("Shortened Pvalue List:\t%d\t%d\n", numSites, pValues.size());
		return pValues.toArray();
	}
	
	// ========================================================================
	private static Boolean classifySites_addToContingencyTables(EventType eventTest, EventType eventTruth, EnumMapSafe<EventType, ContingencyTable> eventTablesMap) {		
		
		if (CompareUtils.isNull(eventTest) || CompareUtils.isNull(eventTruth)) {
			return null;
		}
		
		if (eventTruth == EventType.GainGermline && eventTest == EventType.Ignored) {
			return Boolean.FALSE;
		}

		for (EventType ct : EventType.values()) {
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
	private static Boolean classifySites_addToEventCountTables(EventType eventTruth, EventType eventTest, BucketCounterEnumMatrix<EventType> counter) {
		if (CompareUtils.isNull(eventTest))  return null;
		if (CompareUtils.isNull(eventTruth)) return null;
		
		if (eventTruth == EventType.GainGermline && eventTest == EventType.Ignored) return Boolean.FALSE;
				
		counter.increment(eventTruth, eventTest);
		return Boolean.TRUE;
	}

	// ========================================================================
	public static AllelicBiasTable constructAllelicBiasTable(ArrayList<File> files, AllelicBiasTable biasTable, SNVMap snvMap, FileExtensionAndDelimiter fileExtAndDelim, SeqPlatform platform) {
		biasTable = CompareUtils.isNull(biasTable) ? new AllelicBiasTable() : biasTable;
		for (File inFile : files) {
			ArrayList<String> allLines = readLinesFromFilesAsStringList(inFile, true, true, true);
			for (String line : allLines) { 
				final Chrom chrom  = Chrom.getChrom(   StringUtils.extractNthColumnValue(line, Regions.Col_NAFTAFInput_Chrom,    fileExtAndDelim.mDelimiter) );
				final int position = Integer.parseInt( StringUtils.extractNthColumnValue(line, Regions.Col_NAFTAFInput_Position, fileExtAndDelim.mDelimiter) );
				final float vafNormal = extractVAFNormal(line, platform);
				
				if (AlleleFractionStatsForSample.VAFNormalRange.inRangeBothExclusive(vafNormal)) {
					biasTable.registerSite(chrom, position, vafNormal);
				}
				
				if (snvMap != null) {
					snvMap.registerSNV(chrom, position, 0, Nuc.N, Nuc.N, true, true);
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
	public static void classifySites(LOHcate.Subdirs lohcateDirs, String allelicBiasInFile, SeqPlatform platform, LOHcateSimulator.LOHcateSimulatorParams simulatorParams, String simOutRootFilename) {
		File[] files = (new File(lohcateDirs.getSubDirPath(SubdirsDefault.VAFInputsNormalTumor))).listFiles();

		
		// Get the list of valid files
		ArrayList<File> validFiles = new ArrayList<File>(files.length);		
		for (File file : files) {
			int indexOfSubstring = file.getName().indexOf(Regions.GermlineSuffix);
			if (indexOfSubstring >= 0) {
				validFiles.add(file);
				System.out.println(file.getName().substring(0, indexOfSubstring));
			}
		}
		
		//Arrays.fill(sigPValues, -1);
		
		// We now handle the construction or reading of the allelic bias table.  If we are to correct
		// allelic biases, we either construct it from the input files themselves, or we read pre-calculated
		// allelic biases from an independent file.
		
		// We now calculate allelic biases and read in the genomic map at the same time.  If the user has 
		// selected to provide a pre-defined allelic bias map, we overwrite our map.  However, we still 
		// calculate biases since bias calculation and creation of the genomic map are combined in the
		// same function.
		System.out.println((new Date()).toString() + " Reading sites...");
		SNVMap allSitesMap = new SNVMap();		
		AllelicBiasTable allelicBiasTable = constructAllelicBiasTable(validFiles, null, allSitesMap, StringUtils.FileExtensionTSV, platform);
		System.out.println((new Date()).toString() + " Finished Reading sites...");
				
		if (correctAllelicBias()) {
			if (allelicBiasInFile != null) {
				System.out.println((new Date()).toString() + " Reading Allelic Bias file...");
				allelicBiasTable = AllelicBiasTable.readFileAndConstructTable(allelicBiasInFile, AllelicBiasTable_Col_NumSamples, AllelicBiasTable_Col_AvgVAFNormal);
				System.out.println((new Date()).toString() + " Finished Reading Allelic Bias file...");
			}
		} else {			
			allelicBiasTable = null;   // Reset to since we don't want to handle biases
		}
		
		// Create output directory
		classifySitesHelper_MakeSubDirs(lohcateDirs);
					
		PrintStream simOutputStream = CompareUtils.isNull(simulatorParams) ? System.out : IOUtils.getPrintStream(lohcateDirs.getSubDirPath(SubdirsDefault.Simulation) + File.separator + simOutRootFilename + ".sim.txt");	
		System.out.println(simOutRootFilename);
		
		// Create counters for all events
		EnumMapSafe<EventType, EnumMapSafe<Chrom, DynamicBucketCounter>> eventCounts = new EnumMapSafe<>(EventType.class);
		EnumMapSafe<EventType, ParallelArrayDoubleDynamic> eventCoordinates = new EnumMapSafe<>(EventType.class);
		
		for (EventType event : EventType.values()) {			
			eventCounts.put(event, DynamicBucketCounter.ClassFactory.newEnumMap(Chrom.class));
			eventCoordinates.put(event, new ParallelArrayDoubleDynamic(10000));
		}		
		
		EnumMapSafe<EventType, ContingencyTable> eventTablesTotal = ContingencyTable.ClassFactory.newEnumMap(EventType.class);
		BucketCounterEnumMatrix<EventType> countBreakdownPerEventTotal = new BucketCounterEnumMatrix<>(EventType.class);
		MaxPositionsOnChromosome maxPosOnChrom = new MaxPositionsOnChromosome();				
		
		// Prepare the GISTIC output
		String outFilenameGISTIC = lohcateDirs.getSubDirPath(SubdirsDefault.Regions_GISTIC) + File.separator + "Regions.GISTIC.txt";
		BufferedWriter outGISTIC = IOUtils.getBufferedWriter(outFilenameGISTIC);
		IOUtils.writeToBufferedWriter(outGISTIC, GISTIC_Input_HeaderString, true);
		String outFilenameProbeFile = lohcateDirs.getSubDirPath(SubdirsDefault.Regions_GISTIC) + File.separator + "AllProbes.GISTIC.txt";
		allSitesMap.printMe(true, outFilenameProbeFile, false);		
		
		int fileIndex = 0;		
		for (File file : validFiles) {			
			String filename = file.getName();
			String sampleNameRoot = filename.substring(0, file.getName().indexOf(Regions.GermlineSuffix));  	
			String extension = filename.substring(filename.lastIndexOf(StringUtils.DotStr), filename.length());
			System.out.println("Processing (" + ++fileIndex + " / " + validFiles.size() + "): " + file.getName());				
			classifySitesOneSample(file, sampleNameRoot, extension, fileIndex, allelicBiasTable, allSitesMap, eventCounts, eventCoordinates, eventTablesTotal, countBreakdownPerEventTotal, simulatorParams, simOutputStream, maxPosOnChrom, lohcateDirs, platform, outGISTIC);
		}
				
		// Write out the stats
		printContingencyTableInfo("SIMALL", simOutputStream, 0, eventTablesTotal, countBreakdownPerEventTotal);

		// Close the print stream
		if (CompareUtils.isNotNull(simulatorParams)) {
			IOUtils.closePrintStream(simOutputStream);
		}
		
		IOUtils.closeBufferedWriter(outGISTIC);
		
		// Now plot the LOH recurrence across samples
		ClusteringPlotting.plotRecurrenceGenomeWide(eventCounts, lohcateDirs.getSubDirPath(SubdirsDefault.Plots_Recurrence));
		ClusteringPlotting.plotEventsByCoordinateAcrossSamples(eventCoordinates, maxPosOnChrom, validFiles.size(), lohcateDirs.getSubDirPath(SubdirsDefault.Plots_Recurrence));
		
	
		// Now do the GISTIC processing
		String jisticCovertSEGOutputFilename = lohcateDirs.getSubDirPath(SubdirsDefault.Regions_GISTIC) + File.separator + "Regions.Matrix.JISTIC.txt";
		convertSEG.main(new String[] {outFilenameGISTIC, outFilenameProbeFile, "outputfile=" + jisticCovertSEGOutputFilename});
		String jisticInputPathPrefix = IOUtils.pathConcat(new String[] { "E:", "Research", "Code", "Libraries", "Java", "lib", "JISTIC" }) + File.separator; 
		String jisticlocationsString = "locations=" + jisticInputPathPrefix + "hg18_Gene_Info.txt";
		String jisticCytobandString  = "bands=" + jisticInputPathPrefix + "Human_cytoBand.txt";
		String jisticCopyNumberString = "copynumber=" + jisticCovertSEGOutputFilename;
		
		String jisticInputSpecFile = IOUtils.pathConcat(new String[] { "spec=" + jisticInputPathPrefix + "glioexample", "focal", "GISTICFocal.spec" });		
		//String jisticInputSpecFile = IOUtils.pathConcat(new String[] { "spec=" + jisticInputPathPrefix + "glioexample", "limited", "GISTICLimited.spec" });
		JISTIC.Distribution.main(new String[] { jisticInputSpecFile, jisticCopyNumberString, jisticCytobandString, jisticlocationsString });
	}
	
	// ========================================================================
	public static boolean fillRegionBasedOnVAFMaxLikelihood(RegionRange region, 
														  ClusteringInputOneSample oneSampleData, 
														  ClusteringInputOneSampleMetaData metaData,
														  ArrayList<EventType> events,
														  EventType targetEventType, 
														  PrimitiveWrapper.WDouble probFromMaxLikelihood,
														  boolean fillRegion) {
		
		int indexStart = oneSampleData.getIndex(region.getChromosome(), region.getRangeStart());
		CompareUtils.ensureTrue(indexStart >= 0, "ERROR: Starting index must be > 0");
		int indexEnd   = oneSampleData.getIndex(region.getChromosome(), region.getRangeEnd());
		CompareUtils.ensureTrue(indexEnd >= indexStart, "ERROR: Ending index must be >= starting index!");
		//System.out.println(eventType + "\t" + region);
		
		// We do two iterations.  The first iteration will be based on 							
		int mostLikelyList = 0;		
		boolean isLinearFromRegression = false;
		boolean isLinearFromKMeans = false;
		boolean noDenserPeaks = false;
		
		if (!region.spansOneSite()) {
			String listOfListOfMeansStr = constrctListofListOfMeans();
			ArrayList<DoubleArrayList> listOfListOfMeans = ArrayUtils.getListOfDoubleListsFromStringForm(listOfListOfMeansStr, true);							
			mostLikelyList = determineCopyGainVAFMaxLikelihood(oneSampleData, metaData, region, listOfListOfMeans, TissueType.Tumor, probFromMaxLikelihood);
			System.out.printf("ML:\t%12s\t%d\t%d\t%d\t%d\n", targetEventType, region.getChromosome().ordinal(), indexStart, indexEnd, mostLikelyList);
			//System.out.println("\t" + mostLikelyList);
			
			// Check if still considered an event
			if (mostLikelyList > 0) {									
				DoubleArrayList vafTumorsSorted = new DoubleArrayList();
				for (int row = indexStart; row <= indexEnd; row++) {
					if (metaData.mVAFNormalHetRange.inRangeLowerExclusive(oneSampleData.getSiteAtIndex(row).calcVAFNormal())) {
						vafTumorsSorted.add(metaData.mAdjustedVAFTumor[row]);											
					}											
				}				
				ArrayUtils.sort(vafTumorsSorted);				
				SimpleRegression simpReg = new SimpleRegression(true);
				for (int i = 0; i < vafTumorsSorted.size(); i++) {					
					simpReg.addData(i + 1, vafTumorsSorted.get(i));
				}
				System.out.printf("R:\t%12s\t%d\t%d\t%d\t%g\t%g\n", targetEventType, region.getChromosome().ordinal(), indexStart, indexEnd, simpReg.getRSquare(), simpReg.getSumSquaredErrors()); 		
				if (targetEventType == EventType.cnLOH) {
					//TODO isLinearFromRegression = (simpReg.getSumSquaredErrors() < 0.01);
					//isLinearFromRegression = (simpReg.getSumSquaredErrors() < 0.01);
				}
				
				// Now test via density
				Frequency freqDist = new Frequency();
				double multiplierFreq = 100.0;
				for (int i = 0; i < vafTumorsSorted.size(); i++) {
					int vafTumorMultiplied = Cast.toInt( Math.round(vafTumorsSorted.get(i) * multiplierFreq) );
					freqDist.addValue(vafTumorMultiplied);
				}
				DoubleArrayList mostLikelyMeans = listOfListOfMeans.get(mostLikelyList);
				int mean0 = Cast.toInt(Math.round(mostLikelyMeans.get(0) * multiplierFreq));
				int mean1 = Cast.toInt(Math.round(mostLikelyMeans.get(1) * multiplierFreq));
				int leeway = 4;
				
				double probMean0  = NumberUtils.getCumProbInterval(freqDist, mean0, leeway); 						
				double probMean1  = NumberUtils.getCumProbInterval(freqDist, mean1, leeway);
				double probCenter = NumberUtils.getCumProbInterval(freqDist, 50, leeway);
				boolean denserPeaks = (probMean0 > probCenter) && (probMean1 > probCenter);
				System.out.printf("\tFREQ: %d\t%d\t%g\t%g\t%g\n", mean0, mean1, probMean0, probMean1, probCenter);
				noDenserPeaks = !denserPeaks;
				
				
				// Now do k-means clustering
				if (!isLinearFromRegression) {
					double multiplier = 100.0;
					ArrayList<DataPoint> kMeansPoints = new ArrayList<DataPoint>(vafTumorsSorted.size());
					for (int i = 0; i < vafTumorsSorted.size(); i++) {
						double vafTumorAtSite = vafTumorsSorted.get(i);					
						vafTumorAtSite = NumberUtils.roundToNearest5(Cast.toInt(Math.round(vafTumorAtSite * multiplier))) / multiplier;
						DataPoint dp = new DataPoint(vafTumorAtSite, 0, "");
						kMeansPoints.add(dp);
					}
					JCA kMeansClusterer = new JCA(2, 100, kMeansPoints);
					kMeansClusterer.startAnalysis();
					Cluster c1 = kMeansClusterer.getCluster(0);
					Cluster c2 = kMeansClusterer.getCluster(1);
					double normSumSqC1 = c1.getSumSqr() / c1.getNumDataPoints();
					double normSumSqC2 = c2.getSumSqr() / c2.getNumDataPoints();
					System.out.printf("\tSS:\t%g\t%g\n", normSumSqC1, normSumSqC2);
					double thresh = 0.042;
					//TODO isLinearFromKMeans = (normSumSqC1 >= thresh) && (normSumSqC2 >= thresh);
				}
			}
						
		} 
		
		// Go through VAFs
		boolean resetToHetGermline = (mostLikelyList == 0) || isLinearFromRegression || isLinearFromKMeans || noDenserPeaks;

		if (fillRegion && (events != null)) {
			fillRegion(events, indexStart, indexEnd, resetToHetGermline, targetEventType, oneSampleData, metaData);
		}
		
		return !resetToHetGermline;
	}
	
	// ========================================================================
	private static void fillRegion( ArrayList<EventType> events, 
			int indexStartIncl, int indexEndIncl, boolean resetToHetGermline, 
			EventType targetEventType,  
			ClusteringInputOneSample oneSampleData, 
			ClusteringInputOneSampleMetaData metaData) {
		
		for (int row = indexStartIncl; row <= indexEndIncl; row++) {									
			final EventType event = events.get(row);
			if (resetToHetGermline) {
				if (event == targetEventType) {
					events.set(row, EventType.HETGermline);
				}
			} else {
				if (event == EventType.HETGermline || event == EventType.Noise || event == EventType.Ignored) {
					if (metaData.mVAFNormalHetRange.inRangeLowerExclusive(oneSampleData.getSiteAtIndex(row).calcVAFNormal())) {
						events.set(row, targetEventType);
					}
				}
			}
		}
	}
	
	// ========================================================================
	public static boolean combineTwoRegions(CopyNumberRegionRange regionPrev, CopyNumberRegionRange regionCurr, ClusteringInputOneSample oneSampleData, ClusteringInputOneSampleMetaData metaData) {
		
		// Check that we have the same chromosome
		if (regionPrev.getChromosome() != regionCurr.getChromosome()) {
			return false;
		}
		
		// Test that they have teh same event type
		if (regionPrev.mCopyNumberEventType != regionCurr.mCopyNumberEventType) {
			return false;
		}
		
		int indexRegionPrevStart = oneSampleData.getIndex(regionPrev.getChromosome(), regionPrev.getRangeStart()); 
		CompareUtils.ensureTrue(indexRegionPrevStart >= 0, "ERROR: Starting index must be > 0");
		int indexRegionPrevEnd   = oneSampleData.getIndex(regionPrev.getChromosome(), regionPrev.getRangeEnd());
		CompareUtils.ensureTrue(indexRegionPrevStart >= indexRegionPrevStart, "ERROR: Ending index must be >= starting index!");
		

		int indexRegionCurrStart = oneSampleData.getIndex(regionCurr.getChromosome(), regionCurr.getRangeStart()); 
		CompareUtils.ensureTrue(indexRegionCurrStart >= 0, "ERROR: Starting index must be > 0");
		int indexRegionCurrEnd   = oneSampleData.getIndex(regionCurr.getChromosome(), regionCurr.getRangeEnd());
		CompareUtils.ensureTrue(indexRegionCurrStart >= indexRegionCurrStart, "ERROR: Ending index must be >= starting index!");

		int numHetSitesMid = oneSampleData.getNumHetSitesInRangeByIndex(indexRegionPrevEnd + 1, indexRegionCurrStart - 1, metaData.mVAFNormalHetRange);
		PrimitiveWrapper.WDouble probFromMaxLikelihood = new PrimitiveWrapper.WDouble(0);
		CopyNumberRegionRange midRange = new CopyNumberRegionRange(regionPrev);
		
		// Check how many het sites there are between regions.  If there are none, we can
		// treat as if the regions are adjacent.		
		if (numHetSitesMid == 0) {
			return true;
			
		} else {
			int midRegionStartPos = oneSampleData.getSiteAtIndex(regionPrev.getChromosome(), indexRegionPrevEnd   + 1).getPosition();
			int midRegionEndPos   = oneSampleData.getSiteAtIndex(regionPrev.getChromosome(), indexRegionCurrStart - 1).getPosition();
			InputParameterInteger numHetSitesThreshold = new InputParameterInteger(6, "", JSAP.NO_SHORTFLAG, "", "");
			
			if (numHetSitesMid < numHetSitesThreshold.getValue()) {
				midRegionStartPos = regionPrev.getRangeStart();
				midRegionEndPos   = regionCurr.getRangeEnd();
			}
							
			midRange.set(regionPrev.getChromosome(), midRegionStartPos, midRegionEndPos, false, indexRegionCurrStart - indexRegionPrevEnd - 1);
			CompareUtils.ensureTrue(!midRange.spansOneSite(), "Must always span more than one site!");
			boolean isEventRegion = fillRegionBasedOnVAFMaxLikelihood(midRange, oneSampleData, metaData, null, regionPrev.mCopyNumberEventType, probFromMaxLikelihood, false);
			return isEventRegion;
		}
	}

	
	// ========================================================================
	private static String constrctListofListOfMeans() {
		StringBuilder sb = new StringBuilder(1024);
		sb.append("({0.5}");
		//sb.append(";{0.45,0.55}");
		sb.append(";{0.4,0.6}");
		sb.append(";{0.3333,0.6667}");		
		sb.append(";{0.3,0.7}");
		sb.append(";{0.25,0.75}");
		sb.append(";{0.2,0.8}");
		sb.append(";{0.15,0.85}");		
		sb.append(";{0.1,0.9}");
		sb.append(";{0.05,0.95}");
		sb.append(";{0.005,0.995}");
		sb.append(";{0.0001,0.9999}");
		sb.append(")");
		return sb.toString();
	}
	
	// ========================================================================
	public static void classifySitesOneSample(File file, String sampleNameRoot, String sampleFilenameExtension, 
			int sampleIndex, 
			AllelicBiasTable allelicBiasTable,
			SNVMap allSitesMap,
			EnumMapSafe<EventType, EnumMapSafe<Chrom, DynamicBucketCounter>> eventCountsByCoordinate, 
			EnumMapSafe<EventType, ParallelArrayDoubleDynamic> eventsByCoordinate,
			EnumMapSafe<EventType, ContingencyTable> eventTablesTotal,
			BucketCounterEnumMatrix<EventType> countBreakdownPerEventTotal,
			LOHcateSimulatorParams simulatorParams, PrintStream simOutputStream,
			MaxPositionsOnChromosome maxPosOnChrom,
			LOHcate.Subdirs lohcateDirs,
			//String sitesClassifiedDir, String vafComparisonPlotDir, String vafWaterfallPlotDir, String copyNumberPlotDir, 
			SeqPlatform platform,
			BufferedWriter outGISTIC) {
		
		StringUtils.FileExtensionAndDelimiter fileExtDelim = StringUtils.FileExtensionTSV;	
		StringBuilder sb = new StringBuilder(8192);
		
		// Get the file and data
		File inFile = file;
		if (!inFile.exists()) {
			CompareUtils.ensureTrue(false, "File does not exist: " + inFile.getAbsolutePath());
		}		
		ClusteringInputOneSample oneSampleData = readLinesFromFiles(inFile);
		if (oneSampleData.getNumSites() == 0) return;
		oneSampleData.mSampleNameRoot = sampleNameRoot;

		String[] columnHeaders = new String[] {
			"chr", "coordinate", "refBase", "varBase", "variantBase-N", "variantBase-T", "refEnv-N", "refEnv-T", "Q20_TotCov_N", "Q20_TotCov_T", "Q20_VarCov_N", "Q20_VarCov_T", "Q20_VariantRatio_N", "Q20_VariantRatio_T", "dbsnp", "MutationType", "Hugo_Symbol", 
			"Event", "TumorNormalCopyNumRatio", "TumorCopyNum", "Adjusted_Q20_VariantRatio_T", "Adjusted_Q20_VariantRatio_N", "pVal_Imbalance_Tumor", "pVal_Imbalance_Normal" 
		};
		String headerStr = StringUtils.constructColumnDelimitedString(columnHeaders, fileExtDelim.mDelimiter, sb, true).toString();
		
		// Initialize output file
		String outFilename = sampleNameRoot + ".withCopyNum" + sampleFilenameExtension;
		String outFilenameFullPath = lohcateDirs.getSubDirPath(SubdirsDefault.SitesClassified) + File.separator + outFilename;
		PrintWriter out = new PrintWriter(IOUtils.getBufferedWriter(outFilenameFullPath));
		out.println(headerStr);
		
		// Initialize output file for simulation
		boolean isSimulation = simulatorParams != null;
		String outSimFilename = sampleNameRoot + ".simStats" + sampleFilenameExtension;
		String outSimFilenameFullPath = lohcateDirs.getSubDirPath(SubdirsDefault.SitesClassified) + File.separator + outSimFilename;
		PrintWriter outSim = isSimulation ? new PrintWriter(IOUtils.getBufferedWriter(outSimFilenameFullPath)) : null;
		
		// Initialize tables for simulation-evaluation tallying
		EnumMapSafe<EventType, ContingencyTable> eventTables = ContingencyTable.ClassFactory.newEnumMap(EventType.class);
		BucketCounterEnumMatrix<EventType> testCounts = new BucketCounterEnumMatrix<EventType>(EventType.class); 

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
			ClusteringInputOneSampleMetaData metaData = new ClusteringInputOneSampleMetaData(oneSampleData.getNumSites(), sampleNameRoot);
			classifySitesHelper_preprocessMetaData(oneSampleData, metaData, allelicBiasTable, platform);				

			System.out.printf("FDR Alpha:\t%g\nFDR Tumor p-value: \t%g\nFDR Normal p-value:\t%g\n", ClusteringParams.GlobalClusteringParams.mFDRAlpha.getValue(), metaData.mFDRTumor, metaData.mFDRNormal);
			System.out.println("\tInferred Copy Number and Allelic Imbalance Per Site...");					

			// KEY STEP: Obtain the clusters
			int indexFirstSomaticRowInAllVariants = 0;  // Set as invalid value for now
			ClusteringResults<EventType> clusteringResults = Clustering.getClusters_withPlane(oneSampleData, metaData, metaData.mVAFNormalHetRange, outFilenameFullPath, indexFirstSomaticRowInAllVariants, platform); //get cluster assignments (HET ball, LOH sidelobes, DUP wedge, &c.)
			ArrayList<EventType> events = clusteringResults.mClassificationForPoint;
			System.out.println("\tGot clusters");

			// BEGIN POST-PROCESSING
			boolean doOverride = true;
			if (doOverride) {
				CopyNumberRegionsByChromosome regionsByChrom = Regions.segmentRegionsOneSample(oneSampleData, events, null, null);						
			
				//EnumMapSafe<EventType, CopyNumberRegionsByChromosome> regionsInSamplePerEventType = 
				//		new EnumMapSafe<EventType, CopyNumberRegionsByChromosome>(EventType.class);

				RegionRange midRange = new RegionRange(Chrom.c0, 0);
				PrimitiveWrapper.WDouble probFromMaxLikelihood = new PrimitiveWrapper.WDouble(0);
//				int numMidPossibilities = 2;
//				double[] probMLMidPossibilities = new double[numMidPossibilities];
//				int[] midRegionStart = new int[numMidPossibilities];
//				int[] midRegionEnd   = new int[numMidPossibilities];
				
				for (EventType eventType : EventType.AmpLOHcnLOH) {
					CopyNumberRegionsByChromosome regionsInOneSampleMerged = Regions.mergeRegionsWithConstraints(regionsByChrom, eventType, Regions.REGION_SEGMENTATION_DIST_THRESHOLD, oneSampleData, metaData);
					//regionsInOneSampleMerged.removeSingletonRegions();
					//regionsInSamplePerEventType.put(eventType, regionsInOneSampleMerged);
					//regionsInOneSampleMerged.print(System.out, StringUtils.FileExtensionTSV.mDelimiter);
										
					for (Chrom chrom : Chrom.values()) {
						ArrayList<CopyNumberRegionRange> regionsOnChrom = regionsInOneSampleMerged.getRegions(chrom);
						CopyNumberRegionRange regionPrev = null;
						for (CopyNumberRegionRange region : regionsOnChrom) {							
							CompareUtils.ensureTrue(region.getChromosome() == chrom, "ERROR: Chromosomes should match!");							
							fillRegionBasedOnVAFMaxLikelihood(region, oneSampleData, metaData, events, eventType, probFromMaxLikelihood, true);							
							//System.out.println("CURR\t" + region.toString());
							
							boolean fillMid = false;
							// Fill in the gap between regions if it exists
							if (regionPrev != null && fillMid) {
								int indexEndPrev   = oneSampleData.getIndex(regionPrev.getChromosome(), regionPrev.getRangeEnd());
								int indexStartCurr = oneSampleData.getIndex(    region.getChromosome(),     region.getRangeStart());
								
								// Ensure we aren't comparing index-adjacent regions.  
								if (indexStartCurr > indexEndPrev + 1) {
									
									// Test whether the middle region has enough sites on its own to be checked
									int numHetSitesMid = oneSampleData.getNumHetSitesInRangeByIndex(indexEndPrev + 1, indexStartCurr - 1, metaData.mVAFNormalHetRange);
									int sufficientNumHetSitesMid = 8;
									
									if (numHetSitesMid >= sufficientNumHetSitesMid) {
										int midRegionStartPos = oneSampleData.getSiteAtIndex(chrom, indexEndPrev   + 1).getPosition();
										int midRegionEndPos   = oneSampleData.getSiteAtIndex(chrom, indexStartCurr - 1).getPosition();
										midRange.set(chrom, midRegionStartPos, midRegionEndPos, false, indexStartCurr - indexEndPrev - 1);
										fillRegionBasedOnVAFMaxLikelihood(midRange, oneSampleData, metaData, events, eventType, probFromMaxLikelihood, true);
										
									} else {
										/*
										Arrays.fill(probMLMidPossibilities, 0);
										
										for (int pIndex = 0; pIndex < numMidPossibilities; pIndex++) {
											switch(pIndex) {
											case 0: 
												midRegionStart[pIndex] = regionPrev.getRangeStart();   
												midRegionEnd[pIndex] = region.getRangeEnd();   
												break;
											case 1: 
												midRegionStart[pIndex] = oneSampleData.getSiteAtIndex(chrom, indexEndPrev + 1).getPosition();   
												midRegionEnd[pIndex] = region.getRangeEnd();   
												break;
											default: CompareUtils.ensureTrue(false, "ERROR: Invalid mid-range possibility!"); break;
											}
											
											midRange.set(chrom, midRegionStart[pIndex], midRegionEnd[pIndex], false, indexStartCurr - indexEndPrev - 1);
										}*/
										
										int midRegionStartPos = regionPrev.getRangeStart(); //;								
										int midRegionEndPos   = region.getRangeEnd(); //;
										midRange.set(chrom, midRegionStartPos, midRegionEndPos, false, indexStartCurr - indexEndPrev - 1);
										boolean result = fillRegionBasedOnVAFMaxLikelihood(midRange, oneSampleData, metaData, events, eventType, probFromMaxLikelihood, false);
										if (result) {											
											fillRegionBasedOnVAFMaxLikelihood(midRange, oneSampleData, metaData, events, eventType, probFromMaxLikelihood, true);
										}
											//System.out.println("MID\t" + midRange.toString());										
									}									
								}
							}
							
							
							regionPrev = region;
						}
					}
					
					// Now outout					
					
				}			
				
				// Now, we segment the sample again since now all the sites have been 
				// processed, and any unassigned sites between CNA regions have been filled in
				CopyNumberRegionsByChromosome regionsByChromPostFill = Regions.segmentRegionsOneSample(oneSampleData, events, null, null);
				Regions.ActionerGISTIC actionerGistic = new Regions.ActionerGISTIC(outGISTIC, regionsByChromPostFill, events, oneSampleData, metaData, allSitesMap);				
			}
			
			// Now initialize the data structure needed to plot
			BucketCounterEnum<EventType> eventCounts = ArrayUtils.getEnumTypeCounts(events, EventType.class);			
			EnumMapSafe<EventType, ParallelArrayDouble> clusterCoordinates  = getCoordinateArraysPerClusterType(eventCounts, EventType.class);
			EnumMapSafe<EventType, ParallelArrayDouble> waterfallPlotTumor  = getCoordinateArraysPerClusterType(eventCounts, EventType.class);
			EnumMapSafe<EventType, ParallelArrayDouble> waterfallPlotNormal = getCoordinateArraysPerClusterType(eventCounts, EventType.class);
			EnumMapSafe<EventType, ParallelArrayDouble> copyNumPlot         = getCoordinateArraysPerClusterType(eventCounts, EventType.class);

			// Do the post-processing
			chromPosTracker.clear();
			
			ParallelArrayDoubleDynamic chromBoundaryXY = new ParallelArrayDoubleDynamic();			
			
			int startingRowGermlineOrSomaticOrAll = 0;  // 0 because header line has been stripped away
			int numSites = oneSampleData.getNumSites();
			for (int row = startingRowGermlineOrSomaticOrAll; row < numSites; row++) {
				ClusteringInputOneSite oneSiteInfo = oneSampleData.getSiteAtIndex(row);						
				EventType eventType = events.get(row);

				// Print to output
				if (numLoopIterations == 1) {
					out.print(oneSiteInfo.printToString(sb, true, StringUtils.FileExtensionTSV.mDelimiter).toString());
					out.printf("\t%s\t", eventType.name());
					metaData.printSiteInformation(out, row, true);
					out.flush();
				}
								
				EventType eventTruth = goldStandard.getEvent(row);
				if (metaData.mVAFNormalHetRange.inRangeLowerExclusive(oneSiteInfo.calcVAFNormal())) {
					classifySites_addToContingencyTables(eventType, eventTruth, eventTables);
					classifySites_addToContingencyTables(eventType, eventTruth, eventTablesTotal);
					classifySites_addToEventCountTables(eventTruth, eventType, testCounts);				
					classifySites_addToEventCountTables(eventTruth, eventType, countBreakdownPerEventTotal);
				}

				
				boolean chromCrossed = chromPosTracker.chromCrossedWithCurrentCoordinates(oneSiteInfo.getChrom(), oneSiteInfo.getPosition());						
				for (double d = 0; chromCrossed && (d <= 5.0); d += 0.02) {
					chromBoundaryXY.add(chromPosTracker.getPositionGenomeWide(), d);
				}
				maxPosOnChrom.registerPosition(oneSiteInfo.getChrom(), chromPosTracker.getPositionGenomeWide());
				
				if (eventType.isSomaticCopyNumberEvent() || (eventType == EventType.GainGermline)) {
					eventsByCoordinate.get(eventType).add(chromPosTracker.getPositionGenomeWide(), sampleIndex);
				}

				// Add recurrence count
				eventCountsByCoordinate.get(eventType).get(oneSiteInfo.getChrom()).incrementCount(oneSiteInfo.getPosition());

				clusterCoordinates .get(eventType).add( metaData.mAdjustedVAFTumor[row],         metaData.mAdjustedVAFNormal[row] );
				waterfallPlotTumor .get(eventType).add( chromPosTracker.getPositionGenomeWide(), metaData.mAdjustedVAFTumor[row]);
				waterfallPlotNormal.get(eventType).add( chromPosTracker.getPositionGenomeWide(), metaData.mAdjustedVAFNormal[row]);
				copyNumPlot        .get(eventType).add( chromPosTracker.getPositionGenomeWide(), (metaData.mTumorCopyNumRatiosPerGene[row] * Regions.DefaultDiploidCopyNumber));
			}									

			// Now let's create the datasets needed to
			double[][] boundaryArrays = chromBoundaryXY.toArrays();
			DefaultXYDataset xyDatasetVAFPlot             = classifySitesHelper_createAndFillXYData(clusterCoordinates, null);
			DefaultXYDataset xyDatasetWaterfallPlotTumor  = classifySitesHelper_createAndFillXYData(waterfallPlotTumor,  boundaryArrays);
			DefaultXYDataset xyDatasetWaterfallPlotNormal = classifySitesHelper_createAndFillXYData(waterfallPlotNormal, boundaryArrays);
			DefaultXYDataset xyDatasetCopyNumber          = classifySitesHelper_createAndFillXYData(copyNumPlot, boundaryArrays);

			ClusteringPlotting.plotVAFComparison(xyDatasetVAFPlot,             lohcateDirs.getSubDirPath(SubdirsDefault.Plots_VAF_2D)         + File.separator + sampleNameRoot + ".VAFComparison", sampleNameRoot);
			ClusteringPlotting.plotVAFGenomeWide(xyDatasetWaterfallPlotTumor,  lohcateDirs.getSubDirPath(SubdirsDefault.Plots_VAF_GenomeWide) + File.separator + sampleNameRoot + ".VAF_GenomeWide_Tumor",  sampleNameRoot, true);				
			ClusteringPlotting.plotVAFGenomeWide(xyDatasetWaterfallPlotNormal, lohcateDirs.getSubDirPath(SubdirsDefault.Plots_VAF_GenomeWide) + File.separator + sampleNameRoot + ".VAF_GenomeWide_Normal", sampleNameRoot, false);
			ClusteringPlotting.plotCopyNumGenomeWide(xyDatasetCopyNumber,      lohcateDirs.getSubDirPath(SubdirsDefault.Plots_CopyNumber)     + File.separator + sampleNameRoot + ".CopyNumber_GenomeWide", sampleNameRoot);				
		}
		out.close();
				
		// Print matrix breakdown statistics (distribtion of test events for each truth event)
		printContingencyTableInfo("SIMOUT", simOutputStream, sampleIndex, eventTables, testCounts);
		
		oneSampleData.clear();
	}

	// ========================================================================
	private static void printContingencyTableInfo(String prefix, PrintStream out, int sampleIndex, EnumMapSafe<EventType, ContingencyTable> eventTables, BucketCounterEnumMatrix<EventType> countBreakdownPerEvent) {
		for (EventType eventTruth : EventType.values()) {
			ContingencyTable table = eventTables.get(eventTruth);			

			// Print initial columns
			out.printf("%s:\t%d\t%12s\t|", prefix, sampleIndex, eventTruth.name());

			// Print contingency table raw values
			for (ContingencyTable.ContingencyTableValue ctValue : ContingencyTable.ContingencyTableValue.values()) {
				out.printf("\t%d", table.getCount(ctValue));
			}
			out.printf("\t%d\t%d", table.getTotalTruthPositive(), table.getTotalTruthNegative());

			// Print the contingency table calcuated values
			out.printf("\tSensitivity:\t%g\tSpecificity:\t%g\tF-Measure:\t%g", table.getSensitivity(), table.getSpecificity(), table.getFMeasure());

			// Now print the distributed counts 
			out.printf("\t|");
			for (EventType eventTest : EventType.values()) {
				out.printf("\t%d", countBreakdownPerEvent.getCount(eventTruth, eventTest)); 						
			}
			out.println("");
			out.flush();
		}
	}

	// ========================================================================
	// ========================================================================
	private static DefaultXYDataset classifySitesHelper_createAndFillXYData(
			EnumMapSafe<EventType, ParallelArrayDouble> coordinatesByEvent, 
			double[][] boundaryArrays) {
		
		DefaultXYDataset xyDataset = new DefaultXYDataset();
		
		for (EventType eventType : EventType.values()) {
			xyDataset.addSeries(eventType.name(), coordinatesByEvent.get(eventType).mArray);			
		}
		
		if (boundaryArrays != null) {
			xyDataset.addSeries("Boundary", boundaryArrays);
		}
		
		return xyDataset;
	}
	
	// ========================================================================
	private static<T extends Enum<T>> EnumMapSafe<T, ParallelArrayDouble> getCoordinateArraysPerClusterType(BucketCounterEnum<T> enumCounts, Class<T> enumClass) {
		
		EnumMapSafe<T, ParallelArrayDouble> theMap = new EnumMapSafe<T, ParallelArrayDouble>(enumClass);
		T[] enumConstants = enumClass.getEnumConstants();
		for (T enumConstant : enumConstants) {
			ParallelArrayDouble pArray = new ParallelArrayDouble(enumCounts.getCount(enumConstant));
			theMap.put(enumConstant, pArray);
		}
		
		return theMap;
	}

	// ========================================================================
	private static void calculateAdjustedVAFs(SiteList<ClusteringInputOneSite> sites,
											  ClusteringInputOneSampleMetaData metaData,
											  AllelicBiasTable allelicBiasTable, 											  
											  SeqPlatform platform) {
		
		final float defaultVAFNormal = 0.50f;
			
		for (int row = 0; row < sites.getNumSites(); row++) {
			ClusteringInputOneSite oneSiteInfo = sites.getSiteAtIndex(row);
			
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
	public static void getPValuesImbalance(SiteList<ClusteringInputOneSite> sites, ClusteringInputOneSampleMetaData metaData) { 
		
		for (int row = 0; row < sites.getNumSites(); row++) {
			ClusteringInputOneSite oneSiteInfo = sites.getSiteAtIndex(row);
												
			int varCovgTumor =  (int) Math.round(oneSiteInfo.mCovgTotalTumor  * metaData.mAdjustedVAFTumor[row]);
			int varCovgNormal = (int) Math.round(oneSiteInfo.mCovgTotalNormal * metaData.mAdjustedVAFNormal[row]);
			
			metaData.mImbalancePValuesTumor[row]  = GenotypeUtils.getPValuesImbalanceTissue(oneSiteInfo.mCovgTotalTumor,  varCovgTumor);
			metaData.mImbalancePValuesNormal[row] = GenotypeUtils.getPValuesImbalanceTissue(oneSiteInfo.mCovgTotalNormal, varCovgNormal);		
		}
	}

	// ========================================================================
	private static void smoothCopyNumbers(SiteList<ClusteringInputOneSite> sites, ClusteringInputOneSampleMetaData metaData) {
		
		ArrayList<CopyNumberRegionRange> setOfRegions = new ArrayList<CopyNumberRegionRange>();
		int positionDiffThreshold = 5_000_000;
		
		for (Chrom chrom : Chrom.values()) {
			if (chrom.isInvalid()) continue;
		
			setOfRegions.clear();
			ListIterator<CopyNumberRegionRange> regionsOnChrom = metaData.mGeneRegions.getIteratorForChromosome(chrom);
			
			while (regionsOnChrom.hasNext()) {
				CopyNumberRegionRange currentRegion = regionsOnChrom.next();
				
				if (!setOfRegions.isEmpty()) {
					int positionDiff = sites.getSiteAtIndex(currentRegion.getRangeStart()).getPosition() - 
							           sites.getSiteAtIndex(setOfRegions.get(0).getRangeStart()).getPosition();
					
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
	public static float[] calcRoughCopyNumberRatioPerSite(SiteList<ClusteringInputOneSite> sites, ClusteringInputOneSampleMetaData metaData) {		
		
		ObjectWalkerTracker<String> prevGene = new ObjectWalkerTracker<String>("", String.CASE_INSENSITIVE_ORDER);
		int rowOfFirstInstanceOfGene = -1;
		int positionPrev = 0;
		int numSites = sites.getNumSites();
		
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
		
		Arrays.fill(metaData.mTumorCopyNumRatiosPerGene, Regions.DefaultTumorNormalRatio);  // Initialize
		
		for (int row = 0; row < numSites; row++) {
			ClusteringInputOneSite oneSiteInfo = sites.getSiteAtIndex(row);
			int positionDifference = oneSiteInfo.getPosition() - positionPrev;
			positionPrev = oneSiteInfo.getPosition();

			metaData.mIsSomaticSite[row] = (oneSiteInfo.mCovgTotalNormal <= 0);
			float ratio = metaData.mIsSomaticSite[row] ? Regions.TumorNormalRatioOfSomaticSite : ((float) oneSiteInfo.mCovgTotalTumor / (float) oneSiteInfo.mCovgTotalNormal);						
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
					float averageRatio = Regions.TumorNormalRatioOfSomaticSite;
					
					if (numRowsWithSameGene == 0) {
						CompareUtils.ensureTrue(metaData.mIsSomaticSite[row - 1], "ERROR: Gene must > 0 sites representative of that gene!");
						
					} else if (numRowsWithSameGene == 1 && !AllowGeneWithOneSite) {					
						averageRatio = Regions.DefaultTumorNormalRatio;
						
					} else if (numRowsWithSameGene > 1 || AllowGeneWithOneSite) {
						double readCountAverageTumor  = (double) readCountSumTumor  / (double) numRowsWithSameGene;
						double readCountAverageNormal = (double) readCountSumNormal / (double) numRowsWithSameGene;
						int indexDistTumor  = NumberUtils.getIndexOfPoissonWithHighestProbability((int) Math.round(readCountAverageTumor),  pdTumor);
						int indexDistNormal = NumberUtils.getIndexOfPoissonWithHighestProbability((int) Math.round(readCountAverageNormal), pdNormal);
						//averageRatio = (float) (indexDistTumor + 1) / (float) (indexDistNormal + 1);
						
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
				geneRegion = new CopyNumberRegionRange(EventType.Ignored, oneSiteInfo.getChrom(), row);
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
		performTumorCopyNumRatioOffset(metaData, null);
		return metaData.mTumorCopyNumRatiosPerGene;
	}

	// ========================================================================
	private static void performTumorCopyNumRatioOffset(ClusteringInputOneSampleMetaData metaData, BitSet siteIsHet) {
		int numSites = metaData.mTumorCopyNumRatiosPerGene.length;
		DescriptiveStatistics ds = new DescriptiveStatistics();
		
		for (int row = 0; row < numSites; row++) {
			if (siteIsHet == null || siteIsHet.get(row)) {
				ds.addValue(metaData.mTumorCopyNumRatiosPerGene[row]);
			}
		}
		double median = ds.getPercentile(50);
		System.out.println("Median: " + median);
		double diff = 1.0 - median;
		
		for (int row = 0; row < numSites; row++) {
			metaData.mTumorCopyNumRatiosPerGene[row] += (float) diff;
		}
		
	}
	
	// ========================================================================
	/** Tumor : Normal avg coverage ratio per chromosome. */
	public static float[] calcAvgCoverageRatioPerChrom(ClusteringInputOneSample oneSampleInfo, ClusteringInputOneSampleMetaData metaData) {		
		
		float[] tumorNormalRatioPerChrom = new float[Chrom.values().length];
		int[] totalReadCountPerChromNormal = new int[Chrom.values().length];
		int[] totalReadCountPerChromTumor  = new int[Chrom.values().length];
		
		Arrays.fill(tumorNormalRatioPerChrom, Regions.DefaultTumorNormalRatio);  // Fill with 1.0 since that's an equal ratio
		Arrays.fill(metaData.mCopyNumRatioPerChromNormal, Regions.DefaultNormalRatio);
		
		Arrays.fill(metaData.mNumSitesPerChrom, 0);  // initialize counts
		Arrays.fill(metaData.mAvgReadCountPerChromNormal, 0); // initialize
		Arrays.fill(totalReadCountPerChromNormal, 0);
		Arrays.fill(totalReadCountPerChromTumor,  0);
		
		metaData.mReadCountTalliesTumor.clear();
		metaData.mReadCountTalliesNormal.clear();
		
		int numSites = oneSampleInfo.getNumSites();
		for (int row = 0; row < numSites; row++) {
			ClusteringInputOneSite oneSiteInfo = oneSampleInfo.getSiteAtIndex(row);
			
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
	public static int determineCopyGainVAFMaxLikelihood(ClusteringInputOneSample oneSampleInfo, ClusteringInputOneSampleMetaData metaData, RegionRange range, ArrayList<DoubleArrayList> listOfListOfMeans, TissueType targetTissue,  PrimitiveWrapper.WDouble probFromMaxLikelihood) {
		Chrom chrom = range.getChromosome();
		if (chrom.isInvalid()) return -1;  // Return if not on a valid chromosome

		int indexChromStart = oneSampleInfo.getIndex(chrom, range.getRangeStart());
		int indexChromEnd   = oneSampleInfo.getIndex(chrom, range.getRangeEnd());		

		indexChromStart = Math.max(0, indexChromStart);
		indexChromEnd = Math.min(indexChromEnd, oneSampleInfo.getNumSites() - 1);
		
		// Construct our lists
		int multiplier = 100;
		ArrayList<PoissonDistributionList> pdListofLists = new ArrayList<PoissonDistributionList>();				
		for (DoubleArrayList doubleArray : listOfListOfMeans) {
			PoissonDistributionList pdList = new PoissonDistributionList();
			pdListofLists.add(pdList);
			for (DoubleCursor d : doubleArray) {
				pdList.registerMean(multiplier * d.value);
				//System.out.println("\tPD: " + d.value);
			}			
			//System.out.println("--");			
		}
		
		double[] prob = new double[pdListofLists.size()];		
		Arrays.fill(prob, 0);
		for (int row = indexChromStart; row <= indexChromEnd; row++) {
			ClusteringInputOneSite oneSiteInfo = oneSampleInfo.getSiteAtIndex(row);				
			
			// Skip homozygous or near-homozygous sites
			if (oneSiteInfo.refOrVarHasZeroReadCount()) continue;  // skip homozygous sites
			if (!metaData.mVAFNormalHetRange.inRangeLowerExclusive(oneSiteInfo.calcVAFNormal())) continue;  // skip almost homozygous sites

			// Get the relevant variant allele fraction
			double vafTissue = (targetTissue == TissueType.Normal) ? oneSiteInfo.calcVAFNormal() : metaData.mAdjustedVAFTumor[row];
			
			PrimitiveWrapper.WDouble probVar = new PrimitiveWrapper.WDouble(0);
			for (int indexDist = 0; indexDist < pdListofLists.size(); indexDist++) {
				PoissonDistributionList pdList = pdListofLists.get(indexDist);
				pdList.getMeanOfMostLikelyDistribution((int) Math.round(multiplier * vafTissue), probVar);					
				prob[indexDist] += probVar.mDouble;
			}
		}

		int indexMaxElement = ArrayUtils.getIndexOfMaxElement(prob, 0);
		if (indexMaxElement >= 0) {
			probFromMaxLikelihood.mDouble = prob[indexMaxElement];
		} else {
			probFromMaxLikelihood.mDouble = 0;
			(new Exception("WARNING: Cannot have negative index!")).printStackTrace(); 			
		}
		
		if (indexMaxElement >= 0) {
//			for (int i = 0; i < prob.length; i++) {
//				System.out.println("\t" + prob[i]);
//			}
		}
		
		return indexMaxElement;
	}
	
	
	// ========================================================================
	public static void determineGermlineAneuploidiesByVAF(ClusteringInputOneSample oneSampleInfo, ClusteringInputOneSampleMetaData metaData) {
		String listOfListOfMeansStr = "({0.5};{0.6667,0.3333})";
		ArrayList<DoubleArrayList> listOfListOfMeans = ArrayUtils.getListOfDoubleListsFromStringForm(listOfListOfMeansStr, true);
		PrimitiveWrapper.WDouble probFromMaxLikelihood = new PrimitiveWrapper.WDouble(0);
		
		for (Chrom chrom : Chrom.values()) {
			if (chrom.isInvalid()) continue;  // move on if not on a valid chromosome
			
			int indexChromStart = oneSampleInfo.getIndexChromStart(chrom);
			if (indexChromStart < 0) continue;

			int indexChromEnd = oneSampleInfo.getIndexChromEnd(chrom);
			CompareUtils.ensureTrue(indexChromEnd >= 0, "ERROR: Chromosome must have last position!");
						
			CopyNumberRegionRange range = new CopyNumberRegionRange(EventType.Ignored, chrom, oneSampleInfo.getSiteAtIndex(indexChromStart).getPosition(), oneSampleInfo.getSiteAtIndex(indexChromEnd).getPosition());
			int indexOfMostLikelyList = determineCopyGainVAFMaxLikelihood(oneSampleInfo, metaData, range, listOfListOfMeans, TissueType.Normal, probFromMaxLikelihood);
			
			
			metaData.mChromHasGermlineGain[chrom.ordinal()] = (indexOfMostLikelyList > 0);
			//metaData.mChromHasGermlineGain[chrom.ordinal()] = fillRegionBasedOnVAFMaxLikelihood(range, oneSampleInfo, metaData, null, ClusterType.GainGermline, probFromMaxLikelihood, false);
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
				ClusteringInputOneSite oneSiteInfo = oneSampleInfo.getSiteAtIndex(row);						
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
	public static ClusteringResults<EventType> getClusters_withPlane(
			SiteList<ClusteringInputOneSite> sites, 
			ClusteringInputOneSampleMetaData metaData,
			RangeDouble vafNormalRange,
			String outFilenameFullPath, 
			int startingRowGermlineOrSomaticOrAll, 
			SeqPlatform platform) {

		//apply DBScan to points within NAF frame
		ElementPlaneSplit<Floint> planeSplit = new ElementPlaneSplit<Floint>(sites.getNumSites(), 2);

		getValidPointsForClustering(sites, planeSplit, metaData, ScalingFactor, platform, vafNormalRange);		
		
		int indexLowerPlane = 0, indexUpperPlane = 1;		
		ArrayList<Floint> pointsLowerPlane = planeSplit.getPointsOnPlane(indexLowerPlane);
		ArrayList<Floint> pointsUpperPlane = planeSplit.getPointsOnPlane(indexUpperPlane);
		
		System.out.printf("\tPoints Lower Plane: %d\n", pointsLowerPlane.size());
		System.out.printf("\tPoints Upper Plane: %d\n", pointsUpperPlane.size());
		System.out.println("\tBegin clustering algorithm: " + (new Date()).toString());

		//DBScanFaster dbscanner = new DBScanFaster(points, HET_BALL_EPS + DUP_WEDGE_LASSO, DUP_WEDGE_MINPTS, 0, 0, 1, 1); //parameters for capturing HET ball//KDBSCAN(param, 10, 0.0325f);
					
		ClusteringResults<EventType> clusterResults = new ClusteringResults<EventType>(sites.getNumSites());
		clusterResults.initializeResults(sites.getNumSites(), EventType.Ignored);
		
		// ----------------------- Scan the lower plane
		DBScanFaster dbscannerLowerPlane = new DBScanFaster(pointsLowerPlane, Clustering.HET_BALL_EPS, Clustering.HET_BALL_MINPTS, 0, 0, 1, 1); //parameters for capturing HET ball//KDBSCAN(param, 10, 0.0325f);
		dbscannerLowerPlane.cluster();

		// Now, we find the central heterozygous cluster ID and replace it with a fixed, known number.   
		// In parallel, we track which points in the lower plane were not part of the central cluster.
		// These points are clustered in another pass with lower neighbor density/size parameters.
		ArrayList<Floint> nonHetPoints = new ArrayList<Floint>(pointsLowerPlane.size());
		boolean[] isNonHetPoint                  = new boolean[pointsLowerPlane.size()];				
		int[] clusterAssignmentsLowerPlane = dbscannerLowerPlane.getClustAssignments();  // save and cache	
		final int clusterIDofHetBall = dbscannerLowerPlane.getCentralClusterID();
		
		for (int i = 0; i < clusterAssignmentsLowerPlane.length; i++) {
			int indexInMainList = planeSplit.getIndexOfPlaneElementInMainList(indexLowerPlane, i);
			
			if (clusterAssignmentsLowerPlane[i] == DBSCAN2.ClusterIDOfNoise) {
				clusterResults.setClassification(indexInMainList, EventType.Noise, DBSCAN2.ClusterIDOfNoise);
			} else {			
				isNonHetPoint[i] = (clusterAssignmentsLowerPlane[i] != clusterIDofHetBall); 
				if (isNonHetPoint[i]) {
					clusterResults.setClassification(indexInMainList, EventType.Ignored, clusterAssignmentsLowerPlane[i]);
					//nonHetPoints.add(pointsLowerPlane.get(i));				
				} else {
					clusterResults.setClassification(indexInMainList, EventType.HETGermline, clusterIDofHetBall);					
				}				
			}
		}
		
		System.out.println("\tPoints Lower Plane Non-Het: " + nonHetPoints.size());
		
		// Perform copy number normalization using supposedly diploid hets
		int numSites = sites.getNumSites();
		BitSet isHetSite = new BitSet(numSites);
		for (int row = 0; row < numSites; row++) {
			if (clusterResults.getClassification(row) == EventType.HETGermline) {
				isHetSite.set(row);
			}
		}
		performTumorCopyNumRatioOffset(metaData, isHetSite);
		
		
		//DBScanFaster dbscannerLowerPlaneNonHet = new DBScanFaster(nonHetPoints, Clustering.NON_HET_BALL_EPS, Clustering.NON_HET_BALL_MINPTS, 0, 0, 1, 1);
		//dbscannerLowerPlaneNonHet.cluster();
		//int[] clusterAssignmentsNonHet = dbscannerLowerPlaneNonHet.getClustAssignments();

		// -------------------- UPPER PLANE
		DBScanFaster dbscannerUpperPlane = new DBScanFaster(pointsUpperPlane, Clustering.NON_HET_BALL_EPS, Clustering.NON_HET_BALL_MINPTS, 0, 0, 1, 1);
		dbscannerUpperPlane.cluster();
		int[] clusterAssignmentsUpperPlane = dbscannerUpperPlane.getClustAssignments();
		
		for (int ind = 0; ind < clusterAssignmentsUpperPlane.length; ind++) {	
			int indexIntoOriginalRows = planeSplit.getIndexOfPlaneElementInMainList(indexUpperPlane, ind); 					
			float copyNum = metaData.mTumorCopyNumRatiosPerGene[indexIntoOriginalRows] * Regions.DefaultDiploidCopyNumber;
			
			if (clusterAssignmentsUpperPlane[ind] == DBSCAN2.ClusterIDOfNoise) {
				clusterResults.setClassification(indexIntoOriginalRows, EventType.Noise, DBSCAN2.ClusterIDOfNoise);
			} else {				
				if (isCopyNumAmplified(copyNum)) {
					clusterResults.setClassification(indexIntoOriginalRows, EventType.GainSomatic, clusterAssignmentsUpperPlane[ind]);					
				} else if (isCopyNumInDiploidRange(copyNum)) {					
					clusterResults.setClassification(indexIntoOriginalRows, EventType.cnLOH, clusterAssignmentsUpperPlane[ind]);
				} else {
					clusterResults.setClassification(indexIntoOriginalRows, EventType.LOH, clusterAssignmentsUpperPlane[ind]);
				}
			}
			
			ClusteringInputOneSite oneSiteInfo = sites.getSiteAtIndex(indexIntoOriginalRows);
			Chrom chrom = oneSiteInfo.getChrom();
			if (metaData.chromHasGermlineGain(chrom)) {
				clusterResults.setClassification(indexIntoOriginalRows, EventType.GainGermline, clusterAssignmentsUpperPlane[ind]);				
			} else {
				if (ForcePointsOnDiagonalAsNull) {
					if (pointOnDiagonal(pointsUpperPlane.get(ind), ClusterDiagonalLeeway) && isCopyNumInDiploidRange(copyNum)) {						
						clusterResults.setClassification(indexIntoOriginalRows, EventType.Ignored, clusterAssignmentsUpperPlane[ind]);
					}
				}
				
			}
		}
		
		System.out.println("\tEnd clustering algorithm: " + (new Date()).toString());
		
		assignClustersToNonClusteredSites(sites, platform, vafNormalRange, clusterResults);
		return clusterResults;
	}

	// ========================================================================
	// Convenience private function to break up caller function
	private static void assignClustersToNonClusteredSites(SiteList<ClusteringInputOneSite> sites, 
			SeqPlatform platform,
			//int startingRowSomatic,					
			RangeDouble vafNormalRange,
			ClusteringResults<EventType> clusterResults) {
		
		for (int row = 0; row < sites.getNumSites(); row++) {

			ClusteringInputOneSite oneSiteInfo = sites.getSiteAtIndex(row);
			float vafNormal = oneSiteInfo.calcVAFNormal();
			boolean vafInRangeNormal = vafNormalRange.inRangeLowerExclusive(vafNormal); 								
			
			if (vafInRangeNormal) {
				// points already classified in this case
				continue; 
				
			} else {
				EventType eventAtSite = null;
				
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
						eventAtSite = EventType.HETSomatic;
					} else if (vafTumor <= hetBoundaryLower) {
						// Normal: AA, Tumor: AA [Thus homozygous reference in both, no events]
						eventAtSite = EventType.Ignored;
					} else if (vafTumor > hetBoundaryUpper) {
						// Normal: AA, Tumor: BB or CC [made by: AA -> AB or AC (somatic het mutation) -> B or C (LOH, loss of A)]
						eventAtSite = EventType.LOH;
					} else {
						// Normal: AA, Tumor: AB [made by: AA -> AB (somatic het mutation)
						eventAtSite = EventType.Ignored;
					}					
				} else if (RangeDouble.inRangeLowerExclusive(vafNormal, vafNormalRange.getBoundUpper(), Regions.MaxVariantAlleleFrequency)) {
					// We are above the upper frame boundary
					if (vafNormal >= (1.0 - 0.05) && vafTumor >= 0.10 && vafTumor <= 0.90) {
						eventAtSite = EventType.HETSomatic;
					} else if (vafTumor <= hetBoundaryLower) {
						// Normal: BB, Tumor: AA [made by: BB -> AB (reverse somatic het mutation) -> A (LOH, loss of B)]
						eventAtSite = EventType.LOH;
					} else if (vafTumor > hetBoundaryUpper) {
						// Normal: BB, Tumor: BB or CC (ambiguous until we know exact variant for tumor)
						// TODO - Leave as Null for now, but will need to change later to resolve the
						// ambiguity mentioned above
						eventAtSite = EventType.Ignored;
					} else {
						// Normal: BB, Tumor: AB or CB [made by: BB -> AB (reverse somatic het mutation) or BB -> CB (somatic het mutation)
						eventAtSite = EventType.Ignored;
					}
				} else {
					CompareUtils.throwErrorAndExit("ERROR: Contradiction - variant allele frequency cannot be in and out of bounds simultanteously!" + vafNormal);
				}
				
				// Now finally assign the cluster.
				CompareUtils.ensureNotNull(eventAtSite, "ERROR: Event type for point cannot be unassigned!");
				clusterResults.setClassification(row, eventAtSite, ClusteringResults.InvalidSubClusteringID);				
			}			
		}		
	}

	// ========================================================================	
	private static void getValidPointsForClustering(SiteList<ClusteringInputOneSite> sites, // input
													ElementPlaneSplit<Floint> planeSplit,
													ClusteringInputOneSampleMetaData metaData,
													double scalingFactor,
													SeqPlatform platform, 
													RangeDouble vafNormalRange) {
		
		planeSplit.clear();
		for (int row = 0; row < sites.getNumSites(); row++) {
			ClusteringInputOneSite oneSiteInfo = sites.getSiteAtIndex(row);			
						
			if (vafNormalRange.inRangeLowerExclusive(oneSiteInfo.calcVAFNormal())) {				
				Floint thePoint = new FlointImpl(metaData.mAdjustedVAFTumor[row], metaData.mAdjustedVAFNormal[row], (float) (metaData.mImbalancePValuesTumor[row] * scalingFactor));
				
				// Assume p-value as vertical row factor
				boolean  tumorSigImbalanced = (metaData.mImbalancePValuesTumor[row] <= metaData.mSigPValuesPerSite[row]);
		//      boolean  tumorSigImbalanced = (metaData.mImbalancePValuesTumor[row] <= metaData.mFDRTumor);
				boolean normalSigImbalanced = (metaData.chromHasGermlineGain(oneSiteInfo.getChrom()));
				int planeID = (tumorSigImbalanced || normalSigImbalanced) ? 1 : 0; 
				planeSplit.registerElement(planeID, thePoint, row);
			} else {
				planeSplit.registerElement(ElementPlaneSplit.InvalidPlaneID, null, row);				
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
		float copyNum = copyNumRatio * Regions.DefaultDiploidCopyNumber;
		return isCopyNumInDiploidRange(copyNum);
	}

	// ========================================================================	
	// Tests whether the point lies roughly on the x = y diagonal, given a leeway from a ratio of 1.0
	private static boolean pointOnDiagonal(Floint point, double leeway) {
		double ratio = point.getY() / point.getX();
		double equalRatio = 1.0;
		leeway = Math.min(leeway, equalRatio);
		return ( ((equalRatio - leeway) <= ratio) && (ratio <= (equalRatio + leeway)) );
	}

	// ========================================================================
	// Extracts and calculates the variant allele frequency depending on the platform 
	static float extractVAFNormal(String[] components, SeqPlatform platform) {
		switch(platform) {
		case Illumina: return  Float.parseFloat(components[Regions.Col_NAFTAFInput_VariantRatioNormal]);
		case SOLiD:    return (Float.parseFloat(components[6]) / Float.parseFloat(components[5]));
		}
		return Float.NaN;
	}

	// ========================================================================
	// Extracts and calculates the variant allele frequency in the normal depending on the platform tab-delimited file
	static float extractVAFNormal(String line, SeqPlatform platform) {
		switch(platform) {
		case Illumina: return  Float.parseFloat(StringUtils.extractNthColumnValue(line, Regions.Col_NAFTAFInput_VariantRatioNormal, StringUtils.TabStr));
		case SOLiD:    return (Float.parseFloat(StringUtils.extractNthColumnValue(line, 6, StringUtils.TabStr)) / 
				               Float.parseFloat(StringUtils.extractNthColumnValue(line, 5, StringUtils.TabStr)));
		}
		return Float.NaN;		
	}

	// ========================================================================
	static float extractVAFTumor(String[] components, SeqPlatform platform) {
		switch(platform) {
		case Illumina: return  Float.parseFloat(components[Regions.Col_NAFTAFInput_VariantRatioTumor]);
		case SOLiD:    return (Float.parseFloat(components[4]) / Float.parseFloat(components[3]));
		}
		return Float.NaN;
	}

	// ========================================================================
	// Extracts and calculates the variant allele frequency in the tumor depending on the platform tab-delimited file
	static float extractVAFTumor(String line, SeqPlatform platform) {
		switch(platform) {
		case Illumina: return  Float.parseFloat(StringUtils.extractNthColumnValue(line, Regions.Col_NAFTAFInput_VariantRatioTumor, StringUtils.TabStr));
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
