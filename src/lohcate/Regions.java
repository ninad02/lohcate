package lohcate;
import genomeEnums.Chrom;
import genomeEnums.Nuc;
import genomeEnums.VariantLocation;
import genomeUtils.ChromPositionTracker;
import genomeUtils.GenotypeUtils;
import genomeUtils.RegionRange;
import genomeUtils.RegionRange.RegionRangeOverlap;
import genomeUtils.SNVMap;
import static genomeUtils.RegionRange.RegionRangeOverlap.*;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Date;
import java.util.ListIterator;

import nutils.ArgumentParserUtils;
import nutils.ArrayUtils;
import nutils.CompareUtils;
import nutils.EnumMapSafe;
import nutils.IOUtils;
import nutils.NullaryClassFactory;
import nutils.PrimitiveWrapper;
import nutils.StringUtils;
import nutils.counter.BucketCounterEnum;
import nutils.counter.DynamicBucketCounter;

import com.carrotsearch.hppc.IntArrayList;
import com.martiansoftware.jsap.FlaggedOption;
import com.martiansoftware.jsap.JSAP;
import com.martiansoftware.jsap.JSAPResult;
import com.martiansoftware.jsap.StringParser;

import lohcate.LOHcate.SubdirsDefault;
import lohcate.clustering.Clustering;
import lohcate.clustering.ClusteringInputOneSample;
import lohcate.clustering.ClusteringInputOneSampleMetaData;
import lohcate.clustering.ClusteringInputOneSite;
import lohcate.clustering.ClusteringParams;
import lohcate.clustering.ClusteringPlotting;
import lohcateEnums.EventType;
import lohcateEnums.ColorPastel;
import lohcateEnums.MutationType;
import lohcateEnums.SeqPlatform;

/**
 * LOHcate --- A software tool for LOH calling and visualization in cancer genomes
 * D Wheeler & SG Reddy
 * Human Genome Sequencing Center, Baylor College of Medicine (Houston, TX)
 * 
 * This runner class is the big enchilada. See main method for 'switchboard' / pipeline overview.
 * 
 * @author Ninad Dewal, Siddharth G. Reddy
 *
 */
public class Regions {
	
	public static final int DefaultDiploidCopyNumber = 2;
	
	public static final int REGION_SEGMENTATION_DIST_THRESHOLD = 2000000; //greatest possible distance between 2 'adjacent' points of LOH in a region of 'contiguous' LOH
	
	public static final String GermlineSuffix = "." + VariantLocation.Germline.toLowerCase();	
	public static final String NovelStr  = "novel";	
	public static final String ChromPrefix = Chrom.ChromPrefix_chr;
	public static final String MissingGeneNameValue = ".";
	public static final char   MissingAllele = '.';
	public static final float MaxVariantAlleleFrequency = 1.0f;
	public static final String GenBrowserTrack = ".genBrowserTrack";
	
	public static final boolean EliminateExtremeGCSites = true;
	public static final boolean EliminateHighDensitySNVs = true;
	
	public static final PrintStream LogOutput = IOUtils.getPrintStream("LOHcate.Log." + ((new Date()).toString()).replace(' ', '_').replace(':', '-') + ".txt");
	
	/*
	0 1   refName              chr1
	1 2   coordinate           883918
	2 3   refBase              G
	3 4   varBase              A
	4 5   variantBase-N        A
	5 6   variantBase-T        A
	6 7   refEnv-N             GATCACGGAGAAG
	7 8   refEnv-T             GATCACGGAGAAG
	
	8    Q20_TotCov_N         55
	9  Q20_TotCov_T         74
	10  Q20_VarCov_N         32
	11  Q20_VarCov_T         44
	12  Q20_VariantRatio_N   0.581818
	13  Q20_VariantRatio_T   0.594595
	14  dbsnp                rs139116730,byfrequency,G|0.997;A|0.003,.
	15  MutationType         synonymous_SNV
	16  Hugo_Symbol          NOC2L
	 */
	// Column constants for the curated TSV files (files that have a cluster column)
	public static final int Col_NAFTAFInput_Chrom                 = 0;
	public static final int Col_NAFTAFInput_Position              = 1;
	public static final int Col_NAFTAFInput_AlleleRef             = 2;
	public static final int Col_NAFTAFInput_AlleleVarPop          = 3;
	public static final int Col_NAFTAFInput_AlleleVarN            = 4;
	public static final int Col_NAFTAFInput_AlleleVarT            = 5;
	public static final int Col_NAFTAFInput_FlankingStringNormal  = 6;
	public static final int Col_NAFTAFInput_FlankingStringTumor   = 7;
	public static final int Col_NAFTAFInput_TotalCoverageNormal   = 8;
	public static final int Col_NAFTAFInput_TotalCoverageTumor    = 9;
	public static final int Col_NAFTAFInput_VariantCoverageNormal = 10;
	public static final int Col_NAFTAFInput_VariantCoverageTumor  = 11;
	public static final int Col_NAFTAFInput_VariantRatioNormal    = 12;
	public static final int Col_NAFTAFInput_VariantRatioTumor     = 13;
	public static final int Col_NAFTAFInput_DbSNPString           = 14;
	public static final int Col_NAFTAFInput_MutationType          = 15;
	public static final int Col_NAFTAFInput_HugoSymbol            = 16;
	
	private static final double GCContentThresholdLow  = 0.05;
	private static final double GCContentThresholdHigh = 0.80;
	
	private static final int MaxWindowLength = 1000;
	private static final int MaxSitesInWindowAllowed = 3;
	
	// Column constants for the curated TSV files (files that have a cluster column)
	private static final int ColCuratedTSV_Chrom    = Col_NAFTAFInput_Chrom;
	private static final int ColCuratedTSV_Position = Col_NAFTAFInput_Position;
	//private static final int ColCuratedTSV_VariantBaseTumor = 2;
	private static final int ColCuratedTSV_VafNormal    = Col_NAFTAFInput_VariantRatioNormal;
	private static final int ColCuratedTSV_VafTumor     = Col_NAFTAFInput_VariantRatioTumor;
	private static final int ColCuratedTSV_dbSNP        = Col_NAFTAFInput_DbSNPString;	
	private static final int ColCuratedTSV_MutationType = Col_NAFTAFInput_MutationType;
	private static final int ColCuratedTSV_Gene         = Col_NAFTAFInput_HugoSymbol;
	//private static final int ColCuratedTSV_VariantLocation = 7;
	private static final int ColCuratedTSV_Cluster      = 17; //8;
	
	public static final GenotypeUtils.GenomicCoordinateComparatorInTextFileLine LineComparatorTab = new GenotypeUtils.GenomicCoordinateComparatorInTextFileLine();
	public static final NullaryClassFactory<ArrayList> NullClassFactoryArrayList = new NullaryClassFactory<ArrayList>(ArrayList.class);
	
	// ========================================================================
	/** This removes the header lines and filters out sites. */
	public static ArrayList<String> curateSNPCalls_removeHeaderLinesFromRows(ArrayList<String> rows, boolean removeHighDensitySNVSites, boolean removeExtremeGCSites) {
		int windowPositionStart = 0;
		int windowRowStart = 0;
		Chrom chromPrev = Chrom.c0;
		int numSitesRemoved = 0;
		
		for (int rowIndex = 0; rowIndex < rows.size(); rowIndex++) {
			String row = rows.get(rowIndex);
			if (row.indexOf("refName") >= 0 && row.indexOf("coord") >= 0 || row.indexOf("chrX") >= 0 || row.indexOf("chrY") >= 0 || row.indexOf("chrM") >= 0 /*|| row.indexOf("intergenic") >= 0 */) {
				rows.set(rowIndex, null);
			} else {
				// Do some light GC filtering
				Chrom chrom = Chrom.getChrom(    StringUtils.extractNthColumnValue(row, Regions.Col_NAFTAFInput_Chrom,               StringUtils.FileExtensionTSV.mDelimiter));
				int position = Integer.parseInt( StringUtils.extractNthColumnValue(row, Regions.Col_NAFTAFInput_Position,            StringUtils.FileExtensionTSV.mDelimiter));
				String gcString =                StringUtils.extractNthColumnValue(row,        Col_NAFTAFInput_FlankingStringTumor, StringUtils.FileExtensionTSV.mDelimiter);
				
				int posDiff = position - windowPositionStart;
				
				if ((chrom != chromPrev) || (posDiff >= MaxWindowLength)) {
					int numSitesSpanned = (rowIndex - 1) - windowRowStart + 1;
					if (numSitesSpanned > MaxSitesInWindowAllowed) {
						for (int rowToClean = windowRowStart; rowToClean < rowIndex; rowToClean++) {
							//System.out.println("Removing:\t" + rows.get(rowToClean));
							if (removeHighDensitySNVSites) {
								rows.set(rowToClean, null);
								numSitesRemoved++;
							}							
						}
					}	

					windowPositionStart = position;
					windowRowStart = rowIndex;
					chromPrev = chrom;
				} 
				
				double fractionGCNormal = GenotypeUtils.calcFractionGC(gcString);
				if (removeExtremeGCSites && (fractionGCNormal < GCContentThresholdLow) || (fractionGCNormal >= GCContentThresholdHigh)) {
					rows.set(rowIndex, null);
				}
			}
		}
		System.out.println("\tNum Sites Removed in Windows: " + numSitesRemoved);
		return ArrayUtils.removeNullElements(rows);		
	}
	
	// ========================================================================
	public static final float DefaultTumorNormalRatio = 1.0f;
	public static final float DefaultNormalRatio = 1.0f;
	public static final float TumorNormalRatioOfSomaticSite = -1.0f;

	
	
	
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
		StringUtils.FileExtensionAndDelimiter fileExtDelim = StringUtils.FileExtensionTSV;
		File[] files = (new File(inDir)).listFiles();		
		ArrayList<File> validFilesList = new ArrayList<File>(files.length);
		SNVMap snvMap = new SNVMap();
		
		// Create our output directories
		System.out.println("Creating sub-directories...");
		IOUtils.createDirectoryPath(outDir, false);
		IOUtils.createDirectoryPath(outDirBrowserTracks, false);
		
		// First we determine the regions from each sample
		System.out.println("Constructing regions in samples...");
		ArrayList<CopyNumberRegionsByChromosome> regionsInSamples = new ArrayList<CopyNumberRegionsByChromosome>();
		for (File file : files) {
			CopyNumberRegionsByChromosome regionsInOneSample = segmentRegionsOneFile(file, outDir, snvMap);
			if (regionsInOneSample != null) {
				regionsInSamples.add(regionsInOneSample);
				validFilesList.add(file);
			}
		}
		
		// Print out the map
		snvMap.printMe(true, outDir + File.separator + "AllSites.txt");
		
		// Now we have all the contiguous regions from all the samples.  Find the regions of the cluster type		
		// Declare the maximum stretch of a region for a particular cluster type
		EnumMapSafe<EventType, Integer> maxBasePairsContiguousRegion = ArrayUtils.createEnumMap(EventType.class, Integer.MAX_VALUE);
		maxBasePairsContiguousRegion.put(EventType.GainSomatic, REGION_SEGMENTATION_DIST_THRESHOLD);
		maxBasePairsContiguousRegion.put(EventType.cnLOH,       REGION_SEGMENTATION_DIST_THRESHOLD);
		maxBasePairsContiguousRegion.put(EventType.LOH,         REGION_SEGMENTATION_DIST_THRESHOLD);
		
		EnumMapSafe<EventType, ArrayList<CopyNumberRegionsByChromosome>> regionsInSamplesPerEventType =
				ArrayUtils.createEnumMapOfArrayLists(EventType.class, CopyNumberRegionsByChromosome.class);
		
		for (EventType clusterType : EventType.AmpLOHHetG) { 						
			ArrayList<CopyNumberRegionsByChromosome> regionsInSamplesForOneClusterType = regionsInSamplesPerEventType.get(clusterType);	
			
			for (CopyNumberRegionsByChromosome regionsInOneSample : regionsInSamples) {
				CopyNumberRegionsByChromosome regionsInOneSampleMerged = mergeRegionsWithConstraintsOld(regionsInOneSample, clusterType, maxBasePairsContiguousRegion.get(clusterType));
				regionsInOneSampleMerged.removeSingletonRegions();
				regionsInSamplesForOneClusterType.add(regionsInOneSampleMerged);
				printSegmentedRegionsToFile(outDir, regionsInOneSampleMerged, clusterType, snvMap);
			}
		}

		PrintStream out = IOUtils.getPrintStream(outDir + File.separator + "testOut.txt");
		
		// Now we want to determine the recurrent regions
		EnumMapSafe<EventType, CopyNumberRegionsByChromosome> regionsRecurrentPerEventType = 
				new EnumMapSafe<EventType, CopyNumberRegionsByChromosome>(EventType.class);	
		
		for (EventType clusterType : EventType.AmpLOHHetG) { 					
			CopyNumberRegionsByChromosome recurrentRegionsForOneClusterType = determineRecurrentRegions(regionsInSamplesPerEventType.get(clusterType), clusterType);
			regionsRecurrentPerEventType.put(clusterType, recurrentRegionsForOneClusterType);
			recurrentRegionsForOneClusterType.print(out, fileExtDelim.mDelimiter);
			plotRecurrence(recurrentRegionsForOneClusterType, outDir, clusterType);
		}
		
		
		System.out.println("Scoring regions...");
		// Now, we need to go through the recurrent regions and score them, based on the 
		// contents of the individual samples that make up the recurrent regions.
		for (EventType clusterType : EventType.AmpLOHHetG) {
			ArrayList<CopyNumberRegionsByChromosome> regionsInSamplesForOneClusterType = regionsInSamplesPerEventType.get(clusterType);
			CopyNumberRegionsByChromosome recurrentRegionsForOneClusterType = regionsRecurrentPerEventType.get(clusterType);			
			
			countClusterTypesInRegions(recurrentRegionsForOneClusterType, regionsInSamplesForOneClusterType, validFilesList);
			recurrentRegionsForOneClusterType.print(out, fileExtDelim.mDelimiter);
		}
				
		IOUtils.closePrintStream(out);

		System.out.println("Generating Browser tracks...");
		// Generatebrowser tracks
		for (EventType clusterType : EventType.OnlyLOH) {
			//ClusterType clusterType = ClusterType.OnlyLOH[clusterIndex];
			ArrayList<CopyNumberRegionsByChromosome> regionsInSamplesForOneClusterType = regionsInSamplesPerEventType.get(clusterType);
			CopyNumberRegionsByChromosome recurrentRegionsForOneClusterType = regionsRecurrentPerEventType.get(clusterType);
			
			String outDirBrowserTracksForCluster = outDirBrowserTracks + File.separator + clusterType.name();
			IOUtils.createDirectoryPath(outDirBrowserTracksForCluster, false);
			genBrowserTracks(regionsInSamplesForOneClusterType, validFilesList, recurrentRegionsForOneClusterType, clusterType, outDirBrowserTracksForCluster);
		}				
	}
	
	// ========================================================================
	/** Plots a genome wide recurrence plot. */
	public static void plotRecurrence(CopyNumberRegionsByChromosome recurrentRegionsForOneClusterType, String outDir, EventType eventType) {
		
		EnumMapSafe<Chrom, DynamicBucketCounter> eventCount = DynamicBucketCounter.ClassFactory.newEnumMap(Chrom.class);
		
		float maxRecurrenceScore = -1.0f;
		for (Chrom chrom : Chrom.values()) {
			ArrayList<CopyNumberRegionRange> cnrrOnChrom = recurrentRegionsForOneClusterType.mRegionsByChrom.get(chrom);
			if (cnrrOnChrom.isEmpty()) continue;
			
			for (CopyNumberRegionRange cnrr : cnrrOnChrom) {
				DynamicBucketCounter dbc = eventCount.get(chrom);				
				
				int range = cnrr.getRangeLength();
				int rangeIncr = Math.max(range / 20, 1);
				for (int startPos = cnrr.getRangeStart(); startPos <= cnrr.getRangeEnd(); startPos += rangeIncr) {
					dbc.incrementCount(startPos, Math.round(cnrr.mRecurrenceScore));
				}
			}
		}
		
		ClusteringPlotting.plotRecurrenceGenomeWideByEvent(eventCount, outDir, eventType);
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
				CompareUtils.throwErrorAndExit("ERROR: Samples don't match up: " + regionsOneSample.mSampleName + "\t" + inFile.getName());
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
				
				final Chrom chrom  = Chrom.getChrom                       (StringUtils.extractNthColumnValue(line, ColCuratedTSV_Chrom,    StringUtils.TabStr));					
				final int position = Integer.parseInt                     (StringUtils.extractNthColumnValue(line, ColCuratedTSV_Position, StringUtils.TabStr));
				final EventType clusterType = EventType.getClusterType(StringUtils.extractNthColumnValue(line, ColCuratedTSV_Cluster,  StringUtils.TabStr));
				
				if ((prevChrom == null) || (chrom != prevChrom)) {
					regionsInChr =                regionsGenomeWide.mRegionsByChrom.get(chrom);
					regionsInChrSampleSpecific =   regionsOneSample.mRegionsByChrom.get(chrom);					
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
						region.mClusterTypeCounts.increment(clusterType);
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
	public static CopyNumberRegionsByChromosome determineRecurrentRegions(ArrayList<CopyNumberRegionsByChromosome> regionsInSamples, EventType clusterType) {
		
		// Test for a null parameter
		if (regionsInSamples.isEmpty()) return null;
		
		// Create a copy that serves as the target (or "sink") for all the intersecting regions
		CopyNumberRegionsByChromosome target = regionsInSamples.get(0).getCopy();
		
		// Now traverse the rest of the samples
		for (int i = 1; i < regionsInSamples.size(); i++) {
			takeUnionAndBreakDownIntersectingRegions(target, regionsInSamples.get(i), clusterType);
		}
		
		return target;
	}

	// ========================================================================
	/** Given the merged segmented regions in one sample, this method prints the regions to file. */
	public static void printSegmentedRegionsToFile(String outDir, CopyNumberRegionsByChromosome regionsInSample, EventType clusterType) {
				
		String outFilename = outDir + File.separator + "Regions.GISTIC." + clusterType.name() + "." + regionsInSample.mSampleName + ".txt";
		StringBuilder sb = new StringBuilder(2048);
		String delim = StringUtils.FileExtensionTSV.mDelimiter;
		
		
		double copyNumBase = (clusterType == EventType.GainSomatic) ? 2.5 : (clusterType == EventType.LOH ? 1.5 : 2.0);		
//		double log2RatioGisticCol = log2Ratio - 1.0;
//		double log2CopyNumGistic = log2CopyNum - 1.0;
		
		BufferedWriter out = IOUtils.getBufferedWriter(outFilename);
		// Go chromosome by chromosome
		for (Chrom chrom : Chrom.Autosomes) {							
			
			ArrayList<CopyNumberRegionRange> regionsInChrom = regionsInSample.mRegionsByChrom.get(chrom);			
			for (CopyNumberRegionRange cnrr : regionsInChrom) {
				sb.setLength(0);
				
				double errorFactor = Math.random() / 10;  // Get between 0 and 0.1
				boolean addPositiveError = (Math.random() >= 0.5);
				double copyNum = copyNumBase + ((addPositiveError) ? errorFactor : -errorFactor);
				double copyNumRatio = copyNum / DefaultDiploidCopyNumber;
				double log2CopyNum = Math.log10(copyNum) / Math.log10(2);
				double log2Ratio = Math.log10(copyNumRatio) / Math.log10(2);

				
				sb.append(regionsInSample.mSampleName)
				  .append(delim).append(chrom.ordinal())
				  .append(delim).append(cnrr.getRangeStart())
				  .append(delim).append(cnrr.getRangeEnd())
				  .append(delim).append(cnrr.getNumSitesInterrogated())
				  .append(delim).append(log2Ratio);
				
				IOUtils.writeToBufferedWriter(out, sb.toString(), true);
			}
		}
		
		IOUtils.closeBufferedWriter(out);
	}
	
	
	// ========================================================================
	// Printing JISTIC relevant output
	// ========================================================================
	public static void printSegmentedRegionsToFile_GISTIC(
			BufferedWriter outputWriter, 
			CopyNumberRegionsByChromosome regionsInSample, 
			ArrayList<EventType> events, 
			ClusteringInputOneSample oneSampleInfo, 
			ClusteringInputOneSampleMetaData metaData, 
			SNVMap snvMap) {
		
		
		
		StringBuilder sb = new StringBuilder(2048);
		String delim = StringUtils.FileExtensionTSV.mDelimiter;
		CopyNumberRegionRange dummyRange = new CopyNumberRegionRange(EventType.HETGermline, Chrom.c0, 0);
		
		
	}
	
	
	
	public static void printSegmentedRegionsToFile(String outDir, CopyNumberRegionsByChromosome regionsInSample, EventType clusterType, SNVMap snvMap) {
		
		String outFilename = outDir + File.separator + "Regions.GISTIC." + clusterType.name() + "." + regionsInSample.mSampleName + ".txt";
		StringBuilder sb = new StringBuilder(2048);
		String delim = StringUtils.FileExtensionTSV.mDelimiter;
		CopyNumberRegionRange dummyRange = new CopyNumberRegionRange(EventType.HETGermline, Chrom.c0, 0);		
		double dummyRangeCopyNumberBase = 2.0;
		
		double copyNumBase = (clusterType == EventType.GainSomatic) ? 2.5 : (clusterType == EventType.LOH ? 1.5 : 2.0);		
//		double log2RatioGisticCol = log2Ratio - 1.0;
//		double log2CopyNumGistic = log2CopyNum - 1.0;
		
		BufferedWriter out = IOUtils.getBufferedWriter(outFilename);
		// Go chromosome by chromosome
		for (Chrom chrom : Chrom.Autosomes) {				
			
			int numSitesOnChrom = snvMap.getNumSitesOnChromosome(chrom);
			if (numSitesOnChrom <= 0) continue;
			
			ListIterator<CopyNumberRegionRange> regionsInChromIter = regionsInSample.mRegionsByChrom.get(chrom).listIterator();
			CopyNumberRegionRange currentRegion = null;
			int indexInMap = 0;
			boolean dummyRangeValid = false;
						
			// First scan and see if we have any regions preceding the first chromosomal position in the map
			while (regionsInChromIter.hasNext()) {
				currentRegion = regionsInChromIter.next();
				if (currentRegion.afterRange(chrom, snvMap.getPosition(chrom, indexInMap))) {
					printSegmentedRegionsToFile_Helper(currentRegion, sb, delim, copyNumBase, regionsInSample.mSampleName, out);
					currentRegion = null;
				} else {
					break;
				}
			}
			
			// We could have broken the loop because the current region encloses the
			
			for (indexInMap = 0; indexInMap < numSitesOnChrom; indexInMap++) {
				int mapPosition = snvMap.getPosition(chrom, indexInMap);
				
				if ((currentRegion == null) || currentRegion.beforeRange(chrom, mapPosition)) {
					if (dummyRangeValid) {
						if (dummyRange.isFinalized()) {
							printSegmentedRegionsToFile_Helper(dummyRange, sb, delim, dummyRangeCopyNumberBase, regionsInSample.mSampleName, out);
							dummyRange.set(chrom, mapPosition, mapPosition, false, 1);
						} else {
							boolean extendResult = dummyRange.extendRange(chrom, mapPosition);
							if (!extendResult) {
								String errorString = "ERROR: printSegmentedRegionsToFile(): Could not extend range! " + chrom + "\t" + mapPosition + "\t" + dummyRange;
								CompareUtils.ensureTrue(extendResult, errorString);	
							}
														
						}
					} else {
						dummyRange.set(chrom, mapPosition, mapPosition, false, 1);
						dummyRangeValid = true;
					}
					
				} else if (currentRegion.inRange(chrom, mapPosition)) {
					if (dummyRangeValid) {
						printSegmentedRegionsToFile_Helper(dummyRange, sb, delim, dummyRangeCopyNumberBase, regionsInSample.mSampleName, out);
						dummyRangeValid = false;
					}
					
					// Print this region to file 
					int indexOfRegionStartInMap = snvMap.getIndexOfPositionInMap(chrom, currentRegion.getRangeStart());
					if (indexOfRegionStartInMap < 0) {
						String errorString = "ERROR: printSegmentedRegionsToFile(): Region start (" + currentRegion.getRangeStart() + "must exist in map!";
						CompareUtils.ensureTrue(false, errorString);
					}
										
					// Now change the loop index to move to the map position just after the end of this region
					int indexOfRegionEndInMap = snvMap.getIndexOfPositionInMap(chrom, currentRegion.getRangeEnd());
					if (indexOfRegionEndInMap < 0) {
						String errorString = "ERROR: printSegmentedRegionsToFile(): Region end (" + currentRegion.getRangeEnd() + "must exist in map!";
						CompareUtils.ensureTrue(false, errorString);						
					}
					indexInMap = indexOfRegionEndInMap;  // indexInMap will be incremented at loop end
					
					currentRegion.set(chrom, currentRegion.getRangeStart(), currentRegion.getRangeEnd(), true, indexOfRegionEndInMap - indexOfRegionStartInMap + 1);
					printSegmentedRegionsToFile_Helper(currentRegion, sb, delim, copyNumBase, regionsInSample.mSampleName, out);				
					
					// Now move to the next cna-affected region
					currentRegion = regionsInChromIter.hasNext() ? regionsInChromIter.next() : null;
					
				} else if (currentRegion.afterRange(chrom, mapPosition)) {
					CompareUtils.ensureTrue(false, "ERROR: printSegmentedRegionsToFile(): Impossible for region to precede position!");
				} else {
					CompareUtils.ensureTrue(false, "ERROR: printSegmentedRegionsToFile(): Impossible state!");
				}
			}
			
			if (dummyRangeValid) {
				printSegmentedRegionsToFile_Helper(dummyRange, sb, delim, dummyRangeCopyNumberBase, regionsInSample.mSampleName, out);
			}
		}
		
		IOUtils.closeBufferedWriter(out);
	}
	
	// ========================================================================
	private static void printSegmentedRegionsToFile_Helper(CopyNumberRegionRange cnrr, StringBuilder sb, String delim, double copyNumBase, String sampleName, BufferedWriter out) {
		sb.setLength(0);
		
		double errorFactor = Math.random() / 10;  // Get between 0 and 0.1
		boolean addPositiveError = (Math.random() >= 0.5);
		double copyNum = copyNumBase + ((addPositiveError) ? errorFactor : -errorFactor);
		double copyNumRatio = copyNum / DefaultDiploidCopyNumber;
		double log2CopyNum = Math.log10(copyNum) / Math.log10(2);
		double log2Ratio = Math.log10(copyNumRatio) / Math.log10(2);

		
		sb.append(sampleName)
		  .append(delim).append(cnrr.getChromosome().getCode())
		  .append(delim).append(cnrr.getRangeStart())
		  .append(delim).append(cnrr.getRangeEnd())
		  .append(delim).append(cnrr.getNumSitesInterrogated())
		  .append(delim).append(log2Ratio);
		
		IOUtils.writeToBufferedWriter(out, sb.toString(), true);
	}
	
	// ========================================================================
	/** Given the segmented regions in one sample, this method merges the regions of one type (LOH, Amp, etc) 
	 *  in the following manner.  LOH/Amp regions cannot be longer than maxLengthContiguousRegion long. 
	 */
	public static CopyNumberRegionsByChromosome mergeRegionsWithConstraintsOld(CopyNumberRegionsByChromosome regionsInSample, EventType clusterType, int maxBasePairsContiguousRegion) {
		
		// Create an empty return object
		CopyNumberRegionsByChromosome regionsInSampleMerged = new CopyNumberRegionsByChromosome(regionsInSample.mSampleName);
		CopyNumberRegionRange regionTest = new CopyNumberRegionRange(EventType.Ignored, Chrom.c0, 0);
		
		// Go chromosome by chromosome
		for (Chrom chrom : Chrom.Autosomes) {				
			ArrayList<CopyNumberRegionRange> regionsInChromOriginal =       regionsInSample.mRegionsByChrom.get(chrom);
			ArrayList<CopyNumberRegionRange> regionsInChromMerged   = regionsInSampleMerged.mRegionsByChrom.get(chrom);
						
			// We declare a stored region that can be extended.  Initialize to null for now
			CopyNumberRegionRange regionToExtend = null;
			
			// Iterate through the regions for this chromosome
			for (CopyNumberRegionRange currentRegion : regionsInChromOriginal) {

				if (currentRegion.mCopyNumberEventType == clusterType) {
					// Check if there's a region already waiting for extension.  
					// If not, create a new one (and a copy at that), and add to array
					if (regionToExtend == null) {
						regionToExtend = currentRegion.getCopy();
						regionsInChromMerged.add(regionToExtend);  // add this to the new array
						
					} else {
						if (currentRegion.mCopyNumberEventType != regionToExtend.mCopyNumberEventType) {
							CompareUtils.throwErrorAndExit("ERROR: Must have same cluster type!\t" + currentRegion.mCopyNumberEventType + "\t" + regionToExtend.mCopyNumberEventType);
						}
						
						int maxEndIndexInclusive = regionToExtend.getRangeEnd() + maxBasePairsContiguousRegion - 1;						
						if (currentRegion.getRangeStart() <= maxEndIndexInclusive) {
							regionToExtend.setRangeEnd(currentRegion.getRangeEnd());	
							regionToExtend.incrementSitesInterrogated(currentRegion.getNumSitesInterrogated());
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
	/** Given the segmented regions in one sample, this method merges the regions of one type (LOH, Amp, etc) 
	 *  in the following manner.  LOH/Amp regions cannot be longer than maxLengthContiguousRegion long. 
	 */
	public static CopyNumberRegionsByChromosome mergeRegionsWithConstraints(CopyNumberRegionsByChromosome regionsInSample, EventType eventType, int maxBasePairsContiguousRegion, ClusteringInputOneSample oneSampleData, ClusteringInputOneSampleMetaData metaData) {
		
		// Create an empty return object
		CopyNumberRegionsByChromosome regionsInSampleMerged = new CopyNumberRegionsByChromosome(regionsInSample.mSampleName);
		CopyNumberRegionRange regionTest = new CopyNumberRegionRange(EventType.Ignored, Chrom.c0, 0);
		
		// Go chromosome by chromosome
		for (Chrom chrom : Chrom.Autosomes) {				
			ArrayList<CopyNumberRegionRange> regionsInChromOriginal =       regionsInSample.mRegionsByChrom.get(chrom);
			ArrayList<CopyNumberRegionRange> regionsInChromMerged   = regionsInSampleMerged.mRegionsByChrom.get(chrom);
						
			// We declare a stored region that can be extended.  Initialize to null for now
			CopyNumberRegionRange regionToExtend = null;
			
			// Iterate through the regions for this chromosome
			for (CopyNumberRegionRange currentRegion : regionsInChromOriginal) {

				if (currentRegion.mCopyNumberEventType == eventType) {
					// Check if there's a region already waiting for extension.  
					// If not, create a new one (and a copy at that), and add to array
					if (regionToExtend == null) {
						regionToExtend = currentRegion.getCopy();
						regionsInChromMerged.add(regionToExtend);  // add this to the new array
						
					} else {
						if (currentRegion.mCopyNumberEventType != regionToExtend.mCopyNumberEventType) {
							CompareUtils.throwErrorAndExit("ERROR: Must have same cluster type!\t" + currentRegion.mCopyNumberEventType + "\t" + regionToExtend.mCopyNumberEventType);
						}
						
						boolean shouldCombine = Clustering.combineTwoRegions(regionToExtend, currentRegion, oneSampleData, metaData);					
						
						//int maxEndIndexInclusive = regionToExtend.getRangeEnd() + maxBasePairsContiguousRegion - 1;						
						if (shouldCombine /* currentRegion.getRangeStart() <= maxEndIndexInclusive */) {
							regionToExtend.setRangeEnd(currentRegion.getRangeEnd());	
							regionToExtend.incrementSitesInterrogated(currentRegion.getNumSitesInterrogated());
						} else {
							// The current region is out of bounds.  We simply set
							// the current region as the new region to extend.
							regionToExtend = currentRegion.getCopy();
							regionsInChromMerged.add(regionToExtend);  // add this to the new array
						}
					}	
				} else {
					if (eventType.isCopyNumberChangeGermlineOrSomatic() && currentRegion.mCopyNumberEventType.isCopyNumberChangeGermlineOrSomatic()) {
						// We know that we have a copy number change event that does not match the targeted 
						// event type (by virtue of the if-else chain) , and so we do not allow extension 
						// from a previous region of the targeted event type.
						regionToExtend = null;
					}
				}
			}		
		}
		
		return regionsInSampleMerged;
	}
	
	
	// ========================================================================
	/** The continuation of the @method segmentRegionsAllFiles method, but by individual. */ 
	public static CopyNumberRegionsByChromosome segmentRegionsOneFile(File inFile, String outDir, SNVMap snvMap) {
		StringUtils.FileExtensionAndDelimiter fileExtAndDelim = StringUtils.FileTextTabDelim;		
		
		// First check that the file is a file of the desired extension
		String suffix = ".withCopyNum" + fileExtAndDelim.mExtension;
		if (!inFile.getName().endsWith(suffix)) return null;
		String samplenameRoot = inFile.getName().substring(0, inFile.getName().length() - suffix.length());
		System.out.println("Processing Sample: " + samplenameRoot);
		
		StringBuilder sb = new StringBuilder(4096);
		
		// Load all lines into memory, extract the header row, and then sort by chrom/position
		ArrayList<String> allLines = IOUtils.readAllLinesFromFile(inFile.getAbsolutePath(), false, true, sb);
		String headerString = sb.toString();
		Collections.sort(allLines, LineComparatorTab);	
		
		// Read the event results
		ArrayList<EventType> events = new ArrayList<EventType>(allLines.size());
		for (String line : allLines) {
			EventType event = EventType.getClusterType(StringUtils.extractNthColumnValue(line, Regions.ColCuratedTSV_Cluster, fileExtAndDelim.mDelimiter));
			CompareUtils.ensureTrue(CompareUtils.isNotNull(event), "ERROR: segmentRegionsOneFile(): Event is not a valid enumerated event type!");
			events.add(event);
		}
		
		ClusteringInputOneSample oneSampleInfo = new ClusteringInputOneSample(allLines, samplenameRoot);		
		CopyNumberRegionsByChromosome rV = segmentRegionsOneSample(oneSampleInfo, events, outDir, snvMap);
		allLines.clear();
		oneSampleInfo.clear();
		return rV;
	}
	
	// ========================================================================
	public static CopyNumberRegionsByChromosome segmentRegionsOneSample(ClusteringInputOneSample oneSampleInfo, ArrayList<EventType> eventPerSite, String outDir, SNVMap snvMap) {
		// Have an array of regions for amplifications and LOH
		CopyNumberRegionsByChromosome regionsByChrom = new CopyNumberRegionsByChromosome(oneSampleInfo.getSampleNameRoot());		 	
		CopyNumberRegionRange currentRegion = null;
		ChromPositionTracker chromPosTrack = new ChromPositionTracker();		
		
		// We start at index 0 assuming no header and that the rows are sorted by chrom/position
		int numSites = oneSampleInfo.getNumSites();
		for (int row = 0; row < numSites; row++) {			
			ClusteringInputOneSite oneSiteInfo = oneSampleInfo.getSiteAtIndex(row);
			
			final Chrom chrom  = oneSiteInfo.getChrom();									
			final int position = oneSiteInfo.getPosition();
			EventType eventType = eventPerSite.get(row);
			final int rsId = (oneSiteInfo.getRsID() < 0) ? 0 : oneSiteInfo.getRsID();
			
			// Convert cnLOH cluster type to LOH
			//eventType = eventType.isLOH() ? ClusterType.LOH : eventType;
			
			// Register the site in the map
			if (CompareUtils.isNotNull(snvMap)) {
				snvMap.registerSNV(chrom, position, rsId, Nuc.N, Nuc.N, true, true);
			}
			
			// If the chromosome has changed, we set that we have no current region
			if (chromPosTrack.chromCrossedWithCurrentCoordinates(chrom, position) && (currentRegion != null)) {
				currentRegion.makeFinalized();
				currentRegion = null;
			}
			
			// Now determine whether we create or extend regions
			if (currentRegion == null) {
				if (segmentRegionsOneFile_isValidEvent(eventType)) {
					currentRegion = new CopyNumberRegionRange(eventType, chrom, position);
					regionsByChrom.addRegion(chrom, currentRegion);
				}
			} else {
				// Ensure we're on the same chrom
				CompareUtils.ensureTrue(currentRegion.getChromosome() == chrom, "ERROR: Chromosomes different!");
				
				// Compare the position, make sure it's after the current range end
				if (position < currentRegion.getRangeEnd()) {
					CompareUtils.throwErrorAndExit("ERROR: segmentRegionsOneFile(): Rows not sorted!");					
				} else if (position == currentRegion.getRangeEnd()) {
					System.out.println("WARNING: Duplicate coordinates: " + chrom + "\t" + position);
					continue;  // We ignore duplicate positions
				}
								
				// Now we know the position is after the current range end				
				if ((currentRegion.mCopyNumberEventType == eventType) 
					//&& (position - currentRegion.getRangeEnd() < REGION_SEGMENTATION_DIST_THRESHOLD * 1000) 
				 	) {
					boolean result = currentRegion.extendRange(chrom, position);
					if (!result) {
						CompareUtils.throwErrorAndExit("ERROR: Could not extend range! " + currentRegion.toString() + "\t" + chrom + "\t" + position + "\t" + oneSampleInfo.getSampleNameRoot());
					}
				} else {
					// The cluster types are different, or the chrom didn't match, 
					// or the position was too far away.  Create a new region and add it
					// if it is not null or it is not noise.					
					if (segmentRegionsOneFile_isValidEvent(eventType)) {
						currentRegion.makeFinalized();
						currentRegion = new CopyNumberRegionRange(eventType, chrom, position);
						regionsByChrom.addRegion(chrom, currentRegion);
					}
				}
			}
		}
		
		if (currentRegion != null) currentRegion.makeFinalized();
		
		System.out.println(oneSampleInfo.getSampleNameRoot());				
		return regionsByChrom;

	}
	
	// ========================================================================
	private static boolean segmentRegionsOneFile_isValidEvent(EventType eventType) {
		return (eventType != EventType.Noise && 
				eventType != EventType.Ignored  && 
				eventType != EventType.HETSomatic);
	}

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
		takeUnionAndBreakDownIntersectingRegions(CopyNumberRegionsByChromosome regionsTarget, CopyNumberRegionsByChromosome regionsSource, final EventType clusterType) {
		
		// First iterate over the chromosomes
		for (Chrom chrom : Chrom.Autosomes) {
			ArrayList<CopyNumberRegionRange> regionsChrTarget = regionsTarget.mRegionsByChrom.get(chrom);
			ArrayList<CopyNumberRegionRange> regionsChrSource = regionsSource.mRegionsByChrom.get(chrom);
			takeUnionAndBreakDownIntersectingRegions(regionsChrTarget, regionsChrSource, clusterType);	
		}

		return regionsTarget;
	}
	
	public static ArrayList<CopyNumberRegionRange>
		takeUnionAndBreakDownIntersectingRegions(ArrayList<CopyNumberRegionRange> regionsTarget, ArrayList<CopyNumberRegionRange> regionsSource, final EventType clusterType) {		
		
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
			if (regionTarget.mCopyNumberEventType != clusterType) {
				regionsTarget.remove(indexTarget);
				//indexTarget++;
				continue;
			} else if (regionSource.mCopyNumberEventType != clusterType) {
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
				CompareUtils.throwErrorAndExit("ERROR: Predicted and determined region overlap types don't match!" 
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
				CompareUtils.throwErrorAndExit("ERROR: Should be on different chromosomes!");					
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
				CompareUtils.throwErrorAndExit("ERROR: Invalid option!");
			}				
		}
				

		// Add any remaining regions.  We know that either of the following two loops will execute, but not both
		// since the previous loop was broken by failure of one of the loop conditions.
		// 
		// We do not need to add any regions to the target array, as either the array elements
		// were already traversed, or they already exist in the original array.  However, we need to remove
		// elements that may not match the cluster type desired
		while (indexTarget < regionsTarget.size()) {
			if (regionsTarget.get(indexTarget).mCopyNumberEventType != clusterType) {
				regionsTarget.remove(indexTarget);
			} else {
				indexTarget++;
			}
		}
		
		// We only need to add elements (actually, their copies) if more still exist in the source array. 
		for (; indexSource < regionsSource.size(); indexSource++) {
			CopyNumberRegionRange regionSource = regionsSource.get(indexSource);
			if (regionSource.mCopyNumberEventType == clusterType) {
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
	private static void getBrowserTracksHelper_Helper_findMinMaxPositionsAcrossSamplesPerChromosome(int[] minPerChrom, int[] maxPerChrom, ArrayList<File> sampleFiles, EventType clusterType) {
		Arrays.fill(minPerChrom, Integer.MAX_VALUE);
		Arrays.fill(maxPerChrom, Integer.MIN_VALUE);
		
		for (File inFile : sampleFiles) {
			int lineCounter = -1;
			String line;
			BufferedReader in = IOUtils.getBufferedReader(inFile.getAbsolutePath());
			
			while ((line = IOUtils.getNextLineInBufferedReader(in)) != null) {
				if (++lineCounter == 0) { continue; }
				
				final Chrom chrom  = Chrom.getChrom  (StringUtils.extractNthColumnValue(line, ColCuratedTSV_Chrom,    StringUtils.TabStr));					
				final int position = Integer.parseInt(StringUtils.extractNthColumnValue(line, ColCuratedTSV_Position, StringUtils.TabStr));
				final EventType clusterTypeForSite = EventType.getClusterType(StringUtils.extractNthColumnValue(line, ColCuratedTSV_Cluster,      StringUtils.TabStr));

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
										EventType clusterType, String outDir) {		

		// First, we want to find the minimum and maximum positions of sites across samples for each chromosome
		int numChromArrayElements = Chrom.values().length;		
		int[] minStartPositionOnChromosomeAcrossSamples = ArrayUtils.newIntArray(numChromArrayElements, Integer.MAX_VALUE); 				
		int[] maxStartPositionOnChromosomeAcrossSamples = ArrayUtils.newIntArray(numChromArrayElements, Integer.MIN_VALUE);
		getBrowserTracksHelper_Helper_findMinMaxPositionsAcrossSamplesPerChromosome
			(minStartPositionOnChromosomeAcrossSamples, maxStartPositionOnChromosomeAcrossSamples, sampleFiles, clusterType);

		// Initialize color array and event count 2D array
		int[] rgb = ArrayUtils.newIntArray(3, 0);  				
		int[][] eventCount = ArrayUtils.newIntArray2D(numChromArrayElements, sampleFiles.size(), 0); 
		
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
				CompareUtils.throwErrorAndExit("ERROR: Samples don't match up: " + regionsOneSample.mSampleName + "\t" + inFile.getName());
			}
			String sampleNameRoot = regionsOneSample.mSampleName;
			
			System.out.println("Processing sample (" + sampleIndex + " / " + sampleFiles.size() + "):\t" + sampleNameRoot);

			// This is the list of regions
			int[] rowNum = ArrayUtils.newIntArray(numChromArrayElements, 0);			
			
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
	
				for (CopyNumberRegionRange cnrr : regionsOneSample.mRegionsByChrom.get(chrom)) {
					constructRowString_BEDFormat(chrom, sb, cnrr.getRangeStart(), cnrr.getRangeEnd() + 1, ++rowNum[chrom.ordinal()], snowString);
					outListForChrom.add(sb.toString());
				}
			}
						
			// Go through all the lines
			for (int lineIndex = 0; lineIndex < allLines.size(); lineIndex++) {
				String line = allLines.get(lineIndex);

				final Chrom chrom  = Chrom.getChrom                       (StringUtils.extractNthColumnValue(line, ColCuratedTSV_Chrom,    StringUtils.TabStr));					
				final int position = Integer.parseInt                     (StringUtils.extractNthColumnValue(line, ColCuratedTSV_Position, StringUtils.TabStr));
				final EventType clusterTypeForSite = EventType.getClusterType(StringUtils.extractNthColumnValue(line, ColCuratedTSV_Cluster,      StringUtils.TabStr));
				final MutationType mutationType = MutationType.getSNVType                  (StringUtils.extractNthColumnValue(line, ColCuratedTSV_MutationType, StringUtils.TabStr));
				final double vafNormal = Double.parseDouble(StringUtils.extractNthColumnValue(line, ColCuratedTSV_VafNormal, StringUtils.TabStr));
				final VariantLocation varLoc = Clustering.isVariantInGermline(vafNormal) ? VariantLocation.Germline : VariantLocation.Somatic;
				
				boolean addSite = false;
				ColorPastel theColorMutationType = ColorPastel.Pastel_Cyan;  // set as default
				
				if (clusterTypeForSite == clusterType) {					
					// Increment the event count for this chromosome and sample
					eventCount[chrom.ordinal()][sampleIndex]++;					
										
					if (mutationType != null) {
						switch(mutationType) {
						case NonSynonymous_SNV: theColorMutationType = ColorPastel.RGB_Red;    break;
						case Synonymous_SNV:    theColorMutationType = ColorPastel.Dark_Green; break;		
						default: break;
						}
					}
					
					if (varLoc == VariantLocation.Somatic) {
						theColorMutationType = ColorPastel.Blue;
					}
					
				} else if  (clusterTypeForSite == EventType.HETGermline) {					
					//score = 300; //used with grayscale BEDs (historical artifact)
					
					if (mutationType != null) {
						switch(mutationType) {
						case NonSynonymous_SNV: theColorMutationType = ColorPastel.Yellow_Orange; break;
						case Synonymous_SNV:    theColorMutationType = ColorPastel.Light_Green;   break;	
						default: break;
						}
					}
				} else if (clusterTypeForSite == EventType.HETSomatic) {
					addSite = true;
					theColorMutationType = ColorPastel.Darker_Blue_Violet;
				}

				//if (score == greyscaleBEDScore) {
				if (addSite) {
					theColorMutationType.getRGBValues(rgb);					
					sb.setLength(0);
					sb.append(rgb[0]).append(StringUtils.CommaStr).append(rgb[1]).append(StringUtils.CommaStr).append(rgb[2]);
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
			for (CopyNumberRegionRange cnrr : recurrentRegionsForOneClusterType.mRegionsByChrom.get(chrom)) {
				recurrenceMin = Math.min(recurrenceMin, cnrr.mRecurrenceScore);
				recurrenceMax = Math.max(recurrenceMax, cnrr.mRecurrenceScore);
			}
			float recurrenceMinMaxRange = recurrenceMax - recurrenceMin + 1;
			
			// Scale accordingly -- this idea was inspired from Siddharth Reddy
			for (CopyNumberRegionRange cnrr : recurrentRegionsForOneClusterType.mRegionsByChrom.get(chrom)) {
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
		sb.append("browser position chr").append(chrom.ordinal()).append(StringUtils.ColonString)
		  .append(rangeStart).append("-").append(rangeEnd);
		return sb;
	}
	
	// ========================================================================
	public static StringBuilder constructTrackNameString_BED(String sampleNameRoot, StringBuilder sb, boolean useRGB) {
		sb.setLength(0);				
		sb.append("track name=").append(StringUtils.DoubleQuoteStr).append(sampleNameRoot).append(StringUtils.DoubleQuoteStr).append(StringUtils.SpaceString)
		  .append("description=").append(StringUtils.DoubleQuoteStr).append(StringUtils.SpaceString).append(StringUtils.DoubleQuoteStr).append(StringUtils.SpaceString)
		  .append("visibility=dense").append(StringUtils.SpaceString);
		
		if (useRGB) {
			sb.append("itemRgb=").append(StringUtils.DoubleQuoteStr).append("on").append(StringUtils.DoubleQuoteStr);
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
		   append(StringUtils.SpaceString).append(rangeStart).
		   append(StringUtils.SpaceString).append(rangeEnd).
		   append(StringUtils.SpaceString).append("row").append(rowNum).
		   append(StringUtils.SpaceString).append(score).
		   append(StringUtils.SpaceString).append("+").				   
		   append(StringUtils.SpaceString).append(rangeStart).
		   append(StringUtils.SpaceString).append(rangeEnd);
		
		if (!colorString.isEmpty()) {
			sb.append(StringUtils.SpaceString).append(colorString);
		}
		
		return sb;
	}
	
	
	// ========================================================================
	// Convenience function to write the manifest file
	private static void writeManifestFile(String outDir, EventType clusterType, int[][] eventCount, 
										  ArrayList<CopyNumberRegionsByChromosome> samplesWithRegions,
										  String recurrenceNameRoot) {
		
		//String filePrefix = outDir + File.separator;
		String filePrefix = "http://ron.cs.columbia.edu/ninad/LOHcate/crc/browser_tracks/LOH/toUpload/";
		
		for (Chrom chrom : Chrom.values()) {
			if (chrom.isInvalid()) continue;

			String outFilename = 
					outDir + File.separator + clusterType.name() + StringUtils.DotStr + ChromPrefix + chrom.ordinal() + StringUtils.DotStr + "Manifest.txt";
			
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
	private static String constructGenBrowserTrackFilenamePerSampleAndCluster(String sampleName, EventType clusterType, Chrom chrom) {
		return clusterType.name() + StringUtils.DotStr + ChromPrefix + chrom.ordinal() + StringUtils.DotStr 
				                  + sampleName + GenBrowserTrack + StringUtils.FileExtensionTSV.mExtension;
	}
	
	// ========================================================================

	
	/**
	 * Generate 'master' gene enrichment table (used to generate histograms).
	 * @param inDir curated SNP calls
	 */
	public static void getGeneEnrichment(String inDir, String outDir) {
		String outFilename = outDir + File.separator + "geneEnrichment" + StringUtils.FileExtensionTSV.mExtension;
		IOUtils.createDirectoryPath(outDir, false);
		File[] files = (new File(inDir)).listFiles();
		StringUtils.FileExtensionAndDelimiter fileExtDelim = StringUtils.FileTextTabDelim; 
				
		StringBuilder sb = new StringBuilder(4096);
		BucketCounterEnum<EventType> clusterTypeCountsGene = new BucketCounterEnum<EventType>(EventType.class);
		
		// A list of genes that we will keep binary sorted for efficieny purposes
		ArrayList<Gene> genes = new ArrayList<Gene>(); 
		Gene dummyGene = new Gene("", Chrom.c0);  // We'll use this gene for binary searching
		
		for (File file : files) { //iterate through curated SNP calls
			
			int indexOfDelimiter  = file.getName().indexOf(fileExtDelim.mExtension); 			
			if (indexOfDelimiter == file.getName().length() - fileExtDelim.mExtension.length()) {
				System.out.println(file.getName());
				
				// Load all lines into memory, extract the header row, and then sort by chrom/position
				ArrayList<String> allLines = IOUtils.readAllLinesFromFile(file.getAbsolutePath(), false, true, sb);
				String headerString = sb.toString();
				Collections.sort(allLines, LineComparatorTab);					
				
				int startingRow = 0;  // because the header was stripped away
				for (int row = startingRow; row < allLines.size(); row++) {
				
					//if (row % 10000 == 0) System.out.println(row + " of " + allLines.size());
					
					String components[] = allLines.get(row).split(StringUtils.TabStr); 
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
							resultIndex = ArrayUtils.getInsertPoint(resultIndex);  // calculate the proper insertion point 							 
							genes.add(resultIndex, new Gene(geneName, chrom));							
						}
						Gene currentGene = genes.get(resultIndex);
						clusterTypeCountsGene.clear();
						
						if (position > currentGene.mMaxBasePairPosition) { //get right-bound of gene's range of variants
							currentGene.mMaxBasePairPosition = position; 
						}
						if (position < currentGene.mMinBasePairPosition) { //...left-bound...
							currentGene.mMinBasePairPosition = position;
						}
												
						MutationType mutationType = MutationType.getSNVType(components[ColCuratedTSV_MutationType]);
						currentGene.incrementCount(mutationType);  // increment synonymous or nonsynonymous variant count						
						
						VariantLocation variantLocation = 
								Clustering.isVariantInGermline(Double.parseDouble(components[ColCuratedTSV_VafNormal])) ? VariantLocation.Germline : VariantLocation.Somatic; 						
						currentGene.incrementCount(variantLocation);
																
						EventType clusterType = EventType.getClusterType(components[ColCuratedTSV_Cluster]);
						if (clusterType == null) { CompareUtils.throwErrorAndExit("ERROR: Invalid cluster type: " + components[ColCuratedTSV_Cluster]); }
						currentGene.incrementCount(clusterType);  // increment LOH/DUP/&c. count
						double vafTumor = Double.parseDouble(components[ColCuratedTSV_VafTumor]);
						if ((clusterType == EventType.LOH) && (vafTumor > 0.5)) {
							currentGene.mCountLOHreferenceLost++;
							currentGene.addPatientIfNotAlreadyAdded_LOHRefLost(file.getName(), position);
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
				MutationType.NonSynonymous_SNV.toLowerCase(), 
				MutationType.Synonymous_SNV.toLowerCase(), 
				VariantLocation.Germline.toLowerCase(), 
				VariantLocation.Somatic.toLowerCase(),
				
				EventType.GainSomatic.name(),
				EventType.LOH.name(), 
				EventType.LOH.name() + "_refLost",
				EventType.LOH.name() + "_refLost_Positions", 
				EventType.HETGermline.name(),
				EventType.HETSomatic.name(),
				
				EventType.GainSomatic.name() + logStr,
				EventType.LOH.name() + logStr, 
				EventType.LOH.name() + "_refLost" + logStr, 
				EventType.HETGermline.name() + logStr,
				EventType.HETSomatic.name() + logStr,
				
				EventType.GainSomatic.name()              + densityStr,
				EventType.LOH.name()              + densityStr, 
				EventType.LOH.name() + "_refLost" + densityStr, 
				EventType.HETGermline.name()      + densityStr,
				EventType.HETSomatic.name()       + densityStr,
				
				EventType.GainSomatic.name()         + recurrenceStr, 
				EventType.LOH.name()         + recurrenceStr, 
				EventType.HETGermline.name() + recurrenceStr,
				EventType.HETSomatic.name()  + recurrenceStr
		};
		String headerStr = StringUtils.constructColumnDelimitedString(columnHeaders, fileExtDelim.mDelimiter, sb, true).toString();

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
		
		String eventsByPatientPerGeneFilename = outDir + File.separator + "eventsByPatientPerGene" + StringUtils.FileExtensionTSV.mExtension;
		BufferedWriter outBreakdown = IOUtils.getBufferedWriter(eventsByPatientPerGeneFilename);
				
		// Now write out individual gene outputs
		for (Gene gene : genes) {
			// Write samples for each gene out as well
			ArrayList<String> patientsAllClusters = new ArrayList<String>();
			for (EventType ct : EventType.values()) {
				if ((ct == EventType.Noise) || (ct == EventType.Ignored)) continue;
				
				ArrayList<String> patientsForEvent = gene.getPatientsForClusterType(ct);
				for (String patientForEvent : patientsForEvent) {
					String outString = gene.getName() + "\t" + ct.name() + "\t" + patientForEvent;
					IOUtils.writeToBufferedWriter(outBreakdown, outString, true);
				}
			}
		}
		IOUtils.closeBufferedWriter(outBreakdown);
	}
	
	
	// ========================================================================
	// ENTRY POINT
	// ========================================================================
	public static void main(String[] args) {
				
	
	}

}
