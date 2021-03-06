package lohcate;
import genomeEnums.Chrom;
import genomeEnums.Nuc;
import genomeEnums.VariantLocation;
import genomeUtils.ChromPositionTracker;
import genomeUtils.GenotypeUtils;
import genomeUtils.RegionAndSiteWalker;
import genomeUtils.RegionBreakerAndIntersecter;
import genomeUtils.RegionBreakerAndIntersecter.RegionIntersectTester;
import genomeUtils.SNVMap;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.ListIterator;

import nutils.ArrayUtils;
import nutils.Cast;
import nutils.CompareUtils;
import nutils.ControlFlagBool;
import nutils.EnumMapSafe;
import nutils.IOUtils;
import nutils.NullaryClassFactory;
import nutils.NumberUtils;
import nutils.StringUtils;
import nutils.BitUtils.Compactor.CompactorIntoLong;
import nutils.collectionsSorted.ArrayListSortedComparable;
import nutils.counter.DynamicBucketCounter;
import nutils.primitives.wrapper.PrimitiveWrapper;

import lohcate.clustering.AlleleFractionStatsForSample;
import lohcate.clustering.Clustering;
import lohcate.clustering.ClusteringInputOneSample;
import lohcate.clustering.ClusteringInputOneSampleMetaData;
import lohcate.clustering.ClusteringInputOneSite;
import lohcate.clustering.ClusteringParams;
import lohcate.clustering.ClusteringPlotting;
import lohcateEnums.EventType;
import lohcateEnums.ColorPastel;
import lohcateEnums.MutationType;

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
	public static final int ColCuratedTSV_Chrom    = Col_NAFTAFInput_Chrom;
	public static final int ColCuratedTSV_Position = Col_NAFTAFInput_Position;
	//private static final int ColCuratedTSV_VariantBaseTumor = 2;
	private static final int ColCuratedTSV_VafNormal    = Col_NAFTAFInput_VariantRatioNormal;
	public static final int ColCuratedTSV_VafTumor     = Col_NAFTAFInput_VariantRatioTumor;
	public static final int ColCuratedTSV_dbSNP        = Col_NAFTAFInput_DbSNPString;	
	public static final int ColCuratedTSV_MutationType = Col_NAFTAFInput_MutationType;
	public static final int ColCuratedTSV_Gene         = Col_NAFTAFInput_HugoSymbol;	        
	//private static final int ColCuratedTSV_VariantLocation = 7;
	public static final int ColCuratedTSV_Cluster      = 17; //8;
	public static final int ColCuratedTSV_CopyNumber   = 19;
	
	public static final GenotypeUtils.GenomicCoordinateComparatorInTextFileLine LineComparatorTab = new GenotypeUtils.GenomicCoordinateComparatorInTextFileLine();
	public static final NullaryClassFactory<ArrayList> NullClassFactoryArrayList = new NullaryClassFactory<ArrayList>(ArrayList.class);
	
	// ========================================================================
	/** This removes the header lines and filters out sites. */
	public static ArrayList<String> curateSNPCalls_removeHeaderLinesFromRows(ArrayList<String> rows, boolean removeHighDensitySNVSites, boolean removeExtremeGCSites) {
		int windowPositionStart = 0;
		int windowRowStart = 0;		
		int numSitesRemoved = 0;
		int numSitesRemovedGC = 0;
		int numVariantSitesSpanned = 0;
		ChromPositionTracker chromTracker = new ChromPositionTracker();
		
		for (int rowIndex = 0; rowIndex < rows.size(); rowIndex++) {
			String row = rows.get(rowIndex);
			if (row.indexOf("refName") >= 0 && row.indexOf("coord") >= 0 || /*row.indexOf("chrX") >= 0 ||*/ row.indexOf("chrY") >= 0 || row.indexOf("chrM") >= 0 /*|| row.indexOf("intergenic") >= 0 */) {
				rows.set(rowIndex, null);
			} else {
				// Do some light GC filtering
				Chrom chrom = Chrom.getChrom(    StringUtils.extractNthColumnValue(row, Regions.Col_NAFTAFInput_Chrom,               StringUtils.FileExtensionTSV.mDelimiter));
				int position = Integer.parseInt( StringUtils.extractNthColumnValue(row, Regions.Col_NAFTAFInput_Position,            StringUtils.FileExtensionTSV.mDelimiter));
				String gcString =                StringUtils.extractNthColumnValue(row, Regions.Col_NAFTAFInput_FlankingStringTumor, StringUtils.FileExtensionTSV.mDelimiter);				
				
				int posDiff = position - windowPositionStart;
				
				if (chromTracker.chromCrossedWithCurrentCoordinates(chrom, position) || (posDiff >= MaxWindowLength)) {					
					if (numVariantSitesSpanned > MaxSitesInWindowAllowed) {
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
					numVariantSitesSpanned = 0;
				}
				
				// Now check the vaf normal and see whether it's an outlier.
				double vafNormal = 0;
				if (ClusteringParams.GlobalClusteringParams.isTumorOnly()) {
					short readsVariantTumor = Short.parseShort( StringUtils.extractNthColumnValue(row, Regions.Col_NAFTAFInput_VariantCoverageTumor, StringUtils.FileExtensionTSV.mDelimiter) );
					short readsTotalTumor   = Short.parseShort( StringUtils.extractNthColumnValue(row, Regions.Col_NAFTAFInput_TotalCoverageTumor,   StringUtils.FileExtensionTSV.mDelimiter) );
					vafNormal = (double) readsVariantTumor / (double) readsTotalTumor;					
				} else {
					vafNormal = Double.parseDouble(StringUtils.extractNthColumnValue(row, Regions.Col_NAFTAFInput_VariantRatioNormal, StringUtils.FileExtensionTSV.mDelimiter));					
				}
				
				if (vafNormal > AlleleFractionStatsForSample.VAFNormalRange.getBoundLower()) {						
					numVariantSitesSpanned++;
				}
				
				double fractionGCNormal = GenotypeUtils.calcFractionGC(gcString);
				if (fractionGCNormal >= 0) {
					// We know we have a valid GC fraction
					if (removeExtremeGCSites && (fractionGCNormal < GCContentThresholdLow) || (fractionGCNormal >= GCContentThresholdHigh)) {
						rows.set(rowIndex, null);
						numSitesRemovedGC++;
					}
				}
			}
		}

		System.out.println("\tNum Sites Removed in Windows: " + numSitesRemoved);
		System.out.println("\tNum Sites Removed in GC: " + numSitesRemovedGC);
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
		snvMap.printMe(true, outDir + File.separator + "AllSites.txt", true);
		
		// Now we have all the contiguous regions from all the samples.  Find the regions of the cluster type		
		// Declare the maximum stretch of a region for a particular cluster type
		EnumMapSafe<EventType, Integer> maxBasePairsContiguousRegion = EnumMapSafe.createEnumMap(EventType.class, Integer.MAX_VALUE);
		maxBasePairsContiguousRegion.put(EventType.GainSomatic, REGION_SEGMENTATION_DIST_THRESHOLD);
		maxBasePairsContiguousRegion.put(EventType.cnLOH,       REGION_SEGMENTATION_DIST_THRESHOLD);
		maxBasePairsContiguousRegion.put(EventType.LOH,         REGION_SEGMENTATION_DIST_THRESHOLD);
		
		EnumMapSafe<EventType, ArrayList<CopyNumberRegionsByChromosome>> regionsInSamplesPerEventType =
				EnumMapSafe.createEnumMapOfArrayLists(EventType.class, CopyNumberRegionsByChromosome.class);
		
		for (EventType clusterType : EventType.AmpLOHHetG) { 						
			ArrayList<CopyNumberRegionsByChromosome> regionsInSamplesForOneClusterType = regionsInSamplesPerEventType.get(clusterType);	
			/*
			for (CopyNumberRegionsByChromosome regionsInOneSample : regionsInSamples) {
				CopyNumberRegionsByChromosome regionsInOneSampleMerged = mergeRegionsWithConstraintsOld(regionsInOneSample, events, oneSampleData, metaData) 
						mergeRegionsWithConstraintsOld(regionsInOneSample, clusterType, maxBasePairsContiguousRegion.get(clusterType));
				regionsInOneSampleMerged.removeSingletonRegions();
				regionsInSamplesForOneClusterType.add(regionsInOneSampleMerged);
				printSegmentedRegionsToFile(outDir, regionsInOneSampleMerged, clusterType, snvMap);
			}
			*/
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
			ArrayList<CopyNumberRegionRangeLOHcate> cnrrOnChrom = recurrentRegionsForOneClusterType.getRegions(chrom);
			if (cnrrOnChrom.isEmpty()) continue;
			
			for (CopyNumberRegionRangeLOHcate cnrr : cnrrOnChrom) {
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
		regionsGenomeWide.clearEventCounts();
		
		StringBuilder sb = new StringBuilder(4096);

		int sampleIndex = -1;
		for (File inFile : sampleFiles) {
			
			// Load all lines into memory, extract the header row, and then sort by chrom/position
			ArrayList<String> allLines = IOUtils.readAllLinesFromFile(inFile.getAbsolutePath(), false, true, sb);
			String headerString = sb.toString();
			Collections.sort(allLines, LineComparatorTab);
						
			CopyNumberRegionsByChromosome regionsOneSample = regionsGenomeWideSampleSpecific.get(++sampleIndex);
			if (inFile.getName().indexOf(regionsOneSample.getSampleName()) < 0) {
				CompareUtils.throwErrorAndExit("ERROR: Samples don't match up: " + regionsOneSample.getSampleName() + "\t" + inFile.getName());
			}
			
			// We initialize some indices for efficiency purposes
			Chrom prevChrom = null;
			ArrayList<CopyNumberRegionRangeLOHcate> regionsInChr = null;
			ArrayList<CopyNumberRegionRangeLOHcate> regionsInChrSampleSpecific = null;
			PrimitiveWrapper.WInteger regionIndexInChr               = new PrimitiveWrapper.WInteger(-1);
			PrimitiveWrapper.WInteger regionIndexInChrSampleSpecific = new PrimitiveWrapper.WInteger(-1);
			
			// Go through all the lines
			for (int lineIndex = 0; lineIndex < allLines.size(); lineIndex++) {
				String line = allLines.get(lineIndex);
				
				final Chrom chrom  = Chrom.getChrom                       (StringUtils.extractNthColumnValue(line, ColCuratedTSV_Chrom,    StringUtils.TabStr));					
				final int position = Integer.parseInt                     (StringUtils.extractNthColumnValue(line, ColCuratedTSV_Position, StringUtils.TabStr));
				final EventType clusterType = EventType.getClusterType(StringUtils.extractNthColumnValue(line, ColCuratedTSV_Cluster,  StringUtils.TabStr));
				
				if ((prevChrom == null) || (chrom != prevChrom)) {
					regionsInChr =                regionsGenomeWide.getRegions(chrom);
					regionsInChrSampleSpecific =   regionsOneSample.getRegions(chrom);					
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
						CopyNumberRegionRangeLOHcate region = regionsInChr.get(regionIndexInChr.mInt);
						region.mEventTypeCounts.increment(clusterType);
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
	public static Boolean scanRegionsForPoint(ArrayList<? extends CopyNumberRegionRangeLOHcate> regionsInChromosome, 
											  PrimitiveWrapper.WInteger regionIndex, 
											  Chrom chrom, int position) {
		
		for (; regionIndex.mInt < regionsInChromosome.size(); regionIndex.mInt++) {			
			CopyNumberRegionRangeLOHcate region = regionsInChromosome.get(regionIndex.mInt);
			
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
		CopyNumberRegionsByChromosome target = regionsInSamples.get(0).makeClone();
		
		// Now traverse the rest of the samples
		for (int i = 1; i < regionsInSamples.size(); i++) {
			takeUnionAndBreakDownIntersectingRegions(target, regionsInSamples.get(i), clusterType);
		}
		
		return target;
	}

	// ========================================================================
	/** Given the merged segmented regions in one sample, this method prints the regions to file. */
	public static void printSegmentedRegionsToFile(String outDir, CopyNumberRegionsByChromosome regionsInSample, EventType clusterType) {
				
		String outFilename = outDir + File.separator + "Regions.GISTIC." + clusterType.name() + "." + regionsInSample.getSampleName() + ".txt";
		StringBuilder sb = new StringBuilder(2048);
		String delim = StringUtils.FileExtensionTSV.mDelimiter;
		
		
		double copyNumBase = (clusterType == EventType.GainSomatic) ? 2.5 : (clusterType == EventType.LOH ? 1.5 : 2.0);		
//		double log2RatioGisticCol = log2Ratio - 1.0;
//		double log2CopyNumGistic = log2CopyNum - 1.0;
		
		BufferedWriter out = IOUtils.getBufferedWriter(outFilename);
		// Go chromosome by chromosome
		for (Chrom chrom : Chrom.Autosomes) {							
			
			ArrayList<CopyNumberRegionRangeLOHcate> regionsInChrom = regionsInSample.getRegions(chrom);			
			for (CopyNumberRegionRangeLOHcate cnrr : regionsInChrom) {
				sb.setLength(0);
				
				double errorFactor = Math.random() / 10;  // Get between 0 and 0.1
				boolean addPositiveError = (Math.random() >= 0.5);
				double copyNum = copyNumBase + ((addPositiveError) ? errorFactor : -errorFactor);
				double copyNumRatio = copyNum / DefaultDiploidCopyNumber;
				double log2CopyNum = Math.log10(copyNum) / Math.log10(2);
				double log2Ratio = Math.log10(copyNumRatio) / Math.log10(2);

				
				sb.append(regionsInSample.getSampleName())
				  .append(delim).append(chrom.getName())
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

	public static class ActionerGISTIC implements RegionAndSiteWalker.Actioner<CopyNumberRegionRangeLOHcate> {

		// ====================================================================		
		private static final String Delim = StringUtils.FileExtensionTSV.mDelimiter;
		StringBuilder mSB;
		ClusteringInputOneSample mOneSampleInfo; 
		ClusteringInputOneSampleMetaData mMetaData;
		BufferedWriter mOut;
		
		// Variables used in methods
		ArrayListSortedComparable<String> mHugoSymbols;
		
		public ActionerGISTIC(
				BufferedWriter outputWriter, 
				CopyNumberRegionsByChromosome regionsInSample, 
				ArrayList<EventType> events, 
				ClusteringInputOneSample oneSampleInfo, 
				ClusteringInputOneSampleMetaData metaData, 
				SNVMap snvMap) {
						
			mSB = new StringBuilder(2048);
			RegionAndSiteWalker<CopyNumberRegionRangeLOHcate> regionSiteWalker = new RegionAndSiteWalker<CopyNumberRegionRangeLOHcate>(CopyNumberRegionRangeLOHcate.ClassFactory);
			mOneSampleInfo = oneSampleInfo;
			mMetaData = metaData;
			mOut = outputWriter;
			mHugoSymbols = new ArrayListSortedComparable<String>();
			
			
			for (Chrom chrom : Chrom.Autosomes) {
				ArrayList<CopyNumberRegionRangeLOHcate> reigonsInChrom = new ArrayList<CopyNumberRegionRangeLOHcate>(regionsInSample.getRegions(chrom));				
				regionSiteWalker.walk(reigonsInChrom, snvMap, chrom, this);
			}
		}
		
		// ====================================================================
		@Override
		public void takeAction(ArrayList<Boolean> rangeInTargetSet, ArrayList<CopyNumberRegionRangeLOHcate> regions, CopyNumberRegionRangeLOHcate regionLatest) {
			
			int indexOfStart = mOneSampleInfo.getIndex(regionLatest.getChromosome(), regionLatest.getRangeStart());
			int indexOfEnd   = mOneSampleInfo.getIndex(regionLatest.getChromosome(), regionLatest.getRangeEnd());
			
			boolean regionLatestInTargetSet = rangeInTargetSet.get(rangeInTargetSet.size() - 1).booleanValue();
			
			// Test if the latest range is in the target set of ranges.  In that case, we use
			// the copy number in the target range, otherwise we assume the diploid default
			double log2CopyNumRatio = NumberUtils.MathLog2(Regions.DefaultTumorNormalRatio);
			mHugoSymbols.clear();
						
			if (regionLatestInTargetSet) {
				CompareUtils.ensureTrue(indexOfStart >= 0, "ERROR: Site must exist in sample!");
				double copyNumRatio = mMetaData.getTumorNormalCopyNumRatioAtIndex(indexOfStart);
				log2CopyNumRatio = NumberUtils.MathLog2(copyNumRatio);
				//copyNumber = mMetaData.getCopyNumberAtIndex(indexOfSite);
								
				for (int row = indexOfStart; row <= indexOfEnd; row++) { 
					ClusteringInputOneSite oneSiteInfo = mOneSampleInfo.getSiteAtIndex(regionLatest.getChromosome(), row);
					String hugoSymbol = oneSiteInfo.getHugoSymbol();					
					if (mHugoSymbols.get(hugoSymbol) == null) {
						mHugoSymbols.add(hugoSymbol);
					}
				}
			}
			
			// Reset the string builder
			mSB.setLength(0);
			mSB.append(mOneSampleInfo.getSampleNameRoot())
			  .append(Delim).append(regionLatest.getChromosome().getCode())
			  .append(Delim).append(regionLatest.getRangeStart())
			  .append(Delim).append(regionLatest.getRangeEnd())
			  .append(Delim).append(regionLatest.getNumSitesInterrogated())
			  .append(Delim).append(log2CopyNumRatio)
			  .append(Delim).append(regionLatest.mCopyNumberEventType);
			
			// Append the hugo symbols string
			mSB.append(Delim);			
			for (int i = 0; i < mHugoSymbols.size(); i++) {
				if (i > 0) mSB.append(";");
				mSB.append(mHugoSymbols.get(i));
			}
			
			IOUtils.writeToBufferedWriter(mOut, mSB.toString(), true);			
		}
		
	}
	// ========================================================================
	// ========================================================================
	
	
	
	public static void printSegmentedRegionsToFile(String outDir, CopyNumberRegionsByChromosome regionsInSample, EventType clusterType, SNVMap snvMap) {
		
		String outFilename = outDir + File.separator + "Regions.GISTIC." + clusterType.name() + "." + regionsInSample.getSampleName() + ".txt";
		StringBuilder sb = new StringBuilder(2048);
		String delim = StringUtils.FileExtensionTSV.mDelimiter;
		CopyNumberRegionRangeLOHcate dummyRange = new CopyNumberRegionRangeLOHcate(EventType.HETGermline, Chrom.c0, 0);		
		double dummyRangeCopyNumberBase = 2.0;
		
		double copyNumBase = (clusterType == EventType.GainSomatic) ? 2.5 : (clusterType == EventType.LOH ? 1.5 : 2.0);		
//		double log2RatioGisticCol = log2Ratio - 1.0;
//		double log2CopyNumGistic = log2CopyNum - 1.0;
		
		BufferedWriter out = IOUtils.getBufferedWriter(outFilename);
		// Go chromosome by chromosome
		for (Chrom chrom : Chrom.Autosomes) {				
			
			int numSitesOnChrom = snvMap.getNumSitesOnChromosome(chrom);
			if (numSitesOnChrom <= 0) continue;
			
			ListIterator<CopyNumberRegionRangeLOHcate> regionsInChromIter = regionsInSample.getRegions(chrom).listIterator();
			CopyNumberRegionRangeLOHcate currentRegion = null;
			int indexInMap = 0;
			boolean dummyRangeValid = false;
						
			// First scan and see if we have any regions preceding the first chromosomal position in the map
			while (regionsInChromIter.hasNext()) {
				currentRegion = regionsInChromIter.next();
				if (currentRegion.afterRange(chrom, snvMap.getPosition(chrom, indexInMap))) {
					printSegmentedRegionsToFile_Helper(currentRegion, sb, delim, copyNumBase, regionsInSample.getSampleName(), out);
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
							printSegmentedRegionsToFile_Helper(dummyRange, sb, delim, dummyRangeCopyNumberBase, regionsInSample.getSampleName(), out);
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
						printSegmentedRegionsToFile_Helper(dummyRange, sb, delim, dummyRangeCopyNumberBase, regionsInSample.getSampleName(), out);
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
					printSegmentedRegionsToFile_Helper(currentRegion, sb, delim, copyNumBase, regionsInSample.getSampleName(), out);				
					
					// Now move to the next cna-affected region
					currentRegion = regionsInChromIter.hasNext() ? regionsInChromIter.next() : null;
					
				} else if (currentRegion.afterRange(chrom, mapPosition)) {
					CompareUtils.ensureTrue(false, "ERROR: printSegmentedRegionsToFile(): Impossible for region to precede position!");
				} else {
					CompareUtils.ensureTrue(false, "ERROR: printSegmentedRegionsToFile(): Impossible state!");
				}
			}
			
			if (dummyRangeValid) {
				printSegmentedRegionsToFile_Helper(dummyRange, sb, delim, dummyRangeCopyNumberBase, regionsInSample.getSampleName(), out);
			}
		}
		
		IOUtils.closeBufferedWriter(out);
	}
	
	// ========================================================================
	private static void printSegmentedRegionsToFile_Helper(CopyNumberRegionRangeLOHcate cnrr, StringBuilder sb, String delim, double copyNumBase, String sampleName, BufferedWriter out) {
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
	public static CopyNumberRegionsByChromosome reassignRegionsByCopyNumber(CopyNumberRegionsByChromosome regionsInSample, ArrayList<EventType> events, ClusteringInputOneSample oneSampleData, ClusteringInputOneSampleMetaData metaData) {
		
		// Create an empty return object
		CopyNumberRegionsByChromosome regionsInSampleMerged = new CopyNumberRegionsByChromosome(regionsInSample.getSampleName());
		CopyNumberRegionRangeLOHcate regionTest = new CopyNumberRegionRangeLOHcate(EventType.Ignored, Chrom.c0, 0);
		
		PrimitiveWrapper.WDouble probFromMaxLikelihood = new PrimitiveWrapper.WDouble(0);
		// Go chromosome by chromosome
		for (Chrom chrom : Chrom.Autosomes) {				
			ArrayList<CopyNumberRegionRangeLOHcate> regionsInChromOriginal =       regionsInSample.getRegions(chrom);
			ArrayList<CopyNumberRegionRangeLOHcate> regionsInChromMerged   = regionsInSampleMerged.getRegions(chrom);
			
			// We declare a stored region that can be extended.  Initialize to null for now
			CopyNumberRegionRangeLOHcate regionToExtend = null;
			
			// Iterate through the regions for this chromosome
			for (CopyNumberRegionRangeLOHcate currentRegion : regionsInChromOriginal) {
				
				if (currentRegion.mCopyNumberEventType == EventType.HETGermline) {
					double avgCopyNum = Clustering.calcAverageCopyNumberOverRegion(currentRegion, oneSampleData, metaData);	

					if (Clustering.isCopyNumDepletedHemizygous(avgCopyNum)) {
						currentRegion.mCopyNumberEventType = EventType.LOH;
					} else if (Clustering.isCopyNumDepletedHomozygous(avgCopyNum)) {
						currentRegion.mCopyNumberEventType = EventType.DELHom;
					} else if (Clustering.isCopyNumAmplified(avgCopyNum)) {
						currentRegion.mCopyNumberEventType = EventType.GainSomatic;
					} else if (Clustering.isCopyNumAmplified(avgCopyNum)) {
						currentRegion.mCopyNumberEventType = EventType.cnLOH;						
					}
					
					long compactUnit = oneSampleData.getIndicesForRegion(currentRegion);
					int indexStart = Cast.toInt(CompactorIntoLong.TwoIntsIntoLong.Compactor.getValue(CompactorIntoLong.TwoIntsIntoLong.IntMSB, compactUnit));
					int indexEnd   = Cast.toInt(CompactorIntoLong.TwoIntsIntoLong.Compactor.getValue(CompactorIntoLong.TwoIntsIntoLong.IntLSB, compactUnit));
					
					// First fill the region
					Clustering.fillRegion(events, indexStart, indexEnd, false, currentRegion.mCopyNumberEventType, oneSampleData, metaData);
					
					// Now, if we had a cnLOH event, we test
					if (currentRegion.mCopyNumberEventType == EventType.cnLOH) {
						Clustering.fillRegionBasedOnVAFMaxLikelihood(currentRegion, oneSampleData, metaData, events, currentRegion.mCopyNumberEventType, probFromMaxLikelihood, true);
					}
					
				}
			}
		}
		
		return regionsInSample;
	}
	
	// ========================================================================
	private static boolean mergeRegionsWithConstraints_Helper_ShouldCombine(EventType et, double avgCopyNum, boolean shouldCombineInitialCondition) {
		if (LOHcate.LOHcateSensitivity.isHighOrMore()) {
			switch(et) {
			case LOH:         shouldCombineInitialCondition  = Clustering.isCopyNumDepletedHemizygous(avgCopyNum); break;
			case DELHom:      shouldCombineInitialCondition  = Clustering.isCopyNumDepletedHomozygous(avgCopyNum); break;
			case cnLOH:       shouldCombineInitialCondition &= Clustering.isCopyNumInDiploidRange(avgCopyNum); break;
			case GainSomatic: shouldCombineInitialCondition  = Clustering.isCopyNumAmplified(avgCopyNum); break;
			default: break;
			}
		} else if (LOHcate.LOHcateSensitivity.isLow()) {
			switch(et) {
			case LOH:         shouldCombineInitialCondition &= Clustering.isCopyNumDepletedHemizygous(avgCopyNum); break;
			case DELHom:      shouldCombineInitialCondition &= Clustering.isCopyNumDepletedHomozygous(avgCopyNum); break;
			case cnLOH:       shouldCombineInitialCondition &= Clustering.isCopyNumInDiploidRange(avgCopyNum); break;
			case GainSomatic: shouldCombineInitialCondition &= Clustering.isCopyNumAmplified(avgCopyNum); break;
			default: break;
			}
		}
		return shouldCombineInitialCondition;
	}
	
	// ========================================================================
	/** Given the segmented regions in one sample, this method merges the regions of one type (LOH, Amp, etc) 
	 *  in the following manner.  LOH/Amp regions cannot be longer than maxLengthContiguousRegion long. 
	 */
	public static CopyNumberRegionsByChromosome mergeRegionsWithConstraints(CopyNumberRegionsByChromosome regionsInSample, EventType eventType, int maxBasePairsContiguousRegion, ClusteringInputOneSample oneSampleData, ClusteringInputOneSampleMetaData metaData) {
		
		// Create an empty return object
		CopyNumberRegionsByChromosome regionsInSampleMerged = new CopyNumberRegionsByChromosome(regionsInSample.getSampleName());
		CopyNumberRegionRangeLOHcate regionTest = new CopyNumberRegionRangeLOHcate(EventType.Ignored, Chrom.c0, 0);
		
		// Go chromosome by chromosome
		for (Chrom chrom : Chrom.Autosomes) {				
			ArrayList<CopyNumberRegionRangeLOHcate> regionsInChromOriginal =       regionsInSample.getRegions(chrom);
			ArrayList<CopyNumberRegionRangeLOHcate> regionsInChromMerged   = regionsInSampleMerged.getRegions(chrom);
						
			// We declare a stored region that can be extended.  Initialize to null for now
			CopyNumberRegionRangeLOHcate regionToExtend = null;
			
			// Iterate through the regions for this chromosome
			for (CopyNumberRegionRangeLOHcate currentRegion : regionsInChromOriginal) {

				if (currentRegion.mCopyNumberEventType == eventType) {
					// Check if there's a region already waiting for extension.  
					// If not, create a new one (and a copy at that), and add to array
					if (regionToExtend == null) {
						regionToExtend = currentRegion.makeClone();
						regionsInChromMerged.add(regionToExtend);  // add this to the new array
						
					} else {
						if (currentRegion.mCopyNumberEventType != regionToExtend.mCopyNumberEventType) {
							CompareUtils.throwErrorAndExit("ERROR: Must have same cluster type!\t" + currentRegion.mCopyNumberEventType + "\t" + regionToExtend.mCopyNumberEventType);
						}
						
						boolean shouldCombine = Clustering.combineTwoRegions(regionToExtend, currentRegion, oneSampleData, metaData, true);					
						System.out.printf("RegionsCombined: %b\t%s\t%d\t%d\t%d\t%d\t%d\n", shouldCombine, eventType.name(), regionToExtend.getChromosome().getCode(), regionToExtend.getRangeStart(), regionToExtend.getRangeEnd(), currentRegion.getRangeStart(), currentRegion.getRangeEnd());
						
						CopyNumberRegionRangeLOHcate midRegion = Clustering.getMiddleRegion(regionToExtend, currentRegion, oneSampleData, false);
						if (midRegion != null) {
							double avgCopyNum = Clustering.calcAverageCopyNumberOverRegion(midRegion, oneSampleData, metaData);	
							shouldCombine = mergeRegionsWithConstraints_Helper_ShouldCombine(regionToExtend.mCopyNumberEventType, avgCopyNum, shouldCombine);
						}
						
						if (shouldCombine) {
							regionToExtend.setRangeEnd(currentRegion.getRangeEnd());	
							regionToExtend.incrementSitesInterrogated(currentRegion.getNumSitesInterrogated());
						} else {
							// The current region is out of bounds.  We simply set
							// the current region as the new region to extend.
							regionToExtend = currentRegion.makeClone();
							regionsInChromMerged.add(regionToExtend);  // add this to the new array
						}						
					}	
				} else {
					if (eventType.isCopyNumberChangeGermlineOrSomatic() && currentRegion.mCopyNumberEventType.isCopyNumberChangeGermlineOrSomatic()) {
						// We know that we have a copy number change event that does not match the targeted 
						// event type (by virtue of the if-else chain) , and so we do not allow extension 
						// from a previous region of the targeted event type.
						regionToExtend = null;
						
					} else if (regionToExtend != null) {
						
						// Ensure we have the same targeted event type
						CompareUtils.ensureTrue(regionToExtend.mCopyNumberEventType == eventType, "ERROR: Impossible state!");
					
						// Only do something if we're in a non-cna region
						if (currentRegion.mCopyNumberEventType == EventType.HETGermline) {

							boolean shouldCombine = Clustering.combineTwoRegions(regionToExtend, currentRegion, oneSampleData, metaData, false);
							System.out.printf("RegionsCombinedHET: %b\t%s\t%d\t%d\t%d\t%d\t%d\n", shouldCombine, eventType.name(), regionToExtend.getChromosome().getCode(), regionToExtend.getRangeStart(), regionToExtend.getRangeEnd(), currentRegion.getRangeStart(), currentRegion.getRangeEnd());						

							double avgCopyNum = Clustering.calcAverageCopyNumberOverRegion(currentRegion, oneSampleData, metaData);							
							shouldCombine = mergeRegionsWithConstraints_Helper_ShouldCombine(regionToExtend.mCopyNumberEventType, avgCopyNum, shouldCombine);
							
							ControlFlagBool performExtendedSmoothing = new ControlFlagBool(true);							
							if (performExtendedSmoothing.getValue()) {
								if (shouldCombine) {
									regionToExtend.setRangeEnd(currentRegion.getRangeEnd());	
									regionToExtend.incrementSitesInterrogated(currentRegion.getNumSitesInterrogated());								
								} else {
									//if (regionToExtend.mCopyNumberEventType == EventType.cnLOH) {
									//	regionToExtend = null;
									//}
								}
							}
						}
						
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
		CopyNumberRegionRangeLOHcate currentRegion = null;
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
					currentRegion = new CopyNumberRegionRangeLOHcate(eventType, chrom, position);
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
						currentRegion = new CopyNumberRegionRangeLOHcate(eventType, chrom, position);
						regionsByChrom.addRegion(chrom, currentRegion);
					}
				}
			}
		}
		
		if (currentRegion != null) currentRegion.makeFinalized();
		
		System.out.println("Finished initial segmentation of sites into regions for:\t" + oneSampleInfo.getSampleNameRoot());				
		return regionsByChrom;

	}
	
	// ========================================================================
	private static boolean segmentRegionsOneFile_isValidEvent(EventType eventType) {
		return (eventType != EventType.Noise && 
				eventType != EventType.Ignored  && 
				eventType != EventType.HETSomatic);
	}

	// ========================================================================
	public static class RegionTester implements RegionIntersectTester<CopyNumberRegionRangeLOHcate> { 

		// Member variables
		final EventType mEventType;
				
		public RegionTester(EventType eventType) { mEventType = eventType; }
		
		@Override
		public void takeActionOnEqualRegions(CopyNumberRegionRangeLOHcate region) { region.mRecurrenceScore += 1.0; }
		
		@Override
		public boolean isValidRegion(CopyNumberRegionRangeLOHcate region) { return (region.mCopyNumberEventType == mEventType); }
		
		@Override
		public boolean isInvalidRegion(CopyNumberRegionRangeLOHcate region) { return !isValidRegion(region); }
	};
	
	public static EnumMapSafe<EventType, RegionTester> RegionTesterByEvent = new EnumMapSafe<EventType, RegionTester>(EventType.class);

	// ========================================================================
//	public static<T extends CopyNumberRegionsByChromosome, E extends CopyNumberRegionRange> 
//		void takeUnionAndBreakDownIntersectingRegions(T regionsTarget, T regionsSource, RegionTester regionTester, Class<E> theClass) {
//		
//		for (Chrom chrom : Chrom.Autosomes) {
//			ArrayList<E> regionsChrTarget = regionsTarget.getRegions(chrom);					
//			ArrayList<E> regionsChrSource = regionsSource.getRegions(chrom);
//			RegionBreakerAndIntersecter.takeUnionAndBreakDownIntersectingRegions(regionsChrTarget, regionsChrSource, regionTester, CopyNumberRegionRange.class);	
//		}
//		
//		
//	}
	
	// ========================================================================
	public static CopyNumberRegionsByChromosome
		takeUnionAndBreakDownIntersectingRegions(CopyNumberRegionsByChromosome regionsTarget, CopyNumberRegionsByChromosome regionsSource, final EventType eventType) {
		
		
		RegionTester regionTester = RegionTesterByEvent.get(eventType);
		if (CompareUtils.isNull(regionTester)) {
			regionTester = new RegionTester(eventType);
			RegionTesterByEvent.put(eventType, regionTester);
		}		
		
		// First iterate over the chromosomes
		for (Chrom chrom : Chrom.Autosomes) {
			ArrayList<CopyNumberRegionRangeLOHcate> regionsChrTarget = regionsTarget.getRegions(chrom);
			ArrayList<CopyNumberRegionRangeLOHcate> regionsChrSource = regionsSource.getRegions(chrom);
			RegionBreakerAndIntersecter.takeUnionAndBreakDownIntersectingRegions(regionsChrTarget, regionsChrSource, regionTester);	
		}

		return regionsTarget;
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
			if (inFile.getName().indexOf(regionsOneSample.getSampleName()) < 0) {
				CompareUtils.throwErrorAndExit("ERROR: Samples don't match up: " + regionsOneSample.getSampleName() + "\t" + inFile.getName());
			}
			String sampleNameRoot = regionsOneSample.getSampleName();
			
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
	
				for (CopyNumberRegionRangeLOHcate cnrr : regionsOneSample.getRegions(chrom)) {
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
			for (CopyNumberRegionRangeLOHcate cnrr : recurrentRegionsForOneClusterType.getRegions(chrom)) {
				recurrenceMin = Math.min(recurrenceMin, cnrr.mRecurrenceScore);
				recurrenceMax = Math.max(recurrenceMax, cnrr.mRecurrenceScore);
			}
			float recurrenceMinMaxRange = recurrenceMax - recurrenceMin + 1;
			
			// Scale accordingly -- this idea was inspired from Siddharth Reddy
			for (CopyNumberRegionRangeLOHcate cnrr : recurrentRegionsForOneClusterType.getRegions(chrom)) {
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
				String theSampleName = samplesWithRegions.get(indexSample).getSampleName();
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

	
	// ========================================================================
	// ENTRY POINT
	// ========================================================================
	public static void main(String[] args) {
				
	
	}

}
