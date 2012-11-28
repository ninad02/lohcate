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
import java.util.Comparator;
import java.util.ListIterator;

import org.jfree.chart.renderer.xy.XYDotRenderer;

import com.martiansoftware.jsap.FlaggedOption;
import com.martiansoftware.jsap.JSAP;
import com.martiansoftware.jsap.JSAPResult;

import lohcateEnums.Chrom;
import lohcateEnums.ClusterType;
import lohcateEnums.ColorPastel;
import lohcateEnums.MutationType;
import lohcateEnums.Nuc;
import lohcateEnums.SeqPlatform;
import lohcateEnums.VariantLocation;
import shared.ArgumentParserUtils;
import shared.BucketCounter;
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
	
	public static final int DefaultDiploidCopyNumber = 2;
	
	private static final int REGION_SEGMENTATION_DIST_THRESHOLD = 2000000; //greatest possible distance between 2 'adjacent' points of LOH in a region of 'contiguous' LOH
	
	public static final String GermlineSuffix = "." + VariantLocation.Germline.toLowerCase();	
	public static final String NovelStr  = "novel";	
	public static final String ChromPrefix = "chr";
	public static final String MissingGeneNameValue = ".";
	public static final float MaxVariantAlleleFrequency = 1.0f;
	public static final String GenBrowserTrack = ".genBrowserTrack";

	// Column constants for the curated TSV files (files that have a cluster column)
	public static final int Col_NAFTAFInput_Chrom = 0;
	public static final int Col_NAFTAFInput_Position = 1;
	public static final int Col_NAFTAFInput_FlankingStringNormal = 6;
	public static final int Col_NAFTAFInput_FlankingStringTumor  = 7;
	public static final int Col_NAFTAFInput_TotalCoverageNormal  = 8;
	public static final int Col_NAFTAFInput_TotalCoverageTumor   = 9;
	public static final int Col_NAFTAFInput_VariantCoverageNormal = 10;
	public static final int Col_NAFTAFInput_VariantCoverageTumor  = 11;
	public static final int Col_NAFTAFInput_VariantRatioNormal    = 12;
	public static final int Col_NAFTAFInput_VariantRatioTumor     = 13;
	public static final int Col_NAFTAFInput_DbSNPString           = 14;
	public static final int Col_NAFTAFInput_MutationType          = 15;
	public static final int Col_NAFTAFInput_HugoSymbol            = 16;
	
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
	
	private static final double GCContentThresholdLow  = 0.05;
	private static final double GCContentThresholdHigh = 0.80;
	
	private static final int MaxWindowLength = 1000;
	private static final int MaxSitesInWindowAllowed = 4;
	
	// Column constants for the curated TSV files (files that have a cluster column)
	private static final int ColCuratedTSV_Chrom    = 0;
	private static final int ColCuratedTSV_Position = 1;
	//private static final int ColCuratedTSV_VariantBaseTumor = 2;
	private static final int ColCuratedTSV_VafTumor = 3;
	private static final int ColCuratedTSV_dbSNP    = 4;
	private static final int ColCuratedTSV_Gene     = 5;
	private static final int ColCuratedTSV_MutationType = 6;
	private static final int ColCuratedTSV_VariantLocation = 7;
	private static final int ColCuratedTSV_Cluster = 8;
	
	public static final GenomicCoordinateComparatorInTextFileLine LineComparatorTab = new GenomicCoordinateComparatorInTextFileLine();
	
	static ArrayList<String> curateSNPCalls_removeHeaderLinesFromRows(ArrayList<String> rows) {
		int windowPositionStart = 0;
		int windowRowStart = 0;
		Chrom chromPrev = Chrom.c0;
		int numSitesRemoved = 0;
		
		for (int rowIndex = 0; rowIndex < rows.size(); rowIndex++) {
			String row = rows.get(rowIndex);
			if (row.indexOf("refName") >= 0 && row.indexOf("coord") >= 0 || row.indexOf("chrX") >= 0) {
				rows.set(rowIndex, null);
			} else {
				// Do some light GC filtering
				Chrom chrom = Chrom.getChrom(    Utils.extractNthColumnValue(row, Script.Col_NAFTAFInput_Chrom,               Utils.FileExtensionTSV.mDelimiter));
				int position = Integer.parseInt( Utils.extractNthColumnValue(row, Script.Col_NAFTAFInput_Position,            Utils.FileExtensionTSV.mDelimiter));
				String gcString =                Utils.extractNthColumnValue(row,        Col_NAFTAFInput_FlankingStringTumor, Utils.FileExtensionTSV.mDelimiter);
				
				int posDiff = position - windowPositionStart;
				
				if ((chrom != chromPrev) || (posDiff >= MaxWindowLength)) {
					int numSitesSpanned = (rowIndex - 1) - windowRowStart + 1;
					if (numSitesSpanned > MaxSitesInWindowAllowed) {
						for (int rowToClean = windowRowStart; rowToClean < rowIndex; rowToClean++) {
							//System.out.println("Removing:\t" + rows.get(rowToClean));
							rows.set(rowToClean, null);
							numSitesRemoved++;
						}
					}	

					windowPositionStart = position;
					windowRowStart = rowIndex;
					chromPrev = chrom;
				} 
				
				double fractionGCNormal = Utils.calcFractionGC(gcString);
				if ((fractionGCNormal < GCContentThresholdLow) || (fractionGCNormal >= GCContentThresholdHigh)) {
					rows.set(rowIndex, null);
				}
			}
		}
		System.out.println("\tNum Sites Removed in Windows: " + numSitesRemoved);
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
		Utils.FileExtensionAndDelimiter fileExtDelim = Utils.FileExtensionTSV;
		File[] files = (new File(inDir)).listFiles();
		ArrayList<CopyNumberRegionsByChromosome> regionsInSamples = new ArrayList<CopyNumberRegionsByChromosome>();
		ArrayList<File> validFilesList = new ArrayList<File>(files.length);
		SNVMap snvMap = new SNVMap();
		
		// Create our output directories
		IOUtils.createDirectoryPath(outDir, false);
		IOUtils.createDirectoryPath(outDirBrowserTracks, false);
		
		// First we determine the regions from each sample
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
		int[] maxBasePairsContiguousRegion = new int[] {
				REGION_SEGMENTATION_DIST_THRESHOLD,
				REGION_SEGMENTATION_DIST_THRESHOLD,
				Integer.MAX_VALUE
		};  
		
		ArrayList<ArrayList<CopyNumberRegionsByChromosome>> regionsInSamplesPerClusterType = new ArrayList<ArrayList<CopyNumberRegionsByChromosome>>();
		for (int clusterIndex = 0; clusterIndex < ClusterType.AmpLOHHetG.length; clusterIndex++) {
			ClusterType clusterType = ClusterType.AmpLOHHetG[clusterIndex];
			ArrayList<CopyNumberRegionsByChromosome> regionsInSamplesForOneClusterType = new ArrayList<CopyNumberRegionsByChromosome>();
			regionsInSamplesPerClusterType.add(regionsInSamplesForOneClusterType);			
			int maxBasePairsContiguousRegionForCluster = maxBasePairsContiguousRegion[clusterIndex];
			
			for (CopyNumberRegionsByChromosome regionsInOneSample : regionsInSamples) {
				CopyNumberRegionsByChromosome regionsInOneSampleMerged = mergeRegionsWithConstraints(regionsInOneSample, clusterType, maxBasePairsContiguousRegionForCluster);
				regionsInSamplesForOneClusterType.add(regionsInOneSampleMerged);
				printSegmentedRegionsToFile(outDir, regionsInOneSampleMerged, clusterType, snvMap);
			}
		}

		PrintStream out = IOUtils.getPrintStream(outDir + File.separator + "testOut.txt");
		
		// Now we want to determine the recurrent regions
		ArrayList<CopyNumberRegionsByChromosome> regionsRecurrentPerClusterType = new ArrayList<CopyNumberRegionsByChromosome>();
		for (int clusterIndex = 0; clusterIndex < ClusterType.AmpLOHHetG.length; clusterIndex++) {
			ClusterType clusterType = ClusterType.AmpLOHHetG[clusterIndex];
			
			CopyNumberRegionsByChromosome recurrentRegionsForOneClusterType = 
					determineRecurrentRegions(regionsInSamplesPerClusterType.get(clusterIndex), clusterType);
			regionsRecurrentPerClusterType.add(recurrentRegionsForOneClusterType);
			recurrentRegionsForOneClusterType.print(out, fileExtDelim.mDelimiter);
		}
		
		
		System.out.println("Scoring regions...");
		// Now, we need to go through the recurrent regions and score them, based on the 
		// contents of the individual samples that make up the recurrent regions.
		for (int clusterIndex = 0; clusterIndex < ClusterType.AmpLOHHetG.length; clusterIndex++) {
			ClusterType clusterType = ClusterType.AmpLOHHetG[clusterIndex];

			CopyNumberRegionsByChromosome recurrentRegionsForOneClusterType = regionsRecurrentPerClusterType.get(clusterIndex);
			ArrayList<CopyNumberRegionsByChromosome> regionsInSamplesForOneClusterType = regionsInSamplesPerClusterType.get(clusterIndex);
			
			countClusterTypesInRegions(recurrentRegionsForOneClusterType, regionsInSamplesForOneClusterType, validFilesList);
			recurrentRegionsForOneClusterType.print(out, fileExtDelim.mDelimiter);
		}
				
		IOUtils.closePrintStream(out);

		System.out.println("Generating Browser tracks...");
		// Generatebrowser tracks
		//for (int clusterIndex = 0; clusterIndex < ClusterType.OnlyLOH.length; clusterIndex++) {
		for (ClusterType clusterType : ClusterType.OnlyLOH) {
			//ClusterType clusterType = ClusterType.OnlyLOH[clusterIndex];
			ArrayList<CopyNumberRegionsByChromosome> regionsInSamplesForOneClusterType = regionsInSamplesPerClusterType.get(clusterType.ordinal());
			CopyNumberRegionsByChromosome recurrentRegionsForOneClusterType = regionsRecurrentPerClusterType.get(clusterType.ordinal());
			
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
	/** Given the merged segmented regions in one sample, this method prints the regions to file. */
	public static void printSegmentedRegionsToFile(String outDir, CopyNumberRegionsByChromosome regionsInSample, ClusterType clusterType) {
				
		String outFilename = outDir + File.separator + "Regions.GISTIC." + clusterType.name() + "." + regionsInSample.mSampleName + ".txt";
		StringBuilder sb = new StringBuilder(2048);
		String delim = Utils.FileExtensionTSV.mDelimiter;
		
		
		double copyNumBase = (clusterType == ClusterType.Amp) ? 2.5 : (clusterType == ClusterType.LOH ? 1.5 : 2.0);		
//		double log2RatioGisticCol = log2Ratio - 1.0;
//		double log2CopyNumGistic = log2CopyNum - 1.0;
		
		BufferedWriter out = IOUtils.getBufferedWriter(outFilename);
		// Go chromosome by chromosome
		for (Chrom chrom : Chrom.Autosomes) {							
			
			ArrayList<CopyNumberRegionRange> regionsInChrom = regionsInSample.mRegionsByChrom.get(chrom.getCode());			
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
	
	
	public static void printSegmentedRegionsToFile(String outDir, CopyNumberRegionsByChromosome regionsInSample, ClusterType clusterType, SNVMap snvMap) {
		
		String outFilename = outDir + File.separator + "Regions.GISTIC." + clusterType.name() + "." + regionsInSample.mSampleName + ".txt";
		StringBuilder sb = new StringBuilder(2048);
		String delim = Utils.FileExtensionTSV.mDelimiter;
		CopyNumberRegionRange dummyRange = new CopyNumberRegionRange(ClusterType.HETGermline, Chrom.c0, 0);		
		double dummyRangeCopyNumberBase = 2.0;
		
		double copyNumBase = (clusterType == ClusterType.Amp) ? 2.5 : (clusterType == ClusterType.LOH ? 1.5 : 2.0);		
//		double log2RatioGisticCol = log2Ratio - 1.0;
//		double log2CopyNumGistic = log2CopyNum - 1.0;
		
		BufferedWriter out = IOUtils.getBufferedWriter(outFilename);
		// Go chromosome by chromosome
		for (Chrom chrom : Chrom.Autosomes) {				
			
			int numSitesOnChrom = snvMap.getNumSitesOnChromosome(chrom);
			if (numSitesOnChrom <= 0) continue;
			
			ListIterator<CopyNumberRegionRange> regionsInChromIter = regionsInSample.mRegionsByChrom.get(chrom.getCode()).listIterator();
			CopyNumberRegionRange currentRegion = null;
			int indexInMap = 0;
			boolean dummyRangeValid = false;
						
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
								Utils.ensureTrue(extendResult, errorString);	
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
						Utils.ensureTrue(false, errorString);
					}
										
					// Now change the loop index to move to the map position just after the end of this region
					int indexOfRegionEndInMap = snvMap.getIndexOfPositionInMap(chrom, currentRegion.getRangeEnd());
					if (indexOfRegionEndInMap < 0) {
						String errorString = "ERROR: printSegmentedRegionsToFile(): Region end (" + currentRegion.getRangeEnd() + "must exist in map!";
						Utils.ensureTrue(false, errorString);						
					}
					indexInMap = indexOfRegionEndInMap;  // indexInMap will be incremented at loop end
					
					currentRegion.set(chrom, currentRegion.getRangeStart(), currentRegion.getRangeEnd(), true, indexOfRegionEndInMap - indexOfRegionStartInMap + 1);
					printSegmentedRegionsToFile_Helper(currentRegion, sb, delim, copyNumBase, regionsInSample.mSampleName, out);				
					
					// Now move to the next cna-affected region
					currentRegion = regionsInChromIter.hasNext() ? regionsInChromIter.next() : null;
					
				} else if (currentRegion.afterRange(chrom, mapPosition)) {
					Utils.ensureTrue(false, "ERROR: printSegmentedRegionsToFile(): Impossible for region to precede position!");
				} else {
					Utils.ensureTrue(false, "ERROR: printSegmentedRegionsToFile(): Impossible state!");
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
	/** The continuation of the @method segmentRegionsAllFiles method, but by individual. */ 
	public static CopyNumberRegionsByChromosome segmentRegionsOneFile(File inFile, String outDir, SNVMap snvMap) {
		FileExtensionAndDelimiter fileExtAndDelim = Utils.FileExtensionTSV;		
		
		// First check that the file is a file of the desired extension		
		int indexOfDelimiter = inFile.getName().indexOf(fileExtAndDelim.mExtension); 
		if (indexOfDelimiter < 0) return null;
		if (indexOfDelimiter != inFile.getName().length() - fileExtAndDelim.mExtension.length()) return null;
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
			final String rsColumnValue =                              (Utils.extractNthColumnValue(line, ColCuratedTSV_dbSNP,    Utils.TabStr));
			final int rsId = (rsColumnValue.indexOf(NovelStr) >= 0 || rsColumnValue.indexOf(Utils.NAStr) >= 0) ? 0 : GenotypeUtils.getNumberFromRsId(rsColumnValue);
			
			// Register the site in the map
			snvMap.registerSNV(chrom, position, rsId, Nuc.N, Nuc.N, true, true);
			
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
		return (clusterType != ClusterType.Noise && 
				clusterType != ClusterType.Null  && 
				clusterType != ClusterType.HETSomatic);
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
					sb.append(delimiter).append(cnrr.getChromosome().ordinal());
					sb.append(delimiter).append(cnrr.getRangeStart());
					sb.append(delimiter).append(cnrr.getRangeEnd());
					
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
				final MutationType mutationType = MutationType.getSNVType                  (Utils.extractNthColumnValue(line, ColCuratedTSV_MutationType, Utils.TabStr));
				final VariantLocation varLoc = VariantLocation.getVariantLocation(Utils.extractNthColumnValue(line, ColCuratedTSV_VariantLocation, Utils.TabStr));
				
				boolean addSite = false;
				ColorPastel theColorMutationType = ColorPastel.Pastel_Cyan;  // set as default
				
				if (clusterTypeForSite == clusterType) {					
					// Increment the event count for this chromosome and sample
					eventCount[chrom.ordinal()][sampleIndex]++;					
										
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
					//score = 300; //used with grayscale BEDs (historical artifact)
					
					if (mutationType != null) {
						switch(mutationType) {
						case NonSynonymous_SNV: theColorMutationType = ColorPastel.Yellow_Orange; break;
						case Synonymous_SNV:    theColorMutationType = ColorPastel.Light_Green;   break;						
						}
					}
				} else if (clusterTypeForSite == ClusterType.HETSomatic) {
					addSite = true;
					theColorMutationType = ColorPastel.Darker_Blue_Violet;
				}

				//if (score == greyscaleBEDScore) {
				if (addSite) {
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
		
		//String filePrefix = outDir + File.separator;
		String filePrefix = "http://ron.cs.columbia.edu/ninad/LOHcate/renal/browser_tracks/LOH/toUpload/";
		
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
	public static void getGeneEnrichment(String inDir, String outDir) {
		String outFilename = outDir + File.separator + "geneEnrichment." + Utils.FileExtensionTSV.mExtension;
		IOUtils.createDirectoryPath(outDir, false);
		File[] files = (new File(inDir)).listFiles();
		FileExtensionAndDelimiter fileExtDelim = Utils.FileExtensionTSV; 
				
		StringBuilder sb = new StringBuilder(4096);
		BucketCounter clusterTypeCountsGene = new BucketCounter(ClusterType.values().length, 0);
		
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
						clusterTypeCountsGene.clear();
						
						if (position > currentGene.mMaxBasePairPosition) { //get right-bound of gene's range of variants
							currentGene.mMaxBasePairPosition = position; 
						}
						if (position < currentGene.mMinBasePairPosition) { //...left-bound...
							currentGene.mMinBasePairPosition = position;
						}
												
						MutationType mutationType = MutationType.getSNVType(components[ColCuratedTSV_MutationType]);
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
				MutationType.NonSynonymous_SNV.toLowerCase(), 
				MutationType.Synonymous_SNV.toLowerCase(), 
				VariantLocation.Germline.toLowerCase(), 
				VariantLocation.Somatic.toLowerCase(),
				
				ClusterType.Amp.name(),
				ClusterType.LOH.name(), 
				ClusterType.LOH.name() + "_refLost", 
				ClusterType.HETGermline.name(),
				ClusterType.HETSomatic.name(),
				
				ClusterType.Amp.name() + logStr,
				ClusterType.LOH.name() + logStr, 
				ClusterType.LOH.name() + "_refLost" + logStr, 
				ClusterType.HETGermline.name() + logStr,
				ClusterType.HETSomatic.name() + logStr,
				
				ClusterType.Amp.name()              + densityStr,
				ClusterType.LOH.name()              + densityStr, 
				ClusterType.LOH.name() + "_refLost" + densityStr, 
				ClusterType.HETGermline.name()      + densityStr,
				ClusterType.HETSomatic.name()       + densityStr,
				
				ClusterType.Amp.name()         + recurrenceStr, 
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
		
		// Now write out individual gene outputs
		for (Gene gene : genes) {
			// Write samples for each gene out as well
			ArrayList<String> patientsAllClusters = new ArrayList<String>();
			for (ClusterType ct : ClusterType.values()) {
				if ((ct == ClusterType.Noise) || (ct == ClusterType.Null)) continue;
				patientsAllClusters.add(ct.name());
				patientsAllClusters.addAll(gene.getPatientsForClusterType(ct));
			}
			String geneOutFilename = outDir + File.separator + gene.mLabel + ".samples.txt";
			IOUtils.writeOutputFile(geneOutFilename, patientsAllClusters);				
		}
	}
	
	
	// ========================================================================
	// ENTRY POINT
	// ========================================================================
	public static void main(String[] args) {
		
		String nafTafInptus     = "naf_taf_inputs";
		String classifiedSites  = "sites_classified";
		String vafPlots         = "vafPlots";
		String vafWaterfallPlots = "vafWaterfallPlots";
		String regions          = "regions";
		String browserTracks    = "browser_tracks";
		String geneEnrichment   = "gene_enrichment";
		
		String taskClustering = "clustering";
		String taskRegions = "regions";
		String taskGenes = "genes";
		
		String[] argsMajor = (args.length > 0) ? (new String[] { args[0] }) : (new String[] { " " });
		JSAP jsapTask = new JSAP();
		
		String lohcateTask = "TaskName";
		FlaggedOption task = new FlaggedOption(lohcateTask).setStringParser(JSAP.STRING_PARSER).setRequired(true)
				.setShortFlag('t').setLongFlag("task").setUsageName(lohcateTask);
		task.setHelp("Indicate which task LOHcate should perform: {" + taskClustering + ", " + taskRegions + ", " + taskGenes + "}");
		ArgumentParserUtils.registerJSAPParameter(jsapTask, task);
		JSAPResult jsapResult = ArgumentParserUtils.parseAndCheck(argsMajor, jsapTask, LOHcate.class.getName());		
		String taskName = jsapResult.getString(lohcateTask);		
		
		String rootFolderPath = "RootFolderPath";
		FlaggedOption rootFolder = new FlaggedOption(rootFolderPath).setStringParser(JSAP.STRING_PARSER).setRequired(true)
				.setShortFlag(JSAP.NO_SHORTFLAG).setLongFlag("root").setUsageName(rootFolderPath);
		rootFolder.setHelp("Indicate the root directory of the data");
		ArgumentParserUtils.registerJSAPParameter(jsapTask, rootFolder);
		
		long sys_time_init = System.currentTimeMillis();	
		
		if (taskName.equals(taskClustering)) {
			String allelicBiasFile = "AllelicBiasFile";
			FlaggedOption allelicBias = new FlaggedOption(allelicBiasFile).setStringParser(JSAP.STRING_PARSER).setRequired(false)
					.setShortFlag(JSAP.NO_SHORTFLAG).setLongFlag("allelicBias").setUsageName(allelicBiasFile);
			allelicBias.setHelp("Specify a file that has allelic biases");
			ArgumentParserUtils.registerJSAPParameter(jsapTask, allelicBias);
			Clustering.registerClusteringParameters(jsapTask);
			
			jsapResult = ArgumentParserUtils.parseAndCheck(args, jsapTask, LOHcate.class.getName());			
			String rootFolderName      = jsapResult.getString(rootFolderPath);
			String allelicBiasFilename = jsapResult.getString(allelicBiasFile);
			Clustering.configureParameters(jsapResult);
			
			if (allelicBiasFilename == null) {
				Clustering.ParamsBool.IgnoreAllelicBias.setValue(true);
			}
			
			//SeqPlatform platform = SeqPlatform.getPlatform(Integer.parseInt(args[2]));
			Clustering.classifySites(rootFolderName + File.separator + nafTafInptus,
									 allelicBiasFilename,
									 rootFolderName + File.separator + classifiedSites, 
									 rootFolderName + File.separator + vafPlots, 
									 rootFolderName + File.separator + vafWaterfallPlots,
					                 SeqPlatform.Illumina); //args[2] --> 0::Illumina, 1::SOLiD
					                 
		} else if (taskName.equals(taskRegions)) {
			jsapResult = ArgumentParserUtils.parseAndCheck(args, jsapTask, LOHcate.class.getName());			
			String rootFolderName = jsapResult.getString(rootFolderPath);
			
			segmentRegionsAllFiles(rootFolderName + File.separator + classifiedSites, 
								   rootFolderName + File.separator + regions,
								   rootFolderName + File.separator + browserTracks);
			
		} else if (taskName.equals(taskGenes)) {
			jsapResult = ArgumentParserUtils.parseAndCheck(args, jsapTask, LOHcate.class.getName());			
			String rootFolderName = jsapResult.getString(rootFolderPath);
			
			getGeneEnrichment(rootFolderName + File.separator + classifiedSites, 
							  rootFolderName + File.separator + geneEnrichment);			
		}
		
		/*
		String root = args[0]; //project directory
		switch (Integer.parseInt(args[1])) { //args[1] --> 'switchboard' parameter
			case 0:				
				break;
			case 1:
				
				break;
			case 2:
				
				break;
				
			// Everything below this point is Sidd's original code	
				/*
			case 5:
				for (int i = 0; i<Enrichment.cluster_names.length - 1; i++)
					Enrichment.getPathwayEnrichment(root + "/gene_enrichment.csv", root + "/kegg_pathways_roster.tsv", root + "/kegg/pathway_enrichment/" + Enrichment.cluster_names[i] + ".csv", i);
				break;
			case 6:
				Enrichment.annotatePathways(root + "/gene_enrichment.csv", root + "/kegg_pathways_roster.tsv", root + "/KEGG");
				break;
			case 7:
				Enrichment.getGOTermCounts(root + "/gene_enrichment.csv", root + "/gene_association.goa_human", root + "/GO/counts/go_term_counts.csv");
				for (int i = 0; i<Enrichment.cluster_names.length - 1; i++) {
					Enrichment.getGOTermCounts(root + "/gene_enrichment/gene_enrichment_top_" + Enrichment.cluster_names[i] + ".csv", root + "/gene_association.goa_human", root + "/GO/counts/go_term_counts_top_" + Enrichment.cluster_names[i] + ".csv");
					Enrichment.getGOTermEnrichment(root + "/GO/counts/go_term_counts.csv", root + "/GO/counts/go_term_counts_top_" + Enrichment.cluster_names[i] + ".csv", root + "/GO/enrichment/go_term_enrichment_top_" + Enrichment.cluster_names[i] + ".csv");
					
					Enrichment.getGOTermCounts(root + "/gene_enrichment/g1/gene_enrichment_top_" + Enrichment.cluster_names[i] + ".csv", root + "/gene_association.goa_human", root + "/GO/g1/counts/go_term_counts_top_" + Enrichment.cluster_names[i] + ".csv");
					Enrichment.getGOTermEnrichment(root + "/GO/g1/counts/go_term_counts.csv", root + "/GO/g1/counts/go_term_counts_top_" + Enrichment.cluster_names[i] + ".csv", root + "/GO/g1/enrichment/go_term_enrichment_top_" + Enrichment.cluster_names[i] + ".csv");
				}
				break;
				
		}*/
		
		System.out.println("Time elapsed: " + (System.currentTimeMillis()-sys_time_init)/1000 + " seconds");
	}

}
