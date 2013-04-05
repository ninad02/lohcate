package lohcate;
import genomeUtils.GenomeConstants;
import genomeUtils.GenotypeUtils;
import genomeUtils.InfoOneSiteOneSample;
import genomeUtils.RegionSimulator;
import genomeUtils.SeqReadSimulator;
import genomeEnums.Chrom;
import genomeEnums.Genotype;
import genomeEnums.Nuc;
import genomeEnums.Phase;
import genomeEnums.TissueType;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;

import lohcate.clustering.AlleleFractionStatsForSample;
import lohcate.clustering.ClusteringInputOneSample;
import lohcate.clustering.ClusteringInputOneSite;
import lohcateEnums.ClusterType;

import org.apache.commons.math3.random.JDKRandomGenerator;
import org.apache.commons.math3.random.RandomDataGenerator;

import com.carrotsearch.hppc.IntArrayList;
import com.martiansoftware.jsap.JSAP;

import nutils.ArgumentParserUtils.InputParameter;
import nutils.ArgumentParserUtils.InputParameterDouble;
import nutils.ArgumentParserUtils.InputParameterInteger;
import nutils.ArrayUtils;
import nutils.CompareUtils;
import nutils.IOUtils;
import nutils.NumberUtils;
import nutils.StringUtils;
import nutils.counter.DynamicBucketCounter;

public class LOHcateSimulator {

	// ========================================================================
	// INNER CLASS
	// ========================================================================
	public static class LOHcateSimulatorParams implements SeqReadSimulator.SeqReadSimulationParams {

		// ====================================================================
		// MEMBER VARIABLES
		// ====================================================================
		
		protected InputParameterDouble[] mCoverageGenerated = new InputParameterDouble[] {
			new InputParameterDouble(100.0, "CoverageGeneratedTumor", JSAP.NO_SHORTFLAG, "CoverageGeneratedTumor", JSAP.NO_DEFAULT),
			new InputParameterDouble(100.0, "CoverageGeneratedNormal", JSAP.NO_SHORTFLAG, "CoverageGeneratedNormal", JSAP.NO_DEFAULT)
		};				
				
		protected InputParameterDouble[] mCoverageExpected = new InputParameterDouble[] {
				new InputParameterDouble(100.0, "CoverageExpectedTumor",  JSAP.NO_SHORTFLAG, "CoverageExpectedTumor", JSAP.NO_DEFAULT),
				new InputParameterDouble(100.0, "CoverageExpectedNormal", JSAP.NO_SHORTFLAG, "CoverageExpectedNormal", JSAP.NO_DEFAULT)
		};	
		
		protected InputParameterDouble mTumorPurity = new InputParameterDouble(0.40, "TumorPurity", JSAP.NO_SHORTFLAG, "TumorPurity", JSAP.NO_DEFAULT);		
		protected InputParameterInteger mNumCNARegions = new InputParameterInteger(5, "NumCNARegions", JSAP.NO_SHORTFLAG, "NumCNARegions", JSAP.NO_DEFAULT);
		
		protected InputParameterInteger mNumIterations = new InputParameterInteger(1, "NumIterations", JSAP.NO_SHORTFLAG, "NumIterations", JSAP.NO_DEFAULT);
		
		public ArrayList<InputParameter<?>> mParams;
		
		// ====================================================================
		public LOHcateSimulatorParams() {			
			mParams = new ArrayList<InputParameter<?>>();
			for (TissueType t : TissueType.values()) {
				mParams.add(mCoverageGenerated[t.ordinal()]);
				mParams.add(mCoverageExpected[t.ordinal()]);
			}
			mParams.add(mTumorPurity);
			mParams.add(mNumCNARegions);
			
		}
		
		// ====================================================================
		public LOHcateSimulatorParams(String paramFilename) {
			this();  // call other constructor
			parseParamFileAndSet(paramFilename);
		}
		
		// ====================================================================
		protected void parseParamFileAndSet(String paramFilename) {
			String[][] params = IOUtils.readAllLinesFromFileAsMatrix(paramFilename);
			for (String[] paramAndArgs : params) {
				String paramName = paramAndArgs[0];
				String argValue  = paramAndArgs[1];
				
				for (InputParameter<?> ip : mParams) {
					if (ip.getLongFlag().equals(paramName)) {
						ip.parseValue(argValue);
					}
				}
			}
		}
		
		// ========================================================================
		public void printValues(PrintStream out) {
			out.println("------------------------------------------------");
			out.println("------------ Simulation Parameters -------------");

			for (InputParameter<?> ip : mParams) {
				String line = ip.getNameAndValueAsString(StringUtils.FileExtensionTSV.mDelimiter);
				out.println(line);
			}			
			
			out.println("------------------------------------------------");
			out.flush();
		}
		
		// ====================================================================
		// GETTERS
		public double getCoverageGenerated(TissueType t) { return mCoverageGenerated[t.ordinal()].getValue(); }
		
		public double getPurity() { return mTumorPurity.getValue(); }
		
		public int getNumCNARegions() { return mNumCNARegions.getValue(); } 

		public int getNumIterations() { return mNumIterations.getValue(); }
		
		// ====================================================================		
		public double getNumFractionOfReadForGCAlleles() { return 1.0; }
		
		public double getCopyNumberTotal()               { return 2.5; }	
		public boolean isCopyNumberChangeBiAllelic()     { return false; }
		public boolean allowHemizygousOrHomozygousDeletionGenotypes() { return true; }
	}
	
	// ========================================================================
	// ======================= END INNER CLASS ================================
	// ========================================================================

	
	// ========================================================================
	// INNER CLASS
	// ========================================================================
	public static class LOHcateSimulatorGoldStandard implements SeqReadSimulator.SeqReadSimulationGoldStandard {
		
		ArrayList<ClusterType> mSomaticEvents;
		
		public LOHcateSimulatorGoldStandard(int numSites) {
			mSomaticEvents = new ArrayList<ClusterType>(numSites);
			setNumSites(numSites);
		}
		
		public void append(Nuc alleleAffected, Nuc alleleNotAffected, boolean doubleAmp) {}
		
		public void setNumSites(int numSites) {
			clear();
			mSomaticEvents.ensureCapacity(numSites);
			ArrayUtils.addToCollection(mSomaticEvents, ClusterType.Null, numSites, true);
		}
		
		public void setEventAtSite(int siteIndex, ClusterType eventType) { mSomaticEvents.set(siteIndex, eventType); }
		
		public void clear() {
			mSomaticEvents.clear();
		}
		
		public ClusterType getEvent(int index) { return mSomaticEvents.get(index); }
	}

	// ========================================================================
	RandomDataGenerator mRandomGen;
	LOHcateSimulatorParams mParamsSimulation;
	
	public LOHcateSimulator() {
		constructorCommon(NumberUtils.getRandomLong(1, Long.MAX_VALUE - 10));  // the -10 is just a buffer 
	}
	
	public LOHcateSimulator(long randomGeneratorSeed) {
		constructorCommon(randomGeneratorSeed);
	}

	private void constructorCommon(long randomGeneratorSeed) {
		JDKRandomGenerator jdkRG = new JDKRandomGenerator();
		jdkRG.setSeed(randomGeneratorSeed);
		mRandomGen = new RandomDataGenerator(jdkRG);
		mParamsSimulation = new LOHcateSimulatorParams();
	}

	// ========================================================================
	public void generateSimulatedDataForSample(LOHcateSimulatorParams simParams, ClusteringInputOneSample oneSampleData, LOHcateSimulatorGoldStandard goldStandard) {
		
		// First initialize the sample by simulating reads across the sample
		simulateReadsInRegion(simParams, oneSampleData, null, goldStandard);
		
		// Get the regions
		ArrayList<CopyNumberRegionRange> regionsSelected = deduceRegions(simParams, oneSampleData);
		
		// Go through the regions and generate the data
		for (CopyNumberRegionRange regionSelected : regionsSelected) {
			simulateReadsInRegion(simParams, oneSampleData, regionSelected, goldStandard);
		}
	}
	
	// ========================================================================
	private ArrayList<CopyNumberRegionRange> deduceRegions(LOHcateSimulatorParams simParams, ClusteringInputOneSample oneSampleData) {

		boolean specifyRandomSeed = true;
		SeqReadSimulator readSimulator = specifyRandomSeed ? new SeqReadSimulator(Integer.SIZE) : new SeqReadSimulator();
		
		ArrayList<CopyNumberRegionRange> cnRegions = new ArrayList<CopyNumberRegionRange>(simParams.getNumCNARegions());
		
		// We create an array that tracks which chromosomes were used
		boolean[] chromsUsed = new boolean[Chrom.values().length];
		Arrays.fill(chromsUsed, false);
		
		int meanLengthCNARegion = 80000000;
		for (int regionIndex = 0; regionIndex < simParams.getNumCNARegions(); regionIndex++) {
			int regionLength = (int) readSimulator.getRandomDataGenerator().nextExponential(meanLengthCNARegion);
			//SequenceLogger.outputPrintln("Region Length: " + regionLength);
			// Get random event type
			int randomEventIndex = NumberUtils.getRandomInteger(0, ClusterType.AmpLOHcnLOH.length - 1);
			CopyNumberRegionRange newCNRegion = new CopyNumberRegionRange(ClusterType.AmpLOHcnLOH[randomEventIndex], Chrom.c0, 0);
			
			RegionSimulator.generateRegion(regionLength, newCNRegion, cnRegions, oneSampleData);
			chromsUsed[newCNRegion.getChromosome().ordinal()] = true;  // Set this this chromosome was used
		}
		
		// Now we want to assign an entire chromosome to undergo aneuploidy gain in the germline.  We pick
		// a chromosome that hasn't been selected yet and has 2+ sites to house a CN region.
		Chrom germlineGainChrom = null;
		
		while (true) {
			germlineGainChrom = Chrom.getChrom((byte) NumberUtils.getRandomInteger(Chrom.IndexAutosomalStart, Chrom.IndexAutosomalEnd));
			if (!chromsUsed[germlineGainChrom.ordinal()]) {
				int indexChromStart = oneSampleData.getIndexChromStart(germlineGainChrom);
				if (indexChromStart >= 0) {
					int indexChromEnd = oneSampleData.getIndexChromEnd(germlineGainChrom);
					if (indexChromEnd - indexChromStart > 1) {
						CopyNumberRegionRange newCNRegion = 
								new CopyNumberRegionRange(ClusterType.GainGermline, germlineGainChrom, 
										oneSampleData.getSiteAtIndex(germlineGainChrom, indexChromStart).getPosition(),
										oneSampleData.getSiteAtIndex(germlineGainChrom, indexChromEnd).getPosition());
						cnRegions.add(newCNRegion);
						break;
					}
				}
			}
		} 
		
		
		
		// 
		// Go through the regions and assign 
		//Print out the regions
		for (CopyNumberRegionRange region : cnRegions) {
			System.out.println(region + "\t" + region.mCopyNumberClusterType);
		}
		
		// Generate the filename and the outstream		
		//PrintStream outStream = IOUtils.getPrintStream(sep.mOutFilename);
		//simParams.printValues(outStream);
		return cnRegions;		
	}
	
//private void simulateAndEvaluate(SimEvalParams sep, ParametersSimulation simParams) {
		
				
		//if (ParametersHATS.GlobalParams.shouldLog()) {	
		//	SequenceLogger.finalizeLogger();
		//}
		
		// Flush and close the output stream
		//outStream.flush();
		//outStream.close();
	
	// ========================================================================	
	public void simulateReadsInRegion(LOHcateSimulatorParams simParams, ClusteringInputOneSample infoSample, CopyNumberRegionRange cnRegion, LOHcateSimulatorGoldStandard goldStandard) {
		
		Nuc[] genotype = new Nuc[2];
		SeqReadSimulator seqSim = new SeqReadSimulator();
		InfoOneSiteOneSample iosos = new InfoOneSiteOneSample();
		SeqReadSimulator.SeqReadSimulationAdjustReads readAdjuster = SeqReadSimulator.ReadAdjusterNoAdjustment; 		
		StringBuilder sb = new StringBuilder(2048);
		
		// Set the indices to default values (to cover the entire sample, unless a region is specified
		int indexStart = 0;
		int indexEnd   = infoSample.getNumSites() - 1;
				
		// If a region is specified, this limits the simulation only to the specified region
		if (cnRegion != null) {
			indexStart = infoSample.getIndex(cnRegion.getChromosome(), cnRegion.getRangeStart());
			if (indexStart < 0) {
				CompareUtils.ensureTrue(false, "ERROR: Starting index must exist! Chrom:\t" + cnRegion.getChromosome() + "\tPos:\t" + cnRegion.getRangeStart());
			}

			indexEnd = infoSample.getIndex(cnRegion.getChromosome(), cnRegion.getRangeEnd());
			if (indexEnd < 0) {
				CompareUtils.ensureTrue(false, "ERROR: Starting index must exist! Chrom:\t" + cnRegion.getChromosome() + "\tPos:\t" + cnRegion.getRangeEnd());
			}
		}
		
		boolean allowHemizygousGenotype = true;
		
		// Now go through the relevant region
		for (int row = indexStart; row <= indexEnd; row++) {
			iosos.clear();
			ClusteringInputOneSite infoOneSite = infoSample.getSiteAtIndex(row);
		
			Nuc referenceAllele     = infoOneSite.getReferenceAllele();
			Nuc variantAlleleNormal = infoOneSite.getVariantAlleleNormal();
			Nuc variantAlleleTumor  = infoOneSite.getVariantAlleleTumor();
			
			Genotype genotypeEnum = deduceGenotype(infoOneSite.calcVAFNormal());
			GenotypeUtils.defineGenotypeAlleles(genotypeEnum, referenceAllele, variantAlleleNormal, genotype);
			ClusterType eventTypeToAssign = (genotypeEnum == Genotype.EnumHeterozygous ? ClusterType.HETGermline : ClusterType.Null);			
			
			//if (ParametersHATS.GlobalParams.shouldLog()) { SequenceLogger.debugPrint(infoOneSite.getPosition() + "\t" + refAllele); }
//			if (ParametersHATS.GlobalParams.shouldLog()) {
//				SequenceLogger.debugPrint("Gen: " + genotype[0].getChar() + "/" + genotype[1].getChar());
//			}

			double[][] copyNumber = ArrayUtils.newDoubleArray2D(TissueType.numValid(), Phase.numValid(), GenomeConstants.DefaultHaploidCopyNumber);
			boolean isBiAllelicChangeNormal = false;
			boolean isBiAllelicChangeTumor  = false;
			Phase phaseToChange = NumberUtils.getRandomBit() ? Phase.p0 : Phase.p1;
			Phase otherPhase = (phaseToChange == Phase.p0) ? Phase.p1 : Phase.p0;
			
			if (cnRegion == null) {
				// Do nothing, we're just using as a placeholder.
			
			} else if (cnRegion.mCopyNumberClusterType == ClusterType.GainGermline) {
				copyNumber[TissueType.Normal.mCode][phaseToChange.mCode] = ++(copyNumber[TissueType.Tumor.mCode][phaseToChange.mCode]);
				eventTypeToAssign = ClusterType.GainGermline;
						
			} else if (cnRegion.mCopyNumberClusterType == ClusterType.GainSomatic) {
				copyNumber[TissueType.Tumor.mCode][phaseToChange.mCode] = GenotypeUtils.adjustHaploidCopyNumber(GenomeConstants.DefaultDiploidCopyNumber, simParams.getPurity());				
				eventTypeToAssign = ClusterType.GainSomatic;
				
			} else if (cnRegion.mCopyNumberClusterType == ClusterType.LOH) {				
				copyNumber[TissueType.Tumor.mCode][phaseToChange.mCode] = GenotypeUtils.adjustHaploidCopyNumber(0, simParams.getPurity());				
				eventTypeToAssign = ClusterType.LOH;
				
			} else if (cnRegion.mCopyNumberClusterType == ClusterType.cnLOH) {
				copyNumber[TissueType.Tumor.mCode][phaseToChange.mCode] = GenotypeUtils.adjustHaploidCopyNumber(GenomeConstants.DefaultDiploidCopyNumber, simParams.getPurity());
				copyNumber[TissueType.Tumor.mCode][otherPhase.mCode]    = GenotypeUtils.adjustHaploidCopyNumber(0,                                        simParams.getPurity());
				eventTypeToAssign = ClusterType.cnLOH;
			}
			
			// Finally do the gold standard assignment
			goldStandard.setEventAtSite(row, eventTypeToAssign);
			
			seqSim.calculateReadsTissue(iosos.getTissueInfo(TissueType.Normal), simParams.getCoverageGenerated(TissueType.Normal), 
					copyNumber[TissueType.Normal.mCode][Phase.p0.mCode], copyNumber[TissueType.Normal.mCode][Phase.p1.mCode], 
					genotype, referenceAllele, allowHemizygousGenotype, isBiAllelicChangeNormal, null, readAdjuster);
			
			seqSim.calculateReadsTissue(iosos.getTissueInfo(TissueType.Tumor), simParams.getCoverageGenerated(TissueType.Tumor), 
					copyNumber[TissueType.Tumor.mCode][Phase.p0.mCode], copyNumber[TissueType.Tumor.mCode][Phase.p1.mCode], 
					genotype, referenceAllele, allowHemizygousGenotype, isBiAllelicChangeTumor, goldStandard, readAdjuster);
			
			
			// Finally, register the read counts
			infoOneSite.setCovgVarNormal(iosos.mNormal.mAlleleB.mNumReads);
			infoOneSite.setCovgVarTumor ( iosos.mTumor.mAlleleB.mNumReads);
			infoOneSite.setCovgTotalNormal((short) iosos.mNormal.calcNumReadsTotal());
			infoOneSite.setCovgTotalTumor ((short)  iosos.mTumor.calcNumReadsTotal());
			
			// System.out.println(infoOneSite.printToString(sb, true, StringUtils.FileExtensionTSV.mDelimiter) + "\t" + eventTypeToAssign);
		}
	}
	
	// ========================================================================
	/** Based on the variant allele frequency, this determines the genotype. */
	private Genotype deduceGenotype(double vaf) {
		if (vaf <= AlleleFractionStatsForSample.VAFNormalFrameLower) {
			return Genotype.EnumHomozygous00;
		} else if (vaf > AlleleFractionStatsForSample.VAFNormalFrameUpper) {
			return Genotype.EnumHomozygous11;
		} else {
			return Genotype.EnumHeterozygous;
		}
	}

	// ========================================================================
	// Test function.  Tests whether distributions of reads differ between separate
	// copy number and aggregated copy number
	private static void TestSeparateVersusAggregated() {
		LOHcateSimulatorParams params = new LOHcateSimulatorParams();
		params.mTumorPurity.setValue(0.5);
		double haploidCoverage = 50;
		int numReadIterations = 1000000; 
		
		// First the joint		
		double jointCopyNum = 1.0;
		double jointCopyNumTumor = GenotypeUtils.adjustDiploidCopyNumber(jointCopyNum, params.mTumorPurity.getValue());			
		DynamicBucketCounter jointCounts = new DynamicBucketCounter();
		
		// Now the separate
		double sepCopyNumATumor = 1.0;
		double sepCopyNumBTumor = 0.0;
		double sepCopyNumAStroma = 1.0;
		double sepCopyNumBStroma = 1.0;
		DynamicBucketCounter sepCounts = new DynamicBucketCounter();
		
		LOHcateSimulator simulator =  new LOHcateSimulator();
		
		// For the joint, generate
		for (int iter = 0; iter < numReadIterations; iter++) {
			// First, the joint
			int jointNumReads = (int) simulator.mRandomGen.nextPoisson(jointCopyNumTumor * haploidCoverage);
			jointCounts.incrementCount(jointNumReads);
			
			// Now, do seperate
			int sepReadsATumor =  (int) simulator.mRandomGen.nextPoisson(sepCopyNumATumor  * haploidCoverage);
			int sepReadsBTumor =  0; //(int) simulator.mRandomGen.nextPoisson(sepCopyNumBTumor  * haploidCoverage);
			int sepReadsAStroma = (int) simulator.mRandomGen.nextPoisson(sepCopyNumAStroma * haploidCoverage);
			int sepReadsBStroma = (int) simulator.mRandomGen.nextPoisson(sepCopyNumBStroma * haploidCoverage);
			
			sepReadsATumor  = NumberUtils.numSuccessesInTrials(sepReadsATumor,      params.mTumorPurity.getValue(), 0);
			sepReadsBTumor  = NumberUtils.numSuccessesInTrials(sepReadsBTumor,      params.mTumorPurity.getValue(), 0);
			sepReadsAStroma = NumberUtils.numSuccessesInTrials(sepReadsAStroma, 1 - params.mTumorPurity.getValue(), 0);
			sepReadsBStroma = NumberUtils.numSuccessesInTrials(sepReadsBStroma, 1 - params.mTumorPurity.getValue(), 0);
			int sepTotalReads = sepReadsATumor + sepReadsBTumor + sepReadsAStroma + sepReadsBStroma;
			sepCounts.incrementCount(sepTotalReads);
		}
		
		IntArrayList[] jointCountsArray = jointCounts.toArrayListInt();		
		for (int i = 0; i < jointCountsArray[0].size(); i++) {
			System.out.printf("%d\t%d\n", jointCountsArray[0].get(i), jointCountsArray[1].get(i));
		}
		
		IntArrayList[] sepCountsArray = sepCounts.toArrayListInt();
		for (int i = 0; i < sepCountsArray[0].size(); i++) {
			System.out.printf("%d\t%d\n", sepCountsArray[0].get(i), sepCountsArray[1].get(i));
		}

	}
	
	// ========================================================================
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		TestSeparateVersusAggregated();
	}

}
