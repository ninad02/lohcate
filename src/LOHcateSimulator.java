import genomeUtils.GenomeConstants;
import genomeUtils.GenotypeUtils;
import genomeUtils.InfoOneSiteOneSample;
import genomeUtils.RegionSimulator;
import genomeUtils.SeqReadSimulator;
import genomeUtils.SeqReadSimulator.SeqReadSimulationAdjustReads;
import genomeUtils.SeqReadSimulator.SeqReadSimulationParams;
import genomeEnums.Chrom;
import genomeEnums.Genotype;
import genomeEnums.Nuc;
import genomeEnums.Phase;
import genomeEnums.TissueType;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;

import lohcateEnums.ClusterType;

import org.apache.commons.math3.random.JDKRandomGenerator;
import org.apache.commons.math3.random.RandomDataGenerator;

import nutils.ArrayUtils;
import nutils.CompareUtils;
import nutils.IOUtils;
import nutils.NumberUtils;
import nutils.StringUtils;

public class LOHcateSimulator {

	// ========================================================================
	
	public static class LOHcateSimulatorParams implements SeqReadSimulator.SeqReadSimulationParams {
		public double getCoverageGenerated(TissueType t) { return 100; }
		public double getCopyNumberTotal()               { return 2.5; }			
		public double getNumFractionOfReadForGCAlleles() { return 1.0; }
		
		public boolean isCopyNumberChangeBiAllelic()     { return false; }	
		public boolean allowHemizygousOrHomozygousDeletionGenotypes() { return true; }
		public double getPurity() { return 0.80; }
		
		public int getNumCNARegions() { return 5; }
	}
	
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
	public void generateSimulatedDataForSample(LOHcateSimulatorParams simParams, Clustering.ClusteringInputOneSample oneSampleData, LOHcateSimulatorGoldStandard goldStandard) {
		
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
	private ArrayList<CopyNumberRegionRange> deduceRegions(LOHcateSimulatorParams simParams, Clustering.ClusteringInputOneSample oneSampleData) {

		boolean specifyRandomSeed = true;
		SeqReadSimulator readSimulator = specifyRandomSeed ? new SeqReadSimulator(Integer.SIZE) : new SeqReadSimulator();
		
		ArrayList<CopyNumberRegionRange> cnRegions = new ArrayList<CopyNumberRegionRange>(simParams.getNumCNARegions());
		
		// We create an array that tracks which chromosomes were used
		boolean[] chromsUsed = new boolean[Chrom.values().length];
		Arrays.fill(chromsUsed, false);
		
		int meanLengthCNARegion = 20000000;
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
	public void simulateReadsInRegion(LOHcateSimulatorParams simParams, Clustering.ClusteringInputOneSample infoSample, CopyNumberRegionRange cnRegion, LOHcateSimulatorGoldStandard goldStandard) {
		
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
			Clustering.ClusteringInputOneSite infoOneSite = infoSample.getSiteAtIndex(row);
		
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
		if (vaf <= Clustering.AlleleFrequencyStatsForSample.VAFNormalFrameLower) {
			return Genotype.EnumHomozygous00;
		} else if (vaf > Clustering.AlleleFrequencyStatsForSample.VAFNormalFrameUpper) {
			return Genotype.EnumHomozygous11;
		} else {
			return Genotype.EnumHeterozygous;
		}
	}
	
	
	// ========================================================================
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub

	}

}
