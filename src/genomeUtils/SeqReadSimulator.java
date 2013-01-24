package genomeUtils;

import org.apache.commons.math3.random.JDKRandomGenerator;
import org.apache.commons.math3.random.RandomDataGenerator;

import genomeEnums.Genotype;
import genomeEnums.Nuc;
import genomeEnums.Phase;
import genomeEnums.TissueType;
import genomeUtils.InfoOneSiteOneSample.InfoOneSiteSampleTissue;
import nutils.NumberUtils;


public class SeqReadSimulator {
	
	// ========================================================================
	public static interface SeqReadSimulationParams {
		public double getCoverageGenerated(TissueType t);
		public double getCopyNumberTotal();				
		public double getNumFractionOfReadForGCAlleles();
		
		public boolean isCopyNumberChangeBiAllelic();	
		public boolean allowHemizygousOrHomozygousDeletionGenotypes();
	}
	
	// ========================================================================
	public static interface SeqReadSimulationGoldStandard {
		public void append(Nuc alleleAffected, Nuc alleleNotAffected, boolean doubleAmp);
		public void clear();
	}
	
	// ========================================================================
	public static interface SeqReadSimulationAdjustReads {
		public int adjustNumReadsBias(Nuc allele, int numReadsOriginal);
	}
	
	// ========================================================================
	/** A read adjuster that simply returns the number of reads passed in. */
	public static SeqReadSimulationAdjustReads ReadAdjusterNoAdjustment = new SeqReadSimulationAdjustReads() {
		public int adjustNumReadsBias(Nuc allele, int numReadsOriginal) { return numReadsOriginal; }
	};
	
	// ========================================================================
	// ========================================================================
	public static final int MinReadCountThreshold = 1;
	
	RandomDataGenerator mRandomGen;

	public SeqReadSimulator() {
		constructorCommon(NumberUtils.getRandomLong(1, Long.MAX_VALUE - 10));  // the -10 is just a buffer 
	}
	
	public SeqReadSimulator(long randomGeneratorSeed) {
		constructorCommon(randomGeneratorSeed);
	}

	private void constructorCommon(long randomGeneratorSeed) {
		JDKRandomGenerator jdkRG = new JDKRandomGenerator();
		jdkRG.setSeed(randomGeneratorSeed);
		mRandomGen = new RandomDataGenerator(jdkRG);
	}
	
	public RandomDataGenerator getRandomDataGenerator() { return mRandomGen; }
	
	// ========================================================================
	public static double calcCoverageHaploid(double coverageDiploid) { return coverageDiploid / GenomeConstants.DefaultDiploidCopyNumber; } 
	
	// ========================================================================
	/** A helper function. */
	public void calculateReads(InfoOneSiteOneSample samplePackage, SeqReadSimulationParams simParams, Nuc[] genotype, final Phase phaseToAmp, final Nuc refAllele, SeqReadSimulationGoldStandard ausp) {		
		
		// first set the sample package copy number
		samplePackage.mNormal.mIsAmplified = false;
		samplePackage.mTumor.mIsAmplified  = true;
		
		if (genotype[0] == genotype[1]) {
			calculateReadsHomozygous(samplePackage, simParams, genotype[0], refAllele, ausp);
		} else {
			calculateReadsHeterozygous(samplePackage, simParams, genotype, phaseToAmp, refAllele, ausp);
		}
	}

	// ========================================================================
	private void calculateReadsHeterozygous(InfoOneSiteOneSample samplePackage, SeqReadSimulationParams simParams, Nuc[] genotype, final Phase phaseToAmp, final Nuc refAllele, SeqReadSimulationGoldStandard ausp) {
		samplePackage.mNormal.mGenotypeCode = samplePackage.mTumor.mGenotypeCode = Genotype.EnumHeterozygous;
		
		// First set the normal  
		int indexOfReference = (refAllele == genotype[0]) ? 0 : 1;			
		Nuc otherAllele = genotype[1 - indexOfReference];
		//if (ParametersHATS.GlobalParams.shouldLog()) { SequenceLogger.debugPrint("Index of ref: " + indexOfReference + "\tPhaseToAmplify: " + phaseToAmp + "\tRefAllele:" + refAllele); }
							
		// Enter a loop to ensure that the total read count is not zero
		double normalHaploidCoverage = calcCoverageHaploid(simParams.getCoverageGenerated(TissueType.Normal));
		do {
			int numReadsA = adjustNumReadsBias(simParams, samplePackage.mNormal.mAlleleA.mAllele, Math.round(mRandomGen.nextPoisson(normalHaploidCoverage)));
			int numReadsB = adjustNumReadsBias(simParams, samplePackage.mNormal.mAlleleB.mAllele, Math.round(mRandomGen.nextPoisson(normalHaploidCoverage)));
			samplePackage.mNormal.mAlleleA.set(refAllele,   (short) numReadsA);
			samplePackage.mNormal.mAlleleB.set(otherAllele, (short) numReadsB);				
		} while (samplePackage.mNormal.calcNumReadsTotal() < MinReadCountThreshold);
		samplePackage.mNormal.reexamineNumReads();
		
		// Now set the tumor		
		double copyNumberNonAmp = (simParams.isCopyNumberChangeBiAllelic()) ? (simParams.getCopyNumberTotal() / (double) GenomeConstants.DefaultDiploidCopyNumber) : 1.0;
		double copyNumberAmp    = simParams.getCopyNumberTotal() - copyNumberNonAmp; //GenotypeUtils.calcHigherHaplotypeCopyNumber(simParams.mCopyNumberAmplicon);
		copyNumberAmp = Math.max(copyNumberAmp, 0.000000001);

		double haploidCoverage  = simParams.getCoverageGenerated(TissueType.Tumor) / (double) GenomeConstants.DefaultDiploidCopyNumber;
		double coverageTumorAmp    = haploidCoverage * copyNumberAmp;
		double coverageTumorNonAmp = haploidCoverage * copyNumberNonAmp;
		
		//if (ParametersHATS.GlobalParams.shouldLog()) {
			//SequenceLogger.debugPrint("CN:\t" + simParams.mCopyNumberAmplicon + "\tcoverageTumorAmp:\t" + coverageTumorAmp + "\tcoverageTumorNonAmp:\t" + coverageTumorNonAmp);
		//}
		
		short numReadsNonAmp = 0, numReadsAmp = 0;			
		do {
			numReadsNonAmp = (short) Math.round(mRandomGen.nextPoisson(coverageTumorNonAmp));		
			numReadsAmp    = (short) Math.round(mRandomGen.nextPoisson(coverageTumorAmp));
			
			if (!simParams.allowHemizygousOrHomozygousDeletionGenotypes()) {
				numReadsNonAmp = (short) Math.max(numReadsNonAmp, 1);
				numReadsAmp =    (short) Math.max(numReadsAmp,    1);
			}
		} while ((numReadsAmp + numReadsNonAmp) < MinReadCountThreshold);
		
		if (Phase.getPhase(indexOfReference) == phaseToAmp) {
			samplePackage.mTumor.mAlleleA.set(refAllele,   numReadsAmp);
			samplePackage.mTumor.mAlleleB.set(otherAllele, numReadsNonAmp);
			ausp.append(genotype[indexOfReference], genotype[1 - indexOfReference], simParams.isCopyNumberChangeBiAllelic());				
		} else {				
			samplePackage.mTumor.mAlleleA.set(refAllele,   numReadsNonAmp);
			samplePackage.mTumor.mAlleleB.set(otherAllele, numReadsAmp);
			ausp.append(genotype[1 - indexOfReference], genotype[indexOfReference], simParams.isCopyNumberChangeBiAllelic());
		}				

		// Now adjust for bias
		InfoOneSiteOneSample.InfoOneSiteSampleTissue infoOsst = samplePackage.mTumor;
		infoOsst.mAlleleA.mNumReads = (short) adjustNumReadsBias(simParams, infoOsst.mAlleleA.mAllele, infoOsst.mAlleleA.mNumReads);
		infoOsst.mAlleleB.mNumReads = (short) adjustNumReadsBias(simParams, infoOsst.mAlleleB.mAllele, infoOsst.mAlleleB.mNumReads);
		if (infoOsst.mAlleleA.mNumReads + infoOsst.mAlleleB.mNumReads < MinReadCountThreshold) {
			// We check if an allele is a G or C allele.  If so, then set the other allele with
			// a nonzero read count.  Thus, if the genotype is (G|C),(A|T), then the nonzero read
			// count would go to (A|T).  Note however that for genotype (G,C) or (C,G), this assigns
			// either the C or G allele (respectively) the nonzero count.			
			if (infoOsst.mAlleleA.mAllele.isGC()) {
				infoOsst.mAlleleB.mNumReads = MinReadCountThreshold;				
			} else if (infoOsst.mAlleleB.mAllele.isGC()) {
				infoOsst.mAlleleA.mNumReads = MinReadCountThreshold;
			} 
		}

		samplePackage.mTumor.reexamineNumReads();	
	}	
	
	// ========================================================================
	public void calculateReadsTissue(final InfoOneSiteSampleTissue iosst, double coverageDiploid, double copyNumberPhase0, double copyNumberPhase1, Nuc[] genotype, Nuc referenceAllele, boolean isBiAllelicChange,
			SeqReadSimulationGoldStandard goldStandard, SeqReadSimulationAdjustReads numReadAdjuster) {

		if (genotype[0] == genotype[1]) {
			calculateReadsTissueHomozygous(iosst, coverageDiploid, copyNumberPhase0, copyNumberPhase1, isBiAllelicChange, genotype[0], referenceAllele, goldStandard, numReadAdjuster);			
		} else {
			calculateReadsTissueHeterozygous(iosst, coverageDiploid, copyNumberPhase0, copyNumberPhase1, isBiAllelicChange, genotype, referenceAllele, true, goldStandard, numReadAdjuster);
			//calculateReadsHeterozygous(samplePackage, simParams, genotype, phaseToAmp, refAllele, ausp);
		}
	}
	
	// ========================================================================
	private void calculateReadsTissueHeterozygous(final InfoOneSiteSampleTissue iosst, 
												  final double coverageDiploid, 
												  final double copyNumberPhase0, final double copyNumberPhase1,
												  final boolean biAllelicChange, final Nuc[] genotype, final Nuc referenceAllele,
												  boolean allowHemizygousGenotype,
												  SeqReadSimulationGoldStandard goldStandard,
												  SeqReadSimulationAdjustReads numReadAdjuster) {
		
		iosst.mGenotypeCode = Genotype.EnumHeterozygous;
		
		// First, map the correct parameters to the correct allele
		int indexOfReference = (referenceAllele == genotype[0]) ? 0 : 1;			
		Nuc otherAllele = genotype[1 - indexOfReference];
		double copyNumRef = Math.max(copyNumberPhase0, 0.000000001);
		double copyNumAlt = Math.max(copyNumberPhase1, 0.000000001);
		if (indexOfReference == 1) {
			copyNumRef = Math.max(copyNumberPhase1, 0.000000001);
			copyNumAlt = Math.max(copyNumberPhase0, 0.000000001);
		}
		
		double haploidCoverage = calcCoverageHaploid(coverageDiploid);
		double haploidCoverageRef = haploidCoverage * copyNumRef;
		double haploidCoverageAlt = haploidCoverage * copyNumAlt;
		
		//if (ParametersHATS.GlobalParams.shouldLog()) { SequenceLogger.debugPrint("Index of ref: " + indexOfReference + "\tPhaseToAmplify: " + phaseToAmp + "\tRefAllele:" + refAllele); }
							
		// Enter a loop to ensure that the total read count is not zero		
		do {
			int numReadsA = numReadAdjuster.adjustNumReadsBias(referenceAllele, Math.round(mRandomGen.nextPoisson(haploidCoverageRef)));
			int numReadsB = numReadAdjuster.adjustNumReadsBias(otherAllele,     Math.round(mRandomGen.nextPoisson(haploidCoverageAlt)));
			iosst.mAlleleA.set(referenceAllele,   (short) numReadsA);
			iosst.mAlleleB.set(otherAllele,       (short) numReadsB);		
			
			if (!allowHemizygousGenotype) {
				numReadsA = (short) Math.max(numReadsA, 1);
				numReadsB = (short) Math.max(numReadsB, 1);
			}
			
		} while (iosst.calcNumReadsTotal() < MinReadCountThreshold);
		iosst.reexamineNumReads();
				
		//if (ParametersHATS.GlobalParams.shouldLog()) {
			//SequenceLogger.debugPrint("CN:\t" + simParams.mCopyNumberAmplicon + "\tcoverageTumorAmp:\t" + coverageTumorAmp + "\tcoverageTumorNonAmp:\t" + coverageTumorNonAmp);
		//}
	}
	
	// ========================================================================
	private void calculateReadsTissueHomozygous(final InfoOneSiteSampleTissue iosst, 
												final double coverageDiploid, 
												final double copyNumberPhase0, final double copyNumberPhase1, 
												final boolean biAllelicChange, final Nuc allele, final Nuc referenceAllele, 
												SeqReadSimulationGoldStandard goldStandard,
												SeqReadSimulationAdjustReads numReadAdjuster) {
		
		double poissonMean = calcCoverageHaploid(coverageDiploid) * (copyNumberPhase0 + copyNumberPhase1);
		short numReads = 0;

		do {
			numReads = (short) Math.round(mRandomGen.nextPoisson(poissonMean));
		} while (numReads < MinReadCountThreshold);
		
		// Adjust the reads for bias
		numReads = (short) Math.max(MinReadCountThreshold, numReadAdjuster.adjustNumReadsBias(allele, numReads));				
		
		if (referenceAllele == allele) {			
			iosst.mGenotypeCode = Genotype.EnumHomozygous00;
			iosst.mAlleleA.set(allele, numReads);					
		} else {
			iosst.mGenotypeCode = Genotype.EnumHomozygous11;
			iosst.mAlleleB.set(allele, numReads);												
		}			
		
		if (goldStandard != null) {
			goldStandard.append(allele, allele, biAllelicChange);
		}
		
		iosst.reexamineNumReads();
	}
	
	// ========================================================================
	private void calculateReadsHomozygous(final InfoOneSiteOneSample iosos, SeqReadSimulationParams simParams, final Nuc allele, final Nuc refAllele, SeqReadSimulationGoldStandard goldStandard) {
		short numReadsNormal = 0, numReadsTumor = 0;
		
		do {
			numReadsNormal = (short) Math.round(mRandomGen.nextPoisson(simParams.getCoverageGenerated(TissueType.Normal)));
		} while (numReadsNormal < MinReadCountThreshold);
		
		double meanTumor = calcCoverageHaploid(simParams.getCoverageGenerated(TissueType.Tumor)) * simParams.getCopyNumberTotal(); 
		do {			
			numReadsTumor  = (short) Math.round(mRandomGen.nextPoisson(meanTumor));
		} while (numReadsTumor < MinReadCountThreshold);
		
		// Adjust the reads for bias
		numReadsNormal = (short) Math.max(MinReadCountThreshold, adjustNumReadsBias(simParams, allele, numReadsNormal));
		numReadsTumor  = (short) Math.max(MinReadCountThreshold, adjustNumReadsBias(simParams, allele, numReadsTumor ));
		
		if (refAllele == allele) {					
			iosos.mNormal.mGenotypeCode = iosos.mTumor.mGenotypeCode = Genotype.EnumHomozygous00;
			iosos.mNormal.mAlleleA.set(allele, numReadsNormal);
			iosos.mTumor.mAlleleA.set (allele, numReadsTumor);			
		} else {
			iosos.mNormal.mGenotypeCode = iosos.mTumor.mGenotypeCode = Genotype.EnumHomozygous11;
			iosos.mNormal.mAlleleB.set(allele, numReadsNormal);
			iosos.mTumor.mAlleleB.set (allele, numReadsTumor);										
		}			
		goldStandard.append(allele, allele, simParams.isCopyNumberChangeBiAllelic());

		iosos.mNormal.reexamineNumReads();
		iosos.mTumor.reexamineNumReads();

	}

	// ========================================================================
	private static int adjustNumReadsBias(SeqReadSimulationParams simParams, Nuc allele, int numReadsIntended) {
		
		// First test if there is a need to adjust the read counts
		if ((simParams.getNumFractionOfReadForGCAlleles() >= 0.9999) || !allele.isGC()) {
			return numReadsIntended;
		}
		
		// Now we're guaranteed we have a GC allele and there is a fraction to downsample
		return NumberUtils.numSuccessesInTrials(numReadsIntended, simParams.getNumFractionOfReadForGCAlleles(), 0);	
	}
	
	// ========================================================================
	
	private static void TestReadGeneration() {
		InfoOneSiteSampleTissue iosst = new InfoOneSiteSampleTissue();
		SeqReadSimulator srs = new SeqReadSimulator();
		
		int numIter = 10000;
		double coverageDiploid = 100;
		double copyNumberPhase0 = 1.5;
		double copyNumberPhase1 = 1;		
		Nuc referenceAllele = Nuc.A;
		Nuc[] genotype = new Nuc[] { Nuc.A, Nuc.C };
		boolean isBiAllelicChange = false;
		
		int totalReadA = 0;
		int totalReadB = 0;
		
		for (int i = 0; i < numIter; i++) {
			iosst.clear();
			srs.calculateReadsTissue(iosst, coverageDiploid, copyNumberPhase0, copyNumberPhase1, genotype, referenceAllele, isBiAllelicChange, null, ReadAdjusterNoAdjustment);
			totalReadA += iosst.mAlleleA.mNumReads;
			totalReadB += iosst.mAlleleB.mNumReads;
			System.out.println(iosst.mGenotypeCode);
		}		
		System.out.println((totalReadA / (double) numIter) + "\t" + (totalReadB / (double) numIter));
	}
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		TestReadGeneration();
	}

}
