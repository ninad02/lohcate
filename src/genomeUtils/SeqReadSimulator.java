package genomeUtils;

public class SeqReadSimulator {

	/*
	// ========================================================================
	private void calculateReadsHomozygous(final InfoOneSiteOneSample iosos, ParametersSimulation simParams, final Nuc allele, final Nuc refAllele, AmpUnampStringPair ausp) {
		short numReadsNormal = 0, numReadsTumor = 0;
		
		do {
			numReadsNormal = (short) Math.round(mRandomGen.nextPoisson(simParams.getCoverageGenerated(TissueType.Normal)));
		} while (numReadsNormal < MinReadCountThreshold);
		
		double meanTumor = (simParams.getCoverageGenerated(TissueType.Tumor) / (double) GenomeConstants.DefaultDiploidCopyNumber) * simParams.mCopyNumberAmplicon; 
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
		ausp.append(allele, allele, simParams.mAmplificationBiAllelic);

		iosos.mNormal.reexamineNumReads();
		iosos.mTumor.reexamineNumReads();

	}*/
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub

	}

}
