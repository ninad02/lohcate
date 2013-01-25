package genomeUtils;

/**
 * 
 */
import genomeEnums.Genotype;
import genomeEnums.Nuc;
import genomeEnums.TissueType;

import java.util.ArrayList;

import nutils.NumberUtils;

public class InfoOneSiteOneSample {
	
	// ========================================================================
	public InfoOneSiteSampleTissue mNormal;
	public InfoOneSiteSampleTissue mTumor;
	
	public InfoOneSiteOneSample() {
		mNormal = new InfoOneSiteSampleTissue();
		mTumor  = new InfoOneSiteSampleTissue();
		clear();
	}

	// ========================================================================
	public InfoOneSiteSampleTissue getTissueInfo(TissueType t) {
		switch(t) {
		case Normal: return mNormal; 
		case Tumor:  return mTumor;		
		}
		return null;
	}
	
	// ========================================================================
	public void clear() {
		mNormal.clear();
		mTumor.clear();
	}
	
	// ========================================================================
	public static ArrayList<InfoOneSiteOneSample> createList(int numSamples) {
		ArrayList<InfoOneSiteOneSample> iososList = new ArrayList<InfoOneSiteOneSample>(numSamples);
		for (int s = 0; s < numSamples; s++) { iososList.add(new InfoOneSiteOneSample()); }
		return iososList;
	}
					
	// ========================================================================
	// === INNER CLASS ========================================================
	// ========================================================================
	public static class InfoOneSiteSampleTissue {
		public Genotype mGenotypeCode;
		
		public int mNumReadsTotal;
		public boolean mIsAmplified;
		public InfoOneSiteSampleTissueAllele mAlleleA;
		public InfoOneSiteSampleTissueAllele mAlleleB;
		
		// ========================================================================
		public InfoOneSiteSampleTissue() {
			mAlleleA = new InfoOneSiteSampleTissueAllele();
			mAlleleB = new InfoOneSiteSampleTissueAllele();
			clear();
		}		
		
		// ========================================================================
		public void clear() {
			mGenotypeCode = Genotype.EnumInvalidGenotype;
			mNumReadsTotal = 0;		
			mIsAmplified = false;
			mAlleleA.clear();
			mAlleleB.clear();
		}

		// ========================================================================
		public double calcProportionReadsA() { return (double) mAlleleA.mNumReads / (double) calcNumReadsTotal(); }
		
		// ========================================================================
		public double calcProportionReadsB() { return (double) mAlleleB.mNumReads / (double) calcNumReadsTotal(); }
		
		// ========================================================================
		public int calcNumReadsTotal() { return (mAlleleA.mNumReads + mAlleleB.mNumReads); }
		
		// ========================================================================
		public boolean readsValidate() { return mNumReadsTotal == calcNumReadsTotal(); }
		
		// ========================================================================
		/** This checks for any zero count reads.  If so, the genotype is appropriately assigned. */ 
		public void reexamineNumReads() {
			// Now check for any zero count reads.  If there are, then the genotype will not
			// be detected as heteroyzgous.  We need to switch to a homozygous genotype.
			boolean numReadsAZero = mAlleleA.mNumReads == 0;
			boolean numReadsBZero = mAlleleB.mNumReads == 0;
			if (numReadsAZero && numReadsBZero) {
				mAlleleA.clear();
				mAlleleB.clear();
				mGenotypeCode = Genotype.EnumHomoDeletion;
			} else if (numReadsAZero) {
				mAlleleA.clear();
				mGenotypeCode = Genotype.EnumHomozygous11;
			} else if (numReadsBZero) {
				mAlleleB.clear();
				mGenotypeCode = Genotype.EnumHomozygous00;
			} else {
				mGenotypeCode = Genotype.EnumHeterozygous;
			}
		}
		
		// ========================================================================
		/** This checks which allele has the greater count.  If both alles have the same count, N is returned. */
		public Nuc getGreaterCountAllele() {
			if (mAlleleA.mNumReads == mAlleleB.mNumReads) {
				return Nuc.N;
			} else if (mAlleleA.mNumReads > mAlleleB.mNumReads) {
				return mAlleleA.mAllele;
			} else {
				return mAlleleB.mAllele;
			}
		}
	}

	// ========================================================================
	// === INNER CLASS ========================================================
	// ========================================================================
	public static class InfoOneSiteSampleTissueAllele {
		public Nuc mAllele;
		public short mNumReads;
		
		public void set(Nuc allele, short numReadsForAllele) {
			mAllele = allele;
			mNumReads = numReadsForAllele;
		}
				
		public void clear() {			
			set(Nuc.N, (short) 0);
		}
		
		public boolean hasNoReads() { return (mNumReads == 0); }
		
		/** Given a probability to retain reads, this downsamples the read count with the 
		 *  given argument, which is the probability of retaining a read.  This also enforces
		 *  a minimum read count. 
		 */
		public void downsampleCount(double probabilityToRetainRead, int minReadCount) {
			mNumReads = (short) NumberUtils.numSuccessesInTrials(mNumReads, probabilityToRetainRead, minReadCount);
		}
	}
}