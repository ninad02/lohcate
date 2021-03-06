package lohcate.clustering;

import genomeEnums.Chrom;
import java.io.PrintWriter;
import java.util.Arrays;

import lohcate.CopyNumberRegionsByChromosome;
import lohcate.Regions;
import nutils.Cast;
import nutils.RangeDouble;
import nutils.counter.DynamicBucketCounter;
import nutils.primitives.wrapper.PrimitiveWrapper;

// ========================================================================
// ========================================================================
public class ClusteringInputOneSampleMetaData {

	float[]  mAdjustedVAFNormal;
	float[]  mAdjustedVAFTumor;
	double[] mImbalancePValuesTumor;
	double[] mImbalancePValuesNormal;
	double[] mSigPValuesPerSite;
	
	PrimitiveWrapper.WFloat mCoverageRatioTumorToNormal;
	DynamicBucketCounter mReadCountTalliesTumor;
	DynamicBucketCounter mReadCountTalliesNormal;
	
	float mAvgCoverageNormal;
	float mAvgCoverageTumor;
	
	int[]   mNumSitesPerChrom;
	float[] mAvgReadCountPerChromNormal;
	float[] mCopyNumRatioPerChromNormal;
	private float[] mTumorCopyNumRatiosPerGene;	
	boolean[] mIsSomaticSite;
	boolean[] mChromHasGermlineGain;
	boolean[] mPossibleSampleContamination;
	
	CopyNumberRegionsByChromosome mGeneRegions;
	RangeDouble mVAFNormalHetRange;
	
	double mFDRNormal = 0;
	double mFDRTumor = 0;
	
	public ClusteringInputOneSampleMetaData(int numSites, String sampleName) {
		mAdjustedVAFNormal         = new float[numSites];
		mAdjustedVAFTumor          = new float[numSites];
		mImbalancePValuesTumor    = new double[numSites];
		mImbalancePValuesNormal   = new double[numSites];
		mTumorCopyNumRatiosPerGene = new float[numSites];
		mIsSomaticSite             = new boolean[numSites];
		mSigPValuesPerSite        = new double[numSites];
		
		mCoverageRatioTumorToNormal = new PrimitiveWrapper.WFloat(0);
		mReadCountTalliesTumor  = new DynamicBucketCounter();
		mReadCountTalliesNormal = new DynamicBucketCounter();
		mGeneRegions = new CopyNumberRegionsByChromosome(sampleName);
		
		mNumSitesPerChrom           = new   int[Chrom.values().length];
		mAvgReadCountPerChromNormal = new float[Chrom.values().length];
		mCopyNumRatioPerChromNormal = new float[Chrom.values().length];
		mChromHasGermlineGain       = new boolean[Chrom.values().length];
		mPossibleSampleContamination = new boolean[Chrom.values().length];
		
		Arrays.fill(mChromHasGermlineGain, false);
		
		mVAFNormalHetRange = new RangeDouble(AlleleFractionStatsForSample.VAFNormalFrameLower, AlleleFractionStatsForSample.VAFNormalFrameUpper);
		
		clear();
	}		
	
	public void clear() {		
	}
	
	public boolean chromHasGermlineGain(Chrom chrom) {
		//return (mCopyNumRatioPerChromNormal[chrom.ordinal()] > ClusteringParams.GlobalClusteringParams.mGermlineTrisomyThreshold.getValue());
		return mChromHasGermlineGain[chrom.ordinal()];
	}
	
	public float getTumorNormalCopyNumRatioAtIndex(int row) {
		return mTumorCopyNumRatiosPerGene[row];
	}
	
	public void setTumorNormalCopyNumRatioAtIndex(int row, double copyNumRatio) {
		mTumorCopyNumRatiosPerGene[row] = Cast.toFloat(copyNumRatio);
	}

	
	public float getCopyNumberAtIndex(int row) {
		return mTumorCopyNumRatiosPerGene[row] * Regions.DefaultDiploidCopyNumber;
	}
	
	public void setCopyNumberAtIndex(int row, double copyNum) {
		mTumorCopyNumRatiosPerGene[row] = Cast.toFloat(copyNum / Regions.DefaultDiploidCopyNumber);
	}
	
	public void adjustTumorCopyNumRatio(int row, double diffFactor) {		
		mTumorCopyNumRatiosPerGene[row] += Cast.toFloat(diffFactor);
	}
		
	public void setCopyNumberAllSites(double copyNum) {
		for (int i = 0; i < mTumorCopyNumRatiosPerGene.length; i++) {
			setCopyNumberAtIndex(i, copyNum);
		}
	}
	
	public int getNumSites() { return mTumorCopyNumRatiosPerGene.length; }
	
	public void printSiteInformation(PrintWriter out, int row, boolean printNewLine) {
		out.printf(  "%g\t%g", mTumorCopyNumRatiosPerGene[row], getCopyNumberAtIndex(row));  
		out.printf("\t%g\t%g", mAdjustedVAFTumor[row], mAdjustedVAFNormal[row]);
		out.printf("\t%g\t%g", mImbalancePValuesTumor[row], mImbalancePValuesNormal[row]);
		out.printf("\t%g", mSigPValuesPerSite[row]);
		if (printNewLine) {
			out.println("");
		}			
	}
}