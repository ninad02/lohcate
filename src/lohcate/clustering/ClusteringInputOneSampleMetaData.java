package lohcate.clustering;

import genomeEnums.Chrom;

import java.io.PrintWriter;
import java.util.Arrays;

import lohcate.Script;
import nutils.PrimitiveWrapper;
import nutils.counter.DynamicBucketCounter;

// ========================================================================
// ========================================================================
public class ClusteringInputOneSampleMetaData {

	float[]  mAdjustedVAFNormal;
	float[]  mAdjustedVAFTumor;
	double[] mImbalancePValuesTumor;
	double[] mImbalancePValuesNormal;		
	
	PrimitiveWrapper.WFloat mCoverageRatioTumorToNormal;
	DynamicBucketCounter mReadCountTalliesTumor;
	DynamicBucketCounter mReadCountTalliesNormal;
	
	int[]   mNumSitesPerChrom;
	float[] mAvgReadCountPerChromNormal;
	float[] mCopyNumRatioPerChromNormal;
	float[] mTumorCopyNumRatiosPerGene;	
	boolean[] mIsSomaticSite;
	
	double mFDRNormal = 0;
	double mFDRTumor = 0;
	
	public ClusteringInputOneSampleMetaData(int numSites) {
		mAdjustedVAFNormal         = new float[numSites];
		mAdjustedVAFTumor          = new float[numSites];
		mImbalancePValuesTumor    = new double[numSites];
		mImbalancePValuesNormal   = new double[numSites];
		mTumorCopyNumRatiosPerGene = new float[numSites];
		mIsSomaticSite             = new boolean[numSites];
		
		mCoverageRatioTumorToNormal = new PrimitiveWrapper.WFloat(0);
		mReadCountTalliesTumor  = new DynamicBucketCounter();
		mReadCountTalliesNormal = new DynamicBucketCounter();
		
		mNumSitesPerChrom           = new   int[Chrom.values().length];
		mAvgReadCountPerChromNormal = new float[Chrom.values().length];
		mCopyNumRatioPerChromNormal = new float[Chrom.values().length];		
		clear();
	}		
	
	public void clear() {		
	}
	
	public boolean chromHasGermlineGain(Chrom chrom) {
		return (mCopyNumRatioPerChromNormal[chrom.ordinal()] > ClusteringParams.GlobalClusteringParams.mGermlineTrisomyThreshold.getValue());
	}
	
	public void printSiteInformation(PrintWriter out, int row, boolean printNewLine) {
		out.printf(  "%g\t%g", mTumorCopyNumRatiosPerGene[row], (mTumorCopyNumRatiosPerGene[row] * Script.DefaultDiploidCopyNumber));  
		out.printf("\t%g\t%g", mAdjustedVAFTumor[row], mAdjustedVAFNormal[row]);
		out.printf("\t%g\t%g", mImbalancePValuesTumor[row], mImbalancePValuesNormal[row]);
		if (printNewLine) {
			out.println("");
		}			
	}
}