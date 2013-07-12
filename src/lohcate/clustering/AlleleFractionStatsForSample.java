package lohcate.clustering;

import genomeUtils.SiteList;

import java.util.Arrays;

import lohcateEnums.SeqPlatform;
import nutils.RangeDouble;

// ========================================================================
// INNER CLASS
// ========================================================================
/** An inner class that calculates statistics for a set of sites (with variant 
 *  allele frequencies information for tumor and matched normal). */
public class AlleleFractionStatsForSample {
	
	//NAF {het, loh, dup} FRAME definition via peak detection and parameter-tuned standard deviation expansion
	//we have to avoid the often hugely dense peak of homozygous mutations (AF > 0.8) and the occasionally hugely dense peak of neg. tail noise / somatics / &c. (AF < 0.2)
	public static final float VAFNormalFrameLower = Clustering.correctAllelicBias() ? 0.1f : 0.2f;
	public static final float VAFNormalFrameUpper = Clustering.correctAllelicBias() ? 0.9f : 0.8f; 
	public static final RangeDouble VAFNormalRange = new RangeDouble(VAFNormalFrameLower, VAFNormalFrameUpper);
	public static float BinSize             = 0.025f; //smoothing parameter
	
	int   mNumBins;
	int[]   mBinCount;    // The bins in which counts are binned and stored
	float[] mBinValue;    // The value each bin represents
	
	// Statistics values
	float mCountMean = -1;   
	float mVariance = -1;
	float mStdDev = -1;
	
	public AlleleFractionStatsForSample() {
		mNumBins = (int) ((VAFNormalFrameUpper - VAFNormalFrameLower) / BinSize) + 1;
		mBinCount = new   int[mNumBins];
		mBinValue = new float[mNumBins];
		deduceBinValues();
		Arrays.fill(mBinCount, 0);			
	}
	
	public void tabulateAndPerformStatistics(SiteList<ClusteringInputOneSite> sites, SeqPlatform platform) {			
		tallyVariantAlleleFractionsIntoBins(sites, platform, true);
		calculateSummaryStatistics();
	}
	
	public float getValueNStandardDeviationsAway(float n) { return (mCountMean + (n * mStdDev)); }
	
	private void deduceBinValues() {
		for (int i = 0; i < mBinValue.length; i++) {
			mBinValue[i] = VAFNormalFrameLower + ((i + 1) * BinSize);
		}
	}
	
	private void tallyVariantAlleleFractionsIntoBins(SiteList<ClusteringInputOneSite> sites, SeqPlatform platform, boolean clearBinCountBins) {
		if (clearBinCountBins) {
			Arrays.fill(mBinCount, 0);
		}
		
		// First, tally the variant allele frequencies into bins
		for (int row = 0; row < sites.getNumSites(); row++) {
			float vafNormal = sites.getSiteAtIndex(row).calcVAFNormal(); 
			if (RangeDouble.inRangeLowerExclusive(vafNormal, VAFNormalFrameLower, VAFNormalFrameUpper)) {
				int binNumber = (int) ((vafNormal - VAFNormalFrameLower) / BinSize);
				mBinCount[binNumber]++;
			}
		}
	}
			
	private void calculateSummaryStatistics() {
		// Next, calculate summary statistics (mean, variance, standard deviation) 
		float total = 0;
		int numElementsTotal = 0;
					
		for (int i = 0; i < mBinCount.length; i++) {
			total            += (mBinCount[i] * mBinValue[i]);
			numElementsTotal +=  mBinCount[i];
		}	
				
		mCountMean = (float) total / (float) numElementsTotal;
		
		float varianceTotal = 0f;
		for (int i = 0; i < mBinValue.length; i++) { //calculate std. deviation of # points in each bin
			float diff = mCountMean - mBinValue[i];
			float diffSquared = diff * diff;
			float totalForBin = diffSquared * mBinCount[i];
			varianceTotal += totalForBin;			
		}
		
		mVariance = varianceTotal / (float) (numElementsTotal - 1);  // -1 for degrees of freedom		
		mStdDev = (float) Math.sqrt(mVariance);  //standard deviation in NAF-coord across {0.2 < NAF <= 0.8} variants in VAF plot
	}
}