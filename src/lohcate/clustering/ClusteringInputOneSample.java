package lohcate.clustering;

import genomeEnums.Chrom;
import genomeUtils.RegionSimulator;

import java.util.ArrayList;
import java.util.Collections;

import lohcateEnums.SeqPlatform;
import nutils.ArrayUtils;

// ========================================================================
// ========================================================================
public class ClusteringInputOneSample implements RegionSimulator.SampleInformation<ClusteringInputOneSite> {
	ArrayList<ClusteringInputOneSite> mInfoSites;
	ClusteringInputOneSite mDummySite;
	String mSampleNameRoot;
	
	public ClusteringInputOneSample(int numSitesEstimated) {
		mInfoSites = new ArrayList<ClusteringInputOneSite>(numSitesEstimated);
		constructorCommon();
	}
	
	public ClusteringInputOneSample(ArrayList<String> rows) {
		this(rows, "");
	}
	
	public ClusteringInputOneSample(ArrayList<String> rows, String sampleNameRoot) {
		mSampleNameRoot = sampleNameRoot;
		mInfoSites = new ArrayList<ClusteringInputOneSite>(rows.size());
		constructorCommon();
		parseLines(rows);
	}
	
	public String getSampleNameRoot() { return mSampleNameRoot; }
	
	public void clear() {
		for (ClusteringInputOneSite oneSiteInfo : mInfoSites) {
			oneSiteInfo.clear();
		}
		
		mSampleNameRoot = "";
	}
	
	private void constructorCommon() {
		mDummySite = new ClusteringInputOneSite();
	}
	
	private void parseLines(ArrayList<String> rows) {			
		for (String line : rows) {
			ClusteringInputOneSite cios = new ClusteringInputOneSite(line, SeqPlatform.Illumina);
			mInfoSites.add(cios);
		}
	}
	
	/** Returns the index of the first site (on the chromosome) in this sample.  If
	 *  the chromosome does not exist, a negative index is returned in accordance to
	 *  the definition from Collections.binarySearch()
	 * @param chrom
	 * @return
	 */
	public synchronized int getIndexChromStart(Chrom chrom) {
		int resultIndex = getIndex(chrom, 1);
		if (resultIndex >= 0) {
			return resultIndex;
		} else {				
			int insertPoint = ArrayUtils.getInsertPoint(resultIndex);
			if (insertPoint >= mInfoSites.size()) return resultIndex;
			
			ClusteringInputOneSite oneSiteInfo = mInfoSites.get(insertPoint);
			return (oneSiteInfo.getChrom().equals(chrom) ? insertPoint : resultIndex);
		}
	}
	
	public synchronized int getIndexChromEnd(Chrom chrom) {
		Chrom nextChrom = chrom.getNextChrom();
		if (nextChrom == null) {
			// There is no next chromosome
			return mInfoSites.size() - 1;
		} else {
			// There is a position for a next chromosome
			int resultIndex = getIndexChromStart(nextChrom);
			if (resultIndex < 0) {
				// This means the next chromosome doesn't exist in the dataset.  
				// Try the following chromosome via recursion
				return getIndexChromEnd(nextChrom);  
			} else {
				return resultIndex - 1;
			}
		}
	}
	
	public synchronized int getIndex(Chrom chrom, int position) {
		mDummySite.setChrom(chrom);
		mDummySite.setPosition(position);
		return Collections.binarySearch(mInfoSites, mDummySite, ClusteringInputOneSite.ClusteringInputOneSiteComparator);			
	}
	
	public synchronized ClusteringInputOneSite getSiteAtIndex(Chrom chrom, int index) {
		return getSiteAtIndex(index);
	}
	
	public synchronized ClusteringInputOneSite getSiteAtIndex(int index) {
		return mInfoSites.get(index);
	}
	
	public synchronized int getNumSites() { return mInfoSites.size(); }
}