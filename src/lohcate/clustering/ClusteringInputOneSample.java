package lohcate.clustering;

import genomeUtils.RegionSimulator;
import genomeUtils.SiteList;

import java.util.ArrayList;

import nutils.RangeDouble;
import lohcateEnums.SeqPlatform;

/**
 * 
 * @author Ninad Dewal
 *
 */
public class ClusteringInputOneSample extends SiteList<ClusteringInputOneSite> implements RegionSimulator.SampleInformation<ClusteringInputOneSite> {
	String mSampleNameRoot;

	// ========================================================================
	public ClusteringInputOneSample(int numSitesEstimated) {
		super(numSitesEstimated, new ClusteringInputOneSite(), ClusteringInputOneSite.ClusteringInputOneSiteComparator);
	}
	
	// ========================================================================
	public ClusteringInputOneSample(ArrayList<String> rows) {	
		this(rows, "");
	}
	
	// ========================================================================
	public ClusteringInputOneSample(ArrayList<String> rows, String sampleNameRoot) {
		super(rows.size(), new ClusteringInputOneSite(), ClusteringInputOneSite.ClusteringInputOneSiteComparator);		
		mSampleNameRoot = sampleNameRoot;
		parseLines(rows);
	}
	
	// ========================================================================
	public String getSampleNameRoot() { return mSampleNameRoot; }
	
	// ========================================================================		
	public void clear() {
		for (ClusteringInputOneSite oneSiteInfo : mInfoSites) {
			oneSiteInfo.clear();
		}
		
		super.clear();
		
		mSampleNameRoot = "";
	}
	
	// ========================================================================
	/** Indices inclusive. */
	public int getNumHetSitesInRangeByIndex(int indexStart, int indexEnd, RangeDouble hetVAFRange) {
		int numHetSites = 0;
		for (int index = indexStart; index <= indexEnd; index++) {			
			if (hetVAFRange.inRangeLowerExclusive(mInfoSites.get(index).calcVAFNormal())) {
				++numHetSites;
			}			
		}
		return numHetSites;
	}
	
	// ========================================================================
	private void parseLines(ArrayList<String> rows) {			
		for (String line : rows) {
			ClusteringInputOneSite cios = new ClusteringInputOneSite(line, SeqPlatform.Illumina);
			mInfoSites.add(cios);
		}
	}	
}