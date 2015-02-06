package lohcate.clustering;

import genomeEnums.TissueType;
import genomeUtils.RegionRange;
import genomeUtils.RegionSimulator;
import genomeUtils.SiteList;

import java.util.ArrayList;

import com.google.common.primitives.Primitives;

import nutils.CompareUtils;
import nutils.PrimitiveWrapper;
import nutils.RangeDouble;
import nutils.BitUtils.Compactor.CompactorIntoLong;
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
	public double calcAverageCoverage(TissueType t) {
		long sumReads = 0;
		for (ClusteringInputOneSite oneSiteInfo : mInfoSites) {
			sumReads += (t == TissueType.Normal) ? oneSiteInfo.getCovgTotalNormal() : oneSiteInfo.getCovgTotalTumor();
		}		
		return (double) sumReads / (double) getNumSites();
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
	public long getIndicesForRegion(RegionRange range) {
		int indexStart = getIndex(range.getChromosome(), range.getRangeStart()); 
		CompareUtils.ensureTrue(indexStart >= 0, "ERROR: Starting index must be > 0");
		
		int indexEnd = getIndex(range.getChromosome(), range.getRangeEnd());		
		CompareUtils.ensureTrue(indexEnd >= indexStart, "ERROR: Ending index must be >= starting index!");
		
		long compactUnit = 0;
		compactUnit = CompactorIntoLong.TwoIntsIntoLong.Compactor.setValue(CompactorIntoLong.TwoIntsIntoLong.IntMSB, indexStart, compactUnit);
		compactUnit = CompactorIntoLong.TwoIntsIntoLong.Compactor.setValue(CompactorIntoLong.TwoIntsIntoLong.IntLSB, indexEnd,   compactUnit);
		return compactUnit;
	}
	
	// ========================================================================
	private void parseLines(ArrayList<String> rows) {			
		for (String line : rows) {
			ClusteringInputOneSite cios = new ClusteringInputOneSite(line, SeqPlatform.Illumina);
			mInfoSites.add(cios);
		}
	}	
}