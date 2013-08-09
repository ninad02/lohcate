package genomeUtils;

import genomeEnums.Chrom;

import java.util.Comparator;

import nutils.ArrayUtils;
import nutils.collectionsSorted.ArrayListSortedComparable;

public abstract class SiteList<E extends Comparable<E> & SiteInformation.Writeable> {
	
	// ========================================================================
	protected ArrayListSortedComparable<E> mInfoSites;
	protected E mDummySite;

	// ========================================================================
	public SiteList(int numSitesEstimated, E dummySiteInstance, Comparator<E> comparator) {
		mInfoSites = new ArrayListSortedComparable<E>(numSitesEstimated, comparator);
		mDummySite = dummySiteInstance;
	}
	
	// ========================================================================
	public synchronized void clear() {
		mInfoSites.clear();		
	}
		
	// ========================================================================
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
			
			E oneSiteInfo = mInfoSites.get(insertPoint);
			return (oneSiteInfo.getChrom().equals(chrom) ? insertPoint : resultIndex);
		}
	}
	
	// ========================================================================
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

	// ========================================================================
	public synchronized E getSiteAtPosition(Chrom chrom, int position) {
		int theIndex = getIndex(chrom, position);
		return (theIndex < 0) ? null : getSiteAtIndex(chrom, theIndex);
	}
	
	// ========================================================================
	public synchronized int getIndex(Chrom chrom, int position) {
		mDummySite.set(chrom, position);
		return mInfoSites.getIndex(mDummySite);
	}
	
	// ========================================================================
	public synchronized E getSiteAtIndex(Chrom chrom, int index) {
		return getSiteAtIndex(index);
	}
	
	// ========================================================================
	public synchronized E getSiteAtIndex(int index) {
		return mInfoSites.get(index);
	}
	
	// ========================================================================
	public synchronized int getNumSites() { return mInfoSites.size(); }

	// ========================================================================
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub

	}

}
