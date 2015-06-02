package genomeUtils;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.ListIterator;

import genomeEnums.Chrom;
import nutils.CloneInf;
import nutils.EnumMapSafe;
import nutils.IOUtils;

/** A class to house ranges for a particular sample, grouped by chromosome.  Allows for also
 *  holding metainformation about a particular sample
 * @author Ninad Dewal
 *
 * @param <E> The range class specified
 * @param <T> An object that holds meta information for the sample
 */
public class RegionRangesOverGenome<E extends RegionRange<E>, T extends Serializable> 
	implements CloneInf<RegionRangesOverGenome<E, T>>, Serializable {

	/**
	 * 
	 */
	private static final long serialVersionUID = -6561264663430750082L;
	
	// ========================
	// MEMBER VARIABLES
	// ========================
	protected EnumMapSafe<Chrom, ArrayList<E>> mRegionsByChrom;
	protected T mMetaInfo;
	
	// ========================

	// ========================================================================
	public RegionRangesOverGenome(T metaInfo) {
		mRegionsByChrom = EnumMapSafe.createEnumMapOfArrayLists(Chrom.class);
		setMetaInfo(metaInfo);
	}
			
	// ========================================================================
	public RegionRangesOverGenome() {
		this((T) null);
	}
	
	// ========================================================================
	public RegionRangesOverGenome(RegionRangesOverGenome<E, T> rhs) {
		this(rhs, true);
	}
	
	// ========================================================================
	public RegionRangesOverGenome(RegionRangesOverGenome<E, T> rhs, boolean deepCopy) {
		this((T) null);
		performCopy(rhs, deepCopy);
	}
	
	// ========================================================================
	private void performCopy(RegionRangesOverGenome<E, T> rhs, boolean deepCopy) {
				
		setMetaInfo(deepCopy ? IOUtils.makeCopy(rhs.mMetaInfo) : rhs.mMetaInfo);
		
		// Now iterate over the chromosomes
		for (Chrom chrom : Chrom.values()) {
			ArrayList<E> listRhs = rhs.mRegionsByChrom.get(chrom);
			ArrayList<E> list    =     mRegionsByChrom.get(chrom); 
			list.clear();  // Clear the list out just to make sure there are no elements
			
			if (deepCopy) {
				for (E elementRhs : listRhs) {				
					list.add(elementRhs.makeClone());					
				}
			} else {
				list.addAll(listRhs);
			}
		} 
	}

	// ========================================================================
	@Override
	public RegionRangesOverGenome<E, T> makeClone() {
		return makeClone(true);
	}
	
	// ========================================================================
	@Override
	public RegionRangesOverGenome<E, T> makeClone(boolean deepCopy) {
		return new RegionRangesOverGenome<E, T>(this, deepCopy);
	}

	
	// ========================================================================
	public void setMetaInfo(T metaInfo) { this.mMetaInfo = metaInfo; }
	
	// ========================================================================
	public T getMetaInfo() { return this.mMetaInfo; }
	
	// ========================================================================
	public void addRegion(Chrom chrom, E region) {
		mRegionsByChrom.get(chrom).add(region);
	}
		
	// ========================================================================
	public void addRegion(E region) {
		addRegion(region.getChromosome(), region);
	}

	// ========================================================================
	/** Returns the list of regions for a chromosome. */
	public ArrayList<E> getRegions(Chrom chrom) { return mRegionsByChrom.get(chrom); }
	
	// ========================================================================
	/** Given a chromosome and position, returns the region that includes the coordinate, or 
	 *  null if no region includes the coordinate.
	 */
	public E getRegion(Chrom chrom, int position) {
		int resultIndex = getIndexOfRegion(chrom, position);
		return ((resultIndex < 0) ? null : getRegions(chrom).get(resultIndex));
	}
	
	// ========================================================================
	/** Given a chromosome and position, returns the index of the region that includes the coordinate, or 
	 *  -1 if no region includes the coordinate.
	 */
	public int getIndexOfRegion(Chrom chrom, int position) {
		ArrayList<E> regions = getRegions(chrom);
		for (int i = 0; i < regions.size(); i++) {			
			if (regions.get(i).inRange(chrom, position)) {
				return i; 
			} 
		}
		
		return -1 ;
	}
	
	
	// ========================================================================
	/** Returns an iterator for the regions for a chromosome. */
	public ListIterator<E> getIteratorForChromosome(Chrom chrom) {
		return getRegions(chrom).listIterator();
	}

	// ========================================================================	
	public void iterateOverChromosomes(ActionerPerChromosome<E> actioner) {
		for (Chrom chrom : Chrom.ValidChromosomes) {
			actioner.takeAction(chrom, getRegions(chrom));			
		}
	}
	
	// ========================================================================	
	public static interface ActionerPerChromosome<D extends RegionRange<D>> {
		public void takeAction(Chrom chrom, ArrayList<D> regions);
	}
	
	// ========================================================================
	/** Intersects this object with the argument and saves the results in this object. */
	public RegionRangesOverGenome<E, T> intersect(RegionRangesOverGenome<E, T> rhs, RegionBreakerAndIntersecter.RegionIntersectTester<E> intersecterObject) {
		
		// If null then return
		if (rhs == null) return this;
		
		// Now, iterate over the chromosomes
		for (Chrom chrom : Chrom.ValidChromosomes) {			
			RegionBreakerAndIntersecter.takeUnionAndBreakDownIntersectingRegions(this.getRegions(chrom), rhs.getRegions(chrom), intersecterObject);
		}
		
		return this;
	}

	// ========================================================================
	public static<E extends RegionRange<E>, T extends Serializable> RegionRangesOverGenome<E, T> 
		intersectAll(Collection<RegionRangesOverGenome<E, T>> allSamples, 
					 RegionBreakerAndIntersecter.RegionIntersectTester<E> intersecterObject) {

		// Set a return object.  
		RegionRangesOverGenome<E, T> baseSample = null;
		
		for (Iterator<RegionRangesOverGenome<E, T>> sampleIter = allSamples.iterator(); sampleIter.hasNext(); ) {
			RegionRangesOverGenome<E, T> currentSample = sampleIter.next();
			
			if (baseSample == null) {
				baseSample = currentSample.makeClone();
			} else {
				baseSample.intersect(currentSample, intersecterObject);
			}
		}
		
		// Will return null if input list is empty
		return baseSample;
	}
	
	// ========================================================================
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub
	}



}
