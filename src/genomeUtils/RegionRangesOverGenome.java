package genomeUtils;

import java.util.ArrayList;
import java.util.ListIterator;

import genomeEnums.Chrom;
import nutils.EnumMapSafe;
import nutils.UtilsBasic;

/** A class to house ranges for a particular sample, grouped by chromosome.  Allows for also
 *  holding metainformation about a particular sample
 * @author Ninad Dewal
 *
 * @param <E> The range class specified
 * @param <T> An object that holds meta information for the sample
 */
public class RegionRangesOverGenome<E extends RegionRange, T extends Cloneable> {

	// ========================
	// MEMBER VARIABLES
	// ========================
	EnumMapSafe<Chrom, ArrayList<E>> mRegionsByChrom;
	T mMetaInfo;
	
	// ========================

	// ========================================================================
	public RegionRangesOverGenome(T metaInfo) {
		mRegionsByChrom = EnumMapSafe.createEnumMapOfArrayLists(Chrom.class);
		mMetaInfo = metaInfo;		
	}
			
	// ========================================================================
	public RegionRangesOverGenome() {
		this((T) null);
	}
	
	// ========================================================================
	public RegionRangesOverGenome(RegionRangesOverGenome<E, T> rhs) {
		this(UtilsBasic.getClone(rhs.mMetaInfo));
		
		// Now deep copy the regions
		for (Chrom chrom : Chrom.values()) {
			ArrayList<E> listRhs = rhs.mRegionsByChrom.get(chrom);
			ArrayList<E> list    =     mRegionsByChrom.get(chrom); 
			
			for (E elementRhs : listRhs) {				
				list.add((E) elementRhs.getCopy());					
			}
		} 		
	}
	
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
	



	
	
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub
	}

}
