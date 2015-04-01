/**
 * 
 */
package genomeUtils;

import nutils.CloneInf;
import genomeEnums.Chrom;

/**
 * @author Ninad Dewal
 *
 */
public class RegionRangeWithPayload<E extends CloneInf<E> & Comparable<E>>
	extends RegionRange<RegionRangeWithPayload<E>> {

	protected E mPayload;
	
	// ========================================================================
	public RegionRangeWithPayload() {		
		super();	
		setPayload(null);
	}

	// ========================================================================
	public RegionRangeWithPayload(Chrom chrom, int rangeStart, E payload) {
		super(chrom, rangeStart);
		setPayload(payload);
	}

	// ========================================================================
	public RegionRangeWithPayload(Chrom chrom, int rangeStart, int rangeEnd, E payload) {
		super(chrom, rangeStart, rangeEnd);
		setPayload(payload);
	}

	// ========================================================================
	public RegionRangeWithPayload(RegionRangeWithPayload<E> rhs, boolean deepCopy) {
		super(rhs);		
		setPayload(deepCopy ? rhs.mPayload.makeClone() : rhs.mPayload);
	}

	// ========================================================================
	@Override
	public RegionRangeWithPayload<E> makeClone() {
		return makeClone(true);
	}
	
	// ========================================================================
	@Override
	public RegionRangeWithPayload<E> makeClone(boolean deepCopy) {
		return new RegionRangeWithPayload<E>(this, deepCopy);
	}
	
	// ========================================================================
	public void setPayload(E payload) { mPayload = payload; }
	
	// ========================================================================
	public E getPayload() { return mPayload; }

	// ========================================================================
	@Override
	public int compareTo(RegionRangeWithPayload<E> rhs) {		
		int result = super.compareTo(rhs);
		if (result == 0) {
			result = mPayload.compareTo(rhs.mPayload);
		}
		return result;
	}
	
	// ========================================================================
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub

	}


}
