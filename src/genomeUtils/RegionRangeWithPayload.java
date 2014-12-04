/**
 * 
 */
package genomeUtils;

import nutils.Cast;
import genomeEnums.Chrom;

/**
 * @author Ninad Dewal
 *
 */
public class RegionRangeWithPayload<T extends Comparable<T>> 
	extends RegionRange 
	implements Comparable<RegionRangeWithPayload<T>> {

	protected T mPayload;
	
	// ========================================================================
	public RegionRangeWithPayload() {		
		super();	
		setPayload(null);
	}

	// ========================================================================
	/**
	 * @param chrom
	 * @param rangeStart
	 */
	public RegionRangeWithPayload(Chrom chrom, int rangeStart, T payload) {
		super(chrom, rangeStart);
		setPayload(payload);
	}

	// ========================================================================
	/**
	 * @param chrom
	 * @param rangeStart
	 * @param rangeEnd
	 */
	public RegionRangeWithPayload(Chrom chrom, int rangeStart, int rangeEnd, T payload) {
		super(chrom, rangeStart, rangeEnd);
		setPayload(payload);
	}

	// ========================================================================
	public RegionRangeWithPayload(RegionRangeWithPayload<T> rhs) {
		super(Cast.upcast(rhs, RegionRange.class));
		setPayload(rhs.mPayload);
	}

	// ========================================================================
	public void setPayload(T payload) { mPayload = payload; }
	
	// ========================================================================
	public T getPayload() { return mPayload; }

	// ========================================================================
	//@Override
	public int compareTo(RegionRangeWithPayload<T> rhs) {
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
