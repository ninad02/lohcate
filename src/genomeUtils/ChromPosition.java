package genomeUtils;

import genomeEnums.Chrom;
import nutils.Cast;
import nutils.BitUtils.Compactor.CompactorInf;
import nutils.BitUtils.Compactor.CompactorIntoLong;

/** 
 * 
 * @author Ninad Dewal
 * A simple non-optimized class to represent a chromosome and position.    
 */

public class ChromPosition implements Comparable<ChromPosition >{

	// ========================================================================
	// MEMBER VARIABLES
	// ========================================================================
	private long mCompactUnit;
	
	// ========================================================================
	public ChromPosition() {
		this(Chrom.c0, 0);		
	}
	
	public ChromPosition(Chrom chrom, int position) {
		setChromosome(chrom);
		setPosition(position);
	}
	
	// ========================================================================
	public void setChromosome(Chrom chrom) {
		mCompactUnit = ChromPos.Compactor.setValue(ChromPos.Chrom, chrom.ordinal(), mCompactUnit);
	}
	
	// ========================================================================
	public void setPosition(int position) {
		mCompactUnit = ChromPos.Compactor.setValue(ChromPos.Position, position, mCompactUnit);
	}
	
	// ========================================================================
	public Chrom getChromosome() {			
		return Chrom.getChrom(Cast.toByte( ChromPos.Compactor.getValue(ChromPos.Chrom, mCompactUnit) ));
	}
	
	// ========================================================================
	public int getPosition() {
		return Cast.toInt( ChromPos.Compactor.getValue(ChromPos.Position, mCompactUnit) );
	}
	
	// ========================================================================
	private static enum ChromPos implements CompactorInf<ChromPos> {
		Chrom(5),
		Position(29);
	
		private static final CompactorIntoLong<ChromPos> Compactor = new CompactorIntoLong<>(ChromPos.class, false);
		
		private int mNumBits;
		private ChromPos(int size) {
			mNumBits = size;
		}
		
		@Override
		public int getNumBits() { return mNumBits; }
	}
	
	// ========================================================================
	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + (int) (mCompactUnit ^ (mCompactUnit >>> 32));
		return result;
	}

	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		ChromPosition other = (ChromPosition) obj;
		if (mCompactUnit != other.mCompactUnit)
			return false;
		return true;
	}

	// ========================================================================
	@Override
	public int compareTo(ChromPosition o) { return Long.compare(mCompactUnit, o.mCompactUnit); }
	
	// ========================================================================	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub

	}



}
