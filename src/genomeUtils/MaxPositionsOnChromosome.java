package genomeUtils;

import java.util.Arrays;

import genomeEnums.Chrom;

public class MaxPositionsOnChromosome {

	// ========================================================================
	long[] mPositions;
	
	// ========================================================================
	public MaxPositionsOnChromosome() {
		mPositions = new long[Chrom.values().length];
		Arrays.fill(mPositions, 0);
	}

	// ========================================================================
	public void registerPosition(Chrom chrom, long position) {
		int chromIndex = chrom.ordinal();
		mPositions[chromIndex] = Math.max(position, mPositions[chromIndex]);
	}
	
	// ========================================================================
	public long getMaxPositionForChrom(Chrom chrom) { return mPositions[chrom.ordinal()]; }

	// ========================================================================
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub

	}

}
