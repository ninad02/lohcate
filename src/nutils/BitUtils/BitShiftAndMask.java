package nutils.BitUtils;

// ========================================================================
// ========================================================================
public class BitShiftAndMask implements ValueExtractor {
	public int  mNumBits;
	public long mNumBitsToShift;
	public long mMask;
	public long mMaskShifted;
	public long mClearShifted;
	
	// ========================================================================
	public BitShiftAndMask(int numBits, long numBitsToShift) {
		super();
		mNumBits = numBits;
		mNumBitsToShift = numBitsToShift;
		mMask = BitSetUtils.getMask(numBits); 
        mMaskShifted = mMask << mNumBitsToShift;
        mClearShifted = ~mMaskShifted;
	}	
	
	// ========================================================================
	public long setValueInCompactUnit(long value, long compactUnit) {
		value = Math.max(value, 0);
		value = Math.min(value, mMask);
		return ((compactUnit & mClearShifted) | (value << mNumBitsToShift));			
	}
	
	// ========================================================================
	public long extractValue(long compactUnit) {
		return ((compactUnit >>> mNumBitsToShift) & mMask);
	}
	
	// ========================================================================
	public static BitShiftAndMask createBitShiftAndMaskInChain(int numBits, BitShiftAndMask bsm) {
		return new BitShiftAndMask(numBits, bsm.mNumBitsToShift - numBits);
	}
	
	// ========================================================================
	public void toString(StringBuilder sb, boolean clearStringBuilder, String delimiter) {
		if (clearStringBuilder) {
			sb.setLength(0);
		}
		
		sb.append(mNumBits).append(delimiter)
		  .append(mNumBitsToShift).append(delimiter)
		  .append(Long.toBinaryString(mMask)).append(delimiter)
		  .append(Long.toBinaryString(mMaskShifted)).append(delimiter)
		  .append(Long.toBinaryString(mClearShifted));
	}
}
// ========================================================================
// ========================================================================