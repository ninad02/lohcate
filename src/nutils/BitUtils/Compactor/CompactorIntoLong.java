package nutils.BitUtils.Compactor;

import nutils.CompareUtils;


public class CompactorIntoLong<E extends Enum<E> & CompactorInf<E>> extends Compactor<E>  {

	// ========================================================================
	public CompactorIntoLong(Class<E> enumClass, boolean useMostSignificantBit) {
		super(enumClass, useMostSignificantBit);
	}

	// ========================================================================
	@Override
	protected int getBitCapacity() { return Long.SIZE; }

	// ========================================================================
	public long setValue(E variable, long value, long compactUnit) {
		return super.setValue(variable, value, compactUnit);
	}
	
	// ========================================================================
	public long getValue(E variable, long compactUnit) {
		return super.getValue(variable, compactUnit);		
	}
	
	// ========================================================================
	public static enum TwoIntsIntoLong implements CompactorInf<TwoIntsIntoLong> {
		IntMSB(32),
		IntLSB(32);

		public static CompactorIntoLong<TwoIntsIntoLong> Compactor = new CompactorIntoLong<>(TwoIntsIntoLong.class, true);
		
		private int mNumBits;
		private TwoIntsIntoLong(int numBits) {
			mNumBits = numBits;
		}
		
		@Override
		public int getNumBits() { return mNumBits; }		
	}
	
	// ========================================================================
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		TestStress();
	}
	
	// ========================================================================
	private static void TestStress() {
		int numBitsChrom = TestEnum.Chrom.getNumBits();
		int numBitsPosition = TestEnum.Position.getNumBits();
		int numPossibilitiesChrom = (1 << numBitsChrom);
		int numPossibilitiesPosition = (1 << numBitsPosition);
		System.out.println("Num Possibilities Chrom: " + numPossibilitiesChrom);
		System.out.println("Num Possibilities Position: " + numPossibilitiesPosition);
		
		CompactorIntoLong<TestEnum> compactor = new CompactorIntoLong<>(TestEnum.class, true);
		
		for (int iChrom = 0; iChrom < numPossibilitiesChrom; iChrom++) {
			for (int iPos = 0; iPos < numPossibilitiesPosition; iPos++) {
				long compactUnit = -1L;
				compactUnit = compactor.setValue(TestEnum.Chrom, iChrom, compactUnit);
				compactUnit = compactor.setValue(TestEnum.Position, iPos, compactUnit);
				
				long chromNum = compactor.getValue(TestEnum.Chrom, compactUnit);
				long position = compactor.getValue(TestEnum.Position, compactUnit);
				
				if (chromNum != iChrom) {
					CompareUtils.ensureTrue(false, "Chrom doesn't match: " + iChrom + "\t" + chromNum);					
				} else if (position != iPos) {
					CompareUtils.ensureTrue(false, "Position doesn't match: " + iPos + "\t" + position);					
				}
			}
			System.out.println("Chrom " + iChrom + " done!");
		}
	}

	// ========================================================================
	private static enum TestEnum implements CompactorInf<TestEnum> {
		Chrom(5), 
		Position(28),
		Other(31);
		
		private int mNumBits;
		private TestEnum(int numBits) {
			mNumBits = numBits;
		}
		
		public int getNumBits() { return mNumBits; }
	}


}
