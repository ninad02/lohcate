package nutils.BitUtils.Compactor;

import nutils.CompareUtils;

public class CompactorIntoInt<E extends Enum<E> & CompactorInf<E>> extends Compactor<E> {

	// ========================================================================
	public CompactorIntoInt(Class<E> enumClass, boolean useMostSignificantBit) {
		super(enumClass, useMostSignificantBit);
	}

	// ========================================================================
	@Override
	protected int getBitCapacity() { return Integer.SIZE; }

	// ========================================================================	
	public int setValue(E variable, int value, int compactUnit) {
		return (int) super.setValue(variable, value, compactUnit);
	}
	
	// ========================================================================
	public int getValue(E variable, int compactUnit) {
		return (int) super.getValue(variable, compactUnit);
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

		CompactorIntoInt<TestEnum> compactor = new CompactorIntoInt<TestEnum>(TestEnum.class, true);

		for (int iChrom = 0; iChrom < numPossibilitiesChrom; iChrom++) {
			for (int iPos = 0; iPos < numPossibilitiesPosition; iPos++) {
				int compactUnit = 0;
				compactUnit = compactor.setValue(TestEnum.Chrom, iChrom, compactUnit);
				compactUnit = compactor.setValue(TestEnum.Position, iPos, compactUnit);
						
				int chromNum = compactor.getValue(TestEnum.Chrom, compactUnit);
				int position = compactor.getValue(TestEnum.Position, compactUnit);

				if (chromNum != iChrom) {
					CompareUtils.ensureTrue(false, "Chrom doesn't match: " + iChrom + "\t" + chromNum);					
				} else if (position != iPos) {
					CompareUtils.ensureTrue(false, "Chr/Position doesn't match: " + iChrom + "/" + iPos + "\t" + chromNum + "/" + position);					
				}
			}
			System.out.println("Chrom " + iChrom + " done!");
		}
	}

	// ========================================================================
	private static enum TestEnum implements CompactorInf<TestEnum> {
		Chrom(5), 
		Position(27);

		private int mNumBits;
		private TestEnum(int numBits) {
			mNumBits = numBits;
		}

		public int getNumBits() { return mNumBits; }
	}
}
