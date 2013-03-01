package nutils;

public class RangeDouble {

	// ========================================================================
	protected double mBoundLower;
	protected double mBoundUpper;
	
	// ========================================================================
	public RangeDouble(double rangeLower, double rangeUpper) {
		mBoundLower = rangeLower;
		mBoundUpper = rangeUpper;
	}

	// ========================================================================
	public double getBoundLower() { return mBoundLower; }
	public double getBoundUpper() { return mBoundUpper; }

	// ========================================================================
	public void setBoundLower(double boundLower) { this.mBoundLower = boundLower; }
	public void setBoundUpper(double boundUpper) { this.mBoundUpper = boundUpper; }

	// ========================================================================
	public boolean inRange(double value) {
		return ((mBoundLower <= value) && (value <= mBoundUpper));
	}
	
	// ========================================================================
	public boolean inRangeBothExclusive(double value) {
		return ((mBoundLower < value) && (value < mBoundUpper));
	}

	// ========================================================================
	public boolean inRangeLowerExclusive(double value) {
		return ((mBoundLower < value) && (value <= mBoundUpper));
	}
	
	// ========================================================================
	public boolean inRangeUpperExclusive(double value) {
		return ((mBoundLower <= value) && (value < mBoundUpper));
	}

	// ========================================================================
	/** Returns whether the value is in the range, with an exclusive lower bound and inclusive upper bound. */
	public static boolean inRangeLowerExclusive(double value, double boundLower, double boundUpper) {
		return ((boundLower < value) && (value <= boundUpper));
	}

	// ========================================================================
	private static void TestRangeDouble() {
		RangeDouble rd = new RangeDouble(1, 2);
		
		System.out.println(rd.inRange(1));
		System.out.println(rd.inRange(2));

		System.out.println(rd.inRangeBothExclusive(1));
		System.out.println(rd.inRangeBothExclusive(2));	
		
		System.out.println(rd.inRangeLowerExclusive(1));
		System.out.println(rd.inRangeLowerExclusive(2));	

		System.out.println(rd.inRangeUpperExclusive(1));
		System.out.println(rd.inRangeUpperExclusive(2));	

	}
	
	// ========================================================================
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		TestRangeDouble();
	}

}
