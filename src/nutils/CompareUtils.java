package nutils;

public class CompareUtils {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub

	}

	/** Returns the standard comparisons of two integers. */
	public static int compareInt(int num1, int num2) {
		return ((num1 > num2) ? 1 : ((num1 < num2) ? -1 : 0));
	}

	/** Returns the standard comparisons of two longs. */
	public static int compareLong(long num1, long num2) {
		return ((num1 > num2) ? 1 : ((num1 < num2) ? -1 : 0));
	}
	
	/** Returns whether the double value is within a range. */
	public static boolean inRange(double val, double target, double leeway) {
		return ((val >= target - leeway) && (val <= target + leeway));
	}

	/** Returns the number of differences between the two strings. */
	public static int hammingDistance(String s1, String s2) {
		if (s1.length() != s2.length()) {
			// TODO - handle this case better
			return Integer.MAX_VALUE;
		} else {
			int diffCount = 0;
			int sLength = s1.length();
			for (int i = 0; i < sLength; i++) {
				if (s1.charAt(i) != s2.charAt(i)) {
					diffCount++;
				}
			}
			return diffCount;
		}
	}

	/** Compares only the alpha numeric digits. */
	public static boolean approximateMatch(String s0, String s1, boolean substringMatch) {
		StringBuilder s0Buffer = new StringBuilder(s0.length());
		StringBuilder s1Buffer = new StringBuilder(s1.length());
		CompareUtils.approximateMatch_fillStringBuilder(s0Buffer, s0);
		CompareUtils.approximateMatch_fillStringBuilder(s1Buffer, s1);
		
		String s0_Simp = s0Buffer.toString();
		String s1_Simp = s1Buffer.toString();
		
		return (s0_Simp.indexOf(s1_Simp) >= 0);
	}

	public static void approximateMatch_fillStringBuilder(StringBuilder sBuilder, String s) {
		for (int i = 0; i < s.length(); i++) {
			char ch = s.charAt(i);
			if (Character.isLetterOrDigit(ch)) {
				sBuilder.append(Character.toUpperCase(ch));
			}
		}
	}

	/** Returns whether a number is even.  Is done by a binary AND operation. */
	public static boolean isEven(int num) { return ((num & 1) == 0); }

	/** Returns whether an object is null. */
	public static boolean isNull(Object o) { return (o == null); }

	/** Returns whether an object is not null. */
	public static boolean isNotNull(Object o) { return (o != null); }
	
	/** Given a range minimum and range maximum, this ensures that the input value is within
	 *  the range, inclusive.  In other words, if the input value is greater than the maximum,
	 *  the maximum is returned.  If it is smaller than the minimum, the minimum is returned.
	 *  Otherwise, the number itself is returned.
	 */
	public static long ensureInRange(long value, long minValue, long maxValue) {
		ensureTrue(maxValue >= minValue, "ERROR: CompareUtils.ensureInRange(): maxValue is less than minValue!");
		return Math.max(minValue, Math.min(value, maxValue));
	}

	public static void throwErrorAndExit(String errorString) {
		Exception e = new Exception(errorString);
		e.printStackTrace();
		System.exit(-1);			
	}

	/** Given a boolean expression and error message, this quits the program if the
	 *  boolean expression is false (after which, the error message is printed. */
	public static void ensureTrue(boolean condition, String errorString) {
		if (!condition) throwErrorAndExit(errorString);
	}
	
	/** Given an object, this quits the program if the object is null (after which,
	 *  the error message is printed.
	 */
	public static void ensureNotNull(Object o, String errorString) {
		ensureTrue(isNotNull(o), errorString);
	}

}

