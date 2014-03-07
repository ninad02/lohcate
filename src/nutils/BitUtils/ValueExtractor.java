package nutils.BitUtils;

// ========================================================================
public interface ValueExtractor {

	// ========================================================================
	public static final ValueExtractor LongExtractorWhole = new ValueExtractor() {
		public long extractValue(long compactUnit) { return compactUnit; }
	};
	
	// ========================================================================
	public static final ValueExtractor IntExtractorLSB = new ValueExtractor() { 
		public long extractValue(long compactUnit) { return (compactUnit & 0xFFFFFFFFL); }
	};

	// ========================================================================
	public static final ValueExtractor IntExtractorMSB = new ValueExtractor() { 
		public long extractValue(long compactUnit) { return ((compactUnit >>> Integer.SIZE) & 0xFFFFFFFFL); }
	};

	// ========================================================================
	public long extractValue(long compactUnit);
	
	
}