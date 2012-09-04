package lohcateEnums;

public enum SeqPlatform {
	Illumina,
	SOLiD;
	
	public int getNumValid() { return values().length; }
	
	public static SeqPlatform getPlatform(int code) {
		return (inBounds(code)) ? values()[code] : null;
	}
	
	private static boolean inBounds(int code) {
		return (code >= 0 && code < values().length);
	}
}
