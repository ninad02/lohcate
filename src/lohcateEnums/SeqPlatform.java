package lohcateEnums;

public enum SeqPlatform {
	Illumina,
	SOLiD;
	
	public int getNumValid() { return values().length; }
	
	public static SeqPlatform getPlatform(int option) {		
		for (SeqPlatform sp : values()) {
			if (sp.ordinal() == option) {
				return sp;
			}
		}
		return null;
	}
}
