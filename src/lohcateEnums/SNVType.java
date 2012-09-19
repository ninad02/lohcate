package lohcateEnums;

public enum SNVType {
	NonSynonymous_SNV    ("nonsynonymous_SNV"),
	Synonymous_SNV       ("synonymous_SNV"),
	
	downstream,
	exonicOrfErr         ("exonic-orf-err"),
	intergenic,
	intronic,
	ncRNA_exonic,
	ncRNA_intronic,
	ncRNA_splicing,
	ncRNA_UTR3,
	ncRNA_UTR5,
	splicing,
	stopgain_SNV,
	stoploss_SNV,	
	upstream,
	upstream_downstream  ("upstream;downstream"),
	UTR3,
	UTR5
	
	;
	
	private String mLowerCase;
	private String mUpperCase;
	private String mPrintName;
	
	private SNVType() {
		this(null);
	}
	
	private SNVType(String printName) {
		mLowerCase = name().toLowerCase();
		mUpperCase = name().toUpperCase();
		mPrintName = (printName == null) ? name() : printName;		
	}
	
	public String toUpperCase() { return mUpperCase; }
	public String toLowerCase() { return mLowerCase; }
	
	public static SNVType getSNVType(String snvType) {
		for (SNVType enumType : values()) {
			
			if ( snvType.equalsIgnoreCase(enumType.name()) ||
				 snvType.equalsIgnoreCase(enumType.mPrintName)) {
				 
				return enumType;
			}
		}
		return null;
	}
	
	/*
	 * 

	 */
}
