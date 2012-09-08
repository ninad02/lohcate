package lohcateEnums;

public enum SNVType {
	NonSynonymous,
	Synonymous
	;
	
	private String mLowerCase;
	private String mUpperCase;
	
	private SNVType() {
		mLowerCase = name().toLowerCase();
		mUpperCase = name().toUpperCase(); 
	}
	
	public String toUpperCase() { return mUpperCase; }
	public String toLowerCase() { return mLowerCase; }
	
	public static SNVType getSNVType(String snvType) {
		for (SNVType enumType : values()) {
			if ( snvType.equals(enumType.name()) ||
				 snvType.toLowerCase().equals(enumType.mLowerCase) ||
				 snvType.toUpperCase().equals(enumType.mUpperCase)) {
				return enumType;
			}
		}
		return null;
	}
}
