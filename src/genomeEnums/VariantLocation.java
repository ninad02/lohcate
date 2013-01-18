package genomeEnums;

public enum VariantLocation {
	Germline,
	Somatic;
	
	private String mLowerCase;
	private String mUpperCase;
	
	private VariantLocation() {
		mLowerCase = name().toLowerCase();
		mUpperCase = name().toUpperCase(); 
	}
	
	public String toUpperCase() { return mUpperCase; }
	public String toLowerCase() { return mLowerCase; }
	
	public static VariantLocation getVariantLocation(String variantLocation) {
		for (VariantLocation enumType : values()) {
			if (variantLocation.equals(enumType.name()) ||
				variantLocation.toLowerCase().equals(enumType.mLowerCase) ||
				variantLocation.toUpperCase().equals(enumType.mUpperCase)) {
				return enumType;
			}
		}
		return null;
	}
}
