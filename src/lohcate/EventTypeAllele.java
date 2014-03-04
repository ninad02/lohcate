package lohcate;

public enum EventTypeAllele {
	
	LossReference_VariantCommon(true),
	LossReference_VariantRare  (true),
	LossVariant_VariantCommon  (true),
	LossVariant_VariantRare    (true),
	
	GainReference_VariantCommon (false),
	GainReference_VariantRare   (false),
	GainVariant_VariantCommon   (false),
	GainVariant_VariantRare     (false),
	;
	
	private boolean mIsLossType;
	
	private EventTypeAllele(boolean isLossType) {
		mIsLossType = isLossType;
	}
	
	public boolean getIsLossType() { return mIsLossType; }
}
