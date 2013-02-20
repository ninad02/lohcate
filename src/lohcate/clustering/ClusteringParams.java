package lohcate.clustering;

import genomeEnums.TissueType;

import java.util.ArrayList;

import nutils.ArgumentParserUtils;
import nutils.ArgumentParserUtils.InputParameterBoolean;
import nutils.ArgumentParserUtils.InputParameterDouble;
import nutils.ArgumentParserUtils.InputParameterNumber;

import com.martiansoftware.jsap.FlaggedOption;
import com.martiansoftware.jsap.JSAP;
import com.martiansoftware.jsap.JSAPResult;

public class ClusteringParams {
	
	public static ClusteringParams GlobalClusteringParams = new ClusteringParams();

	// ====================================================================
	// MEMBER VARIABLES
	// ====================================================================
	
	protected InputParameterDouble mAmplificationThreshold   = new InputParameterDouble(2.3,  "AmplificationThreshold",   JSAP.NO_SHORTFLAG, "ampThresh", "Copy_Number_Threshold_for_Amplification");
	protected InputParameterDouble mDeletionThreshold        = new InputParameterDouble(1.7,  "DeletionThreshold",        JSAP.NO_SHORTFLAG, "delThresh", "Copy_Number_Threshold_for_Deletion");
	protected InputParameterDouble mGermlineTrisomyThreshold = new InputParameterDouble(1.4, "GermlineTrisomyThreshold", JSAP.NO_SHORTFLAG, "germlineAneuploidyGainThreshold", "Copy_Number_Threshold_for_Germline_Chromosomal_Gain");
	protected InputParameterDouble mFDRAlpha                 = new InputParameterDouble(0.01, "mFDRAlpha",                JSAP.NO_SHORTFLAG, "FDR_Alpha", "FDR_Alpha_Value_for_Allelic_Fraction_Imbalance");

	protected InputParameterBoolean mIgnoreAllelicBias     = new InputParameterBoolean(false, "IgnoreAllelicBias",     JSAP.NO_SHORTFLAG, "ignoreAllelicBias", "Specifies whether allelic bias is to be ignored and not corrected");
	protected InputParameterBoolean mIgnoreMultipleTesting = new InputParameterBoolean(false, "IgnoreMultipleTesting", JSAP.NO_SHORTFLAG, "ignoreMultipleTesting", "Specifies whether multiple testing correction is not to be done"); 
		
	protected ArrayList<InputParameterDouble>  mParamsDouble;
	protected ArrayList<InputParameterBoolean> mParamsBoolean;

	public ClusteringParams() {
		mParamsDouble = new ArrayList<InputParameterDouble>();
		mParamsDouble.add(mAmplificationThreshold);
		mParamsDouble.add(mDeletionThreshold);
		mParamsDouble.add(mGermlineTrisomyThreshold);
		mParamsDouble.add(mFDRAlpha);
		
		mParamsBoolean = new ArrayList<InputParameterBoolean>();
		mParamsBoolean.add(mIgnoreAllelicBias);
		mParamsBoolean.add(mIgnoreMultipleTesting);
	}

	// ========================================================================
	public void setIgnoreAllelicBias(boolean shouldIgnore) {
		mIgnoreAllelicBias.setValue(shouldIgnore);
	}
	
	// ========================================================================
	public JSAP registerClusteringParameters(JSAP jsap) {
		jsap = (jsap == null) ? new JSAP() : jsap;
		
		//ArgumentParserUtils.createSwitch(mIgnoreAllelicBias, jsap);
		ArgumentParserUtils.createSwitch(mIgnoreMultipleTesting, jsap);
		
		for (InputParameterNumber<?> pn : mParamsDouble) {
			FlaggedOption fo = new FlaggedOption(pn.getName())
					.setStringParser(JSAP.DOUBLE_PARSER)
					.setDefault(pn.getDefaultValue().toString())				
					.setShortFlag(pn.getShortFlag())
					.setLongFlag(pn.getLongFlag())
					.setUsageName(pn.getUsageName());
			ArgumentParserUtils.registerJSAPParameter(jsap, fo);
		}
				
		return jsap;
	}
	
	// ========================================================================
	public void configureParameters(JSAPResult config) {
		for (InputParameterDouble pd : mParamsDouble) {			
			pd.setValue(config.getDouble(pd.getName()));
		}		
	}
	
	// ========================================================================
}
	