package lohcate.clustering;

import java.util.ArrayList;

import nutils.ArgumentParserUtils;
import nutils.ArgumentParserUtils.InputParameterBoolean;
import nutils.ArgumentParserUtils.InputParameterDouble;
import nutils.ArgumentParserUtils.InputParameterInteger;
import nutils.ArgumentParserUtils.InputParameterNumber;

import com.martiansoftware.jsap.FlaggedOption;
import com.martiansoftware.jsap.JSAP;
import com.martiansoftware.jsap.JSAPResult;
import com.martiansoftware.jsap.StringParser;

public class ClusteringParams {
	
	public static ClusteringParams GlobalClusteringParams = new ClusteringParams();

	// ====================================================================
	// MEMBER VARIABLES
	// ====================================================================
	
	protected InputParameterDouble mAmplificationThreshold   = new InputParameterDouble(2.2 /*2.85*/,  "AmplificationThreshold",   JSAP.NO_SHORTFLAG, "ampThresh", "Copy_Number_Threshold_for_Amplification");
	protected InputParameterDouble mDeletionThreshold        = new InputParameterDouble(1.85 /*1.15*/, "DeletionThreshold",        JSAP.NO_SHORTFLAG, "delThresh", "Copy_Number_Threshold_for_Deletion");
	protected InputParameterDouble mDeletionHomozygousThreshold = new InputParameterDouble(1.1 /*0.5*/, "HomozygousDeletionThreshold",        JSAP.NO_SHORTFLAG, "homDelThresh", "Copy_Number_Threshold_for_Homozygous_Deletion");
	//protected InputParameterDouble mAmplificationThreshold   = new InputParameterDouble(2.5,  "AmplificationThreshold",   JSAP.NO_SHORTFLAG, "ampThresh", "Copy_Number_Threshold_for_Amplification");
	//protected InputParameterDouble mDeletionThreshold        = new InputParameterDouble(1.5,  "DeletionThreshold",        JSAP.NO_SHORTFLAG, "delThresh", "Copy_Number_Threshold_for_Deletion");
	protected InputParameterDouble mGermlineTrisomyThreshold = new InputParameterDouble(1.4, "GermlineTrisomyThreshold", JSAP.NO_SHORTFLAG, "germlineAneuploidyGainThreshold", "Copy_Number_Threshold_for_Germline_Chromosomal_Gain");
	protected InputParameterDouble mFDRAlpha                 = new InputParameterDouble(0.01, "mFDRAlpha",                JSAP.NO_SHORTFLAG, "FDR_Alpha", "FDR_Alpha_Value_for_Allelic_Fraction_Imbalance");
	protected InputParameterInteger mAllelicBiasMinNumSamples = new InputParameterInteger(2, "AllelicBiasMinNumSamples", JSAP.NO_SHORTFLAG, "allelicBiasNumSites", "Minimum_Number_of_Samples_at_Site_for_Allelic_Bias_Correction");
	protected InputParameterInteger mGermlineCopyNumMinNumSamples = new InputParameterInteger(2, "GermlineCopyNumMinNumSamples", JSAP.NO_SHORTFLAG, "excessiveGermlineCopyNumber", "Minimum_Number_of_Samples_at_Site_for_Excessive_Germline_Copy_Number");
	protected InputParameterInteger mNumReadsMinNormal        = new InputParameterInteger(15, "NumReadsMinForNormal", JSAP.NO_SHORTFLAG, "numReadsMinNormal", "Minimum_Number_of_Reads_in_Matched_Normal_to_Include_Site");
	protected InputParameterInteger mNumReadsMinTumor         = new InputParameterInteger( 0, "NumReadsMinForTumor", JSAP.NO_SHORTFLAG, "numReadsMinTumor", "Minimum_Number_of_Reads_in_Tumor_to_Include_Site");

	protected InputParameterBoolean mIgnoreAllelicBias     = new InputParameterBoolean(false, "IgnoreAllelicBias",     JSAP.NO_SHORTFLAG, "ignoreAllelicBias", "Specifies whether allelic bias is to be ignored and not corrected");
	protected InputParameterBoolean mIgnoreMultipleTesting = new InputParameterBoolean(false, "IgnoreMultipleTesting", JSAP.NO_SHORTFLAG, "ignoreMultipleTesting", "Specifies whether multiple testing correction is not to be done");
	protected InputParameterBoolean mIsTumorOnly           = new InputParameterBoolean(false, "IsTumorOnly", JSAP.NO_SHORTFLAG, "isTumorOnly", "Indicates that only tumor data is present and that it should be compared with a simulated normal");
	protected InputParameterBoolean mSimulation            = new InputParameterBoolean(false, "simulation", 'S', "simulation", "Specifies whether to run internal testing by simulation");
	
	
	
	protected InputParameterInteger mNumPasses = new InputParameterInteger(3, "NumPasses", JSAP.NO_SHORTFLAG, "numPasses", "The number of passes to go over a region during segmentation.");
	
	protected InputParameterInteger mSmoothBlockLength = new InputParameterInteger(1_000_000, "SmoothingBlockLength", JSAP.NO_SHORTFLAG, "smoothingBlockLength", "The maximum size of a block (in bases) for which copy number smoothing is to be applied.");
	
	protected double mVAFUpperFrame;
	protected double mVAFLowerFrame; 
		
	protected ArrayList<InputParameterDouble>  mParamsDouble;
	protected ArrayList<InputParameterInteger> mParamsInteger;
	protected ArrayList<InputParameterBoolean> mParamsBoolean;

	public ClusteringParams() {
		mParamsDouble = new ArrayList<InputParameterDouble>();
		mParamsDouble.add(mAmplificationThreshold);
		mParamsDouble.add(mDeletionThreshold);
		mParamsDouble.add(mGermlineTrisomyThreshold);
		mParamsDouble.add(mFDRAlpha);
		
		mParamsInteger = new ArrayList<InputParameterInteger>();
		mParamsInteger.add(mAllelicBiasMinNumSamples);
		mParamsInteger.add(mNumReadsMinNormal);
		mParamsInteger.add(mNumReadsMinTumor);
		
		mParamsBoolean = new ArrayList<InputParameterBoolean>();
		mParamsBoolean.add(mIgnoreAllelicBias);
		mParamsBoolean.add(mIgnoreMultipleTesting);
		mParamsBoolean.add(mIsTumorOnly);
		//mParamsBoolean.add(mSimulation);
	}

	// ========================================================================
	public void setIsSimulation(boolean isSimulation) {
		mSimulation.setValue(isSimulation);
	}
	
	// ========================================================================
	public void setIgnoreAllelicBias(boolean shouldIgnore) {
		mIgnoreAllelicBias.setValue(shouldIgnore);
	}
	
	// ========================================================================
	public boolean isTumorOnly() { return mIsTumorOnly.getValue(); }
	
	// ========================================================================
	public JSAP registerClusteringParameters(JSAP jsap) {
		jsap = (jsap == null) ? new JSAP() : jsap;
		
		for (InputParameterBoolean pb : mParamsBoolean) {
			ArgumentParserUtils.createSwitch(pb, jsap);
		}
		
		for (InputParameterNumber<?> pn : mParamsDouble) {
			registerParameterHelper(jsap, pn, JSAP.DOUBLE_PARSER);
		}
		
		for (InputParameterNumber<?> pn : mParamsInteger) {
			registerParameterHelper(jsap, pn, JSAP.INTEGER_PARSER);
		}
				
		return jsap;
	}

	// ========================================================================
	protected static void registerParameterHelper(JSAP jsap, InputParameterNumber<?> parameterNumber, StringParser parser) {		
		FlaggedOption fo = new FlaggedOption(parameterNumber.getName())
				.setStringParser(parser)
				.setDefault(  parameterNumber.getDefaultValue().toString())				
				.setShortFlag(parameterNumber.getShortFlag())
				.setLongFlag( parameterNumber.getLongFlag())
				.setUsageName(parameterNumber.getUsageName());
		ArgumentParserUtils.registerJSAPParameter(jsap, fo);
	}
	
	// ========================================================================
	public void configureParameters(JSAPResult config) {
		for (InputParameterDouble pd : mParamsDouble) {			
			pd.setValue(config.getDouble(pd.getName()));
		}		
		
		for (InputParameterInteger pi : mParamsInteger) {
			pi.setValue(config.getInt(pi.getName()));
		}
		
		for (InputParameterBoolean pb : mParamsBoolean) {
			pb.setValue(config.getBoolean(pb.getName()));
		}
	}
	
	// ========================================================================
}
	