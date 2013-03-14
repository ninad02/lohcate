package lohcate.clustering;

import genomeEnums.Chrom;
import genomeEnums.Nuc;
import genomeUtils.GenotypeUtils;
import genomeUtils.RegionSimulator;

import java.util.Comparator;

import lohcate.Script;
import lohcateEnums.MutationType;
import lohcateEnums.SeqPlatform;
import nutils.BitSetUtils;
import nutils.CompareUtils;
import nutils.NumberUtils;
import nutils.StringUtils;
import nutils.BitSetUtils.BitShiftAndMask;
import shared.Utils;

public class ClusteringInputOneSite implements Comparable<ClusteringInputOneSite>, RegionSimulator.SiteInformation {

	// ========================================================================
	public static Comparator<ClusteringInputOneSite> ClusteringInputOneSiteComparator = new Comparator<ClusteringInputOneSite>() {
		public int compare(ClusteringInputOneSite site1, ClusteringInputOneSite site2) {
			//System.out.println(site1.getChrom() + "\t" + site1.getPosition() + "\t" + site2.getChrom() + "\t" + site2.getPosition());
			int result = site1.getChrom().compareTo(site2.getChrom()); 
			if (result == 0) {
				result = Integer.compare(site1.getPosition(), site2.getPosition());
			}
			return result;
		}
	};
	
	// ========================================================================
	// MEMBER VARIABLES
	// ========================================================================
	private long mDataUnit_ChromProsRevVarAllelesMutType;
	
	private static final BitSetUtils.BitShiftAndMask bsmChrom = new BitShiftAndMask(5, 58);
	private static final BitSetUtils.BitShiftAndMask bsmPos   = BitShiftAndMask.createBitShiftAndMaskInChain(28, bsmChrom);
	private static final BitSetUtils.BitShiftAndMask bsmAlleleRef  = BitShiftAndMask.createBitShiftAndMaskInChain(3, bsmPos);
	private static final BitSetUtils.BitShiftAndMask bsmAlleleVarN = BitShiftAndMask.createBitShiftAndMaskInChain(3, bsmAlleleRef);
	private static final BitSetUtils.BitShiftAndMask bsmAlleleVarT = BitShiftAndMask.createBitShiftAndMaskInChain(3, bsmAlleleVarN);
	private static final BitSetUtils.BitShiftAndMask bsmAlleleVarPop = BitShiftAndMask.createBitShiftAndMaskInChain(3, bsmAlleleVarT);		
	private static final BitSetUtils.BitShiftAndMask bsmMutType    = BitShiftAndMask.createBitShiftAndMaskInChain(5, bsmAlleleVarPop);
	
	public String mFlankingNormal;
	public String mFlankingTumor;
	public String mHugoSymbol;
	
	public short mCovgTotalNormal;
	public short mCovgTotalTumor;
	public short mCovgVarNormal;
	public short mCovgVarTumor;
	public int mRsID; 
	
	// ========================================================================
	public ClusteringInputOneSite() {
		clear();
	}			
	
	// ========================================================================
	public ClusteringInputOneSite(String line, SeqPlatform platform) {
		mFlankingNormal = StringUtils.extractNthColumnValue(line, Script.Col_NAFTAFInput_FlankingStringNormal,  StringUtils.FileExtensionTSV.mDelimiter);
		mFlankingTumor  = StringUtils.extractNthColumnValue(line, Script.Col_NAFTAFInput_FlankingStringTumor,   StringUtils.FileExtensionTSV.mDelimiter);
		
		mCovgTotalNormal = Short.parseShort(StringUtils.extractNthColumnValue(line, Script.Col_NAFTAFInput_TotalCoverageNormal,   StringUtils.FileExtensionTSV.mDelimiter));
		mCovgTotalTumor  = Short.parseShort(StringUtils.extractNthColumnValue(line, Script.Col_NAFTAFInput_TotalCoverageTumor,    StringUtils.FileExtensionTSV.mDelimiter));
		mCovgVarNormal   = Short.parseShort(StringUtils.extractNthColumnValue(line, Script.Col_NAFTAFInput_VariantCoverageNormal, StringUtils.FileExtensionTSV.mDelimiter));
		mCovgVarTumor    = Short.parseShort(StringUtils.extractNthColumnValue(line, Script.Col_NAFTAFInput_VariantCoverageTumor,  StringUtils.FileExtensionTSV.mDelimiter));

		// Now for rsID
		if (platform == SeqPlatform.SOLiD) {
			//variantAnnotation = (germCols[7].indexOf(Script.NovelStr) >= 0) ? Script.NovelStr : Utils.NAStr;
			// Sidd: For the NAStr case, strangely, the variant base is n/a in the SOLiD naf-taf-inputs 
			// (and there's not much point in looking up the reference base's allele frequency)
		} else if (platform == SeqPlatform.Illumina) {
			String dbsnpStr = StringUtils.extractNthColumnValue(line, Script.Col_NAFTAFInput_DbSNPString, StringUtils.FileExtensionTSV.mDelimiter); 
			if (dbsnpStr.indexOf(Script.NovelStr) >= 0) {
				mRsID = GenotypeUtils.RsID_Novel;
			} else {
				String rsNumTemp = GenotypeUtils.extractRsNumberFromLine(dbsnpStr);
				
				if ((rsNumTemp == null) || 
					 rsNumTemp.equals(Utils.rsPrefix) || 
					 rsNumTemp.equals(Utils.rsNegative)) {
					
					mRsID = GenotypeUtils.RsID_Unknown;
				} else {
					mRsID = GenotypeUtils.getNumberFromRsId(rsNumTemp);
				}
			}
		} else {
			mRsID = GenotypeUtils.RsID_Unknown;  // Could change if another platform added later
		}
	
		// Now for Hugo symbol and mutation type		
		if (platform == SeqPlatform.SOLiD) {
			String mutationTypeStr = StringUtils.extractNthColumnValue(line, 9, StringUtils.FileExtensionTSV.mDelimiter);
			mutationTypeStr.replace("syn", "synonymous").replace("nsynonymous", "nonsynonymous");
			mutationTypeStr = (mutationTypeStr == "") ? Utils.NAStr : mutationTypeStr;			
			mHugoSymbol = StringUtils.extractNthColumnValue(line, 8, StringUtils.FileExtensionTSV.mDelimiter);
			mHugoSymbol = (mHugoSymbol == "") ? Utils.NAStr : mHugoSymbol;

		} else if (platform == SeqPlatform.Illumina) {			
			String mutationTypeStr = StringUtils.extractNthColumnValue(line, Script.Col_NAFTAFInput_MutationType, StringUtils.FileExtensionTSV.mDelimiter);
			int mutTypeCode = MutationType.getSNVType(mutationTypeStr).ordinal();
			mDataUnit_ChromProsRevVarAllelesMutType = bsmMutType.setValueInCompactUnit(mutTypeCode, mDataUnit_ChromProsRevVarAllelesMutType);
			
			mHugoSymbol     = StringUtils.extractNthColumnValue(line, Script.Col_NAFTAFInput_HugoSymbol, StringUtils.FileExtensionTSV.mDelimiter);
			
		}
	
		// Now set the chrom, position, and alleles
		Chrom chrom  = Chrom.getChrom(   StringUtils.extractNthColumnValue(line, Script.Col_NAFTAFInput_Chrom,     StringUtils.FileExtensionTSV.mDelimiter) );
		int position = Integer.parseInt( StringUtils.extractNthColumnValue(line, Script.Col_NAFTAFInput_Position,  StringUtils.FileExtensionTSV.mDelimiter) );
		
		char nucChar = StringUtils.extractNthColumnValue(line, Script.Col_NAFTAFInput_AlleleRef,  StringUtils.FileExtensionTSV.mDelimiter).charAt(0);
		Nuc aRef     = nucChar == Script.MissingAllele ? Nuc.N : Nuc.getNuc(nucChar);
		
		nucChar      = StringUtils.extractNthColumnValue(line, Script.Col_NAFTAFInput_AlleleVarN, StringUtils.FileExtensionTSV.mDelimiter).charAt(0);
		Nuc aVarN    = nucChar == Script.MissingAllele ? Nuc.N : Nuc.getNuc(nucChar);
		
		nucChar      = StringUtils.extractNthColumnValue(line, Script.Col_NAFTAFInput_AlleleVarT, StringUtils.FileExtensionTSV.mDelimiter).charAt(0);
		Nuc aVarT    = nucChar == Script.MissingAllele ? Nuc.N : Nuc.getNuc(nucChar);
		
		nucChar      = StringUtils.extractNthColumnValue(line, Script.Col_NAFTAFInput_AlleleVarPop, StringUtils.FileExtensionTSV.mDelimiter).charAt(0);
		Nuc aVarPop  = nucChar == Script.MissingAllele ? Nuc.N : Nuc.getNuc(nucChar);
		
		setChrom(chrom);
		setPosition(position);			
		mDataUnit_ChromProsRevVarAllelesMutType = bsmAlleleRef.setValueInCompactUnit( aRef.getCode(),  mDataUnit_ChromProsRevVarAllelesMutType);
		mDataUnit_ChromProsRevVarAllelesMutType = bsmAlleleVarN.setValueInCompactUnit(aVarN.getCode(), mDataUnit_ChromProsRevVarAllelesMutType);
		mDataUnit_ChromProsRevVarAllelesMutType = bsmAlleleVarT.setValueInCompactUnit(aVarT.getCode(), mDataUnit_ChromProsRevVarAllelesMutType);
		mDataUnit_ChromProsRevVarAllelesMutType = bsmAlleleVarPop.setValueInCompactUnit(aVarPop.getCode(), mDataUnit_ChromProsRevVarAllelesMutType);
	}

	// ====================================================================
	public void setChrom(Chrom chrom) {
		mDataUnit_ChromProsRevVarAllelesMutType = bsmChrom.setValueInCompactUnit(chrom.ordinal(), mDataUnit_ChromProsRevVarAllelesMutType);
	}
	
	// ====================================================================
	public void setPosition(int position) {			
		mDataUnit_ChromProsRevVarAllelesMutType = bsmPos.setValueInCompactUnit(position,  mDataUnit_ChromProsRevVarAllelesMutType);
	}
	
	// ====================================================================
	public void setCovgTotalNormal(short mCovgTotalNormal) {
		this.mCovgTotalNormal = mCovgTotalNormal;
	}

	// ====================================================================
	public void setCovgTotalTumor(short mCovgTotalTumor) {
		this.mCovgTotalTumor = mCovgTotalTumor;
	}

	// ====================================================================
	public void setCovgVarNormal(short mCovgVarNormal) {
		this.mCovgVarNormal = mCovgVarNormal;
	}

	// ====================================================================
	public void setCovgVarTumor(short mCovgVarTumor) {
		this.mCovgVarTumor = mCovgVarTumor;
	}
	
	// ====================================================================
	public boolean refOrVarHasZeroReadCount() {
		return ((this.mCovgVarNormal == this.mCovgTotalNormal) || (this.mCovgVarNormal == 0));
	}
	
	// ====================================================================
	public float calcVAFNormal() { return (float) mCovgVarNormal / (float) mCovgTotalNormal; }
	public float calcVAFTumor()  { return (float) mCovgVarTumor  / (float) mCovgTotalTumor;  }
	public int getRsID()         { return mRsID; }
	
	public Chrom getChrom()  { return Chrom.getChrom((byte) bsmChrom.getValueInCompactUnit(mDataUnit_ChromProsRevVarAllelesMutType)); }
	public int getPosition() { return                ((int)   bsmPos.getValueInCompactUnit(mDataUnit_ChromProsRevVarAllelesMutType)); }
	
	public Nuc getReferenceAllele()     { return Nuc.getAllele(bsmAlleleRef.getValueInCompactUnit(mDataUnit_ChromProsRevVarAllelesMutType));  }
	public Nuc getVariantAlleleNormal() { return Nuc.getAllele(bsmAlleleVarN.getValueInCompactUnit(mDataUnit_ChromProsRevVarAllelesMutType)); }
	public Nuc getVariantAlleleTumor()  { return Nuc.getAllele(bsmAlleleVarT.getValueInCompactUnit(mDataUnit_ChromProsRevVarAllelesMutType)); }
	public Nuc getVariantAllelePopulation() { return Nuc.getAllele(bsmAlleleVarPop.getValueInCompactUnit(mDataUnit_ChromProsRevVarAllelesMutType)); }
	
	public MutationType getMutationType() {
		return MutationType.getMutationType((int) bsmMutType.getValueInCompactUnit(mDataUnit_ChromProsRevVarAllelesMutType));
	}

	// ====================================================================
	public void clear() {
		mDataUnit_ChromProsRevVarAllelesMutType = 0;
		mFlankingTumor = mFlankingNormal = mHugoSymbol = "";
		mCovgTotalNormal = mCovgTotalTumor = mCovgVarNormal = mCovgVarTumor = 0;
		mRsID = GenotypeUtils.RsID_Unknown;
	}
	
	// ====================================================================
	public int compareTo(ClusteringInputOneSite rhs) {
		// The compact unit variable has chromosome and position ordered such that it's easy to sort
		return Long.compare(mDataUnit_ChromProsRevVarAllelesMutType, rhs.mDataUnit_ChromProsRevVarAllelesMutType);
	}
	
	// ====================================================================
	
	public StringBuilder printToString(StringBuilder sb, boolean clearStringBuilder, String delimiter) {
		if (clearStringBuilder) {
			sb.setLength(0);
		}
		
		sb.append("chr").append(getChrom().ordinal())
		  .append(delimiter).append(getPosition())
		  .append(delimiter).append(getReferenceAllele())
		  .append(delimiter).append(getVariantAllelePopulation())
		  .append(delimiter).append(getVariantAlleleNormal())
		  .append(delimiter).append(getVariantAlleleTumor())
		  .append(delimiter).append(mFlankingNormal)
		  .append(delimiter).append(mFlankingTumor)
		  .append(delimiter).append(mCovgTotalNormal)
		  .append(delimiter).append(mCovgTotalTumor)
		  .append(delimiter).append(mCovgVarNormal)
		  .append(delimiter).append(mCovgVarTumor)
		  .append(delimiter).append(calcVAFNormal())
		  .append(delimiter).append(calcVAFTumor());
		
		if (mRsID == GenotypeUtils.RsID_Novel) {
			sb.append(delimiter).append(Script.NovelStr);
		} else {
			sb.append(delimiter).append(Utils.rsPrefix).append(mRsID);
		}
		
		sb.append(delimiter).append(getMutationType().getPrintName());
		sb.append(delimiter).append(mHugoSymbol);
		return sb;
	}
	
	// ====================================================================		
	public static void TestClusteringInputOneSite_Robust() {
		int numTrials = 10000000;
		
		ClusteringInputOneSite oneSiteInfo = new ClusteringInputOneSite();
		for (int trial = 0; trial < numTrials; trial++) {
			Chrom chrom = Chrom.getChrom((byte) NumberUtils.getRandomInteger(0, Chrom.values().length - 1));
			int position = NumberUtils.getRandomInteger(0, (int) BitSetUtils.getMask(28));
			oneSiteInfo.setChrom(chrom);
			oneSiteInfo.setPosition(position);
			
			//System.out.println(trial);
			CompareUtils.ensureTrue(oneSiteInfo.getChrom().equals(chrom), "Chrom not equal!");
			CompareUtils.ensureTrue(position == oneSiteInfo.getPosition(), "Position not equal!");
		}
	}
}