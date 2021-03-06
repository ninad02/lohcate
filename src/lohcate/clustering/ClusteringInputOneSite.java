package lohcate.clustering;

import genomeEnums.Chrom;
import genomeEnums.Nuc;
import genomeEnums.VariantFrequency;
import genomeUtils.GenotypeUtils;
import genomeUtils.SiteInformation;

import java.util.Comparator;
import java.util.HashSet;

import lohcate.Regions;
import lohcateEnums.MutationType;
import lohcateEnums.SeqPlatform;
import nutils.CompareUtils;
import nutils.NumberUtils;
import nutils.StringUtils;
import nutils.BitUtils.BitSetUtils;
import nutils.BitUtils.BitShiftAndMask;
import nutils.collectionsSorted.ArrayListSortedComparable;
import shared.Utils;

public class ClusteringInputOneSite implements Comparable<ClusteringInputOneSite>, SiteInformation.Writeable {

	protected static ArrayListSortedComparable<String> HugoSymbolArray = new ArrayListSortedComparable<String>();
	
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
	
	private static final BitShiftAndMask bsmChrom = new BitShiftAndMask(5, 58);
	private static final BitShiftAndMask bsmPos   = BitShiftAndMask.createBitShiftAndMaskInChain(28, bsmChrom);
	private static final BitShiftAndMask bsmAlleleRef  = BitShiftAndMask.createBitShiftAndMaskInChain(3, bsmPos);
	private static final BitShiftAndMask bsmAlleleVarN = BitShiftAndMask.createBitShiftAndMaskInChain(3, bsmAlleleRef);
	private static final BitShiftAndMask bsmAlleleVarT = BitShiftAndMask.createBitShiftAndMaskInChain(3, bsmAlleleVarN);
	private static final BitShiftAndMask bsmAlleleVarPop = BitShiftAndMask.createBitShiftAndMaskInChain(3, bsmAlleleVarT);		
	private static final BitShiftAndMask bsmMutType    = BitShiftAndMask.createBitShiftAndMaskInChain(5, bsmAlleleVarPop);
	
	public String mFlankingNormal;
	public String mFlankingTumor;
	public String mHugoSymbol;
	
	public short mCovgTotalNormal;
	public short mCovgTotalTumor;
	public short mCovgVarNormal;
	public short mCovgVarTumor;
	public int mRsID; 
	public VariantFrequency mVarFrequency;

	// ========================================================================
	/** Switches the tumor data with the normal data. */
	public void switchTumorAndNormal() {
		
		// Swap the coverages		
		short tempCovg = mCovgTotalNormal;
		mCovgTotalNormal = mCovgTotalTumor;
		mCovgTotalTumor = tempCovg;
		
		tempCovg = mCovgVarNormal;
		mCovgVarNormal = mCovgVarTumor;
		mCovgVarTumor = tempCovg;
		
		String tempFlanking = mFlankingNormal;
		mFlankingNormal = mFlankingTumor;
		mFlankingTumor = tempFlanking;
		
		Nuc tempNucVarNormal = getVariantAlleleNormal();
		Nuc tempNucVarTumor  = getVariantAlleleTumor();
		setVariantAlleleNormal(tempNucVarTumor);
		setVariantAlleleTumor(tempNucVarNormal);
	}
	
	public void copyTumorIntoNormal() {
		mCovgTotalNormal = mCovgTotalTumor;
		mCovgVarNormal   = mCovgVarTumor;
		mFlankingNormal  = mFlankingTumor;
		setVariantAlleleNormal(getVariantAlleleTumor());
	}
	
	public void copyNormalIntoTumor() {
		mCovgTotalTumor = mCovgTotalNormal;
		mCovgVarTumor   = mCovgVarNormal;
		mFlankingTumor  = mFlankingNormal;
		
		setVariantAlleleTumor(getVariantAlleleNormal());
	}
	
	// ========================================================================
	public ClusteringInputOneSite() {
		clear();
	}			

	// ========================================================================
	private static short parseShortHelper(String theString) {		
		if (theString.charAt(0) == Regions.MissingAllele || theString.isEmpty()) {
			return 0;
		} else {
			return Short.parseShort(theString);
		}
	}
	
	// ========================================================================
	public ClusteringInputOneSite(String line, SeqPlatform platform) {
		String delimiter = StringUtils.FileExtensionTSV.mDelimiter;
		
		mFlankingNormal = StringUtils.extractNthColumnValue(line, Regions.Col_NAFTAFInput_FlankingStringNormal, delimiter);
		mFlankingTumor  = StringUtils.extractNthColumnValue(line, Regions.Col_NAFTAFInput_FlankingStringTumor,  delimiter);
		
		mCovgTotalNormal = parseShortHelper(StringUtils.extractNthColumnValue(line, Regions.Col_NAFTAFInput_TotalCoverageNormal,   delimiter));
		mCovgTotalTumor  = parseShortHelper(StringUtils.extractNthColumnValue(line, Regions.Col_NAFTAFInput_TotalCoverageTumor,    delimiter));
		mCovgVarNormal   = parseShortHelper(StringUtils.extractNthColumnValue(line, Regions.Col_NAFTAFInput_VariantCoverageNormal, delimiter));
		mCovgVarTumor    = parseShortHelper(StringUtils.extractNthColumnValue(line, Regions.Col_NAFTAFInput_VariantCoverageTumor,  delimiter));

		// Now set the chrom, position, and alleles
		Chrom chrom  = Chrom.getChrom(   StringUtils.extractNthColumnValue(line, Regions.Col_NAFTAFInput_Chrom,     delimiter) );
		int position = Integer.parseInt( StringUtils.extractNthColumnValue(line, Regions.Col_NAFTAFInput_Position,  delimiter) );
		
		char nucChar = StringUtils.extractNthColumnValue(line, Regions.Col_NAFTAFInput_AlleleRef,  delimiter).charAt(0);
		Nuc aRef     = nucChar == Regions.MissingAllele ? Nuc.N : Nuc.getNuc(nucChar);
		
		nucChar      = StringUtils.extractNthColumnValue(line, Regions.Col_NAFTAFInput_AlleleVarN, delimiter).charAt(0);
		Nuc aVarN    = nucChar == Regions.MissingAllele ? Nuc.N : Nuc.getNuc(nucChar);
		
		nucChar      = StringUtils.extractNthColumnValue(line, Regions.Col_NAFTAFInput_AlleleVarT, delimiter).charAt(0);
		Nuc aVarT    = nucChar == Regions.MissingAllele ? Nuc.N : Nuc.getNuc(nucChar);
		
		nucChar      = StringUtils.extractNthColumnValue(line, Regions.Col_NAFTAFInput_AlleleVarPop, delimiter).charAt(0);
		Nuc aVarPop  = nucChar == Regions.MissingAllele ? Nuc.N : Nuc.getNuc(nucChar);
		
		setChrom(chrom);
		setPosition(position);			
		setReferenceAllele(aRef);
		setVariantAlleleNormal(aVarN);
		setVariantAlleleTumor(aVarT);
		setVariantAllelePopulation(aVarPop);
		
		// Now for rsID
		mVarFrequency = VariantFrequency.Unknown;
		if (platform == SeqPlatform.SOLiD) {
			//variantAnnotation = (germCols[7].indexOf(Script.NovelStr) >= 0) ? Script.NovelStr : Utils.NAStr;
			// Sidd: For the NAStr case, strangely, the variant base is n/a in the SOLiD naf-taf-inputs 
			// (and there's not much point in looking up the reference base's allele frequency)
		} else if (platform == SeqPlatform.Illumina) {
			String dbsnpStr = StringUtils.extractNthColumnValue(line, Regions.Col_NAFTAFInput_DbSNPString, StringUtils.FileExtensionTSV.mDelimiter); 
			if (dbsnpStr.indexOf(Regions.NovelStr) >= 0) {
				mRsID = GenotypeUtils.RsID_Novel;
				mVarFrequency = VariantFrequency.Novel;
			} else {
				String rsNumTemp = GenotypeUtils.extractRsNumberFromLine(dbsnpStr);
				
				if ((rsNumTemp == null) || 
					 rsNumTemp.equals(Utils.rsPrefix) || 
					 rsNumTemp.equals(Utils.rsNegative)) {
					
					mRsID = GenotypeUtils.RsID_Unknown;
				} else {
					mRsID = GenotypeUtils.getNumberFromRsId(rsNumTemp);
										
					String[] items = dbsnpStr.split(",");  // Split on the comma
					String[] alleles = items[2].split(";|\\|");
					int numElementsTogether = 2;
					if (alleles.length < numElementsTogether) {
						mVarFrequency = VariantFrequency.Unknown;
					} else if (alleles.length == numElementsTogether) {
						mVarFrequency = VariantFrequency.Unknown;
					} else if (alleles.length >= numElementsTogether) {
						mVarFrequency = VariantFrequency.Common;
						
						double varFreqLowest = 1.0;
						Nuc varFreqLowestNuc = Nuc.N;
						Nuc tempNuc = Nuc.N;
						for (int i = 0; i < alleles.length; i++) {
							String allele = alleles[i];
							if (NumberUtils.isEven(i)) {  // is at allele position
								if (allele.length() > 0 && Character.isLetter(allele.charAt(0))) {
									tempNuc = Nuc.getNuc(allele.charAt(0));
								}
							} else {
								if (allele.length() > 0) {
									double varFreq = Double.parseDouble(allele); 
									if (varFreq < varFreqLowest) {
										varFreqLowest = varFreq;
										varFreqLowestNuc = tempNuc;
									}
								}
								tempNuc = Nuc.N; // reset the allele
							}
						}
						
						// check that we're not the reference
						if (varFreqLowestNuc == aRef) {
							varFreqLowest = 1.0 - varFreqLowest;  // get variant allele fraction in population 
						}
						
						if (varFreqLowest < 0.01) {
							mVarFrequency = VariantFrequency.Rare;
						}
					}
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
			String mutationTypeStr = StringUtils.extractNthColumnValue(line, Regions.Col_NAFTAFInput_MutationType, StringUtils.FileExtensionTSV.mDelimiter);
			int mutTypeCode = MutationType.getSNVType(mutationTypeStr).ordinal();
			mDataUnit_ChromProsRevVarAllelesMutType = bsmMutType.setValueInCompactUnit(mutTypeCode, mDataUnit_ChromProsRevVarAllelesMutType);
			
			mHugoSymbol = StringUtils.extractNthColumnValue(line, Regions.Col_NAFTAFInput_HugoSymbol, StringUtils.FileExtensionTSV.mDelimiter);
			mHugoSymbol = findAndAddHugoSymbol(mHugoSymbol); 
			
		}
	
	}

	// ====================================================================
	public void set(Chrom chrom, int position) {
		setChrom(chrom);
		setPosition(position);
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
	public short getCovgTotalNormal() { return this.mCovgTotalNormal; }
	
	// ====================================================================
	public short getCovgTotalTumor()  { return this.mCovgTotalTumor; }
	
	// ====================================================================
	public void setCovgTotalNormal(short covgTotalNormal) {
		this.mCovgTotalNormal = covgTotalNormal;
	}

	// ====================================================================
	public void setCovgTotalTumor(short mCovgTotalTumor) {
		this.mCovgTotalTumor = mCovgTotalTumor;
	}

	// ====================================================================
	public void setCovgVarNormal(short covgVarNormal) {
		this.mCovgVarNormal = covgVarNormal;
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
	
	public Chrom getChrom()  { return Chrom.getChrom((byte) bsmChrom.extractValue(mDataUnit_ChromProsRevVarAllelesMutType)); }
	public int getPosition() { return                ((int)   bsmPos.extractValue(mDataUnit_ChromProsRevVarAllelesMutType)); }
	
	public String getHugoSymbol() { return mHugoSymbol; }
	
	public Nuc getReferenceAllele()     { return Nuc.getAllele(bsmAlleleRef.extractValue(mDataUnit_ChromProsRevVarAllelesMutType));  }
	public Nuc getVariantAlleleNormal() { return Nuc.getAllele(bsmAlleleVarN.extractValue(mDataUnit_ChromProsRevVarAllelesMutType)); }
	public Nuc getVariantAlleleTumor()  { return Nuc.getAllele(bsmAlleleVarT.extractValue(mDataUnit_ChromProsRevVarAllelesMutType)); }
	public Nuc getVariantAllelePopulation() { return Nuc.getAllele(bsmAlleleVarPop.extractValue(mDataUnit_ChromProsRevVarAllelesMutType)); }
	
	public void setReferenceAllele(Nuc aRef) {
		mDataUnit_ChromProsRevVarAllelesMutType = bsmAlleleRef.setValueInCompactUnit(aRef.getCode(),  mDataUnit_ChromProsRevVarAllelesMutType);
	}
	
	public void setVariantAlleleNormal(Nuc aVarN) {
		mDataUnit_ChromProsRevVarAllelesMutType = bsmAlleleVarN.setValueInCompactUnit(aVarN.getCode(), mDataUnit_ChromProsRevVarAllelesMutType);
	}
	
	public void setVariantAlleleTumor(Nuc aVarT) {
		mDataUnit_ChromProsRevVarAllelesMutType = bsmAlleleVarT.setValueInCompactUnit(aVarT.getCode(), mDataUnit_ChromProsRevVarAllelesMutType);
	}
	
	public void setVariantAllelePopulation(Nuc aVarPop) {
		mDataUnit_ChromProsRevVarAllelesMutType = bsmAlleleVarPop.setValueInCompactUnit(aVarPop.getCode(), mDataUnit_ChromProsRevVarAllelesMutType);
	}
	
	public MutationType getMutationType() {
		return MutationType.getMutationType((int) bsmMutType.extractValue(mDataUnit_ChromProsRevVarAllelesMutType));
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
	public String toString() {		
		return printToString(new StringBuilder(2048), true, StringUtils.FileExtensionTSV.mDelimiter).toString();
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
			sb.append(delimiter).append(Regions.NovelStr);
		} else {
			sb.append(delimiter).append(Utils.rsPrefix).append(mRsID);
		}
		sb.append(";").append(mVarFrequency);
		
		sb.append(delimiter).append(getMutationType().getPrintName());
		sb.append(delimiter).append(mHugoSymbol);
		return sb;
	}
	
	// ====================================================================
	/** Checks if hugo symbol exists in table.  If it does, it returns the copy in the table.  If not, it adds 
	 *  the symbol and returns the same input argument hugosymbol string. 
	 * @param hugoSymbol
	 * @return
	 */
	protected synchronized String findAndAddHugoSymbol(String hugoSymbol) {
		String result = HugoSymbolArray.get(hugoSymbol);
		if (result == null) {
			HugoSymbolArray.add(hugoSymbol);
			result = hugoSymbol;
		}
		return result;
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

	
	// ====================================================================
	// ====================================================================		
	// INTERFACE IMPLEMENTATIONS: Floint
	// ====================================================================
	// ====================================================================
	/*
	@Override
	public float getX() { return calcVAFTumor(); }

	// ====================================================================
	@Override
	public float getY() { return calcVAFNormal(); }

	@Override
	public float getZ() { return 0; }

	@Override
	public double getCartesianDistanceSquared(Floint rhs) {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	public double getCartesianDistance(Floint rhs) { return Math.sqrt(getCartesianDistanceSquared(rhs)); }

	@Override
	public double getTheoreticalFurthestDifferenceXOrYWithinDistance(double distance) { return distance; }
	*/
}