package lohcate.geneAnalysis;
import genomeEnums.Chrom;
import genomeEnums.VariantLocation;

import java.util.ArrayList;

import org.apache.commons.math3.util.MathArrays;

import nutils.ArrayUtils;
import nutils.EnumMapSafe;
import nutils.NumberUtils;
import nutils.counter.BucketCounterEnum;
import nutils.counter.DynamicBucketCounter;

import lohcate.EventTypeAllele;
import lohcateEnums.EventType;
import lohcateEnums.MutationType;

import shared.Utils;
/**
 * LOHcate --- A software tool for LOH calling and visualization in cancer genomes 
 * Human Genome Sequencing Center, Baylor College of Medicine (Houston, TX)
 * 
 * @author Ninad Dewal & Siddharth G. Reddy
 *
 */
public class GeneCounter implements Comparable<GeneCounter> {	
	
	public String mLabel;
	public Chrom mChrom;
	public int mNumericID;
	public EnumMapSafe<EventType, ArrayList<String>> mPatients; 	
	
	public BucketCounterEnum<MutationType> mMutationTypeCounts;	
	public BucketCounterEnum<EventType> mEventTypeCounts;	
	public BucketCounterEnum<VariantLocation> mVariantLocationCounts;
	
	public EnumMapSafe<EventTypeAllele, ArrayList<String>> mAlleleEventType;
	
	public int mCountLOHreferenceLost;
	public int mMinBasePairPosition, mMaxBasePairPosition;
	
	EnumMapSafe<EventTypeAllele, DynamicBucketCounter> mVariantAlleleCounts;
	
	
	
	public GeneCounter(String name, Chrom chrom, int numericID) {		
		mLabel = name;
		mChrom = chrom;
		mNumericID = numericID;
		initializePatients();
				
		mMutationTypeCounts    = new BucketCounterEnum<MutationType>(MutationType.class);  // stores hit counts for synonymous, nonsynonymous
		mEventTypeCounts       = new BucketCounterEnum<EventType>(EventType.class);				
		mVariantLocationCounts = new BucketCounterEnum<VariantLocation>(VariantLocation.class);  // stores hit counts for germline, somatic
		
		mAlleleEventType = EnumMapSafe.createEnumMapOfArrayLists(EventTypeAllele.class, String.class); 				
		
		clearCounts();
				
		mVariantAlleleCounts = DynamicBucketCounter.ClassFactory.newEnumMap(EventTypeAllele.class);	
				
		mMinBasePairPosition = Integer.MAX_VALUE;
		mMaxBasePairPosition = Integer.MIN_VALUE;
	}
	
	public void clearCounts() {
		mMutationTypeCounts.clear();
		mEventTypeCounts.clear();
		mVariantLocationCounts.clear();		
		//mAlleleEventTypeCounts.clear();
		mCountLOHreferenceLost = 0;
	}
	
	public int getNumericID() { return mNumericID; }
	
	private void initializePatients() {
		mPatients = EnumMapSafe.createEnumMapOfArrayLists(EventType.class, String.class);		
	}
	
	/** Adds a patient to the reference-lost list in the case of LOH. */
	public boolean addPatientIfNotAlreadyAdded_byAlleleEvent(EventTypeAllele eventForAllele, String patientName, int position) {
		mVariantAlleleCounts.get(eventForAllele).incrementCount(position);
		return ArrayUtils.checkInListAndInsertIfMissing(mAlleleEventType.get(eventForAllele), patientName);
	}
	
	/** Returns true if the patient was already existing for the clustertype, false otherwise. */
	public boolean addPatientIfNotAlreadyAdded(String patientName, EventType eventType) {		
		ArrayList<String> patientList = getPatientsForEventType(eventType);		
		return ArrayUtils.checkInListAndInsertIfMissing(patientList, patientName);
	}
	
	public int getNumPatientsForEventType(EventType eventType) {
		ArrayList<String> patientList = getPatientsForEventType(eventType);
		return patientList.size();
	}
	
	public ArrayList<String> getPatientsForEventType(EventType ct) {
		return mPatients.get(ct);
	}
	
	public void incrementCount(EventType eventType)    { mEventTypeCounts.increment(eventType); }	
	public void incrementCount(MutationType mutationType) { mMutationTypeCounts.increment(mutationType); } 			
	public void incrementCount(VariantLocation varLoc) { mVariantLocationCounts.increment(varLoc); }
	
	public int getCount(EventType eventType)    { return mEventTypeCounts.getCount(eventType); }
	public int getCount(MutationType mutationType) { return mMutationTypeCounts.getCount(mutationType); }
	public int getCount(VariantLocation varLoc) { return mVariantLocationCounts.getCount(varLoc); }
	
		
	
	public float getRecurrence(EventType ct, int total) { return (float) mPatients.get(ct).size() / (float) total; }
	
	public int getRangeLength() { return (mMaxBasePairPosition - mMinBasePairPosition + 1); }
	
	public String getName() { return mLabel; }
	
	public float getDensityVariantLocation(VariantLocation varLoc) {
		return (float) getCount(varLoc) / (float) getRangeLength();
	}
	
	public float getDensityMutationType(MutationType mutationType) {
		return (float) getCount(mutationType) / (float) getRangeLength();
	}
	
	public float getDensityEventType(EventType eventType) {
		return (float) getCount(eventType) / (float) getRangeLength(); 
	}
	
	public static String getHeaders(String delim) {
		//String logStr = "_log";
		//String recurrenceStr = "_recurrence";
		//String densityStr = "_density";
		String inGeneSuffix = "_in_gene";
		String numberSignPrefix = "#_";
		String numberSignPtsPrefix = numberSignPrefix + "pts_w_";
		String pValPrefix = "pval_";
		String locationsSuffix = "_locations";
		String eventsTotalSuffix = "_events_total";
		
		StringBuilder sb = new StringBuilder(2048);
		
		sb.append(numberSignPrefix).append(MutationType.NonSynonymous_SNV.toLowerCase()).append(locationsSuffix)
		  .append(delim).append(numberSignPrefix).append(MutationType.Synonymous_SNV.toLowerCase()).append(locationsSuffix);
		  
		  for (EventType et : EventType.AmpLOHcnLOHhetSomaticDelHom) {
			  sb.append(delim).append(numberSignPrefix).append(et.name()).append(eventsTotalSuffix);
		  }
		  
		  for (EventType et : EventType.AmpLOHcnLOHhetSomaticDelHom) {
			  sb.append(delim).append("ln(").append(numberSignPrefix).append(et.name()).append(eventsTotalSuffix).append(")");
		  }
		  
		  for (EventType et : EventType.AmpLOHcnLOHhetSomaticDelHom) {
			  sb.append(delim).append(et.name()).append("_density");
		  }

		  for (EventType et : EventType.AmpLOHcnLOHhetSomaticDelHom) {
			  sb.append(delim).append(numberSignPtsPrefix).append(et.name()).append(inGeneSuffix);			  
		  }
		  
		  for (EventTypeAllele eta : EventTypeAllele.values()) {
			  sb.append(delim).append(numberSignPtsPrefix).append(eta.name()).append(inGeneSuffix);			  
		  }
		  
		  sb.append(delim).append(numberSignPtsPrefix).append("LossReference").append(inGeneSuffix);
		  sb.append(delim).append(numberSignPtsPrefix).append("LossVariant").append(inGeneSuffix);
		  
		  sb.append(delim).append(pValPrefix).append("LossReference").append(inGeneSuffix);
		  sb.append(delim).append(pValPrefix).append("LossReference_when_VariantRare").append(inGeneSuffix);
		  sb.append(delim).append(pValPrefix).append("LossReference_when_VariantCommon").append(inGeneSuffix);
		  sb.append(delim).append(pValPrefix).append("VariantRare_when_LossReference").append(inGeneSuffix);
		  sb.append(delim).append(pValPrefix).append("VariantCommon_when_LossVariant").append(inGeneSuffix);
		   
		  
		  for (EventTypeAllele eta : EventTypeAllele.values()) {
			  sb.append(delim).append(eta.name()).append("_Positions_in_").append(EventType.LOH.name()).append("_or_").append(EventType.cnLOH.name());			  
		  }		  

		//EventType.LOH.name() + "_refLost" + logStr, 		
		//EventType.LOH.name() + "_refLost" + densityStr, 
		  return sb.toString();		
	}
		
	public String countsToString(String delim) {
		StringBuilder sb = new StringBuilder(8192);

		int refLostTotal = mAlleleEventType.get(EventTypeAllele.LossReference_VariantRare).size() + mAlleleEventType.get(EventTypeAllele.LossReference_VariantCommon).size();
		int varLostTotal = mAlleleEventType.get(EventTypeAllele.LossVariant_VariantRare).size()   + mAlleleEventType.get(EventTypeAllele.LossVariant_VariantCommon).size();
		int varRareTotal = mAlleleEventType.get(EventTypeAllele.LossReference_VariantRare).size() + mAlleleEventType.get(EventTypeAllele.LossVariant_VariantRare).size();
		int varCommonTotal = mAlleleEventType.get(EventTypeAllele.LossReference_VariantCommon).size() + mAlleleEventType.get(EventTypeAllele.LossVariant_VariantCommon).size();		
		int lostTotal = refLostTotal + varLostTotal;
		
		sb.append(getCount(MutationType.NonSynonymous_SNV))
		  .append(delim).append(getCount(MutationType.Synonymous_SNV));
		  
		for (EventType et : EventType.AmpLOHcnLOHhetSomaticDelHom) {
			sb.append(delim).append(getCount(et));
		}

		for (EventType et : EventType.AmpLOHcnLOHhetSomaticDelHom) {
			sb.append(delim).append(Utils.log(getCount(et)));					  
		}
		
		for (EventType et : EventType.AmpLOHcnLOHhetSomaticDelHom) {
			sb.append(delim).append(getDensityEventType(et));					  
		}
		
		for (EventType et : EventType.AmpLOHcnLOHhetSomaticDelHom) {
			sb.append(delim).append(getNumPatientsForEventType(et));					  
		}

		// Report the counts from the event type allele options
		for (EventTypeAllele eta : EventTypeAllele.values()) {
			sb.append(delim).append(mAlleleEventType.get(eta).size());					  
		}
		  
		// Report the totals
		sb.append(delim).append(refLostTotal);
		sb.append(delim).append(varLostTotal);
		  
		// Comparison is
		sb.append(delim).append(NumberUtils.cumulativeProbabilitySuccessBinomial(lostTotal, refLostTotal, 0.5));
		sb.append(delim).append(NumberUtils.cumulativeProbabilitySuccessBinomial(varRareTotal,   mAlleleEventType.get(EventTypeAllele.LossReference_VariantRare).size(),   0.5));
		sb.append(delim).append(NumberUtils.cumulativeProbabilitySuccessBinomial(varCommonTotal, mAlleleEventType.get(EventTypeAllele.LossReference_VariantCommon).size(), 0.5));
		sb.append(delim).append(NumberUtils.cumulativeProbabilitySuccessBinomial(refLostTotal,   mAlleleEventType.get(EventTypeAllele.LossReference_VariantRare).size(),   0.5));
		sb.append(delim).append(NumberUtils.cumulativeProbabilitySuccessBinomial(varLostTotal,   mAlleleEventType.get(EventTypeAllele.LossVariant_VariantCommon).size(),   0.5));
		  		  		
		
		for (EventTypeAllele eta : EventTypeAllele.values()) {
			sb.append(delim).append(mVariantAlleleCounts.get(eta).toString());			  
		}
		 // .append(delim).append(mPatientsWithRefLost.size() / (float) getRangeLength())
		  		
		return sb.toString();
	}
	
	public int compareTo(GeneCounter rhs) {
		return mLabel.compareToIgnoreCase(rhs.mLabel);
	}
	
	/* (non-Javadoc)
	 * @see java.lang.Object#hashCode()
	 */
	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + ((mLabel == null) ? 0 : mLabel.hashCode());
		return result;
	}

	/* (non-Javadoc)
	 * @see java.lang.Object#equals(java.lang.Object)
	 */
	@Override
	public boolean equals(Object obj) {
		if (this == obj) return true;
		if (obj == null) return false;
		if (getClass() != obj.getClass()) return false;
		
		GeneCounter other = (GeneCounter) obj;
		if (mLabel == null) {
			if (other.mLabel != null)
				return false;
		} else if (mLabel.compareToIgnoreCase(other.mLabel) != 0) {				
			return false;
		} 
		return true;
	}
	
}