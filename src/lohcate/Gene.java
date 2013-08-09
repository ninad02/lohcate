package lohcate;
import genomeEnums.Chrom;
import genomeEnums.VariantLocation;

import java.util.ArrayList;
import nutils.ArrayUtils;
import nutils.EnumMapSafe;
import nutils.counter.BucketCounterEnum;
import nutils.counter.DynamicBucketCounter;

import lohcateEnums.EventType;
import lohcateEnums.MutationType;

import shared.Utils;
/**
 * LOHcate --- A software tool for LOH calling and visualization in cancer genomes
 * D Wheeler & SG Reddy
 * Human Genome Sequencing Center, Baylor College of Medicine (Houston, TX)
 * 
 * @author Siddharth G. Reddy
 *
 */
public class Gene implements Comparable<Gene> {	
	
	public String mLabel;
	public Chrom mChrom;
	public EnumMapSafe<EventType, ArrayList<String>> mPatients; 
	public ArrayList<String> mPatientsWithRefLost;
	
	public BucketCounterEnum<MutationType> mMutationTypeCounts;	
	public BucketCounterEnum<EventType> mClusterTypeCounts;	
	public BucketCounterEnum<VariantLocation> mVariantLocationCounts;
	public int mCountLOHreferenceLost;
	public int mMinBasePairPosition, mMaxBasePairPosition;
	
	DynamicBucketCounter mVariantAlleleCounts;
	
	public Gene(String name, Chrom chrom) {		
		mLabel = name;
		mChrom = chrom;
		initializePatients();
				
		mMutationTypeCounts    = new BucketCounterEnum<MutationType>(MutationType.class);  // stores hit counts for synonymous, nonsynonymous
		mClusterTypeCounts     = new BucketCounterEnum<EventType>(EventType.class);				
		mVariantLocationCounts = new BucketCounterEnum<VariantLocation>(VariantLocation.class);  // stores hit counts for germline, somatic
		clearCounts();
		
		mVariantAlleleCounts = new DynamicBucketCounter();
				
		mMinBasePairPosition = Integer.MAX_VALUE;
		mMaxBasePairPosition = Integer.MIN_VALUE;
	}
	
	public void clearCounts() {
		mMutationTypeCounts.clear();
		mClusterTypeCounts.clear();
		mVariantLocationCounts.clear();		
		mCountLOHreferenceLost = 0;
	}
	
	private void initializePatients() {
		mPatients = EnumMapSafe.createEnumMapOfArrayLists(EventType.class, String.class);
		mPatientsWithRefLost = new ArrayList<String>();
	}
	
	/** Adds a patient to the reference-lost list in the case of LOH. */
	public boolean addPatientIfNotAlreadyAdded_LOHRefLost(String patientName, int position) {
		mVariantAlleleCounts.incrementCount(position);
		return ArrayUtils.checkInListAndInsertIfMissing(mPatientsWithRefLost, patientName);
	}
	
	/** Returns true if the patient was already existing for the clustertype, false otherwise. */
	public boolean addPatientIfNotAlreadyAdded(String patientName, EventType clusterType) {		
		ArrayList<String> patientList = getPatientsForClusterType(clusterType);		
		return ArrayUtils.checkInListAndInsertIfMissing(patientList, patientName);
	}
	
	public int getNumPatientsForClusterType(EventType clusterType) {
		ArrayList<String> patientList = getPatientsForClusterType(clusterType);
		return patientList.size();
	}
	
	public ArrayList<String> getPatientsForClusterType(EventType ct) {
		return mPatients.get(ct);
	}
	
	public void incrementCount(EventType clusterType)    { mClusterTypeCounts.increment(clusterType); }	
	public void incrementCount(MutationType mutationType) { mMutationTypeCounts.increment(mutationType); } 			
	public void incrementCount(VariantLocation varLoc) { mVariantLocationCounts.increment(varLoc); }
	
	public int getCount(EventType clusterType)    { return mClusterTypeCounts.getCount(clusterType); }
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
	
	public float getDensityClusterType(EventType clusterType) {
		return (float) getCount(clusterType) / (float) getRangeLength(); 
	}
		
	public String countsToString(String delim) {
		StringBuilder sb = new StringBuilder(8192);
		
		sb.append(getCount(MutationType.NonSynonymous_SNV))
		  .append(delim).append(getCount(MutationType.Synonymous_SNV))
		  .append(delim).append(getCount(VariantLocation.Germline))
		  .append(delim).append(getCount(VariantLocation.Somatic))
		  
		  .append(delim).append(getCount(EventType.GainSomatic))
		  .append(delim).append(getCount(EventType.LOH))
		  .append(delim).append(mPatientsWithRefLost.size())
		  .append(delim).append(mVariantAlleleCounts.toString())
		  .append(delim).append(getCount(EventType.HETGermline))
		  .append(delim).append(getCount(EventType.HETSomatic))
		  
		  .append(delim).append(Utils.log(getCount(EventType.GainSomatic)))
		  .append(delim).append(Utils.log(getCount(EventType.LOH)))
		  .append(delim).append(Utils.log(mPatientsWithRefLost.size()))
		  .append(delim).append(Utils.log(getCount(EventType.HETGermline)))
		  .append(delim).append(Utils.log(getCount(EventType.HETSomatic)))
		  
		  .append(delim).append(getDensityClusterType(EventType.GainSomatic))
		  .append(delim).append(getDensityClusterType(EventType.LOH))
		  .append(delim).append(mPatientsWithRefLost.size() / (float) getRangeLength())
		  .append(delim).append(getDensityClusterType(EventType.HETGermline))
		  .append(delim).append(getDensityClusterType(EventType.HETSomatic))
		  
		  .append(delim).append(getNumPatientsForClusterType(EventType.GainSomatic))
		  .append(delim).append(getNumPatientsForClusterType(EventType.LOH))
		  .append(delim).append(getNumPatientsForClusterType(EventType.HETGermline))
		  .append(delim).append(getNumPatientsForClusterType(EventType.HETSomatic));
		
		return sb.toString();
	}
	
	public int compareTo(Gene rhs) {
		return mLabel.compareToIgnoreCase(rhs.mLabel);
	}
}