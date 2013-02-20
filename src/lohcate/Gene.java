package lohcate;
import genomeEnums.Chrom;
import genomeEnums.VariantLocation;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.EnumMap;

import nutils.ArrayUtils;
import nutils.EnumMapSafe;
import nutils.counter.BucketCounterCore;
import nutils.counter.BucketCounterEnum;

import lohcateEnums.ClusterType;
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
	public EnumMapSafe<ClusterType, ArrayList<String>> mPatients;  	
	
	public BucketCounterEnum<MutationType> mMutationTypeCounts;	
	public BucketCounterEnum<ClusterType> mClusterTypeCounts;	
	public BucketCounterEnum<VariantLocation> mVariantLocationCounts;
	public int mCountLOHreferenceLost;
	public int mMinBasePairPosition, mMaxBasePairPosition;
	
	public Gene(String name, Chrom chrom) {		
		mLabel = name;
		mChrom = chrom;
		initializePatients();
				
		mMutationTypeCounts    = new BucketCounterEnum<MutationType>(MutationType.class);  // stores hit counts for synonymous, nonsynonymous
		mClusterTypeCounts     = new BucketCounterEnum<ClusterType>(ClusterType.class);				
		mVariantLocationCounts = new BucketCounterEnum<VariantLocation>(VariantLocation.class);  // stores hit counts for germline, somatic
		clearCounts();
				
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
		mPatients = ArrayUtils.createEnumMapOfArrayLists(ClusterType.class, String.class);
	}
	
	/** Returns true if the patient was already existing for the clustertype, false otherwise. */
	public boolean addPatientIfNotAlreadyAdded(String patientName, ClusterType clusterType) {		
		ArrayList<String> patientList = getPatientsForClusterType(clusterType);		
		return ArrayUtils.checkInListAndInsertIfMissing(patientList, patientName);
	}
	
	public int getNumPatientsForClusterType(ClusterType clusterType) {
		ArrayList<String> patientList = getPatientsForClusterType(clusterType);
		return patientList.size();
	}
	
	public ArrayList<String> getPatientsForClusterType(ClusterType ct) {
		return mPatients.get(ct);
	}
	
	public void incrementCountForClusterType(ClusterType clusterType)    { mClusterTypeCounts.increment(clusterType); }	
	public void incrementCountForMutationType(MutationType mutationType) { mMutationTypeCounts.increment(mutationType); } 			
	public void incrementCountForVariantLocation(VariantLocation varLoc) { mVariantLocationCounts.increment(varLoc); }
	
	public int getCountVariantLocation(VariantLocation varLoc) { return mVariantLocationCounts.getCount(varLoc); }
	public int getCountMutationType(MutationType mutationType) { return mMutationTypeCounts.getCount(mutationType); }
	public int getCountClusterType(ClusterType clusterType)    { return mClusterTypeCounts.getCount(clusterType); }	
	
	public float getRecurrence(ClusterType ct, int total) { return (float) mPatients.get(ct).size() / (float) total; }
	
	public int getRangeLength() { return (mMaxBasePairPosition - mMinBasePairPosition + 1); }
	
	public float getDensityVariantLocation(VariantLocation varLoc) {
		return (float) getCountVariantLocation(varLoc) / (float) getRangeLength();
	}
	
	public float getDensityMutationType(MutationType mutationType) {
		return (float) getCountMutationType(mutationType) / (float) getRangeLength();
	}
	
	public float getDensityClusterType(ClusterType clusterType) {
		return (float) getCountClusterType(clusterType) / (float) getRangeLength(); 
	}
		
	public String countsToString(String delim) {
		StringBuilder sb = new StringBuilder(8192);
		
		sb.append(getCountMutationType(MutationType.NonSynonymous_SNV))
		  .append(delim).append(getCountMutationType(MutationType.Synonymous_SNV))
		  .append(delim).append(getCountVariantLocation(VariantLocation.Germline))
		  .append(delim).append(getCountVariantLocation(VariantLocation.Somatic))
		  
		  .append(delim).append(getCountClusterType(ClusterType.GainSomatic))
		  .append(delim).append(getCountClusterType(ClusterType.LOH))
		  .append(delim).append(mCountLOHreferenceLost)
		  .append(delim).append(getCountClusterType(ClusterType.HETGermline))
		  .append(delim).append(getCountClusterType(ClusterType.HETSomatic))
		  
		  .append(delim).append(Utils.log(getCountClusterType(ClusterType.GainSomatic)))
		  .append(delim).append(Utils.log(getCountClusterType(ClusterType.LOH)))
		  .append(delim).append(Utils.log(mCountLOHreferenceLost))
		  .append(delim).append(Utils.log(getCountClusterType(ClusterType.HETGermline)))
		  .append(delim).append(Utils.log(getCountClusterType(ClusterType.HETSomatic)))
		  
		  .append(delim).append(getDensityClusterType(ClusterType.GainSomatic))
		  .append(delim).append(getDensityClusterType(ClusterType.LOH))
		  .append(delim).append(mCountLOHreferenceLost / (float) getRangeLength())
		  .append(delim).append(getDensityClusterType(ClusterType.HETGermline))
		  .append(delim).append(getDensityClusterType(ClusterType.HETSomatic))
		  
		  .append(delim).append(getNumPatientsForClusterType(ClusterType.GainSomatic))
		  .append(delim).append(getNumPatientsForClusterType(ClusterType.LOH))
		  .append(delim).append(getNumPatientsForClusterType(ClusterType.HETGermline))
		  .append(delim).append(getNumPatientsForClusterType(ClusterType.HETSomatic));
		
		return sb.toString();
	}
	
	public int compareTo(Gene rhs) {
		return mLabel.compareToIgnoreCase(rhs.mLabel);
	}
}