import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;

import lohcateEnums.Chrom;
import lohcateEnums.ClusterType;
import lohcateEnums.SNVType;
import lohcateEnums.VariantLocation;

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
	public ArrayList<ArrayList<String>> mPatients;
	
	public int[] mSNVTypeCounts;
	public int[] mClusterTypeCounts;	
	public int[] mVariantLocationCounts;
	public int   mCountLOHreferenceLost;
	public int mMinBasePairPosition, mMaxBasePairPosition;
	
	public Gene(String name, Chrom chrom) {		
		mLabel = name;
		mChrom = chrom;
		initializePatients();
				
		mSNVTypeCounts         = new int[        SNVType.values().length];  // stores hit counts for synonymous, nonsynonymous
		mClusterTypeCounts     = new int[    ClusterType.values().length];  //stores hit counts for LOH, dup, &c. variants
		mVariantLocationCounts = new int[VariantLocation.values().length];  // stores hit counts for germline, somatic
		clearCounts();
				
		mMinBasePairPosition = Integer.MAX_VALUE;
		mMaxBasePairPosition = Integer.MIN_VALUE;
	}
	
	public void clearCounts() {
		Arrays.fill(mSNVTypeCounts,         0);
		Arrays.fill(mClusterTypeCounts,     0);
		Arrays.fill(mVariantLocationCounts, 0);
		mCountLOHreferenceLost = 0;
	}
	
	private void initializePatients() {
		mPatients = new ArrayList<ArrayList<String>>();
		Utils.addNewEmptyArrayLists(mPatients, ClusterType.values().length);
	}
	
	/** Returns true if the patient was already existing for the clustertype, false otherwise. */
	public boolean addPatientIfNotAlreadyAdded(String patientName, ClusterType clusterType) {
		ArrayList<String> patientList = mPatients.get(clusterType.ordinal());
		int resultIndex = Collections.binarySearch(patientList, patientName);
		if (resultIndex < 0) {
			int insertionIndex = -(resultIndex + 1);
			patientList.add(insertionIndex, patientName);
			return false;
		} else {
			return true;
		}		
	}
	
	public int getNumPatientsForClusterType(ClusterType clusterType) {
		ArrayList<String> patientList = mPatients.get(clusterType.ordinal());
		return patientList.size();
	}
	
	public void incrementCountForClusterType(ClusterType clusterType) {  mClusterTypeCounts[clusterType.ordinal()]++; }
	
	public void incrementCountForMutationType(SNVType mutationType) { mSNVTypeCounts[mutationType.ordinal()]++; }
	
	public void incrementCountForVariantLocation(VariantLocation varLoc) { mVariantLocationCounts[varLoc.ordinal()]++; }
	
	public int getCountVariantLocation(VariantLocation varLoc) { return mVariantLocationCounts[varLoc.ordinal()]; }
	public int getCountMutationType(SNVType mutationType)      { return mSNVTypeCounts[mutationType.ordinal()]; }
	public int getCountClusterType(ClusterType clusterType)    { return mClusterTypeCounts[clusterType.ordinal()]; }
	
	
	public float getRecurrence(int ind, int total) { return (float)mPatients.get(ind).size() / (float)total; }
	
	public int getRangeLength() { return (mMaxBasePairPosition - mMinBasePairPosition + 1); }
	
	public float getDensityVariantLocation(VariantLocation varLoc) {
		return (float) getCountVariantLocation(varLoc) / (float) getRangeLength();
	}
	
	public float getDensityMutationType(SNVType mutationType) {
		return (float) getCountMutationType(mutationType) / (float) getRangeLength();
	}
	
	public float getDensityClusterType(ClusterType clusterType) {
		return (float) getCountClusterType(clusterType) / (float) getRangeLength(); 
	}
		
	public String countsToString(String delim) {
		StringBuilder sb = new StringBuilder(8192);
		
		sb.append(getCountMutationType(SNVType.NonSynonymous))
		  .append(delim).append(getCountMutationType(SNVType.Synonymous))
		  .append(delim).append(getCountVariantLocation(VariantLocation.Germline))
		  .append(delim).append(getCountVariantLocation(VariantLocation.Somatic))
		  
		  .append(delim).append(getCountClusterType(ClusterType.Dup))
		  .append(delim).append(getCountClusterType(ClusterType.LOH))
		  .append(delim).append(mCountLOHreferenceLost)
		  .append(delim).append(getCountClusterType(ClusterType.HET))
		  
		  .append(delim).append(Utils.log(getCountClusterType(ClusterType.Dup)))
		  .append(delim).append(Utils.log(getCountClusterType(ClusterType.LOH)))
		  .append(delim).append(Utils.log(mCountLOHreferenceLost))
		  .append(delim).append(Utils.log(getCountClusterType(ClusterType.HET)))
		  
		  .append(delim).append(getDensityClusterType(ClusterType.Dup))
		  .append(delim).append(getDensityClusterType(ClusterType.LOH))
		  .append(delim).append(mCountLOHreferenceLost / (float) getRangeLength())
		  .append(delim).append(getDensityClusterType(ClusterType.HET))
		  
		  .append(delim).append(getNumPatientsForClusterType(ClusterType.Dup))
		  .append(delim).append(getNumPatientsForClusterType(ClusterType.LOH))
		  .append(delim).append(getNumPatientsForClusterType(ClusterType.HET));
		
		return sb.toString();
	}
	
	public int compareTo(Gene rhs) {
		return mLabel.compareToIgnoreCase(rhs.mLabel);
	}
}