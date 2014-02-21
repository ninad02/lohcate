package lohcate.geneAnalysis;

import genomeEnums.Chrom;
import genomeEnums.VariantFrequency;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;

import com.google.common.collect.ArrayListMultimap;

import nutils.ArrayUtils;
import nutils.CompareUtils;
import nutils.StringUtils;

import lohcate.EventTypeAllele;
import lohcate.Regions;
import lohcateEnums.EventType;
import lohcateEnums.MutationType;

/** 
 * A wrapper class that holds tables, etc. when tallying over genes from multiple patients
 * @author Ninad Dewal
 *
 */
public class GenesTallier {

	// ========================================================================
	// MEMBER VARIABLES
	// ________________________________________________________________________
	
	// A list of genes that we will keep binary sorted for efficieny purposes
	public ArrayList<GeneCounter> mGenes;	
	protected ArrayList<GeneCounter> mAllGeneCountersUnsorted;  // an unsorted array to map numeric geneIDs to genes 	
	
	// ========================================================================	
	public GenesTallier() {
		mGenes                   = new ArrayList<GeneCounter>();
		mAllGeneCountersUnsorted = new ArrayList<GeneCounter>();		
	}

	// ========================================================================
	public ArrayList<GeneCounter> getGenes() { return mGenes; }
	
	// ========================================================================
	public GeneCounter getGene(int geneID) { return mAllGeneCountersUnsorted.get(geneID); }
	
	// ========================================================================
	public void performBasicGeneTallyingForPatient(
		String patientName,
		ArrayList<String> allLines,
		boolean hasHeader,
		ArrayListMultimap<GeneCounter, EventType> eventsForGene) {
		
		// Create a dummy gene
		GeneCounter dummyGene = new GeneCounter("", Chrom.c0, 0);  // We'll use this gene for binary searching
		
		// Clear out the hashmap to be filled
		eventsForGene.clear();
		
		// Set the starting row depending on whether there's a header line or not
		int startingRow = hasHeader ? 1 : 0;  
		
		// Loop through the rows
		for (int row = startingRow; row < allLines.size(); row++) {

			//ClusteringInputOneSite oneSiteInfo = new ClusteringInputOneSite(allLines.get(row), SeqPlatform.Illumina);					

			String components[] = allLines.get(row).split(StringUtils.TabStr); 
			String geneName =               components[Regions.ColCuratedTSV_Gene];
			Chrom chrom  = Chrom.getChrom(  components[Regions.ColCuratedTSV_Chrom]);
			int position = Integer.parseInt(components[Regions.ColCuratedTSV_Position]);					
			String[] dbsnpStrSplit = components[Regions.ColCuratedTSV_dbSNP].split(";");
			VariantFrequency varFreq = VariantFrequency.valueOf(dbsnpStrSplit[1]);					

			// Avoid invalid gene names (e.g. "", ".", or whitespace				
			if (geneName.equals(Regions.MissingGeneNameValue) || geneName.trim().isEmpty()) { 
				continue;
			}

			//gene names can sometimes come with an uninteresting/irrelevant prefix
			if (geneName.indexOf("dist") >= 0) { 
				geneName = geneName.split("\\(")[0];
			}

			// Set the dummy gene for binary search
			dummyGene.mLabel = geneName;
			dummyGene.mChrom = chrom;

			// Find the gene.
			int resultIndex = Collections.binarySearch(mGenes, dummyGene);
			if (resultIndex < 0) {  
				resultIndex = ArrayUtils.getInsertPoint(resultIndex);  // calculate the proper insertion point
				GeneCounter geneCounter = new GeneCounter(geneName, chrom, mAllGeneCountersUnsorted.size());
				mAllGeneCountersUnsorted.add(geneCounter);
				mGenes.add(resultIndex, geneCounter);							
			}
			GeneCounter currentGene = mGenes.get(resultIndex);						

			currentGene.mMaxBasePairPosition = Math.max(position, currentGene.mMaxBasePairPosition); // Set right-bound
			currentGene.mMinBasePairPosition = Math.min(position, currentGene.mMinBasePairPosition); // Set left-bound

			MutationType mutationType = MutationType.getSNVType(components[Regions.ColCuratedTSV_MutationType]);
			currentGene.incrementCount(mutationType);  // increment synonymous or nonsynonymous variant count						

			//VariantLocation variantLocation = 
			//		Clustering.isVariantInGermline(Double.parseDouble(components[ColCuratedTSV_VafNormal])) ? VariantLocation.Germline : VariantLocation.Somatic; 						
			//currentGene.incrementCount(variantLocation);

			EventType eventType = EventType.getClusterType(components[Regions.ColCuratedTSV_Cluster]);
			if (eventType == null) { CompareUtils.throwErrorAndExit("ERROR: Invalid cluster type: " + components[Regions.ColCuratedTSV_Cluster]); }
			currentGene.incrementCount(eventType); 

			// Check if the gene-event combo already exists for this
			if (!GeneEnrichment.isEventToIgnore(eventType)) {
				List<EventType> values = eventsForGene.get(currentGene);
				if (values.isEmpty()) {
					eventsForGene.put(currentGene, eventType);
				} else {
					if (!values.contains(eventType)) {
						eventsForGene.put(currentGene, eventType);
					}
				}				
			}

			double vafTumor = Double.parseDouble(components[Regions.ColCuratedTSV_VafTumor]);
			tallyAlleleSpecificBreakdown(currentGene, eventType, vafTumor, varFreq, patientName, position);
			currentGene.addPatientIfNotAlreadyAdded(patientName, eventType);
		}		
	}
	
	// ========================================================================
	private static void tallyAlleleSpecificBreakdown(GeneCounter currentGene, EventType eventType, double vafTumor, VariantFrequency varFreq, String patientName, int position) {
		boolean vafTumorIsEnriched = (vafTumor > 0.5);
		boolean isLOHType = (eventType == EventType.LOH || eventType == EventType.cnLOH);
		boolean isAmpType = (eventType == EventType.GainSomatic);
		boolean varAlleleNotCommon = (varFreq == VariantFrequency.Rare || varFreq == VariantFrequency.Novel);
		if (vafTumorIsEnriched) {
			if (isLOHType) {
				if (varAlleleNotCommon) {
					currentGene.mCountLOHreferenceLost++;
					currentGene.addPatientIfNotAlreadyAdded_byAlleleEvent(EventTypeAllele.LossReference_VariantRare, patientName, position);
				} else {
					currentGene.addPatientIfNotAlreadyAdded_byAlleleEvent(EventTypeAllele.LossReference_VariantCommon, patientName, position);
				}
			} else if (isAmpType) {
				if (varAlleleNotCommon) {
					currentGene.addPatientIfNotAlreadyAdded_byAlleleEvent(EventTypeAllele.GainVariant_VariantRare, patientName, position);
				} else {
					currentGene.addPatientIfNotAlreadyAdded_byAlleleEvent(EventTypeAllele.GainVariant_VariantCommon, patientName, position);
				}								
			} 
		} else {
			// Variant allele was lost
			if (isLOHType) {
				if (varAlleleNotCommon) {
					currentGene.addPatientIfNotAlreadyAdded_byAlleleEvent(EventTypeAllele.LossVariant_VariantRare, patientName, position);
				} else {
					currentGene.addPatientIfNotAlreadyAdded_byAlleleEvent(EventTypeAllele.LossVariant_VariantCommon, patientName, position);
				}
			} else if (isAmpType) {
				if (varAlleleNotCommon) {
					currentGene.addPatientIfNotAlreadyAdded_byAlleleEvent(EventTypeAllele.GainReference_VariantRare, patientName, position);
				} else {
					currentGene.addPatientIfNotAlreadyAdded_byAlleleEvent(EventTypeAllele.GainReference_VariantCommon, patientName, position);
				}
			}							
		}
	}
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub

	}

}
