package lohcate;

import genomeEnums.Chrom;
import genomeEnums.VariantFrequency;

import java.io.BufferedWriter;
import java.io.File;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.Collections;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;

import lohcate.clustering.Clustering;
import lohcateEnums.EventType;
import lohcateEnums.MutationType;
import nutils.ArrayUtils;
import nutils.Cast;
import nutils.CompareUtils;
import nutils.EnumMatrix;
import nutils.IOUtils;
import nutils.NullaryClassFactory;
import nutils.NumberUtils;
import nutils.PrimitiveWrapper;
import nutils.StringUtils;
import nutils.BitUtils.Compactor.CompactorInf;
import nutils.BitUtils.Compactor.CompactorIntoInt;
import nutils.BitUtils.Compactor.CompactorIntoLong;
import nutils.StringUtils.FileExtensionAndDelimiter;
import nutils.collectionsSorted.ArrayListSortedComparable;
import nutils.counter.BucketCounterEnum;

public class GeneEnrichment {

	public GeneEnrichment() {
		// TODO Auto-generated constructor stub
	}

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		TestCompactUnitRobust();

	}
	
	private static void TestCompactUnitRobust() {
		int numIter = 100000000;
		int milestone = numIter / 10;
		for (int i = 0; i < numIter; i++) {
			if (i % milestone == 0) System.out.println("Reached iteration: " + i);
			int randomGeneID1 = NumberUtils.getRandomInteger(0, 65535);
			int randomGeneID2;
			do {
				randomGeneID2 = NumberUtils.getRandomInteger(0, 65535);
			} while (randomGeneID2 == randomGeneID1);
			
			EventType e1 = EventType.values()[NumberUtils.getRandomInteger(0, EventType.values().length - 1)];
			EventType e2 = EventType.values()[NumberUtils.getRandomInteger(0, EventType.values().length - 1)];
			long unit1 = compactIntoUnit(randomGeneID1, randomGeneID2, e1, e2);
			long unit2 = compactIntoUnit(randomGeneID2, randomGeneID1, e2, e1);
			if (unit1 != unit2) {
				System.err.println("ERROR: Different compact units: \t" + unit1 + "\t" + unit2 + "\t" + randomGeneID1 + "\t" + randomGeneID2 + "\t" + e1 + "\t" + e2);
			}
		}
	}
	
	private static boolean isEventToIgnore(EventType eventType) {
		return (eventType == EventType.Noise) 
		    || (eventType == EventType.Ignored)
		    || (eventType == EventType.HETGermline)
		; 
	}

	/**
	 * Generate 'master' gene enrichment table (used to generate histograms).
	 * @param inDir curated SNP calls
	 */
	public static void getGeneEnrichment(String inDir, String outDir) {
		String outFilename = outDir + File.separator + "geneEnrichment" + StringUtils.FileExtensionTSV.mExtension;
		IOUtils.createDirectoryPath(outDir, false);
		File[] files = (new File(inDir)).listFiles();
		StringUtils.FileExtensionAndDelimiter fileExtDelim = StringUtils.FileTextTabDelim; 
				
		StringBuilder sb = new StringBuilder(4096);		
		
		// A list of genes that we will keep binary sorted for efficieny purposes
		ArrayList<GeneCounter> genes = new ArrayList<GeneCounter>(); 
		GeneCounter dummyGene = new GeneCounter("", Chrom.c0, 0);  // We'll use this gene for binary searching
		
		HashMap<GeneCounter, EventType> eventForGene = new HashMap<GeneCounter, EventType>();
				
		//TwoGenesAndEvents2 coOcurrenceWithPatientsDummy = new TwoGenesAndEvents2();
		
		HashMap<PrimitiveWrapper.WLong, ArrayList<String>> coOcurrenceWithPatientsCompact = 
				new HashMap<PrimitiveWrapper.WLong, ArrayList<String>>();
		PrimitiveWrapper.WLong twoGenesTwoEventsCompactDummy = new PrimitiveWrapper.WLong(0);
		HashMap<PrimitiveWrapper.WLong, BitSet> coOcurrenceWithPatients = new HashMap<PrimitiveWrapper.WLong, BitSet>();
		
		NullaryClassFactory<DualEventsWithPatientsList> listClassFactory = new NullaryClassFactory<>(DualEventsWithPatientsList.class);
		EnumMatrix<EventType, DualEventsWithPatientsList> countMatrix = new EnumMatrix<>(EventType.class, listClassFactory);
		DualEventsWithPatients dualEventDummy = new DualEventsWithPatients();
				
		int geneCounterID = 0;
		ArrayList<GeneCounter> allGeneCountersUnsorted = new ArrayList<GeneCounter>();  // keep unsorted to map to index
		ArrayList<String> allPatients = new ArrayList<String>();
		
		for (File file : files) { //iterate through curated SNP calls
			int suffixIndex = -1;
			
			if (file.getName().endsWith(fileExtDelim.mExtension) 
			    && ((suffixIndex = file.getName().indexOf(Clustering.SitesClassifiedOutputSuffix)) >= 0)) {
				String patientName = file.getName().substring(0, suffixIndex);
				allPatients.add(patientName);
				int patientIndex = allPatients.size() - 1;
				System.out.println(allPatients.size() + "\t" + file.getName() + "\t" + patientName);
				
				// Load all lines into memory, extract the header row, and then sort by chrom/position
				ArrayList<String> allLines = IOUtils.readAllLinesFromFile(file.getAbsolutePath(), false, true, sb);
				String headerString = sb.toString();
				Collections.sort(allLines, Regions.LineComparatorTab);					
		
				eventForGene.clear();
				int startingRow = 0;  // because the header was stripped away
				for (int row = startingRow; row < allLines.size(); row++) {
				
					//if (row % 10000 == 0) System.out.println(row + " of " + allLines.size());
					
					//ClusteringInputOneSite oneSiteInfo = new ClusteringInputOneSite(allLines.get(row), SeqPlatform.Illumina);					
										
					
					String components[] = allLines.get(row).split(StringUtils.TabStr); 
					String geneName =               components[Regions.ColCuratedTSV_Gene];
					Chrom chrom  = Chrom.getChrom(  components[Regions.ColCuratedTSV_Chrom]);
					int position = Integer.parseInt(components[Regions.ColCuratedTSV_Position]);					
					String[] dbsnpStrSplit = components[Regions.ColCuratedTSV_dbSNP].split(";");
					VariantFrequency varFreq = VariantFrequency.valueOf(dbsnpStrSplit[1]);					
					
					if (!geneName.equals(Regions.MissingGeneNameValue)) { //avoid ".", which crops up a lot
						if (geneName.indexOf("dist") >= 0) { //gene names can sometimes come with an uninteresting/irrelevant prefix
							geneName = geneName.split("\\(")[0];
						}
					
						// Set the dummy gene for binary search
						dummyGene.mLabel = geneName;
						dummyGene.mChrom = chrom;
						
						int resultIndex = Collections.binarySearch(genes, dummyGene);
						if (resultIndex < 0) {  
							// we haven't seen this gene before
							resultIndex = ArrayUtils.getInsertPoint(resultIndex);  // calculate the proper insertion point
							GeneCounter geneCounter = new GeneCounter(geneName, chrom, geneCounterID++);
							allGeneCountersUnsorted.add(geneCounter);
							genes.add(resultIndex, geneCounter);							
						}
						GeneCounter currentGene = genes.get(resultIndex);						
						
						if (position > currentGene.mMaxBasePairPosition) { //get right-bound of gene's range of variants
							currentGene.mMaxBasePairPosition = position; 
						}
						if (position < currentGene.mMinBasePairPosition) { //...left-bound...
							currentGene.mMinBasePairPosition = position;
						}
												
						MutationType mutationType = MutationType.getSNVType(components[Regions.ColCuratedTSV_MutationType]);
						currentGene.incrementCount(mutationType);  // increment synonymous or nonsynonymous variant count						
						
						//VariantLocation variantLocation = 
						//		Clustering.isVariantInGermline(Double.parseDouble(components[ColCuratedTSV_VafNormal])) ? VariantLocation.Germline : VariantLocation.Somatic; 						
						//currentGene.incrementCount(variantLocation);
																
						EventType eventType = EventType.getClusterType(components[Regions.ColCuratedTSV_Cluster]);
						if (eventType == null) { CompareUtils.throwErrorAndExit("ERROR: Invalid cluster type: " + components[Regions.ColCuratedTSV_Cluster]); }
						
						// Check if the gene-event combo already exists for this
						if (!isEventToIgnore(eventType) && eventType == EventType.LOH) {
							EventType tableValue = eventForGene.get(currentGene);
							if (tableValue == null) {
								eventForGene.put(currentGene, eventType);
							} else if (eventType != tableValue) {
								//CompareUtils.throwErrorAndExit("ERROR: Different events for gene: \t" + tableValue + "\t" + eventType);
								System.out.println("WARNING: Different events for gene:\t" + geneName + "\t" + tableValue + "\t" + eventType);
							}
						}
						
						currentGene.incrementCount(eventType); 
						double vafTumor = Double.parseDouble(components[Regions.ColCuratedTSV_VafTumor]);
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
		
						currentGene.addPatientIfNotAlreadyAdded(file.getName(), eventType);
					}
				}
				
				// Now that we've gotten all the gene-event combinations for a particular
				// sample, we store all the combinations.  We do an O(n^2) operation
				System.out.println("Num Events for Sample: " + eventForGene.size());
				int numCombos = 0;			
				int numNewCombos = 0;
				Set<Map.Entry<GeneCounter, EventType>> keysAndValuesMap = eventForGene.entrySet();
				for (Map.Entry<GeneCounter, EventType> keyValueIter1 : keysAndValuesMap) {
					for (Map.Entry<GeneCounter, EventType> keyValueIter2: keysAndValuesMap) {
						GeneCounter gene1 = keyValueIter1.getKey();
						GeneCounter gene2 = keyValueIter2.getKey();
						EventType event1 = keyValueIter1.getValue();
						EventType event2 = keyValueIter2.getValue();
						
						if ( (keyValueIter1.getValue() == keyValueIter2.getValue()) && (gene1.compareTo(gene2) == 0) ) {
							// No need to continue the loop since combinations will be checked anyway with outer loop eventually
							break;
						} else {				
							/*
							twoGenesTwoEventsCompactDummy.mLong = compactIntoUnit(gene1.getNumericID(), gene2.getNumericID(), keyValueIter1.getValue(), keyValueIter2.getValue());
							numCombos++;
							BitSet patients = coOcurrenceWithPatients.get(twoGenesTwoEventsCompactDummy);
							if (patients == null) {
								patients = new BitSet(256); 								
								//coOcurrenceWithPatients.put(twoGenesTwoEventsCompactDummy, patients);
								twoGenesTwoEventsCompactDummy = new PrimitiveWrapper.WLong(0);
								numNewCombos++;
							}
							patients.set(patientIndex);
							*/
							
							// Swap if necessary
							if (gene1.getNumericID() > gene2.getNumericID()) {
								GeneCounter geneTemp = gene1;
								gene1 = gene2;
								gene2 = geneTemp;
								
								EventType eventTemp = event1;
								event1 = event2;
								event2 = eventTemp;
							}
							numCombos++;
							dualEventDummy.setEvents(gene1.getNumericID(), gene2.getNumericID(), event1, event2);
							DualEventsWithPatientsList dualEventList = countMatrix.get(event1, event2);
							DualEventsWithPatients dualEvent = dualEventList.mList.get(dualEventDummy);
							if (dualEvent == null) {
								dualEventDummy.mPatients.clear();
								dualEventDummy.mPatients.set(patientIndex);
								dualEventList.mList.add(dualEventDummy);
								dualEventDummy = new DualEventsWithPatients();		
								numNewCombos++;
							} else {
								dualEvent.mPatients.set(patientIndex);
							}							
						}					
					}
				}
				System.out.println("Num Combos: " + numCombos);
				System.out.println("Num New Combos: " + numNewCombos);
				eventForGene.clear();
				System.gc();
				System.gc();
			}
		}
		
		String[] columnHeaders = new String[] { 
				"chr", "bp_start", "bp_end", "length", "gene", GeneCounter.getHeaders(fileExtDelim.mDelimiter)
		};
		String headerStr = StringUtils.constructColumnDelimitedString(columnHeaders, fileExtDelim.mDelimiter, sb, true).toString();
	
		BufferedWriter out = IOUtils.getBufferedWriter(outFilename);
		IOUtils.writeToBufferedWriter(out, headerStr, true);
		for (GeneCounter gene : genes) {			
			sb.setLength(0);
			sb.append(gene.mChrom.getCode())
				.append(fileExtDelim.mDelimiter).append(gene.mMinBasePairPosition)
				.append(fileExtDelim.mDelimiter).append(gene.mMaxBasePairPosition)
				.append(fileExtDelim.mDelimiter).append(gene.getRangeLength())
				.append(fileExtDelim.mDelimiter).append(gene.mLabel)
				.append(fileExtDelim.mDelimiter).append(gene.countsToString(fileExtDelim.mDelimiter));
			IOUtils.writeToBufferedWriter(out, sb.toString(), true);
			IOUtils.flushBufferedWriter(out);
		}
		IOUtils.closeBufferedWriter(out);
		
		String eventsByPatientPerGeneFilename = outDir + File.separator + "eventsByPatientPerGene" + StringUtils.FileExtensionTSV.mExtension;
		BufferedWriter outBreakdown = IOUtils.getBufferedWriter(eventsByPatientPerGeneFilename);
				
		// Now write out individual gene outputs
		for (GeneCounter gene : genes) {
			// Write samples for each gene out as well			
			for (EventType eventType : EventType.values()) {
				if (isEventToIgnore(eventType)) continue;
				
				ArrayList<String> patientsForEvent = gene.getPatientsForEventType(eventType);
				for (String patientForEvent : patientsForEvent) {
					String outString = gene.getName() + "\t" + eventType.name() + "\t" + patientForEvent;
					IOUtils.writeToBufferedWriter(outBreakdown, outString, true);
				}
			}
		}
		IOUtils.closeBufferedWriter(outBreakdown);
		
		// Store co-occurrence counts	
		
		String coOccurrenceFilename = outDir + File.separator + "coOccurrence" + StringUtils.FileExtensionTSV.mExtension;
		BufferedWriter outCoOccurrence = IOUtils.getBufferedWriter(coOccurrenceFilename);		
		for (Map.Entry<PrimitiveWrapper.WLong, BitSet> entry : coOcurrenceWithPatients.entrySet()) {
			
			PrimitiveWrapper.WLong compactUnit = entry.getKey();			
			int geneID1 = Cast.toInt(TwoGenesWithEvents.Compactor.getValue(TwoGenesWithEvents.ID_Gene1, compactUnit.mLong));
			int geneID2 = Cast.toInt(TwoGenesWithEvents.Compactor.getValue(TwoGenesWithEvents.ID_Gene2, compactUnit.mLong));
			GeneCounter gene1 = allGeneCountersUnsorted.get(geneID1);
			GeneCounter gene2 = allGeneCountersUnsorted.get(geneID2);
			
			int eventIndex1 = Cast.toInt(TwoGenesWithEvents.Compactor.getValue(TwoGenesWithEvents.Event_Gene1, compactUnit.mLong));
			EventType eventGene1 = EventType.values()[eventIndex1];
			
			int eventIndex2 = Cast.toInt(TwoGenesWithEvents.Compactor.getValue(TwoGenesWithEvents.Event_Gene2, compactUnit.mLong));
			EventType eventGene2 = EventType.values()[eventIndex2];
			
			BitSet patients = entry.getValue();
			sb.setLength(0);
			sb.append(gene1.getName())
			  .append(StringUtils.FileTextTabDelim.mDelimiter)
			  .append(gene2.getName())
			  .append(StringUtils.FileTextTabDelim.mDelimiter)
			  .append(eventGene1.name())
			  .append(StringUtils.FileTextTabDelim.mDelimiter)			  
			  .append(eventGene2.name())
			  .append(StringUtils.FileTextTabDelim.mDelimiter)
			  .append(patients.size())
			  .append(StringUtils.FileTextTabDelim.mDelimiter)
			  .append(patients.toString());			  
			IOUtils.writeToBufferedWriter(outCoOccurrence, sb.toString(), true);
		}
		
		IOUtils.closeBufferedWriter(outCoOccurrence);
	}

	// ========================================================================
	private static long compactIntoUnit(int numIDGene1, int numIDGene2, EventType eventGene1, EventType eventGene2) {		
		if (numIDGene1 > numIDGene2) {
			return compactIntoUnit(numIDGene2, numIDGene1, eventGene2, eventGene1);  // Call recursively
		} else {
			long unit = 0;
			unit = TwoGenesWithEvents.Compactor.setValue(TwoGenesWithEvents.ID_Gene1, numIDGene1, unit);
			unit = TwoGenesWithEvents.Compactor.setValue(TwoGenesWithEvents.ID_Gene2, numIDGene2, unit);
			unit = TwoGenesWithEvents.Compactor.setValue(TwoGenesWithEvents.Event_Gene1, eventGene1.ordinal(), unit);
			unit = TwoGenesWithEvents.Compactor.setValue(TwoGenesWithEvents.Event_Gene2, eventGene2.ordinal(), unit);
			return unit;
		}		
	}

	// ========================================================================
	public static class DualEventsWithPatientsList {
		protected ArrayListSortedComparable<DualEventsWithPatients> mList;
		public DualEventsWithPatientsList() {
			mList = new ArrayListSortedComparable<>(1000);
		}
	}
	
	// ========================================================================
	public static class DualEventsWithPatients implements Comparable<DualEventsWithPatients> {
		protected long mDualEventCompact;
		protected BitSet mPatients;
		
		public DualEventsWithPatients() {
			mDualEventCompact = 0;
			mPatients = new BitSet();
		}
		
		public void setEvents(int numIDGene1, int numIDGene2, EventType eventGene1, EventType eventGene2) {
			mDualEventCompact = compactIntoUnit(numIDGene1, numIDGene2, eventGene1, eventGene2);
		}

		@Override
		public int compareTo(DualEventsWithPatients o) {			
			return Long.compare(this.mDualEventCompact, o.mDualEventCompact);
		}
		
		
	}
	
	// ========================================================================
	public static enum TwoGenesWithEvents implements CompactorInf<TwoGenesWithEvents> {
		ID_Gene1(16),
		ID_Gene2(16),
		Event_Gene1(4),
		Event_Gene2(4)
		;

		public static CompactorIntoLong<TwoGenesWithEvents> Compactor = new CompactorIntoLong<TwoGenesWithEvents>(TwoGenesWithEvents.class, false);
		
		private short mNumBits;
		
		private TwoGenesWithEvents(int numBits) {
			mNumBits = Cast.toShort(numBits);
		}

		@Override
		public int getNumBits() { return mNumBits; }	
	}
	
	// ========================================================================
	public static class TwoGenesAndEvents2 implements Comparable<TwoGenesAndEvents2> {
		private String mGeneLabel1, mGeneLabel2;
		private EventType mEventForGene1, mEventForGene2;
		
		public TwoGenesAndEvents2() { clear(); }
		
				
		/* (non-Javadoc)
		 * @see java.lang.Object#hashCode()
		 */
		@Override
		public int hashCode() {
			final int prime = 31;
			int result = 1;
			result = prime
					* result
					+ ((mEventForGene1 == null) ? 0 : mEventForGene1.hashCode());
			result = prime
					* result
					+ ((mEventForGene2 == null) ? 0 : mEventForGene2.hashCode());
			result = prime * result
					+ ((mGeneLabel1 == null) ? 0 : mGeneLabel1.hashCode());
			result = prime * result
					+ ((mGeneLabel2 == null) ? 0 : mGeneLabel2.hashCode());
			return result;
		}



		/* (non-Javadoc)
		 * @see java.lang.Object#equals(java.lang.Object)
		 */
		@Override
		public boolean equals(Object obj) {
			if (this == obj)
				return true;
			if (obj == null)
				return false;
			if (getClass() != obj.getClass())
				return false;
			TwoGenesAndEvents2 other = (TwoGenesAndEvents2) obj;
			if (mEventForGene1 != other.mEventForGene1)
				return false;
			if (mEventForGene2 != other.mEventForGene2)
				return false;
			if (mGeneLabel1 == null) {
				if (other.mGeneLabel1 != null)
					return false;
			} else if (!mGeneLabel1.equals(other.mGeneLabel1))
				return false;
			if (mGeneLabel2 == null) {
				if (other.mGeneLabel2 != null)
					return false;
			} else if (!mGeneLabel2.equals(other.mGeneLabel2))
				return false;
			return true;
		}

		public void set(String geneLabel1, String geneLabel2, EventType eventForGene1, EventType eventForGene2) {			
			if (geneLabel1.compareToIgnoreCase(geneLabel2) > 0) {
				setInternal(geneLabel2, geneLabel1, eventForGene2, eventForGene1);
			} else {
				setInternal(geneLabel1, geneLabel2, eventForGene1, eventForGene2);
			}
		}
		
		private void setInternal(String geneLabel1, String geneLabel2, EventType eventForGene1, EventType eventForGene2) {			
			mGeneLabel1 = geneLabel1;
			mGeneLabel2 = geneLabel2;
			mEventForGene1 = eventForGene1;
			mEventForGene2 = eventForGene2;			
		}

		public void clear() {
			mGeneLabel1 = mGeneLabel2 = "";
			mEventForGene1 = mEventForGene2 = null;
		}
		
		public String toString() {
			StringBuilder sb = new StringBuilder(2048);
			sb.append(mGeneLabel1).append(StringUtils.FileTextTabDelim.mDelimiter)			  
			  .append(mEventForGene1).append(StringUtils.FileTextTabDelim.mDelimiter)
			  .append(mGeneLabel2).append(StringUtils.FileTextTabDelim.mDelimiter)
			  .append(mEventForGene2).append(StringUtils.FileTextTabDelim.mDelimiter);
			return sb.toString();
		}

		@Override
		public int compareTo(TwoGenesAndEvents2 rhs) {
			// Do the fast checks first (by event)
			int result = compareTo_Helper(this.mEventForGene1, rhs.mEventForGene1);
			if (result == 0) {
				result = compareTo_Helper(this.mEventForGene2, rhs.mEventForGene2);
				if (result == 0) {
					result = this.mGeneLabel1.compareToIgnoreCase(rhs.mGeneLabel1);
					if (result == 0) {
						result = this.mGeneLabel2.compareToIgnoreCase(rhs.mGeneLabel2);
					}
				}
			}
			return result;
		}
		
		private int compareTo_Helper(EventType eventForGene1, EventType eventForGene2) {
			if (eventForGene1 == eventForGene2) {
				return 0;
			} else if (eventForGene1 == null && eventForGene2 != null) {
				return -1;
			} else if (eventForGene1 != null && eventForGene2 == null) {
				return 1;
			} else {
				return eventForGene1.compareTo(eventForGene2);
			}
		}
		
		
	}
}
