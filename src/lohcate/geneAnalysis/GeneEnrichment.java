package lohcate.geneAnalysis;

import genomeEnums.Chrom;
import genomeEnums.VariantFrequency;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.Collection;
import java.util.Collections;
import java.util.Date;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.carrotsearch.hppc.DoubleArrayList;
import com.carrotsearch.hppc.IntArrayList;
import com.carrotsearch.hppc.LongArrayList;
import com.google.common.collect.ArrayListMultimap;

import lohcate.EventTypeAllele;
import lohcate.Regions;
import lohcate.clustering.AlleleFractionStatsForSample;
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
import nutils.BitUtils.BitSetUtils;
import nutils.BitUtils.ValueExtractor;
import nutils.BitUtils.Compactor.CompactorInf;
import nutils.BitUtils.Compactor.CompactorIntoInt;
import nutils.BitUtils.Compactor.CompactorIntoLong;
import nutils.StringUtils.FileExtensionAndDelimiter;
import nutils.array.FastArrayListLong;
import nutils.array.FastArrayListLong.FastArrayListIndex;
import nutils.array.FastArrayListLong.FastArrayListIterator;
import nutils.collectionsSorted.ArrayListSortedComparable;
import nutils.counter.BucketCounterEnum;

public class GeneEnrichment {

	private static final int InitialDualEventCount = 1;
	
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
			long unit1 = compactIntoUnitLong(randomGeneID1, randomGeneID2, e1, e2, Cast.toShort(InitialDualEventCount));
			long unit2 = compactIntoUnitLong(randomGeneID2, randomGeneID1, e2, e1, Cast.toShort(InitialDualEventCount));
			if (unit1 != unit2) {
				System.err.println("ERROR: Different compact units: \t" + unit1 + "\t" + unit2 + "\t" + randomGeneID1 + "\t" + randomGeneID2 + "\t" + e1 + "\t" + e2);
			}
		}
	}
	
	// ========================================================================
	protected static boolean isEventToIgnore(EventType eventType) {
		return (eventType == EventType.Noise) 
		    || (eventType == EventType.Ignored)
		    || (eventType == EventType.HETGermline)
		; 
	}

	// ========================================================================
	/**
	 * Generate 'master' gene enrichment table (used to generate histograms).
	 * @param inDir curated SNP calls
	 */
	public static void getGeneEnrichment(String inDir, String outDir) {
		IOUtils.createDirectoryPath(outDir, false);
		File[] files = (new File(inDir)).listFiles();
		StringUtils.FileExtensionAndDelimiter fileExtDelim = StringUtils.FileTextTabDelim;
		
		String outFilename                  = outDir + File.separator + "geneEnrichment"    + StringUtils.FileExtensionTSV.mExtension;		
		String outFilenameDualEventsOneGene = outDir + File.separator + "dualEventsOneGene" + StringUtils.FileExtensionTSV.mExtension;
		BufferedWriter outDualEventsOneGene = IOUtils.getBufferedWriter(outFilenameDualEventsOneGene);
		String outfilenameDualEventsRecurrent = outDir + File.separator + "dualEventsRecurrent" + StringUtils.FileExtensionTSV.mExtension;
		BufferedWriter outDualEventsRecurrent = IOUtils.getBufferedWriter(outfilenameDualEventsRecurrent);
				
		StringBuilder sb = new StringBuilder(4096);		
		//TwoGenesAndEvents2 coOcurrenceWithPatientsDummy = new TwoGenesAndEvents2();		
		HashMap<PrimitiveWrapper.WLong, ArrayList<String>> coOcurrenceWithPatientsCompact = 
				new HashMap<PrimitiveWrapper.WLong, ArrayList<String>>();
		PrimitiveWrapper.WLong twoGenesTwoEventsCompactDummy = new PrimitiveWrapper.WLong(0);
		HashMap<PrimitiveWrapper.WLong, BitSet> coOcurrenceWithPatients = new HashMap<PrimitiveWrapper.WLong, BitSet>();
		
		NullaryClassFactory<DualEventsWithPatientsList> listClassFactory = new NullaryClassFactory<>(DualEventsWithPatientsList.class);
		EnumMatrix<EventType, DualEventsWithPatientsList> countMatrix = new EnumMatrix<>(EventType.class, listClassFactory);
		DualEventsWithPatients dualEventDummy = new DualEventsWithPatients();
				
		// Per-patient Map for associating gene with event(s)		
		ArrayListMultimap<GeneCounter, EventType> eventsForGene = ArrayListMultimap.create(20000, 2); 				
		
		// Define an object to keep the tallying results over the patients
		GenesTallier tallier = new GenesTallier();
		
		// Create the co-occurrence cache file
		System.out.println(outDir);
		String cacheFilename = outDir + File.separator + "GeneCacheFile." + NumberUtils.getRandomInteger(1, Integer.MAX_VALUE) + ".lohcate.txt";
		System.out.println(cacheFilename);
		File cacheFile = IOUtils.createBlankFile(cacheFilename);		
		LongArrayList allDualEvents = new LongArrayList();
		FastArrayListLong allDualEventsFast = new FastArrayListLong();
		
		ArrayList<String> allPatients = new ArrayList<String>(files.length);
		ArrayList<File> mValidFiles = new ArrayList<File>(files.length);		
		String targetSuffix = Clustering.SitesClassifiedOutputSuffix + ".lohcateInput_txt" ; // fileExtDelim.mExtension;
		System.out.println("Reading Patients...");
		for (File file : files) {  // Iterate through clustering/smoothing results			
			if (file.getName().endsWith(targetSuffix)) {
				String patientName = file.getName().substring(0, file.getName().indexOf(targetSuffix));
				allPatients.add(patientName);
				mValidFiles.add(file);					
				System.out.println(allPatients.size() + "\t" + file.getName() + "\t" + patientName);
			}
		}
		
		for (int patientIndex = 0; patientIndex < allPatients.size(); patientIndex++) {
			File patientFile = mValidFiles.get(patientIndex);
			String patientName = allPatients.get(patientIndex);
			
			System.out.printf("Patient: %d\t%s\t%s\n", patientIndex, patientFile.getName(), patientName);
			
			// Load all lines into memory, extract the header row, and then sort by chrom/position
			ArrayList<String> allLines = IOUtils.readAllLinesFromFile(patientFile.getAbsolutePath(), false, true, sb);
			String headerString = sb.toString();
			Collections.sort(allLines, Regions.LineComparatorTab);					
		
			tallier.performBasicGeneTallyingForPatient(patientName, allLines, false, eventsForGene);
			System.out.println("Num Events for Sample: " + eventsForGene.size());
						
			// Now that we've gotten all the gene-event combinations for a particular
			// sample, we store all the combinations.  We do an O(n^2) operation
			int numCombos = 0;			
			Collection<Map.Entry<GeneCounter, EventType>> keysAndValuesMap = eventsForGene.entries();
				
			for (GeneCounter gene : eventsForGene.keySet()) {
				List<EventType> events = eventsForGene.get(gene);
				int numEvents = events.size();
				if (numEvents > 1) {
					sb.setLength(0);
					sb.append(patientName);
					appendGeneCoordinates(gene, sb);
					sb.append(StringUtils.FileExtensionTSV.mDelimiter).append(numEvents);
					for (EventType event : events) {
						sb.append(StringUtils.FileExtensionTSV.mDelimiter)
						  .append(event.name());
					}
					IOUtils.writeToBufferedWriter(outDualEventsOneGene, sb.toString(), true);
				}
			}
			int numPairs = (keysAndValuesMap.size() * (keysAndValuesMap.size() - 1)) / 2;
			LongArrayList allPairsOnePatient = new LongArrayList(numPairs);
			
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

						twoGenesTwoEventsCompactDummy.mLong = compactIntoUnitLong(gene1.getNumericID(), gene2.getNumericID(), keyValueIter1.getValue(), keyValueIter2.getValue(), Cast.toShort(InitialDualEventCount));
						numCombos++;
						allPairsOnePatient.add(twoGenesTwoEventsCompactDummy.mLong);
						/*
							BitSet patients = coOcurrenceWithPatients.get(twoGenesTwoEventsCompactDummy);
							if (patients == null) {
								patients = new BitSet(256); 								
								coOcurrenceWithPatients.put(twoGenesTwoEventsCompactDummy, patients);
								twoGenesTwoEventsCompactDummy = new PrimitiveWrapper.WLong(0);
								numNewCombos++;
							}
							patients.set(patientIndex);
						 */

						/*
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
							int resultIndex = dualEventList.mList.getIndex(dualEventDummy);
							if (resultIndex < 0) {
								dualEventDummy.mPatients.clear();
								dualEventDummy.mPatients.set(patientIndex);
								dualEventList.mList.add(ArrayUtils.getInsertPoint(resultIndex), dualEventDummy);
								dualEventDummy = new DualEventsWithPatients();								
								numNewCombos++;
							} else {
								DualEventsWithPatients dualEvent = dualEventList.mList.get(resultIndex);
								CompareUtils.ensureTrue(dualEvent.compareTo(dualEventDummy) == 0, "ERROR: Impossible state with Dual Events array!");
								dualEvent.mPatients.set(patientIndex);
							}	
						 */

					}					
				}
			}

			System.out.println("Num Combos: " + numCombos);
			System.out.println("Sorting...");
			ArrayUtils.sort(allPairsOnePatient);
			System.out.println("Sorting done; now reading file...");
			//readFromCoOcurrenceFileAndUpdate(allPairsOnePatient, patientIndex, cacheFile.getAbsolutePath(), System.out);
			//readFromCoOccurrenceListAndUpdate(allPairsOnePatient, allDualEvents, System.out);
			//readFromCoOccurrenceListAndUpdate(allPairsOnePatient, allDualEventsFast, System.out);

			//System.out.println("Num New Combos: " + numNewCombos);
			eventsForGene.clear();
			System.gc();
			System.gc();
		}
		
		// Now count all the dual events
		int[] numPatientsTalliesForDualEvent = new int[allPatients.size() + 1];
		Arrays.fill(numPatientsTalliesForDualEvent, 0);
		for (FastArrayListIterator iter = allDualEventsFast.iterator(); iter.hasNextAndGet(twoGenesTwoEventsCompactDummy); ) {
			long compactUnit = twoGenesTwoEventsCompactDummy.mLong;
			int patientCount = Cast.toInt(TwoGenesWithEvents.Compactor.getValue(TwoGenesWithEvents.PatientCount, compactUnit));
			numPatientsTalliesForDualEvent[patientCount]++;			
		}
		
		int numEventsTotal = ArrayUtils.arraySum(numPatientsTalliesForDualEvent);
		System.out.printf("\n-----\nNumEventsTotal:\t%d\n", numEventsTotal);
		for (int i = 0; i < numPatientsTalliesForDualEvent.length; i++) {
			System.out.printf("Num Dual Events with\t%d\tpatients observing dual event:\t%d\n", i, numPatientsTalliesForDualEvent[i]);
		}
		
		int numEventsTop = Cast.toInt(numEventsTotal * 0.01);
		int numPatientsLowerBound = 1;
		int sum = 0;
		for (int i = numPatientsTalliesForDualEvent.length - 1; i >= 0; i--) {
			sum += numPatientsTalliesForDualEvent[i];
			if (sum >= numEventsTop) {
				numPatientsLowerBound = i;
				break;
			}
		}
		
		for (FastArrayListIterator iter = allDualEventsFast.iterator(); iter.hasNextAndGet(twoGenesTwoEventsCompactDummy); ) {
			long compactUnit = twoGenesTwoEventsCompactDummy.mLong;
			int patientCount = Cast.toInt(TwoGenesWithEvents.Compactor.getValue(TwoGenesWithEvents.PatientCount, compactUnit));
			if (patientCount >= numPatientsLowerBound && patientCount > 1) {
				int gene1ID = Cast.toInt(TwoGenesWithEvents.Compactor.getValue(TwoGenesWithEvents.ID_Gene1, compactUnit));
				int gene2ID = Cast.toInt(TwoGenesWithEvents.Compactor.getValue(TwoGenesWithEvents.ID_Gene2, compactUnit));
				GeneCounter gene1 = tallier.getGene(gene1ID);
				GeneCounter gene2 = tallier.getGene(gene2ID);
				EventType event1 = EventType.getClusterType( Cast.toInt(TwoGenesWithEvents.Compactor.getValue(TwoGenesWithEvents.Event_Gene1, compactUnit)) );
				EventType event2 = EventType.getClusterType( Cast.toInt(TwoGenesWithEvents.Compactor.getValue(TwoGenesWithEvents.Event_Gene2, compactUnit)) );
				double distanceScore = calcDistanceScore(gene1, gene2);
				double score = distanceScore * patientCount;
				if (score >= 0) {
				sb.setLength(0);
				sb.append(gene1.mChrom.getCode())
				  .append(StringUtils.FileExtensionTSV.mDelimiter).append(gene1.getName())
				  .append(StringUtils.FileExtensionTSV.mDelimiter).append(event1.name())
				  .append(StringUtils.FileExtensionTSV.mDelimiter).append(gene2.mChrom.getCode())
				  .append(StringUtils.FileExtensionTSV.mDelimiter).append(gene2.getName())
				  .append(StringUtils.FileExtensionTSV.mDelimiter).append(event2.name())
				  .append(StringUtils.FileExtensionTSV.mDelimiter).append(patientCount)
				  .append(StringUtils.FileExtensionTSV.mDelimiter).append((gene1.mChrom == gene2.mChrom ? 0 : 1))
				  .append(StringUtils.FileExtensionTSV.mDelimiter).append(score)
				  ;
				
				}
				IOUtils.writeToBufferedWriter(outDualEventsRecurrent, sb.toString(), true);
			}
		}		
		
		
		// Close the dual events file
		IOUtils.closeBufferedWriter(outDualEventsRecurrent);
		IOUtils.closeBufferedWriter(outDualEventsOneGene);
		
		// Print out the basic summary statistics of the genes
		printBasicGeneSummary(outFilename, fileExtDelim, tallier.getGenes());
		
		// Print out the breakdown of (Gene, Event, Patient)
		String eventsByPatientPerGeneFilename = outDir + File.separator + "eventsByPatientPerGene" + StringUtils.FileExtensionTSV.mExtension;
		printGenePatientBreakdown(eventsByPatientPerGeneFilename, tallier.getGenes());
		
		
		// Store co-occurrence counts			
		String coOccurrenceFilename = outDir + File.separator + "coOccurrence" + StringUtils.FileExtensionTSV.mExtension;
		BufferedWriter outCoOccurrence = IOUtils.getBufferedWriter(coOccurrenceFilename);		
		for (Map.Entry<PrimitiveWrapper.WLong, BitSet> entry : coOcurrenceWithPatients.entrySet()) {
			
			PrimitiveWrapper.WLong compactUnit = entry.getKey();			
			int geneID1 = Cast.toInt(TwoGenesWithEvents.Compactor.getValue(TwoGenesWithEvents.ID_Gene1, compactUnit.mLong));
			int geneID2 = Cast.toInt(TwoGenesWithEvents.Compactor.getValue(TwoGenesWithEvents.ID_Gene2, compactUnit.mLong));
			GeneCounter gene1 = tallier.getGene(geneID1);
			GeneCounter gene2 = tallier.getGene(geneID2);
			
			
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
		
		cacheFile.delete();
		IOUtils.closeBufferedWriter(outCoOccurrence);
	}

	// ========================================================================
	private static double calcDistanceScore(GeneCounter gene1, GeneCounter gene2) {
		if (gene1.getChrom() == gene2.getChrom()) {
			if (gene1.mMinBasePairPosition > gene2.mMinBasePairPosition) {
				return calcDistanceScore(gene2, gene1);
			} else {
				int diffDistance = gene2.mMinBasePairPosition - gene1.mMaxBasePairPosition;
				return ((double) diffDistance) / (double) gene1.getChrom().getLength();
			}
		} else {
			return 1;
		}
	}

	// ========================================================================
	private static StringBuilder appendGeneCoordinates(GeneCounter gene, StringBuilder sb) {
		sb.append(StringUtils.FileExtensionTSV.mDelimiter).append(gene.getChrom().ordinal())
		  .append(StringUtils.FileExtensionTSV.mDelimiter).append(gene.getChrom().getArm(gene.mMaxBasePairPosition).name())							  
		  .append(StringUtils.FileExtensionTSV.mDelimiter).append(gene.mMinBasePairPosition)
		  .append(StringUtils.FileExtensionTSV.mDelimiter).append(gene.mMaxBasePairPosition)
		  .append(StringUtils.FileExtensionTSV.mDelimiter).append(gene.getName())
		  ;
		return sb;
	}
	
	// ========================================================================
	private static void printGenePatientBreakdown(String outFilename, ArrayList<GeneCounter> genes) {
		BufferedWriter outBreakdown = IOUtils.getBufferedWriter(outFilename);
		
		StringBuilder sb = new StringBuilder(1024);
		// Now write out individual gene outputs
		for (GeneCounter gene : genes) {
			// Write samples for each gene out as well			
			for (EventType eventType : EventType.values()) {
				if (isEventToIgnore(eventType)) continue;
				
				ArrayList<String> patientsForEvent = gene.getPatientsForEventType(eventType);
				for (String patientForEvent : patientsForEvent) {					
					sb.setLength(0);
					appendGeneCoordinates(gene, sb);
					sb.append(StringUtils.FileExtensionTSV.mDelimiter).append(eventType.name());
					sb.append(StringUtils.FileExtensionTSV.mDelimiter).append(patientForEvent);
					IOUtils.writeToBufferedWriter(outBreakdown, sb.toString(), true);
				}
			}
			
			for (EventTypeAllele eta : EventTypeAllele.values()) {
				ArrayList<String> patientsForEvent = gene.getPatientsForAlleleEventType(eta);
				for (String patientForEvent : patientsForEvent) {
					sb.setLength(0);
					appendGeneCoordinates(gene, sb);
					sb.append(StringUtils.FileExtensionTSV.mDelimiter).append(eta.name());					
					sb.append(StringUtils.FileExtensionTSV.mDelimiter).append(patientForEvent);
					IOUtils.writeToBufferedWriter(outBreakdown, sb.toString(), true);
				}				
			}
		}
		IOUtils.closeBufferedWriter(outBreakdown);
	}
	
	// ========================================================================
	private static void printBasicGeneSummary(String outFilename, StringUtils.FileExtensionAndDelimiter fileExtDelim, ArrayList<GeneCounter> genes) {
		String[] columnHeaders = new String[] { 
				"chr", "arm", "bp_start", "bp_end", "length", "gene", GeneCounter.getHeaders(fileExtDelim.mDelimiter)
		};		
		StringBuilder sb = new StringBuilder(2048);
		String headerStr = StringUtils.constructColumnDelimitedString(columnHeaders, fileExtDelim.mDelimiter, sb, true).toString();
	
		BufferedWriter out = IOUtils.getBufferedWriter(outFilename);
		IOUtils.writeToBufferedWriter(out, headerStr, true);
		for (GeneCounter gene : genes) {			
			sb.setLength(0);
			sb.append(gene.mChrom.getCode())
			    .append(fileExtDelim.mDelimiter).append(gene.getChrom().getArm(gene.mMaxBasePairPosition))
				.append(fileExtDelim.mDelimiter).append(gene.mMinBasePairPosition)
				.append(fileExtDelim.mDelimiter).append(gene.mMaxBasePairPosition)
				.append(fileExtDelim.mDelimiter).append(gene.getRangeLength())
				.append(fileExtDelim.mDelimiter).append(gene.mLabel)
				.append(fileExtDelim.mDelimiter).append(gene.countsToString(fileExtDelim.mDelimiter));
			IOUtils.writeToBufferedWriter(out, sb.toString(), true);
			IOUtils.flushBufferedWriter(out);
		}
		IOUtils.closeBufferedWriter(out);
	}

	// ========================================================================
	private static void coOccurrenceAnalysis() {
		
	}

	// ========================================================================
	private static void readFromCoOccurrenceListAndUpdate(LongArrayList coOccurrencesOnePatient, FastArrayListLong allDualEvents, PrintStream outStream) {
		
		allDualEvents.setValueExtractor(TwoGenesWithEvents.Extractor);
		
		int numCoOccurrencesAlreadyExisted = 0;
		int[] numXTons = new int[92]; 		
		
		// Iterate through the dual events in the patient
		FastArrayListIndex indexAndValue = FastArrayListLong.getNewIndexToken();
		for (int indexInPatient = 0; indexInPatient < coOccurrencesOnePatient.size(); indexInPatient++) {
			long dualEventInPatient = coOccurrencesOnePatient.get(indexInPatient);
			FastArrayListIndex exists = allDualEvents.getIndex(dualEventInPatient, indexAndValue);
			long patientCount = 1;
			
			if (exists == null) {
				// Does not exist.  Just add
				allDualEvents.add(dualEventInPatient);				
			} else {				
				patientCount = TwoGenesWithEvents.Compactor.getValue(TwoGenesWithEvents.PatientCount, indexAndValue.mValue);
				long replaceUnit  = TwoGenesWithEvents.Compactor.setValue(TwoGenesWithEvents.PatientCount, ++patientCount, indexAndValue.mValue);
				allDualEvents.replace(indexAndValue, replaceUnit);
				++numCoOccurrencesAlreadyExisted;
			}	
			
			++numXTons[Cast.toInt(Math.min(patientCount, numXTons.length - 1))];
		}
		
		outStream.printf("Num Co-occurrences Visited, New, Total:\t%d\t%d\t%d\n", 
				numCoOccurrencesAlreadyExisted, coOccurrencesOnePatient.size() - numCoOccurrencesAlreadyExisted, coOccurrencesOnePatient.size());

		for (int j = 1; j < numXTons.length; j++) {
			outStream.printf("Num\t%d\t-ton:\t%d\n", j, numXTons[j]);					
		}

	}
	
	// ========================================================================
	private static void readFromCoOccurrenceListAndUpdate(LongArrayList coOccurrencesOnePatient, LongArrayList allDualEvents, PrintStream outStream) {
		
		// Set a bitset to see which entries in the patient list were in the master list
		BitSet coOccurenceNotVisited = new BitSet(coOccurrencesOnePatient.size());
		coOccurenceNotVisited.set(0, coOccurrencesOnePatient.size());  // set all as not visited
		
		int[] numXTons = new int[92]; 		
		int numCoOccurrencesAlreadyExisted = 0;
		
		int arraySize = allDualEvents.size();
		for (int i = 0; i < arraySize; i++) {
			long masterUnit = allDualEvents.get(i);
			long patientCount = TwoGenesWithEvents.Compactor.getValue(TwoGenesWithEvents.PatientCount, masterUnit);
			
			// Reset the count in the master so that it can be searched within the patient list
			long searchUnit = TwoGenesWithEvents.Compactor.setValue(TwoGenesWithEvents.PatientCount, InitialDualEventCount, masterUnit);
			
			// Now search for the unit in the co-occurrences
			int resultIndex = ArrayUtils.binarySearchValue(searchUnit, coOccurrencesOnePatient, ValueExtractor.LongExtractorWhole);
			if (resultIndex >= 0) {				
				long masterUnitNew = TwoGenesWithEvents.Compactor.setValue(TwoGenesWithEvents.PatientCount, ++patientCount, masterUnit);
				allDualEvents.set(i, masterUnitNew);
				coOccurenceNotVisited.set(resultIndex, false);
				++numCoOccurrencesAlreadyExisted;
			}
			
			++numXTons[Cast.toInt(Math.min(patientCount, numXTons.length - 1))];
		}
		
		outStream.printf("Num Co-occurrences Visited, New, Total:\t%d\t%d\t%d\n", 
				numCoOccurrencesAlreadyExisted, coOccurrencesOnePatient.size() - numCoOccurrencesAlreadyExisted, coOccurrencesOnePatient.size());
		
		for (int j = 1; j < numXTons.length; j++) {
			outStream.printf("Num %d-ton:\t%d\n", j, numXTons[j]);					
		}
		
		
		allDualEvents.ensureCapacity(allDualEvents.size() + coOccurenceNotVisited.cardinality());
		for (int bitIndex = coOccurenceNotVisited.nextSetBit(0); bitIndex >= 0; bitIndex = coOccurenceNotVisited.nextSetBit(bitIndex + 1)) {
			long compactUnit = coOccurrencesOnePatient.get(bitIndex);
			allDualEvents.add(compactUnit);
		}
		
		outStream.println("Sorting All Events...");
		//ArrayUtils.sort(allDualEvents);
		outStream.println("Sorting All Events Finished...");
	}
	
	// ========================================================================
	private static void readFromCoOcurrenceFileAndUpdate(LongArrayList coOccurrencesOnePatient, int patientIndex, String filename, PrintStream outStream) {
		File inFile = new File(filename);		
		BufferedReader in = IOUtils.getBufferedReader(filename);
		
		File outFile = new File(filename + ".temp.txt");
		BufferedWriter out = IOUtils.getBufferedWriter(outFile.getAbsolutePath());
		ArrayUtils.sort(coOccurrencesOnePatient);
		
		String line;		
		StringBuilder sb = new StringBuilder(2048);
		BitSet coOccurenceNotVisited = new BitSet(coOccurrencesOnePatient.size());
		coOccurenceNotVisited.set(0, coOccurrencesOnePatient.size());  // set all as not visited
		
		int numSingletons = 0;
		int numCoOccurrencesAlreadyExisted = 0;
		boolean timeToPrune = ((patientIndex % 8) == 0);
		
		while ((line = IOUtils.getNextLineInBufferedReader(in)) != null) {
			String compactValueStr = StringUtils.extractNthColumnValue(line, 0, StringUtils.TabStr);
			String countStr        = StringUtils.extractNthColumnValue(line, 1, StringUtils.TabStr);
			String bitString =       StringUtils.extractNthColumnValue(line, 2, StringUtils.TabStr);
			
			long compactValue = Long.parseLong(compactValueStr, 16);
			int eventCount = Integer.parseInt(countStr);
			String newBitString = bitString; // default
						
			int resultIndex = ArrayUtils.binarySearchValue(compactValue, coOccurrencesOnePatient, ValueExtractor.LongExtractorWhole);
			if (resultIndex >= 0) {
				newBitString = BitSetUtils.setBit(bitString, patientIndex, true);
				coOccurenceNotVisited.set(resultIndex, false);
				++numCoOccurrencesAlreadyExisted;
				++eventCount;
			} 
			int x = Character.MAX_RADIX;
			if (eventCount == 1) {
				++numSingletons;
			}
			
			boolean shouldPrint = !(timeToPrune && ((eventCount == 1) || (patientIndex >= 40 && eventCount <= 2)));
			
			// Now, write to the output
			if (shouldPrint) {
				sb.setLength(0);
				sb.append(compactValueStr).append(StringUtils.TabStr).append(eventCount).append(StringUtils.TabStr).append(newBitString);
				IOUtils.writeToBufferedWriter(out, sb.toString(), true);
			}
		}
		
		// Close the input file
		IOUtils.closeBufferedReader(in);
		
		outStream.printf("Num Co-occurrences Visited, New, Total:\t%d\t%d\t%d\n", 
				numCoOccurrencesAlreadyExisted, coOccurrencesOnePatient.size() - numCoOccurrencesAlreadyExisted, coOccurrencesOnePatient.size());
		outStream.printf("Num Singletons: %d\n", numSingletons);
		
		// Now, there are entries in the input list that were not in the iput file.  We need to write these entries.
		for (int bitIndex = coOccurenceNotVisited.nextSetBit(0); bitIndex >= 0; bitIndex = coOccurenceNotVisited.nextSetBit(bitIndex + 1)) {
			long compactUnit = coOccurrencesOnePatient.get(bitIndex);
			String compactUnitStr = Long.toHexString(compactUnit);
			String newBitString = BitSetUtils.setBit("", patientIndex, true);
			
			sb.setLength(0);
			sb.append(compactUnitStr).append(StringUtils.TabStr).append(1).append(StringUtils.TabStr).append(newBitString);			
			IOUtils.writeToBufferedWriter(out, sb.toString(), true);			
		}
		
		IOUtils.closeBufferedWriter(out);
		
		// Now rename
		inFile.delete();
		outFile.renameTo(inFile);
	}

	// ========================================================================
	private static int compactIntoUnit(int numIDGene1, int numIDGene2) {
		if (numIDGene1 > numIDGene2) {
			return compactIntoUnit(numIDGene2, numIDGene1);
		} else {
			int unit = 0;
			unit = TwoGenes.Compactor.setValue(TwoGenes.ID_Gene1, numIDGene1, unit);
			unit = TwoGenes.Compactor.setValue(TwoGenes.ID_Gene2, numIDGene2, unit);
			return unit;
		}
	}
	
	// ========================================================================
	private static long compactIntoUnitLong(int numIDGene1, int numIDGene2, EventType eventGene1, EventType eventGene2, short count) {		
		if (numIDGene1 > numIDGene2) {
			return compactIntoUnitLong(numIDGene2, numIDGene1, eventGene2, eventGene1, count);  // Call recursively
		} else {
			long unit = 0;
			unit = TwoGenesWithEvents.Compactor.setValue(TwoGenesWithEvents.ID_Gene1, numIDGene1, unit);
			unit = TwoGenesWithEvents.Compactor.setValue(TwoGenesWithEvents.ID_Gene2, numIDGene2, unit);
			unit = TwoGenesWithEvents.Compactor.setValue(TwoGenesWithEvents.Event_Gene1, eventGene1.ordinal(), unit);
			unit = TwoGenesWithEvents.Compactor.setValue(TwoGenesWithEvents.Event_Gene2, eventGene2.ordinal(), unit);
			unit = TwoGenesWithEvents.Compactor.setValue(TwoGenesWithEvents.PatientCount, count, unit);
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
	public static class TwoGenesListWithPatients {
		
		IntArrayList mGenePairs;
		BitSet mPatients;
		int mNumPatients;
		
		NullaryClassFactory<TwoGenesListWithPatients> ClassFactory = new NullaryClassFactory<>(TwoGenesListWithPatients.class);
		
		public TwoGenesListWithPatients(int numPatients) {
			mNumPatients = numPatients;
			mGenePairs = new IntArrayList(10000000);
			mPatients = new BitSet(100000000);
		}
		
		public void setPatient(int pairIndex, int patientIndex, boolean value) {
			int trueIndex = (pairIndex * mNumPatients) + patientIndex;
			mPatients.set(trueIndex, value);
		}
		
	}
	
	// ========================================================================
	public static class DualEventsWithPatients implements Comparable<DualEventsWithPatients> {
		protected long mDualEventCompact;
		protected BitSet mPatients;
		
		public DualEventsWithPatients() {
			mDualEventCompact = 0;
			mPatients = new BitSet(256);
		}
		
		public void setEvents(int numIDGene1, int numIDGene2, EventType eventGene1, EventType eventGene2) {
			mDualEventCompact = compactIntoUnitLong(numIDGene1, numIDGene2, eventGene1, eventGene2, Cast.toShort(0));
		}

		@Override
		public int compareTo(DualEventsWithPatients o) {			
			return Long.compare(this.mDualEventCompact, o.mDualEventCompact);
		}
		
		
	}

	// ========================================================================
	public static enum TwoGenes implements CompactorInf<TwoGenes> {
		ID_Gene1(16),
		ID_Gene2(16)
		;
		
		public static CompactorIntoInt<TwoGenes> Compactor = new CompactorIntoInt<TwoGenes>(TwoGenes.class, true);
		
		private short mNumBits;
		private TwoGenes(int numBits) {
			mNumBits = Cast.toShort(numBits);
		}
		
		@Override
		public int getNumBits() { return mNumBits; }
	}
	
	// ========================================================================
	public static enum TwoGenesWithEvents implements CompactorInf<TwoGenesWithEvents> {
		ID_Gene1(16),
		ID_Gene2(16),
		Event_Gene1(4),
		Event_Gene2(4),
		PatientCount(16);
		;

		public static CompactorIntoLong<TwoGenesWithEvents> Compactor = new CompactorIntoLong<TwoGenesWithEvents>(TwoGenesWithEvents.class, false);
		
		private short mNumBits;
		
		private TwoGenesWithEvents(int numBits) {
			mNumBits = Cast.toShort(numBits);
		}

		@Override
		public int getNumBits() { return mNumBits; }
		
		public int compare_ExcludingCounts(long dualEvent1, long dualEvent2) {
			return Long.compare(dualEvent1 >>> PatientCount.mNumBits, dualEvent2 >>> PatientCount.mNumBits);
		}
		
		public static ValueExtractor Extractor = new ValueExtractor() {			
			@Override
			public long extractValue(long compactUnit) { return (compactUnit >>> PatientCount.mNumBits); }
		};
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
