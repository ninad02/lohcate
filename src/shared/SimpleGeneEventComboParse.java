package shared;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;

import lohcate.geneAnalysis.GeneEnrichment.TwoGenesAndEvents2;
import lohcateEnums.EventType;

import nutils.IOUtils;

public class SimpleGeneEventComboParse {

	public SimpleGeneEventComboParse(String inFilename) {
		// TODO Auto-generated constructor stub
		ArrayList<String> allLines = IOUtils.readAllLinesFromFile(inFilename);

		HashMap<String, ArrayList<GeneEvent>> geneEventByPatient = new HashMap<String, ArrayList<GeneEvent>>();
		HashMap<TwoGenesAndEvents2, ArrayList<String>> doubleEventCounts = new HashMap<>();
		TwoGenesAndEvents2 doubleEventDummy = new TwoGenesAndEvents2();
		
		for (int i = 0; i < allLines.size(); i++) {
			String line = allLines.get(i);
			String[] cols = line.split("\\s");
						
			String gene = cols[0];
			EventType event = EventType.valueOf(cols[1]);
			String patient = cols[2];			
					
			ArrayList<GeneEvent> eventsForPatient = geneEventByPatient.get(patient); 
			if (eventsForPatient == null) {
				eventsForPatient = new ArrayList<GeneEvent>();
				geneEventByPatient.put(patient, eventsForPatient);				
			}
			eventsForPatient.add(new GeneEvent(gene, event));
		}
		
		// Now, we have the events by patient
		for (Map.Entry<String, ArrayList<GeneEvent>> entry : geneEventByPatient.entrySet()) {
			String patient = entry.getKey();
			ArrayList<GeneEvent> geneEvents = entry.getValue();
			
			for (GeneEvent ge1 : geneEvents) {
				for (GeneEvent ge2 : geneEvents) {
					if (ge1.equals(ge2)) break;
					
					doubleEventDummy.set(ge1.mGene, ge2.mGene, ge1.mEvent, ge2.mEvent);
					ArrayList<String> patientsForDoubleEvent = doubleEventCounts.get(doubleEventDummy);
					if (patientsForDoubleEvent == null) {
						patientsForDoubleEvent = new ArrayList<>();
						doubleEventCounts.put(doubleEventDummy, patientsForDoubleEvent);
						doubleEventDummy = new TwoGenesAndEvents2();
					}
					
					if (patientsForDoubleEvent.contains(patient)) {
						System.err.println("ERROR: Patient added twice for double-event!\t" + patient);
					} else {
						patientsForDoubleEvent.add(patient);
					}						
				}
			}			
		}
		
		for (Map.Entry<TwoGenesAndEvents2, ArrayList<String>> entry : doubleEventCounts.entrySet()) {
			TwoGenesAndEvents2 doubleEvent = entry.getKey();
			ArrayList<String> patients = entry.getValue();
			System.out.println(doubleEvent.toString() + "\t" + patients.size() + "\t" + patients.toString());
		}
	}

	// =====
	public static class GeneEvent {
		String mGene;
		EventType mEvent;
		
		public GeneEvent(String gene, EventType event) {
			mGene = gene;
			mEvent = event;
		}

		/* (non-Javadoc)
		 * @see java.lang.Object#hashCode()
		 */
		@Override
		public int hashCode() {
			final int prime = 31;
			int result = 1;
			result = prime * result
					+ ((mEvent == null) ? 0 : mEvent.hashCode());
			result = prime * result + ((mGene == null) ? 0 : mGene.hashCode());
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
			GeneEvent other = (GeneEvent) obj;
			if (mEvent != other.mEvent)
				return false;
			if (mGene == null) {
				if (other.mGene != null)
					return false;
			} else if (mGene.compareToIgnoreCase(other.mGene) != 0) 					
				return false;
			return true;
		}
		

	}
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		new SimpleGeneEventComboParse(args[0]);
	}

}
