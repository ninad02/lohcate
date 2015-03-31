package lohcate.clustering;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

import java.io.File;
import java.util.Map;

import nutils.PrimitiveWrapper;
import nutils.counter.ObjectCounter;

public class Reader {

	public static SamReaderFactory SamReaderFactoryInstance = SamReaderFactory.makeDefault();
	
	// ========================================================================
	public static ClusteringInputOneSample readBAMPair(String inFile1, String inFile2) {
		return readBAMPair(new File(inFile1), new File(inFile2));
	}
	
	// ========================================================================
	public static ClusteringInputOneSample readBAMPair(File inFile1, File inFile2) {
		SamReader reader1 = SamReaderFactoryInstance.open(inFile1);
		SamReader reader2 = SamReaderFactoryInstance.open(inFile2);
		
		SAMRecordIterator iter1 = reader1.iterator();
		SAMRecordIterator iter2 = reader2.iterator();
		
		ObjectCounter<String> oc = new ObjectCounter<>();
		
		while (iter1.hasNext()) {
			SAMRecord record1 = iter1.next();
			String refName1 = record1.getReferenceName();
			
//			record1.get
//			
//			oc.increment(refName1);
			
			//System.out.println(refName1);
			//System.out.println("\t" + record1.getReferenceIndex());
		}
		
		for (Map.Entry<String, PrimitiveWrapper.WInteger> entry : oc.getAllCounts()) {
			System.out.println(entry.getKey() + "\t" + entry.getValue().mInt);
		}
	
		return null;
	}

	// ========================================================================
	public Reader() {
		// TODO Auto-generated constructor stub
	}

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		Reader.readBAMPair(args[0], args[0]);
	}

}
