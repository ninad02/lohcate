package shared;

import java.io.File;
import java.io.IOException;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.InvalidPathException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

public class TCGASampleListAnalysis {

	// ========================================================================
	public static final String PrefixTCGA = "TCGA";
	public static final String DelimiterTCGABarcode = "-";
	
	public static final int NormalThreshold = 10;
	public static final int Col_TCGAPrefix  = 0;
	public static final int Col_CenterName  = 1;
	public static final int Col_PatientName = 2;
	public static final int Col_SampleRole  = 3;
	public static final int NumDigitsInRole = 2;
	public static final int NumColsInBarcode = 7;
	
	// ========================================================================
	
	
	// ========================================================================
	public TCGASampleListAnalysis() {
		// TODO Auto-generated constructor stub
	}
	
	public static void performAnalysis(String inFilePath) {
		ArrayList<String> allLines = readAllLinesFromFile(inFilePath);
		HashMap<String, TCGAPatient> patientMap = new HashMap<String, TCGAPatient>();
		TCGASample dummySample = new TCGASample();
		ArrayList<TCGASample> successfulSamples = new ArrayList<TCGASample>();
		
		for (String line : allLines) {
			String[] cols = line.split(DelimiterTCGABarcode);			
			
			if (!cols[Col_TCGAPrefix].equalsIgnoreCase(PrefixTCGA)) {
				System.out.println("Does not contain TCGA Prefix, Ignoring: " + line);
			} else if (cols.length < NumColsInBarcode) {
				System.out.printf("Number of columns (%d) does not meet the required number of columns (%d) in TCGA Barcode, Ignoring: (%s)\n", cols.length, NumColsInBarcode, line);						
			} else {
				  TCGAPatient patient = patientMap.get(cols[Col_PatientName]);
				  if (patient == null) {
					  patient = new TCGAPatient(cols[Col_CenterName], cols[Col_PatientName]);
					  patientMap.put( cols[Col_PatientName], patient );
				  }
				  
				  // Check the role
				  dummySample.mFullBarcode = line;
				  dummySample.mSampleRoleName = cols[Col_SampleRole]; 				  
				  int role = new Integer(dummySample.mSampleRoleName.substring(0, NumDigitsInRole));
				  dummySample.mTissue = (role >= NormalThreshold) ? Tissue.Normal : Tissue.Tumor;
				  
				  
				  ArrayList<TCGASample> samplesInPatient = (dummySample.mTissue == Tissue.Normal) ? patient.mNormals : patient.mTumors;
				  
				  
				  if (samplesInPatient.contains(dummySample)) {
					  System.out.println("Duplicate entry, ignoring: " + line);
				  } else {
					  samplesInPatient.add(dummySample);
					  successfulSamples.add(dummySample);
					  dummySample = new TCGASample();					  
				  }
			}
		}

		// Now finally do the output
		System.out.printf("Successfully added (%d) samples!\n", successfulSamples.size());
		
		//ArrayList<String> pairs 
		for (TCGAPatient patient : patientMap.values()) {
			
		}
	}
	
	// ========================================================================
	public static ArrayList<String> readAllLinesFromFile(String inFilePath) {
		Charset utf8 = Charset.forName("UTF-8");
		List<String> allLines = null;
		try {
			allLines = Files.readAllLines((new File(inFilePath)).toPath(), utf8);			
		} catch (IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}
		
		return new ArrayList<String>(allLines);
	}

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub

	}

	// ========================================================================
	public static class TCGAPatient {
		String mCenterName;
		String mPatientName;				
		ArrayList<TCGASample> mNormals;
		ArrayList<TCGASample> mTumors;
		
		public TCGAPatient (String centerName, String patientName) {
			mCenterName = centerName;
			mPatientName = patientName;			
			mNormals = new ArrayList<TCGASample>();
			mTumors  = new ArrayList<TCGASample>();			
		}

		/* (non-Javadoc)
		 * @see java.lang.Object#hashCode()
		 */
		@Override
		public int hashCode() {
			final int prime = 31;
			int result = 1;
			result = prime * result
					+ ((mCenterName == null) ? 0 : mCenterName.hashCode());
			result = prime * result
					+ ((mPatientName == null) ? 0 : mPatientName.hashCode());
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
			TCGAPatient other = (TCGAPatient) obj;
			if (mCenterName == null) {
				if (other.mCenterName != null)
					return false;
			} else if (!mCenterName.equalsIgnoreCase(other.mCenterName))
				return false;
			if (mPatientName == null) {
				if (other.mPatientName != null)
					return false;
			} else if (!mPatientName.equalsIgnoreCase(other.mPatientName))
				return false;
			return true;
		}
		
		
	}

	// ========================================================================
	public static enum Tissue { 
		Normal, Tumor;
	}
	
	// ========================================================================
	public static class TCGASample {
		String mSampleRoleName;  // 01A, 11A, 10
		String mFullBarcode;
		Tissue mTissue;
		
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
			TCGASample other = (TCGASample) obj;
			if (mFullBarcode == null) {
				if (other.mFullBarcode != null)
					return false;
			} else if (!mFullBarcode.equalsIgnoreCase(other.mFullBarcode))
				return false;
			return true;
		}
	}
}
