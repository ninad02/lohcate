package lohcate.clustering;

import java.io.File;

import nutils.IOUtils;

public class JISTICWrapper {

	public JISTICWrapper() {
		// TODO Auto-generated constructor stub
	}

	public static void callJISTIC(String inputCopyNumberFile, boolean doFocal) {
		String jisticInputPathPrefix = IOUtils.pathConcat(new String[] { "E:", "Research", "Code", "Libraries", "Java", "lib", "JISTIC" }) + File.separator; 
		String jisticlocationsString = "locations=" + jisticInputPathPrefix + "hg18_Gene_Info.txt";
		String jisticCytobandString  = "bands=" + jisticInputPathPrefix + "Human_cytoBand.txt";
		String jisticCopyNumberString = "copynumber=" + inputCopyNumberFile;
		

		String jisticInputSpecFile = doFocal ?
				IOUtils.pathConcat(new String[] { "spec=" + jisticInputPathPrefix + "glioexample", "focal", "GISTICFocal.spec" }) :		
				IOUtils.pathConcat(new String[] { "spec=" + jisticInputPathPrefix + "glioexample", "limited", "GISTICLimited.spec" });
				
		JISTIC.Distribution.main(new String[] { jisticInputSpecFile, jisticCopyNumberString, jisticCytobandString, jisticlocationsString });
	}
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		callJISTIC(args[0], Boolean.parseBoolean(args[1]));
	}

}
