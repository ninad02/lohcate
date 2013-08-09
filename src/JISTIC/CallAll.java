package JISTIC;

import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;

public class CallAll {

	/**
	 * @param args
	 */
	public static void main(String[] args) throws Exception {
		
		String specfileLoc=args[0].split("=")[1];
		File specFile=new File(specfileLoc);
		System.out.println();
		String folder=specFile.getParent();
		if(folder==null)
			folder=".";
		String[] args2=new String[] {folder};
		Distribution.main(args);
		Convert2IGV.main(args2);
		CombineMatricesNoParam.main(args2);
		
	}
	
	

}
