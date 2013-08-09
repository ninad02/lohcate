package JISTIC;

import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;

public class CallAllByChrom {

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
		for (int i = 1; i <= 22; i++) 
		{					
			String chromFolder=folder+"/chrom"+i;			
			File newFolder=new File(chromFolder);
			newFolder.mkdir();						
			String[] args2=new String[] {chromFolder};
			String[] newargs=new String[args.length+1];
			for (int j = 1; j < args.length; j++) {
				newargs[j]=args[j];
			}
			newargs[0]="spec="+newFolder+"/"+specFile.getName();
			copyFile(specfileLoc,newFolder+"/"+specFile.getName());
			newargs[args.length]="singlechrom="+i;
			System.out.println("Running on chromosome "+i);
			Distribution.main(newargs);
			Convert2IGV.main(args2);
			CombineMatricesNoParam.main(args2);
		}				
	}
	
	public static void copyFile(String infile,String outfile) throws Exception
	{
		FileReader in = new FileReader(infile);
	    FileWriter out = new FileWriter(outfile);
	    int c;
	    while ((c = in.read()) != -1)
	      out.write(c);
	    in.close();
	    out.close();
	}

}
