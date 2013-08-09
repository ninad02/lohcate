package JISTIC;

import java.io.FileOutputStream;
import java.io.PrintStream;
import java.util.Iterator;

public class CreateIGVSessionFile {

	/**
	 * @param args
	 */
	public static void main(String[] args) throws Exception{
		String[] folders=args;
		String outputFile=null;
		String segFile=null;
		String[] peakFiles=null;
		String[] gscoreFiles=null;
		try
		{
			outputFile=args[0];
			segFile=args[1];
			peakFiles=args[2].split(",");
			gscoreFiles=args[3].split(",");
		}catch(Exception e)
		{
			System.err.println("Format: outputfile segfile commaseparatedpeakfiles commaseparatedgscorefiles");
			throw new RuntimeException("Problems parsing the input for CreateIGVSessionFile");
		}
		
		PrintStream fout = new PrintStream(new FileOutputStream(outputFile));
		fout.println("<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\"?>");
		fout.println("<Global genome=\"hg18\" version=\"2\">");
		fout.println("\t<Files>");
		// input files here
		fout.println("\t\t<DataFile name=\""+segFile+"\" relativePath=\"false\"/>");
		for (int i = 0; i < peakFiles.length; i++) {
			fout.println("\t\t<DataFile name=\""+peakFiles[i]+"\" relativePath=\"false\"/>");
		}
		for (int i = 0; i < peakFiles.length; i++) {
			fout.println("\t\t<DataFile name=\""+gscoreFiles[i]+"\" relativePath=\"false\"/>");
		}		
		//fout.println("\t\t");
		//fout.println("\t\t");
		fout.println("\t</Files>");
		fout.println("\t<Panel name=\"DataPanel\">");
		// all samples here, it works if you leave it empty
		//fout.println("\t\t");
		fout.println("\t</Panel>");
		fout.println("\t<Panel name=\"FeaturePanel\">");
		// gene track first
		fout.println("\t\t<Track color=\"0,0,178\" displayName=\"Gene\" expand=\"false\" height=\"59\" id=\"Gene\" renderer=\"GENE_TRACK\" scale=\"0.0,0.0,10.0\" visible=\"true\" windowFunction=\"count\"/>");
		// jistic tracks here, it works if you leave it empty 
		fout.println("\t</Panel>");		
		/*
		// filter tracks you are interested (seg samples and gistic results)
		fout.println("\t<Filter match=\"all\" name=\"\" showTracks=\"false\">");
		fout.println("\t\t<FilterElement booleanOperator=\"OR\" item=\"NAME\" operator=\"contains\" value=\"PEAK\"/>");		
		fout.println("\t</Filter>");
		*/		
		fout.println("</Global>");
		fout.close();
	}

}
