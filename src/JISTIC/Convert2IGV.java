package JISTIC;

import java.io.BufferedReader;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.PrintStream;
import java.util.HashSet;

public class Convert2IGV {

	
	public Convert2IGV()
	{
		
	}
	/**
	 * @param args
	 */
	public static void main(String[] args) throws Exception {
		String folder=".";
		if(args.length>0)
		{
			folder=args[0];
		}
		Convert2IGV c2i=new Convert2IGV();
		c2i.convertMarkerFile(folder, "markers.txt");
		c2i.convertPeaksFile(folder, "peaks.txt");
	}
	
	public void convertPeaksFile(String folder,String fileName) throws Exception
	{
		PrintStream fout = new PrintStream(new FileOutputStream(folder+"/"+"peaks.gp_gistic.gistic.txt"));
		fout.print("Type	Chromosome	Start	End	-log10(q-value)	G-score	average amplitude	frequency");
				
		BufferedReader fin=new BufferedReader(new FileReader(folder+"/"+fileName));		
		String line=fin.readLine();
		line=fin.readLine();		
		while(line!=null)
		{
			String[] tokens=line.split("\t");			
			String[] loc=tokens[0].split(":");
			int chromosome=Integer.parseInt(loc[2].substring(3));
			String[] locs=loc[3].split("\\(")[0].split("-");
			int startlocation=Integer.parseInt(locs[0]);
			int endlocation=startlocation;
			if(locs.length>1)
			{
				endlocation=Integer.parseInt(locs[1]);
			}
			String type="Amp";
			float score=1.0f;
			if(loc[1].equals("DEL"))
			{
				type="Del";
				score=1.1f;
			}
			
			// need to enter 0 before peak
			printline(fout,type,chromosome, startlocation-2, startlocation-1, 0.0,0);						
			// peak
			printline(fout,type,chromosome, startlocation, endlocation, score,(int)score);
			// need to enter 0 also after peak
			printline(fout,type,chromosome, endlocation+1, endlocation+2, 0.0,0);
			
			line=fin.readLine();	
		}
		
		fin.close();
		fout.close();
	}
	
	
	public void convertMarkerFile(String folder,String fileName) throws Exception
	{
		PrintStream foutgscore = new PrintStream(new FileOutputStream(folder+"/"+"gscore.gp_gistic.gistic.txt"));
		foutgscore.print("Type	Chromosome	Start	End	-log10(q-value)	G-score	average amplitude	frequency");
		PrintStream foutpval = new PrintStream(new FileOutputStream(folder+"/"+"pval.gp_gistic.gistic.txt"));
		foutpval.print("Type	Chromosome	Start	End	-log10(q-value)	G-score	average amplitude	frequency");
		dealWithFileMarkerFile(folder+"/"+ fileName,0,foutgscore,foutpval);		
		dealWithFileMarkerFile(folder+"/"+ fileName,1,foutgscore,foutpval);
		foutgscore.close();
		foutpval.close();	
	}
	
	public void dealWithFileMarkerFile(String filePath,int type,PrintStream foutgscore,PrintStream foutpval) throws Exception
	{
		String typeString=type==0?"Amp":"Del";
		BufferedReader fin=new BufferedReader(new FileReader(filePath));
		int currentchromosome=0;
		int startlocation=0;
		int currentgscore=0;
		float currentqval=0;
		int lastlocation=0;
		float qval=0;
	    int gscore=0;
		String line=fin.readLine();
		line=fin.readLine();		
		while(line!=null)
		{
			String[] tokens=line.split("\t");			
			String[] loc=tokens[1].split(":");
			line=fin.readLine();
			int chromosome=Integer.parseInt(loc[0].substring(3));
			int location=Integer.parseInt(loc[1]);
			if(type==0)
			{
				qval=Float.parseFloat(tokens[2]);
			    gscore=Integer.parseInt(tokens[3]);
			}else
			{
				if(tokens.length<5)
				{
					line=fin.readLine();
					continue;
				}
				qval=Float.parseFloat(tokens[4]);
				gscore=Integer.parseInt(tokens[5]);				
			}
			
			
			if(currentchromosome==0)
			{
				currentchromosome=chromosome;
				startlocation=location;
				lastlocation=location;
				currentgscore=gscore;
				currentqval=qval;
			}
			
			if(currentchromosome==chromosome && currentgscore==gscore)
			{
				lastlocation=location;
			}			
			else
			{				
				printline(foutgscore,typeString,currentchromosome, startlocation, lastlocation, currentgscore,currentgscore);
				printline(foutpval,typeString,currentchromosome, startlocation, lastlocation, Math.log10(currentqval),currentgscore);
				
				currentchromosome=chromosome;
				startlocation=location;
				lastlocation=location;
				currentgscore=gscore;
				currentqval=qval;				
			}								
			if(line==null)
			{				
				printline(foutgscore,typeString,currentchromosome, startlocation, lastlocation, currentgscore,currentgscore);
				printline(foutpval,typeString,currentchromosome, startlocation, lastlocation, -Math.log10(currentqval),currentgscore);
			}			
		}
		fin.close();
	}
	
	void printline(PrintStream fout,String type,int chrom, int startloc, int endloc, double qval,int gscore)
	{
		fout.print("\n");
		fout.print(type+"\t");
		fout.print(chrom+"\t");
		fout.print(startloc+"\t");
		fout.print(endloc+"\t");
		fout.print(qval+"\t");
		fout.print(gscore+"\t");
		fout.print(0+"\t");							
		fout.print(0);
	}

}
