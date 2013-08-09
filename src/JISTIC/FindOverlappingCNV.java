package JISTIC;

import java.io.BufferedReader;
import java.io.FileReader;

public class FindOverlappingCNV {
	public static void main(String[] args) throws Exception {
		boolean mydata=false;
		Genome gene_genome = null;
		int chrom=-1;
		int start=-1;
		int end=-1;		
		String name="";
		String regFilesString=null;
		String locFilesString=null;		
		
		if(args.length==1 && args[0].equalsIgnoreCase("local"))
			mydata=true;
		
		if(!mydata)
		{
			if(args.length==4)
			{
				chrom=Integer.parseInt(args[0]);
				start=Integer.parseInt(args[1]);
				end=Integer.parseInt(args[2]);
				regFilesString=args[3];
			}else if(args.length==3)
			{
				locFilesString=args[0];
				regFilesString=args[2];
				gene_genome=new Genome(args[1], null);
			}else if(args.length==2)
			{
				locFilesString=args[0];
				regFilesString=args[1];
			}else
			{
				System.out.println("Parameters option 1: chromomsome start end comma_sep_cnv_files");
				System.out.println("Parameters option 2: inputfile comma_sep_cnv_files");
				System.out.println("input file should have format 'geneid\tchrom\tstart\tend' with no header");
				System.out.println("Parameters option 3: inputfilegeneids genomefile comma_sep_cnv_files");							
				System.exit(0);
			}
			
		}else
		{			
			chrom=6;
			start=19943714;
			end=19951655;
			regFilesString="./ovarian/CNV.verified_080606.combined.mod.txt,./ovarian/CNV.XY.txt,./ovarian/10-01-10_Ovarian_CNV_Dylan.txt,./ovarian/May10OvarianBatchCVN1.txt,./ovarian/May10OvarianUpDownCVN1.txt,./ovarian/May10OvarianNormalsCVN1.txt";
			//regFilesString="./ovarian/CNV.verified_080606.combined.mod.txt,./ovarian/CNV.XY.txt,./ovarian/10-01-10_Ovarian_CNV_Dylan.txt,./ovarian/May10OvarianBatchCVN1.txt,./ovarian/May10OvarianUpDownCVN1.txt,./ovarian/May10OvarianNormalsCVNTest.txt";						
			//locFilesString="C:/MyDocs/Docs/Columbia/subject/CompBio/Ovarian/keygeneslocations.txt";
			//locFilesString="C:/MyDocs/Docs/Columbia/subject/CompBio/Ovarian/keygenes.txt";
			//locFilesString="C:/Documents and Settings/fs2282/My Documents/subject/CompBio/Ovarian/keygenes.txt";
			//gene_genome=new Genome("./may10_hg18_Gene_Info.txt", null);
		}
		
		BufferedReader fin1=null;
		String line1=null;
		if(locFilesString!=null)
		{
			fin1=new BufferedReader(new FileReader(locFilesString));
			line1=fin1.readLine();
		}
		
		while(line1!=null || chrom!=-1)
		{
			if(line1!=null)
			{
				String[] tokens=line1.split("\t");
				if(gene_genome==null)
				{
					name=tokens[0];
					chrom=Integer.parseInt(tokens[1]);
					start=Integer.parseInt(tokens[2]);
					end=Integer.parseInt(tokens[3]);
				}else
				{
					name=tokens[0];
					Gene gene=gene_genome.byName(name);
					chrom=gene.chromosome;
					start=gene.start;
					end=gene.end;	
					name=gene.symbol+"("+name+")";					
				}
			}
			System.out.println("\n"+name+"\t"+chrom+":"+start+"-"+end);
			String[] regFiles=regFilesString.split(",");
			for (int i = 0; i < regFiles.length; i++) {
				int j=0;
				BufferedReader fin=new BufferedReader(new FileReader(regFiles[i]));
				String line=fin.readLine();
				line=fin.readLine();
				while(line!=null)
				{
					j++;
					String[] tokens=line.split("\t");
					if(tokens.length>3)
					{
						int cnvchrom=-1;
						if(tokens[1].equalsIgnoreCase("X"))
						{
							cnvchrom=23;
						}else if(tokens[1].equalsIgnoreCase("Y"))
						{
							cnvchrom=24;
						}else
						{
							cnvchrom=Integer.parseInt(tokens[1]);
						}
						
						if(chrom==cnvchrom)
						{
							int cnvstart=Integer.parseInt(tokens[2]);
							int cnvend=Integer.parseInt(tokens[3]);
							if((cnvstart<=start && cnvend>=start) || // empieza antes
									(cnvstart<=end && cnvend>=end) ||  // empieza despues
									(cnvstart>=start && cnvstart<=end) ) // contenida dentro
							{
								System.out.println("\t"+regFiles[i]+" line "+j+" :");
								System.out.println("\t"+line);
							}
						}
					}
					line=fin.readLine();
				}			
				fin.close();
			}
			chrom=-1;
			if(line1!=null)
				line1=fin1.readLine();
		}
	}

}
