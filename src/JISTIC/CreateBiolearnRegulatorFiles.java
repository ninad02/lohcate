package JISTIC;

import java.io.BufferedReader;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Iterator;

public class CreateBiolearnRegulatorFiles {
	public static void main(String[] args) throws Exception {
		// input are folders
		String[] abbType={"AMP","DEL"};
		String[] abbLevel={"Broad","NonBroad","Peak"};
		String[] files={"AMP.genes.All.matrix","DEL.genes.All.matrix"};
		String outputfolder=null;
		String[] folders=null;
		try
		{
			outputfolder=args[0];
			folders=args[1].split(",");			
		}catch(Exception e)
		{		
			System.err.println("Format: outputfolder commaseparatedjisticfolders");
			throw new RuntimeException("Problems parsing the input for CreateIGVSessionFile");
		}
		HashSet<String>[][] regs=(HashSet<String>[][]) new HashSet[2][3];
		
		for (int i = 0; i < abbType.length; i++) {
			for (int j = 0; j < abbLevel.length; j++) {						
				regs[i][j]=new HashSet<String>();				
			}
		}
		
		for (int i = 0; i < abbType.length; i++) {
			for (int j = 0; j < folders.length; j++) {
				BufferedReader fin=new BufferedReader(new FileReader(folders[j]+"/"+files[i]));
				String line=fin.readLine();
				// process header
				String[] tokens=line.split("\t");
				int[] regiontype=new int[tokens.length];
				for (int k = 1; k < regiontype.length; k++) {
					if(tokens[k].startsWith("PEAK"))
					{
						regiontype[k]=3;
					}else if(tokens[k].startsWith("BROAD"))
					{
						regiontype[k]=1;
					}else
					{
						regiontype[k]=2;
					}
				}								
				// process each line
				line=fin.readLine();
				while(line!=null)
				{						
					tokens=line.split("\t");
					int maxtype=0;
					for (int k = 1; k < tokens.length; k++) {
						if(tokens[k].equals("1"))
						{
							maxtype=Math.max(maxtype, regiontype[k]);
							if(maxtype==3)
								break;
						}
					}
					for (int k = 0; k < maxtype; k++) {
						regs[i][k].add(tokens[0]);
					}
					line=fin.readLine();
				}																													
				fin.close();
			}
		}
		
		for (int i = 0; i < abbLevel.length; i++) {
			HashSet<String> combinedRegs=new HashSet<String>();
			for (int j = 0; j < abbType.length; j++) {
				// output the list for aberration type j and abberration level i
				PrintStream fout = new PrintStream(new FileOutputStream(outputfolder+"/"+abbLevel[i]+"Genes."+abbType[j]));
				Iterator iterator = regs[j][i].iterator();
				if(iterator.hasNext())
				{
					fout.print(iterator.next());
					while(iterator.hasNext()) 
						fout.print("\n"+iterator.next());					
				}
				fout.close();
				combinedRegs.addAll(regs[j][i]);
			}
			
			// now output the combined list
			PrintStream fout = new PrintStream(new FileOutputStream(outputfolder+"/"+abbLevel[i]+"Genes.ALL"));
			Iterator iterator = combinedRegs.iterator();
			if(iterator.hasNext())
			{
				fout.print(iterator.next());
				while(iterator.hasNext())
					fout.print("\n"+iterator.next());				
			}			
			fout.close();
			
		}
		
	}
}
