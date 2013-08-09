package JISTIC;
import java.util.*;
import java.io.*;

public class combine_matrices {

    public static void main(String[] args){
      try {
	List<String> column_names = new Vector<String>();
	List<String> row_names = new Vector<String>();
	List< Set<Integer> > ones = new Vector< Set<Integer> >();
	for(int file = 0; file < args.length; file++){
		File f=new File(args[file]);
		if(!f.exists())
		{
			System.err.println("Waiting for file to be available: "+args[file]);
			while(!f.exists());
			System.err.println("File is available: "+args[file]);
		}
	    BufferedReader infile = new BufferedReader(new FileReader(f));
	    String[] columns = infile.readLine().split("\\t");
	    if(columns.length == 0)
		continue;
	    for(String next = infile.readLine(); next != null; next = infile.readLine()){
		String[] tokens = next.split("\\t");
		row_names.add(tokens[0]);
		Set<Integer> row_ones = new HashSet<Integer>();
		for(int i=1; i < tokens.length; i++)
		  if(tokens[i].equals("1"))
		    row_ones.add(i + column_names.size() - 1);
		ones.add(row_ones);
	    }
	    column_names.addAll(Arrays.asList(columns).subList(1, columns.length));
	    infile.close();
	}
	String header_line = "";
	for(Iterator<String> iter = column_names.iterator(); iter.hasNext(); )
	    header_line += "\t" + iter.next();
	System.out.println(header_line);
	for(ListIterator<String> iter = row_names.listIterator(); iter.hasNext(); )
	{
	    String line = iter.next();
	    for(int i=0; i < column_names.size(); i++)
	      if(ones.get(iter.previousIndex()).contains(i))
		line += "\t1";
	      else
		line += "\t0";
	    System.out.println(line);
	}
      } catch (Exception e) {
	e.printStackTrace();
      }
    }
}
