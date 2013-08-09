package JISTIC;
import java.util.*;
import java.io.*;

public class Genome extends Vector<Gene> {

	static Map<String, Gene> genes_by_name=null;
	static class Band extends GenomeLocation {
		String ID;	

		boolean fake;
		Band(String line) {
			super(Gene.translateChr(line.split("\\t")[0].substring(3)),
					Integer.parseInt(line.split("\\t")[1])+1,
					Integer.parseInt(line.split("\\t")[2]));
			ID = line.split("\\t")[3].replace('.',':');	    
			fake = false;
		}

		/** constructor creating a fake band, covering locations not included
		 *  in the bands file
		 */
		Band(int chr, int start, int end){
			super(chr, start, end);
			ID = start == 0 ? "pter" : "qter";
			fake = true;
		}


		int Arm() { return Character.getNumericValue(ID.charAt(0))
				- Character.getNumericValue('p'); }

		public String toString()
		{ return String.valueOf(chromosome) + ID; }

		boolean SameArm(Band o) {
			return chromosome == o.chromosome && ID.charAt(0) == o.ID.charAt(0);
		}
	}

	public Genome(String locfile, String bands_file)
			throws IOException
			{
		TreeSet<Gene> genes = new TreeSet<Gene>();
		genes_by_name = new HashMap<String, Gene>();
		if(locfile != null){
			boolean parsingError=false;
			String fileErrorLogName=Distribution.dir + (new File(locfile)).getName() +"_error.log";
			PrintStream fileErrorLog = new PrintStream(new FileOutputStream(fileErrorLogName));
			BufferedReader infile = new BufferedReader(new FileReader(locfile));
			for(String next = infile.readLine(); next != null; next = infile.readLine())
			{
				try {
					if(next.startsWith("Centrosome"))
						continue;
					Gene gene = new Gene(next);
					if(genes.contains(gene))
						System.err.println("Duplicate genes in " + locfile + ": " + gene);
					else if(genes_by_name.containsKey(gene.ID)){
						System.err.println("Duplicate genes with inconsistent location in " + locfile + ": "
								+ gene + " and " + genes_by_name.get(gene.ID));
						System.exit(1);
					} else {
						// If we have gene in the exact same location but
						// with a different name, print a warning but
						// still include both genes
						if(!genes.tailSet(gene).isEmpty() && genes.tailSet(gene).first().same_location(gene))
							System.err.println("Warning: duplicate genes with inconsistent location in " + locfile + ": "
									+ gene + " and " + genes.tailSet(gene).first()
									+ "; both will be included in the output");
						if(!genes.headSet(gene).isEmpty() && genes.headSet(gene).last().same_location(gene))
							System.err.println("Warning: duplicate genes with inconsistent location in " + locfile + ": "
									+ gene + " and " + genes.headSet(gene).last()
									+ "; both will be included in the output");
						genes.add(gene);
						genes_by_name.put(gene.ID, gene);
					}
				} catch (NumberFormatException e) {
					parsingError=true;
					fileErrorLog.println(next);	    	
				}
			}
			fileErrorLog.close();
			if(parsingError)
			{
				System.err.println("Some genome locations could not be parsed in " +locfile);
				System.err.println("For details a log file with unparsed locations has been created: " +fileErrorLogName);	    			
			}
		}
		addAll(genes);

		if(bands_file != null){
			SignificantRegion.broad_threshold=new int[3][Marker.numchrom];
			TreeSet<Band> band_set = new TreeSet<Band>();
			// Track, for each chromosome, the lowest and highest location covered in the bands file
			int[] lowest_covered = new int[num_chromosomes], highest_covered = new int[num_chromosomes];
			Arrays.fill(lowest_covered, Integer.MAX_VALUE);
			Arrays.fill(highest_covered, -1);
			BufferedReader infile = new BufferedReader(new FileReader(bands_file));
			for(String next = infile.readLine(); next != null; next = infile.readLine())
			{
				Genome.Band band = new Band(next);
				band_set.add(band);
				lowest_covered[band.chromosome] = Math.min(lowest_covered[band.chromosome], band.start);
				highest_covered[band.chromosome] = Math.max(highest_covered[band.chromosome], band.end);
			}
			bands = new Vector<Band>(band_set);
			for(ListIterator<Band> iter = bands.listIterator(); iter.hasNext(); ){
				Band band = iter.next();
				if(!iter.hasNext() || band.Arm() != bands.get(iter.nextIndex()).Arm())
					if(band.Arm() == 0)
						SignificantRegion.broad_threshold[1][band.chromosome-1] = band.end / 2;
					else {
						SignificantRegion.broad_threshold[2][band.chromosome-1] =
								band.end / 2 - SignificantRegion.broad_threshold[1][band.chromosome-1];
						SignificantRegion.broad_threshold[0][band.chromosome-1] =band.end / 4;
					}
			}

			// For each chromosome, add a qter band covering everything
			// after the last band in the file., anda pter band if needed
			for(int i=0; i < num_chromosomes; i++)
				if(highest_covered[i] >= 0){
					band_set.add(new Band(i, highest_covered[i]+1, Integer.MAX_VALUE));
					if(lowest_covered[i] > 0)
						band_set.add(new Band(i, 0, lowest_covered[i]-1));
				}
			bands = new Vector<Band>(band_set);
		}
			}

	int[] ArmLimits(ListIterator<Marker> iter)
	{
		int[] results = new int[2];
		int marker_band = Collections.binarySearch(bands, iter.previous());
		if(marker_band < 0)
			marker_band = - marker_band - 2;
		int start_band = marker_band; 
		while(start_band > 0 && bands.get(start_band-1).SameArm(bands.get(marker_band)))
			start_band--;
		GenomeLocation arm_start = bands.get(start_band);
		while(iter.hasPrevious() && iter.previous().compareTo(arm_start) >= 0);
		if(iter.hasNext() && iter.next().compareTo(arm_start) >= 0)
			iter.previous();
		results[0] = iter.nextIndex();
		int end_band = marker_band; 
		while(end_band < bands.size()-1 && bands.get(end_band+1).SameArm(bands.get(marker_band)))
			end_band++;
		GenomeLocation arm_end = bands.get(end_band).end_location();
		while(iter.hasNext() && iter.next().compareTo(arm_end) <= 0);
		if(iter.hasPrevious() && iter.previous().compareTo(arm_end) <= 0)
			iter.next();
		results[1] = iter.previousIndex();
		return results;
	}

	List<Band> bands = null;

	// Just in case we ever run this for other organisms that have more chromosomes
	static final int num_chromosomes = 50;

	public Gene byName(String name)
	{
		return genes_by_name.get(name);
	}
}
