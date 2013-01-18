import java.io.File;
import java.util.ArrayList;

import shared.FileOps;
import shared.Utils;

/**
 * LOHcate --- A software tool for LOH calling and visualization in cancer genomes
 * D Wheeler & SG Reddy
 * Human Genome Sequencing Center, Baylor College of Medicine (Houston, TX)
 * 
 * Here's some garbage code that can be used to parse .copynum & .rawcnv files into BED files (that can be uploaded to UCSC genome browser)
 * 
 * CNV calling from Affy .CEL files --> http://www.openbioinformatics.org/penncnv/penncnv_tutorial_affy_gw6.html
 * 
 * @author Siddharth G. Reddy
 *
 */
public class CopyNum {

	public static void addCNVTrack(String inDir, String outDir) { //[other: hepato, ovarian, burkitts]
		String[] split;
		String toWrite = "", chr, load;
		ArrayList<Point> regions;
		int loc_min = Integer.MAX_VALUE, loc_max = Integer.MIN_VALUE;
		float min, max;
		File[] files = (new File(inDir)).listFiles();
		for (int i = 1; i<=22; i++) { //iterate through chromosomes
			if (i<=21)
				chr = Integer.toString(i);
			else
				chr = "X";
			System.out.println("chr" + chr);
			for (File file : files) {
				if (file.getName().indexOf(".txt")!=-1) {
					System.out.print(file.getName().split(".txt")[0] + " ");
					load = FileOps.loadFromFile(file.getAbsolutePath());
					split = load.split("\tchr" + chr + "\t"); //sans-"chr" for hepato
					try {
						toWrite = "";
						regions = new ArrayList<Point>();
						min = Float.MAX_VALUE; max = -1 * Float.MAX_VALUE;
						loc_min = Integer.MAX_VALUE; loc_max = 0;
						
						for (int j = 1; j<split.length; j++) { //iterate through regions -- grab bounds
							if (Float.parseFloat(split[j].split("\t")[1].split("\n")[0]) > max)
								max = Float.parseFloat(split[j].split("\t")[1].split("\n")[0]);
							if (Float.parseFloat(split[j].split("\t")[1].split("\n")[0]) < min)
								min = Float.parseFloat(split[j].split("\t")[1].split("\n")[0]);
							
							if (Integer.parseInt(split[j].split("\t")[0]) < loc_min)
								loc_min = Integer.parseInt(split[j].split("\t")[0]);
							if (Integer.parseInt(split[j].split("\t")[0]) > loc_max)
								loc_max = Integer.parseInt(split[j].split("\t")[0]);
						}
						
						for (int j = 1; j<split.length; j++) { //iterate through regions -- collect/merge
							if (regions.size()==0)
								regions.add(new Point(Integer.parseInt(split[j].split("\t")[0]), Integer.parseInt(split[j].split("\t")[0])
										, Float.parseFloat(split[j].split("\t")[1].split("\n")[0]) - min));
							else
								regions = cnvSplice(new Point(Integer.parseInt(split[j].split("\t")[0]), Integer.parseInt(split[j].split("\t")[0])
										, Float.parseFloat(split[j].split("\t")[1].split("\n")[0]) - min), regions);
						}
						//System.out.println(split.length + " | " + regions.size());
						for (Point region : regions) //iterate through regions -- dump
							FileOps.appendToFile(outDir + "/chr" + chr + ".txt", "chr" + chr + " " + region.x + " " + region.y + " cn" + (region.score / (float)region.z) + " " + Utils.normalize((float)Math.pow(2, Math.abs(region.score / (float)region.z)), (float)Math.pow(2, 0), (float)Math.pow(2, Math.abs(max - min)), 200, 900) + " + " + region.x + " " + region.y + "\n");
						
					} catch (Exception e) { e.printStackTrace(); System.out.println("throwing chr" + chr); }
				}
			}
			System.out.println();
			toWrite = "browser position chr" + chr + ":" + loc_min + "-" + loc_max + "\nbrowser hide all\n"
				+ "track name=\"CopyNumVar\" description=\" \" visibility=dense useScore=1\n";
			FileOps.writeToFile(outDir + "/chr" + chr + ".txt", toWrite + FileOps.loadFromFile(outDir + "/chr" + chr + ".txt"));
		}
	}
	
	/*public static void addCNVTrack(String inDir, String outDir) { //varscan-cnv
		String[] split;
		String toWrite = "", chr, load;
		ArrayList<Point> regions;
		int loc_min = Integer.MAX_VALUE, loc_max = Integer.MIN_VALUE;
		float min, max;
		File[] files = (new File(inDir)).listFiles();
		for (int i = 1; i<=22; i++) { //iterate through chromosomes
			if (i<=21)
				chr = Integer.toString(i);
			else
				chr = "X";
			System.out.println("chr" + chr);
			for (File file : files) {
				if (file.getName().indexOf(".copynumber")!=-1) {
					System.out.print(file.getName().split(".copynumber")[0] + " ");
					load = FileOps.loadFromFile(file.getAbsolutePath());
					split = load.split("\n" + chr + "\t");
					try {
						toWrite = "";
						regions = new ArrayList<Point>();
						min = Float.MAX_VALUE; max = -1 * Float.MAX_VALUE;
						loc_min = Integer.MAX_VALUE; loc_max = 0;
						
						for (int j = 1; j<split.length; j++) { //iterate through regions -- grab bounds
							if (Float.parseFloat(split[j].split("\t")[5]) > max)
								max = Float.parseFloat(split[j].split("\t")[5]);
							if (Float.parseFloat(split[j].split("\t")[5]) < min)
								min = Float.parseFloat(split[j].split("\t")[5]);
							
							if (Integer.parseInt(split[j].split("\t")[0]) < loc_min)
								loc_min = Integer.parseInt(split[j].split("\t")[0]);
							if (Integer.parseInt(split[j].split("\t")[1]) > loc_max)
								loc_max = Integer.parseInt(split[j].split("\t")[1]);
						}
						
						for (int j = 1; j<split.length; j++) { //iterate through regions -- collect/merge
							if (regions.size()==0)
								regions.add(new Point(Integer.parseInt(split[j].split("\t")[0]), Integer.parseInt(split[j].split("\t")[1])
										, Float.parseFloat(split[j].split("\t")[5]) - min));
							else
								regions = cnvSplice(new Point(Integer.parseInt(split[j].split("\t")[0]), Integer.parseInt(split[j].split("\t")[1])
										, Float.parseFloat(split[j].split("\t")[5]) - min), regions);
						}
						//System.out.println(split.length + " | " + regions.size());
						for (Point region : regions) //iterate through regions -- dump
							FileOps.appendToFile(outDir + "/chr" + chr + ".txt", "chr" + chr + " " + region.x + " " + region.y + " cn" + (region.score / (float)region.z) + " " + normalize((float)Math.pow(2, Math.abs(region.score / (float)region.z)), (float)Math.pow(2, 0), (float)Math.pow(2, Math.abs(max - min)), 200, 900) + " + " + region.x + " " + region.y + "\n");
						
					} catch (Exception e) { e.printStackTrace(); System.out.println("throwing chr" + chr); }
				}
			}
			System.out.println();
			toWrite = "browser position chr" + chr + ":" + loc_min + "-" + loc_max + "\nbrowser hide all\n"
				+ "track name=\"CopyNumVar\" description=\" \" visibility=dense useScore=1\n";
			FileOps.writeToFile(outDir + "/chr" + chr + ".txt", toWrite + FileOps.loadFromFile(outDir + "/chr" + chr + ".txt"));
		}
	}*/
	
	/*public static void addCNVTrack(String inDir, String outDir) { //penncnv
		String[] split;
		String toWrite, chr, load = FileOps.loadFromFile(inDir);
		ArrayList<Point> regions;
		int min, max, loc_min, loc_max;
		for (int i = 1; i<=22; i++) { //iterate through chromosomes
			if (i<=21)
				chr = Integer.toString(i);
			else
				chr = "X";
			split = load.split("chr" + chr + ":");
			System.out.println("chr" + chr);
			try {
				toWrite = "";
				regions = new ArrayList<Point>();
				min = Integer.MAX_VALUE; max = 0;
				loc_min = Integer.MAX_VALUE; loc_max = 0;
				for (int j = 1; j<split.length; j++) { //iterate through regions -- collect/merge, grab min/max
					if (Integer.parseInt(split[j].split("cn=")[1].split(" ")[0]) > max)
						max = Integer.parseInt(split[j].split("cn=")[1].split(" ")[0]);
					if (Integer.parseInt(split[j].split("cn=")[1].split(" ")[0]) < min)
						min = Integer.parseInt(split[j].split("cn=")[1].split(" ")[0]);
					
					if (Integer.parseInt(split[j].split("-")[0]) < loc_min)
						loc_min = Integer.parseInt(split[j].split("-")[0]);
					if (Integer.parseInt(rmTrailingSpaces(split[j].split("numsnp=")[0].split("-")[1])) > loc_max)
						loc_max = Integer.parseInt(rmTrailingSpaces(split[j].split("numsnp=")[0].split("-")[1]));
					
					if (regions.size()==0)
						regions.add(new Point(Integer.parseInt(split[j].split("numsnp=")[0].split("-")[0]), Integer.parseInt(rmTrailingSpaces(split[j].split("numsnp=")[0].split("-")[1]))
								, Float.parseFloat(split[j].split("cn=")[1].split(" ")[0])));
					else
						regions = cnvSplice(new Point(Integer.parseInt(split[j].split("numsnp=")[0].split("-")[0]), Integer.parseInt(rmTrailingSpaces(split[j].split("numsnp=")[0].split("-")[1]))
								, Float.parseFloat(split[j].split("cn=")[1].split(" ")[0])), regions);
				}
				//System.out.println(split.length + " | " + regions.size());
				for (Point region : regions) //iterate through regions -- dump
					toWrite += "chr" + chr + " " + region.x + " " + region.y + " cn" + (int)(region.score / (float)region.z) + " " + normalize(region.score / (float)region.z, min, max, 200, 900) + " + " + region.x + " " + region.y + "\n";
				
				if (toWrite.length() > 0) {
					toWrite = "browser hide all\n"
						+ "track name=\"CopyNumVar\" description=\" \" visibility=dense useScore=1\n" + toWrite;
					toWrite = "browser position chr" + chr + ":" + loc_min + "-" + loc_max + "\n" + toWrite;
					FileOps.writeToFile(outDir + "/chr" + chr + ".txt", toWrite);
				}
				
			} catch (Exception e) { e.printStackTrace(); System.out.println("throwing chr" + chr); }
		}
	}*/
	
	public static ArrayList<Point> cnvSplit(Point r, ArrayList<Point> arr) { //splice, w/out insertion at end
		ArrayList<Point> rtn = new ArrayList<Point>();
		for (Point elem : arr) {
			if (!OldCode.intersecting(r, elem))
				rtn.add(elem);
			else {
				if (r.x - elem.x > 0)
					rtn.add(new Point(elem.x, r.x, elem.score)); //add non-zero left region
				if (elem.y - r.y > 0)
					rtn.add(new Point(r.y, elem.y, elem.score)); //add non-zero right region
				rtn.add(new Point((int)Math.max(elem.x, r.x), (int)Math.min(elem.y, r.y), elem.z + 1, elem.score + r.score));
			}
		}
		return rtn;
	}
	
	public static ArrayList<Point> cnvSplice(Point r, ArrayList<Point> arr) {
		ArrayList<Point> rtn = new ArrayList<Point>(), toAdd = new ArrayList<Point>();
		toAdd.add(r);
		for (Point elem : arr) { //iterate through 'template' arr
			if (OldCode.intersecting(elem, r)) { //if takeover
				if (r.x - elem.x > 0)
					rtn.add(new Point(elem.x, r.x, elem.score)); //add non-zero left region
				if (elem.y - r.y > 0)
					rtn.add(new Point(r.y, elem.y, elem.score)); //add non-zero right region
				cnvSplit(elem, toAdd);
			}
			else
				rtn.add(elem);
		}
		for (Point elem : toAdd)
			rtn.add(elem);
		return rtn;
	}
	
}
