package enrichmentanalysis;

import shared.*;
/**
 * LOHcate --- A software tool for LOH calling and visualization in cancer genomes
 * D Wheeler & SG Reddy
 * Human Genome Sequencing Center, Baylor College of Medicine (Houston, TX)
 * 
 * @author Siddharth G. Reddy
 *
 */
public class Gene {
	
	public float rank, rank_2;
	private String name;
	private PathwayElement pathElem;
	
	public Gene(String name) {
		rank = 0;
		rank_2 = 0;
		this.name = Utils.gClean(name);
		pathElem = new PathwayElement("default");
	}
	
	public void setPathElem(Pathway param, String inDir) {
		String source = "";
		try {
			source = FileOps.loadFromFile(inDir);
		} catch (Exception e) { e.printStackTrace(); }
		pathElem = new PathwayElement(param.getDiagram(), source, name);
	}
	
	public void setName(String temp) { name = temp; }
	public String getName() { return name; }
	
	public String toString() { return getName(); }
	public PathwayElement getPathElem() { return pathElem; }
	
}