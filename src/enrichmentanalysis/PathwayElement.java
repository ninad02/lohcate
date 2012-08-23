package enrichmentanalysis;

import java.awt.Color;
import java.awt.Font;
import java.awt.Graphics2D;
//import java.awt.Rectangle;
import java.awt.image.BufferedImage;
import java.util.ArrayList;

import shared.FileOps;
import shared.Utils;
/**
 * LOHcate --- A software tool for LOH calling and visualization in cancer genomes
 * D Wheeler & SG Reddy
 * Human Genome Sequencing Center, Baylor College of Medicine (Houston, TX)
 * 
 * @author Siddharth G. Reddy
 *
 */
public class PathwayElement {
	
	public final int WIDTH = 46, HEIGHT = 17; //default in KEGG Pathways
	
	private String name;
	private BufferedImage diagram; //blank pathway wiring diagram
	private ArrayList<Rectangle> boxes; //multiple genes/aliases can be assigned to a block in a KEGG Pathways wiring diagram
	
	public PathwayElement(BufferedImage diagram, String source, String gene_name) {
		this.diagram = diagram;
		boxes = new ArrayList<Rectangle>();
		setLoc(source, gene_name);
	}
	public PathwayElement(String name) {
		this.name = Utils.gClean(name);
	}
	/**
	 * Grab coordinates of pathway element location(s) from HTML image-map annotations. 
	 */
	public void setLoc(String source, String gene_name) {
		for (int i = 0; i<source.split("\\(" + gene_name + "\\)").length-1; i++)
			boxes.add(new Rectangle(Integer.parseInt(source.split("\\(" + gene_name + "\\)")[i].split("coords=")[source.split("\\(" + gene_name + "\\)")[i].split("coords=").length - 1].split(",")[0])
							, Integer.parseInt(source.split("\\(" + gene_name + "\\)")[i].split("coords=")[source.split("\\(" + gene_name + "\\)")[i].split("coords=").length - 1].split(",")[1])
							, WIDTH, HEIGHT));
	}

	public ArrayList<Rectangle> getBoxes() { return boxes; }
	public void setDiagram(BufferedImage temp) { diagram = temp; }
	public String getName() { return name; }
	public BufferedImage getDiagram() { return diagram; }
	
	public boolean isClone(PathwayElement potClone) {
		boolean rtn = false;
		for (Rectangle box1 : boxes)
			for (Rectangle box2 : potClone.getBoxes())
				if (box1.equals(box2))
					return true;
		return rtn;
	}
}