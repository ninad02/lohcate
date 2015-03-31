package enrichmentanalysis;

import java.awt.Color;
import java.awt.image.BufferedImage;
import shared.FileOps;
/**
 * LOHcate --- A software tool for LOH calling and visualization in cancer genomes
 * D Wheeler & SG Reddy
 * Human Genome Sequencing Center, Baylor College of Medicine (Houston, TX)
 * 
 * @author Siddharth G. Reddy
 *
 */
public class Pathway {
	
	private String name, hsaID;
	private BufferedImage diagram;
	
	public Pathway(String hsaID) {
		this("", hsaID);
	}
	public Pathway(String name, String hsaID) {
		this.name = name;
		this.hsaID = hsaID;
	}
	
	/**
	 * Remove default coloring
	 */
	public static BufferedImage processRawPathway(BufferedImage img) {
		for (int row = 0; row<img.getHeight(); row++)
			for (int col = 0; col<img.getWidth(); col++)
				if (!(new Color(img.getRGB(col, row))).equals(Color.BLACK) && !(new Color(img.getRGB(col, row))).equals(Color.RED) && !(new Color(img.getRGB(col, row))).equals(Color.WHITE))
					img.setRGB(col, row, (Color.WHITE).getRGB());
		return img;
	}
	
	/**
	 * Grab pathway diagram+annotations from KEGG Pathways
	 */
	public void setDiagram(String filePath, String annot_path) {	
		diagram = processRawPathway(FileOps.getImage("http://www.genome.jp/kegg/pathway/hsa/hsa" + hsaID + ".png"));
		FileOps.writeImageToFile(diagram, filePath);
		try {
			FileOps.writeToFile(annot_path, FileOps.getHTML("http://www.genome.jp/kegg/pathway/hsa/hsa" + hsaID + ".html"));
		} catch (Exception e) { e.printStackTrace(); }
	}
	
	public BufferedImage getDiagram() { return diagram; }
	public String getHSAID() { return hsaID; }
	public String getName() { return name; }
	
}