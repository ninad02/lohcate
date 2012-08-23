package enrichmentanalysis;

import java.awt.Color;
import java.awt.Graphics2D;
import java.awt.image.BufferedImage;
/**
 * LOHcate --- A software tool for LOH calling and visualization in cancer genomes
 * D Wheeler & SG Reddy
 * Human Genome Sequencing Center, Baylor College of Medicine (Houston, TX)
 * 
 * @author Siddharth G. Reddy
 *
 */
public class Rectangle {
	public int x, y, width, height;
	private Float rank, rank_2;
	private Color rank_col, rank_2_col; 
	
	public Rectangle() {
		x = 0; y = 0; width = 0; height = 0;
	}
	public Rectangle(int x, int y, int width, int height) {
		this.x = x; this.y = y; this.width = width; this.height = height;
	}
	
	/**
	 * Fills object-defined bounds with rank indicator color
	 */
	public BufferedImage writeToDiagram(BufferedImage img) {
		Graphics2D obj = img.createGraphics();
		Color fill = new Color(img.getRGB(this.x + 2, this.y + 2));
		for (int x = this.x; x<this.x + width; x++) { //iterate through pixels within bounds
			for (int y = this.y; y<this.y + height; y++) {
				if (img.getRGB(x, y)==(fill).getRGB())// && x < this.x + width/2)
					img.setRGB(x, y, (rank_col).getRGB()); //change to rank indicator color
				//else if (img.getRGB(x, y)==(fill).getRGB()) { // && x >= this.x + width/2)
				
				//}
			}
		}
		
		return img;
	}
	public boolean equals(Rectangle param) { return (param.x==x && param.y==y && param.width==width && param.height==height); }
	public void setRankColor(Color param) { rank_col = param; }
	public Color getRankColor() { return rank_col; }
	public void setRank2Color(Color param) { rank_2_col = param; }
	public void setRank(Float param) { rank = param; }
	public void setRank2(Float param) { rank_2 = param; }
	public Float getRank() { return rank; }
	public Float getRank2() { return rank_2; }
}
