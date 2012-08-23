package enrichmentanalysis;


import java.awt.Color;

/**
 * LOHcate --- A software tool for LOH calling and visualization in cancer genomes
 * D Wheeler & SG Reddy
 * Human Genome Sequencing Center, Baylor College of Medicine (Houston, TX)
 * 
 * @author Siddharth G. Reddy
 *
 */
public class Cell {
	
	private Color color;
	private float value;
	private int height, width;
	
	public Cell(float value, int width, int height) {
		this.value = value;
		this.width = width;
		this.height = height;
	}
	
	public Cell(int width, int height) {
		this(0f, width, height);
	}
	
	public Cell(Color color, float value) {
		this.color = color;
		this.value = value;
	}
	
	public void setColor(Color param) { color = param; }
	public void setValue(float param) { value = param; }
	
	public Color getColor() { return color; }
	public float getValue() { return value; }
	public int getWidth() { return width; }
	public int getHeight() { return height; }
}
