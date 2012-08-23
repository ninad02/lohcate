/**
 * LOHcate --- A software tool for LOH calling and visualization in cancer genomes
 * D Wheeler & SG Reddy
 * Human Genome Sequencing Center, Baylor College of Medicine (Houston, TX)
 * 
 * 'Point', but with float coordinates + an extra integer
 * 
 * @author Siddharth G. Reddy
 *
 */
public class Floint {
	public float x, y;
	public int loc;
	public Floint(float x, float y) { this(x, y, -1); }
	public Floint(float x, float y, int loc) { this.x = x; this.y = y; this.loc = loc; }
}