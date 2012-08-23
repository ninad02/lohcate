/**
 * LOHcate --- A software tool for LOH calling and visualization in cancer genomes
 * D Wheeler & SG Reddy
 * Human Genome Sequencing Center, Baylor College of Medicine (Houston, TX)
 * 
 * A pretty straightforward mod of the traditional Point class, which only has fields x,y (int)
 * 
 * @author Siddharth G. Reddy
 *
 */
public class Point {
	public int x, y, z;
	public float score;
	public Point(int x, int y) {
		this(x, y, 0);
	}
	public Point(int x, int y, float score) {
		this(x, y, 1, score);
	}
	public Point(int x, int y, int z, float score) {
		this.x = x; this.y = y; this.score = score; this.z = z;
	}
}
