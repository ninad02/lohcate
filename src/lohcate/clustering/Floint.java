package lohcate.clustering;
import java.util.Comparator;
import java.util.GregorianCalendar;


/**
 * LOHcate --- A software tool for LOH calling and visualization in cancer genomes
 * D Wheeler & SG Reddy
 * Human Genome Sequencing Center, Baylor College of Medicine (Houston, TX)
 * 
 * 'Point', but with float coordinates + an extra integer
 * 
 * @author Siddharth G. Reddy, Ninad Dewal
 *
 */
public class Floint implements Comparable<Floint> {
	public float mX, mY, mZ;
	
	public Floint(float x, float y, float z) { 
		mX = x; 
		mY = y; 				
		mZ = z;
	}
	
	public double getCartesianDistanceSquared(Floint rhs) {
		return (this.mX - rhs.mX) * 
			   (this.mX - rhs.mX) +
			   
			   (this.mY - rhs.mY) * 
			   (this.mY - rhs.mY) +
			   
			   (Clustering.Doing3D ?
					   (this.mZ - rhs.mZ) *
					   (this.mZ - rhs.mZ) 
					   : 0)
			   ;
	}
	
	public double getCartesianDistance(Floint rhs) {
		return Math.sqrt(getCartesianDistanceSquared(rhs));
	}
	
	// Cartesian distance is given by:
	// dist = sqrt( (x2 - x1)^2 + (y2 - y1)^2 )
	// Thus, the maximum possible difference on x (given a fixed distance) must assume that (y2 - y1) = 0
	// The equation for calculating this difference is then:
	// dist = sqrt ( (x2 - x1)^2 + 0^2 ) = sqrt ((x2 - x1)^2) = x2 - x1
	// The same applies for y: dist = y2 - y1
	public double getTheoreticalFurthestDifferenceXOrYWithinDistance(double distance) {
		return distance;
	}
	
	public int compareTo(Floint rhs) {		
		int rV = Float.compare(mX, rhs.mX);
		if (rV == 0) {
			rV = Float.compare(mY, rhs.mY);
		}
		return rV;
	}
	
	// ========================================================================
	// Global Comparator object for public use	
	public static final FlointCompareX  FlointCompareXObj  = new FlointCompareX();
	public static final FlointCompareY  FlointCompareYObj  = new FlointCompareY();
	
	public static final FlointDistanceX FlointDistanceXObj = new FlointDistanceX();
	public static final FlointDistanceY FlointDistanceYObj = new FlointDistanceY();
	
	// Comparator
	public static class FlointCompareY implements Comparator<Floint> {
		public int compare(Floint n1, Floint n2) {
			int rV = Float.compare(n1.mY, n2.mY);
			if (rV == 0) {
				rV = Float.compare(n1.mX, n2.mX);
			}
			return rV;
		}		
	}
	
	public static class FlointCompareX implements Comparator<Floint> {
		public int compare(Floint n1, Floint n2) {		
			int rV = Float.compare(n1.mX, n2.mX);
			if (rV == 0) {
				rV = Float.compare(n1.mY, n2.mY);
			}
			return rV;
		}
	}
	
	public static abstract class FlointDistance {
		public abstract double distance(Floint n1, Floint n2); 
	}
	
	public static class FlointDistanceX extends FlointDistance {
		public double distance(Floint n1, Floint n2) { return n1.mX - n2.mX; }
	}
	
	public static class FlointDistanceY extends FlointDistance {
		public double distance(Floint n1, Floint n2) { return n1.mY - n2.mY; }
	}
	
	public static void TestPowers() {
		GregorianCalendar gc1a = new GregorianCalendar();  
		int numIter = 100;
		
		System.out.println();
		for (int i = 0; i < numIter; i++) {
			for (int j = 0;  j < numIter; j++) {
				for (int k = 0; k < numIter; k++) {
					for (int m = 0; m < numIter; m++) {
						int diff1 = i - j;
						int diff2 = k - m;
						double distance = Math.sqrt( diff1 * diff1 + diff2 * diff2 );
					}
				}
			}
		}
		GregorianCalendar gc1b = new GregorianCalendar();
				
		for (int i = 0; i < numIter; i++) {
			for (int j = 0;  j < numIter; j++) {
				for (int k = 0; k < numIter; k++) {
					for (int m = 0; m < numIter; m++) {
						double distance = Math.sqrt( Math.pow((i - j), 2) + Math.pow((k - m), 2) );
					}
				}
			}
		}
		GregorianCalendar gc1c = new GregorianCalendar();
		
		System.out.println(gc1a.toString() + "\t" + gc1b.toString() + "\t" + (gc1b.getTimeInMillis() - gc1a.getTimeInMillis()));
		System.out.println(gc1b.toString() + "\t" + gc1c.toString() + "\t" + (gc1c.getTimeInMillis() - gc1b.getTimeInMillis()));
	}
	
	public static void main(String[] args) {
		TestPowers();
	}
}