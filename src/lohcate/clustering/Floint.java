package lohcate.clustering;

import java.util.Comparator;

/**
 * An interface that represents a point.  This allows for any class to represent a 
 * point, enabling the direct use of the class for geometric purposes rather than
 * creating new point objects.  Furthermore, this prevents the disconnect between
 * a geometric point and its source object.
 * 
 * @author Ninad Dewal
 * 
 */
public interface Floint {
	
	// ========================================================================
	public static final Comparator<Floint> XYComparator = new Comparator<Floint>() {
		public int compare(Floint f1, Floint f2) {
			int result = Float.compare(f1.getX(), f2.getX());
			if (result == 0) {
				result = Float.compare(f1.getY(), f2.getY());
			}
			return result;
		}
	};	
	
	// ========================================================================
	public static final Comparator<Floint> YXComparator = new Comparator<Floint>() {
		public int compare(Floint f1, Floint f2) {
			int result = Float.compare(f1.getY(), f2.getY());
			if (result == 0) {
				result = Float.compare(f1.getX(), f2.getX());
			}
			return result;
		}
	};	
	
	// ========================================================================
	public static final Comparator<Floint> XYZComparator = new Comparator<Floint>() {
		public int compare(Floint f1, Floint f2) {
			int result = XYComparator.compare(f1, f2);
			if (result == 0) {
				result = Float.compare(f1.getZ(), f2.getZ());
			}
			return result;
		}
	};	

	// ========================================================================
	public float getX();
	public float getY();
	public float getZ();
	
	// ========================================================================
	public <E extends Floint> double getCartesianDistanceSquared(E rhs);

	// ========================================================================
	public <E extends Floint> double getCartesianDistance(E rhs);
	
	// Cartesian distance is given by:
	// dist = sqrt( (x2 - x1)^2 + (y2 - y1)^2 )
	// Thus, the maximum possible difference on x (given a fixed distance) must assume that (y2 - y1) = 0
	// The equation for calculating this difference is then:
	// dist = sqrt ( (x2 - x1)^2 + 0^2 ) = sqrt ((x2 - x1)^2) = x2 - x1
	// The same applies for y: dist = y2 - y1
	public double getTheoreticalFurthestDifferenceXOrYWithinDistance(double distance);
	
	// ========================================================================
	// Global Comparator object for public use	
	public static interface FlointDistance {
		public double distance(Floint n1, Floint n2); 
	}
	
	public static final FlointDistance FlointDistanceX = new FlointDistance() {
		public double distance(Floint n1, Floint n2) { return n1.getX() - n2.getX(); }
	};
	
	public static final FlointDistance FlointDistanceY = new FlointDistance() {
		public double distance(Floint n1, Floint n2) { return n1.getY() - n2.getY(); }
	};
		
	
}
