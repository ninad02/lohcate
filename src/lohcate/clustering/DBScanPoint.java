package lohcate.clustering;

import java.util.Comparator;

/** This class contains meta information for each point in DBScan in order
 *  to decrease run-time for the DBScan algorithm.
 *  
 *  @author Ninad Dewal
 */
public class DBScanPoint /*implements Comparable<DBScanPoint>*/ {
	
	// ========================================================================
	private boolean mVisited;
	public boolean mAdded;
	public int mClusterAssigned;

	// ========================================================================
	public int mIndexSortedX;
	public int mIndexSortedY;		
	
	// ========================================================================
	protected Floint mFloint; 
	
	// ========================================================================
	protected DBScanPoint(Floint floint) {		
		mFloint = floint;
		reset();
	}
	
	// ========================================================================
	public boolean getVisited() { return mVisited; }
	public void    setVisited() { mVisited = true; }	
	public void  resetVisited() { mVisited = false; }
	
	// ========================================================================
	public void reset() {
		resetForClustering();
		mIndexSortedX = mIndexSortedY = -1;
	}
	
	// ========================================================================
	public void resetForClustering() {
		resetVisited();
		mClusterAssigned = DBSCAN2.ClusterIDOfNoise;
		mAdded = false;
	}
	
	// ========================================================================
	public double getCartesianDistance(DBScanPoint rhs) {
		return mFloint.getCartesianDistance(rhs.mFloint);
	}
	
	// ========================================================================
	public double getCartesianDistanceSquared(DBScanPoint rhs) {
		return mFloint.getCartesianDistanceSquared(rhs.mFloint);
	}
	
	// ========================================================================
	/** Default sorting by x-coordinate */
	/*public int compareTo(DBScanPoint rhs) {
		return mFloint.compareTo(rhs.mFloint);
	}*/

	// ========================================================================
	// Global Comparator object for public use	
	public static final Comparator<DBScanPoint> DBScanPointCompareXY = new Comparator<DBScanPoint>() {
		public int compare(DBScanPoint p1, DBScanPoint p2) {
			return Floint.XYComparator.compare(p1.mFloint, p2.mFloint);
		}
	};
	
	public static final Comparator<DBScanPoint> DBScanPointCompareYX  = new Comparator<DBScanPoint>() {
		public int compare(DBScanPoint p1, DBScanPoint p2) {
			return Floint.YXComparator.compare(p1.mFloint, p2.mFloint);
		}
	};
	
	
	// ====================================================================
	public static interface DBScanPointDistance {
		public double distance(DBScanPoint n1, DBScanPoint n2); 
	}
	
	public static final DBScanPointDistance DBScanPointDistanceX = new DBScanPointDistance() {
		public double distance(DBScanPoint n1, DBScanPoint n2) { 
			return Floint.FlointDistanceX.distance(n1.mFloint, n2.mFloint); 
		} 
	};
	
	public static final DBScanPointDistance DBScanPointDistanceY = new DBScanPointDistance() {
		public double distance(DBScanPoint n1, DBScanPoint n2) { 
			return Floint.FlointDistanceY.distance(n1.mFloint, n2.mFloint); 
		}
	};
}