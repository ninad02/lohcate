import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.ListIterator;

/**
 * LOHcate --- A software tool for LOH calling and visualization in cancer genomes
 * @author Ninad Dewal, Siddharth G. Reddy, David Wheeler
 * Human Genome Sequencing Center, Baylor College of Medicine (Houston, TX)
 * 
 * The superclass for faster versions of DBSCAN.  Subclasses can override this method's 
 * getNeighbors() method, which is the speed bottleneck.  
 * 
 * A summary of the DBSCAN algorithm (density-based spatial clustering of applications with noise) is given at: 
 * Link --> https://en.wikipedia.org/wiki/DBSCAN
 *
 */
public class DBSCAN2 {
	
	public static final int ClusterIDOfNoise = 0;
	
	protected ArrayList<DBScanPoint> mPoints; //data points
	protected float mEpsilon; //neighborhood parameter
	protected float mEpsilonSquared; // used for optimization
	protected int mMinPts; //neighborhood density parameter	
	protected int mClusterIndex;
	
	// ========================================================================
	public static int getClusterIDOfNoise() { return ClusterIDOfNoise; }
	
	// ========================================================================
	public DBSCAN2(ArrayList<Floint> points, float epsilon, int minPoints) {
		this.mPoints = addPoints( points, new ArrayList<DBScanPoint>(points.size()) );
		changeParams(epsilon, minPoints);
	}
	
	// ========================================================================
	/** Creates DBScanPoint object wrappers for the points. */
	protected static ArrayList<DBScanPoint> addPoints(ArrayList<Floint> points, ArrayList<DBScanPoint> targetList) {
		for (Floint floint : points) {
			DBScanPoint dbsp = new DBScanPoint(floint);
			targetList.add(dbsp);
		}
		return targetList;
	}
	
	// ========================================================================
	/** Changes the parameter values
	 *  @param newEpsilon   The new distance desired.  Enter a negative number for no change
	 *  @param newMinPoints The new number of points required for density.  Enter a negative number for no change.
	 */
	public void changeParams(float newEpsilon, int newMinPoints) {
		if (newEpsilon >= 0) {
			mEpsilon = newEpsilon;
			mEpsilonSquared = mEpsilon * mEpsilon; // don't use Math.pow(), as it's slow
		}
		
		if (newMinPoints >= 0) {
			mMinPts = newMinPoints;
		}
		
		// Now we need to reset the clustering meta-information per point
		// since the paramters changed.  cluster() will need to be run again
		for (DBScanPoint point : mPoints) {
			point.resetForClustering();
		}
		
		mClusterIndex = ClusterIDOfNoise;  // reset our cluster index
	}
	
	// ========================================================================
	/**
	 * Clusters points[] using DBSCAN
	 */
	public void cluster() {
		mClusterIndex = ClusterIDOfNoise; // start at default
		ArrayList<DBScanPoint> neighbors = new ArrayList<DBScanPoint>(3000000);
		
		for (DBScanPoint thePoint : mPoints) {				
			if (!thePoint.getVisited()) { 
				//System.out.println(i + " of " + points.length);
				thePoint.setVisited();
				neighbors = getNeighbors(thePoint, neighbors, true); //grab neighbors for points[i]
				//System.out.println(neighbors.size());
				if (neighbors.size() >= mMinPts) { //if neighborhood is dense enough
					//System.out.println("expanding cluster...");
					++mClusterIndex;
					expandCluster(thePoint, neighbors, mClusterIndex); //start expanding a cluster around points[i]					
				}
			}
		}
	}
	
	// ========================================================================
	protected void expandCluster(DBScanPoint point, ArrayList<DBScanPoint> neighbors, int clusterIndex) {
		point.mClusterAssigned = clusterIndex;   // assign our 'core' point to cluster c
		
		ArrayList<DBScanPoint> neighborsOfNeighbor = new ArrayList<DBScanPoint>(3000000);		
		
		// We perform a loop on the neighbors list, which can grow 
		// in size during the body of the loop.  Thus, neighbors list becomes
		// neighbors' neighbors, neighbors' neighbors' neighbors
		for (int i = 0; i < neighbors.size(); i++) {
			DBScanPoint theNeighbor = neighbors.get(i);

			if (i % 10000 == 0) {
				System.out.printf("\tNum Neighbors Processed (%d) out of %d for cluster: (%d)\n", i, mPoints.size(), clusterIndex);
			}
			
			if (!theNeighbor.getVisited()) {
				//count++;
				//System.out.println(count);
				theNeighbor.setVisited();
				
				neighborsOfNeighbor = getNeighbors(theNeighbor, neighborsOfNeighbor, true); //grab neighbor's neighbors
				if (neighborsOfNeighbor.size() >= mMinPts) { //if neighbor's neighborhood is dense enough
					//then let's iterate through them as well (they might be 'eligible' for inclusion in cluster c)
					//System.out.printf("Neighbors: %d, NeighborsToAdd: %d, Total: %d\n", neighbors.size(), neighborsOfNeighbor.size(), neighbors.size() + neighborsOfNeighbor.size());
					addAllElementsNotAlreadyAdded(neighbors, neighborsOfNeighbor);
				}
			}

			if (theNeighbor.mClusterAssigned == ClusterIDOfNoise) {
				theNeighbor.mClusterAssigned = clusterIndex;
			}
		}
	}

	// ========================================================================
	// We create this function because it is ironically more efficient than the ArrayList.addAll() method,
	// which stupidly allocates memory to create an extra and needless temporary array in its implementaion.
	// This adds a point if it wasn't already added.
	protected static void addAllElementsNotAlreadyAdded(Collection<DBScanPoint> listToWhichToAdd, Collection<DBScanPoint> elementsToAdd) {
		for (DBScanPoint element : elementsToAdd) {
			if (!element.mAdded) {
				listToWhichToAdd.add(element);
				element.mAdded = true;
			}
		}	
	}
	
	
	// ========================================================================
	// We create this function because it is ironically more efficient than the ArrayList.addAll() method,
	// which stupidly allocates memory to create an extra and needless temporary array in its implementaion.
	public static void addAll(Collection<DBScanPoint> listToWhichToAdd, Collection<DBScanPoint> elementsToAdd) {
		for (DBScanPoint element : elementsToAdd) {
			listToWhichToAdd.add(element);
		}
	}
	
	// ========================================================================
	// We must add elements to a list this way if the iterator for that list has already been invoked.
	// If we do not add elements to the list via the iterator, the iterator will throw an exception.
	public static void addAllToNextPosition(ListIterator<DBScanPoint> iter, Collection<DBScanPoint> elementsToAdd) {
		for (DBScanPoint element : elementsToAdd) {
			iter.add(element);
		}
	}
	
	// ========================================================================
	protected ArrayList<DBScanPoint> getNeighbors(DBScanPoint point, ArrayList<DBScanPoint> neighbors, boolean clearList) {
		neighbors = (neighbors == null) ? new ArrayList<DBScanPoint>(20000) : neighbors; 
		if (clearList) { 
			neighbors.clear(); 
		}
		
		for (DBScanPoint elem : mPoints) {
			if (point.mFloint.getCartesianDistance(elem.mFloint) < mEpsilon) { //if point is within parameter-defined neighborhood
				neighbors.add(elem);
			}				
		}
		return neighbors;
	}
		
	// ========================================================================
	public int[] getClustAssignments() {
		int[] thePoints = new int[mPoints.size()];
		for (int i = 0; i < thePoints.length; i++) {
			thePoints[i] = mPoints.get(i).mClusterAssigned; 
		}
		return thePoints;
	}
	
	// ========================================================================
	/**
	 * Returns the value in clust_assignments[] that occurs most often
	 */
	public int getLargestCluster() {
		if (mClusterIndex == ClusterIDOfNoise) return ClusterIDOfNoise;  // means we have nothing but noise
		
		int[] counts = new int[mClusterIndex + 1];
		for (DBScanPoint thePoint : mPoints) {
			counts[thePoint.mClusterAssigned]++;
		}
		
		int startIndex = ClusterIDOfNoise + 1;  // we start at 1 so we only consider points with valid cluster assignments
		int indexOfMaxCount = startIndex;
		for (int i = startIndex + 1; i < counts.length; i++) {
			if (counts[i] > counts[indexOfMaxCount]) {
				indexOfMaxCount = i;
			}
		}		

		//System.out.println("partition " + rtn + " :: " + counts.get(rtn) + " of " + points.length);
		return indexOfMaxCount;
	}
	
	// ========================================================================
	// INNER CLASS: DBScanPoint
	// ========================================================================
	/** This class contains meta information for each point in DBScan in order
	 *  to decrease run-time for the DBScan algorithm.
	 */
	protected static class DBScanPoint implements Comparable<DBScanPoint> {
		
		private boolean mVisited;
		public boolean mAdded;
		public int mClusterAssigned;
		
		public int mIndexSortedX;
		public int mIndexSortedY;		
		
		protected Floint mFloint; 
		
		protected DBScanPoint(Floint floint) {
			mFloint = floint;
			reset();
		}
		
		public boolean getVisited() { return mVisited; }
		public void    setVisited() { mVisited = true; }	
		public void  resetVisited() { mVisited = false; }
		
		public void reset() {
			resetForClustering();
			mIndexSortedX = mIndexSortedY = -1;
		}
		
		public void resetForClustering() {
			resetVisited();
			mClusterAssigned = DBSCAN2.ClusterIDOfNoise;
			mAdded = false;
		}
		
		public double getCartesianDistance(DBScanPoint rhs) {
			return mFloint.getCartesianDistance(rhs.mFloint);
		}
		
		public double getCartesianDistanceSquared(DBScanPoint rhs) {
			return mFloint.getCartesianDistanceSquared(rhs.mFloint);
		}
		
		/** Default sorting by x-coordinate */
		public int compareTo(DBScanPoint rhs) {
			return mFloint.compareTo(rhs.mFloint);
		}

		// ========================================================================
		// Global Comparator object for public use	
		public static final DBScanPointCompareX  DBScanPointCompareXObj  = new DBScanPointCompareX();
		public static final DBScanPointCompareY  DBScanPointCompareYObj  = new DBScanPointCompareY();
		
		public static final DBScanPointDistanceX DBScanPointDistanceXObj = new DBScanPointDistanceX();
		public static final DBScanPointDistanceY DBScanPointDistanceYObj = new DBScanPointDistanceY();
		
		// ====================================================================
		/** Comparator. */
		public static class DBScanPointCompareX implements Comparator<DBScanPoint> {
			public int compare(DBScanPoint p1, DBScanPoint p2) {
				return Floint.FlointCompareXObj.compare(p1.mFloint, p2.mFloint);
			}
		}
		
		public static class DBScanPointCompareY implements Comparator<DBScanPoint> {
			public int compare(DBScanPoint p1, DBScanPoint p2) {
				return Floint.FlointCompareYObj.compare(p1.mFloint, p2.mFloint);
			}
		}
		
		// ====================================================================
		public static abstract class DBScanPointDistance {
			public abstract double distance(DBScanPoint n1, DBScanPoint n2); 
		}
		
		public static class DBScanPointDistanceX extends DBScanPointDistance {
			public double distance(DBScanPoint n1, DBScanPoint n2) { 
				return Floint.FlointDistanceXObj.distance(n1.mFloint, n2.mFloint); 
			} 
		}
		
		public static class DBScanPointDistanceY extends DBScanPointDistance {
			public double distance(DBScanPoint n1, DBScanPoint n2) { 
				return Floint.FlointDistanceYObj.distance(n1.mFloint, n2.mFloint); 
			}
		}		
	}
}
