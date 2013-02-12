package lohcate.clustering;
import java.util.ArrayList;

/**
 * LOHcate --- A software tool for LOH calling and visualization in cancer genomes
 * D Wheeler & SG Reddy
 * Human Genome Sequencing Center, Baylor College of Medicine (Houston, TX)
 * 
 * An implementation of the DBSCAN (density-based spatial clustering of applications with noise) algorithm.
 * Link --> https://en.wikipedia.org/wiki/DBSCAN
 * 
 * @author Siddharth G. Reddy
 *
 */
public class DBSCAN {
	
	protected ArrayList<DBSCAN2.DBScanPoint> mPoints; //data points
	protected float mEps; //neighborhood parameter
	protected int mMinPts; //neighborhood density parameter
	
	protected int mClusterIndex = 0;
	
	public DBSCAN(ArrayList<Floint> points, float eps, int minPts) {
		//this.mPoints = points;
		this.mEps = eps;
		this.mMinPts = minPts;
	}
	
	/**
	 * Clusters points[] using DBSCAN
	 */
	public void cluster() {
		mClusterIndex = 0; // start at default
		ArrayList<DBSCAN2.DBScanPoint> neighbors;
		
		for (DBSCAN2.DBScanPoint thePoint : mPoints) {				
			if (!thePoint.getVisited()) { 
				//System.out.println(i + " of " + points.length);
				thePoint.setVisited();
				neighbors = getNeighbors(thePoint); //grab neighbors for points[i]
				//System.out.println(neighbors.size());
				if (neighbors.size() >= mMinPts) { //if neighborhood is dense enough
					//System.out.println("expanding cluster...");
					++mClusterIndex;
					expandCluster(thePoint, neighbors, mClusterIndex); //start expanding a cluster around points[i]					
				}
			}
		}
	}
	
	private void expandCluster(DBSCAN2.DBScanPoint point, ArrayList<DBSCAN2.DBScanPoint> neighbors, int c) {
		point.mClusterAssigned = c;   // assign our 'core' point to cluster c
		
		ArrayList<DBSCAN2.DBScanPoint> neighborsOfNeighbors;		
		//int count = 0;
		
		// We perform a loop on the neighbors list, which can grow 
		// in size during the body of the loop. 
		for (int i = 0; i < neighbors.size(); i++) { //iterate through neighbors (later becomes neighbors' neighbors', neighbors' neighbors' neighbors', &c.) of core point
			DBSCAN2.DBScanPoint theNeighbor = neighbors.get(i);
			
			if (!theNeighbor.getVisited()) {
				//count++;
				//System.out.println(count);
				theNeighbor.setVisited();								
				
				neighborsOfNeighbors = getNeighbors(theNeighbor); //grab neighbor's neighbors			
				if (neighborsOfNeighbors.size() >= mMinPts) { //if neighbor's neighborhood is dense enough
					//then let's iterate through them as well (they might be 'eligible' for inclusion in cluster c)
					neighbors.addAll(neighborsOfNeighbors);
				}
			}

			// Assign neighbor to current cluster if it hasn't already been assigned to a cluster
			if (theNeighbor.mClusterAssigned == 0) {
				theNeighbor.mClusterAssigned = c;
			}
		}
	}
	
	
	private ArrayList<DBSCAN2.DBScanPoint> getNeighbors(DBSCAN2.DBScanPoint point) {
		ArrayList<DBSCAN2.DBScanPoint> pre_rtn = new ArrayList<DBSCAN2.DBScanPoint>();
		
		for (DBSCAN2.DBScanPoint elem : mPoints)
			if (Math.sqrt(Math.pow(elem.mFloint.mX - point.mFloint.mX, 2) + Math.pow(elem.mFloint.mY - point.mFloint.mY, 2)) < mEps) //if point is within parameter-defined neighborhood
				pre_rtn.add(elem);
		return pre_rtn;
	}
	
	public int[] getClustAssignments() {
		int[] thePoints = new int[mPoints.size()];
		for (int i = 0; i < thePoints.length; i++) {
			thePoints[i] = mPoints.get(i).mClusterAssigned; 
		}
		return thePoints;
	}
	
	/**
	 * Returns the value in clust_assignments[] that occurs most often
	 */
	public int getLargestCluster() {
		if (mClusterIndex == 0) return -1;
		
		int[] counts = new int[mClusterIndex + 1];
		for (DBSCAN2.DBScanPoint thePoint : mPoints) {
			counts[thePoint.mClusterAssigned]++;
		}
		
		int startIndex = 1;  // we start at 1 so we only consider points with valid cluster assignments
		int indexOfMaxCount = startIndex;
		for (int i = startIndex + 1; i < counts.length; i++) {
			if (counts[i] > counts[indexOfMaxCount]) {
				indexOfMaxCount = i;
			}
		}

		//System.out.println("partition " + rtn + " :: " + counts.get(rtn) + " of " + points.length);
		return indexOfMaxCount;
	}
}
