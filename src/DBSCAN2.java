import java.util.ArrayList;
import java.util.Collection;
import java.util.LinkedList;
import java.util.ListIterator;
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
public class DBSCAN2 {
	
	protected ArrayList<Floint> mPoints; //data points
	protected float mEps; //neighborhood parameter
	protected int mMinPts; //neighborhood density parameter
	
	protected int mClusterIndex = 0;
	
	public DBSCAN2(ArrayList<Floint> points, float eps, int minPts) {
		this.mPoints = points;
		this.mEps = eps;
		this.mMinPts = minPts;
	}
	
	/**
	 * Clusters points[] using DBSCAN
	 */
	public void cluster() {
		mClusterIndex = 0; // start at default
		ArrayList<Floint> neighbors = new ArrayList<Floint>(30000);
		
		for (Floint thePoint : mPoints) {				
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
	
	protected void expandCluster(Floint point, ArrayList<Floint> neighbors, int c) {
		point.mClusterAssigned = c;   // assign our 'core' point to cluster c
		
		ArrayList<Floint> neighborsOfNeighbor = new ArrayList<Floint>(30000);		
		
		// We perform a loop on the neighbors list, which can grow 
		// in size during the body of the loop.  Thus, neighbors list becomes
		// neighbors' neighbors, neighbors' neighbors' neighbors
		//for (ListIterator<Floint> iter = neighbors.listIterator(); iter.hasNext(); ) {
			//Floint theNeighbor = iter.next();
		for (int i = 0; i < neighbors.size(); i++) {
			Floint theNeighbor = neighbors.get(i);

			if (!theNeighbor.getVisited()) {
				//count++;
				//System.out.println(count);
				theNeighbor.setVisited();
				
				neighborsOfNeighbor = getNeighbors(theNeighbor, neighborsOfNeighbor, true); //grab neighbor's neighbors
				if (neighborsOfNeighbor.size() >= mMinPts) { //if neighbor's neighborhood is dense enough
					//then let's iterate through them as well (they might be 'eligible' for inclusion in cluster c)
					//addAllToNextPosition(iter, neighborsOfNeighbor);
					addAll(neighbors, neighborsOfNeighbor);
				}
			}

			if (theNeighbor.mClusterAssigned == 0) {
				theNeighbor.mClusterAssigned = c;
			}
		}
	}
	
	// We create this function because it is ironically more efficient than the ArrayList.addAll() method,
	// which stupidly allocates memory to create an extra and needless temporary array in its implementaion.
	public static void addAll(Collection<Floint> listToWhichToAdd, Collection<Floint> elementsToAdd) {
		for (Floint element : elementsToAdd) {
			listToWhichToAdd.add(element);
		}
	}
	
	// We must add elements to a list this way if the iterator for that list has already been invoked.
	// If we do not add elements to the list via the iterator, the iterator will throw an exception.
	public static void addAllToNextPosition(ListIterator<Floint> iter, Collection<Floint> elementsToAdd) {
		for (Floint element : elementsToAdd) {
			iter.add(element);
		}
	}
	
	protected ArrayList<Floint> getNeighbors(Floint point, ArrayList<Floint> neighbors, boolean clearList) {
		neighbors = (neighbors == null) ? new ArrayList<Floint>(10000) : neighbors; 
		if (clearList) { 
			neighbors.clear(); 
		}
		
		for (Floint elem : mPoints) {
			if (point.getCartesianDistance(elem) < mEps) { //if point is within parameter-defined neighborhood
				neighbors.add(elem);
			}				
		}
		return neighbors;
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
		for (Floint thePoint : mPoints) {
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
