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
	
	private int[] visited; //0 --> unvisited, 1 --> visited
	private Floint[] points; //data points
	private float eps; //neighborhood parameter
	private int minPts; //neighborhood density parameter
	private int[] clust_assignments; //cluster assignments
	
	public DBSCAN(Floint[] points, float eps, int minPts) {
		this.points = points;
		this.eps = eps;
		this.minPts = minPts;
		clust_assignments = new int[points.length];
		visited = new int[points.length];
		for (int i = 0; i<points.length; i++) {
			points[i].loc = i;
			clust_assignments[i] = 0;
			visited[i] = 0;
		}
	}
	
	/**
	 * Clusters points[] using DBSCAN
	 */
	public void cluster() {
		int c = 1; //need to start at an integer other than default (== 0)
		ArrayList<Floint> neighbors;
		for (int i = 0; i<points.length; i++) { //iterate through data points
			if (visited[i]==0) { //if points[i] has not been visited yet
				//System.out.println(i + " of " + points.length);
				visited[i] = 1; //we don't need to visit points[i] again
				neighbors = getNeighbors(points[i]); //grab neighbors for points[i]
				//System.out.println(neighbors.size());
				if (neighbors.size() >= minPts) { //if neighborhood is dense enough
					//System.out.println("expanding cluster...");
					expandCluster(points[i], neighbors, c); //start expanding a cluster around points[i]
					c++;
				}
			}
		}
	}
	
	private void expandCluster(Floint point, ArrayList<Floint> neighbors, int c) {
		clust_assignments[point.loc] = c; //assign our 'core' point to cluster c
		ArrayList<Floint> nneighbors;
		//int count = 0;
		for (int i = 0; i<neighbors.size(); i++) { //iterate through neighbors (later becomes neighbors' neighbors', neighbors' neighbors' neighbors', &c.) of core point
			if (visited[neighbors.get(i).loc]==0) { //if extended neighbor has not been visited yet
				//System.out.println(count);
				visited[neighbors.get(i).loc] = 1; //mark extended neighbor as visited
				//count++;
				nneighbors = getNeighbors(neighbors.get(i)); //grab neighbor's neighbors
				if (nneighbors.size() >= minPts) //if neighbor's neighborhood is dense enough
					neighbors = join(neighbors, nneighbors); //then let's iterate through them as well (they might be 'eligible' for inclusion in cluster c)
			}
			if (clust_assignments[neighbors.get(i).loc]==0) //if the neighbor has not already been assigned to a cluster
				clust_assignments[neighbors.get(i).loc] = c; //neighbor has become part of the current cluster
		}
	}
	
	private ArrayList<Floint> join(ArrayList<Floint> one, ArrayList<Floint> two) {
		for (Floint point : two)
			one.add(point);
		return one;
	}
	
	private ArrayList<Floint> getNeighbors(Floint point) {
		ArrayList<Floint> pre_rtn = new ArrayList<Floint>();
		for (Floint elem : points)
			if (Math.sqrt(Math.pow(elem.x - point.x, 2) + Math.pow(elem.y - point.y, 2)) < eps) //if point is within parameter-defined neighborhood
				pre_rtn.add(elem);
		return pre_rtn;
	}
	
	public int[] getClustAssignments() { return clust_assignments; }
	
	/**
	 * Returns the value in clust_assignments[] that occurs most often
	 */
	public int getLargestCluster() {
		int rtn = 0, ind;
		ArrayList<Integer> counts = new ArrayList<Integer>();
		int max = 0;
		for (int i = 0; i<points.length; i++)
			if (clust_assignments[i] > max)
				max = clust_assignments[i];
		for (int i = 0; i<=max; i++)
			counts.add(0);
		
		for (int i = 0; i<points.length; i++)
			counts.set(clust_assignments[i], counts.get(clust_assignments[i]) + 1);
		
		for (int i = 0; i<counts.size(); i++)
			if (counts.get(i) > counts.get(rtn))
				rtn = i;
		//System.out.println("partition " + rtn + " :: " + counts.get(rtn) + " of " + points.length);
		return rtn;
	}
}
