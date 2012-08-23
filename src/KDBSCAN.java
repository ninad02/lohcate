import java.util.ArrayList;
/**
 * LOHcate --- A software tool for LOH calling and visualization in cancer genomes
 * D Wheeler & SG Reddy
 * Human Genome Sequencing Center, Baylor College of Medicine (Houston, TX)
 * 
 * @author Siddharth G. Reddy
 *
 */

import java.util.ArrayList;
/**
 * LOHcate --- A software tool for LOH calling and visualization in cancer genomes
 * D Wheeler & SG Reddy
 * Human Genome Sequencing Center, Baylor College of Medicine (Houston, TX)
 * 
 * I tried implementing a variation on DBSCAN, but it was too slow and ended up requiring the same # of parameters
 * Link --> http://www.scribd.com/doc/21053830/Improved-DBSCAN-for-spatial-databases-with-noise-and-different-densities
 * 
 * @author Siddharth G. Reddy
 *
 */
public class KDBSCAN {
	
	private int[] visited;
	private Floint[] points;
	private float dr;
	private int k;
	private int[] clust_assignments;
	private Neighborhood[] big_neighbs;
	
	public KDBSCAN(Floint[] points, int k, float dr) {
		this.points = points;
		this.dr = dr;
		this.k = k;
		clust_assignments = new int[points.length];
		visited = new int[points.length];
		big_neighbs = new Neighborhood[points.length];
		for (int i = 0; i<points.length; i++) {
			points[i].loc = i;
			clust_assignments[i] = 0;
			visited[i] = 0;
		}
	}
	
	private void calcNeighbs() {
		for (int i = 0; i<big_neighbs.length; i++) {
			//System.out.println(i + " of " + big_neighbs.length);
			big_neighbs[i] = getNeighbors(i);
		}
	}
	
	private Neighborhood getNeighbors(int param) {
		ArrayList<Floint> neighbors = new ArrayList<Floint>();//new Floint[k];
		for (int i = 0; i<k; i++)
			neighbors.add(new Floint(0, 0, 0));
		float[] dists = new float[k];
		float temp, prev = 0;
		for (int i = 0; i<neighbors.size(); i++) {
			dists[i] = Float.MAX_VALUE;
			for (int j = 0; j<points.length; j++) {
				if (j!=param) {
					temp = (float) Math.sqrt(Math.pow(points[j].x - points[param].y, 2) + Math.pow(points[j].y - points[param].y, 2));
					if (prev < temp && temp < dists[i]) {
						dists[i] = temp;
						neighbors.set(i, points[j]);
					}
				}
			}
			prev = dists[i];
		}
		return new Neighborhood(neighbors, dists);
	}
	
	public void cluster() {
		calcNeighbs();
		int c = 1;
		ArrayList<Floint> neighbors;
		for (int i = 0; i<points.length; i++) {
			if (visited[i]==0) {
				//System.out.println(i + " of " + points.length);
				visited[i] = 1;
				if (big_neighbs[i].density >= dr) {
					//System.out.println("expanding cluster...");
					expandCluster(points[i], big_neighbs[i].members, c);
					c++;
				}
			}
		}
	}
	
	private void expandCluster(Floint point, ArrayList<Floint> neighbors, int c) {
		clust_assignments[point.loc] = c;
		ArrayList<Floint> nneighbors;
		//int count = 0;
		for (int i = 0; i<neighbors.size(); i++) {
			if (visited[neighbors.get(i).loc]==0) {
				//System.out.println(count);
				visited[neighbors.get(i).loc] = 1;
				//count++;
				if (big_neighbs[neighbors.get(i).loc].density >= dr)
					neighbors = join(neighbors, big_neighbs[neighbors.get(i).loc].members);
			}
			if (clust_assignments[neighbors.get(i).loc]==0)
				clust_assignments[neighbors.get(i).loc] = c;
		}
	}
	
	private ArrayList<Floint> join(ArrayList<Floint> one, ArrayList<Floint> two) {
		for (Floint point : two)
			one.add(point);
		return one;
	}
	
	private ArrayList<Floint> getNeighbors(Floint point) {
		return big_neighbs[point.loc].members;
	}
	
	public int[] getClustAssignments() { return clust_assignments; }
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
		System.out.println("partition " + rtn + " :: " + counts.get(rtn) + " of " + points.length);
		return rtn;
	}
}

class Neighborhood {
	public ArrayList<Floint> members;
	public float[] dists;
	public float density;
	public Neighborhood(ArrayList<Floint> members, float[] dists) {
		this.members = members;
		this.dists = dists;
		calcDensity();
		//System.out.println(density);
	}
	public void calcDensity() {
		density = 0f;
		for (int i = 0; i<dists.length; i++)
			density += dists[i];
	}
}