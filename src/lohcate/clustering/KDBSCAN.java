package lohcate.clustering;
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
	/*
	private Floint[] points;
	private float dr;
	private int k;
	private Neighborhood[] big_neighbs;
	
	public KDBSCAN(Floint[] points, int k, float dr) {
		this.points = points;
		this.dr = dr;
		this.k = k;
		
		big_neighbs = new Neighborhood[points.length];
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
			neighbors.add(new Floint(0, 0));
		float[] dists = new float[k];
		float temp, prev = 0;
		for (int i = 0; i<neighbors.size(); i++) {
			dists[i] = Float.MAX_VALUE;
			for (int j = 0; j<points.length; j++) {
				if (j!=param) {
					temp = (float) Math.sqrt(Math.pow(points[j].mX - points[param].mY, 2) + Math.pow(points[j].mY - points[param].mY, 2));
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
		clust_assignments[point.mIndex] = c;
		ArrayList<Floint> nneighbors;
		//int count = 0;
		for (int i = 0; i<neighbors.size(); i++) {
			if (visited[neighbors.get(i).mIndex]==0) {
				//System.out.println(count);
				visited[neighbors.get(i).mIndex] = 1;
				//count++;
				if (big_neighbs[neighbors.get(i).mIndex].density >= dr)
					neighbors.addAll(big_neighbs[neighbors.get(i).mIndex].members);					
			}
			if (clust_assignments[neighbors.get(i).mIndex]==0)
				clust_assignments[neighbors.get(i).mIndex] = c;
		}
	}

	
	private ArrayList<Floint> getNeighbors(Floint point) {
		return big_neighbs[point.mIndex].members;
	}
	
	public int[] getClustAssignments() { return clust_assignments; }
	
	
	public int getLargestCluster() {
		if (mClusterIndex == 0) return -1;
		
		int[] counts = new int[mClusterIndex + 1];
		for (Floint thePoint : points) {
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

		System.out.println("partition " + indexOfMaxCount + " :: " + counts[indexOfMaxCount] + " of " + points.length);
		return indexOfMaxCount;
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
	}*/
}