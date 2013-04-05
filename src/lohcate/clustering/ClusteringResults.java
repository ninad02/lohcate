package lohcate.clustering;

import java.util.ArrayList;

import com.carrotsearch.hppc.IntArrayList;

public class ClusteringResults<E> {
	
	public static final int InvalidSubClusteringID = Integer.MIN_VALUE;

	// ========================================================================
	public ArrayList<E> mClassificationForPoint;
	public IntArrayList mSubClusterID;
		
	// ========================================================================
	public ClusteringResults(int numExpectedPoints) {
		mClassificationForPoint = new ArrayList<E>(numExpectedPoints);
		mSubClusterID = new IntArrayList(numExpectedPoints);
	}
	
	// ========================================================================
	public ClusteringResults() {
		this(10);
	}

	// ========================================================================
	public void initializeResults(int numPoints) {
		initializeResults(numPoints, null);
	}
	
	// ========================================================================
	public void initializeResults(int numPoints, E initialValue) {
		clear();
		for (int i = 0; i < numPoints; i++) {
			mClassificationForPoint.add(initialValue);
			mSubClusterID.add(InvalidSubClusteringID);
		}
	}
	
	// ========================================================================
	public void clear() {
		mClassificationForPoint.clear();
		mSubClusterID.clear();
	}
	
	// ========================================================================
	public E getClassification(int index) { return mClassificationForPoint.get(index); }
	
	// ========================================================================
	public int getSubClusterID(int index) { return mSubClusterID.get(index); }
	
	// ========================================================================
	public void setClassification(int index, E cluster, int subClusterID) {
		mClassificationForPoint.set(index, cluster);
		mSubClusterID.set(index, subClusterID);
	}
	
	// ========================================================================
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub

	}

}
