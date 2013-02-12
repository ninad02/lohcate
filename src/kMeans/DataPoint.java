package kMeans;

/**
 * From: http://www.codecodex.com/wiki/K-means_cluster_analysis_algorithm
This class represents a candidate for Cluster analysis. A candidate must have
a name and two independent variables on the basis of which it is to be clustered.
A Data Point must have two variables and a name. A List of  Data Point object
is fed into the constructor of the JCA class. JCA and DataPoint are the only
classes which may be available from other packages.
@author Shyam Sivaraman
@version 1.0
@see JCA
@see Cluster
 */

public class DataPoint {
	private double mX,mY;
	private String mObjName;
	private Cluster mCluster;
	private double mEuDt;

	public DataPoint(double x, double y, String name) {
		mX = x;
		mY = y;
		mObjName = name;
	}

	public void setCluster(Cluster cluster) {
		mCluster = cluster;
		calcEuclideanDistance();
	}

	public void calcEuclideanDistance() { 

		//called when DP is added to a cluster or when a Centroid is recalculated.
		mEuDt = Math.hypot(mX - mCluster.getCentroid().getCx(),
				mY - mCluster.getCentroid().getCy());
	}

	public double testEuclideanDistance(Centroid c) {
		return Math.sqrt(Math.pow((mX - c.getCx()), 2) + Math.pow((mY - c.getCy()), 2));
	}

	public double getX() {
		return mX;
	}

	public double getY() {
		return mY;
	}

	public Cluster getCluster() {
		return mCluster;
	}

	public double getCurrentEuDt() {
		return mEuDt;
	}

	public String getObjName() {
		return mObjName;
	}

}

