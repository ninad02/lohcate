package kMeans;

/**
 * From: http://www.codecodex.com/wiki/K-means_cluster_analysis_algorithm
 * This class represents the Centroid for a Cluster. The initial centroid is calculated
 * using a equation which divides the sample space for each dimension into equal parts
 * depending upon the value of k.
 * @author Shyam Sivaraman
 * @version 1.0
 * @see Cluster
 */

public class Centroid implements Comparable<Centroid> {
    private double mCx, mCy;
    private Cluster mCluster;

    public Centroid(double cx, double cy) {
        mCx = cx;
        mCy = cy;
    }

    public void calcCentroid() { //only called by CAInstance
        int numDP = mCluster.getNumDataPoints();
        double tempX = 0, tempY = 0;
       // int i;
        //caluclating the new Centroid
        for (DataPoint dp : mCluster.getDataPoints()) {
            tempX += dp.getX(); 
            //total for x
            tempY += dp.getY(); 
            //total for y
        }
        mCx = tempX / numDP;
        mCy = tempY / numDP;
        //calculating the new Euclidean Distance for each Data Point
        for (DataPoint dp : mCluster.getDataPoints()) {
            dp.calcEuclideanDistance();
        }
        //calculate the new Sum of Squares for the Cluster
        mCluster.calcSumOfSquares();
    }

    public void setCluster(Cluster c) {
        mCluster = c;
    }

    public double getCx() {
        return mCx;
    }

    public double getCy() {
        return mCy;
    }

    public Cluster getCluster() {
        return mCluster;
    }

    public int compareTo(Centroid rhs) {
    	if (mCx > rhs.mCx) {
    		return 1;
    	} else if (mCx < rhs.mCx) {
    		return -1;
    	} else {
    		if (mCy > rhs.mCy) {
    			return 1;
    		} else if (mCy < rhs.mCy) {
    			return -1;
    		} else {
    			return 0;
    		}
    	}
    }
}

