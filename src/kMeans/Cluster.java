package kMeans;

import java.util.ArrayList;

/**
 * From: http://www.codecodex.com/wiki/K-means_cluster_analysis_algorithm
 * This class represents a Cluster in a Cluster Analysis Instance. A Cluster is associated
 * with one and only one JCA Instance. A Cluster is related to more than one DataPoints and
 * one centroid.
 * @author Shyam Sivaraman
 * @version 1.1
 * @see DataPoint
 * @see Centroid
 */



public class Cluster implements Comparable<Cluster> {
    private String mName;
    private Centroid mCentroid; //will be set by calling setCentroid()
    private double mSumSqr;
    private ArrayList<DataPoint> mDataPoints = new ArrayList<DataPoint>();

    public Cluster(String name) {
        mName = name;
    }

    public void setCentroid(Centroid c) {
        mCentroid = c;
    }

    public Centroid getCentroid() {
        return mCentroid;
    }

    public void addDataPoint(DataPoint dp) { //called from CAInstance
        dp.setCluster(this); //initiates a inner call to calcEuclideanDistance() in DP.
        mDataPoints.add(dp);
        calcSumOfSquares();
    }

    public void removeDataPoint(DataPoint dp) {
        mDataPoints.remove(dp);
        calcSumOfSquares();
    }

    public int getNumDataPoints() {
        return mDataPoints.size();
    }

    public DataPoint getDataPoint(int pos) {
        return mDataPoints.get(pos);
    }

    public void calcSumOfSquares() { //called from Centroid
        double temp = 0;
        for (DataPoint dp : mDataPoints) {
            temp += dp.getCurrentEuDt();
        }
        mSumSqr = temp;
    }

    public double getSumSqr() {
        return mSumSqr;
    }

    public String getName() {
        return mName;
    }

    public ArrayList<DataPoint> getDataPoints() {
        return mDataPoints;
    }
    
    public int compareTo(Cluster rhs) {
    	return mCentroid.compareTo(rhs.mCentroid);
    }

}
