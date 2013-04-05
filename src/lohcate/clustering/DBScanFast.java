package lohcate.clustering;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.Date;


public class DBScanFast extends DBSCAN2 {
	
	private ArrayList<DBScanPoint> mPointsSortedX;
	private ArrayList<DBScanPoint> mPointsSortedY;	
	
	public DBScanFast(ArrayList<Floint> points, float eps, int minPts) {
		super(points, eps, minPts);
		mPointsSortedX = new ArrayList<DBScanPoint>(points.size());
		mPointsSortedY = new ArrayList<DBScanPoint>(points.size());
		duplicateAndSortPoints();
	}
	
	protected void duplicateAndSortPoints() {
		mPointsSortedX.clear();
		mPointsSortedY.clear();
		
		mPointsSortedX.addAll(mPoints);
		mPointsSortedY.addAll(mPoints);
		
		Collections.sort(mPointsSortedX, DBScanPoint.DBScanPointCompareXY);
		Collections.sort(mPointsSortedY, DBScanPoint.DBScanPointCompareYX);		
		
		System.out.println("\tAssignStart: " + (new Date()).toString());
		assignIndexToPoints(true,  mPointsSortedX, DBScanPoint.DBScanPointCompareXY);
		assignIndexToPoints(false, mPointsSortedY, DBScanPoint.DBScanPointCompareYX);
		System.out.println("\tAssignEnd: " + (new Date()).toString());
	}
	
	private void assignIndexToPoints(boolean doX, ArrayList<DBScanPoint> pointsSorted, Comparator<DBScanPoint> theComparator) {
		boolean[] assigned = new boolean[mPoints.size()];
		Arrays.fill(assigned, false);
		
		for (DBScanPoint point : mPoints) {
			int resultIndex = Collections.binarySearch(pointsSorted, point, theComparator);
			
			if (resultIndex >= 0) {
				int assignedIndex = resultIndex;
				if (assigned[resultIndex]) {
					assignedIndex = assignIndexHelper(point, assigned, pointsSorted, theComparator, resultIndex, true);
					if (assignedIndex < 0) {
						assignedIndex = assignIndexHelper(point, assigned, pointsSorted, theComparator, resultIndex, false);
						if (assignedIndex < 0) {
							System.err.println("The point cannot be found!");
							System.exit(-1);
						}
					}					
				} 
				
				assigned[assignedIndex] = true;
				if (doX) {
					point.mIndexSortedX = assignedIndex;
				} else {
					point.mIndexSortedY = assignedIndex;
				}
			} else {
				System.err.println("BAAD!");
				System.exit(-1);
			}
		}
	}
	
	private int assignIndexHelper(DBScanPoint point, boolean[] assigned, ArrayList<DBScanPoint> pointsSorted, Comparator<DBScanPoint> theComparator, int centralIndex, boolean forward) {
		int increment = forward ? 1 : -1;		
		
		for (int i = centralIndex + increment; 
				(i >= 0 && i < pointsSorted.size()) 
				&& (theComparator.compare(point, pointsSorted.get(i)) == 0); 				
				i += increment) {
			
			if (!assigned[i]) {
				 return i;				
			}
		}
		return -1;
	}

	/** @Override 
	 * 	This procedure will try to optimize determining the neighbors, examining only
	 *  those points within a bounded search space, in which the bounds of the space
	 *  represent the maximum theoretical distances as given 
	 */
	protected ArrayList<DBScanPoint> getNeighbors(DBScanPoint point, ArrayList<DBScanPoint> neighbors, boolean clearList) {
		neighbors = (neighbors == null) ? new ArrayList<DBScanPoint>(20000) : neighbors; 
		if (clearList) { 
			neighbors.clear(); 
		}
		
		getNeighborsHelper(point, neighbors, mPointsSortedX, DBScanPoint.DBScanPointCompareXY, DBScanPoint.DBScanPointDistanceX, true);
		getNeighborsHelper(point, neighbors, mPointsSortedY, DBScanPoint.DBScanPointCompareYX, DBScanPoint.DBScanPointDistanceY, false);
		return neighbors;
	}
	
	private void getNeighborsHelper(DBScanPoint point, ArrayList<DBScanPoint> neighbors, ArrayList<DBScanPoint> sortedPoints, Comparator<DBScanPoint> theComparator, DBScanPoint.DBScanPointDistance distanceMethodOneAxis, boolean doX) {
		getNeighborsHelperHelper(point, neighbors, sortedPoints, theComparator, distanceMethodOneAxis, true,  doX);
		getNeighborsHelperHelper(point, neighbors, sortedPoints, theComparator, distanceMethodOneAxis, false, doX);				
	}

	private void getNeighborsHelperHelper(DBScanPoint point, ArrayList<DBScanPoint> neighbors, ArrayList<DBScanPoint> sortedPoints, Comparator<DBScanPoint> theComparator, DBScanPoint.DBScanPointDistance distanceMethodOneAxis, boolean moveRight, boolean doX) {
		double theoreticalMaxDistance = point.mFloint.getTheoreticalFurthestDifferenceXOrYWithinDistance(mEpsilon);
		int increment = moveRight ? 1 : -1;
		DBScanPoint candidateNeighbor = null;
		int centralIndex = doX ? point.mIndexSortedX : point.mIndexSortedY;
				
		for (int i = centralIndex + increment; 
				((i >= 0) && (i < sortedPoints.size())) 
				&& ((candidateNeighbor = sortedPoints.get(i)) != null) 
				&& (Math.abs(distanceMethodOneAxis.distance(point, candidateNeighbor)) <= theoreticalMaxDistance); 
				i += increment) {
			
			if (point.getCartesianDistanceSquared(candidateNeighbor) <= mEpsilonSquared) {
				if (doX) {
					neighbors.add(candidateNeighbor);
					candidateNeighbor.mAdded = true;
				} else {
					if (!candidateNeighbor.mAdded) {
						neighbors.add(candidateNeighbor);
					}
					candidateNeighbor.mAdded = false;  // reset the status since no need for this field anymore
				}
			}				
		}
	}

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub

	}

}
