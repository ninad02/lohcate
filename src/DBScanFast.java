import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.Date;


public class DBScanFast extends DBSCAN2 {
	
	private ArrayList<Floint> mPointsSortedX;
	private ArrayList<Floint> mPointsSortedY;	
	
	public DBScanFast(ArrayList<Floint> points, float eps, int minPts) {
		super(points, eps, minPts);
		mPointsSortedX = new ArrayList<Floint>();
		mPointsSortedY = new ArrayList<Floint>();
		duplicateAndSortPoints();
	}
	
	protected void duplicateAndSortPoints() {
		mPointsSortedX.clear();
		mPointsSortedY.clear();
		
		mPointsSortedX.addAll(mPoints);
		mPointsSortedY.addAll(mPoints);
		
		Collections.sort(mPointsSortedX, Floint.FlointCompareXObj);
		Collections.sort(mPointsSortedY, Floint.FlointCompareYObj);		
		
		System.out.println("\tAssignStart: " + (new Date()).toString());
		assignIndexToPoints(true,  mPointsSortedX, Floint.FlointCompareXObj);
		assignIndexToPoints(false, mPointsSortedY, Floint.FlointCompareYObj);
		System.out.println("\tAssignEnd: " + (new Date()).toString());
	}
	
	private void assignIndexToPoints(boolean doX, ArrayList<Floint> pointsSorted, Comparator<Floint> theComparator) {
		boolean[] assigned = new boolean[mPoints.size()];
		Arrays.fill(assigned, false);
		
		for (Floint point : mPoints) {
			int resultIndex = Collections.binarySearch(pointsSorted, point, theComparator);
			
			if (resultIndex >= 0) {
				int assignedIndex = resultIndex;
				if (assigned[resultIndex]) {
					assignedIndex = assignIndexHelper(point, assigned, pointsSorted, resultIndex, true);
					if (assignedIndex < 0) {
						assignedIndex = assignIndexHelper(point, assigned, pointsSorted, resultIndex, false);
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
	
	private int assignIndexHelper(Floint point, boolean[] assigned, ArrayList<Floint> pointsSorted, int centralIndex, boolean forward) {
		int increment = forward ? 1 : -1;		
		
		for (int i = centralIndex + increment; 
				(i >= 0 && i < pointsSorted.size()) 
				&& point.compareTo(pointsSorted.get(i)) == 0; 
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
	protected ArrayList<Floint> getNeighbors(Floint point, ArrayList<Floint> neighbors, boolean clearList) {
		neighbors = (neighbors == null) ? new ArrayList<Floint>(10000) : neighbors; 
		if (clearList) { 
			neighbors.clear(); 
		}
		
		getNeighborsHelper(point, neighbors, mPointsSortedX, Floint.FlointCompareXObj, Floint.FlointDistanceXObj, true);
		getNeighborsHelper(point, neighbors, mPointsSortedY, Floint.FlointCompareYObj, Floint.FlointDistanceYObj, false);
		return neighbors;
	}
	
	private void getNeighborsHelper(Floint point, ArrayList<Floint> neighbors, ArrayList<Floint> sortedPoints, Comparator<Floint> theComparator, Floint.FlointDistance distanceMethodOneAxis, boolean doX) {
		getNeighborsHelperHelper(point, neighbors, sortedPoints, theComparator, distanceMethodOneAxis, true,  doX);
		getNeighborsHelperHelper(point, neighbors, sortedPoints, theComparator, distanceMethodOneAxis, false, doX);				
	}

	private void getNeighborsHelperHelper(Floint point, ArrayList<Floint> neighbors, ArrayList<Floint> sortedPoints, Comparator<Floint> theComparator, Floint.FlointDistance distanceMethodOneAxis, boolean goForward, boolean doX) {
		double theoreticalMaxDistance = point.getTheoreticalFurthestDifferenceXOrYWithinDistance(mEps);
		int increment = goForward ? 1 : -1;
		Floint candidateNeighbor = null;
		int centralIndex = doX ? point.mIndexSortedX : point.mIndexSortedY;
				
		for (int i = centralIndex + increment; 
				((i >= 0) && (i < sortedPoints.size())) 
				&& ((candidateNeighbor = sortedPoints.get(i)) != null) 
				&& (Math.abs(distanceMethodOneAxis.distance(point, candidateNeighbor)) <= theoreticalMaxDistance); 
				i += increment) {
			
			if (point.getCartesianDistance(candidateNeighbor) < mEps) {
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
