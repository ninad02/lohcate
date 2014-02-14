package genomeUtils;

import java.util.ArrayList;

import nutils.PrimitiveWrapper;

import com.carrotsearch.hppc.IntArrayList;

public class ElementPlaneSplit<E> {

	public static final int InvalidPlaneID = -1;
	
	protected ArrayList<PlaneElements<E>> mPlanes;
	
	protected IntArrayList mPlaneID;
	protected IntArrayList mIndexInPlane;

	// ========================================================================
	public ElementPlaneSplit(int numPoints, int numPlanes) {
		mPlaneID      = new IntArrayList(numPoints);
		mIndexInPlane = new IntArrayList(numPoints);
		
		mPlanes = new ArrayList<PlaneElements<E>>();
		for (int i = 0; i < numPlanes; i++) {
			mPlanes.add(new PlaneElements<E>(numPoints));
		}
	}

	// ========================================================================
	public int getNumPlanes() { return mPlanes.size(); }
	
	// ========================================================================
	public void registerElement(int planeID, E element, int indexOfElementInMainList) {
		if (planeID < 0) {
			mPlaneID.add(InvalidPlaneID);
			mIndexInPlane.add(Integer.MAX_VALUE);
		} else {
			int indexInPlane = mPlanes.get(planeID).registerElement(element, indexOfElementInMainList);
			mPlaneID.add(planeID);
			mIndexInPlane.add(indexInPlane);
		}
	}
	
	// ========================================================================
	public ArrayList<E> getPointsOnPlane(int planeID) {
		return mPlanes.get(planeID).mList;
	}
	
	// ========================================================================
	public void clear() {
		for (PlaneElements<E> plane : mPlanes) {
			plane.clear();
		}
	}
	
	// ========================================================================
	public int getIndexOfPlaneElementInMainList(int planeID, int elementIndexInPlane) {
		return mPlanes.get(planeID).mIndexToMainList.get(elementIndexInPlane);
	}
	
	// ========================================================================
	public int getPlane(int indexOfElementInMainList, PrimitiveWrapper.WInteger indexInPlane) {
		indexInPlane.mInt = mIndexInPlane.get(indexOfElementInMainList);
		return mPlaneID.get(indexOfElementInMainList);
	}

	// ========================================================================
	// INNER CLASS
	// ========================================================================
	public static class PlaneElements<E> {
		
		ArrayList<E> mList;
		IntArrayList mIndexToMainList;
		
		public PlaneElements(int numEstimatedElements) {
			mList = new ArrayList<E>(numEstimatedElements);
			mIndexToMainList = new IntArrayList(numEstimatedElements);
		}
		
		public int registerElement(E element, int indexOfElementInMainList) {
			mList.add(element);
			mIndexToMainList.add(indexOfElementInMainList);
			return mList.size() - 1; // Returns the index of the added element
		}
		
		public void clear() {
			mList.clear();
			mIndexToMainList.clear();
		}
	}
	
	// ========================================================================
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub

	}

}
