package nutils.counter;

import java.util.HashMap;
import java.util.Map.Entry;
import java.util.Set;

import nutils.PrimitiveWrapper;
import nutils.PrimitiveWrapper.WInteger;

public class ObjectCounter<E> {
	
	HashMap<E, PrimitiveWrapper.WInteger> mCounterMap;
	
	// ========================================================================
	public ObjectCounter() {
		mCounterMap = new HashMap<>();
	}

	// ========================================================================
	/** Looks up the object in the counter.  If it exists, then the count is 
	 * incremented, and the new count is returned.  If it does not exist, 1
	 * is returned.
	 * @return (existing count + 1) if exists, 1 if not
	 */
	public int increment(E theObject) {
		PrimitiveWrapper.WInteger theCount = mCounterMap.get(theObject);
		if (theCount == null) {
			theCount = new PrimitiveWrapper.WInteger(1);			
			mCounterMap.put(theObject, theCount);			
		} else {
			++theCount.mInt;
		}
		return theCount.mInt;
	}
	
	// ========================================================================
	/** Returns the objects and counts as a set.
	 * @return 
	 * 
	 */
	public Set<Entry<E, PrimitiveWrapper.WInteger>> getAllCounts() {
		return mCounterMap.entrySet();
	}
	

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub

	}

}
