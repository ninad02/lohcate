package nutils.counter;

import java.io.PrintStream;
import java.util.HashMap;
import java.util.Map.Entry;
import java.util.Set;

import nutils.primitives.wrapper.PrimitiveWrapper;

public class ObjectCounter<E> {
	
	HashMap<E, PrimitiveWrapper.WInteger> mCounterMap;
	
	// ========================================================================
	public ObjectCounter() {
		mCounterMap = new HashMap<>();
	}

	// ========================================================================
	public void clear() { 
		mCounterMap.clear();
	}
	
	// ========================================================================
	/** Looks up the object in the counter.  If it exists, then the count is 
	 * incremented, and the new count is returned.  If it does not exist, 1
	 * is returned.
	 * @return (existing count + 1) if exists, 1 if not
	 */
	public int increment(E theObject) {
		return increment(theObject, 1);
	}
	
	// ========================================================================
	/** Looks up the object in the counter.  If it exists, then the count is increment
	 *  by the amount specified, and the new count is returned.  If not, then the 
	 *  increment amount is returned
	 *  @ return (existing count + increment) if exists, increment if not
	 */
	public int increment(E theObject, int incrementAmount) {
		PrimitiveWrapper.WInteger theCount = mCounterMap.get(theObject);
		if (theCount == null) {
			theCount = new PrimitiveWrapper.WInteger(incrementAmount);			
			mCounterMap.put(theObject, theCount);			
		} else {
			theCount.mInt += incrementAmount;
		}
		return theCount.mInt;
	}
	
	// ========================================================================
	/** Given another counter, this adds the counts from that counter to this one. */
	public void increment(ObjectCounter<E> rhs) {
		for (E key : rhs.mCounterMap.keySet()) {
			this.increment(key, rhs.getCount(key));
		}
	}

	// ========================================================================
	public boolean contains(E theObject) { return (getCount(theObject) > 0); }
	
	// ========================================================================
	public int getCount(E theObject) {
		PrimitiveWrapper.WInteger theCount = mCounterMap.get(theObject);
		return (theCount == null) ? 0 : theCount.mInt; 
	}
	
	// ========================================================================
	public void print(PrintStream out, String delimiter) {
		for (E key : mCounterMap.keySet()) {
			PrimitiveWrapper.WInteger theCount = mCounterMap.get(key);
			out.printf("%s%s%d\n", key.toString(), delimiter, theCount.mInt);			
		}
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
