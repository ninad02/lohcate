package nutils.objectPool;

import java.util.concurrent.ArrayBlockingQueue;

/**
 * This abstract class offers a thread-safe framework for a pool that contains a finite
 * number of elements.  
 * 
 * @author Ninad Dewal
 */

public abstract class ObjectPoolLimited<E> extends ObjectPoolAbstract<E> {
	
	// ========================================================================
	protected ArrayBlockingQueue<E> mPool;	
	
	// ========================================================================
	/** Protected Constructor
	 * @param numObjectsInPool Fixed number of objects that the pool should be set to. 
	 * @param poolName Name of the pool.  Can be set to null
	 */

	protected ObjectPoolLimited(int numObjectsInPool, String poolName) {
		super(poolName);
		mPool = new ArrayBlockingQueue<E>(numObjectsInPool, true);
		fillPool();
	}
	
	// ========================================================================
	// Private method for filling the pool with actual elements.  Should only
	// be called by the constructor.
	private void fillPool() {		
		synchronized(mPool) {
			//System.out.println(mPoolName + " Pool Remaining Capacity: " + mPool.remainingCapacity());
			int remainingCapacity = mPool.remainingCapacity();
			for (int i = 0; i < remainingCapacity; i++) {
				mPool.add(instantiateNewObject());
			}
			//System.out.println(mPoolName + " Pool Filled Size: " + mPool.size());
		}
	}

	// ========================================================================
	/** Gets an element from the pool.  Returns null if no elements left in pool. */
	public E borrowElement() { return mPool.poll(); }
	
	// ========================================================================
	/** Returns an element back to the pool.  Function fails if the pool already has
	 *  the specified number of elements. Returns whether success or failure as boolean. */
	public boolean returnElement(E e) { 		
		synchronized(mPool) {
			if (mPool.remainingCapacity() == 0) return false;  // avoids exception being thrown
			return mPool.add(e);
		} 
	}

	// ========================================================================
	/** Returns the current number of elements contained by the pool. */
	public int size() { return mPool.size(); }
}
