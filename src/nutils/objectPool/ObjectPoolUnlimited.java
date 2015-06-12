package nutils.objectPool;

import java.util.concurrent.LinkedBlockingQueue;

/**
 * This abstract class offers a thread-safe framework for a pool that contains a limit of
 * Integer.MAX_VALUE number of elements, practically no limit.
 * 
 * The idea is that elements can be recycled into the pool, which
 * can then be borrowed.  If the pool is empty at the time of a borrowing request, the
 * pool will simply create a new instance and return it.  In other words, this framework
 * offers a balance between element creation and element recycling.
 * 
 * @author Ninad Dewal
 */


public abstract class ObjectPoolUnlimited<E> extends ObjectPoolAbstract<E> {
	
	// ========================================================================
	protected LinkedBlockingQueue<E> mPool;

	// ========================================================================
	/** Protected Constructor 
	 * @param poolName Name of the pool.  Can be set to null
	 */
	protected ObjectPoolUnlimited(String poolName) {
		super(poolName);
		mPool = new LinkedBlockingQueue<E>();		
	}

	// ========================================================================
	/** Gets an element from the pool.  Returns null if no elements left in pool. */
	public E borrowElement() {
		synchronized(mPool) {
			return (mPool.size() == 0) ? instantiateNewObject() : mPool.poll(); 
		}
	}

	// ========================================================================
	/** Returns an element back to the pool.  Function fails if the pool already has
	 *  the specified number of elements. Returns whether success or failure. */
	public boolean returnElement(E e) { return mPool.offer(e); }
	
	// ========================================================================
	/** Returns the current number of elements contained by the pool. */
	public int size() { return mPool.size(); }
	
	// ========================================================================
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub

	}

}
