package nutils.objectPool;

/**
 * This is the interface that determines the basic methods for an object pool design pattern.  This does not necessitate thread-safety and
 * leaves thread-safety to sublcasses, thus allowing flexibility of avoiding threads
 * if the use-case is a single threaded program.
 *    
 * @author Ninad Dewal
 * @version 0.0.2 
 */

public interface ObjectPool<E> {

	// ========================================================================
	/** Gets an element from the pool.  Returns null if no elements left in pool. */
	public E borrowElement();

	// ========================================================================
	/** Returns an element back to the pool.  Function fails if the pool already has
	 *  the specified number of elements. Returns whether success or failure. */
	public boolean returnElement(E e);
	
	// ========================================================================
	/** Returns the current number of elements contained by the pool. */
	public int size();	
	
}
