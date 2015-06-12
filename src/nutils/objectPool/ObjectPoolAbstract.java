package nutils.objectPool;

/**
 * This is the abstract class that implements ObjectPool basic methods for the 
 * object pool design pattern.  This does not necessitate thread safety and 
 * leaves thread-safety to sublcasses, thus allowing flexibility of avoiding 
 * threads if the use-case is a single threaded program.
 *    
 * @author Ninad Dewal
 * @version 0.0.2 
 */

public abstract class ObjectPoolAbstract<E> implements ObjectPool<E> {

	// ========================================================================
	protected String mPoolName = "";

	// ========================================================================
	protected ObjectPoolAbstract() {
		this("");
	}
	
	// ========================================================================
	protected ObjectPoolAbstract(String poolName) {
		setPoolName(poolName);
	}
	
	// ========================================================================
	/** Sets the name of the pool. If null string is passed, the pool name is a blank string. */
	public void setPoolName(String poolName) {
		mPoolName = (poolName == null) ? mPoolName : poolName;
	}
	
	// ========================================================================
	/** Returns the name of the pool. */
	public String getPoolName() { return mPoolName; }
	
	// ========================================================================
	/** This forces a subclass to implement a method to instantiate a new element object. 
	 *  This was designed such that class E doesn't need to implement a nullary constructor. 
	 */
	protected abstract E instantiateNewObject();

}
