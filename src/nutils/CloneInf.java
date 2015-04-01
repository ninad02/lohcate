package nutils;

/**
 * This class is meant for enforcing a class to be able to clone itself. 
 * @author Ninad Dewal
 *
 *
 * @param <T> Any type
 */

public interface CloneInf<T> {
	
	/** Returns a deep copy. */
	public T makeClone();
	
	/** Returns either a deep or shallow copy, depending on the argument. */
	public T makeClone(boolean deepCopy);
}
