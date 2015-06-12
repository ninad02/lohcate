package nutils.objectPool;

import java.util.ArrayList;

/**
 * This is a package-private toy test class that is meant purely for testing purposes.
 *    
 * @author Ninad Dewal
 * @version 0.0.2 
 */

// ========================================================================
class ObjectPoolLimitedString_TestClass extends ObjectPoolLimited<String> {

	// ========================================================================
	int mCurrentID = 0;
	
	// ========================================================================
	protected ObjectPoolLimitedString_TestClass(int numObjectsInPool, String poolName) {
		super(numObjectsInPool, poolName);
		mCurrentID = 0;
	}

	// ========================================================================
	@Override
	protected String instantiateNewObject() {
		synchronized(this) {
			return new String("Str_" + ++mCurrentID);
		}		
	}
	
	// ========================================================================
	public static void main(String[] args) {
		TestPool();
	}
	
	// ========================================================================
	// ========================================================================
	private static void TestPool() {
		int numElementsTotal = 100;
		
		ObjectPoolLimitedString_TestClass testPool = new ObjectPoolLimitedString_TestClass(numElementsTotal, "String Test Pool");
		
		// Print the pool size
		System.out.printf("Expected Pool Size:\t%d\tActual:\t%d\n", numElementsTotal, testPool.size());
		
		// Now, borrow all elements from the pool
		ArrayList<String> tempList = new ArrayList<String>();
		for (int i = 0; i < numElementsTotal; i++) {
			tempList.add(testPool.borrowElement());			
		}
		
		// Print the pool size
		System.out.printf("Expected Pool Size:\t%d\tActual:\t%d\n", 0, testPool.size());		
		
		// Try to borrow an element from the bool
		String extraElement = testPool.borrowElement(); 
		System.out.printf("Expected Element Value:\t%s\tActual:\t%s\n", "null", extraElement);
		
		// Cache one of the elements
		extraElement = tempList.get(0);
		
		// Return all elements
		for (String s : tempList) {
			testPool.returnElement(s);
		}
		
		// Print the pool size
		System.out.printf("Expected Pool Size:\t%d\tActual:\t%d\n", numElementsTotal, testPool.size());
		
		// Try to add the extra element
		boolean isSuccess = testPool.returnElement(extraElement);
		System.out.printf("Expected Adding Value:\t%b\tActual:\t%b\n", false, isSuccess);
	}
}
