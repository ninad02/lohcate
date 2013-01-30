package nutils;

import java.util.ArrayList;
import java.util.Collection;

public class NullaryClassFactory<E> {
	
	Class<E> mClass;
	
	// ========================================================================
	public NullaryClassFactory(Class<E> theClass) {
		mClass = theClass;
	}
	
	// ========================================================================
	public E newInstanceRaw() throws InstantiationException, IllegalAccessException { 		
		return mClass.newInstance();
	}

	// ========================================================================
	/** The same as newInstanceRaw(), but the exceptions are caught internally. */
	public E newInstance() {
		E newObject = null;
		
		try {
			newObject = newInstanceRaw();
			
		} catch (IllegalAccessException e) {			
			e.printStackTrace();
			System.exit(-1);
			
		} catch (InstantiationException e) {
			e.printStackTrace();
			System.exit(-1);		
		}
		
		return newObject;
	}

	// ========================================================================
	public static<V, T extends Collection<V>> T newList(Class<T> listClass, NullaryClassFactory<V> elementFactory, int numElements) {
		NullaryClassFactory<T> classFactory = new NullaryClassFactory<T>(listClass);
		T newCollection = classFactory.newInstance();
				
		for (int i = 0; i < numElements; i++) {
			V newObject = elementFactory.newInstance();
			newCollection.add(newObject);
		}
		
		return newCollection;
	}
	
	// ========================================================================
	public ArrayList<E> newArrayList(int numElements) {
		ArrayList<E> newList = new ArrayList<E>(numElements);
		for (int i = 0; i < numElements; i++) {
			E newObject = newInstance();
			newList.add(newObject);
		}
		return newList;
	}

	// ========================================================================
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub

	}

}
