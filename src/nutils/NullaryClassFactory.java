package nutils;

import java.util.ArrayList;
import java.util.Collection;

public class NullaryClassFactory<V> {
	
	Class<V> mClass;
	
	// ========================================================================
	public NullaryClassFactory(Class<V> theClass) {
		mClass = theClass;
	}
	
	// ========================================================================
	public Class<V> getClassOfTarget() { return mClass; }
	
	// ========================================================================
	public V newInstanceRaw() throws InstantiationException, IllegalAccessException { 		
		return mClass.newInstance();
	}

	// ========================================================================
	/** The same as newInstanceRaw(), but the exceptions are caught internally. */
	public V newInstance() {
		return NullaryClassFactory.getNewInstance(mClass);
	}

	// ========================================================================
	public <Key extends Enum<Key>> EnumMapSafe<Key, V> newEnumMap(Class<Key> keyClass) {		
		EnumMapSafe<Key, V> enumMap = new EnumMapSafe<Key, V>(keyClass);
		
		Key[] keys = keyClass.getEnumConstants();
		for (Key theKey : keys) {
			enumMap.put(theKey, newInstance());
		}
		
		return enumMap;
	}
	
	// ========================================================================
	public <T extends Collection<V>> T newList(Class<T> listClass, int numElements) {
		NullaryClassFactory<T> classFactory = new NullaryClassFactory<T>(listClass);
		T newCollection = classFactory.newInstance();
				
		for (int i = 0; i < numElements; i++) {			
			newCollection.add(newInstance());
		}
		
		return newCollection;
	}
	
	// ========================================================================
	public ArrayList<V> newArrayList(int numElements) {
		ArrayList<V> newList = new ArrayList<V>(numElements);
		for (int i = 0; i < numElements; i++) {
			V newObject = newInstance();
			newList.add(newObject);
		}
		return newList;
	}

	// ========================================================================
	public static<V> V getNewInstance(Class<V> theClass) {
		V newObject = null;
		
		try {
			newObject = theClass.newInstance();
			
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
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub

	}

}
