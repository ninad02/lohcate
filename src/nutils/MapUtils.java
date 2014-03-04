package nutils;

import java.util.ArrayList;
import java.util.Map;

public class MapUtils {
	
	// ========================================================================
	public static<K, V> V getOrCreate(K theKey, NullaryClassFactory<V> classFactory, Map<K, V> theMap) {
		V resultValue = theMap.get(theKey);
		if (resultValue == null) {
			resultValue = classFactory.newInstance();
			theMap.put(theKey, resultValue);
		}
		return resultValue;
	}
	
	// ========================================================================
	public static<K, V> ArrayList<V> getOrCreateArrayList(K theKey, Map<K, ArrayList<V>> theMap) {
		ArrayList<V> resultList = theMap.get(theKey);
		if (resultList == null) {
			resultList = new ArrayList<V>();
			theMap.put(theKey, resultList);
		}
		return resultList;
	}
}
