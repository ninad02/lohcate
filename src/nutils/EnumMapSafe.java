package nutils;

import java.util.ArrayList;
import java.util.EnumMap;
import java.util.Map;

public class EnumMapSafe<K extends Enum<K>, V> extends EnumMap<K, V> {
	
	/**
	 * 
	 */
	private static final long serialVersionUID = 6320062537394900382L;

	// ========================================================================
	/** Creates an empty enum map with the specified key type. */
	public EnumMapSafe(Class<K> keyType) {
		super(keyType);
	}

	// ========================================================================
	/** Creates an enum map with the same key type as the specified enum map, initially containing the same mappings (if any). */
	public EnumMapSafe(EnumMapSafe<K, ? extends V> m) {
		super(m);
	}
	
	// ========================================================================
	/** Creates an enum map initialized from the specified map. */
	public EnumMapSafe(Map<K, ? extends V> m) {
		super(m);
	}
	
	// ========================================================================
	public V get(K key) {
		return super.get(key);
	}
	
	// ========================================================================
	@Deprecated	
	public V get(Object o) throws UnsupportedOperationException {
		((new Exception("ERROR: Unsupported EnumMapSafe operation: get(Object)!"))).printStackTrace();
		System.exit(-1);
		return null;
	}

	// ========================================================================
	public static <E extends Enum<E>, V> EnumMapSafe<E, ArrayList<V>> createEnumMapOfArrayLists(Class<E> enumClass) {
		EnumMapSafe<E, ArrayList<V>> newMap = new EnumMapSafe<E, ArrayList<V>>(enumClass);
		E[] enumValues = enumClass.getEnumConstants();
		
		for (E enumValue : enumValues) {
			newMap.put(enumValue, new ArrayList<V>());
		}
		
		return newMap;
	}
	
	// ========================================================================
	public static <E extends Enum<E>, V> EnumMapSafe<E, ArrayList<V>> createEnumMapOfArrayLists(Class<E> enumClass, Class<V> arrayTargetObject) {
		return createEnumMapOfArrayLists(enumClass);
	}

	// ========================================================================
	/** Creates an EnumMap and initializes it with (maps all keys with) the given value. */
	public static <E extends Enum<E>, V> EnumMapSafe<E, V> createEnumMap(Class<E> enumClass, V initialValue) {
		EnumMapSafe<E, V> newMap = new EnumMapSafe<E, V>(enumClass);
		E[] enumValues = enumClass.getEnumConstants();
		
		for (E enumValue : enumValues) {
			newMap.put(enumValue, initialValue);
		}
		
		return newMap;
	}

	// ========================================================================
	/** An inner static enum class for testing purposes. */
	private static enum MyNuc {
		A, B, C, D;
	}
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		EnumMapSafe<MyNuc, Integer> map = new EnumMapSafe<>(MyNuc.class);
		for (MyNuc nuc : MyNuc.values()) {
			map.put(nuc, nuc.ordinal() * 2);
		}
		
		for (MyNuc nuc : MyNuc.values()) {
			System.out.println(map.get(nuc));
		}

		System.out.println(map.equals(map));
		//EnumMapSafe<MyNuc, Integer> map2 = map.clone();


	}

}
