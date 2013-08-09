package nutils;

import java.util.Set;

public class EnumMatrix<E extends Enum<E>, V> {

	//=========================================================================
	public static interface MatrixTraverseAction<E extends Enum<E>, V> {
		public void takeAction(EnumMatrix<E, V> matrix, E row, E column, V value);
	}

	//=========================================================================
	private MatrixTraverseAction<E, V> mTraversePut = new MatrixTraverseAction<E, V>() {
		public void takeAction(EnumMatrix<E, V> matrix, E row, E column, V value) {
			matrix.put(row, column, value);
		}
	};
	
	//=========================================================================
	protected EnumMapSafe<E, EnumMapSafe<E, V>> mMatrix;

	//=========================================================================
	public EnumMatrix(Class<E> keyType) {
		this(keyType, null);
	}
	
	//=========================================================================
	public EnumMatrix(Class<E> keyType, V initialValue) {
		mMatrix = new EnumMapSafe<E, EnumMapSafe<E, V>>(keyType);
		allocateRowElements(keyType);
		clear(initialValue);		
	}
	
	//=========================================================================
	private void allocateRowElements(Class<E> keyType) {
		E[] enumConstants = keyType.getEnumConstants();
		for (E enumConstant : enumConstants) {
			mMatrix.put(enumConstant, new EnumMapSafe<E, V>(keyType));
		}				
	}
	
	//=========================================================================
	public Boolean put(E row, E column, V value) {
		if (row    == null) return null;
		if (column == null) return Boolean.FALSE;
		
		mMatrix.get(row).put(column, value);
		return Boolean.TRUE;
	}
	
	//=========================================================================
	public V get(E row, E column) {
		return mMatrix.get(row).get(column);
	}

	//=========================================================================
	/** Resets all elements to null. */
	public void clear() { clear(null); }
	
	//=========================================================================
	/** Resets all elements to the specified value. */
	public void clear(V value) {
		Set<E> enumKeys = mMatrix.keySet();
		for (E enumKeyRow : enumKeys) {
			for (E enumKeyCol : enumKeys) {
				put(enumKeyRow, enumKeyCol, value);
			}
		}
	}
	
	//=========================================================================
	protected void traverse(MatrixTraverseAction<E, V> traverseAction) {
		
	}

	//=========================================================================
	//public static
	
	//=========================================================================
	//=========================================================================
	//=========================================================================
	//=========================================================================
	//=========================================================================
	//=========================================================================
	//=========================================================================
	private static enum TestEnums { A, B, C, D; }
	
	//=========================================================================
	private static void TestEnumMatrix() {
		EnumMatrix<TestEnums, Integer> matrix = new EnumMatrix<>(TestEnums.class, -1);
		int counter = 10;
		for (TestEnums theEnumRow : TestEnums.values()) {
			for (TestEnums theEnumCol : TestEnums.values()) {
				matrix.put(theEnumRow, theEnumCol, ++counter);
			}
		}
		

		for (TestEnums theEnumRow : TestEnums.values()) {
			for (TestEnums theEnumCol : TestEnums.values()) {
				System.out.printf("%s\t", matrix.get(theEnumRow, theEnumCol));
			}
			System.out.println();
		}

	}
	
	//=========================================================================
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		TestEnumMatrix();
	}

}
