package nutils;

public class UtilsBasic {

	public UtilsBasic() {
		// TODO Auto-generated constructor stub
	}
	
	// ========================================================================
	public static<T> boolean isNull(T t)    { return (t == null); }
		
	// ========================================================================
	public static<T> boolean isNotNull(T t) { return (t != null); }

	// ========================================================================
	public static<E extends Enum<E>> E getEnumValue(String enumName, Class<E> enumType) {
		E theEnum = null;
		try {
			theEnum = Enum.valueOf(enumType, enumName);
		} catch (IllegalArgumentException e) {}
		
		return theEnum;
	}
	
	// ========================================================================
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub

	}

}
