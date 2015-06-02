package nutils;

import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;

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
	public static interface ProcedureToBeTimed {
		public void takeAction();
	}
	
	// ========================================================================
	public static long getTimeTaken(ProcedureToBeTimed p) {
		long timeStart = System.currentTimeMillis();
		p.takeAction();
		long timeEnd = System.currentTimeMillis();
		return (timeEnd - timeStart);		
	}
	
	// ========================================================================
	@SuppressWarnings("unchecked")
	public static<E extends Cloneable> E getClone(E obj) {
		try {
			Method m = obj.getClass().getMethod("clone");
			return (E) m.invoke(obj);
		} catch (NoSuchMethodException | SecurityException | IllegalAccessException | IllegalArgumentException | InvocationTargetException e) {
			e.printStackTrace();
			System.exit(-1);
		}
		return null;
	}
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub

	}

}
