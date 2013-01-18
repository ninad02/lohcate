/**
 * 
 */
package nutils;

/**
 * @author ninad
 *
 */
public abstract class PrimitiveWrapper<T extends PrimitiveWrapper<T>> implements Comparable<T> {
		
	//=========================================================================
	/** Provides a wrapper for a primitive int value that can be changed.  Useful
	 *  when you want to have a changeable integer within a container. */
	public static class WInteger extends PrimitiveWrapper<WInteger> {
		public int mInt;
		public WInteger(int i) { mInt = i; }
		
		public int compareTo(WInteger rhs) {
			return Integer.compare(mInt, rhs.mInt);
		}
	}
	
	//=========================================================================
	/** A class that sets an initial value, with an option of setting it to a new
	 *  value.  If rejects any attempt to set the value to a second new value.
	 */
	public static class WIntegerSetOnce extends WInteger {
		private boolean mAlreadySetOnce;
		
		public WIntegerSetOnce(int i) {
			super(i);
			mAlreadySetOnce = false;
		}
		
		public boolean set(final int i) {			
			if (mAlreadySetOnce) {
				return (i == mInt);
			} else {
				mInt = i;
				return (mAlreadySetOnce = true);
			}
		}
		
		public boolean isAlreadySetOnce() { return mAlreadySetOnce; }
	}

	//=========================================================================
	public static class WDouble extends PrimitiveWrapper<WDouble> {
		public double mDouble;
		public WDouble(double d) { mDouble = d; }
		
		public boolean equals(Object o) {			
			return mDouble == ((WDouble) o).mDouble;
		}
		
		public int compareTo(WDouble rhs) {
			return Double.compare(mDouble, rhs.mDouble);
		}
	}

	//=========================================================================
	public static class WFloat extends PrimitiveWrapper<WFloat> {
		public float mFloat;
		public WFloat(float f) { mFloat = f; }
		
		public boolean equals(Object o) {
			return mFloat == ((WFloat) o).mFloat;
		}
		
		public int compareTo(WFloat rhs) {
			return Float.compare(mFloat, rhs.mFloat);
		}
	}

	//=========================================================================
	public static class WByte extends PrimitiveWrapper<WByte> {
		public byte mByte;
		public WByte(byte b) { mByte = b; }
		
		public int compareTo(WByte rhs) {
			return Byte.compare(mByte, rhs.mByte);
		}
	}

	//=========================================================================
	public static class WBoolean extends PrimitiveWrapper<WBoolean> {
		public boolean mBoolean;		
		public WBoolean(boolean b) { mBoolean = b; }
		
		public int compareTo(WBoolean rhs) {
			return Boolean.compare(mBoolean, rhs.mBoolean);
		}
	}

	//=========================================================================
	public static class WLong extends PrimitiveWrapper<WLong> {
		public long mLong;
		public WLong(long l) { mLong = l; }
		
		public int compareTo(WLong rhs) {
			return Long.compare(mLong, rhs.mLong);
		}
	}

	//=========================================================================
	public static class WChar extends PrimitiveWrapper<WChar> {
		public char mChar;
		public WChar(char ch) { mChar = ch; }
		
		public int compareTo(WChar rhs) {
			return Character.compare(mChar, rhs.mChar);
		}
	}

	//=========================================================================
	public static class WString extends PrimitiveWrapper<WString> {
		public String mStr;
		public WString(String s) { mStr = s; }
		
		public int compareTo(WString rhs) {
			return mStr.compareTo(rhs.mStr);
		}
	}

	//=========================================================================
	public static class WObject<E extends Object & Comparable<? super E>> extends PrimitiveWrapper<WObject<E>> {
		public E mObj;
		public WObject(E o) { mObj = o; }
		
		public int compareTo(WObject<E> rhs) {
			return mObj.compareTo(rhs.mObj);
		}
	}
	
	//=========================================================================
	public static void TestWIntSetOnce() {
		WIntegerSetOnce wiso = new WIntegerSetOnce(-5);
		System.out.println(wiso.mInt);
		System.out.println(wiso.set(2) + "\t" + wiso.mInt);
		System.out.println(wiso.set(3) + "\t" + wiso.mInt);
		System.out.println(wiso.set(2) + "\t" + wiso.mInt);
	}
	
	//=========================================================================
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		TestWIntSetOnce();
	}
}
