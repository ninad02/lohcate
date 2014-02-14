/**
 * 
 */
package nutils;

/**
 * @author ninad
 *
 */
public abstract class PrimitiveWrapper<T extends PrimitiveWrapper<T>> implements Comparable<T> {
		
	protected static final byte DefaultValue = 0;

	//=========================================================================
	public abstract boolean equals(Object o);
	public abstract int hashCode();
	
	//=========================================================================
	/** Provides a wrapper for a primitive int value that can be changed.  Useful
	 *  when you want to have a changeable integer within a container. */
	public static class WInteger extends PrimitiveWrapper<WInteger> {
		public static final NullaryClassFactory<WInteger> ClassFactory = new NullaryClassFactory<WInteger>(WInteger.class);
		
		public int mInt;
		public WInteger(int i) { mInt = i; }
		public WInteger()      { this(DefaultValue); }
		
		public int compareTo(WInteger rhs) {
			return Integer.compare(mInt, rhs.mInt);
		}
		/* (non-Javadoc)
		 * @see java.lang.Object#hashCode()
		 */
		@Override
		public int hashCode() {
			final int prime = 31;
			int result = 1;
			result = prime * result + mInt;
			return result;
		}
		
		/* (non-Javadoc)
		 * @see java.lang.Object#equals(java.lang.Object)
		 */
		@Override
		public boolean equals(Object obj) {
			if (this == obj)
				return true;
			if (obj == null)
				return false;
			if (getClass() != obj.getClass())
				return false;
			WInteger other = (WInteger) obj;
			if (mInt != other.mInt)
				return false;
			return true;
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
		
		/* (non-Javadoc)
		 * @see java.lang.Object#hashCode()
		 */
		@Override
		public int hashCode() {
			final int prime = 31;
			int result = 1;
			long temp;
			temp = Double.doubleToLongBits(mDouble);
			result = prime * result + (int) (temp ^ (temp >>> 32));
			return result;
		}

		/* (non-Javadoc)
		 * @see java.lang.Object#equals(java.lang.Object)
		 */
		@Override
		public boolean equals(Object obj) {
			if (this == obj)
				return true;
			if (obj == null)
				return false;
			if (getClass() != obj.getClass())
				return false;
			WDouble other = (WDouble) obj;
			if (Double.doubleToLongBits(mDouble) != Double
					.doubleToLongBits(other.mDouble))
				return false;
			return true;
		}
		
		public int compareTo(WDouble rhs) {
			return Double.compare(mDouble, rhs.mDouble);
		}
		
		
	}

	//=========================================================================
	public static class WFloat extends PrimitiveWrapper<WFloat> {
		public float mFloat;
		public WFloat(float f) { mFloat = f; }
		
		/* (non-Javadoc)
		 * @see java.lang.Object#hashCode()
		 */
		@Override
		public int hashCode() {
			final int prime = 31;
			int result = 1;
			result = prime * result + Float.floatToIntBits(mFloat);
			return result;
		}

		/* (non-Javadoc)
		 * @see java.lang.Object#equals(java.lang.Object)
		 */
		@Override
		public boolean equals(Object obj) {
			if (this == obj)
				return true;
			if (obj == null)
				return false;
			if (getClass() != obj.getClass())
				return false;
			WFloat other = (WFloat) obj;
			if (Float.floatToIntBits(mFloat) != Float
					.floatToIntBits(other.mFloat))
				return false;
			return true;
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

		/* (non-Javadoc)
		 * @see java.lang.Object#hashCode()
		 */
		@Override
		public int hashCode() {
			final int prime = 31;
			int result = 1;
			result = prime * result + mByte;
			return result;
		}

		/* (non-Javadoc)
		 * @see java.lang.Object#equals(java.lang.Object)
		 */
		@Override
		public boolean equals(Object obj) {
			if (this == obj)
				return true;
			if (obj == null)
				return false;
			if (getClass() != obj.getClass())
				return false;
			WByte other = (WByte) obj;
			if (mByte != other.mByte)
				return false;
			return true;
		}
		
		
	}

	//=========================================================================
	public static class WBoolean extends PrimitiveWrapper<WBoolean> {
		public boolean mBoolean;		
		public WBoolean(boolean b) { mBoolean = b; }
		
		public int compareTo(WBoolean rhs) {
			return Boolean.compare(mBoolean, rhs.mBoolean);
		}

		/* (non-Javadoc)
		 * @see java.lang.Object#hashCode()
		 */
		@Override
		public int hashCode() {
			final int prime = 31;
			int result = 1;
			result = prime * result + (mBoolean ? 1231 : 1237);
			return result;
		}

		/* (non-Javadoc)
		 * @see java.lang.Object#equals(java.lang.Object)
		 */
		@Override
		public boolean equals(Object obj) {
			if (this == obj)
				return true;
			if (obj == null)
				return false;
			if (getClass() != obj.getClass())
				return false;
			WBoolean other = (WBoolean) obj;
			if (mBoolean != other.mBoolean)
				return false;
			return true;
		}
		
		
	}

	//=========================================================================
	public static class WLong extends PrimitiveWrapper<WLong> {
		public long mLong;
		public WLong(long l) { mLong = l; }
		
		public int compareTo(WLong rhs) {
			return Long.compare(mLong, rhs.mLong);
		}

		/* (non-Javadoc)
		 * @see java.lang.Object#hashCode()
		 */
		@Override
		public int hashCode() {
			final int prime = 31;
			int result = 1;
			result = prime * result + (int) (mLong ^ (mLong >>> 32));
			return result;
		}

		/* (non-Javadoc)
		 * @see java.lang.Object#equals(java.lang.Object)
		 */
		@Override
		public boolean equals(Object obj) {
			if (this == obj)
				return true;
			if (obj == null)
				return false;
			if (getClass() != obj.getClass())
				return false;
			WLong other = (WLong) obj;
			if (mLong != other.mLong)
				return false;
			return true;
		}
		
		
	}

	//=========================================================================
	public static class WChar extends PrimitiveWrapper<WChar> {
		public char mChar;
		public WChar(char ch) { mChar = ch; }
		
		public int compareTo(WChar rhs) {
			return Character.compare(mChar, rhs.mChar);
		}

		/* (non-Javadoc)
		 * @see java.lang.Object#hashCode()
		 */
		@Override
		public int hashCode() {
			final int prime = 31;
			int result = 1;
			result = prime * result + mChar;
			return result;
		}

		/* (non-Javadoc)
		 * @see java.lang.Object#equals(java.lang.Object)
		 */
		@Override
		public boolean equals(Object obj) {
			if (this == obj)
				return true;
			if (obj == null)
				return false;
			if (getClass() != obj.getClass())
				return false;
			WChar other = (WChar) obj;
			if (mChar != other.mChar)
				return false;
			return true;
		}
		
		
	}

	//=========================================================================
	public static class WString extends PrimitiveWrapper<WString> {
		public String mStr;
		public WString(String s) { mStr = s; }
		
		public int compareTo(WString rhs) {
			return mStr.compareTo(rhs.mStr);
		}

		/* (non-Javadoc)
		 * @see java.lang.Object#hashCode()
		 */
		@Override
		public int hashCode() {
			final int prime = 31;
			int result = 1;
			result = prime * result + ((mStr == null) ? 0 : mStr.hashCode());
			return result;
		}

		/* (non-Javadoc)
		 * @see java.lang.Object#equals(java.lang.Object)
		 */
		@Override
		public boolean equals(Object obj) {
			if (this == obj)
				return true;
			if (obj == null)
				return false;
			if (getClass() != obj.getClass())
				return false;
			WString other = (WString) obj;
			if (mStr == null) {
				if (other.mStr != null)
					return false;
			} else if (this.compareTo(other) != 0)
				return false;
			return true;
		}
		
		
	}

	//=========================================================================
	public static class WObject<E extends Object & Comparable<? super E>> extends PrimitiveWrapper<WObject<E>> {
		public E mObj;
		public WObject(E o) { mObj = o; }
		
		public int compareTo(WObject<E> rhs) {
			return mObj.compareTo(rhs.mObj);
		}

		/* (non-Javadoc)
		 * @see java.lang.Object#hashCode()
		 */
		@Override
		public int hashCode() {
			final int prime = 31;
			int result = 1;
			result = prime * result + ((mObj == null) ? 0 : mObj.hashCode());
			return result;
		}

		/* (non-Javadoc)
		 * @see java.lang.Object#equals(java.lang.Object)
		 */
		@Override
		public boolean equals(Object obj) {
			if (this == obj)
				return true;
			if (obj == null)
				return false;
			if (getClass() != obj.getClass())
				return false;
			WObject other = (WObject) obj;
			if (mObj == null) {
				if (other.mObj != null)
					return false;
			} else if (!mObj.equals(other.mObj))
				return false;
			return true;
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
