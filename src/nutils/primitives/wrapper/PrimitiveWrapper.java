/**
 * 
 */
package nutils.primitives.wrapper;

import java.io.Serializable;

import nutils.CloneInf;
import nutils.NullaryClassFactory;


/**
 * @author ninad
 *
 */
public abstract class PrimitiveWrapper<T extends PrimitiveWrapper<T>> implements Comparable<T>, Serializable {
		
	/**
	 * 
	 */
	private static final long serialVersionUID = 487200172011561389L;
	protected static final byte DefaultValue = 0;

	//=========================================================================
	public abstract boolean equals(Object o);
	public abstract int hashCode();
	
	//=========================================================================
	/** Provides a wrapper for a primitive int value that can be changed.  Useful
	 *  when you want to have a changeable integer within a container. */
	public static class WInteger extends PrimitiveWrapper<WInteger> implements CloneInf<WInteger> {
		/**
		 * 
		 */
		private static final long serialVersionUID = 5364611144472487047L;

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
		
		@Override
		public WInteger makeClone() { return makeClone(true); }
		
		@Override
		public WInteger makeClone(boolean deepCopy) { return new WInteger(this.mInt); }
		
		
	}
	
	//=========================================================================
	/** A class that sets an initial value, with an option of setting it to a new
	 *  value.  If rejects any attempt to set the value to a second new value.
	 */
	public static class WIntegerSetOnce extends WInteger {
		
		/**
		 * 
		 */
		private static final long serialVersionUID = 3513153186232452102L;

		public static final NullaryClassFactory<WIntegerSetOnce> ClassFactory = new NullaryClassFactory<WIntegerSetOnce>(WIntegerSetOnce.class);
		
		private boolean mAlreadySetOnce;
		
		public WIntegerSetOnce() { this(DefaultValue); }
		
		public WIntegerSetOnce(int i) {
			super(i);
			reset();
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
		
		/** Method the resets the state of whether the integer has been set or not.  
		 *  This method allows the object to be reused, while enforcing the one-set 
		 *  restriction, as the user needs to make an explicitly call to this function
		 *  to reset the one-set restriction.
		 */
		public void reset() { mAlreadySetOnce = false; }	
		
		public void resetAndSet(final int i) {
			reset();
			set(i);
		}
	}

	//=========================================================================
	public static class WDouble extends PrimitiveWrapper<WDouble> {
		/**
		 * 
		 */
		private static final long serialVersionUID = -1695323296329779788L;
		public double mDouble;
		public WDouble(double d) { mDouble = d; }
		public WDouble() { this(DefaultValue); }
		
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
	/** A class that sets an initial value, with an option of setting it to a new
	 *  value.  If rejects any attempt to set the value to a second new value.
	 */
	public static class WDoubleSetOnce extends WDouble {
		/**
		 * 
		 */
		private static final long serialVersionUID = -7701062470079569861L;
		private boolean mAlreadySetOnce;
		
		public WDoubleSetOnce() {
			this(DefaultValue);
		}
		
		public WDoubleSetOnce(final double d) {
			super(d);
			reset();
		}
		
		public boolean set(final double d) {			
			if (mAlreadySetOnce) {
				return (d == mDouble);
			} else {
				mDouble = d;
				return (mAlreadySetOnce = true);
			}
		}
		
		public boolean isAlreadySetOnce() { return mAlreadySetOnce; }
		
		/** Method the resets the state of whether the integer has been set or not.  
		 *  This method allows the object to be reused, while enforcing the one-set 
		 *  restriction, as the user needs to make an explicitly call to this function
		 *  to reset the one-set restriction.
		 */
		public void reset() { mAlreadySetOnce = false; }	
		
		public void resetAndSet(final double d) {
			reset();
			set(d);
		}
	}

	//=========================================================================
	public static class WFloat extends PrimitiveWrapper<WFloat> {
		/**
		 * 
		 */
		private static final long serialVersionUID = -5959443692165371017L;
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
		/**
		 * 
		 */
		private static final long serialVersionUID = 1371871778317109403L;
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
		/**
		 * 
		 */
		private static final long serialVersionUID = 4062487221340800789L;
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
		/**
		 * 
		 */
		private static final long serialVersionUID = -6008319569123390431L;
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
		/**
		 * 
		 */
		private static final long serialVersionUID = -7489775918672747877L;
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
		/**
		 * 
		 */
		private static final long serialVersionUID = 3766531879207274819L;
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
	public static class WObject<E extends Object & Comparable<E>> extends PrimitiveWrapper<WObject<E>> {
		/**
		 * 
		 */
		private static final long serialVersionUID = 4222977275080522829L;
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
			WObject<?> other = (WObject<?>) obj;
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
