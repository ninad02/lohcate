package genomeUtils;

/**
 * This is a convenience class used when traversing a list.  The 
 * purpose is to implement a pattern of testing the current object
 * against a previous object, and indicating whether there's a change.
 * If there's a change, the previous object is set to the current object
 * thereafter.
 *  
 * @author Ninad Dewal
 *
 */
public class ObjectWalkerTracker<E extends Comparable<E>> {
	
	// ========================================================================
	// ENUM: Represents the type of change or no change
	// ========================================================================
	public static enum ChangeType {
		Changed_Null_To_NotNull(true),
		Changed_NotNull_To_Null(true),
		Changed_BothNotNull_IncreasingOrder(true),
		Changed_BothNotNull_DecreasingOrder(true),
		Unchanged_BothNotNull(false),
		Unchanged_BothNull(false)
		;
		
		boolean mIsChanged = false;
		
		private ChangeType(boolean isChangedIntrinsic) {
			mIsChanged = isChangedIntrinsic;
		}
		
		public boolean hasChanged() { return mIsChanged; }
	}
	
	
	// ========================================================================
	protected E mObjPrev;
	protected E mObjInitial;

	// ========================================================================
	public ObjectWalkerTracker() {
		this(null);
	}
	
	// ========================================================================
	public ObjectWalkerTracker(E initialObject) {
		mObjPrev = mObjInitial = initialObject;
	}
	
	// ========================================================================
	public void clear() {
		mObjPrev = mObjInitial;
	}
	
	// ========================================================================
	public E getPrevious() { return mObjPrev; }
	
	// ========================================================================
	public boolean isPreviousUninitialized() { return (mObjPrev == null); }
	
	// ========================================================================
	public ChangeType hasChanged(E nextObj, boolean setPreviousObjectToNewOne) {
		boolean prevNull = isPreviousUninitialized();
		boolean currNull = (nextObj == null);
		ChangeType changeType = null;
		
		if (prevNull && currNull) {
			changeType = ChangeType.Unchanged_BothNull;
		} else if (prevNull && !currNull) {
			changeType = ChangeType.Changed_Null_To_NotNull;
		} else if (!prevNull && currNull) {
			changeType = ChangeType.Changed_NotNull_To_Null;
		} else {
			// Both are not null.  Must check now if they are the same
			int result = mObjPrev.compareTo(nextObj);
			if (result == 0) {
				changeType = ChangeType.Unchanged_BothNotNull;
			} else if (result < 0) {
				changeType = ChangeType.Changed_BothNotNull_IncreasingOrder;
			} else {
				changeType = ChangeType.Changed_BothNotNull_DecreasingOrder;
			}
		}
		
		// Now, finally set the previous object to the current one
		if (setPreviousObjectToNewOne) {
			mObjPrev = nextObj;
		}
		
		// Return the result
		return changeType;
	}

	// ========================================================================
	public static enum TestEnum { A, B, C, D, E };
	
	// ========================================================================
	protected static void TestThisClass() {
		ObjectWalkerTracker<TestEnum> tracker = new ObjectWalkerTracker<TestEnum>();
		
		for (TestEnum testEnum : TestEnum.values()) {
			ChangeType ct = tracker.hasChanged(testEnum, true);
			System.out.println(ct);
		}
		
		System.out.println(tracker.hasChanged(TestEnum.D, true));
		System.out.println(tracker.hasChanged(TestEnum.D, true));
		System.out.println(tracker.hasChanged(null, true));
		System.out.println(tracker.hasChanged(null, true));
		System.out.println(tracker.hasChanged(TestEnum.B, true));
		
	}
	
	// ========================================================================
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		TestThisClass();
	}

}
