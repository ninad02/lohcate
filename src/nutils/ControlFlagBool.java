package nutils;

/** 
 * A convenience class to place control flags inside a program.  The intent is for
 * the developer to use this class instead of a boolean variable; the benefit is that
 * it's easier for an IDE to track this class rather than a boolean primitive.  This
 * way, the developer can track all the internal control flags in the program.
 * 
 * @author Ninad Dewal
 *
 */
public class ControlFlagBool {

	protected Boolean mFlag;
	
	// ========================================================================
	public ControlFlagBool(boolean initialValue) {
		setValue(initialValue);
	}

	// =======================================================================
	public boolean getValue() { return mFlag.booleanValue(); }
	
	// ========================================================================
	public void setValue(boolean value) {
		mFlag = Boolean.valueOf(value);
	}

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub

	}

}
