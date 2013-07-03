package nutils;

import java.util.Collection;

import com.martiansoftware.jsap.JSAP;
import com.martiansoftware.jsap.JSAPException;
import com.martiansoftware.jsap.JSAPResult;
import com.martiansoftware.jsap.Parameter;
import com.martiansoftware.jsap.Switch;

public class ArgumentParserUtils {

	public static final char  NoShortFlag = JSAP.NO_SHORTFLAG;
	public static final String NoLongFlag = JSAP.NO_LONGFLAG;
	public static final String NoUsageName = "";
	
	// ========================================================================
	public static abstract class InputParameter<T> {
		protected char   mFlagShort;
		protected String mFlagLong;
		protected String mUsageName;
		protected String mName;
		
		protected T      mDefaultValue;
		protected T      mActualValue;
		
		public InputParameter(T defaultValue, String name, char flagShort, String flagLong, String usageName) {
			mFlagShort = flagShort;
			mFlagLong  = flagLong;
			mUsageName = usageName;
			mName      = name;
			mActualValue = mDefaultValue = defaultValue;
		}
		
		public char getShortFlag() { return mFlagShort; }
		public String getLongFlag() { return mFlagLong; }
		public String getUsageName() { return mUsageName; }
		public String getName() { return mName; }
		
		public T getValue()        { return mActualValue; } 
		public T getDefaultValue() { return mDefaultValue; }
		
		public void setValue(T newValue) { mActualValue = newValue; } 
		
		public String getNameAndValueAsString(String delimiter) { 
			return (mFlagLong + delimiter + getValue());
		}
		
		public abstract void parseValue(String newValueStr);
	}

	// ========================================================================
	public static abstract class InputParameterNumber<T extends Number> extends InputParameter<T> {
		
		public InputParameterNumber(T defaultValue, String name, char flagShort, String flagLong, String usageName) {
			super(defaultValue, name, flagShort, flagLong, usageName);
		}
	}
	
	// ========================================================================
	public static class InputParameterInteger extends InputParameterNumber<Integer> {
		
		public InputParameterInteger(Integer defaultValue, String name, char flagShort, String flagLong, String usageName) {
			super(defaultValue, name, flagShort, flagLong, usageName);
		}
		
		public void parseValue(String newValueStr) {
			mActualValue = Integer.parseInt(newValueStr);
		}		
	}
	
	// ========================================================================
	public static class InputParameterDouble extends InputParameterNumber<Double> {
		
		public InputParameterDouble(Double defaultValue, String name, char flagShort, String flagLong, String usageName) {
			super(defaultValue, name, flagShort, flagLong, usageName);
		}
		
		public void parseValue(String newValueStr) {
			mActualValue = Double.parseDouble(newValueStr);
		}
	}
	
	// ========================================================================
	public static class InputParameterBoolean extends InputParameter<Boolean> {
		
		public InputParameterBoolean(Boolean defaultValue, String name, char flagShort, String flagLong, String usageName) {
			super(defaultValue, name, flagShort, flagLong, usageName);
		}
		
		public InputParameterBoolean(Boolean defaultValue, String name) {
			super(defaultValue, name, NoShortFlag, NoLongFlag, NoUsageName);
		}
		
		public void parseValue(String newValueStr) {
			mActualValue = Boolean.parseBoolean(newValueStr);
		}
	}
	
	// ========================================================================
	public static class InputParameterString extends InputParameter<String> {

		public InputParameterString(String defaultValue, String name, char flagShort, String flagLong, String usageName) {
			super(defaultValue, name, flagShort, flagLong, usageName);
		}
		
		public void parseValue(String newValueStr) {
			mActualValue = newValueStr;
		}
	}
	
	// ========================================================================
	public static void registerJSAPParameters(JSAP jsap, Collection<Parameter> params) {
		for (Parameter p : params) { registerJSAPParameter(jsap, p); }
	}
	
	// ========================================================================
	public static void registerJSAPParameters(JSAP jsap, Parameter[] params) {
		for (Parameter p : params) { registerJSAPParameter(jsap, p); }
	}
	
	// ========================================================================
	public static void registerJSAPParameter(JSAP jsap, Parameter p) {
		try { 
			jsap.registerParameter(p); 
		} catch (JSAPException e) { 
			e.printStackTrace(); 
			System.exit(-1); 
		}	
	}

	// ========================================================================
	/** Creates a switch and registers it with the JSAP object.  Returns the switch. */
	public static<T> Switch createSwitch(InputParameter<T> param, JSAP jsap) {
		return createSwitch(param.getName(), param.getShortFlag(), param.getLongFlag(), param.getUsageName(), jsap);
	}
	
	// ========================================================================
	/** Creates a switch and registers it with the JSAP object. Returns the switch. 
	 *  @param name The name of the switch
	 *  @param shortFlag The short flag of the switch
	 *  @param longFlag The long flag of the switch.  If null, then sets it to the name. 
	 *  @param the JSAP object
	 *  @return the switch. */
	public static Switch createSwitch(String name, char shortFlag, String longFlag, String help, JSAP jsap) {
		longFlag = (longFlag == null) ? name : longFlag; 		
		help     = (help == null) ? JSAP.NO_HELP : help;
		Switch sw = new Switch(name, shortFlag, longFlag, help);		
		registerJSAPParameter(jsap, sw);		
		return sw;
	}

	// ========================================================================
	/** Creates a flagged option and registers it with the JSAP object.  Returns the flagged option
	 *  @param name The name of the flagged option
	 *  @param valueParser The parser used to parse the value of the object
	 *  @param defaultValue
	 */
	
	// ========================================================================
	/** Parses the argument list with the JSAP object and returns the result.  If there is a failure,
	 *  this prints a failure message and exits.
	 */
	public static JSAPResult parseAndCheck(String[] args, JSAP jsap, String className) {
		JSAPResult config = jsap.parse(args);
		printJSAPParseFail(config, jsap, className);
		return config;
	}
		
	// ========================================================================
	/* The source code for this method is taken from the JSAP website examples. */
	public static void printJSAPParseFail(JSAPResult config, JSAP jsap, String className) {
		if (!config.success()) {
            
            System.err.println();

            // print out specific error messages describing the problems
            // with the command line, THEN print usage, THEN print full
            // help.  This is called "beating the user with a clue stick."
            for (java.util.Iterator errs = config.getErrorMessageIterator();
                    errs.hasNext();) {
                System.err.println("Error: " + errs.next());
            }
            
            System.err.println();
            System.err.println("Usage: java " + className);
            System.err.println("                " + jsap.getUsage());
            System.err.println();
            System.err.println(jsap.getHelp());
            System.exit(1);
        }
	}
}
