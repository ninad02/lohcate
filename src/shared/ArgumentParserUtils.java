package shared;

import java.util.Collection;

import com.martiansoftware.jsap.JSAP;
import com.martiansoftware.jsap.JSAPException;
import com.martiansoftware.jsap.JSAPResult;
import com.martiansoftware.jsap.Parameter;
import com.martiansoftware.jsap.Switch;

public class ArgumentParserUtils {

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
