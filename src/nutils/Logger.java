package nutils;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.OutputStream;
import java.io.PrintStream;
import java.io.UnsupportedEncodingException;

public class Logger extends PrintStream {

	public Logger(File arg0) throws FileNotFoundException {
		super(arg0);
		// TODO Auto-generated constructor stub
	}

	public Logger(File file, String csn) throws FileNotFoundException,
			UnsupportedEncodingException {
		super(file, csn);
		// TODO Auto-generated constructor stub
	}

	public Logger(OutputStream out, boolean autoFlush, String encoding)
			throws UnsupportedEncodingException {
		super(out, autoFlush, encoding);
		// TODO Auto-generated constructor stub
	}

	public Logger(OutputStream out, boolean autoFlush) {
		super(out, autoFlush);
		// TODO Auto-generated constructor stub
	}

	public Logger(OutputStream out) {
		super(out);
		// TODO Auto-generated constructor stub
	}

	public Logger(String fileName, String csn) throws FileNotFoundException,
			UnsupportedEncodingException {
		super(fileName, csn);
		// TODO Auto-generated constructor stub
	}

	public Logger(String fileName) throws FileNotFoundException {
		super(fileName);
		// TODO Auto-generated constructor stub
	}

	public PrintStream printf(String format, Object ... args) {
		return super.format(format, args);
	}
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub		
		
	}

}
