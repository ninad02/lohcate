package shared;

import java.awt.image.BufferedImage;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.net.MalformedURLException;
import java.net.URL;
import java.util.HashMap;

import javax.imageio.ImageIO;

/**
 * @author Siddharth G. Reddy
 *
 */
public class FileOps {
	
	public static final String rootDir = "";
	
	public static BufferedImage loadImage(String fileName) {
		BufferedImage img = null;
		try {
		    img = ImageIO.read(new File(fileName));
		} catch (IOException e) { }
		return img;
	}
	
	public static BufferedImage getImage(String url) {
		URL page;
		try {
			page = new URL(url);
		} catch (MalformedURLException e1) {
			e1.printStackTrace();
		}
		BufferedImage img = null;
		try {
			page = new URL(url);
			img = ImageIO.read(page);
		} catch (Exception e) { e.printStackTrace(); }
		return img;
	}
	
	public static void writeImageToFile(BufferedImage img, String fileName) {
		try {
		    File f = new File(fileName);
		    File fi = new File(fileName.replace(fileName.split("/")[fileName.split("/").length-1], ""));
			if (!fi.isDirectory())
				fi.mkdirs();
		    ImageIO.write(img, "png", f);
		} catch (IOException e) { e.printStackTrace(); }
	}
	
	public static String getHTML(String url) throws IOException {
		URL page = new URL(url);
		BufferedReader in = new BufferedReader(new InputStreamReader(page.openStream()));
		String contents = "", inputLine = "";
		while ((inputLine = in.readLine()) != null)
			contents += inputLine + "\n";
		in.close();
		return contents;
	}
	
	public static String loadFromFile(String fileName) {
		File file = new File(fileName);
        StringBuffer contents = new StringBuffer();
        BufferedReader reader = null;
        

        try {
            reader = new BufferedReader(new FileReader(file));
            String text = null;

            while ((text = reader.readLine()) != null) {
                contents.append(text)
                    .append(System.getProperty(
                        "line.separator"));
            }
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        } finally {
            try {
                if (reader != null) {
                    reader.close();
                }
            } catch (IOException e) {
                e.printStackTrace();
            }
        }
		
		return contents.toString();
	}
	
	public static String t_loadFromFile(String fileName) throws Exception {
		File file = new File(fileName);
        StringBuffer contents = new StringBuffer();
        BufferedReader reader = null;

        try {
            reader = new BufferedReader(new FileReader(file));
            String text = null;

            while ((text = reader.readLine()) != null) {
                contents.append(text)
                    .append(System.getProperty(
                        "line.separator"));
            }
        } finally {
        	if (reader != null)
        		reader.close();
        }
        
		return contents.toString();
	}
	
	
	public static void writeToFile(String filename, String toWrite) {
		writeToFile(filename, toWrite, false, false);
	}
	
	
	
	public static void writeToFile(String filename, String toWrite, boolean writeNewLine, boolean appendTextToFile) {
		File f = new File(filename);		
		
		File parentDir = new File(f.getParent());
		if (!parentDir.isDirectory())
			parentDir.mkdirs();
			    
	    try {	    	
	    	BufferedWriter out = new BufferedWriter(new FileWriter(f, appendTextToFile));
	    	out.write(toWrite);
	    	if (writeNewLine) {
	    		out.newLine();
	    	}
	    	out.close();
	    	
	    } catch (Exception e) {
	    	e.printStackTrace();
	    	System.exit(-1);
	    }
	}
	
	public static void t_writeToFile(String filename, String toWrite) throws Exception {
		File f = new File(filename);
		File fi = new File(filename.replace(filename.split("/")[filename.split("/").length-1], ""));
		if (!fi.isDirectory())
			fi.mkdirs();
	    FileOutputStream fop = null;
	    fop = new FileOutputStream(f);	    
	    fop.write(toWrite.getBytes());
		fop.close();
	}
	

	public static void appendToFile(String fileName, String toWrite) {
		try { 
			BufferedWriter out = new BufferedWriter(new FileWriter(fileName, true)); 
			out.write(toWrite); 
			out.close(); 
		} catch (IOException e) { 
			try {
				writeToFile(fileName, toWrite);
			} catch (Exception e_2) { e.printStackTrace(); }
		} 
	}

}
