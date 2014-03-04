package nutils.image;

import java.awt.GridBagLayout;
import java.awt.Point;
import java.awt.geom.Point2D;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;

import javax.imageio.ImageIO;

import nutils.ArrayUtils;
import nutils.CompareUtils;
import nutils.IOUtils;

/**
 * A class for appending images together
 * @author Ninad Dewal
 *
 */
public class ImageAppender {

	private static final int IndexNumRows = 0;
	private static final int IndexNumCols = 1;
	private static final int IndexDirection = 2;
	private static final int IndexFirstImageFilename = IndexDirection + 1;
	
	// ========================================================================
	public static void appendImages(int numRows, int numCols, boolean directionDownFirst, String outFilename, String[] inFilenames) {
		
		ArrayList<BufferedImage> imagesIn = new ArrayList<BufferedImage>(inFilenames.length);		
		
		// Read all the images		
		for (String inFilename : inFilenames) {			
			try {
				BufferedImage imageIn = ImageIO.read(new File(inFilename));
				imagesIn.add(imageIn);
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
				System.exit(-1);
			}
			
		}
		
		// Ensure that the number of cells exceeds the number of images
		numRows = Math.max(numRows, 1);
		numCols = Math.max(numCols, 1);
		
		if (numRows == 1) {
			numCols = Math.max(numCols, imagesIn.size());
		} else if (numCols == 1) {
			numRows = Math.max(numRows, imagesIn.size());
		}
		
		// Ensure enough cells
		while ((numRows * numCols) < imagesIn.size()) {
			if (directionDownFirst) {
				++numRows;
			} else {
				++numCols;
			}
		}
			
		// Create matrix
		BufferedImage[][] images = new BufferedImage[numRows][numCols];
		ArrayUtils.arrayFill(images, null);
		
		// Now distribute the images
		int indexRow = 0;
		int indexCol = 0;
		for (BufferedImage theImage : imagesIn) {
			images[indexRow][indexCol] = theImage;
			if (directionDownFirst) {		
				indexRow++;
				if (indexRow >= numRows) {
					// Move to the next col
					indexRow = 0;
					indexCol++;
				}
			} else {
				indexCol++;
				if (indexCol >= numCols) {
					indexCol = 0;
					indexRow++;
				}
			}
		}
						
		// Get total max height and width
		Point widthHeight = calcTotalMaxWidthHeight(images, null);		
		
		// Now create new dimensions for image
		BufferedImage img = new BufferedImage(widthHeight.x, widthHeight.y, BufferedImage.TYPE_INT_RGB);

		int numImagesDrawn = 0;
		Point coordinate = new Point(0, 0);		
		
		// Now draw from left to right, since we already know the bounds and ordering
		for (int rowIndex = 0; rowIndex < images.length; rowIndex++) {
			BufferedImage[] allImagesInRow = images[rowIndex];			
			for (int colIndex = 0; colIndex < allImagesInRow.length; colIndex++) {
				BufferedImage theImage = allImagesInRow[colIndex];
				if (theImage != null) {
					boolean imageDrawn = img.createGraphics().drawImage(theImage, coordinate.x, coordinate.y, null);
					if (!imageDrawn) {
						CompareUtils.ensureTrue(false, "ERROR: Image could not be appended!");
					}
					++numImagesDrawn;
					
					// Move the coordinate over horizontally
					coordinate.x += theImage.getWidth();					
				}								
			}
			coordinate.x = 0;
			coordinate.y += allImagesInRow[0].getHeight();			
		}
		
		// horizontally
		File outFile = new File(outFilename);		  		  		
		try {
			boolean drewImageFile = ImageIO.write(img, "png", outFile);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}		
	}

	// ========================================================================
	private static Point calcTotalMaxWidthHeight(BufferedImage[][] images, Point widthHeight) {		
		
		// create a new return value if not already passed in
		widthHeight = (widthHeight == null) ? new Point(0, 0) : widthHeight;
		
		int[] heightsInCols = ArrayUtils.newIntArray(images[0].length, 0);  // first row will always have most # of filled cols
		int maxWidthTotal = Integer.MIN_VALUE;
		
		for (BufferedImage[] allImagesInRow : images) {
			int widthTotal = 0;
			
			for (int colIndex = 0; colIndex < allImagesInRow.length; colIndex++) {
				BufferedImage imageInRow = allImagesInRow[colIndex]; 
				
				if (imageInRow != null) {
					heightsInCols[colIndex] += imageInRow.getHeight(); 
					widthTotal              += imageInRow.getWidth();
				}
			}
			
			// Now that we finished a row, figure out the total max width
			maxWidthTotal = Math.max(maxWidthTotal, widthTotal);
		}
		
		int indexMaxElement = ArrayUtils.getIndexOfMaxElement(heightsInCols, 0);
		CompareUtils.ensureTrue(indexMaxElement >= 0, "ERROR: No Columns present!");
		widthHeight.y = heightsInCols[indexMaxElement];
		widthHeight.x = maxWidthTotal;
		return widthHeight;
	}
		
	// ========================================================================
	public static void main(String[] args) {
		//ImageAppender.appendImages(args);
		String[] inFilenames = new String[] {
				IOUtils.pathConcat(new String[] { "E:", "Research", "Data", "tcga-acc", "plots", "copyNumber",     "TCGA-OR-A5J2-01A-11D-A29I-10.01A-vs-10A.CopyNumber_GenomeWide.png"}),
				IOUtils.pathConcat(new String[] { "E:", "Research", "Data", "tcga-acc", "plots", "VAF_Genomewide", "TCGA-OR-A5J2-01A-11D-A29I-10.01A-vs-10A.VAF_GenomeWide_Normal.png"}),
				IOUtils.pathConcat(new String[] { "E:", "Research", "Data", "tcga-acc", "plots", "VAF_Genomewide", "TCGA-OR-A5J2-01A-11D-A29I-10.01A-vs-10A.VAF_GenomeWide_Tumor.png"}),
		};
		String outFilename = IOUtils.pathConcat(new String[] { "E:", "Research", "Data", "tcga-acc", "plots", "VAF_Genomewide", "TCGA-OR-A5J2-01A-11D-A29I-10.01A-vs-10A.VAF_GenomeWide_Tumor.CONCAT.png"});
		ImageAppender.appendImages(3, 1, true, outFilename, inFilenames);
	}
}
