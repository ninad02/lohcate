package shared;

import genomeEnums.Chrom;
import genomeUtils.RegionRange;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.concurrent.ArrayBlockingQueue;

import org.apache.commons.math3.analysis.interpolation.SplineInterpolator;
import org.apache.commons.math3.analysis.polynomials.PolynomialSplineFunction;
import org.bcm.hgsc.cancer.bed.BedGraph;

import com.carrotsearch.hppc.DoubleArrayList;
import com.google.common.io.Files;

import nutils.ArrayUtils;
import nutils.EnumMapSafe;
import nutils.IOUtils;
import nutils.NumberUtils;
import nutils.StringUtils;

public class OriginReplicationBEDReaderTester {

	private static final int SplineIncrement = 2;
	
	public static void doWork(File inFile, String outDir, int windowSize) {
		
		IOUtils.createDirectoryPath(outDir, false);
		
		StringBuilder sb = new StringBuilder(1024);
		ArrayList<String> allLines = IOUtils.readAllLinesFromFile(inFile.getAbsolutePath(), true, true, sb);
		String header = sb.toString();
		
		RegionRange range = new RegionRange();
		EnumMapSafe<Chrom, ArrayList<MutBEDLine>> infoPerChrom = EnumMapSafe.createEnumMapOfArrayLists(Chrom.class, MutBEDLine.class);
						
		for (int row = 0; row < allLines.size(); row++) {			
			String[] cols = allLines.get(row).split(StringUtils.TabPatternStr);
			
			MutBEDLine mbl = new MutBEDLine();
			mbl.mChrom = Chrom.getChrom(cols[0]);
			mbl.mPosStart = Integer.parseInt(cols[1]);
			mbl.mPosEnd   = Integer.parseInt(cols[2]);			
			int mutValue = Integer.parseInt(cols[4]);
			mbl.mMutType = (mutValue > 0) ? MutPOLEType.TCTtoTAT : MutPOLEType.AGAtoATA; 			
			
			infoPerChrom.get(mbl.mChrom).add(mbl);			
		}
		
		// Clear the strings
		allLines.clear();

		// The blocking queue
		ArrayBlockingQueue<MutBEDLine> abq = new ArrayBlockingQueue<MutBEDLine>(windowSize);
				
		String outFilenameRoot = outDir + File.separator + inFile.getName();
		
		// The writers
		BedGraph bgPvalue = new BedGraph(inFile.getName() + " P-Val", "ORI Window P-values for " + inFile.getName());		
		BedGraph bgScaled = new BedGraph(inFile.getName() + " -log_10(P-Val)", "ORI Window -log_10(P-values) for " + inFile.getName());	
		BedGraph bgORI    = new BedGraph(inFile.getName() + " ORI/Collision", "ORI(-1)/Collision(1) Locations for " + inFile.getName());		
		BedGraph bgSpline = new BedGraph(inFile.getName() + " Spline", "ORI Window Splined (" + SplineIncrement + ") P-values for " + inFile.getName());
		
		// Now go thorugh all chroms		
		for (Chrom chrom : Chrom.values()) {
			
			if (chrom.isInvalid()) continue;
			
			ArrayList<MutBEDLine> bedItemList = infoPerChrom.get(chrom);
			ArrayList<BedGraphEntry> bedGraphEntries = new ArrayList<>(bedItemList.size());
			
			// Clear the blocking queue since on new chromosome
			abq.clear();
			int[] countMutType = new int[MutPOLEType.values().length];
			Arrays.fill(countMutType, 0);
			
			boolean pValPrevSet = false;
			double pValPrev = Double.NaN;
			
			int lastIndex = bedItemList.size() - 1;			
			for (int row = 0; row <= lastIndex; row++) {
				
				MutBEDLine bedItem = bedItemList.get(row);
				
				if (abq.remainingCapacity() == 0) {
					// Remove the head element and decrease the count
					MutBEDLine headElement = null;
					try {
						headElement = abq.take();
					} catch (InterruptedException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
					--countMutType[headElement.mMutType.ordinal()];
				}
				
				boolean added = abq.offer(bedItem);				
				if (!added) {
					System.err.println("Impossible state!");
					System.exit(-1);
				}
				
				++countMutType[bedItem.mMutType.ordinal()];
				
				// Now we see if the queue is full again or if we've reached the last element.
				// If so, we calculate the statistical test
				if ((abq.remainingCapacity() == 0) || (row == lastIndex)) {
					int totalCount = ArrayUtils.arraySum(countMutType);
					assert((row < lastIndex && totalCount == windowSize) || (row == lastIndex));
					int indexMaxElement = ArrayUtils.getIndexOfMaxElement(countMutType, 0);							
					int multiplier = (MutPOLEType.values()[indexMaxElement] == MutPOLEType.AGAtoATA) ? -1 : 1;
					
					int maxCount = countMutType[indexMaxElement];
					double pVal = NumberUtils.cumulativeProbabilitySuccessBinomial(totalCount, maxCount, 0.5) ;
					double pValMultiplied = multiplier * pVal;
					
					double scaledPVal = -NumberUtils.MathLog10Safe(pVal);					
					int midPoint = getMidPoint(abq);					
					//System.out.println(chrom.ordinal() + "\t" + midPoint + "\t" + scaledPVal);
					BedGraphEntry bge = new BedGraphEntry(chrom, midPoint, pValMultiplied);
					bedGraphEntries.add(bge);					
					
					if (pValPrevSet) {
						Boolean resultORI = isORI(pValPrev, pValMultiplied);
						if (resultORI != null) {
							double value = resultORI.booleanValue() ? -1 : 1; 								
							try {
								bgORI.add("chr" + chrom.getName(), midPoint, midPoint, value);
							} catch (IOException e) {
								// TODO Auto-generated catch block
								e.printStackTrace();
							}
						}
					}
					pValPrev = pValMultiplied;
					pValPrevSet = true;

					
					try {
						bgPvalue.add("chr" + chrom.getName(), midPoint, midPoint, pValMultiplied);
						bgScaled.add("chr" + chrom.getName(), midPoint, midPoint, pValMultiplied);						
					} catch (IOException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}					
				}
			}
						
			int incr = SplineIncrement;			
			
			System.out.println(chrom);
			DoubleArrayList splineX = new DoubleArrayList(bedGraphEntries.size());
			DoubleArrayList splineY = new DoubleArrayList(bedGraphEntries.size());
			int indexAdded = 0;
			for (int i = 0; i < bedGraphEntries.size(); i += incr) {
				BedGraphEntry bge = bedGraphEntries.get(i);
				splineX.add(bge.mPos);
				splineY.add(bge.mValue);
				indexAdded = i;					
			}				
			if (indexAdded < (bedGraphEntries.size() - 1)) {
				splineX.add(bedGraphEntries.get(bedGraphEntries.size() - 1).mPos);
				splineY.add(bedGraphEntries.get(bedGraphEntries.size() - 1).mValue);
			}
			//System.out.println(splineX);

			if (splineX.size() > 10) {
				//Spline spline = new Spline(splineX, splineY);
				SplineInterpolator splineInterpolator = new SplineInterpolator();
				PolynomialSplineFunction psf = splineInterpolator.interpolate(splineX.toArray(), splineY.toArray());
				
				for (BedGraphEntry bge : bedGraphEntries) {
					try {
						bgSpline.add("chr" + chrom.getName(), bge.mPos, bge.mPos, psf.value(bge.mPos));
					} catch (IOException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
				}
			}

		}
		
		String bedGraphPval           = outFilenameRoot + ".ORI.windowSize." + windowSize + ".pval.txt";
		String bedGraphPvalMinusLog10 = outFilenameRoot + ".ORI.windowSize." + windowSize + ".pval.minusLog10.txt";
		String bedGraphSpline         = outFilenameRoot + ".ORI.windowSize." + windowSize + ".pval.spline." + SplineIncrement + ".txt";
		String bedGraphORI            = outFilenameRoot + ".ORI.windowSize." + windowSize + ".pval.ORI_Collisions.txt";
		
		try {
			bgPvalue.write(new File(bedGraphPval));
			bgPvalue.close();
			
			bgScaled.write(new File(bedGraphPvalMinusLog10));
			bgScaled.close();
			
			bgORI.write(new File(bedGraphORI));
			bgORI.close();
			
			bgSpline.write(new File(bedGraphSpline));
			bgSpline.close();
			
		} catch (IOException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}
		
		
		//track name="ColorByStrandDemo" description="Color by strand demonstration" visibility=2 colorByStrand="255,0,0 0,0,255"
		
		
	}
	
	public static void convertHeader() {
		
	}

	// ===========================
	public static Boolean isORI(double pValScaledPrev, double pValScaledCurr) {
		if (pValScaledPrev >= 0 && pValScaledCurr >= 0) {
			return null;
		} else if (pValScaledPrev < 0 && pValScaledCurr < 0) {
			return null;
		} else if (pValScaledPrev > 0 && pValScaledCurr < 0) {
			// We have a collision 
			return Boolean.FALSE;
		} else if (pValScaledPrev < 0 && pValScaledCurr > 0) {
			return Boolean.TRUE;
		} else {
			System.out.println("Huha!");
			return null;
		}
	}
	
	// ===========================
	public static int getMidPoint(ArrayBlockingQueue<MutBEDLine> abq) {
		long sum = 0;
		
		for (Iterator<MutBEDLine> iter = abq.iterator(); iter.hasNext(); ) {
			sum += iter.next().mPosStart;
		}
		
		return (int) Math.round((double) sum / ((double) abq.size())); 
	}


	
	// ===========================
	private static enum MutPOLEType {
		TCTtoTAT,
		AGAtoATA;
	}

	// ===========================
	public static class BedGraphEntry {
		Chrom mChrom;
		int mPos;
		double mValue;

		public BedGraphEntry(Chrom chrom, int position, double value) {
			mChrom = chrom;
			mPos = position;
			mValue = value;
		}
	}
	
	// ===========================
	public static class MutBEDLine {
		Chrom mChrom;
		int mPosStart, mPosEnd;
		MutPOLEType mMutType;		
	}

	// ===========================
	public static void workOnDirectory(String inDir, String outDir, int[] windowSizes) {
		File[] files = (new File(inDir)).listFiles();
		ArrayList<File> validFiles = new ArrayList<>();
		
		for (File file : files) {
			if (file.getName().endsWith(".BED.TXT")) {
				validFiles.add(file);
			}
		}
		
		IOUtils.createDirectoryPath(outDir, false);		
		
		for (File file : validFiles) {
			for (int windowSize : windowSizes) {
				doWork(file, outDir, windowSize);
			}
		}
	}
	
	// ===========================
	public static void main(String[] args) {
		int[] windowSizes = new int[] { 11, 21, 31, 41 };
		workOnDirectory(args[0], args[1], windowSizes);
	}
}



