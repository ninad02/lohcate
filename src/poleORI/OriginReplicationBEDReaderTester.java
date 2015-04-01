package poleORI;

import genomeEnums.Chrom;
import genomeUtils.RegionRange;
import genomeUtils.RegionRangeWithPayload;
import htsjdk.tribble.util.MathUtils;

import java.awt.Color;
import java.awt.Font;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.Iterator;
import java.util.concurrent.ArrayBlockingQueue;

import lohcate.clustering.ClusteringPlotting;
import lohcateEnums.ColorPastel;
import lohcateEnums.EventType;

import org.apache.commons.math3.analysis.interpolation.SplineInterpolator;
import org.apache.commons.math3.analysis.polynomials.PolynomialSplineFunction;
import org.bcm.hgsc.cancer.bed.BedGraph;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.LogAxis;
import org.jfree.chart.axis.LogarithmicAxis;
import org.jfree.chart.axis.NumberTickUnit;
import org.jfree.chart.axis.TickUnit;
import org.jfree.chart.axis.TickUnits;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYItemRenderer;
import org.jfree.chart.renderer.xy.XYSplineRenderer;
import org.jfree.chart.title.LegendTitle;
import org.jfree.data.xy.DefaultXYDataset;
import org.jfree.data.xy.XYSeries;

import com.carrotsearch.hppc.FloatArrayList;
import com.carrotsearch.hppc.IntArrayList;
import com.carrotsearch.hppc.DoubleArrayList;
import com.google.common.io.Files;

import nutils.ArrayUtils;
import nutils.ArrayUtils.ParallelArrayDouble;
import nutils.Cast;
import nutils.CloneInf;
import nutils.CompareUtils;
import nutils.EnumMapSafe;
import nutils.GraphUtils;
import nutils.IOUtils;
import nutils.NullaryClassFactory;
import nutils.NumberUtils;
import nutils.PrimitiveWrapper;
import nutils.StringUtils;
import nutils.BitUtils.BitSetUtils;
import nutils.counter.DynamicBucketCounter;

public class OriginReplicationBEDReaderTester {

	private static final int SplineIncrement = 3;
	
	// ===========================
	public static enum MutPOLEType implements CloneInf<MutPOLEType> {
		TCTtoTAT,
		AGAtoATA;

		@Override
		public MutPOLEType makeClone() { return this; }

		@Override
		public MutPOLEType makeClone(boolean deepCopy) { return this; }
	}
	
	// ========================================================================
	public static EnumMapSafe<Chrom, ArrayList<MutBEDLine>> readMutBEDFile(File inFile, StringBuilder headerStringBuffer) {
		ArrayList<String> allLines = IOUtils.readAllLinesFromFile(inFile.getAbsolutePath(), true, true, headerStringBuffer);
		return readMutBEDAllLines(allLines);
	}
	
	// ========================================================================
	public static EnumMapSafe<Chrom, ArrayList<MutBEDLine>> readMutBEDAllLines(ArrayList<String> allLines) {
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
		
		return infoPerChrom;
	}

	// ========================================================================
	public static void doWork(File inFile, String outDir, String transcriptStartSiteFile, int windowSize,
			DefaultXYDataset xyDatasetCumCount, DefaultXYDataset xyDatasetCumFraction,
			DefaultXYDataset xyDatasetCount, DefaultXYDataset xyDatasetFraction) {
		
		IOUtils.createDirectoryPath(outDir, false);
		
		String sampleName = inFile.getName().substring(0, inFile.getName().indexOf("."));
		
		// Read in the data from the file
		StringBuilder sb = new StringBuilder(1024);
		EnumMapSafe<Chrom, ArrayList<MutBEDLine>> infoPerChrom = readMutBEDFile(inFile, sb);
		String header = sb.toString();
			
		// The blocking queue for the sliding window
		ArrayBlockingQueue<MutBEDLine> abq = new ArrayBlockingQueue<MutBEDLine>(windowSize);
		
		// The counter for the inter-ORI distances
		DynamicBucketCounter oriDistCount = new DynamicBucketCounter();		
		DynamicBucketCounter oriDistCountRaw = new DynamicBucketCounter();
		DynamicBucketCounter oriDistCountPermutation = new DynamicBucketCounter();
		
		// Map for keeping all the windows over the chromosome		
		EnumMapSafe<Chrom, ArrayList<BinTestWindow>> binTestWindowsOnePatient = EnumMapSafe.createEnumMapOfArrayLists(Chrom.class, BinTestWindow.class);
		EnumMapSafe<Chrom, BitSet> oriIndicesOnePatient = new EnumMapSafe<Chrom, BitSet>(Chrom.class); 
		
		InterORIDistancePermutationTask interORIpermTask = new InterORIDistancePermutationTask(null);
		
		String outFilenameRoot = outDir + File.separator + inFile.getName();
		
		// The writers
		BedGraph bgPvalue = new BedGraph(inFile.getName() + " P-Val", "ORI Window P-values for " + inFile.getName());		
		BedGraph bgScaled = new BedGraph(inFile.getName() + " -log_10(P-Val)", "ORI Window -log_10(P-values) for " + inFile.getName());	
		//BedGraph bgORI    = new BedGraph(inFile.getName() + " ORI/Collision", "ORI(1)/Collision(-1) Locations for " + inFile.getName());		
		BedGraph bgORI    = new BedGraph(inFile.getName() + " ORI", "ORI Locations for " + inFile.getName());
		BedGraph bgSpline = new BedGraph(inFile.getName() + " Spline", "ORI Window Splined (" + SplineIncrement + ") P-values for " + inFile.getName());
		
		// Now go thorugh all chroms		
		for (Chrom chrom : Chrom.values()) {
			
			if (chrom.isInvalid()) continue;
			
			ArrayList<MutBEDLine> bedItemList = infoPerChrom.get(chrom);
			ArrayList<BinTestWindow> binTestWindows = binTestWindowsOnePatient.get(chrom);  // corresponds to the windows
			ArrayList<MutBEDLine> midPoints = new ArrayList<MutBEDLine>(); 
			
			// Prepare the spline points
			BitSet binTestWindowIsORI       = new BitSet(bedItemList.size());
			BitSet binTestWindowIsCollision = new BitSet(bedItemList.size());
			oriIndicesOnePatient.put(chrom, binTestWindowIsORI);
			
			// Clear the blocking queue since on new chromosome
			abq.clear();
			
			// Make an arraylist of ORI indices
			IntArrayList oriIndices = new IntArrayList();
			
			// A buffer to store mutation type counts 
			int[] countMutType = ArrayUtils.newIntArray(MutPOLEType.values().length, 0);
						
			// A holder for the previous window's p-value, with capabilities to test if the value was set or not
			PrimitiveWrapper.WDoubleSetOnce pValPrev = new PrimitiveWrapper.WDoubleSetOnce(Double.NaN);
			
			// The final index of the list of items
			int lastIndex = bedItemList.size() - 1;
			
			// Create a map of the first indices of a mutation type
			EnumMapSafe<MutPOLEType, PrimitiveWrapper.WIntegerSetOnce> firstIndices = PrimitiveWrapper.WIntegerSetOnce.ClassFactory.newEnumMap(MutPOLEType.class);

			// Now, loop through the rows
			for (int row = 0; row <= lastIndex; row++) {
				
				MutBEDLine bedItem = bedItemList.get(row);
				
				// Now remove the head of the queue and adjust the count
				if (abq.remainingCapacity() == 0) {
					// Remove the head element and decrease the count					
					try {
						MutBEDLine headElement = abq.take();
						--countMutType[headElement.mMutType.ordinal()];
					} catch (InterruptedException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}					
				}
				
				// Now add to the tail of the queue and adjust the count
				boolean added = abq.offer(bedItem);				
				if (!added) {
					System.err.println("Impossible state!");
					System.exit(-1);
				}				
				++countMutType[bedItem.mMutType.ordinal()];
				boolean firstPointWasSet = firstIndices.get(bedItem.mMutType).set(row);
				
				// Now we see if the queue is full again or if we've reached the last element.
				// If so, we calculate the statistical test
				if ((abq.remainingCapacity() == 0) || (row == lastIndex)) {
					int totalCount = ArrayUtils.arraySum(countMutType);
					assert((row < lastIndex && totalCount == windowSize) || (row == lastIndex));
					int indexMaxElement = ArrayUtils.getIndexOfMaxElement(countMutType, 0);							
					
					int maxCount = countMutType[indexMaxElement];					
					double pVal = NumberUtils.cumulativeProbabilitySuccessBinomial(totalCount, maxCount, 0.5) ;
					int multiplier = getMultiplierPvalue(MutPOLEType.values()[indexMaxElement]);
					double pValMultiplied = multiplier * pVal;
					
					double scaledPVal = -NumberUtils.MathLog10Safe(pVal);					
					int midPoint = getMidPoint(abq, midPoints);					
					//System.out.println(chrom.ordinal() + "\t" + midPoint + "\t" + scaledPVal);
					BinTestWindow windowCurrent = new BinTestWindow(chrom, midPoints.get(0).mPosStart, midPoints.get(1).mPosStart, pValMultiplied);
					binTestWindows.add(windowCurrent);					
					
					if (pValPrev.isAlreadySetOnce()) {
						PvalueTrackState trackState = isORI(pValPrev.mDouble, pValMultiplied);
						
						if (trackState.equals(PvalueTrackState.Positive)) {
							// We are not in a easy-to-detect ORI region
							if (firstIndices.get(MutPOLEType.AGAtoATA).isAlreadySetOnce()) {
								//System.out.println("First index AGA: " + firstIndices.get(MutPOLEType.AGAtoATA).mInt + "\t" + row);
								int diff = row - firstIndices.get(MutPOLEType.AGAtoATA).mInt + 1;
								if (diff == windowSize) {
									int binTestWindowEarlierIndex = binTestWindows.size() - (windowSize >>> 1);
									binTestWindowEarlierIndex = Math.max(binTestWindowEarlierIndex, 0);
									binTestWindowIsORI.set(binTestWindowEarlierIndex);
									//BedGraphEntry bgePrevious = bedGraphEntries.get(bgeEarlierIndex); 									
									firstIndices.get(MutPOLEType.AGAtoATA).reset();  // finally reset the first index
									firstIndices.get(MutPOLEType.TCTtoTAT).reset();
								}
							}
							
						} else if (trackState.equals(PvalueTrackState.Negative)) {
							// We are not in a easy-to-detect ORI region
							if (firstIndices.get(MutPOLEType.TCTtoTAT).isAlreadySetOnce()) {
								if (midPoints.get(0).mMutType == MutPOLEType.AGAtoATA) {
									if (pValMultiplied >= -0.1) {
										// We backtrack windows until we reach an ORI
										int minIndex = -1;
										double minPVal = Double.MAX_VALUE;
										boolean foundUpwardSlope = false;
										for (int currentIndex = binTestWindows.size() - 2; (currentIndex >= 0) && !binTestWindowIsORI.get(currentIndex); currentIndex--) {
											double pValOther = binTestWindows.get(currentIndex).mValueWindow;
											if (pValOther < minPVal) {
												minPVal = pValOther;
												minIndex = currentIndex;
												foundUpwardSlope = false;
											} else if ((minIndex >= 0) && (pValOther > minPVal)) {
												foundUpwardSlope = true;
												break;
											}
										}
																				
										if (pValMultiplied >= minPVal) {
											if (foundUpwardSlope) {
												binTestWindowIsORI.set(minIndex);
											
												firstIndices.get(MutPOLEType.AGAtoATA).reset();  // finally reset the first index
												firstIndices.get(MutPOLEType.TCTtoTAT).reset();
												
												
												for (int rowTemp = row - 1, currentIndex = binTestWindows.size() - 2; 
													 rowTemp >= 0 && !binTestWindowIsORI.get(currentIndex); 
													 rowTemp--, currentIndex--) {
													if (bedItemList.get(rowTemp).mMutType == MutPOLEType.TCTtoTAT) {
														firstIndices.get(MutPOLEType.TCTtoTAT).resetAndSet(rowTemp);
														break;
													}
												}
												
											}
										}
									}
								}
							}
						} else {
							// The result is an easy to detect ORI 								
							if (trackState.equals(PvalueTrackState.ORI)) {
								binTestWindowIsORI.set(binTestWindows.size() - 1);
							} else {
								binTestWindowIsCollision.set(binTestWindows.size() - 1);
							}
							firstIndices.get(MutPOLEType.AGAtoATA).reset();  // finally reset the first index
							firstIndices.get(MutPOLEType.TCTtoTAT).reset();
							
						}
					}
					pValPrev.resetAndSet(pValMultiplied);
					
					try {
						bgPvalue.add("chr" + chrom.getName(), midPoint, midPoint, pValMultiplied);
						bgScaled.add("chr" + chrom.getName(), midPoint, midPoint, pValMultiplied);						
					} catch (IOException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}					
				}
			}
			
			
			// Ensure that we don't have overlapping ORI calls and collisions
			// Create a block so that we can deallocate
			{
				BitSet oriCollisionIntersectTest = (BitSet) binTestWindowIsORI.clone();
				oriCollisionIntersectTest.and(binTestWindowIsCollision);
				CompareUtils.ensureTrue(oriCollisionIntersectTest.cardinality() == 0, "ERROR: Have overlap between ORI and collision bitsets!");
			}
			
			// Filter the back and forth oscillation between ORI and collision points
			filterORICollisionBackAndForth(binTestWindowIsORI, binTestWindowIsCollision, binTestWindows);
			
			// Set up the ORI output
			BitSet bedGraphEntriesORIorCollision = new BitSet();
			bedGraphEntriesORIorCollision.or(binTestWindowIsORI);
			bedGraphEntriesORIorCollision.or(binTestWindowIsCollision);
			

			// Write ORI track
			BitSet bitsetToUse = binTestWindowIsORI;
			for (int i = bitsetToUse.nextSetBit(0); i >= 0; i = bitsetToUse.nextSetBit(i + 1)) {
				BinTestWindow bgeCurrent = binTestWindows.get(i);				
				try {
					double tickValue = binTestWindowIsORI.get(i) ? getTickValue(PvalueTrackState.ORI) : getTickValue(PvalueTrackState.Collision);
					bgORI.add("chr" + chrom.getName(), bgeCurrent.mRegionWindow.getRangeStart(), bgeCurrent.mRegionWindow.getRangeEnd(), tickValue);
				} catch (IOException e) {
					e.printStackTrace();
				}
			}
			
			// Calculate the inter-ORI distances
			PrimitiveWrapper.WIntegerSetOnce prevORIIndex = new PrimitiveWrapper.WIntegerSetOnce(); 
			for (int i = binTestWindowIsORI.nextSetBit(0); i >= 0; i = binTestWindowIsORI.nextSetBit(i + 1)) {
				if (prevORIIndex.isAlreadySetOnce()) {
					BinTestWindow bgeCurrent = binTestWindows.get(i);
					BinTestWindow bgeEarlier = binTestWindows.get(prevORIIndex.mInt);
					int distance = bgeCurrent.mRegionWindow.getRangeStart() - bgeEarlier.mRegionWindow.getRangeEnd();
					int distanceDiscretized = discretizeORIDistance(distance);
					oriDistCount.incrementCount(distanceDiscretized);
					oriDistCountRaw.incrementCount(distance);

				}
				prevORIIndex.resetAndSet(i);
			}
			
			// Perform permutation for ORI distance
			interORIpermTask.mBedGrahEntries = binTestWindows;
			performPermutation(binTestWindows.size(), binTestWindowIsORI.cardinality(), 50, interORIpermTask);
			
			// Do spline
			System.out.println(chrom);
			DoubleArrayList splineX = new DoubleArrayList(binTestWindows.size());
			DoubleArrayList splineY = new DoubleArrayList(binTestWindows.size());
			int indexAdded = 0;
			
			// Add known ORI and collision points for sure to spline list.  Remainder of points
			// are subject to spline increment
			for (int i = 0; i < binTestWindows.size(); i++) {
				boolean shouldAdd = false;
				if (binTestWindowIsORI.get(i) || binTestWindowIsCollision.get(i)) {
					shouldAdd = true;
				} else if (i == 0) {
					shouldAdd = true;
				} else if ((i - indexAdded) == SplineIncrement) {
					shouldAdd = true;
				} else if ((i - indexAdded) > SplineIncrement) {
					CompareUtils.ensureTrue(false, "ERROR: Impossible state!");
				}
				
				if (shouldAdd) {
					BinTestWindow bge = binTestWindows.get(i);
					splineX.add(bge.mRegionWindow.getRangeStart());
					splineY.add(bge.mValueWindow);
					indexAdded = i;
				}
			}
					
			// Add final point
			if (indexAdded < (binTestWindows.size() - 1)) {
				splineX.add(binTestWindows.get(binTestWindows.size() - 1).mRegionWindow.getRangeStart());
				splineY.add(binTestWindows.get(binTestWindows.size() - 1).mValueWindow);
			}
			//System.out.println(splineX);

			if (splineX.size() > 10) {				
				//Spline spline = new Spline(splineX, splineY);
				SplineInterpolator splineInterpolator = new SplineInterpolator();
				PolynomialSplineFunction psf = splineInterpolator.interpolate(splineX.toArray(), splineY.toArray());
				
				for (BinTestWindow bge : binTestWindows) {
					try {
						double value = psf.value(bge.mRegionWindow.getRangeStart());
						bgSpline.add("chr" + chrom.getName(), bge.mRegionWindow.getRangeStart(), bge.mRegionWindow.getRangeStart(), Math.max(Math.min(value, 0.5), -0.5));
					} catch (IOException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
				}
			}
		}
		
		// We finished processing all the chromosomes.  Now write the output
		String bedGraphPval           = outFilenameRoot + ".ORI.windowSize." + windowSize + ".pval.txt";
		String bedGraphPvalMinusLog10 = outFilenameRoot + ".ORI.windowSize." + windowSize + ".pval.minusLog10.txt";
		String bedGraphSpline         = outFilenameRoot + ".ORI.windowSize." + windowSize + ".pval.spline." + SplineIncrement + ".txt";
		String bedGraphORI            = outFilenameRoot + ".ORI.windowSize." + windowSize + ".pval.ORIs.txt";
		
		String interORIDistance       = outFilenameRoot + ".ORI.windowSize." + windowSize + ".interORIDistance.txt";
		
		String transcriptPositions    = outFilenameRoot + ".ORI.windowSize." + windowSize + ".transcriptSites.txt";
		parseTranscriptAnnotationsByChrom(transcriptStartSiteFile, transcriptPositions, binTestWindowsOnePatient, oriIndicesOnePatient);
		
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

		IntArrayList   cumulativeKeysBuffer     = oriDistCountRaw.getKeys(null);		
		IntArrayList   cumulativeCountsBuffer   = oriDistCountRaw.getCumulativeSumsForward(null);
		FloatArrayList cumulativeFractionBuffer = oriDistCountRaw.getCumulativeProportionsForward(null);
		
		xyDatasetCumCount.addSeries(sampleName,    ArrayUtils.combineTwoDyanamicArraysIntoOneStatisDouble(cumulativeKeysBuffer, cumulativeCountsBuffer));
		xyDatasetCumFraction.addSeries(sampleName, ArrayUtils.combineTwoDyanamicArraysIntoOneStatisDouble(cumulativeKeysBuffer, cumulativeFractionBuffer));
		
		// Now do for the bucketed counts
		IntArrayList keysBuffer                  = oriDistCount.getKeys(null);
		IntArrayList   countsBufferDiscretized   = oriDistCount.getCounts(null);
		FloatArrayList fractionBufferDiscretized = oriDistCount.getProportions(null);
		
		xyDatasetCount.addSeries(sampleName,    ArrayUtils.combineTwoDyanamicArraysIntoOneStatisDouble(keysBuffer, countsBufferDiscretized));
		xyDatasetFraction.addSeries(sampleName, ArrayUtils.combineTwoDyanamicArraysIntoOneStatisDouble(keysBuffer, fractionBufferDiscretized));
		
		// Now do for permutation
		IntArrayList       keysBufferPerm           = interORIpermTask.mCounterDiscretized.getKeys(null);
		FloatArrayList fractionBufferPerm           = interORIpermTask.mCounterDiscretized.getProportions(null);
		
		IntArrayList       keysBufferPermCumulative = interORIpermTask.mCounterRaw.getKeys(null);
		FloatArrayList fractionBufferPermCumulative = interORIpermTask.mCounterRaw.getCumulativeProportionsForward(null);
		
		DefaultXYDataset xyDatasetPermCumulative = new DefaultXYDataset();
		xyDatasetPermCumulative.addSeries(sampleName,    ArrayUtils.combineTwoDyanamicArraysIntoOneStatisDouble(cumulativeKeysBuffer, cumulativeFractionBuffer));		
		xyDatasetPermCumulative.addSeries("Permutation", ArrayUtils.combineTwoDyanamicArraysIntoOneStatisDouble(keysBufferPermCumulative, fractionBufferPermCumulative));
		String outFilenamePerm = outDir + File.separator + sampleName + ".Permutation.Cumulative";
		createGraphs(xyDatasetPermCumulative, outFilenamePerm, false, false);
		
		DefaultXYDataset xyDatasetPerm = new DefaultXYDataset();
		xyDatasetPerm.addSeries(sampleName,    ArrayUtils.combineTwoDyanamicArraysIntoOneStatisDouble(keysBuffer, fractionBufferDiscretized));		
		xyDatasetPerm.addSeries("Permutation", ArrayUtils.combineTwoDyanamicArraysIntoOneStatisDouble(keysBufferPerm, fractionBufferPerm));
		outFilenamePerm = outDir + File.separator + sampleName + ".Permutation";
		createGraphs(xyDatasetPerm, outFilenamePerm, false, false);
		
		
		//xyDatasetPerm.addSeries("Permutation", ArrayUtils.combineTwoDyanamicArraysIntoOneStatisDouble(keysBufferPerm, fractionBufferPerm));
		
		BufferedWriter outORIDistance = IOUtils.getBufferedWriter(interORIDistance);
		//IntArrayList keysBuffer = new IntArrayList();
		//oriDistCount.getKeys(keysBuffer);
		//for (int i = 0; i < keysBuffer.size(); i++) {
		for (int i = 0; i < 9; i++) {
			int key = Cast.toInt(Math.round(Math.pow(10, i)));
			//int key = keysBuffer.get(i);
			int count = oriDistCount.getCount(key);
			IOUtils.writeToBufferedWriter(outORIDistance, key + "\t" + count, true);
		}
		IOUtils.closeBufferedWriter(outORIDistance);
		
		//track name="ColorByStrandDemo" description="Color by strand demonstration" visibility=2 colorByStrand="255,0,0 0,0,255"
		
		
	}

	// ===========================
	private static void filterORICollisionBackAndForth(BitSet bsORI, BitSet bsCollision, ArrayList<BinTestWindow> bedGraphEntries) {
		BitSet bedGraphEntriesORIorCollision = new BitSet();
		bedGraphEntriesORIorCollision.or(bsORI);
		bedGraphEntriesORIorCollision.or(bsCollision);
		
		double leeway = 0.05;
		int gapThreshold = 2;
		PrimitiveWrapper.WDoubleSetOnce pValPrev = new PrimitiveWrapper.WDoubleSetOnce();
		
		// Filter out the back-and-forth jumping that can occur after a collision track
		for (int indexCollision = bsCollision.nextSetBit(0); indexCollision >= 0; indexCollision = bsCollision.nextSetBit(indexCollision + 1)) {
			// We found our collision index.  Now we scan across the ORI indices to make sure we don't
			// have a concentration of back-and-forth collision and ORI indices
			
			// We keep a "previous" index
			int prevORIOrCollisionIndex = indexCollision;
			pValPrev.resetAndSet(bedGraphEntries.get(indexCollision).mValueWindow);
			
			for (int nextIndex = prevORIOrCollisionIndex + 1; (nextIndex - prevORIOrCollisionIndex <= gapThreshold) && (nextIndex < bedGraphEntries.size()); nextIndex++) {
				double pValCurr = bedGraphEntries.get(nextIndex).mValueWindow;
				if (bsORI.get(nextIndex) || bsCollision.get(nextIndex) || CompareUtils.inRange(pValCurr, pValPrev.mDouble, leeway)) {
					prevORIOrCollisionIndex = nextIndex;
					pValPrev.resetAndSet(pValCurr);
				}
			}
			
			// Take action if we have semi-contiguous ORI and collision signals
			if (prevORIOrCollisionIndex > indexCollision) {
				bsORI.clear      (indexCollision + 1, prevORIOrCollisionIndex + 1);
				bsCollision.clear(indexCollision + 1, prevORIOrCollisionIndex + 1);				
			}
			
			// Increment
			indexCollision = prevORIOrCollisionIndex;
		}
	}
	
	// ===========================
	private static class InterORIDistancePermutationTask implements BitSetPermutationTask {

		DynamicBucketCounter mCounterDiscretized;
		DynamicBucketCounter mCounterRaw;
		ArrayList<BinTestWindow> mBedGrahEntries;
		PrimitiveWrapper.WIntegerSetOnce mPrevIndex;
		
		public InterORIDistancePermutationTask(ArrayList<BinTestWindow> bedGraphEntries) {
			mCounterDiscretized = new DynamicBucketCounter();
			mCounterRaw         = new DynamicBucketCounter();
			mBedGrahEntries = bedGraphEntries;
			mPrevIndex = new PrimitiveWrapper.WIntegerSetOnce(); 
		}
		
		@Override
		public void takeAction(boolean[] bitArray, int permutationNumber) {
			mPrevIndex.reset();
			for (int i = BitSetUtils.nextTrueBoolean(bitArray, 0); i >= 0; i = BitSetUtils.nextTrueBoolean(bitArray, i + 1)) {
				if (mPrevIndex.isAlreadySetOnce()) {
					BinTestWindow bgePrev = mBedGrahEntries.get(mPrevIndex.mInt);
					BinTestWindow bgeCurr = mBedGrahEntries.get(i);
					int distance = bgeCurr.mRegionWindow.getRangeStart() - bgePrev.mRegionWindow.getRangeStart();
					int distanceDiscretized = discretizeORIDistance(distance);
					mCounterDiscretized.incrementCount(distanceDiscretized);
					mCounterRaw.incrementCount(distance);
				}
				mPrevIndex.resetAndSet(i);
			}
			
		}
	}

	// ===========================
	public static interface BitSetPermutationTask {
		public void takeAction(boolean[] bitArray, int permutationNumber);
	}
	
	// ===========================
	private static void performPermutation(int numSlots, int numSlotsDesiredTrue, int numIterations, BitSetPermutationTask permutationTask) {
		boolean[] bitArray = new boolean[numSlots];
		Arrays.fill(bitArray, false);
		
		for (int iter = 0; iter < numIterations; iter++) {
			BitSetUtils.setBooleanArrayWithRandomValues(bitArray, numSlotsDesiredTrue);
			permutationTask.takeAction(bitArray, iter);
		}
	}
	
	// ===========================
	private static int discretizeORIDistance(int distance) {
		// We round to the nearest multiple of nearest log10
		int floorBase = Cast.toInt(NumberUtils.MathLog10Safe(distance));  // get the floor
		int log10Rounded = NumberUtils.roundToNearestNthPositivePowerOf10(distance, Math.min(floorBase, 4));
		//System.out.println(distance + "\t" + log10Rounded);
							
		//int log10 = (int) Math.round(Math.pow(10, Cast.toInt(NumberUtils.MathLog10Safe(distance))));
		//int log10 = Cast.toInt( Math.round( Math.pow(10, Math.round(NumberUtils.MathLog10Safe(distance)) ) ) );
		return log10Rounded;
	}
	
	public static void convertHeader() {
		
	}
	
	// ===========================
	public static int getMultiplierPvalue(MutPOLEType mpt) {
		return (mpt == MutPOLEType.AGAtoATA) ? -1 : 1;
	}

	// ===========================
	public static PvalueTrackState isORI(double pValScaledPrev, double pValScaledCurr) {
		if (pValScaledPrev >= 0 && pValScaledCurr >= 0) {
			return PvalueTrackState.Positive;
		} else if (pValScaledPrev < 0 && pValScaledCurr < 0) {
			return PvalueTrackState.Negative;
		} else if (pValScaledPrev > 0 && pValScaledCurr < 0) {
			return PvalueTrackState.ORI;
		} else if (pValScaledPrev < 0 && pValScaledCurr > 0) {
			return PvalueTrackState.Collision;
		} else {
			System.out.println("Huha!");
			return null;
		}
	}
	
	
	// ===========================
	public static int getMidPoint(ArrayBlockingQueue<MutBEDLine> abq, ArrayList<MutBEDLine> midPoints) {
		/*
		long sum = 0;
		
		for (Iterator<MutBEDLine> iter = abq.iterator(); iter.hasNext(); ) {
			sum += iter.next().mPosStart;
		}
		
		return (int) Math.round((double) sum / ((double) abq.size()));
		*/
		int abqSize = abq.size();
		boolean isOdd = (abqSize & 1) == 1;		
		int midPoint = abqSize >>> 1;
		
		int index = 0;
		MutBEDLine mblPrev = null;
		
		if (midPoints != null) { midPoints.clear(); }
		
		for (Iterator<MutBEDLine> iter = abq.iterator(); iter.hasNext(); index++) {
			MutBEDLine mbl = iter.next();
			
			if (index == midPoint) {
				if (isOdd) {
					if (midPoints != null) { 
						midPoints.add(mbl); 
					}
					return mbl.mPosStart;
				} else {
					if (midPoints != null) {
						midPoints.add(mblPrev);
						midPoints.add(mbl);						
					}
					return ((mblPrev.mPosStart + mbl.mPosStart) >>> 1);
				}
			}
			mblPrev = mbl;
		}
		
		CompareUtils.ensureTrue(false, "ERROR: Cannot reach this point!");
		return -1;
	}


	// ===========================
	private static double getTickValue(PvalueTrackState trackState) {
		switch(trackState) {
		case ORI:       return -1;
		case Collision: return  1;
		default:        return  0;
		}		
	}
	
	// ===========================
	public static class BinTestWindow {
		RegionRangeWithPayload<PrimitiveWrapper.WInteger> mRegionWindow;
		double mValueWindow;

		public BinTestWindow(Chrom chrom, int position, double value) {
			this(chrom, position, position, value);
		}
		
		public BinTestWindow(Chrom chrom, int positionStart, int positionEnd, double value) {
			mRegionWindow = new RegionRangeWithPayload<>(chrom, positionStart, positionEnd, new PrimitiveWrapper.WInteger(0));
			mValueWindow = value;			
		}
	}

	// ===========================
	public static enum PvalueTrackState { 
		Positive,
		Negative,
		ORI,
		Collision;
	}
	
	// ===========================
	public static class MutBEDLine extends RegionRangeWithPayload<MutPOLEType> {
		Chrom mChrom;
		int mPosStart, mPosEnd;
		MutPOLEType mMutType;		
	}

	// ===========================
	public static void createGraphs(DefaultXYDataset xyDatasetCumulCounts, String outFilename, boolean isCount, boolean isCumulative) {
		
		String yAxisLabel = isCount ? "Count" : "Fraction";
		String title = (isCumulative ? "Cumulative " : "") + yAxisLabel + " Per Sample"; 
		JFreeChart xyLineChart = ChartFactory.createXYLineChart(title, "Inter-ORI Distance", yAxisLabel, xyDatasetCumulCounts, PlotOrientation.VERTICAL, true, true, false);
		
		XYPlot xyPlot = (XYPlot) xyLineChart.getPlot();
		
		//XYItemRenderer itemRenderer = ClusteringPlotting.getXYItemRendererHelper(1);		
		//setSeriesPaintPerCluster(itemRenderer, xyDataset, EventType.class, ColorPastelForEnumEventType);
		//xyPlot.setRenderer(itemRenderer);		
		
		xyPlot.setBackgroundPaint(ColorPastel.Gray_05.getColor());
		xyPlot.setDomainGridlinePaint(Color.white);
		xyPlot.setRangeGridlinePaint(Color.white);		
		
		LegendTitle legendTitle = xyLineChart.getLegend();
		legendTitle.setID("Event");
		legendTitle.setItemFont(new Font("Arial", Font.BOLD, 20));
	
		LogAxis logAxis = new LogAxis("Inter-ORI Distance (log10)");
		//LogarithmicAxis logAxis = new LogarithmicAxis("Inter-ORI Distance (log10)");
		DecimalFormat dc = new DecimalFormat();
		//dc.setGroupingSize(3);
		dc.setGroupingUsed(true);					
		logAxis.setNumberFormatOverride(dc);
		//logAxis.setExpTickLabelsFlag(true);
		
		xyPlot.setDomainAxis(logAxis);
		
		XYSplineRenderer splineRenderer = new XYSplineRenderer(2);
		//xyPlot.setRenderer(splineRenderer);
		
		if (!isCount) {			
			//xyPlot.getRangeAxis().setRange(0, 1.02);
		}
		
		Font rangeAxisLabelFont = new Font("Arial", Font.BOLD, 20);
		Font rangeAxisTickFont = new Font("Arial", Font.BOLD, 20);
		xyPlot.getRangeAxis().setLabelFont(rangeAxisLabelFont);		
		xyPlot.getRangeAxis().setTickLabelFont(rangeAxisTickFont);
		//TickUnit tickUnit =  new NumberTickUnit(1.0);
		//TickUnits tickUnits = new TickUnits();
		//tickUnits.add(tickUnit);		
		//xyPlot.getRangeAxis().setStandardTickUnits(tickUnits);
		xyPlot.getDomainAxis().setLabelFont(rangeAxisLabelFont);
		xyPlot.getDomainAxis().setTickLabelFont(new Font("Arial", Font.BOLD, 24));	
		GraphUtils.saveChartAsPNG(outFilename, xyLineChart, 1200, 450);
	}

	// ===========================
	private static void parseTranscriptAnnotationsByChrom_ParseLine(String line, AnnotationLine annot) {
		int leeway = 1000;
		String chromString = StringUtils.extractNthColumnValue(line, 1, StringUtils.TabStr);
		Chrom chrom = Chrom.getChrom(chromString);
		
		int posStart = Integer.parseInt(StringUtils.extractNthColumnValue(line, 4, StringUtils.TabStr));
		int posEnd   = Integer.parseInt(StringUtils.extractNthColumnValue(line, 5, StringUtils.TabStr));
		
		annot.mRangeStrict.set(  chrom, posStart,          posEnd,          true, 0);
		annot.mRangeFlexible.set(chrom, posStart - leeway, posEnd + leeway, true, 0);
		
		String targetGeneNamePrefixString = "gene_name \"";
		int indexGeneName = line.indexOf(targetGeneNamePrefixString) + targetGeneNamePrefixString.length() + 1;
		int indexEndQuote = line.indexOf("\"", indexGeneName);
		annot.mGeneName = line.substring(indexGeneName, indexEndQuote);
	}

	// ===========================
	private static class AnnotationLine {
		RegionRange<?> mRangeStrict = new RegionRange.Default();
		RegionRange<?> mRangeFlexible = new RegionRange.Default();
		String mGeneName;		
	}
	
	// ===========================
	private static void parseTranscriptAnnotationsByChrom(
			String annotationsFile, 
			String outFileName, 
			EnumMapSafe<Chrom, ArrayList<BinTestWindow>> bedGraphEntries,
			EnumMapSafe<Chrom, BitSet> oriIndices) {
		
		BufferedReader inAnnot = IOUtils.getBufferedReader(annotationsFile);
		BufferedWriter out = IOUtils.getBufferedWriter(outFileName, false);		
		AnnotationLine annotLine = new AnnotationLine();
		String line = IOUtils.getNextLineInBufferedReader(inAnnot);
		parseTranscriptAnnotationsByChrom_ParseLine(line, annotLine);
		
		for (Chrom chrom : Chrom.values()) {
			ArrayList<BinTestWindow> bedGraphEntriesChrom = bedGraphEntries.get(chrom);
			BitSet oriIndicesChrom = oriIndices.get(chrom);
			
			for (int i = 0; i < bedGraphEntriesChrom.size(); ) {
				BinTestWindow bge = bedGraphEntriesChrom.get(i);
				
				if (annotLine.mRangeFlexible.beforeRange(chrom, bge.mRegionWindow.getRangeStart())) {
					i++;  // increment to next point
				} else if (annotLine.mRangeFlexible.afterRange(chrom, bge.mRegionWindow.getRangeStart())) {
					line = IOUtils.getNextLineInBufferedReader(inAnnot);					
					if (line == null) break;					
					parseTranscriptAnnotationsByChrom_ParseLine(line, annotLine);
				} else if (annotLine.mRangeFlexible.inRange(chrom, bge.mRegionWindow.getRangeStart())) {
					if (oriIndicesChrom.get(i)) {
						String outLine = chrom.ordinal() + "\t" + bge.mRegionWindow.getRangeStart() + "\t" + annotLine.mGeneName + "\t" + line;
						IOUtils.writeToBufferedWriter(out, outLine, true);						
					}
					i++;
				} else {
					CompareUtils.ensureTrue(false, "Impossible State!");
				}
			}
		}
		
		IOUtils.closeBufferedReader(inAnnot);
		IOUtils.closeBufferedWriter(out);
	}
	
	
	// ===========================
	public static void workOnDirectory(String inDir, String outDirRoot, String transcriptStartSiteFile, int[] windowSizes) {
		File[] files = (new File(inDir)).listFiles();
		ArrayList<File> validFiles = new ArrayList<>();
		
		for (File file : files) {
			if (file.getName().endsWith(".BED.TXT")) {
				validFiles.add(file);
			}
		}
		
		IOUtils.createDirectoryPath(outDirRoot, false);		

		
		for (int windowSize : windowSizes) {
			String outDir = outDirRoot + File.separator + "Window_" + windowSize;
			IOUtils.createDirectoryPath(outDir, false);
			
			DefaultXYDataset xyDatasetCumCount = new DefaultXYDataset();
			DefaultXYDataset xyDatasetCumFraction = new DefaultXYDataset();
			DefaultXYDataset xyDatasetCount = new DefaultXYDataset();
			DefaultXYDataset xyDatasetFraction = new DefaultXYDataset();
			
			for (File file : validFiles) {								
				doWork(file, outDir, transcriptStartSiteFile, windowSize, xyDatasetCumCount, xyDatasetCumFraction, xyDatasetCount, xyDatasetFraction);
			}
			
			String outFilenameRoot = outDir + File.separator + "Inter_ORI.Cumulative.Counts.window." + windowSize;
			createGraphs(xyDatasetCumCount, outFilenameRoot, true, true);
			
			outFilenameRoot = outDir + File.separator + "Inter_ORI.Cumulative.Fractions.window." + windowSize;
			createGraphs(xyDatasetCumFraction, outFilenameRoot, false, true);
			
			outFilenameRoot = outDir + File.separator + "Inter_ORI.Counts.window." + windowSize;
			createGraphs(xyDatasetCount, outFilenameRoot, true, false);
			
			outFilenameRoot = outDir + File.separator + "Inter_ORI.Fractions.window." + windowSize;
			createGraphs(xyDatasetFraction, outFilenameRoot, false, false);

		}
	}
	
	// ===========================
	public static void main(String[] args) {
		int[] windowSizes = new int[] { /*11, */ 20 /*, 31, 41*/ };
		workOnDirectory(args[0], args[1], args[2], windowSizes);
	}
}



