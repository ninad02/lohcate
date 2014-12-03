package shared;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Font;
import java.io.BufferedWriter;
import java.io.File;
import java.util.ArrayList;
import java.util.BitSet;

import lohcate.clustering.ClusteringPlotting;
import lohcateEnums.ColorPastel;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.ValueMarker;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYItemRenderer;
import org.jfree.chart.title.LegendTitle;
import org.jfree.data.xy.DefaultXYDataset;
import org.jfree.data.xy.XYDataset;

import com.carrotsearch.hppc.FloatArrayList;

import cern.colt.list.DoubleArrayList;

import nutils.GraphUtils;
import nutils.IOUtils;
import nutils.NumberUtils;
import nutils.PrimitiveWrapper;
import nutils.StringUtils;
import nutils.counter.DynamicRoundedDoubleCounter;

public class GeneMutationPermutation {

	// ========================================================================
	private static enum SamplingMode {
		ByLength,
		ByCount;
	}
	
	// ========================================================================	
	public static final boolean CollapseGenesIntoOne = true;
	public static final boolean UseCaptureLengths = true;
	
	public static void doWork(String inFileTotal, String inFileTarget) {
		System.out.println(inFileTotal);
		System.out.println(inFileTarget);
		
		ArrayList<GeneMutationLengthAndCount> totalList = readGeneMutationLengthAndCountsFromFile(inFileTotal);
		ArrayList<GeneMutationLengthAndCount> targetList = readGeneMutationLengthAndCountsFromFile(inFileTarget);
		
		PrimitiveWrapper.WInteger totalLengthTarget   = new PrimitiveWrapper.WInteger();
		PrimitiveWrapper.WInteger totalMutCountTarget = new PrimitiveWrapper.WInteger();
		double mutFreqTarget = calcMutFrequency(targetList, totalLengthTarget, totalMutCountTarget);
		
		PrimitiveWrapper.WInteger totalLengthTotal   = new PrimitiveWrapper.WInteger();
		PrimitiveWrapper.WInteger totalMutCountTotal = new PrimitiveWrapper.WInteger();
		double mutFreqTotal = calcMutFrequency(totalList, totalLengthTotal, totalMutCountTotal);
		
		System.out.printf("Mut Freq Total: %g\n", mutFreqTotal);
		System.out.printf("Mut Freq Target: %g\n", mutFreqTarget);
		System.out.printf("Mut Count Target: %d\n", totalMutCountTarget.mInt);
		System.out.printf("Mut Length Target: %d\n", totalLengthTarget.mInt);
		System.out.printf("Gene Count: %d\n", targetList.size());
		
		DefaultXYDataset xyDataset = new DefaultXYDataset();
		double maxPval = -Double.MAX_VALUE;
		//double mutFreqTotal  = calcMutFrequency(totalList);		
		
		for (SamplingMode samplingMode : SamplingMode.values()) {
			System.out.printf("--------------\nSampling Mode: %s\n--------------\n", samplingMode);

			int numIter = 1_000_000;
			DoubleArrayList mutFreqList = new DoubleArrayList(numIter);
			ArrayList<GeneMutationLengthAndCount> listForTrial = new ArrayList<GeneMutationLengthAndCount>();
			PrimitiveWrapper.WInteger totalLengthForTrial   = new PrimitiveWrapper.WInteger();
			PrimitiveWrapper.WInteger totalMutCountForTrial = new PrimitiveWrapper.WInteger();

			int criterion = (samplingMode == SamplingMode.ByLength) ? totalLengthTarget.mInt : targetList.size();
			for (int iter = 0; iter < numIter; iter++) {							
				//if (iter % 100000 == 0) { System.out.printf("Iter: %d\t%d\n", iter, criterion); }
				performRandomization(totalList, criterion, listForTrial, samplingMode);
				double mutFreq = calcMutFrequency(listForTrial, totalLengthForTrial, totalMutCountForTrial);
				mutFreqList.add(mutFreq);
			}

			// Now sort the mutation frequencies
			mutFreqList.sort();

			//
			BufferedWriter out = IOUtils.getBufferedWriter(inFileTarget + ".mutFreqs.txt");
			for (int i = 0; i < mutFreqList.size(); i++) {
				IOUtils.writeToBufferedWriter(out, "" + mutFreqList.get(i), true);			
			}
			IOUtils.closeBufferedWriter(out);

			// Now go from end to the beginning and find the index whose element is less than the target mutation frequency
			int theIndex = mutFreqList.size() - 1; 
			for (; theIndex >= 0; theIndex--) {
				if (mutFreqList.get(theIndex) < mutFreqTarget) {
					break;
				}
			}

			double pVal = 1;
			if (theIndex >= 0) {
				int numElem = mutFreqList.size() - theIndex;
				pVal = (double) numElem / (double) mutFreqList.size();
			}
			maxPval = Math.max(maxPval, pVal);			
			System.out.printf("P-Value: %g\n", pVal);

			DynamicRoundedDoubleCounter roundedCounter = new DynamicRoundedDoubleCounter(2);
			for (int i = 0; i < mutFreqList.size(); i++) {
				roundedCounter.register(mutFreqList.get(i));
			}

			com.carrotsearch.hppc.DoubleArrayList mutFreqUnique = new com.carrotsearch.hppc.DoubleArrayList();
			FloatArrayList mutFreqUniqueProportions = new FloatArrayList();
			roundedCounter.getValuesWithProportions(mutFreqUnique, mutFreqUniqueProportions);


			double[][] data = new double[2][mutFreqUnique.size()];			
			for (int i = 0; i < mutFreqUnique.size(); i++) {
				//System.out.printf("%g\t%g\n", mutFreqUnique.get(i), mutFreqUniqueProportions.get(i));
				data[0][i] = mutFreqUnique.get(i);
				data[1][i] = mutFreqUniqueProportions.get(i);
			}
			String suffix = " (";
			switch(samplingMode) {
			case ByLength: suffix += Math.round(criterion / 1000) + " kb)"; break;
			case ByCount:  suffix += criterion + " genes)"; break;
			}
			xyDataset.addSeries(samplingMode + suffix, data);
		}
		
		String theName = (new File(inFileTarget)).getName();
		String thePrefix = theName.substring(0, theName.indexOf("-mutated"));
		double mutFreqRounded = DynamicRoundedDoubleCounter.roundToNumSigFigs(mutFreqTarget, 2);		
		plotPermutationResults(xyDataset, thePrefix, inFileTarget + ".permuted", mutFreqTarget, "Target MF = " + mutFreqRounded + ", p = " + maxPval);
	}

	// ========================================================================
	private static BitSet ElementDone = new BitSet(65536); 
	public static synchronized void performRandomization(ArrayList<GeneMutationLengthAndCount> theList, final int totalRegionLength, ArrayList<GeneMutationLengthAndCount> returnList, SamplingMode samplingMode) {
		ElementDone.clear();
		returnList.clear();
		
		int lastIndex = theList.size() - 1;
		for (int lengthSoFar = 0; lengthSoFar < totalRegionLength; ) {
			int randomIndex = NumberUtils.getRandomInteger(0, lastIndex);
			if (!ElementDone.get(randomIndex)) {				
				ElementDone.set(randomIndex);
				GeneMutationLengthAndCount glmc = theList.get(randomIndex);
				returnList.add(glmc);
				lengthSoFar += ((samplingMode == SamplingMode.ByLength) ? glmc.mLength : 1);
			}			
		}		
	}
	
	// ========================================================================
	/** Plots null permutation distribution, along with the x-axis value. */
	public static void plotPermutationResults(XYDataset xyDataset, String title, String outFilenameRoot, double valueOfInterest, String valueOfInterestLabel) {
		String xAxisLabel = "Mutation Frequency in Permutation Trial i (MF_i)";
		String yAxisLabel = "Fraction of Iterations with Given MF_i";
		
		JFreeChart theChart = 
				//ChartFactory.createXYLineChart(title, xAxisLabel, yAxisLabel, xyDataset, PlotOrientation.VERTICAL, true, false, false);
				ChartFactory.createScatterPlot(title, xAxisLabel, yAxisLabel, xyDataset, PlotOrientation.VERTICAL, true, false, false);				
		XYPlot xyPlot = (XYPlot) theChart.getPlot();
		
		// itemRenderer = ClusteringPlotting.getXYItemRendererHelper(7);		
		//setSeriesPaintPerCluster(itemRenderer);
		//xyPlot.setRenderer(itemRenderer);
		ValueMarker valueOfInterestMarker = new ValueMarker(valueOfInterest, ColorPastel.Gray_90.getColor(), new BasicStroke());
		valueOfInterestMarker.setLabel(valueOfInterestLabel);
		valueOfInterestMarker.setLabelAnchor(org.jfree.ui.RectangleAnchor.CENTER);
		valueOfInterestMarker.setLabelFont(new Font("Arial", Font.BOLD, 14));		
		xyPlot.addDomainMarker(valueOfInterestMarker);
		
		xyPlot.setBackgroundPaint(ColorPastel.Gray_15.getColor());
		xyPlot.setDomainGridlinePaint(Color.white);
		xyPlot.setRangeGridlinePaint(Color.white);		
		
		LegendTitle legendTitle = theChart.getLegend();
		legendTitle.setID("Permutation Mode");
		legendTitle.setItemFont(new Font("Arial", Font.BOLD, 20));
		
		Font rangeAxisLabelFont = new Font("Arial", Font.BOLD, 20);
		Font rangeAxisTickFont = new Font("Arial", Font.BOLD, 20);
		//xyPlot.getRangeAxis().setRange(0, 1.02);
		//xyPlot.getDomainAxis().setRange(0.0002, 0.0007);
		xyPlot.getRangeAxis().setLabelFont(rangeAxisLabelFont);		
		xyPlot.getRangeAxis().setTickLabelFont(rangeAxisTickFont);
		xyPlot.getDomainAxis().setLabelFont(rangeAxisLabelFont);
		xyPlot.getDomainAxis().setTickLabelFont(rangeAxisTickFont);		
	
		// Now write the plot		
		int width  = 1000;
		int height = 1000;
		GraphUtils.saveChartAsPNG(outFilenameRoot, theChart, width, height);
	}
	
	// ========================================================================
	public static double calcMutFrequency(ArrayList<GeneMutationLengthAndCount> gmlcList, PrimitiveWrapper.WInteger totalLength, PrimitiveWrapper.WInteger totalMutCount) {
		totalMutCount.mInt = 0;
		totalLength.mInt = 0;
		
		for (GeneMutationLengthAndCount gmlc : gmlcList) {
			totalMutCount.mInt += gmlc.mCount;
			totalLength.mInt   += gmlc.mLength;
		}
		
		return (double) totalMutCount.mInt / (double) totalLength.mInt;
	}
	
	// ========================================================================
	public static ArrayList<GeneMutationLengthAndCount> readGeneMutationLengthAndCountsFromFile(String inFilename) {
		ArrayList<String> allLines = IOUtils.readAllLinesFromFile(inFilename);
		
		ArrayList<GeneMutationLengthAndCount> gmlcList = new ArrayList<GeneMutationLengthAndCount>(allLines.size() * 5);
		for (String line : allLines) {
			String[] cols = line.split(StringUtils.TabPatternStr);
			//System.out.println(line);
			//if (cols.length < 4) continue;
			
			int theCount = Integer.parseInt(cols[0]);
			int colLength = UseCaptureLengths ? 4 : 3; 
					
			if (CollapseGenesIntoOne) {
				GeneMutationLengthAndCount gmlc = new GeneMutationLengthAndCount(cols[1], cols[2], Integer.parseInt(cols[colLength]), theCount);
				gmlcList.add(gmlc);
			} else {
				for (int i = 0; i < theCount; i++) {
					GeneMutationLengthAndCount gmlc = new GeneMutationLengthAndCount(cols[1], cols[2], Integer.parseInt(cols[colLength]), 1);
					gmlcList.add(gmlc);
				}		
			}
		}
		return gmlcList;
	}
	
	// ========================================================================
	private static class GeneMutationLengthAndCount {
		public String mGeneName;
		public String mTranscriptID;
		public int mCount;
		public int mLength;
		
		public GeneMutationLengthAndCount(String geneName, String transcriptID, int length, int count) {
			mCount = count;
			mLength = length;
			mGeneName = geneName;
			mTranscriptID = transcriptID;
		}
	}

	// ========================================================================
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		for (int i = 1; i < args.length; i++) {
			doWork(args[0], args[i]);
		}
		
	}

}
