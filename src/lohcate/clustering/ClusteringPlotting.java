package lohcate.clustering;

import genomeEnums.Chrom;
import genomeUtils.MaxPositionsOnChromosome;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Font;
import java.awt.Rectangle;
import java.awt.geom.Ellipse2D;
import java.awt.geom.Line2D;
import java.awt.geom.Rectangle2D;
import java.awt.geom.Ellipse2D.Double;
import java.io.File;
import java.util.ArrayList;
import java.util.EnumMap;
import java.util.Iterator;
import java.util.Set;

import lohcateEnums.EventType;
import lohcateEnums.ColorPastel;
import nutils.ArrayUtils;
import nutils.ArgumentParserUtils.InputParameterBoolean;
import nutils.ArrayUtils.ParallelArrayDoubleDynamic;
import nutils.ArgumentParserUtils;
import nutils.ControlFlagBool;
import nutils.EnumMapSafe;
import nutils.GraphUtils;
import nutils.PrimitiveWrapper;
import nutils.PrimitiveWrapper.WInteger;
import nutils.counter.DynamicBucketCounter;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.LegendItem;
import org.jfree.chart.LegendItemCollection;
import org.jfree.chart.annotations.XYBoxAnnotation;
import org.jfree.chart.axis.NumberTickUnit;
import org.jfree.chart.axis.TickUnit;
import org.jfree.chart.axis.TickUnits;
import org.jfree.chart.plot.IntervalMarker;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.ValueMarker;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.StackedXYAreaRenderer;
import org.jfree.chart.renderer.xy.StackedXYAreaRenderer2;
import org.jfree.chart.renderer.xy.StackedXYBarRenderer;
import org.jfree.chart.renderer.xy.XYItemRenderer;
import org.jfree.chart.renderer.xy.XYShapeRenderer;
import org.jfree.chart.title.LegendTitle;
import org.jfree.chart.title.TextTitle;
import org.jfree.data.xy.DefaultXYDataset;
import org.jfree.data.xy.XYBarDataset;
import org.jfree.data.xy.XYDataset;
import org.jfree.ui.HorizontalAlignment;
import org.jfree.ui.Layer;

import com.carrotsearch.hppc.DoubleArrayList;

public class ClusteringPlotting {

	private static final String ChromTitleString = buildChromString();
	
	// ========================================================================
	/** Plots the VAF (variant allele frequency) of the normal tissue comapred to the tumor tissue. */
	public static void plotVAFComparison(XYDataset xyDataset, String outFilenameRoot, String sampleName) {
		String xAxisLabel = "VAF Tumor";
		String yAxisLabel = "VAF Normal";
		String title = "VAF Comparison: " + sampleName;
		
		JFreeChart theChart = 
				ChartFactory.createScatterPlot(title, xAxisLabel, yAxisLabel, xyDataset, PlotOrientation.VERTICAL, true, false, false);
		XYPlot xyPlot = (XYPlot) theChart.getPlot();
		
		XYItemRenderer itemRenderer = getXYItemRendererHelper(5);		
		setSeriesPaintPerCluster(itemRenderer);
		xyPlot.setRenderer(itemRenderer);		
		
		xyPlot.setBackgroundPaint(ColorPastel.Gray_15.getColor());
		xyPlot.setDomainGridlinePaint(Color.white);
		xyPlot.setRangeGridlinePaint(Color.white);		
		
		LegendTitle legendTitle = theChart.getLegend();
		legendTitle.setID("Clusters");
		legendTitle.setItemFont(new Font("Arial", Font.BOLD, 20));
		
		Font rangeAxisLabelFont = new Font("Arial", Font.BOLD, 20);
		Font rangeAxisTickFont = new Font("Arial", Font.BOLD, 20);
		xyPlot.getRangeAxis().setRange(0, 1.02);
		xyPlot.getDomainAxis().setRange(0, 1.02);
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
	public static<E extends Enum<E>> void plotEventsByCoordinateAcrossSamples(EnumMapSafe<E, ParallelArrayDoubleDynamic> eventsByCoordinate, MaxPositionsOnChromosome maxPosByChrom, int numSamples, String outDir) {
		DefaultXYDataset xyDataset = new DefaultXYDataset();
		
		Set<E> enumKeys = eventsByCoordinate.keySet();		
		for (E enumKey : enumKeys) {			
			double[][] eventsXY = eventsByCoordinate.get(enumKey).toArrays();
			xyDataset.addSeries(enumKey.toString(), eventsXY);
		}
		 
		ParallelArrayDoubleDynamic chromBoundaryXY = new ParallelArrayDoubleDynamic();
		double incr = numSamples / 200.0;
		for (Chrom chrom : Chrom.Autosomes) {
			long maxPos = maxPosByChrom.getMaxPositionForChrom(chrom);
			System.out.println("MaxPos:\t" + chrom + "\t" + maxPos);
			for (double yCoordinate = 0; yCoordinate <= numSamples; yCoordinate += incr) {
				chromBoundaryXY.add(maxPos, yCoordinate);
			}			
		}
		xyDataset.addSeries("Boundary", chromBoundaryXY.toArrays());
		
		String xAxisLabel = "Position";
		String yAxisLabel = "Sample ID with Event";
		String title = "HeatMap for Events Per Sample";
		
		JFreeChart theChart = 
				ChartFactory.createScatterPlot(title, xAxisLabel, yAxisLabel, xyDataset, PlotOrientation.VERTICAL, true, false, false);
		XYPlot xyPlot = (XYPlot) theChart.getPlot();
		
		XYItemRenderer itemRenderer = ClusteringPlotting.getXYItemRendererHelper(1);		
		setSeriesPaintPerCluster(itemRenderer);
		xyPlot.setRenderer(itemRenderer);		
		
		xyPlot.setBackgroundPaint(ColorPastel.White.getColor());
		xyPlot.setDomainGridlinePaint(Color.white);
		xyPlot.setRangeGridlinePaint(Color.white);		
		
		LegendTitle legendTitle = theChart.getLegend();
		legendTitle.setID("Event");
		legendTitle.setItemFont(new Font("Arial", Font.BOLD, 20));
	
		//xyPlot.getRangeAxis().setRange(0, 1.02);		
		
		Font rangeAxisLabelFont = new Font("Arial", Font.BOLD, 20);
		Font rangeAxisTickFont = new Font("Arial", Font.BOLD, 20);
		xyPlot.getRangeAxis().setLabelFont(rangeAxisLabelFont);		
		xyPlot.getRangeAxis().setTickLabelFont(rangeAxisTickFont);
		TickUnit tickUnit =  new NumberTickUnit(1.0);
		TickUnits tickUnits = new TickUnits();
		tickUnits.add(tickUnit);		
		xyPlot.getRangeAxis().setStandardTickUnits(tickUnits);
		xyPlot.getDomainAxis().setLabelFont(rangeAxisLabelFont);
		xyPlot.getDomainAxis().setTickLabelFont(rangeAxisTickFont);	
		GraphUtils.saveChartAsPNG(outDir + File.separator + "HeapMap.AllEvents.AllSamples", theChart, 2400, 1000);


	}

	// ========================================================================
	/** Plots the event recurrence genome-wide. */
	public static void plotRecurrenceGenomeWideByEvent(EnumMapSafe<Chrom, DynamicBucketCounter> eventCount, String outDir, EventType eventType) {
		DoubleArrayList positionsGenomeWideAllSamples = new DoubleArrayList();
		DoubleArrayList eventCountAllSamples          = new DoubleArrayList();		
		
		DoubleArrayList chromBoundaryXValue = new DoubleArrayList();
		DoubleArrayList chromBoundaryYValue = new DoubleArrayList();
		PrimitiveWrapper.WInteger lastKey = new PrimitiveWrapper.WInteger(0);		
		
		// First get the max count for the chromosome boundaries
		int maxCountAllChromosomes = 1;
		for (Chrom chrom : Chrom.Autosomes) {
			maxCountAllChromosomes = Math.max(eventCount.get(chrom).getCountMax(), maxCountAllChromosomes);
		}
		
		for (Chrom chrom : Chrom.values()) {
			if (chrom.isAutosomal()) {
				long chromPosStartGenomeWide = chrom.calculateGenomeWidePositionStart();
				DoubleArrayList[] positionAndCount = eventCount.get(chrom).toArrayListDouble();
				for (int i = 0; i < positionAndCount[0].size(); i++) {
					positionAndCount[0].set(i, positionAndCount[0].get(i) + chromPosStartGenomeWide);
				}
				
				positionsGenomeWideAllSamples.addAll(positionAndCount[0]);
				eventCountAllSamples.addAll         (positionAndCount[1]);
				
				boolean lastKeyExists = eventCount.get(chrom).getKeyLast(lastKey);							
				
				double increment = (maxCountAllChromosomes) / 50.0;
				for (double d = 0; d <= maxCountAllChromosomes; d += increment) {
					chromBoundaryXValue.add(chromPosStartGenomeWide);
					chromBoundaryYValue.add(d);
				}
			}
		}
		DefaultXYDataset posAndEventCountDataset = new DefaultXYDataset();
		posAndEventCountDataset.addSeries(eventType.name(), ArrayUtils.combineTwoDynamicArraysIntoOneStatic(positionsGenomeWideAllSamples, eventCountAllSamples));
		posAndEventCountDataset.addSeries("Boundary",       ArrayUtils.combineTwoDynamicArraysIntoOneStatic(chromBoundaryXValue, chromBoundaryYValue));
		ClusteringPlotting.plotEventSampleRecurrence(posAndEventCountDataset, outDir + File.separator + "All_Samples." + eventType.name(), getColorForEvent(eventType).getColor(), ColorPastel.Gray_70.getColor());
	}
	
	// ========================================================================
	public static void plotRecurrenceGenomeWide(EnumMapSafe<EventType, EnumMapSafe<Chrom, DynamicBucketCounter>> eventCounts, String outDir) {
		for (EventType event : EventType.values()) {
			plotRecurrenceGenomeWideByEvent(eventCounts.get(event), outDir, event);
		}
	}

	// ========================================================================
	public static void plotEventSampleRecurrence(XYDataset xyDataset, String outFilenameRoot, Color seriesColor, Color boundaryColor) {
		String xAxisLabel = "Position";
		String yAxisLabel = "# Samples with Event";
		String title = "Recurrence over Samples for Event: " + xyDataset.getSeriesKey(0);				
		
		JFreeChart theChart = 
				ChartFactory.createScatterPlot(title, xAxisLabel, yAxisLabel, xyDataset, PlotOrientation.VERTICAL, true, false, false);
		XYPlot xyPlot = (XYPlot) theChart.getPlot();
		theChart.addSubtitle(getChromSubtitle());
		
		XYItemRenderer itemRenderer = ClusteringPlotting.getXYItemRendererHelper(5);		
		//setSeriesPaintPerCluster(itemRenderer);
		itemRenderer.setSeriesPaint(0, seriesColor);
		itemRenderer.setSeriesPaint(1, boundaryColor);
		xyPlot.setRenderer(itemRenderer);		
		
		xyPlot.setBackgroundPaint(ColorPastel.White.getColor());
		xyPlot.setDomainGridlinePaint(Color.white);
		xyPlot.setRangeGridlinePaint(Color.white);		
		
		LegendTitle legendTitle = theChart.getLegend();
		legendTitle.setID("Event");
		legendTitle.setItemFont(new Font("Arial", Font.BOLD, 20));
	
		//xyPlot.getRangeAxis().setRange(0, 1.02);		
		
		TickUnits tunits = new TickUnits();
		tunits.add(new NumberTickUnit(2));
		//xyPlot.getRangeAxis().setRange(0, 62);
		xyPlot.getRangeAxis().setStandardTickUnits(tunits);
		
		Font rangeAxisLabelFont = new Font("Arial", Font.BOLD, 20);
		Font rangeAxisTickFont = new Font("Arial", Font.BOLD, 20);
		xyPlot.getRangeAxis().setLabelFont(rangeAxisLabelFont);		
		xyPlot.getRangeAxis().setTickLabelFont(rangeAxisTickFont);
		xyPlot.getDomainAxis().setLabelFont(rangeAxisLabelFont);
		xyPlot.getDomainAxis().setTickLabelFont(rangeAxisTickFont);	
		GraphUtils.saveChartAsPNG(outFilenameRoot, theChart, 2400, 800);
		
	}

	// ========================================================================
	/** Plots the VAF (variant allele frequency) of a tissue vs the genomic position. */
	public static void plotCopyNumGenomeWide(XYDataset xyDataset, String outFilenameRoot, String sampleName) {
		String xAxisLabel = "Position";
		String yAxisLabel = "Copy Number";
		String title = "Copy Number: " + sampleName;				
		
		JFreeChart theChart = 
				ChartFactory.createScatterPlot(title, xAxisLabel, yAxisLabel, xyDataset, PlotOrientation.VERTICAL, true, false, false);
		XYPlot xyPlot = (XYPlot) theChart.getPlot();
		theChart.addSubtitle(getChromSubtitle());
		
		XYItemRenderer itemRenderer = ClusteringPlotting.getXYItemRendererHelper(5);		
		setSeriesPaintPerCluster(itemRenderer);
		xyPlot.setRenderer(itemRenderer);		
		
		xyPlot.setBackgroundPaint(ColorPastel.White.getColor());
		xyPlot.setDomainGridlinePaint(Color.white);
		xyPlot.setRangeGridlinePaint(Color.white);		
		
		LegendTitle legendTitle = theChart.getLegend();
		legendTitle.setID("Clusters");
		legendTitle.setItemFont(new Font("Arial", Font.BOLD, 20));
	
		//xyPlot.getRangeAxis().setRange(0, 1.02);		
		
		Font rangeAxisLabelFont = new Font("Arial", Font.BOLD, 20);
		Font rangeAxisTickFont = new Font("Arial", Font.BOLD, 20);
		xyPlot.getRangeAxis().setLabelFont(rangeAxisLabelFont);		
		xyPlot.getRangeAxis().setTickLabelFont(rangeAxisTickFont);
		xyPlot.getDomainAxis().setLabelFont(rangeAxisLabelFont);
		xyPlot.getDomainAxis().setTickLabelFont(rangeAxisTickFont);	
		GraphUtils.saveChartAsPNG(outFilenameRoot, theChart, 2400, 800);
	}

	// ========================================================================
	public static String buildChromString() {
		StringBuilder chromSB = new StringBuilder(2048);
		
		/* For font size 20 
		int[] numSpaces = new int[] {
				27,  //1 
				27,  //2 
				23,  //3 
				23,  //4
				20,  //5 
				19,  //6 
				17,  //7 
				16,  //8
				14,  //9 
				14, //10 
				15, //11
				13, //12
				10, //13 
				9,  //14 
				7,  //15 
				7,  //16
				7,  //17 
				5,  //18 
				4,  //19 
				4,  //20
				3,  //21 
				2
				
		};
		*/
		
		// For font size 24
		int[] numSpaces = new int[] {
				25,  //1 
				21,  //2 
				19,  //3 
				18,  //4
				18,  //5 
				15,  //6 
				15,  //7 
				14,  //8
				12,  //9 
				11, //10 
				11, //11
				11, //12
				8, //13 
				7,  //14 
				6,  //15 
				5,  //16
				5,  //17 
				5,  //18 
				3,  //19 
				3,  //20
				3,  //21 
				2
				
		};
		
		
		for (Chrom chrom : Chrom.Autosomes) {
			chromSB.append(chrom.ordinal());
			for (int i = 0; i < numSpaces[chrom.ordinal() - 1]; i++) {
				chromSB.append(" ");
			}
		}		
		return chromSB.toString();
	}

	// ========================================================================
	private static TextTitle getChromSubtitle() {
		TextTitle chromSubtitle = new TextTitle(ChromTitleString);
		chromSubtitle.setHorizontalAlignment(HorizontalAlignment.LEFT);
		chromSubtitle.setFont(new Font("Arial", Font.BOLD, 24));
		chromSubtitle.setMargin(10, 190, 0, 0);
		return chromSubtitle;
	}
	
	// ========================================================================
	/** Plots the VAF (variant allele frequency) of a tissue vs the genomic position. */
	public static void plotVAFGenomeWide(XYDataset xyDataset, String outFilenameRoot, String sampleName, boolean isTumor) {
		String xAxisLabel = "Position";
		String yAxisLabel = "VAF";
		String title = "VAF GenomeWide: " + sampleName + (isTumor ? " [Tumor]" : " [Normal]");
		
		//+ "\n" + ChromTitleString;
		
		JFreeChart theChart = 
				ChartFactory.createScatterPlot(title, xAxisLabel, yAxisLabel, xyDataset, PlotOrientation.VERTICAL, true, false, false);
		XYPlot xyPlot = (XYPlot) theChart.getPlot();
		theChart.addSubtitle(getChromSubtitle());
		
		XYItemRenderer itemRenderer = ClusteringPlotting.getXYItemRendererHelper(5);		
		setSeriesPaintPerCluster(itemRenderer);
		xyPlot.setRenderer(itemRenderer);		
		
		xyPlot.setBackgroundPaint(ColorPastel.White.getColor());
		xyPlot.setDomainGridlinePaint(Color.white);
		xyPlot.setRangeGridlinePaint(Color.white);
				
		Color backgroundGreyArmP = new Color(235, 235, 235, 220);
	
		//xyPlot.addAnnotation(new XYBoxAnnotation(0, 0, Chrom.c1.calculateGenomeWidePositionEndOfArmP(), 1.0, null, null, backgroundGreyArmP));
		
		for (Chrom chrom : Chrom.Autosomes) {
			xyPlot.addDomainMarker(new IntervalMarker(chrom.calculateGenomeWidePositionStart(), chrom.calculateGenomeWidePositionEndOfArmP(), backgroundGreyArmP), Layer.BACKGROUND);
			xyPlot.addDomainMarker(new ValueMarker(chrom.calculateGenomeWidePositionStart() + 1, ColorPastel.Gray_90.getColor(), new BasicStroke()));
		}
		
		LegendItemCollection legendItemCollection = xyPlot.getLegendItems();
		int numItems = legendItemCollection.getItemCount();
		for (int i = 0; i < numItems; i++) {
			LegendItem lItem = legendItemCollection.get(i);
			lItem.setShape(new Rectangle(25, 25));
		}
		xyPlot.setFixedLegendItems(legendItemCollection);
		
		LegendTitle legendTitle = theChart.getLegend();
		legendTitle.setID("Clusters");
		legendTitle.setItemFont(new Font("Arial", Font.BOLD, 32));
		
	
		xyPlot.getRangeAxis().setRange(0, 1.02);		
		
		Font rangeAxisLabelFont = new Font("Arial", Font.BOLD, 32);
		Font rangeAxisTickFont = new Font("Arial", Font.BOLD, 20);

		xyPlot.getRangeAxis().setLabelFont(rangeAxisLabelFont);		
		xyPlot.getRangeAxis().setTickLabelFont(rangeAxisTickFont);
		xyPlot.getDomainAxis().setLabelFont(rangeAxisLabelFont);
		xyPlot.getDomainAxis().setTickLabelFont(rangeAxisTickFont);	
		GraphUtils.saveChartAsPNG(outFilenameRoot, theChart, 2400, 800);
	}

	// ========================================================================
	public static XYItemRenderer getXYItemRendererHelper(int size) {
		XYShapeRenderer xyShapeRend = new XYShapeRenderer();				
		
		//XYDotRenderer xyDotRend = new XYDotRenderer();
		Ellipse2D.Double theEllipse = new Ellipse2D.Double(0, 0, size, size);
		Rectangle2D.Double theRect = new Rectangle2D.Double(0, 0, size, size);  
		//xyDotRend.setBaseShape(theEllipse, true);
		//xyShapeRend.setBaseShape(theEllipse);
		xyShapeRend.setBaseShape(theRect);
		//xyDotRend.setDotWidth(size);
		//xyDotRend.setDotHeight(size);
		//return xyDotRend;
		return xyShapeRend;
	}

	// ========================================================================
	static void setSeriesPaintPerCluster(XYItemRenderer itemRenderer) {
		
		ControlFlagBool allGray                = new ControlFlagBool(false);
		ControlFlagBool makeNonEventsInvisible = new ControlFlagBool(false);
			
		if (allGray.getValue()) {
			for (EventType eventType : EventType.values()) {
				itemRenderer.setSeriesPaint(eventType.ordinal(), ColorPastel.Gray_50.getColor());
			}
			if (makeNonEventsInvisible.getValue()) {
				itemRenderer.setSeriesPaint(EventType.HETGermline.ordinal(), ColorPastel.Gray_15.getColor());
				itemRenderer.setSeriesPaint(EventType.HETSomatic.ordinal(),  ColorPastel.Gray_15.getColor());
				itemRenderer.setSeriesPaint(EventType.Noise.ordinal(),       ColorPastel.Gray_15.getColor());
				itemRenderer.setSeriesPaint(EventType.Ignored.ordinal(),     ColorPastel.Gray_15.getColor());
			}
		} else {
			for (EventType eventType : EventType.values()) {
				itemRenderer.setSeriesPaint(eventType.ordinal(), getColorForEvent(eventType).getColor());
			}
		}
		//itemRenderer.setSeriesPaint(EventType.HETSomatic.ordinal(), getColorForEvent(EventType.HETSomatic).getColor());
		itemRenderer.setSeriesPaint(EventType.values().length, ColorPastel.Gray_70.getColor());
		itemRenderer.setSeriesShape(EventType.values().length, new Rectangle(1, 20)); 				
	}

	// ========================================================================
	public static ColorPastel getColorForEvent(EventType eventType) {
		System.out.println("EV:\t" + eventType);
		switch(eventType) {
		case GainGermline: return ColorPastel.Violet;
		case GainSomatic:  return ColorPastel.Dark_Red;
		case LOH:          return ColorPastel.RGB_Blue;
		case cnLOH:        return ColorPastel.Dark_Yellow;
		case HETGermline:  return ColorPastel.Gray_60;
		case HETSomatic:   return ColorPastel.Red_Orange;
		case DELHom:       return ColorPastel.Light_Cyan_Blue;
		case Noise:        return ColorPastel.Dark_Pea_Green;
		case Ignored:      return ColorPastel.Gray_30;			
		default:           return null;
		}
	}
}
