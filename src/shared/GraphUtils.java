package shared;

import java.io.File;
import java.io.IOException;

import org.jfree.chart.ChartUtilities;
import org.jfree.chart.JFreeChart;

public class GraphUtils {

	// ========================================================================
	/** Saves a JFreeChart as an image, given certain parameters. */
	public static void saveChartAsPNG(String outFilenameRoot, JFreeChart theChart, int chartWidth, int chartHeight) {
		try {
			//ChartUtilities.saveChartAsJPEG(theFile, theChart, chartWidth, chartHeight);
			File theFilePNG = new File(outFilenameRoot + ".png");
			ChartUtilities.saveChartAsPNG(theFilePNG, theChart, chartWidth, chartHeight);
		} catch (IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}
	}
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub

	}

}
