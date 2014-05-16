package shared;

import java.util.ArrayList;
import java.util.Collections;
import java.util.EnumSet;
import java.util.ListIterator;

import com.carrotsearch.hppc.DoubleArrayList;
import com.carrotsearch.hppc.FloatArrayList;

import nutils.CompareUtils;
import nutils.EnumMapSafe;
import nutils.IOUtils;
import nutils.NullaryClassFactory;
import nutils.NumberUtils;
import nutils.PrimitiveWrapper;
import nutils.StringUtils;

public class SortMultiCenterMutationCalling {

	private static enum Center {
		BC("bcgsc.ca"),
		Broad("broad.mit.edu"),
		UCSC("ucsc.edu"),
		HGSC("hgsc.bcm.edu"),
		MDA_HGSC("mdanderson.org")
		;
		
		String mAbbreviation;
		private Center(String abbreviation) {
			mAbbreviation = abbreviation;
		}
		
		public static Center getCenter(String abbreviation) {
			for (Center c : Center.values()) {
				if (abbreviation.equalsIgnoreCase(c.mAbbreviation)) {
					return c;
				}
			}
			return null;
		}
		
		public String getAbbreviation() { return mAbbreviation; }	
	}
	
	// ========================================================================
	public static double calcFMeasure(double recall, double precision) {
		return (2 * ((recall * precision) / (recall + precision)));
	}
	
	// ========================================================================
	public static void printOverallScores(ArrayList<FloatArrayList> overallScore, String title, int maxCardinalityToInclude) {
		
		System.out.println("-------- " + title + " Overall ----------");
		System.out.print(title);
		for (Center c : Center.values()) {
			System.out.print("\t" + c);
		}
		System.out.println("");
		ListIterator<FloatArrayList> listIter = overallScore.listIterator();
		for (int minCardinality = 1; minCardinality <= maxCardinalityToInclude; minCardinality++) {
			
			FloatArrayList theList = listIter.next();
		
			for (int j = 0; j < theList.size(); j++) {
				if (j == 0) {
					System.out.print(minCardinality);
				}
				System.out.print("\t" + theList.get(j));
			}
			System.out.println("");
		}
	}

	// ========================================================================
	private static int doScalingScore(int scalingFactor, boolean doScaling) {
		return doScaling ? scalingFactor : 1;
	}
	
	// ========================================================================
	public static void doWork(String inFilename) {
		ArrayList<String> allLines = IOUtils.readAllLinesFromFile(inFilename);		
		
		// Determine the power set
		ArrayList<EnumSet<Center>> centerPowerSet = NumberUtils.powerSetEnum(Center.class, false);
		
		// Now set the master list, broken by cardinality
		ArrayList<ArrayList<RowCount>> allRowsCardinality = new ArrayList<>();
		for (int i = 0; i < Center.values().length; i++) {
			allRowsCardinality.add(new ArrayList<RowCount>());
		}
		
		// Now iterate through the powerset and add to the appropriate cardinality row
		for (EnumSet<Center> centerSet : centerPowerSet) {
			int cardinality = centerSet.size();			
			allRowsCardinality.get(cardinality - 1).add(new RowCount(0, centerSet));
		}
		
		// Now sort the sub-lists
		for (ArrayList<RowCount> rowForCardinality : allRowsCardinality) {
			Collections.sort(rowForCardinality);
		}
		
		// Now go through the data lines
		for (String line: allLines) {
			
			String[] colsInitial = line.split(StringUtils.TabPatternStr);
			System.out.print(colsInitial[0]);
			int count = Integer.parseInt(colsInitial[0]);
			String[] theCenters = colsInitial[1].split(";");
			
			EnumSet<Center> tempSet = EnumSet.noneOf(Center.class);			
			for (String theCenterStr : theCenters) {
				Center c = Center.getCenter(theCenterStr);
				CompareUtils.ensureTrue(c != null, "Cannot be null: " + theCenterStr);
				tempSet.add(c);								
			}			
						
			RowCount rowCount = new RowCount(count, tempSet);			
			int cardinality = rowCount.cardinality();
			int resultIndex = Collections.binarySearch(allRowsCardinality.get(cardinality - 1), rowCount);
			CompareUtils.ensureTrue(resultIndex >= 0, "Must exist: " + resultIndex);
			allRowsCardinality.get(cardinality - 1).get(resultIndex).mCount = count;
		}
		System.out.println("----\nAdded All Rows\n----");
		
		for (ArrayList<RowCount> rowForCardinality : allRowsCardinality) {
			for (RowCount rc : rowForCardinality) {
				System.out.printf("%d\t%s\n", rc.mCount, rc.mSet);
			}
		}
		
		String delim = "\t";
		StringBuilder sb = new StringBuilder(1024);
		
		StringBuilder sbForR = new StringBuilder(1024);
		sbForR.append("(");
		
		StringBuilder sbForRVenneuler = new StringBuilder(1024);
		sbForRVenneuler.append("venneuler(c(");
		
		/*
		for (int i = allRowsCardinality.size() - 2; i >= 0; i--) {
			ArrayList<RowCount> rowsCurrent  = allRowsCardinality.get(i);
			ArrayList<RowCount> rowsOneLevel = allRowsCardinality.get(i + 1);
			
			for (RowCount rowCurrent : rowsCurrent) {
				for (RowCount rowOneLevel : rowsOneLevel) {
					if (rowCurrent.contains(rowOneLevel)) {
						rowCurrent.mCount += rowOneLevel.mCount;
					}
				}
			}
		}
		*/
		
		// Create an enum map of scores
		
		int maxCardinality = Center.values().length;		
		int maxCardinalityToInclude = Math.min(5, maxCardinality);
		boolean doScaling = true;

		ArrayList<FloatArrayList> scoresSensOverall = new ArrayList<FloatArrayList>();
		ArrayList<FloatArrayList> scoresPPVOverall = new ArrayList<FloatArrayList>();
		ArrayList<FloatArrayList> scoresF1Overall = new ArrayList<FloatArrayList>();
		
		for (int minCardinality = 1; minCardinality <= maxCardinalityToInclude; minCardinality++) {

			EnumMapSafe<Center, PrimitiveWrapper.WInteger> mapScorePPV = PrimitiveWrapper.WInteger.ClassFactory.newEnumMap(Center.class);
			EnumMapSafe<Center, PrimitiveWrapper.WInteger> mapTotalPPV = PrimitiveWrapper.WInteger.ClassFactory.newEnumMap(Center.class);
			EnumMapSafe<Center, PrimitiveWrapper.WInteger> mapTotalSens = PrimitiveWrapper.WInteger.ClassFactory.newEnumMap(Center.class);
			EnumMapSafe<Center, PrimitiveWrapper.WInteger> mapScoreSens = PrimitiveWrapper.WInteger.ClassFactory.newEnumMap(Center.class);

			for (int i = 0; i < allRowsCardinality.size(); i++) {
				int cardinality = i + 1;
				for (RowCount rowCount : allRowsCardinality.get(i)) {				
					for (Center c : Center.values()) {
						
						if (cardinality >= minCardinality && cardinality <= maxCardinalityToInclude) {
							PrimitiveWrapper.WInteger totalSens = mapTotalSens.get(c);
							totalSens.mInt += (doScalingScore(cardinality, doScaling) * rowCount.mCount);
						}

						if (rowCount.contains(c)) {
							PrimitiveWrapper.WInteger totalCount = mapTotalPPV.get(c);
							totalCount.mInt += (doScalingScore(maxCardinalityToInclude, doScaling) * rowCount.mCount);

							if (cardinality >= minCardinality && cardinality <= maxCardinalityToInclude) {
								
								totalCount = mapScorePPV.get(c);
								totalCount.mInt += (doScalingScore(cardinality, doScaling) * rowCount.mCount);
								
								totalCount = mapScoreSens.get(c);
								totalCount.mInt += (doScalingScore(cardinality, doScaling) * rowCount.mCount);
							}
						}
					}
				}
			}

			// Now print the scores
			FloatArrayList scoresPPV = new FloatArrayList(maxCardinality);
			scoresPPVOverall.add(scoresPPV);
			System.out.println("\n---------Prediction Scores---------");
			System.out.println("Center\tTotal Calls\tPrediction Score\tMax Predication Score\tPrediction Score Fraction");
			for (Center c : Center.values()) {
				int score = mapScorePPV.get(c).mInt;
				int total = mapTotalPPV.get(c).mInt;
				float fraction = (float) score / (float) total;
				scoresPPV.add(fraction);
				System.out.printf("%s\t%d\t%d\t%d\t%f\t\n", c, total / maxCardinality, score, total, fraction);			
			}

			FloatArrayList scoresSens = new FloatArrayList(maxCardinality);
			scoresSensOverall.add(scoresSens);
			System.out.println("\n---------Prediction Scores---------");
			System.out.println("Center\tCall Score\tTotal Call Score\tCall Score Fraction");
			for (Center c : Center.values()) {
				int score = mapScoreSens.get(c).mInt;
				int total = mapTotalSens.get(c).mInt;
				float fraction = (float) score / (float) total;
				scoresSens.add(fraction);
				System.out.printf("%s\t%d\t%d\t%f\t\n", c, score, total, fraction);
			}

			FloatArrayList scoresF = new FloatArrayList(maxCardinality);
			scoresF1Overall.add(scoresF);
			System.out.println("\n---------F Scores---------");
			int index = 0;
			for (Center c : Center.values()) {
				double fMeasure = calcFMeasure(scoresSens.get(index), scoresPPV.get(index));
				scoresF.add((float) fMeasure);
				System.out.printf("%s\t%g\n", c, fMeasure);
				index++;
			}
			System.out.println("\n------------------");


			StringBuilder alphaSB = new StringBuilder(1024);
			alphaSB.append("alpha=c(");
			for (int i = 0; i < scoresPPV.size(); i++) {
				if (i > 0) alphaSB.append(", ");
				alphaSB.append(scoresPPV.get(i));
			}
			alphaSB.append(")");
			System.out.println(alphaSB.toString());

			System.out.println("\n------------------");

			System.out.println("");
		}
		
		// Print matrices
		printOverallScores(scoresSensOverall, "Sens_Score", maxCardinalityToInclude);
		printOverallScores(scoresPPVOverall, "PPV_Score", maxCardinalityToInclude);
		printOverallScores(scoresF1Overall, "F1_Score", maxCardinalityToInclude);
		
		
		for (int i = 0; i < allRowsCardinality.size() - 1; i++) {
			ArrayList<RowCount> rowsCurrent  = allRowsCardinality.get(i);			
			for (int j = i + 1; j < allRowsCardinality.size(); j++) {
				ArrayList<RowCount> rowsFutureLevel = allRowsCardinality.get(j);
				
				for (RowCount rowCurrent : rowsCurrent) {
					for (RowCount rowFutureLevel : rowsFutureLevel) {
						if (rowCurrent.contains(rowFutureLevel)) {
							rowCurrent.mCount += rowFutureLevel.mCount;
						}
					}
				}
			}
		}
		

		
		
		for (ArrayList<RowCount> allRows : allRowsCardinality) {
			Collections.sort(allRows);
			
			for (RowCount rc : allRows) {
				sb.setLength(0);
				sb.append(rc.mCount);
				for (Center c : Center.values()) {
					if (rc.contains(c)) {
						sb.append(delim).append(c);
					}
				}
				sb.append(delim).append(rc.cardinalityString("&"));
				System.out.println(sb.toString());
				if (sbForR.length() > 1) {
					sbForR.append(", ");
				}
				sbForR.append(rc.mCount);
				sbForRVenneuler.append("\"").append(rc.cardinalityString("&")).append("\"").append("=").append(rc.mCount).append(", ");
			}
		}
		
		sbForR.append(", c(");
		int centerIndex = 0;
		for (Center c : Center.values()) {
			sbForR.append(centerIndex++ > 0 ? ", " : "").append("\"").append(c.name()).append("\"");
		}
		sbForR.append(")");		
		sbForR.append(", cat.dist=c(.08,.08,.08,.08,.08)");
		sbForR.append(", fill=c(\"blue\", \"red\", \"green\", \"purple\", \"orange\")");
		sbForR.append(", alpha=c(0.5, 0.5, 0.5, 0.5, 0.5)");
		
		sbForR.append(")");
		sbForR.insert(0, "draw.quintuple.venn");
		sbForRVenneuler.append("))");
		System.out.println(sbForR.toString());
		System.out.println(sbForRVenneuler);
	}
	
	
	// ========================================================================
	private static class RowCount implements Comparable<RowCount> { 
		EnumSet<Center> mSet;		
		int mCount;
		
		public RowCount(int count, EnumSet<Center> theSet) {
			mSet = theSet;
			mCount = count;			
		}
		
		public void setCenters(ArrayList<Center> theCenters) { mSet.addAll(theCenters); }
		
		public int cardinality() { return mSet.size(); }
		
		/** Returns true if this object contains the parameter-argument. */
		public boolean contains(RowCount rhs) {
			// String rowThis = cardinalityString("&");
			// String rowRhs  = rhs.cardinalityString("&");
			// System.out.println("\n----\n" + rowThis + "\n" + rowRhs + "\n----");
			int matchCount = 0;
			int mismatchCount = 0;
			
			for (Center c : Center.values()) {
				boolean resultThis = this.mSet.contains(c);
				boolean resultRhs  =  rhs.mSet.contains(c);				
				if (resultThis && !resultRhs) {					
					//System.out.println("FALSE");
					return false;
				} else if (resultThis && resultRhs) {
					matchCount++;
				} else if (!resultThis && resultRhs) {
					mismatchCount++;
				}
			}
			if (matchCount == cardinality() && (mismatchCount >= 1)) {
				//System.out.println("TRUE");
				return true;
			} else return false;
		}
		
		public boolean contains(Center c) { return mSet.contains(c); }
		
		public String cardinalityString(String delimiter) {
			StringBuilder sb = new StringBuilder(1024);
			for (Center c : Center.values()) {
				if (mSet.contains(c)) {				
					if (sb.length() > 0) {
						sb.append(delimiter);
					}
					sb.append(c);
				}				
			}
			return sb.toString();
		}
		
		public void clear() {
			mCount = 0;
		}
		
		public int compareTo(RowCount rhs) {
			for (Center c : Center.values()) {
				boolean resultThis = this.mSet.contains(c);
				boolean resultRhs  =  rhs.mSet.contains(c);
				int result = Boolean.compare(resultThis, resultRhs);
				if (result != 0) return -result;
			}
			return 0;
		}
	}
	
	public static void main(String[] args) {
		doWork(args[0]);
	}
}
