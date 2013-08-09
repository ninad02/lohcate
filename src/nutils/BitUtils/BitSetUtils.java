package nutils.BitUtils;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.Collection;
import java.util.Iterator;
import java.util.Random;

import nutils.CompareUtils;
import nutils.NumberUtils;

public class BitSetUtils {
	
	public static final int NumBitsInLong = Long.SIZE;
	
	public static final long[] MasksUnsignedLong = BitSetUtils.createUnsignedLongMasks();
	
	public static Random randomGenerator = new Random();

	// ========================================================================
	/** Given a number of items, this returns a bitset of length: number of items.  
	 * Given a probability, each bit is set to 1 with that given probability.
	 */
	public static BitSet getBitsetWithRandomBitsSet(int numBitsDesired, double probabilitySet) {
		BitSet bs = new BitSet(numBitsDesired);
		for (int i = 0; i < numBitsDesired; i++) {
			if (Math.random() <= probabilitySet) {
				bs.set(i);
			}
		}
		return bs;
	}

	// ========================================================================
	public static BitSet getBitsetWithRandomBitsSetExact(int numBitsDesired, double probabilitySet) {
		boolean[] theArray = new boolean[numBitsDesired];
		BitSetUtils.setBooleanArrayWithRandomValues(theArray, probabilitySet, true);				
		BitSet bs = new BitSet(numBitsDesired);
		for (int i = 0; i < numBitsDesired; i++) {
			if (theArray[i]) bs.set(i);
		}
		return bs;
	}

	// ========================================================================
	public static void fillBitsetPlusWithRandomBitsSetExact(BitSetPlus bs, int numBitsDesired, int numBitsToSet) {				
		boolean invertMode = (numBitsToSet > (numBitsDesired / 2));
		
		if (invertMode) {
			bs.set(0, numBitsDesired);
		} else {
			bs.clear();
			for (int i = 0; i < numBitsToSet; i++) { 				
				
			}
		}
	}

	// ========================================================================
	public static BitSetPlus getBitsetPlusWithRandomBitsSetExact(int numBitsDesired, double probabilitySet) {
		boolean[] theArray = new boolean[numBitsDesired];
		BitSetUtils.setBooleanArrayWithRandomValues(theArray, probabilitySet, true);				
		BitSetPlus bs = new BitSetPlus(numBitsDesired);
		for (int i = 0; i < numBitsDesired; i++) {
			if (theArray[i]) bs.set(i);
		}
		return bs;
	}

	// ========================================================================
	public static ArrayList<BitSetPlus> getListOfBitSetPlus(int numBitSetPlus, int numSamplesInOneBitSetPlus, double probabilitySet) {
		ArrayList<BitSetPlus> bsPermList = new ArrayList<BitSetPlus>(numBitSetPlus);
		for (int i = 0; i < numBitSetPlus; i++) {
			BitSetPlus bsFlip = getBitsetPlusWithRandomBitsSetExact(numSamplesInOneBitSetPlus, probabilitySet);
			bsPermList.add(bsFlip); 
			//System.out.println("Cardinality: " + bsFlip.cardinality());
			//System.out.println(bsFlip.toString());
		}
		return bsPermList;
	}

	// ========================================================================
	/** Given a boolean array, this returns the number of true values in the array. */
	public static int getCardinality(boolean[] theArray) {
		int sum = 0;
		for (int i = 0; i < theArray.length; i++) {
			if (theArray[i]) { sum++; }
		}
		return sum;
	}

	// ========================================================================
	/** Given a number of items, this returns a boolean array of length: number of 
	 *  items.  Given a probability, each element is set to true with that given probability.
	 *  Returns th enumber of true values in the array
	 */
	public static int setBooleanArrayWithRandomValues(boolean[] theArray, double probabilityTrue, boolean beStrict) {
		for (int i = 0; i < theArray.length; i++) {
			theArray[i] = (Math.random() <= probabilityTrue);
		}
		
		if (beStrict) {
			BitSetUtils.flipValuesInBooleanArray(theArray, probabilityTrue);
		}
		
		return getCardinality(theArray);
	}

	// ========================================================================
	public static int setBooleanArrayWithRandomValues(boolean[] theArray, int numTrueValuesDesired) {
		double probabilityTrue = (double) numTrueValuesDesired / (double) theArray.length;
		for (int i = 0; i < theArray.length; i++) {
			theArray[i] = (randomGenerator.nextFloat() <= probabilityTrue);
		}
		BitSetUtils.flipValuesInBooleanArray(theArray, numTrueValuesDesired);
		return getCardinality(theArray);
	}

	// ========================================================================
	public static void setBooleanArrayWithRandomValues2(boolean[] theArray, int numTrueValuesDesired) {
			boolean defaultValue = false;
			boolean newValue = true;
			int numTargetValuesDesired = numTrueValuesDesired;
			if (numTrueValuesDesired > theArray.length / 2) {
				defaultValue = true;
				newValue = false;	
				numTargetValuesDesired = theArray.length - numTrueValuesDesired;
			}
					
			int index;
			Arrays.fill(theArray, defaultValue);
			int arrayLen = theArray.length - 1;
			for (int i = 0; i < numTargetValuesDesired; i++) {
				/*
				do {
					index = CopyNumberMath.randomGenerator.nextInt(theArray.length); 
				} while (theArray[index] == newValue);
				theArray[index] = newValue;
				*/
				
				//boolean success = false;
				index = randomGenerator.nextInt(theArray.length);
				
				/*
				if (theArray[index] != newValue) {
					theArray[index] = newValue;
				} else {
					for (int j = 0; j < theArray.length; j++) {
						if ((index + j < theArray.length) && (theArray[index + j] != newValue)) {
							theArray[index + j] = newValue;
							break;
						}
						
						if ((index - j >= 0) && (theArray[index - j] != newValue)) {
							theArray[index - j] = newValue;
							break;
						}
					}
				}
				*/
							
				
				while (theArray[index] == newValue) {
					index++;
					if (index == theArray.length) { index = 0; }
				}
				if (theArray[index] != newValue) {
					theArray[index] = newValue;
				}
				
	//			if (theArray[index] != newValue) {
	//				theArray[index] = newValue;
	//			} else if (theArray[theArray.length - 1 - index] != newValue) {
	//				theArray[theArray.length - 1 - index] = newValue;
	//			} else {
	//				i--;
	//			}
			}
		}

	// ========================================================================
	public static void setBooleanArrayWithRandomValues4(boolean[] theArray, int[] costs, int numTrueValuesDesired) {
		boolean defaultValue = false;
		boolean newValue = true;
		int numTargetValuesDesired = numTrueValuesDesired;
		if (numTrueValuesDesired > theArray.length / 2) {
			defaultValue = true;
			newValue = false;	
			numTargetValuesDesired = theArray.length - numTrueValuesDesired;
		}
				
		int index;
		Arrays.fill(theArray, defaultValue);
		Arrays.fill(costs, 0);
		
		for (int i = 0; i < numTargetValuesDesired; i++) {			
			index = randomGenerator.nextInt(theArray.length);					
			if (theArray[index] != newValue) {
				theArray[index] = newValue;				
			} else {
				index = index + costs[index];
				if (index < theArray.length) {
					theArray[index] = newValue;
				} else {
					i--;
					continue;
				}
			}
			
			for (int j = index; (j >= 0) && (theArray[j] == newValue); j--) {					
				costs[j]++;
			}
		}
	}
	
	// ========================================================================
	private static int[] indicesOfExistingValuesThatCanBeFlippedFastTrueToFalse = new int[1000];
	private static int[] indicesOfExistingValuesThatCanBeFlippedFastFalseToTrue = new int[1000];
	public static void setBooleanArrayWithRandomValues3(boolean[] theArray, int numTrueValuesDesired) {
		int numTrue = 0;
		int numFalse = 0;
		double probabilityTrue = (double) numTrueValuesDesired / (double) theArray.length;
		for (int i = 0; i < theArray.length; i++) {
			theArray[i] = (randomGenerator.nextDouble() <= probabilityTrue);
			if (theArray[i]) {
				indicesOfExistingValuesThatCanBeFlippedFastTrueToFalse[numTrue] = i;
				numTrue++;
			} else {
				indicesOfExistingValuesThatCanBeFlippedFastFalseToTrue[numFalse] = i;
				numFalse++;
			}
		}		
		
		// Now determine how we did
		if (numTrue == numTrueValuesDesired) {
			return;
		} else if (numTrue > numTrueValuesDesired) {
			int diff = numTrue - numTrueValuesDesired;
			BitSetUtils.setBooleanArrayWithRandomValues3_Helper(theArray, indicesOfExistingValuesThatCanBeFlippedFastTrueToFalse, numTrue, true, diff);			
		} else {
			int diff = numTrueValuesDesired - numTrue;
			BitSetUtils.setBooleanArrayWithRandomValues3_Helper(theArray, indicesOfExistingValuesThatCanBeFlippedFastFalseToTrue, numFalse, false, diff);
		}
	}

	// ========================================================================
	public static void setBooleanArrayWithRandomValues3_Helper(boolean[] theArray, int[] valueToFlipArray, int valueToFlipArrayLength, boolean valueToFlip, int diff) {
		for (int i = diff; i > 0; i--) {
			int randomIndex = NumberUtils.getRandomInteger(0, valueToFlipArrayLength - 1);
			int trueIndex = valueToFlipArray[randomIndex];
			if (trueIndex == -1) { i++; continue; }
			CompareUtils.ensureTrue(theArray[trueIndex] == valueToFlip, "ERROR: Value is already flipped!");
			theArray[trueIndex] = !valueToFlip;			
			valueToFlipArray[randomIndex] = -1;
		}	
	}

	// ========================================================================
	/** Given a desired probability of an array element to be true, this flips the values in
	 *  the boolean array to meet the expected number of values to be true.    
	 */
	public static void flipValuesInBooleanArray(boolean[] theArray, double probabilityTrue) {
		int numExpectedTrue = (int) Math.round(theArray.length * probabilityTrue);
		int numObservedTrue = getCardinality(theArray);
		int diff = numObservedTrue - numExpectedTrue;			
		BitSetUtils.flipValuesRandomlyInBooleanArray(theArray, Math.abs(diff), (diff > 0));		
	}

	// ========================================================================
	public static void flipValuesInBooleanArray(boolean[] theArray, int numTrueValuesDesired) {
		int diff = getCardinality(theArray) - numTrueValuesDesired;			
		BitSetUtils.flipValuesRandomlyInBooleanArray(theArray, Math.abs(diff), (diff > 0));		
	}
	
	// ========================================================================
	private static int[] indicesOfExistingValuesThatCanBeFlipped = new int[1000];
	/** This flips a given number of values to flip in an array from valueToFlip -> !valueToFlip.  If there
	 *  are fewer existing instances of (valueToFlip) than numValuesToFlip, it flips all the (valueToFlip) 
	 *  instances it can find.  This method is done naively.
	 */	
	public static void flipValuesRandomlyInBooleanArray(boolean[] theArray, int numValuesToFlip, boolean valueToFlip) {
		if (numValuesToFlip == 0) return;
		numValuesToFlip = Math.abs(numValuesToFlip);
		
		int numTrueValues = getCardinality(theArray);
		int numExistingValuesThatCanBeFlipped = (valueToFlip) ? numTrueValues : theArray.length - numTrueValues;
		
		// Initialize the arraylist
		if (numExistingValuesThatCanBeFlipped > indicesOfExistingValuesThatCanBeFlipped.length) {
			indicesOfExistingValuesThatCanBeFlipped = new int[numExistingValuesThatCanBeFlipped];
		}
		int indicesUsed = 0;
		
		for (int i = 0; i < theArray.length; i++) {
			if (theArray[i] == valueToFlip) {
				indicesOfExistingValuesThatCanBeFlipped[indicesUsed] = i;
				indicesUsed++;				
			}
		}
		
		// Adjust the amount to flip and perform the flipping
		numValuesToFlip = (numExistingValuesThatCanBeFlipped < numValuesToFlip) ? numExistingValuesThatCanBeFlipped : numValuesToFlip;		
		for (int i = numValuesToFlip; i > 0; i--) {
			int randomIndex = randomGenerator.nextInt(indicesUsed); 				
			int trueIndex = indicesOfExistingValuesThatCanBeFlipped[randomIndex];
			CompareUtils.ensureTrue(theArray[trueIndex] == valueToFlip, "ERROR: Value is already flipped!");
			theArray[trueIndex] = !valueToFlip;
			indicesOfExistingValuesThatCanBeFlipped[randomIndex] = indicesOfExistingValuesThatCanBeFlipped[indicesUsed - 1];
			--indicesUsed;			
		}		
	}

	// ========================================================================
	/** This flips a given number of values to flip in an array from valueToFlip -> !valueToFlip.  If there
	 *  are fewer existing instances of (valueToFlip) than numValuesToFlip, it flips all the (valueToFlip) 
	 *  instances it can find.  This method is done naively.
	 */
	public static void flipValuesRandomlyInBooleanArrayOld(boolean[] theArray, int numValuesToFlip, boolean valueToFlip) {
		if (numValuesToFlip == 0) return;
		numValuesToFlip = Math.abs(numValuesToFlip);
		
		int numTrueValues = getCardinality(theArray);
		int numExistingValuesThatCanBeFlipped = (valueToFlip) ? numTrueValues : theArray.length - numTrueValues;
		
		// Initialize the arraylist
		ArrayList<Integer> indicesOfExistingValuesThatCanBeFlipped = new ArrayList<Integer>(numExistingValuesThatCanBeFlipped);
		for (int i = 0; i < theArray.length; i++) {
			if (theArray[i] == valueToFlip) {
				indicesOfExistingValuesThatCanBeFlipped.add(i);
			}
		}
		
		// Adjust the amount to flip and perform the flipping
		numValuesToFlip = (numExistingValuesThatCanBeFlipped < numValuesToFlip) ? numExistingValuesThatCanBeFlipped : numValuesToFlip;		
		for (int i = numValuesToFlip; i > 0; i--) {
			int randomIndex = NumberUtils.getRandomInteger(0, indicesOfExistingValuesThatCanBeFlipped.size() - 1);
			int trueIndex = indicesOfExistingValuesThatCanBeFlipped.get(randomIndex);
			CompareUtils.ensureTrue(theArray[trueIndex] == valueToFlip, "ERROR: Value is already flipped!");
			theArray[trueIndex] = !valueToFlip;
			indicesOfExistingValuesThatCanBeFlipped.remove(randomIndex);
		}		
	}

	// ========================================================================
	/** Given a boolean array, this scans the array for the true and false values.  It retains
	 *  the true values with a probability of probabilityRetain, turning them to false if they
	 *  fail this probability.  False values remain false.
	 */
	public static void retainRandomValuesInBooleanArray(boolean[] theArray, double probabilityRetain) {
		for (int i = 0; i < theArray.length; i++) {
			if (theArray[i] && (Math.random() > probabilityRetain)) {
				theArray[i] = false;
			}
		}
	}

	// ========================================================================
	/** This is a helper function that takes in a bitset queue as well as a bitset to
	 *  store the intersection of the bitsets in the queue.  This returns the bitset
	 *  passed in as an argument. 
	 */ 
	public static BitSet intersectionOfBitSets(Collection<BitSet> bitSetQueue, BitSet bsToUse) {
		bsToUse.clear();     // First, clear the bitset to use.		
		if (bitSetQueue.size() == 0) {			
			return bsToUse;
		} else {
			boolean firstRun = true;
			for (Iterator<BitSet> iter = bitSetQueue.iterator(); iter.hasNext(); ) {
				if (firstRun) {
					bsToUse.or(iter.next());  // Create a clone of the first bitset in the queue
					firstRun = false;
				} else {
					bsToUse.and(iter.next());
				}
			}
			return bsToUse;
		}		
	}

	// ========================================================================
	/** This is a helper function that takes in a bitset queue as well as a bitset to
	 *  store the union of the bitsets in the queue.  This returns the bitset passed in
	 *  as an argument.
	 */
	public static BitSet unionOfBitSets(Collection<BitSet> bitSetQueue, BitSet bsToUse) {
		bsToUse.clear();
		if (bitSetQueue.size() > 0) {
			for (Iterator<BitSet> iter = bitSetQueue.iterator(); iter.hasNext(); ) {
				bsToUse.or(iter.next());
			}			
		}	
		return bsToUse;
	}

	// ========================================================================
	/** Given a BitSet, this returns a string form of the bitset as 1s and 0s. */
	public static String getBitSetAsString(BitSet bs) {
		StringBuilder sb = new StringBuilder(bs.length());
		for (int i = 0; i < bs.length(); i++) {			
			sb.append(booleanToChar(bs.get(i)));
		}
		return sb.toString();
	}

	// ========================================================================
	public static long[] createUnsignedLongMasks() {
		long[] masks = new long[65];
		
		masks[0] = 0;
		for (int i = 1; i < masks.length; i++) {			
			masks[i] = (masks[i - 1] << 1) | 1L; 
		}
		
		return masks;
	}
	
	// ========================================================================
	/** Given a bitset and an index, this returns the previously set bit starting at that index.  If
	 *  no such bit can be found, -1 is returned. */
	public static int prevSetBit(BitSet bs, int index) {		
		if (index < 0) return -1;
		
		for (int i = index; i >= 0; i--) {
			if (bs.get(i)) return i;
		}
		
		return -1;
	}

	// ========================================================================
	public static boolean intToBoolean(int i) { return (i > 0); }
	
	// ========================================================================
	public static int booleanToInt(boolean b) { return (b ? 1 : 0); }
	
	// ========================================================================
	public static char booleanToChar(boolean b) { return (b ? '1' : '0'); }	
	
	// ========================================================================
	/**
	 * 
	 */
	public static void main(String[] args) {
		//ArrayUtils.printLongArray(MasksUns64>>ignedLong);
		//System.out.println((-1L >>> 63) >> 1);
		for (int i = 0; i <= NumBitsInLong; i++) { System.out.println(getMask(i)); }		
	}

	// ========================================================================
	/** Given a number of bits n, this returns a mask of length n up to 64 bits. */
	public static long getMask(int numBitsDesiredInMask) {		
		//return ((1L << numBitsDesiredInMask) - 1);
		return (numBitsDesiredInMask == 0 ? 0 : -1L >>> (NumBitsInLong - numBitsDesiredInMask));
	}

	// ========================================================================
	public interface ValueExtractor {
		public long extractValue(long compactUnit);		
	}

	public static BitSetUtils.ValueExtractor IntExtractorLSB = new BitSetUtils.ValueExtractor() { 
		public long extractValue(long compactUnit) { return (compactUnit & 0xFFFFFFFFL); }
	};

	public static BitSetUtils.ValueExtractor IntExtractorMSB = new BitSetUtils.ValueExtractor() { 
		public long extractValue(long compactUnit) { return ((compactUnit >>> Integer.SIZE) & 0xFFFFFFFFL); }
	};
	

}
