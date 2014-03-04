package nutils.array;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;

import nutils.ArrayUtils;
import nutils.CompareUtils;
import nutils.PrimitiveWrapper;
import nutils.BitUtils.BitSetUtils;

import com.carrotsearch.hppc.LongArrayList;

/**
 * A class to store primitive long variables with fast addition and deletion.
 * Very useful when storing and modifying lists of large size.
 * 
 * @author Ninad Dewal
 *
 */

public class FastArrayListLong {
	
	private static final int SubListMaxLength = 200;
	private static final int FirstSubListIndex = 0;
	
	/***********************
	 * Member Variables
	 ***********************/
	
	protected ArrayList<SubList> mTopLevel;
	protected int mNumElements;
	protected PrimitiveWrapper.WInteger mSublistIndex;
	protected PrimitiveWrapper.WBoolean mWasFound;
	protected SubList mDummySubList;
	protected SortedSubListComparator mComparator;
	
	// ========================================================================
	public FastArrayListLong() {
		mTopLevel = new ArrayList<>();
		mTopLevel.add(new SubList());
		clear();
		mSublistIndex = new PrimitiveWrapper.WInteger(0);
		mDummySubList = new SubList();
		mComparator = new SortedSubListComparator();
	}
	
	// ========================================================================
	public void clear() {
		for (SubList longList : mTopLevel) {
			longList.clear();
		}
		mNumElements = 0;		
	}

	// ========================================================================
	public void setValueExtractor(BitSetUtils.ValueExtractor valueExtractor) {
		mComparator.setValueExtractor(valueExtractor);
	}

	// ========================================================================
	public FastArrayListIndex getIndex(long value, FastArrayListIndex returnIndex) {
		//printDebugAllSubLists();
		int indexSubList = getSubListIndex(value);
		if (indexSubList < 0) {
			return null;
		}
		
		SubList subList = mTopLevel.get(indexSubList);
		int indexWithinSubList = subList.getIndexOfValue(value, mComparator.mValueExtractor); 				
		if (indexWithinSubList < 0) {
			return null;
		}
		
		if (returnIndex == null) {
			returnIndex = new FastArrayListIndex();
		}
		returnIndex.set(indexSubList, indexWithinSubList, subList.mList.get(indexWithinSubList));
		return returnIndex;
	}

	// ========================================================================
	public long get(FastArrayListIndex index) {
		return mTopLevel.get(index.mIndexSubList).mList.get(index.mIndexWithinSubList);
	}
	
	// ========================================================================
	public FastArrayListIterator iterator() { return new FastArrayListIterator(this); }
	
	// ========================================================================
	public boolean contains(long value) {
		return (getIndex(value, null) != null);
	}

	// ========================================================================
	public boolean replace(FastArrayListIndex arrayIndex, long value) {
		LongArrayList targetList = mTopLevel.get(arrayIndex.mIndexSubList).mList;
		
		if (targetList.isEmpty()) {
			CompareUtils.ensureTrue(false, "ERROR: Target List changed!");
		}
		
		int nextIndex = arrayIndex.mIndexWithinSubList + 1;
		int prevIndex = arrayIndex.mIndexWithinSubList - 1;
		
		if (arrayIndex.mIndexWithinSubList == 0) {			
			if (value > targetList.get(nextIndex)) return false;  // return because value to be set will break sorted state	
			
		} else if (arrayIndex.mIndexWithinSubList == targetList.size() - 1) {
			if (value < targetList.get(prevIndex)) return false;  // return because value to be set will break sorted state
				
		} else {
			if ((value > targetList.get(nextIndex)) || (value < targetList.get(prevIndex))) return false; 
		}
		
		// Now we are guaranteed that we will not break the ordering
		targetList.set(arrayIndex.mIndexWithinSubList, value);
		return true;
	}
	
	// ========================================================================
	public void add(long value) {		
		int subListIndex = getSubListIndex(value);
		subListIndex = ArrayUtils.getInsertPoint(subListIndex);
		
		// Now add to the sublist		
		SubList subListCurr = mTopLevel.get(subListIndex);
		SubList potientialNewSubList = subListCurr.add(value, mComparator.mValueExtractor);
		if (potientialNewSubList != null) {
			mTopLevel.add(subListIndex + 1, potientialNewSubList);
		}
	}
	
	// ========================================================================
	protected int getSubListIndex(final long targetValue) {		
		
		// First see whether we have an empty list
		if (mTopLevel.get(FirstSubListIndex).mList.isEmpty()) {
			return FirstSubListIndex;
		} 
		
		// Initialize the dummy
		mDummySubList.mList.clear();
		mDummySubList.mList.add(targetValue);
		
		int subListResultIndex = Collections.binarySearch(mTopLevel, mDummySubList, mComparator);
		if (subListResultIndex >= 0) {
			// A suitable sub-list was found.
			return subListResultIndex;
		} else {
			int newIndex = ArrayUtils.getInsertPoint(subListResultIndex);
			if (newIndex == 0) {
				// Means that the value needs to be inserted in the first sublist.
				return subListResultIndex;
			} else if (newIndex == mTopLevel.size()) {
				// Means that it's after the last element in the last sublist.
				// We still return the index of the last sublist so that the sublist
				// can later be split if the element is to be added.
				int lastSubListIndex = mTopLevel.size() - 1;
				return ArrayUtils.getFailPoint(lastSubListIndex); 						
				//return mTopLevel.get(lastSubListIndex).isFull() ? newIndex : lastSubListIndex; 
			} else {
				// At this point, newIndex is guaranteed to point to sublist 1 ... last and not sublist 0
				// We see whether the previous or the current sublist has more room and insert in that one.				
				int prevIndex = newIndex - 1;					
				
				return ((mTopLevel.get(prevIndex).mList.size() > mTopLevel.get(newIndex).mList.size()) ?
						ArrayUtils.getFailPoint(newIndex) :  ArrayUtils.getFailPoint(prevIndex));
			}
		}		
	}
	
	
	// ========================================================================	
	protected static class SortedSubListComparator implements Comparator<SubList> {

		BitSetUtils.ValueExtractor mValueExtractor = BitSetUtils.LongExtractorWhole;

		public void setValueExtractor(final BitSetUtils.ValueExtractor valueExtractor) {
			mValueExtractor = valueExtractor;
		}
		
		@Override
		public int compare(SubList list1, SubList list2) {
			
			if (list1.mList.size() == 1) {
				
				if (list2.mList.isEmpty()) {
					CompareUtils.ensureTrue(false, "ERROR: Impossible state.  Sublist 2 cannot be empty!");
				}
				
				final long list2FirstElem = mValueExtractor.extractValue(list2.getFirstElem()); 
				final long list2LastElem  = mValueExtractor.extractValue(list2.getLastElem());		
				
				final long list1Elem = mValueExtractor.extractValue(list1.getFirstElem());
				
				if (list1Elem < list2FirstElem) {
					return -1;
				} else if (list1Elem > list2LastElem) {
					return 1;
				} else {
					return 0;
				}
				
			} else if (list2.mList.size() == 1) {
				return (-compare(list2, list1)); // reverse
				
			} else {
				CompareUtils.ensureTrue(false, "ERROR: Impossible state.  Both compared lists have != 1 element!");
			}
			
			return 0;
		}
		
	}

	// ========================================================================
	public static class FastArrayListIterator {
		private FastArrayListIndex mIndex;
		private FastArrayListLong mList;
		
		protected FastArrayListIterator(FastArrayListLong theList) {
			mIndex = new FastArrayListIndex();
			mList = theList;
			reset();
		}
		
		public void reset() {
			mIndex.mIndexSubList = 0;
			mIndex.mIndexWithinSubList = -1;
		}
		
		/** Checks if there is a subsequent value in this iteration.  If so, this returns true,
		 *  iterates to the next element, and writes the value in the primitive-holder argument.
		 *  If not, this returns false and does not modify the argument.  If the passed in argument 
		 *  is null, this method only returns true or false without iterating to the next element. 
		 * @param value -- A primitive holder that will contain the next value, if applicable.  If 
		 *  set to null, no value will be returned, though true or false still will be. 
		 * @return
		 */
		public boolean hasNextAndGet(PrimitiveWrapper.WLong value) {
			
			// Check whether we have a valid index
			if (mIndex.mIndexSubList >= mList.mTopLevel.size()) {
				CompareUtils.ensureTrue(false, "ERROR: Impossible State!  Cannot be residing in invalid sublist!");
			}
			
			// We're in a valid sublist
			SubList subList = mList.mTopLevel.get(mIndex.mIndexSubList);
			
			int nextIndex = mIndex.mIndexWithinSubList + 1;
			if (nextIndex < subList.mList.size()) {

				// If not null, we write it in 
				if (value != null) {
					mIndex.mIndexWithinSubList = nextIndex;
					value.mLong = mIndex.mValue = mList.get(mIndex);
				}
				
				return true;
				
			} else {
				// We've reached the current sublist's max # elements.
				// We test whether there is another sublist
				for (int nextSubListIndex = mIndex.mIndexSubList + 1; nextSubListIndex < mList.mTopLevel.size(); nextSubListIndex++) {
					SubList subListNext = mList.mTopLevel.get(nextSubListIndex);
					if (!subListNext.mList.isEmpty()) {
						if (value != null) {
							mIndex.mIndexSubList = nextSubListIndex;
							mIndex.mIndexWithinSubList = 0;
							value.mLong = mIndex.mValue = mList.get(mIndex);
						}
						return true;
					}
				}

				return false;
			}
		}		
	}
	
	// ========================================================================
	public static class FastArrayListIndex {
		private int mIndexSubList;
		private int mIndexWithinSubList;
		public long mValue;
		
		protected FastArrayListIndex() {
			mIndexSubList = -1;
			mIndexWithinSubList = -1;
		}
		
		protected void set(int indexSubList, int indexWithinSubList, long value) {
			mIndexSubList = indexSubList;
			mIndexWithinSubList = indexWithinSubList;
			mValue = value;
		}
		
		public long getValue() { return mValue; }
	}
	
	// ========================================================================
	public static FastArrayListIndex getNewIndexToken() { return new FastArrayListIndex(); }
	
	// ========================================================================
	private void printDebugAllSubLists() {
		int index = 0;
		for (SubList subList : mTopLevel) {		
			++index;
			System.out.printf("SubList\t%d\t%d\t%d\t%d\n", index, subList.mList.size(), subList.getFirstElem(), subList.getLastElem());
		}
	}

	// ========================================================================
	// INNER CLASS: SubList
	// ========================================================================
	protected static class SubList {
		LongArrayList mList;

		// ====================================================================
		public SubList() {
			mList = new LongArrayList(SubListMaxLength + 1);
		}
		
		// ====================================================================
		public boolean isFull() { return mList.size() >= SubListMaxLength; }

		// ====================================================================
		public void clear() { mList.clear(); }

		// ====================================================================
		public long getFirstElem() { return mList.get(0); }
		
		// ====================================================================
		public long getLastElem() { return mList.get(mList.size() - 1); }

		// ====================================================================
		public int getIndexOfValue(long value, BitSetUtils.ValueExtractor valueExtractor) {
			value = valueExtractor.extractValue(value);
			return ArrayUtils.binarySearchValue(value, mList, valueExtractor);
		}
		
		// ====================================================================
		/** Adds a value.  If adding the value exceeds the max list size, this 
		 * splits the SubList and returns the new sublist created, that follows
		 * this one.  
		 * @param value
		 * @return A new sublist (which follows this one), if one is to be created; else @null is returned otherwise. 
		 */
		public SubList add(long value, BitSetUtils.ValueExtractor valueExtractor) {
			boolean listIsAlreadyFull = isFull();
			
			int resultIndex = getIndexOfValue(value, valueExtractor);
			int newIndex = ArrayUtils.getInsertPoint(resultIndex);
			mList.insert(newIndex, value);
			
			// Now if we've exceeded the max size, we split.
			if (listIsAlreadyFull) {
				SubList newSubList = new SubList();
				int firstIndexToTransfer = mList.size() >>> 1; // Divide by 2				
				int numElementsToTransfer = mList.size() - firstIndexToTransfer;
				
				// Transfer the elements
				newSubList.mList.add(this.mList.buffer, firstIndexToTransfer, numElementsToTransfer);
				
				// Now remove the elements from this sublist
				this.mList.removeRange(firstIndexToTransfer, mList.size());
				return newSubList;
			}
			
			return null;
		}
	}

	// ============================================================================
	private static void Test2() {
		FastArrayListLong theList = new FastArrayListLong();
		long memoryStart = Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory();
		long timeStart = System.currentTimeMillis();
		
		int target = 2_000_000_000;
		int incr = 20;
		for (int even = 0; even < target; even += incr) {
			theList.add(even);
		}
		
		for (int even = 10; even < target; even += incr) {
			theList.add(even);
		}
		
		int newOffset = 4;
		FastArrayListIndex indexToken = FastArrayListLong.getNewIndexToken();
		for (int even = 0; even < target; even += (incr >>> 1)) {
			theList.getIndex(even, indexToken);
			boolean result = theList.replace(indexToken, indexToken.mValue + newOffset);
			if (!result) {
				System.out.printf("Could not replace %d with %d\n", even, indexToken.mValue + newOffset);
			}
		}
		
		for (int even = 0; even < target; even += (incr >>> 1)) {
			if (!theList.contains(even + newOffset)) {
				CompareUtils.ensureTrue(false, "Must contain: " + (even + newOffset));
			}	
		}
		
		long timeEnd = System.currentTimeMillis();
		long memoryEnd = Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory();
		System.out.println("Success! " + ((timeEnd - timeStart) / 1000.0) + "\t" + ((memoryEnd - memoryStart) / 1048576.0));
	}
	
	// ============================================================================
	private static void Test1() {
		FastArrayListLong theList = new FastArrayListLong();
		
		long memoryStart = Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory();
		long timeStart = System.currentTimeMillis();
		int target = 200_000_000;
		for (int even = 0; even < target; even += 2) {
			theList.add(even);		
			//System.out.println(even);
		}

		for (int odd = 1; odd < target; odd += 2) {
			theList.add(odd);		
			//System.out.println(even);
		}

		for (int i = 0; i < target; i++) {
			if (!theList.contains(i)) {
				CompareUtils.ensureTrue(false, "Must contain: " + i);
			}
		}
		long timeEnd = System.currentTimeMillis();
		/*
		for (int even = 0; even < target; even += 2) {
			//System.out.println("Trying: " + even);
			if (!theList.contains(even)) {
				CompareUtils.ensureTrue(false, "Must contain: " + even);
			} else if (theList.contains(even + 1)) {
				CompareUtils.ensureTrue(false, "Must not contain: " + (even + 1));
			}
		}
		*/
		long memoryEnd = Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory();
		System.out.println("Success! " + (timeEnd - timeStart) + "\t" + (memoryStart - memoryEnd));
	}
	
	// ============================================================================
	public static void main(String[] args) {
		Test2();
	}
	
}

