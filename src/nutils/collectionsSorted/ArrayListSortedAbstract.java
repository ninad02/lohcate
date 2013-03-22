package nutils.collectionsSorted;

import java.lang.reflect.Array;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Comparator;
import java.util.Iterator;
import java.util.List;
import java.util.ListIterator;
import java.util.NoSuchElementException;

public abstract class ArrayListSortedAbstract<E> extends ArrayList<E> {

	/**
	 * 
	 */
	private static final long serialVersionUID = 5176003946276199793L;
	
	protected static int DefaultNumInitialCapacity = 10;
	
	// ========================================================================
	// MEMBER VARIABLES
	protected boolean       mIsSorted; 
	protected Comparator<E> mComparator;
	
	// ========================================================================
	// Constructors
	/** Constructs an empty list with an initial capacity of ten. */
	protected ArrayListSortedAbstract() { this(DefaultNumInitialCapacity, null); }

	// ========================================================================
	/** Constructs an empty list with an initial capacity of ten and with the given comparator. */
	protected ArrayListSortedAbstract(Comparator<E> comparator) { this(DefaultNumInitialCapacity, comparator); }
	
	// ========================================================================
	/** Constructs an empty list with the specified initial capacity. */
	protected ArrayListSortedAbstract(int initialCapacity) { this(initialCapacity, null); }

	// ========================================================================
	/** Constructs an empty list with the specified initial capacity and a comparator. */
	protected ArrayListSortedAbstract(int initialCapacity, Comparator<E> comparator) {
		super(initialCapacity);		
		mComparator = comparator;
		mIsSorted = true;
	}
	
	// ========================================================================
	/** Constructs a list containing the elements of the specified collection, in the order they are returned by the collection's iterator. */
	protected ArrayListSortedAbstract(Collection<? extends E> c, Comparator<E> comparator) {
		super(c);
		mComparator = comparator;
		mIsSorted = false;
		ifNotSortedThenSort(true);
	}
	
	// ========================================================================
	protected ArrayListSortedAbstract(ArrayListSortedAbstract<E> rhs) {
		super(rhs);
		mComparator = rhs.mComparator;
		mIsSorted   = rhs.mIsSorted;
	}

	// ========================================================================
	/** Sets the comparator and re-sorts the list if the comparator argument differs from the one already stored internally. 
	 * @return this object, which allows for compounded statements: setComparator(c).get(0) */
	public ArrayListSortedAbstract<E> setComparator(final Comparator<E> comparator) {
		if (mComparator != comparator) {
			mComparator = comparator;
			ifNotSortedThenSort(true);
		}
		return this;
	}
	
	// ========================================================================
	/** Returns the index for a specified element. */
	public abstract int getIndex(final E element);
	
	// ========================================================================
	/** Returns the object at the specified index.  Ensures the array is sorted beforehand. */
	public E get(int index) {
		ifNotSortedThenSort();
		return super.get(index);		
	}
	
	// ========================================================================
	protected E getUnsorted(int index) { return super.get(index); }
		
	// ========================================================================
	/** Returns the object that is stored in this array. 
	 * @return The element if it exists, or null if not. */	
	public E get(E element) {
		int resultIndex = getIndex(element);
		return (resultIndex < 0) ? null : super.get(resultIndex);
	}	

	// ========================================================================
	/** Returns the index of the last element. */
	public int getLastIndex() { return size() - 1; }
	
	// ========================================================================
	/** Removes all of the elements from this list. */
	public void clear() {
		super.clear();
		mIsSorted = true;
	}
	
	// ========================================================================
	/** Adds the specified element to the list. */
	public boolean add(E element) {
		this.add(size(), element);
		return true;
	}
	
	// =======================================================================-
	/** Inserts the specified element at the specified position in this list. */
	@Override
	public void add(int index, E element) {
		super.add(index, element);
		if (!wasInsertedInSortedManner(index)) {
			mIsSorted = false;
			System.out.println("Not Sorted: " + index + "\t" + element);
		}
	}
	
	// ========================================================================
	/** Adds all of the elements in the specified collection to the list, in the order that they are returned by the specified collection's Iterator.
	 */
	@Override
	public boolean addAll(Collection<? extends E> c) {
		for (E element : c) {
			this.add(element);
		}
		return true;
	}
	
	// ========================================================================
	/** This adds all of the elements in the specified collection into this list.  It ignores the specified position,
	 *  as the elements are to be sorted.  This provides identical functionality as the other addAll() method.  
	 */
	@Override
	public boolean addAll(int index, Collection<? extends E> c) {
		return this.addAll(c);		
	}

	
	// ========================================================================
	/** Replaces the element at the specified position in this list with the specified element. */
	@Override
	public E set(int index, E element) {
		E rV = get(index);
		super.set(index, element);
		if (!wasInsertedInSortedManner(index)) {
			mIsSorted = false;
		}
		return rV;		
	}
	
	// ========================================================================
	protected E setUnsorted(int index, E element) { return super.set(index, element); }
	
	// ========================================================================
	/** Returns a shallow copy of this ArrayList instance. */
	@Override
	public abstract Object clone();

	// ========================================================================
	/** Returns the index of the first occurrence of the specified element in this list, or -1 if this list does not contain the element. */
	@Override
	public int indexOf(Object o) {
		ifNotSortedThenSort();
		return super.indexOf(o);
	}

	// ========================================================================
	/** Returns an iterator over the elements in this list in proper sequence. */
	@Override
	public Iterator<E> iterator() {
		return new ArrayListSortedIterator<E>(this, 0);
	}

	// ========================================================================
	/** Returns the index of the last occurrence of the specified element in this list, or -1 if this list does not contain the element. */
	@Override
	public int lastIndexOf(Object o) {
		ifNotSortedThenSort();
		return super.lastIndexOf(o);
	}

	// ========================================================================
	/** Returns a list iterator over the elements in this list (in proper sequence). */
	@Override
	public ListIterator<E> listIterator() {
		return new ArrayListSortedIterator<E>(this, 0);
	}

	// ========================================================================
	/** Returns a list iterator over the elements in this list (in proper sequence), starting at the specified position in the list. */
	@Override
	public ListIterator<E> listIterator(int index) {
		return new ArrayListSortedIterator<E>(this, index);
	}

	// ========================================================================
	/** Removes the element at the specified position in this list. */
	@Override
	public E remove(int index) {
		ifNotSortedThenSort();
		return super.remove(index);
	}

	// ========================================================================
	/** Removes the first occurrence of the specified element from this list, if it is present. */
	@Override
	public boolean remove(Object o) {
		ifNotSortedThenSort();
		return super.remove(o);
	}
	
	// ========================================================================
	/* (non-Javadoc)
	 * @see java.util.AbstractCollection#containsAll(java.util.Collection)
	 */
	@Override
	public boolean containsAll(Collection<?> c) {
		ifNotSortedThenSort();
		return super.containsAll(c);
	}

	// ========================================================================
	/** Returns a string representation of this collection. The string representation consists of a list of the collection's 
	 *  elements in the order they are returned by its iterator, enclosed in square brackets ("[]"). Adjacent elements are 
	 *  separated by the characters ", " (comma and space). Elements are converted to strings as by String.valueOf(Object). 
	 */
	@Override
	public String toString() {
		ifNotSortedThenSort();
		return super.toString();
	}

	// ========================================================================
	/** Removes from this list all of its elements that are contained in the specified collection. */
	@Override
	public boolean removeAll(Collection<?> c) {
		ifNotSortedThenSort();
		return super.removeAll(c);
	}

	// ========================================================================
	/**	Removes from this list all of the elements whose index is between fromIndex, inclusive, and toIndex, exclusive. */
	@Override
	protected void removeRange(int fromIndex, int toIndex) {
		ifNotSortedThenSort();
		super.removeRange(fromIndex, toIndex);
	}

	// ========================================================================
	/** Returns a view of the portion of this list between the specified fromIndex, inclusive, and toIndex, exclusive. */
	/* (non-Javadoc)
	 * @see java.util.ArrayList#subList(int, int)
	 */
	@Override
	public List<E> subList(int fromIndex, int toIndex) {
		ifNotSortedThenSort();
		return super.subList(fromIndex, toIndex);
	}
	
	// ========================================================================
	/** Returns an array containing all of the elements in this list in proper sequence (from first to last element). 
	 * @override */
	@Override
	public Object[]	toArray() {
		ifNotSortedThenSort();
		return super.toArray();
	}

	// ========================================================================
	/** Returns an array containing all of the elements in this list in proper sequence (from first to last element); the runtime type of the returned array is that of the specified array. */
	@Override
	public<T> T[] toArray(T[] a) {
		ifNotSortedThenSort();
		return super.toArray(a);
	}
	
	// ========================================================================
	
	// ========================================================================
	protected boolean wasInsertedInSortedManner(final int indexOfElement) {
		if (size() <= 1) return true;  // return true for empty or singleton lists
		
		E element = super.get(indexOfElement);
		
		if (indexOfElement == getLastIndex()) {			
			E elementPrev = super.get(indexOfElement - 1);
			return (compareHelper(element, elementPrev) >= 0);
			
		} else if (indexOfElement == 0) {
			E elementNext = super.get(indexOfElement + 1);
			return (compareHelper(element, elementNext) <= 0);
			
		} else {
			E elementPrev = super.get(indexOfElement - 1);
			E elementNext = super.get(indexOfElement + 1);		
			return (compareHelper(element, elementPrev) >= 0) && (compareHelper(element, elementNext) <= 0); 
		}
	}
	
	// ========================================================================
	protected void ifNotSortedThenSort() { 
		ifNotSortedThenSort(false);
	}
	
	// ========================================================================	

	// ========================================================================
	protected abstract int compareHelper(E elementCurrent, E elementPrevOrNext);

	// ========================================================================
	public void sortArray(E[] theArray) {
		if (mComparator == null) {
			Arrays.sort(theArray);			
		} else {				
			Arrays.sort(theArray, mComparator);
		}		
	}
		
	// ========================================================================
	protected void ifNotSortedThenSort(boolean forceSort) {		
		if ((forceSort || !mIsSorted) && (size() > 1)) {
			//System.out.println("Sorting");
			E firstElement = getUnsorted(0);
			E[] newArray = (E[]) Array.newInstance(firstElement.getClass(), size());
			
			for (int index = 0; index < size(); index++) {
				newArray[index] = getUnsorted(index);
			}
			
			sortArray(newArray);
			
			for (int index = 0; index < size(); index++) {
				setUnsorted(index, newArray[index]);
			}						
		} else {
			//System.out.println("Not sorting");
		}
		mIsSorted = true;
	}

	
	// ========================================================================
	// INNER CLASS
	// ========================================================================
	public static class ArrayListSortedIterator<E> implements ListIterator<E> {

		int mCurrentIndex = -1;
		ArrayListSortedAbstract<E> mList;
		
		// ========================================================================
		
		protected ArrayListSortedIterator(ArrayListSortedAbstract<E> theList, int startingIndex) {
			mList = theList;
			mCurrentIndex = startingIndex - 1;
		}
		
		// ========================================================================
		/** Inserts the specified element into the list (optional operation). */
		@Override
		public void add(E e) {
			mList.add(e);
		}

		// ========================================================================
		/** Returns true if this list iterator has more elements when traversing the list in the forward direction. */	
		@Override
		public boolean hasNext() {
			return (mCurrentIndex < mList.size() - 1);			
		}

		// ========================================================================
		/** Returns true if this list iterator has more elements when traversing the list in the reverse direction. */
		@Override
		public boolean hasPrevious() {			
			return (mCurrentIndex > 0);
		}

		// ========================================================================
		
		/** Returns the next element in the list and advances the cursor position. */		
		/** Returns the index of the element that would be returned by a subsequent call to next(). */		
		/** Returns the previous element in the list and moves the cursor position backwards. */

		@Override
		public E next() {
			mList.ifNotSortedThenSort();
			++mCurrentIndex;
			
			if (mCurrentIndex < mList.size()) {
				return mList.get(mCurrentIndex);
			} else {
				System.err.println("Invalid Index: " + mCurrentIndex + "\tList Size: " + mList.size());
				((new NoSuchElementException())).printStackTrace();
				System.exit(-1);
			}
			return null;
		}

		// ========================================================================
		/** Returns the index of the element that would be returned by a subsequent call to next(). */
		@Override
		public int nextIndex() {			
			return mCurrentIndex + 1;
		}

		// ========================================================================
		@Override
		public E previous() {
			mList.ifNotSortedThenSort();
			--mCurrentIndex;
			
			if (mCurrentIndex >= 0) {
				return mList.get(mCurrentIndex);
			} else {
				((new NoSuchElementException())).printStackTrace();
				System.exit(-1);
			}
			return null;
		}

		// ========================================================================
		/** Returns the index of the element that would be returned by a subsequent call to previous(). */
		@Override
		public int previousIndex() {
			return mCurrentIndex - 1;			
		}

		// ========================================================================
		/** Removes from the list the last element that was returned by next() or previous() (optional operation). */
		@Override
		public void remove() {
			mList.remove(mCurrentIndex);			
		}

		// ========================================================================
		/** Replaces the last element returned by next() or previous() with the specified element (optional operation). */
		@Override
		public void set(E e) {
			mList.set(mCurrentIndex, e);
		}
	}
	
	
	// ========================================================================
	// ========================================================================
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub

	}

}
