package nutils.collectionsSorted;

import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;


public class ArrayListSorted<E> extends ArrayListSortedAbstract<E> {
	
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	
	// ========================================================================
	/** Constructs an empty list with an initial capacity of ten and with the given comparator. */
	public ArrayListSorted(final Comparator<E> comparator) { this(DefaultNumInitialCapacity, comparator); }
	
	// ========================================================================
	/** Constructs an empty list with the specified initial capacity and a comparator. */
	public ArrayListSorted(final int initialCapacity, final Comparator<E> comparator) {
		super(initialCapacity, comparator);
		checkComparatorForNull(comparator);
	}
	
	// ========================================================================
	/** Constructs a list containing the elements of the specified collection, in the order they are returned by the collection's iterator. */
	public ArrayListSorted(Collection<? extends E> c, final Comparator<E> comparator) {
		super(c, comparator);
		checkComparatorForNull(comparator);
		ifNotSortedThenSort(true);
	}
	
	// ========================================================================
	public ArrayListSorted(ArrayListSorted<E> rhs) {
		super(rhs);		
	}

	// ========================================================================
	protected void checkComparatorForNull(final Comparator<E> comparator) {
		if (comparator == null) {
			throw new IllegalArgumentException("Comparator cannot be null!");
		}
	}
	
	// ========================================================================
	/** Sets the comparator and re-sorts the list if the comparator argument differs from the one already stored internally. 
	 * @return this object, which allows for compounded statements: setComparator(c).get(0) */
	public ArrayListSorted<E> setComparator(final Comparator<E> comparator) {
		checkComparatorForNull(comparator);
		return (ArrayListSorted<E>) super.setComparator(comparator);
	}
	
	// ========================================================================
	/** Returns the index for a specified element. */
	public int getIndex(final E element) {
		ifNotSortedThenSort();		
		return Collections.binarySearch(this, element, mComparator);
	}
		
	// ========================================================================
	/** Returns a shallow copy of this ArrayList instance. */
	@Override
	public Object clone() {
		return new ArrayListSorted<E>(this);
	}
		
	// ========================================================================
	protected int compareHelper(E elementCurrent, E elementPrevOrNext) {
		return (mComparator.compare(elementCurrent, elementPrevOrNext));
	}	
	
	// ========================================================================
	// ========================================================================
	private static void TestSortedList() {		
		Comparator<Integer> intComparatorReverse = new Comparator<Integer>() {
			public int compare(Integer num1, Integer num2) {
				return -num1.compareTo(num2);
			}
		};
		
		Comparator<Integer> theComparator = null;
		ArrayListSorted<Integer> theList = new ArrayListSorted<Integer>(intComparatorReverse);
			
		//for (int i = 4; i >= 0; i--) {
		for (int i = 0; i <= 40; i++) {
			theList.add(i);
			System.out.println(theList.getIndex(i));
			System.out.println(theList);
		}
		
		theList.add(5);		
		System.out.println(theList);
			
		theList.remove(new Integer(3));		
		System.out.println(theList);
		
		theList.add(3);		
		System.out.println(theList);
		
		theList.remove(0);		
		System.out.println(theList);
		
		//theList.setComparator(null);
		//System.out.println(theList);
		
		//theList.setComparator(intComparatorReverse);
		//System.out.println(theList);

		theList.set(4,100);		
		System.out.println(theList);
		System.out.println(theList.clone());
	}
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		TestSortedList();
	}

}
