package nutils.collectionsSorted;

import java.lang.reflect.Array;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.Iterator;
import java.util.List;
import java.util.ListIterator;
import java.util.NoSuchElementException;
import java.util.RandomAccess;

public class ArrayListSortedComparable<E extends Comparable<E>> extends ArrayListSortedAbstract<E> {
	
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	protected static int DefaultNumInitialCapacity = 10;
		
	// ========================================================================
	// Constructors
	/** Constructs an empty list with an initial capacity of ten. */
	public ArrayListSortedComparable() { this(DefaultNumInitialCapacity, null); }

	// ========================================================================
	/** Constructs an empty list with an initial capacity of ten and with the given comparator. */
	public ArrayListSortedComparable(Comparator<E> comparator) { this(DefaultNumInitialCapacity, comparator); }
	
	// ========================================================================
	/** Constructs an empty list with the specified initial capacity. */
	public ArrayListSortedComparable(int initialCapacity) { this(initialCapacity, null); }

	// ========================================================================
	/** Constructs an empty list with the specified initial capacity and a comparator. */
	public ArrayListSortedComparable(int initialCapacity, Comparator<E> comparator) {
		super(initialCapacity, comparator);
	}
	
	// ========================================================================
	/** Constructs a list containing the elements of the specified collection, in the order they are returned by the collection's iterator. */
	public ArrayListSortedComparable(Collection<? extends E> c) {
		super(c, null);
	}
	
	// ========================================================================
	public ArrayListSortedComparable(ArrayListSortedComparable<E> rhs) {
		super(rhs);
	}
	
	// ========================================================================
	/** Returns the index for a specified element. */
	public int getIndex(final E element) {
		ifNotSortedThenSort();
		
		if (mComparator == null) {
			return Collections.binarySearch(this, element);
		} else {
			return Collections.binarySearch(this, element, mComparator);
		}
	}
			
	// ========================================================================
	/** Returns a shallow copy of this ArrayList instance. */
	@Override
	public Object clone() {
		return new ArrayListSortedComparable<E>(this);
	}
	
	// ========================================================================
	protected int compareHelper(E elementCurrent, E elementPrevOrNext) {
		if (mComparator == null) {
			return (elementCurrent.compareTo(elementPrevOrNext));					
		} else {
			return (mComparator.compare(elementCurrent, elementPrevOrNext));
		}
	}
	
	
	// ========================================================================
	private static void TestSortedList() {		
		Comparator<Integer> intComparatorReverse = new Comparator<Integer>() {
			public int compare(Integer num1, Integer num2) {
				return -num1.compareTo(num2);
			}
		};
		
		ArrayListSortedComparable<Integer> theList = new ArrayListSortedComparable<Integer>(intComparatorReverse);
			
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
		
		theList.setComparator(null);
		System.out.println(theList);
		
		theList.setComparator(intComparatorReverse);
		System.out.println(theList);

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
