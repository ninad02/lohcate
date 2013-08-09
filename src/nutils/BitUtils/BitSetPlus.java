package nutils.BitUtils;
/*
 * @(#)BitSetPlus.java	1.67 06/04/07
 *
 * Copyright 2006 Sun Microsystems, Inc. All rights reserved.
 * SUN PROPRIETARY/CONFIDENTIAL. Use is subject to license terms.
 */


import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.ObjectStreamField;

import nutils.ArrayUtils;
import nutils.CompareUtils;

/**
 * This class implements a vector of bits that grows as needed. Each
 * component of the bit set has a <code>boolean</code> value. The
 * bits of a <code>BitSetPlus</code> are indexed by nonnegative integers.
 * Individual indexed bits can be examined, set, or cleared. One
 * <code>BitSetPlus</code> may be used to modify the contents of another
 * <code>BitSetPlus</code> through logical AND, logical inclusive OR, and
 * logical exclusive OR operations.
 * <p>
 * By default, all bits in the set initially have the value
 * <code>false</code>.
 * <p>
 * Every bit set has a current size, which is the number of bits
 * of space currently in use by the bit set. Note that the size is
 * related to the implementation of a bit set, so it may change with
 * implementation. The length of a bit set relates to logical length
 * of a bit set and is defined independently of implementation.
 * <p>
 * Unless otherwise noted, passing a null parameter to any of the
 * methods in a <code>BitSetPlus</code> will result in a
 * <code>NullPointerException</code>.
 *
 * <p>A <code>BitSetPlus</code> is not safe for multithreaded use without
 * external synchronization.
 *
 * @author  Arthur van Hoff
 * @author  Michael McCloskey
 * @author  Martin Buchholz
 * @version 1.67, 04/07/06
 * @since   JDK1.0
 */
public class BitSetPlus implements Cloneable, java.io.Serializable {
	/*
	 * BitSetPluss are packed into arrays of "words."  Currently a word is
	 * a long, which consists of 64 bits, requiring 6 address bits.
	 * The choice of word size is determined purely by performance concerns.
	 */
	private final static int ADDRESS_BITS_PER_WORD = 6;
	private final static int BITS_PER_WORD = 1 << ADDRESS_BITS_PER_WORD;
	private final static int BIT_INDEX_MASK = BITS_PER_WORD - 1;

	/* Used to shift left or right for a partial word mask */
	private static final long WORD_MASK = 0xffffffffffffffffL;

	/**
	 * @serialField bits long[]
	 *
	 * The bits in this BitSetPlus.  The ith bit is stored in bits[i/64] at
	 * bit position i % 64 (where bit position 0 refers to the least
	 * significant bit and 63 refers to the most significant bit).
	 */
	private static final ObjectStreamField[] serialPersistentFields = {
		new ObjectStreamField("bits", long[].class),
	};

	/**
	 * The internal field corresponding to the serialField "bits".
	 */
	private long[] mWords;

	/**
	 * The number of words in the logical size of this BitSetPlus.
	 */
	private transient int mWordsInUse = 0;

	/**
	 * Whether the size of "words" is user-specified.  If so, we assume
	 * the user knows what he's doing and try harder to preserve it.
	 */
	private transient boolean mSizeIsSticky = false;

	/* use serialVersionUID from JDK 1.0.2 for interoperability */
	private static final long serialVersionUID = 7997698588986878753L;

	/**
	 * Given a bit index, return word index containing it.
	 */
	private static int wordIndex(int bitIndex) {
		return bitIndex >> ADDRESS_BITS_PER_WORD;
	}

	/**
	 * Every public method must preserve these invariants.
	 */
	private void checkInvariants() {
		assert(mWordsInUse == 0 || mWords[mWordsInUse - 1] != 0);
		assert(mWordsInUse >= 0 && mWordsInUse <= mWords.length);
		assert(mWordsInUse == mWords.length || mWords[mWordsInUse] == 0);
	}

	/**
	 * Set the field wordsInUse with the logical size in words of the bit
	 * set.  WARNING:This method assumes that the number of words actually
	 * in use is less than or equal to the current value of wordsInUse!
	 */
	private void recalculateWordsInUse() {
		// Traverse the BitSetPlus until a used word is found
		int i;
		for (i = mWordsInUse-1; i >= 0; i--)
			if (mWords[i] != 0)
				break;

		mWordsInUse = i+1; // The new logical size
	}

	/**
	 * Creates a new bit set. All bits are initially <code>false</code>.
	 */
	public BitSetPlus() {
		initWords(BITS_PER_WORD);
		mSizeIsSticky = false;
	}

	/**
	 * Creates a bit set whose initial size is large enough to explicitly
	 * represent bits with indices in the range <code>0</code> through
	 * <code>nbits-1</code>. All bits are initially <code>false</code>.
	 *
	 * @param     nbits   the initial size of the bit set.
	 * @exception NegativeArraySizeException if the specified initial size
	 *               is negative.
	 */
	public BitSetPlus(int nbits) {
		// nbits can't be negative; size 0 is OK
		if (nbits < 0)
			throw new NegativeArraySizeException("nbits < 0: " + nbits);

		initWords(nbits);
		mSizeIsSticky = true;
	}

	private void initWords(int nbits) {
		mWords = new long[wordIndex(nbits-1) + 1];
	}

	/**
	 * Ensures that the BitSetPlus can hold enough words.
	 * @param wordsRequired the minimum acceptable number of words.
	 */
	private void ensureCapacity(int wordsRequired) {
		if (mWords.length < wordsRequired) {
			// Allocate larger of doubled size or required size
			int request = Math.max(2 * mWords.length, wordsRequired);
			mWords = copyOf(mWords, request); 
			mSizeIsSticky = false;
		}
	}

	/**
	 * Ensures that the BitSetPlus can accommodate a given wordIndex,
	 * temporarily violating the invariants.  The caller must
	 * restore the invariants before returning to the user,
	 * possibly using recalculateWordsInUse().
	 * @param	wordIndex the index to be accommodated.
	 */
	private void expandTo(int wordIndex) {
		int wordsRequired = wordIndex+1;
		if (mWordsInUse < wordsRequired) {
			ensureCapacity(wordsRequired);
			mWordsInUse = wordsRequired;
		}
	}

	/**
	 * Checks that fromIndex ... toIndex is a valid range of bit indices.
	 */
	private static void checkRange(int fromIndex, int toIndex) {
		if (fromIndex < 0)
			throw new IndexOutOfBoundsException("fromIndex < 0: " + fromIndex);
		if (toIndex < 0)
			throw new IndexOutOfBoundsException("toIndex < 0: " + toIndex);
		if (fromIndex > toIndex)
			throw new IndexOutOfBoundsException("fromIndex: " + fromIndex +
					" > toIndex: " + toIndex);
	}

	/**
	 * Sets the bit at the specified index to the complement of its
	 * current value.
	 *
	 * @param   bitIndex the index of the bit to flip.
	 * @exception IndexOutOfBoundsException if the specified index is negative.
	 * @since   1.4
	 */
	public void flip(int bitIndex) {
		if (bitIndex < 0)
			throw new IndexOutOfBoundsException("bitIndex < 0: " + bitIndex);

		int wordIndex = wordIndex(bitIndex);
		expandTo(wordIndex);

		mWords[wordIndex] ^= (1L << bitIndex);

		recalculateWordsInUse();
		checkInvariants();
	}

	/**
	 * Sets each bit from the specified <tt>fromIndex</tt> (inclusive) to the
	 * specified <tt>toIndex</tt> (exclusive) to the complement of its current
	 * value.
	 *
	 * @param     fromIndex   index of the first bit to flip.
	 * @param     toIndex index after the last bit to flip.
	 * @exception IndexOutOfBoundsException if <tt>fromIndex</tt> is negative,
	 *            or <tt>toIndex</tt> is negative, or <tt>fromIndex</tt> is
	 *            larger than <tt>toIndex</tt>.
	 * @since   1.4
	 */
	public void flip(int fromIndex, int toIndex) {
		checkRange(fromIndex, toIndex);

		if (fromIndex == toIndex)
			return;

		int startWordIndex = wordIndex(fromIndex);
		int endWordIndex   = wordIndex(toIndex - 1);
		expandTo(endWordIndex);

		long firstWordMask = WORD_MASK << fromIndex;
		long lastWordMask  = WORD_MASK >>> -toIndex;
		if (startWordIndex == endWordIndex) {
			// Case 1: One word
			mWords[startWordIndex] ^= (firstWordMask & lastWordMask);
		} else {
			// Case 2: Multiple words
			// Handle first word
			mWords[startWordIndex] ^= firstWordMask;

			// Handle intermediate words, if any
			for (int i = startWordIndex+1; i < endWordIndex; i++)
				mWords[i] ^= WORD_MASK;

			// Handle last word
			mWords[endWordIndex] ^= lastWordMask;
		}

		recalculateWordsInUse();
		checkInvariants();
	}

	/**
	 * Sets the bit at the specified index to <code>true</code>.
	 *
	 * @param     bitIndex   a bit index.
	 * @exception IndexOutOfBoundsException if the specified index is negative.
	 * @since     JDK1.0
	 */
	public void set(int bitIndex) {
		if (bitIndex < 0)
			throw new IndexOutOfBoundsException("bitIndex < 0: " + bitIndex);

		int wordIndex = wordIndex(bitIndex);
		expandTo(wordIndex);

		mWords[wordIndex] |= (1L << bitIndex); // Restores invariants

		checkInvariants();
	}

	/**
	 * Sets the bit at the specified index to the specified value.
	 *
	 * @param     bitIndex   a bit index.
	 * @param     value a boolean value to set.
	 * @exception IndexOutOfBoundsException if the specified index is negative.
	 * @since     1.4
	 */
	public void set(int bitIndex, boolean value) {
		if (value)
			set(bitIndex);
		else
			clear(bitIndex);
	}

	/**
	 * Sets the bits from the specified <tt>fromIndex</tt> (inclusive) to the
	 * specified <tt>toIndex</tt> (exclusive) to <code>true</code>.
	 *
	 * @param     fromIndex   index of the first bit to be set.
	 * @param     toIndex index after the last bit to be set.
	 * @exception IndexOutOfBoundsException if <tt>fromIndex</tt> is negative,
	 *            or <tt>toIndex</tt> is negative, or <tt>fromIndex</tt> is
	 *            larger than <tt>toIndex</tt>.
	 * @since     1.4
	 */
	public void set(int fromIndex, int toIndex) {
		checkRange(fromIndex, toIndex);

		if (fromIndex == toIndex)
			return;

		// Increase capacity if necessary
		int startWordIndex = wordIndex(fromIndex);
		int endWordIndex   = wordIndex(toIndex - 1);
		expandTo(endWordIndex);

		long firstWordMask = WORD_MASK << fromIndex;
		long lastWordMask  = WORD_MASK >>> -toIndex;
		if (startWordIndex == endWordIndex) {
			// Case 1: One word
			mWords[startWordIndex] |= (firstWordMask & lastWordMask);
		} else {
			// Case 2: Multiple words
			// Handle first word
			mWords[startWordIndex] |= firstWordMask;

			// Handle intermediate words, if any
			for (int i = startWordIndex+1; i < endWordIndex; i++)
				mWords[i] = WORD_MASK;

			// Handle last word (restores invariants)
			mWords[endWordIndex] |= lastWordMask;
		}

		checkInvariants();
	}

	/**
	 * Sets the bits from the specified <tt>fromIndex</tt> (inclusive) to the
	 * specified <tt>toIndex</tt> (exclusive) to the specified value.
	 *
	 * @param     fromIndex   index of the first bit to be set.
	 * @param     toIndex index after the last bit to be set
	 * @param     value value to set the selected bits to
	 * @exception IndexOutOfBoundsException if <tt>fromIndex</tt> is negative,
	 *            or <tt>toIndex</tt> is negative, or <tt>fromIndex</tt> is
	 *            larger than <tt>toIndex</tt>.
	 * @since     1.4
	 */
	public void set(int fromIndex, int toIndex, boolean value) {
		if (value)
			set(fromIndex, toIndex);
		else
			clear(fromIndex, toIndex);
	}

	/**
	 * Sets the bit specified by the index to <code>false</code>.
	 *
	 * @param     bitIndex   the index of the bit to be cleared.
	 * @exception IndexOutOfBoundsException if the specified index is negative.
	 * @since     JDK1.0
	 */
	public void clear(int bitIndex) {
		if (bitIndex < 0)
			throw new IndexOutOfBoundsException("bitIndex < 0: " + bitIndex);

		int wordIndex = wordIndex(bitIndex);
		if (wordIndex >= mWordsInUse)
			return;

		mWords[wordIndex] &= ~(1L << bitIndex);

		recalculateWordsInUse();
		checkInvariants();
	}

	/**
	 * Sets the bits from the specified <tt>fromIndex</tt> (inclusive) to the
	 * specified <tt>toIndex</tt> (exclusive) to <code>false</code>.
	 *
	 * @param     fromIndex   index of the first bit to be cleared.
	 * @param     toIndex index after the last bit to be cleared.
	 * @exception IndexOutOfBoundsException if <tt>fromIndex</tt> is negative,
	 *            or <tt>toIndex</tt> is negative, or <tt>fromIndex</tt> is
	 *            larger than <tt>toIndex</tt>.
	 * @since     1.4
	 */
	public void clear(int fromIndex, int toIndex) {
		checkRange(fromIndex, toIndex);

		if (fromIndex == toIndex)
			return;

		int startWordIndex = wordIndex(fromIndex);
		if (startWordIndex >= mWordsInUse)
			return;

		int endWordIndex = wordIndex(toIndex - 1);
		if (endWordIndex >= mWordsInUse) {
			toIndex = length();
			endWordIndex = mWordsInUse - 1;
		}

		long firstWordMask = WORD_MASK << fromIndex;
		long lastWordMask  = WORD_MASK >>> -toIndex;
		if (startWordIndex == endWordIndex) {
			// Case 1: One word
			mWords[startWordIndex] &= ~(firstWordMask & lastWordMask);
		} else {
			// Case 2: Multiple words
			// Handle first word
			mWords[startWordIndex] &= ~firstWordMask;

			// Handle intermediate words, if any
			for (int i = startWordIndex+1; i < endWordIndex; i++)
				mWords[i] = 0;

			// Handle last word
			mWords[endWordIndex] &= ~lastWordMask;
		}

		recalculateWordsInUse();
		checkInvariants();
	}

	/**
	 * Sets all of the bits in this BitSetPlus to <code>false</code>.
	 *
	 * @since   1.4
	 */
	public void clear() {
		while (mWordsInUse > 0)
			mWords[--mWordsInUse] = 0;
	}

	/**
	 * Returns the value of the bit with the specified index. The value
	 * is <code>true</code> if the bit with the index <code>bitIndex</code>
	 * is currently set in this <code>BitSetPlus</code>; otherwise, the result
	 * is <code>false</code>.
	 *
	 * @param     bitIndex   the bit index.
	 * @return    the value of the bit with the specified index.
	 * @exception IndexOutOfBoundsException if the specified index is negative.
	 */
	public boolean get(int bitIndex) {
		if (bitIndex < 0)
			throw new IndexOutOfBoundsException("bitIndex < 0: " + bitIndex);

		checkInvariants();

		int wordIndex = wordIndex(bitIndex);
		return (wordIndex < mWordsInUse)
		&& ((mWords[wordIndex] & (1L << bitIndex)) != 0);
	}

	/**
	 * Returns a new <tt>BitSetPlus</tt> composed of bits from this <tt>BitSetPlus</tt>
	 * from <tt>fromIndex</tt> (inclusive) to <tt>toIndex</tt> (exclusive).
	 *
	 * @param     fromIndex   index of the first bit to include.
	 * @param     toIndex     index after the last bit to include.
	 * @return    a new <tt>BitSetPlus</tt> from a range of this <tt>BitSetPlus</tt>.
	 * @exception IndexOutOfBoundsException if <tt>fromIndex</tt> is negative,
	 *            or <tt>toIndex</tt> is negative, or <tt>fromIndex</tt> is
	 *            larger than <tt>toIndex</tt>.
	 * @since   1.4
	 */
	public BitSetPlus get(int fromIndex, int toIndex) {
		checkRange(fromIndex, toIndex);

		checkInvariants();

		int len = length();

		// If no set bits in range return empty BitSetPlus
		if (len <= fromIndex || fromIndex == toIndex)
			return new BitSetPlus(0);

		// An optimization
		if (toIndex > len)
			toIndex = len;

		BitSetPlus result = new BitSetPlus(toIndex - fromIndex);
		int targetWords = wordIndex(toIndex - fromIndex - 1) + 1;
		int sourceIndex = wordIndex(fromIndex);
		boolean wordAligned = ((fromIndex & BIT_INDEX_MASK) == 0);

		// Process all words but the last word
		for (int i = 0; i < targetWords - 1; i++, sourceIndex++)
			result.mWords[i] = wordAligned ? mWords[sourceIndex] :
				(mWords[sourceIndex] >>> fromIndex) |
				(mWords[sourceIndex+1] << -fromIndex);

			// Process the last word
			long lastWordMask = WORD_MASK >>> -toIndex;
				result.mWords[targetWords - 1] =
					((toIndex-1) & BIT_INDEX_MASK) < (fromIndex & BIT_INDEX_MASK)
					? /* straddles source words */
							((mWords[sourceIndex] >>> fromIndex) |
									(mWords[sourceIndex+1] & lastWordMask) << -fromIndex)
									:
										((mWords[sourceIndex] & lastWordMask) >>> fromIndex);

							// Set wordsInUse correctly
							result.mWordsInUse = targetWords;
							result.recalculateWordsInUse();
							result.checkInvariants();

							return result;
	}

	/**
	 * Returns the index of the first bit that is set to <code>true</code>
	 * that occurs on or after the specified starting index. If no such
	 * bit exists then -1 is returned.
	 *
	 * To iterate over the <code>true</code> bits in a <code>BitSetPlus</code>,
	 * use the following loop:
	 *
	 * <pre>
	 * for (int i = bs.nextSetBit(0); i >= 0; i = bs.nextSetBit(i+1)) {
	 *     // operate on index i here
	 * }</pre>
	 *
	 * @param   fromIndex the index to start checking from (inclusive).
	 * @return  the index of the next set bit.
	 * @throws  IndexOutOfBoundsException if the specified index is negative.
	 * @since   1.4
	 */
	public int nextSetBit(int fromIndex) {
		if (fromIndex < 0)
			throw new IndexOutOfBoundsException("fromIndex < 0: " + fromIndex);

		checkInvariants();

		int u = wordIndex(fromIndex);
		if (u >= mWordsInUse)
			return -1;

		long word = mWords[u] & (WORD_MASK << fromIndex);

		while (true) {
			if (word != 0)
				return (u * BITS_PER_WORD) + Long.numberOfTrailingZeros(word);
			if (++u == mWordsInUse)
				return -1;
			word = mWords[u];
		}
	}

	/**
	 * Returns the index of the first bit that is set to <code>false</code>
	 * that occurs on or after the specified starting index.
	 *
	 * @param   fromIndex the index to start checking from (inclusive).
	 * @return  the index of the next clear bit.
	 * @throws  IndexOutOfBoundsException if the specified index is negative.
	 * @since   1.4
	 */
	public int nextClearBit(int fromIndex) {
		// Neither spec nor implementation handle BitSetPluss of maximal length.
		// See 4816253.
		if (fromIndex < 0)
			throw new IndexOutOfBoundsException("fromIndex < 0: " + fromIndex);

		checkInvariants();

		int u = wordIndex(fromIndex);
		if (u >= mWordsInUse)
			return fromIndex;

		long word = ~mWords[u] & (WORD_MASK << fromIndex);

		while (true) {
			if (word != 0)
				return (u * BITS_PER_WORD) + Long.numberOfTrailingZeros(word);
			if (++u == mWordsInUse)
				return mWordsInUse * BITS_PER_WORD;
			word = ~mWords[u];
		}
	}

	/**
	 * Returns the "logical size" of this <code>BitSetPlus</code>: the index of
	 * the highest set bit in the <code>BitSetPlus</code> plus one. Returns zero
	 * if the <code>BitSetPlus</code> contains no set bits.
	 *
	 * @return  the logical size of this <code>BitSetPlus</code>.
	 * @since   1.2
	 */
	public int length() {
		if (mWordsInUse == 0)
			return 0;

		return BITS_PER_WORD * (mWordsInUse - 1) +
		(BITS_PER_WORD - Long.numberOfLeadingZeros(mWords[mWordsInUse - 1]));
	}

	/**
	 * Returns true if this <code>BitSetPlus</code> contains no bits that are set
	 * to <code>true</code>.
	 *
	 * @return    boolean indicating whether this <code>BitSetPlus</code> is empty.
	 * @since     1.4
	 */
	public boolean isEmpty() {
		return mWordsInUse == 0;
	}

	/**
	 * Returns true if the specified <code>BitSetPlus</code> has any bits set to
	 * <code>true</code> that are also set to <code>true</code> in this
	 * <code>BitSetPlus</code>.
	 *
	 * @param	set <code>BitSetPlus</code> to intersect with
	 * @return  boolean indicating whether this <code>BitSetPlus</code> intersects
	 *          the specified <code>BitSetPlus</code>.
	 * @since   1.4
	 */
	public boolean intersects(BitSetPlus set) {
		for (int i = Math.min(mWordsInUse, set.mWordsInUse) - 1; i >= 0; i--)
			if ((mWords[i] & set.mWords[i]) != 0)
				return true;
		return false;
	}

	/**
	 * Returns the number of bits set to <tt>true</tt> in this
	 * <code>BitSetPlus</code>.
	 *
	 * @return  the number of bits set to <tt>true</tt> in this
	 *          <code>BitSetPlus</code>.
	 * @since   1.4
	 */
	public int cardinality() {
		int sum = 0;
		for (int i = 0; i < mWordsInUse; i++)
			sum += Long.bitCount(mWords[i]);
		return sum;
	}

	/**
	 * Performs a logical <b>AND</b> of this target bit set with the
	 * argument bit set. This bit set is modified so that each bit in it
	 * has the value <code>true</code> if and only if it both initially
	 * had the value <code>true</code> and the corresponding bit in the
	 * bit set argument also had the value <code>true</code>.
	 *
	 * @param   set   a bit set.
	 */
	public void and(BitSetPlus set) {
		if (this == set)
			return;

		while (mWordsInUse > set.mWordsInUse)
			mWords[--mWordsInUse] = 0;

		// Perform logical AND on words in common
		for (int i = 0; i < mWordsInUse; i++)
			mWords[i] &= set.mWords[i];

		recalculateWordsInUse();
		checkInvariants();
	}
	
	public void andFast(BitSetPlus set) {
		CompareUtils.ensureTrue(mWords.length == set.mWords.length, "ERROR: andFast(): Inconsistent size!");
		for (int i = 0; i < mWords.length; i++) {
			mWords[i] &= set.mWords[i];
		}
	}

	public boolean set(BitSetPlus rhs) {
		CompareUtils.ensureTrue(mWords.length == rhs.mWords.length, "ERROR: Inconsistent lengths!");
		if (mWords.length == rhs.mWords.length) {
			for (int i = 0; i < mWords.length; i++) {
				mWords[i] = rhs.mWords[i];
			}
			return true;
		}
		return false;
	}
	
	/** Given two bitsets to and together, this stores the result in this bitset. */
	public boolean andTwoBitsets(BitSetPlus bs1, BitSetPlus bs2) {
		CompareUtils.ensureTrue(bs1.mWords.length == bs2.mWords.length, "ERROR: Inconsistent lengths in argument bitsets!");
		CompareUtils.ensureTrue(bs1.mWords.length == this.mWords.length, "ERROR: Inconsistent length in this bitset!");
	
		int wordsLength = mWords.length;		
		for (int i = 0; i < wordsLength; i++) {
			mWords[i] = (long) (bs1.mWords[i] & bs2.mWords[i]);
		}
		mWordsInUse = wordsLength;
		return true;
	}
	
	/**
	 * Performs a logical <b>OR</b> of this bit set with the bit set
	 * argument. This bit set is modified so that a bit in it has the
	 * value <code>true</code> if and only if it either already had the
	 * value <code>true</code> or the corresponding bit in the bit set
	 * argument has the value <code>true</code>.
	 *
	 * @param   set   a bit set.
	 */
	public void or(BitSetPlus set) {
		if (this == set)
			return;

		int wordsInCommon = Math.min(mWordsInUse, set.mWordsInUse);

		if (mWordsInUse < set.mWordsInUse) {
			ensureCapacity(set.mWordsInUse);
			mWordsInUse = set.mWordsInUse;
		}

		// Perform logical OR on words in common
		for (int i = 0; i < wordsInCommon; i++)
			mWords[i] |= set.mWords[i];

		// Copy any remaining words
		if (wordsInCommon < set.mWordsInUse)
			System.arraycopy(set.mWords, wordsInCommon,
					mWords, wordsInCommon,
					mWordsInUse - wordsInCommon);

		// recalculateWordsInUse() is unnecessary
		checkInvariants();
	}

	/**
	 * Performs a logical <b>XOR</b> of this bit set with the bit set
	 * argument. This bit set is modified so that a bit in it has the
	 * value <code>true</code> if and only if one of the following
	 * statements holds:
	 * <ul>
	 * <li>The bit initially has the value <code>true</code>, and the
	 *     corresponding bit in the argument has the value <code>false</code>.
	 * <li>The bit initially has the value <code>false</code>, and the
	 *     corresponding bit in the argument has the value <code>true</code>.
	 * </ul>
	 *
	 * @param   set   a bit set.
	 */
	public void xor(BitSetPlus set) {
		int wordsInCommon = Math.min(mWordsInUse, set.mWordsInUse);

		if (mWordsInUse < set.mWordsInUse) {
			ensureCapacity(set.mWordsInUse);
			mWordsInUse = set.mWordsInUse;
		}

		// Perform logical XOR on words in common
		for (int i = 0; i < wordsInCommon; i++)
			mWords[i] ^= set.mWords[i];

		// Copy any remaining words
		if (wordsInCommon < set.mWordsInUse)
			System.arraycopy(set.mWords, wordsInCommon,
					mWords, wordsInCommon,
					set.mWordsInUse - wordsInCommon);

		recalculateWordsInUse();
		checkInvariants();
	}

	/**
	 * Clears all of the bits in this <code>BitSetPlus</code> whose corresponding
	 * bit is set in the specified <code>BitSetPlus</code>.
	 *
	 * @param     set the <code>BitSetPlus</code> with which to mask this
	 *            <code>BitSetPlus</code>.
	 * @since     1.2
	 */
	public void andNot(BitSetPlus set) {
		// Perform logical (a & !b) on words in common
		for (int i = Math.min(mWordsInUse, set.mWordsInUse) - 1; i >= 0; i--)
			mWords[i] &= ~set.mWords[i];

		recalculateWordsInUse();
		checkInvariants();
	}

	/**
	 * Returns a hash code value for this bit set. The hash code
	 * depends only on which bits have been set within this
	 * <code>BitSetPlus</code>. The algorithm used to compute it may
	 * be described as follows.<p>
	 * Suppose the bits in the <code>BitSetPlus</code> were to be stored
	 * in an array of <code>long</code> integers called, say,
	 * <code>words</code>, in such a manner that bit <code>k</code> is
	 * set in the <code>BitSetPlus</code> (for nonnegative values of
	 * <code>k</code>) if and only if the expression
	 * <pre>((k&gt;&gt;6) &lt; words.length) && ((words[k&gt;&gt;6] & (1L &lt;&lt; (bit & 0x3F))) != 0)</pre>
	 * is true. Then the following definition of the <code>hashCode</code>
	 * method would be a correct implementation of the actual algorithm:
	 * <pre>
	 * public int hashCode() {
	 *      long h = 1234;
	 *      for (int i = words.length; --i &gt;= 0; ) {
	 *           h ^= words[i] * (i + 1);
	 *      }
	 *      return (int)((h &gt;&gt; 32) ^ h);
	 * }</pre>
	 * Note that the hash code values change if the set of bits is altered.
	 * <p>Overrides the <code>hashCode</code> method of <code>Object</code>.
	 *
	 * @return  a hash code value for this bit set.
	 */
	public int hashCode() {
		long h = 1234;
		for (int i = mWordsInUse; --i >= 0; )
			h ^= mWords[i] * (i + 1);

		return (int)((h >> 32) ^ h);
	}

	/**
	 * Returns the number of bits of space actually in use by this
	 * <code>BitSetPlus</code> to represent bit values.
	 * The maximum element in the set is the size - 1st element.
	 *
	 * @return  the number of bits currently in this bit set.
	 */
	public int size() {
		return mWords.length * BITS_PER_WORD;
	}

	/**
	 * Compares this object against the specified object.
	 * The result is <code>true</code> if and only if the argument is
	 * not <code>null</code> and is a <code>BitSetPlus</code> object that has
	 * exactly the same set of bits set to <code>true</code> as this bit
	 * set. That is, for every nonnegative <code>int</code> index <code>k</code>,
	 * <pre>((BitSetPlus)obj).get(k) == this.get(k)</pre>
	 * must be true. The current sizes of the two bit sets are not compared.
	 * <p>Overrides the <code>equals</code> method of <code>Object</code>.
	 *
	 * @param   obj   the object to compare with.
	 * @return  <code>true</code> if the objects are the same;
	 *          <code>false</code> otherwise.
	 * @see     nutils.BitUtils.nutils.BitSetPlus#size()
	 */
	public boolean equals(Object obj) {
		if (!(obj instanceof BitSetPlus))
			return false;
		if (this == obj)
			return true;

		BitSetPlus set = (BitSetPlus) obj;

		checkInvariants();
		set.checkInvariants();

		if (mWordsInUse != set.mWordsInUse)
			return false;

		// Check words in use by both BitSetPluss
		for (int i = 0; i < mWordsInUse; i++)
			if (mWords[i] != set.mWords[i])
				return false;

		return true;
	}

	/**
	 * Cloning this <code>BitSetPlus</code> produces a new <code>BitSetPlus</code>
	 * that is equal to it.
	 * The clone of the bit set is another bit set that has exactly the
	 * same bits set to <code>true</code> as this bit set.
	 *
	 * <p>Overrides the <code>clone</code> method of <code>Object</code>.
	 *
	 * @return  a clone of this bit set.
	 * @see     nutils.BitUtils.nutils.BitSetPlus#size()
	 */
	public Object clone() {
		if (! mSizeIsSticky)
			trimToSize();

		try {
			BitSetPlus result = (BitSetPlus) super.clone();
			result.mWords = mWords.clone();
			result.checkInvariants();
			return result;
		} catch (CloneNotSupportedException e) {
			throw new InternalError();
		}
	}

	/**
	 * Attempts to reduce internal storage used for the bits in this bit set.
	 * Calling this method may, but is not required to, affect the value
	 * returned by a subsequent call to the {@link #size()} method.
	 */
	private void trimToSize() {
		if (mWordsInUse != mWords.length) {
			mWords = copyOf(mWords, mWordsInUse);
			checkInvariants();
		}
	}

	/**
	 * Save the state of the <tt>BitSetPlus</tt> instance to a stream (i.e.,
	 * serialize it).
	 */
	private void writeObject(ObjectOutputStream s)
	throws IOException {

		checkInvariants();

		if (! mSizeIsSticky)
			trimToSize();

		ObjectOutputStream.PutField fields = s.putFields();
		fields.put("bits", mWords);
		s.writeFields();
	}

	/**
	 * Reconstitute the <tt>BitSetPlus</tt> instance from a stream (i.e.,
	 * deserialize it).
	 */
	private void readObject(ObjectInputStream s)
	throws IOException, ClassNotFoundException {

		ObjectInputStream.GetField fields = s.readFields();
		mWords = (long[]) fields.get("bits", null);

		// Assume maximum length then find real length
		// because recalculateWordsInUse assumes maintenance
		// or reduction in logical size
		mWordsInUse = mWords.length;
		recalculateWordsInUse();
		mSizeIsSticky = (mWords.length > 0 && mWords[mWords.length-1] == 0L); // heuristic
		checkInvariants();
	}

	/**
	 * Returns a string representation of this bit set. For every index
	 * for which this <code>BitSetPlus</code> contains a bit in the set
	 * state, the decimal representation of that index is included in
	 * the result. Such indices are listed in order from lowest to
	 * highest, separated by ",&nbsp;" (a comma and a space) and
	 * surrounded by braces, resulting in the usual mathematical
	 * notation for a set of integers.<p>
	 * Overrides the <code>toString</code> method of <code>Object</code>.
	 * <p>Example:
	 * <pre>
	 * BitSetPlus drPepper = new BitSetPlus();</pre>
	 * Now <code>drPepper.toString()</code> returns "<code>{}</code>".<p>
	 * <pre>
	 * drPepper.set(2);</pre>
	 * Now <code>drPepper.toString()</code> returns "<code>{2}</code>".<p>
	 * <pre>
	 * drPepper.set(4);
	 * drPepper.set(10);</pre>
	 * Now <code>drPepper.toString()</code> returns "<code>{2, 4, 10}</code>".
	 *
	 * @return  a string representation of this bit set.
	 */
	public String toString() {
		checkInvariants();

		int numBits = (mWordsInUse > 128) ?
				cardinality() : mWordsInUse * BITS_PER_WORD;
				StringBuilder b = new StringBuilder(6*numBits + 2);
				b.append('{');

				int i = nextSetBit(0);
				if (i != -1) {
					b.append(i);
					for (i = nextSetBit(i+1); i >= 0; i = nextSetBit(i+1)) {
						int endOfRun = nextClearBit(i);
						do { b.append(", ").append(i); }
						while (++i < endOfRun);
					}
				}

				b.append('}');
				return b.toString();
	}
	
    public static long[] copyOf(long[] original, int newLength) {
        long[] copy = new long[newLength];
        ArrayUtils.arrayCopy(copy, original, Math.min(original.length, newLength));
        return copy;
    }
}
