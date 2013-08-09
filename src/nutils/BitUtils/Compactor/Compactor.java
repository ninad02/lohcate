package nutils.BitUtils.Compactor;

import nutils.CompareUtils;
import nutils.EnumMapSafe;
import nutils.BitUtils.BitShiftAndMask;


public abstract class Compactor<E extends Enum<E> & CompactorInf<E>> {
	
	protected Class<E> mClass;
	protected EnumMapSafe<E, BitShiftAndMask> mShiftsAndMasks;
	
	// ========================================================================
	public Compactor(Class<E> enumClass, boolean useMostSignificantBit) {
		mClass = enumClass;
		mShiftsAndMasks = new EnumMapSafe<E, BitShiftAndMask>(enumClass);
		initialize(useMostSignificantBit);
	}

	// ========================================================================
	protected abstract int getBitCapacity();
	
	// ========================================================================
	private void initialize(boolean useMostSignificantBit) {
		// Get all the constants in the enum type
		E[] enumValues = mClass.getEnumConstants();
		
		// First, calculate the total number of bits required across the enumerated values
		final int numBitsRequired = totalNumBits(enumValues);
		
		// If the sum exceeds the capacity, then we throw an error
		int effectiveBitCapacity = useMostSignificantBit ? getBitCapacity() : getBitCapacity() - 1; 
		if (numBitsRequired > effectiveBitCapacity) {
			CompareUtils.ensureTrue(false, "ERROR: Cannot fit " + mClass.getName() + " elements, which require " + numBitsRequired + " bits, within " + effectiveBitCapacity + " bits!");			
		}
		
		// Now, we know we can fit the enums within the bits.  The enums are listed 
		// in order from MSB to LSB, so we traverse the array backwards.		
		BitShiftAndMask bsam = null;
		//StringBuilder sb = new StringBuilder(1024);
		for (E enumValue : enumValues) {
			if (bsam == null) {
				// We're at the first enum value
				bsam = new BitShiftAndMask(enumValue.getNumBits(), effectiveBitCapacity - enumValue.getNumBits());
			} else {
				// Create and overwrite bit-shift & mask object
				bsam = BitShiftAndMask.createBitShiftAndMaskInChain(enumValue.getNumBits(), bsam);
			}
			
			// Insert into the map
			mShiftsAndMasks.put(enumValue, bsam);
			
			//bsam.toString(sb, false, "\t");
			//sb.append("\n");
		}		
		
		//System.out.println(sb.toString());
	}

	// ========================================================================
	protected long setValue(E variable, long value, long compactUnit) {
		return mShiftsAndMasks.get(variable).setValueInCompactUnit(value, compactUnit);
	}
	
	// ========================================================================
	protected long getValue(E variable, long compactUnit) {
		return mShiftsAndMasks.get(variable).getValueInCompactUnit(compactUnit);		
	}
	
	// ========================================================================
	private int totalNumBits(E[] enumValues) {
		int numBits = 0;
		for (E enumValue : enumValues) {
			numBits += enumValue.getNumBits();
		}
		return numBits;
	}
}