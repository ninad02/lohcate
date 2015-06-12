package nutils.objectPool;

import static org.junit.Assert.*;

import java.util.ArrayList;

import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

public class ObjectPoolLimitedTest {

	private static final int NumElementsToTest = 100;
	
	// ========================================================================
	ObjectPoolLimitedString_TestClass mTestClass;
	ArrayList<String>          mElementBuffer;
	
	@BeforeClass
	public static void setUpBeforeClass() throws Exception {
		
	}

	@AfterClass
	public static void tearDownAfterClass() throws Exception {
	}

	@Before
	public void setUp() throws Exception {
		mTestClass = new ObjectPoolLimitedString_TestClass(NumElementsToTest, "StringPool");
		mElementBuffer = new ArrayList<String>();
	}

	@After
	public void tearDown() throws Exception {
	}

	@Test
	public void testObjectPoolLimited() {
		System.out.println("---Testing Constructor---");
		mTestClass = new ObjectPoolLimitedString_TestClass(NumElementsToTest, "StringPool");		
		org.junit.Assert.assertTrue(mTestClass.size() == NumElementsToTest);
	}

	@Test
	public void testBorrowElement() {
		System.out.println("---Testing Borrowing---");
		for (int i = 0; i < NumElementsToTest; i++) {
			mElementBuffer.add(mTestClass.borrowElement());
		}
		System.out.println(mElementBuffer);
		
		org.junit.Assert.assertTrue(mTestClass.size() == 0);
		org.junit.Assert.assertNull(mTestClass.borrowElement());		
	}

	@Test
	public void testReturnElement() {
		
		System.out.println("---Testing Returning---");
		
		// Borrowing elements
		for (int i = 0; i < NumElementsToTest; i++) {
			mElementBuffer.add(mTestClass.borrowElement());
		}
		
		// Check the size
		org.junit.Assert.assertTrue(mTestClass.size() == 0);
		
		// Cache an element
		String tempElement = mElementBuffer.get(0);

		// Return to pool
		for (String s : mElementBuffer) {
			org.junit.Assert.assertTrue(mTestClass.returnElement(s));
		}
		
		// Check the size
		org.junit.Assert.assertTrue(mTestClass.size() == NumElementsToTest);
		
		// Try to return extra element
		org.junit.Assert.assertFalse(mTestClass.returnElement(tempElement));
	}

	@Test
	public void testSize() {
		System.out.println("---Testing Size---");
		assertTrue(mTestClass.size() == NumElementsToTest);
	}

}
