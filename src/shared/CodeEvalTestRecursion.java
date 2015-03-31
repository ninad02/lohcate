package shared;

import java.util.ArrayList;
import java.util.Collections;

import nutils.ArrayUtils;
import nutils.NumberUtils;

public class CodeEvalTestRecursion {

	public static int blockCount = 0;
	public static ArrayList<Point> points = new ArrayList<>(500000);	
	public static Point dummy = new Point(0, 0);
	
	public static void PerformTest() {
		Point p = new Point(0, 0);
		points.add(p);
		move(p);
		System.out.println("# Num Blocks Visted: " + blockCount);
	}

	public static void move(Point p) {
		boolean pass = conditionPasses(p, 19); // || (sumOfDigits(p) != 8);
		if (p.visited || !pass) {
			return;
		} else { 
			p.visited = true;
			blockCount++;
			
			// Move left
			testAndAdd(p.x - 1, p.y);
			testAndAdd(p.x + 1, p.y);
			testAndAdd(p.x, p.y - 1);
			testAndAdd(p.x, p.y + 1);
		}
	}
	
	public static void testAndAdd(int newX, int newY) {
		dummy.x = newX;
		dummy.y = newY;
		int resultIndex = Collections.binarySearch(points, dummy);
		if (resultIndex < 0) {
			Point pNew = new Point(dummy.x, dummy.y);
			points.add(ArrayUtils.getInsertPoint(resultIndex), pNew);
			move(pNew);
		}		
	}
	
	
	
	public static boolean conditionPasses(Point p, int thresh) {		
		return (sumOfDigits(p) <= thresh);
	}
	
	public static int sumOfDigits(Point p) {
		return sumOfDigits(Math.abs(p.x)) + sumOfDigits(Math.abs(p.y));
	}
	
	public static int sumOfDigits(int n) {
		int sum = 0;
		while (n > 0) {
			sum += (n % 10);
			n /= 10;
		}
		return sum;
	}
	
	public static class Point implements Comparable<Point> {
		int x = 0, y = 0;
		boolean visited;
		
		public Point(int x, int y) {
			this.x = x;
			this.y = y;
			this.visited = false;
		}
		
		public int compareTo(Point rhs) {			
			int result = Integer.compare(x, rhs.x);
			if (result == 0) {
				result = Integer.compare(y, rhs.y);
			}
			return result;
		}
	}
	
	public static void main(String[] args) {
		//Test3();
		TestGenericTypes();
		//System.out.println(NumberUtils.MathLog2(1));
		//dbsnptest();
		
		
		/*
		//System.out.println(sumOfDigits(55));
		long time1 = System.currentTimeMillis();
		PerformTest();
		long time2 = System.currentTimeMillis();
		System.out.println(time2 - time1);
		for (Point p : points) {
		//	System.out.println(p.x + "\t" + p.y);
		}
		*/
		
	}
	
	
	private static void dbsnptest() {
		String s = "rs141081392,byfrequency,A|0.000;G|0.999;T|0.000,multiAllele";
		s = "rs148291270,byfrequency;alternate_allele,G|0.998;C|0.002,.";
		s = "rs74900103,unknown,.,.";
		String[] cols = s.split(",");
		//System.out.println(cols);
		for (String col : cols) {
			System.out.println(col);
		}
		
		String[] allelesAndFrequencies = cols[2].split(";|\\|");
		for (String alleleAndFreq : allelesAndFrequencies) {
			System.out.println(alleleAndFreq);
		}
		
		
	}

	
	private static class MyInt {
		int mValue;		
		
		public MyInt(int value) {
			mValue = value;
		}
	}
	
	private static void Test3() {
		
		ArrayList myList = new ArrayList();
		for (int i = 0; i < 10; i++) {
			myList.add(new MyInt(i));
		}
		
		for (int i = 0; i < 10; i++) {
			MyInt myInt = (MyInt) myList.get(i);
			//System.out.println(myInt.mValue);
		}
		
		
		ArrayList<MyInt> myList2 = new ArrayList<MyInt>();
		for (int i = 0; i < 10; i++) {
			myList2.add(new MyInt(i));
		}
		
		for (int i = 0; i < 10; i++) {
			MyInt myInt = myList2.get(i);
			System.out.println(myInt.mValue);
		}
		
		
		
	}
	
	public static interface TestInterface<E> {
		public E foo(E e);
	}
	
	public static class TestA implements TestInterface<TestA> {
		public TestA foo(TestA a) { return this; }
	}
	
	// =====
	public static interface CloneInf<T> {
		public T clone(); 
	}
	
	public static abstract class A<T extends A<T> & CloneInf<T>> {
		public abstract T foo();
		public abstract T clone();
	}
	
	public static class B extends A<B> implements CloneInf<B> {

		@Override
		public B foo() {
			System.out.println("Running foo " + this.getClass().getName());
			return this;
		}		
		
		public B clone() { return null; }
	}
	
	public static class C extends B  {
		public C foo() { 
			System.out.println("Running foo " + this.getClass().getName());
			return this;
		}
		
		@Override
		public C clone() { return null; }
	}
	

	public static void TestGenericTypes_Help(A<?> objA) {
		ArrayList<C> listC = new ArrayList<C>();
		listC.add((C) objA);
	}
	
	public static void TestGenericTypes() {
		B objB = new C();
		A<?> objA = objB;
		
		objB.foo();
		objA.foo();
		
		
	}
	
//	public static class TestB extends TestA implements TestInterface<TestB> {		
//		public TestB foo(TestB b) { return this; }
//	}
}
