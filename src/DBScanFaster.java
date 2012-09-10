import java.util.ArrayList;

public class DBScanFaster extends DBSCAN2 {

	private static final float AdjustedMinX = 0;
	private static final float AdjustedMinY = 0;
	private static final int EpsilonDivider = 3;   // keep at 2 for optimum performance
	
	protected Block[][] mBlocks;
	protected float mOffsetX;
	protected float mOffsetY;
	protected float mEpsilonDivided;
	
	public DBScanFaster(ArrayList<Floint> points, float epsilon, int minPoints, float minValueX, float minValueY, float maxValueX, float maxValueY) {
		super(points, epsilon, minPoints);
		calculateOffsets(minValueX, minValueY);
		initializeBlocks(minValueX, minValueY, maxValueX, maxValueY);
	}
	
	protected void initializeBlocks(float minValueX, float minValueY, float maxValueX, float maxValueY) {
		mEpsilonDivided = mEpsilon / (float) EpsilonDivider;		
		int numBlocksAlongXAxis = getNumBlocksAlongAxis(minValueX, maxValueX, mEpsilonDivided);
		int numBlocksAlongYAxis = getNumBlocksAlongAxis(minValueY, maxValueY, mEpsilonDivided);
		mBlocks = new Block[numBlocksAlongXAxis][numBlocksAlongYAxis];
		
		for (DBScanPoint point : mPoints) {
			int blockIndexX = getBlockIndexAlongXAxis(point.mFloint.mX);
			int blockIndexY = getBlockIndexAlongYAxis(point.mFloint.mY);
			if (mBlocks[blockIndexX][blockIndexY] == null) {
				mBlocks[blockIndexX][blockIndexY] = new Block(getBlockValueX(blockIndexX), getBlockValueY(blockIndexY), mEpsilonDivided);
			}
			mBlocks[blockIndexX][blockIndexY].mPointsInBlock.add(point);
		}
	}
	
	protected float getBlockValueX(int indexX) { return ((indexX * mEpsilonDivided) + mOffsetX); }
	protected float getBlockValueY(int indexY) { return ((indexY * mEpsilonDivided) + mOffsetY); }
	
	public int getBlockIndexAlongXAxis(float x) { return (int) ((x - mOffsetX) / mEpsilonDivided); }
	public int getBlockIndexAlongYAxis(float y) { return (int) ((y - mOffsetY) / mEpsilonDivided); }
	
	public static int getNumBlocksAlongAxis(float minValueOnAxis, float maxValueOnAxis, float blockLength) {
		return ((int) ((maxValueOnAxis - minValueOnAxis) / blockLength) + 1);
	}
	
	/** Calculates the offsets for x and y. */
	private void calculateOffsets(float minValueX, float minValueY) {
		mOffsetX = minValueX - AdjustedMinX;
		mOffsetY = minValueY - AdjustedMinY;
	}
	
	/** @Override 
	 * 	This procedure will try to optimize determining the neighbors, examining only
	 *  those points within a bounded search space, in which the bounds of the space
	 *  represent the maximum theoretical distances as given 
	 */
	protected ArrayList<DBScanPoint> getNeighbors(DBScanPoint point, ArrayList<DBScanPoint> neighbors, boolean clearList) {
		neighbors = (neighbors == null) ? new ArrayList<DBScanPoint>(20000) : neighbors; 
		if (clearList) { 
			neighbors.clear(); 
		}

		int blockIndexX = getBlockIndexAlongXAxis(point.mFloint.mX);
		int blockIndexY = getBlockIndexAlongYAxis(point.mFloint.mY);
		
		int maxNumAdditionalBlocksThatCouldBeWithinRnage = DBScanFaster.EpsilonDivider;
		int blockIndexMinX = Math.max(0, blockIndexX - maxNumAdditionalBlocksThatCouldBeWithinRnage);
		int blockIndexMinY = Math.max(0, blockIndexY - maxNumAdditionalBlocksThatCouldBeWithinRnage);
		
		int blockIndexMaxX = Math.min(mBlocks.length    - 1, blockIndexX + maxNumAdditionalBlocksThatCouldBeWithinRnage);
		int blockIndexMaxY = Math.min(mBlocks[0].length - 1, blockIndexY + maxNumAdditionalBlocksThatCouldBeWithinRnage);
		
		float pointX = point.mFloint.mX;  // cache values for efficiency
		float pointY = point.mFloint.mY;  
		
		int allCorners = 4;
		for (int xIndex = blockIndexMinX; xIndex <= blockIndexMaxX; xIndex++) {
			for (int yIndex = blockIndexMinY; yIndex <= blockIndexMaxY; yIndex++) {
				Block theBlock = mBlocks[xIndex][yIndex];
				if (theBlock == null) {  // no points exist in the block anyway
					continue;  
				}
				
				int numCornersInRange = theBlock.numCornersWithinDistance(pointX, pointY, mEpsilonSquared);				
				if (numCornersInRange == allCorners) { 
					// All corners of the block fall within the range of the point.  Thus, all block.points
					// within the block will then be within range of the point.  There is no need to test the
					// block.points within the block.
					addAll(neighbors, theBlock.mPointsInBlock);

				} else if (numCornersInRange >= 1 && numCornersInRange < allCorners) {
					// We must test each point within the block
					for (DBScanPoint pointInBlock : theBlock.mPointsInBlock) {
						if (point.getCartesianDistanceSquared(pointInBlock) <= mEpsilonSquared) {
							neighbors.add(pointInBlock);
						}						
					}
				}
			}
		}
		
		return neighbors;
	}
	

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub

	}
	
	// ========================================================================
	// INNER CLASS
	// ========================================================================
	public static class Block {
		public float mLowerLeftCornerX;
		public float mLowerLeftCornerY;
		public float mBlockLength;
		ArrayList<DBScanPoint> mPointsInBlock;
		
		public Block(float lowerLeftCornerX, float lowerLeftCornerY, float blockLength) {
			mLowerLeftCornerX = lowerLeftCornerX;
			mLowerLeftCornerY = lowerLeftCornerY;
			mBlockLength = blockLength;
			mPointsInBlock = new ArrayList<DBScanPoint>();
		}
		
		/** Given a point, this returns how many corners of the block lie within the
		 *  passed-in distance^2.  
		 */
		public int numCornersWithinDistance(float otherX, float otherY, float distanceSquared) {
			int loopMax = 2;  // We hard-set this value because it forces each loop below to run twice
			
			// Cycle through the corners.  We don't use a loop in order for runtime efficiency
			int numCornersWithinRange = 0;
			if (getCartesianDistanceSquared(mLowerLeftCornerX, mLowerLeftCornerY, otherX, otherY) <= distanceSquared) {
				numCornersWithinRange++;
			}
			
			if (getCartesianDistanceSquared(mLowerLeftCornerX + mBlockLength, mLowerLeftCornerY, otherX, otherY) <= distanceSquared) {
				numCornersWithinRange++;
			}

			if (getCartesianDistanceSquared(mLowerLeftCornerX, mLowerLeftCornerY + mBlockLength, otherX, otherY) <= distanceSquared) {
				numCornersWithinRange++;
			}

			if (getCartesianDistanceSquared(mLowerLeftCornerX + mBlockLength, mLowerLeftCornerY + mBlockLength, otherX, otherY) <= distanceSquared) {
				numCornersWithinRange++;
			}

			
			
			/*
			for (int xIndex = 0; xIndex < loopMax; xIndex++) {
				float cornerX = mLowerLeftCornerX + (xIndex * mBlockLength);
				
				for (int yIndex = 0; yIndex < loopMax;  yIndex++) {
					
					float cornerY = mLowerLeftCornerY + (yIndex * mBlockLength);
					float distanceActual = getCartesianDistanceSquared(cornerX, cornerY, otherX, otherY);
					if (distanceActual <= distanceSquared) {
						numCornersWithinRange++;
					}
				}
			}*/
			
			return numCornersWithinRange;
		}
		
		public float getCartesianDistanceSquared(float thisX, float thisY, float otherX, float otherY) {
			return (thisX - otherX) * 
				   (thisX - otherX) +
				   
				   (thisY - otherY) *
				   (thisY - otherY);
		}
	}
}