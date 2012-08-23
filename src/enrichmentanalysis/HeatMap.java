package enrichmentanalysis;

import shared.FileOps;

import java.awt.Color;
import java.awt.Font;
import java.awt.Graphics2D;
import java.awt.image.BufferedImage;
/**
 * LOHcate --- A software tool for LOH calling and visualization in cancer genomes
 * D Wheeler & SG Reddy
 * Human Genome Sequencing Center, Baylor College of Medicine (Houston, TX)
 * 
 * @author Siddharth G. Reddy
 *
 */
public class HeatMap {
	
	private Cell[][] matrix; //[row][col]
	private Cell[] gradient_pos, gradient_neg;
	private BufferedImage image;
	private String[] row_labels;
	
	private Color zero_color = Color.BLACK;
	
	private float mid_val;
	
	public HeatMap() { this(1,1); }
	
	public HeatMap(int num_rows, int num_cols) {
		matrix = new Cell[num_rows][num_cols];
		gradient_pos = new Cell[num_cols];
		gradient_neg = new Cell[num_cols];
	}
	
	public void setZeroColor(Color param) { zero_color = param; }
	public void setGradientPos(Color param, float[] value) { 
		for (int i = 0; i<matrix[0].length; i++)
			gradient_pos[i] = new Cell(param, value[i]);
	}
	public void setGradientNeg(Color param, float[] value) { 
		for (int i = 0; i<matrix[0].length; i++)
			gradient_neg[i] = new Cell(param, value[i]);
	}
	public void setMidVal(float param) { mid_val = param; }
	/*public void setMatrixDims(int num_rows, int num_cols) { matrix = new Cell[num_rows][num_cols]; }
	public void setRowLabels(String[] param) { row_labels = param; }
	public void setCell(Cell param, int row, int col) { matrix[row][col] = param; }
	
	public int getNumCols() { return matrix[0].length; }
	public int getNumRows() { return matrix.length; }*/
	public Cell getGradientPos(int col) { return gradient_pos[col]; }
	public Cell getGradientNeg(int col) { return gradient_neg[col]; }
	public BufferedImage getImage() { return image; }
	//public Cell[][] getMatrix() { return matrix; }
	
	/**
	 * Returns color in the gradient from {gradient_neg.color -> zero_color -> gradient_pos.color} that corresponds to the position of 'value' in the gradient {gradient_neg.value -> mid_val -> gradient_pos.value}
	 */
	public Color getColor(float value, int col) {		
		if (value>mid_val) {
			int[] pos_col = {gradient_pos[col].getColor().getRed(), gradient_pos[col].getColor().getGreen(), gradient_pos[col].getColor().getBlue()};
			int[] neg_col = {zero_color.getRed(), zero_color.getGreen(), zero_color.getBlue()};
			return new Color((int)(pos_col[0] + ( (value - gradient_pos[col].getValue()) / (gradient_neg[col].getValue() - gradient_pos[col].getValue()) ) * (neg_col[0] - pos_col[0]))
					,(int)(pos_col[1] + ( (value - gradient_pos[col].getValue()) / (gradient_neg[col].getValue() - gradient_pos[col].getValue()) ) * (neg_col[1] - pos_col[1]))
					,(int)(pos_col[2] + ( (value - gradient_pos[col].getValue()) / (gradient_neg[col].getValue() - gradient_pos[col].getValue()) ) * (neg_col[2] - pos_col[2]))
					);
		}
		else {
			value = Math.abs(value);
			int[] neg_col = {zero_color.getRed(), zero_color.getGreen(), zero_color.getBlue()};
			int[] pos_col = {gradient_neg[col].getColor().getRed(), gradient_neg[col].getColor().getGreen(), gradient_neg[col].getColor().getBlue()};
			return new Color((int)(pos_col[0] + ( (value - gradient_pos[col].getValue()) / (gradient_neg[col].getValue() - gradient_pos[col].getValue()) ) * (neg_col[0] - pos_col[0]))
					,(int)(pos_col[1] + ( (value - gradient_pos[col].getValue()) / (gradient_neg[col].getValue() - gradient_pos[col].getValue()) ) * (neg_col[1] - pos_col[1]))
					,(int)(pos_col[2] + ( (value - gradient_pos[col].getValue()) / (gradient_neg[col].getValue() - gradient_pos[col].getValue()) ) * (neg_col[2] - pos_col[2]))
					);
		}
	}
	
	/*public void setCellColors() {
		for (int row = 0; row<getNumRows(); row++)
			for (int col = 0; col<getNumCols(); col++)
				matrix[row][col].setColor(getColor(matrix[row][col].getValue(), col));
	}
	
	public void genImage(Cell template) {
		int max_len = 0;
		for (int row = 0; row<getNumRows(); row++)
			if (row_labels[row].length() > max_len)
				max_len = row_labels[row].length();
		
		int row_lbl_buffer = 10 + max_len * 50, bott_buff = 10;
		
		image = new BufferedImage(getNumCols() * template.getWidth() + row_lbl_buffer, getNumRows() * template.getHeight() + bott_buff,BufferedImage.TYPE_INT_ARGB);
		
		Graphics2D obj = image.createGraphics();
		for (int row = 0; row<getNumRows(); row++) {
			obj.setColor(Color.BLACK);
			obj.setFont(new Font("Helvetica", Font.PLAIN, matrix[0][0].getHeight()));
			obj.drawString(row_labels[row], 0, (row+1) * matrix[0][0].getHeight());
			for (int col = 0; col<getNumCols(); col++) {
				obj.setColor(matrix[row][col].getColor());
				obj.fillRect(col * matrix[row][col].getWidth() + row_lbl_buffer + 1, row * matrix[row][col].getHeight() + 1, matrix[row][col].getWidth() - 1, matrix[row][col].getHeight() - 1);
			}
		}
	}
	
	public void exportImage(String outPath) { FileOps.writeImageToFile(image, outPath); }*/
	
}
