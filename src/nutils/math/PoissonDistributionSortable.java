package nutils.math;

import org.apache.commons.math3.distribution.PoissonDistribution;
import org.apache.commons.math3.exception.NotStrictlyPositiveException;

public class PoissonDistributionSortable extends PoissonDistribution implements Comparable<PoissonDistributionSortable> {

	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;

	public PoissonDistributionSortable(double p)
			throws NotStrictlyPositiveException {
		super(p);
		// TODO Auto-generated constructor stub
	}

	@Override
	public int compareTo(PoissonDistributionSortable rhs) {
		return Double.compare(this.getMean(), rhs.getMean());
	}
	
	
}