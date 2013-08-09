package JISTIC;
import java.util.*;

public class PeakRegion extends Region {

    double peak_qvalue;
    boolean has_focal = true;

    PeakRegion(Marker.aberration type,
		List<Marker> markers, Genome genome, Genome miRNA_genome,
		double peak_qvalue)
    {
	super(type, markers, genome, miRNA_genome);
	this.peak_qvalue = peak_qvalue;
    }

    public String toString() {
	return "PEAK:" + super.toString();
    }
}
