package JISTIC;

public class GenomeLocation implements Comparable<GenomeLocation> {
    public int chromosome;
    int start;
    int end;

    protected GenomeLocation(int c, int pos1, int pos2)
    {
	chromosome = c;
	start = Math.min(pos1, pos2);
	end = Math.max(pos1, pos2);
    }

    public int compareTo(GenomeLocation o)
    {
	return chromosome == o.chromosome ?
			 (start == o.start ? end - o.end : start - o.start)
			 : chromosome - o.chromosome;
    }

    public String toString()
    {
	String result = "chr" + chromosome + ':' + start;
	if(end != start)
	    result += "-" + end;
	return result;
    }

    public boolean equals(Object o)
    {
	return (o instanceof GenomeLocation) &&
		same_location((GenomeLocation)o);
    }

    public boolean same_location(GenomeLocation o)
    {
	    return o.chromosome == chromosome
		&& o.start == start
		&& o.end == end;
    }

    public boolean overlaps(GenomeLocation o)
    {
	return chromosome == o.chromosome && start <= o.end && o.start <= end;
    }

    public int inCommon(GenomeLocation o)
    {
	return Math.min(Math.min(o.end - start, end - o.start),
			 Math.min(o.end - o.start, end - start)) + 1;
    }

    public int distance(GenomeLocation o)
    {
	int result;
	if(o == null || chromosome != o.chromosome)
	    result = Integer.MAX_VALUE;
	else if(overlaps(o))
	    result = 0;
	else
	    result = Math.max(start - o.end, o.start - end);
	return result;
    }

    public GenomeLocation end_location()
    {
	return new GenomeLocation(chromosome, end, end);
    }
}
