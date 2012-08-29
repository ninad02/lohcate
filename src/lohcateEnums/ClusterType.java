package lohcateEnums;

public enum ClusterType {
	Dup,
	LOHvar,
	LOHref,  // right-of-center LOH
	
	HET,
	Somatic,
	Noise,
	Null;
	
	public int getNumAberrantClusterTypes() { return values().length - 4; }
	
	public ClusterType getClusterType(int code) {
		return (code >= values().length) ? null : values()[code];
	}
}
