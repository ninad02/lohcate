package lohcateEnums;

public enum ClusterType {
	Dup,
	LOHvar,
	LOHref,  // right-of-center LOH
	
	HET,
	Null;
	
	public int getNumAberrantClusterTypes() { return values().length - 2; }
}
