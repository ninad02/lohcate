package lohcateEnums;

public enum ClusterType {
	Dup,
	LOH,
	rocLOH,  // right-of-center LOH
	
	HET,
	Null;
	
	public int getNumAberrantClusterTypes() { return values().length - 2; }
}
