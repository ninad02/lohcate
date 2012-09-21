package lohcateEnums;

import java.util.ArrayList;

public enum ClusterType {
	Amp,
	LOH,	
	HETGermline,
	
	HETSomatic,
	Noise,
	Null;
	
	public static final ArrayList<String> namesLowercase = getClusterTypesLowercase();
	public static final ClusterType[] AmpLOHHetG = new ClusterType[] { Amp, LOH, HETGermline };
	public static final ClusterType[] AmpLOH = new ClusterType[]     { Amp, LOH };
	public static final ClusterType[] OnlyLOH = new ClusterType[]         { LOH };
	
	public static int getNumAberrantClusterTypes() { return 2; }
	
	public static int getNumClusterTypes() { return values().length; }
	
	
	public static ClusterType getClusterType(int code) {
		return (code >= values().length) ? null : values()[code];
	}
	
	/** Compare with lowercase to increase compatibility. */
	public static ClusterType getClusterType(String clusterTypeName) {
		clusterTypeName = clusterTypeName.toLowerCase();
		for (int i = 0; i < namesLowercase.size(); i++) {
			if (clusterTypeName.equals(namesLowercase.get(i))) {
				return getClusterType(i);
			}
		}
		return null;
	}
	
	private static ArrayList<String> getClusterTypesLowercase() {
		ArrayList<String> clusterTypeNamesLowerCase = new ArrayList<String>(values().length);
		for (ClusterType ct : values()) {
			clusterTypeNamesLowerCase.add(ct.name().toLowerCase());
		}
		return clusterTypeNamesLowerCase;
	}
}
