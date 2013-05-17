package lohcateEnums;

import java.util.ArrayList;

public enum ClusterType {
	GainGermline,	
	GainSomatic,
	LOH,	
	cnLOH,
	HETGermline,
	
	HETSomatic,
	Noise,
	Null;
		
	public static final ClusterType[] AmpLOHHetG = new ClusterType[] { GainSomatic, LOH, HETGermline };
	public static final ClusterType[] AmpLOH = new ClusterType[]     { GainSomatic, LOH };
	public static final ClusterType[] OnlyLOH = new ClusterType[]    { LOH };
	// TODO public static final ClusterType[] AmpLOHcnLOH = new ClusterType[] { GainSomatic, LOH, cnLOH };
	public static final ClusterType[] AmpLOHcnLOH = new ClusterType[] { GainSomatic, cnLOH, LOH };
	
	//public static int getNumAberrantClusterTypes() { return 2; }
	
	public static int getNumClusterTypes() { return values().length; }
	
	
	public static ClusterType getClusterType(int code) {
		return (code >= values().length) ? null : values()[code];
	}
	
	public boolean isLOH() { 
		//return (this == LOH || this == cnLOH);
		return (this == LOH);
	}
	
	public boolean isSomaticEvent() {
		return (this == LOH 
				|| this == cnLOH 
				|| this == GainSomatic
				);	
	}
	
	/** Compare with lowercase to increase compatibility. */
	public static ClusterType getClusterType(String clusterTypeName) {
		for (ClusterType ct : values()) {
			if (ct.name().equalsIgnoreCase(clusterTypeName)) {
				return ct;
			}
		}

		return null;
	}
		
}
