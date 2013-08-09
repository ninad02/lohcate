package lohcateEnums;

import java.util.ArrayList;

public enum EventType {
	GainGermline,	
	GainSomatic,
	LOH,	
	cnLOH,
	HETGermline,
	
	HETSomatic,
	DELHom,
	//Homozygous,
	Noise,
	Ignored;
		
	public static final EventType[] AmpLOHHetG = new EventType[] { GainSomatic, LOH, HETGermline };
	public static final EventType[] AmpLOH = new EventType[]     { GainSomatic, LOH };
	public static final EventType[] OnlyLOH = new EventType[]    { LOH };
	public static final EventType[] AmpLOHcnLOH = new EventType[] { GainSomatic, LOH, cnLOH };
	//public static final ClusterType[] AmpLOHcnLOH = new ClusterType[] { GainSomatic, cnLOH, LOH };
	
	//public static int getNumAberrantClusterTypes() { return 2; }
	
	public static int getNumClusterTypes() { return values().length; }
	
	
	public static EventType getClusterType(int code) {
		return (code >= values().length) ? null : values()[code];
	}
	
	public boolean isLOH() { 
		//return (this == LOH || this == cnLOH);
		return (this == LOH);
	}
	
	public boolean isSomaticCopyNumberEvent() {
		return (   this == LOH 
				|| this == cnLOH 
				|| this == GainSomatic
				);	
	}
	
	public boolean isCopyNumberChangeGermlineOrSomatic() {
		return (isSomaticCopyNumberEvent() || (this == GainGermline)); 
	}
	
	/** Compare with lowercase to increase compatibility. */
	public static EventType getClusterType(String clusterTypeName) {
		for (EventType ct : values()) {
			if (ct.name().equalsIgnoreCase(clusterTypeName)) {
				return ct;
			}
		}

		return null;
	}
		
}
