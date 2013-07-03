package lohcate;

import java.io.File;

import nutils.EnumMapSafe;

/** 
 * For now, this is just a wrapper class that serves as an entry point to LOHcate.
 * 
 * @author Ninad Dewal
 *
 */
public class LOHcate {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		Script.main(args);
	}
	
	// ========================================================================
	public static enum SubdirsDefault {
		VAFInputsNormalTumor("naf_taf_inputs"),
		SitesClassified("sites_classified"),		
		Plots_VAF_2D        ("plots_VAF_2D"),
		Plots_VAF_GenomeWide("plots_VAF_Genomewide"),
		Plots_CopyNumber    ("plots_copyNumber"),
		Plots_Recurrence    ("plots_recurrence"),
		Regions("regions"),
		BrowserTracks("browser_tracks"),
		GeneEnrichment("gene_enrichment"),
		Simulation("simulation");
		
		String mDirName;
		
		private SubdirsDefault(String dirName) {
			this.mDirName = dirName;
		}
		
		public String getSubDirName() { return mDirName; }
	}
	
	// ========================================================================
	public static class Subdirs {
		private String mDirRoot;
		private EnumMapSafe<SubdirsDefault, String> mSubDirNames;
		
		public Subdirs(String dirRoot) {
			this.mDirRoot = dirRoot;
			this.mSubDirNames = new EnumMapSafe<SubdirsDefault, String>(SubdirsDefault.class);
			initializeDefaultSubDirNames();
		}
		
		public String getSubDirPath(SubdirsDefault sd) { return mSubDirNames.get(sd); }
		
		public String getDirRoot() { return mDirRoot; }
		
		public void setSubDirPath(SubdirsDefault sd, String customSubDirName) {
			mSubDirNames.put(sd, mDirRoot + File.separator + customSubDirName);
		}
		
		private void initializeDefaultSubDirNames() {
			for (SubdirsDefault sd : SubdirsDefault.values()) {
				setSubDirPath(sd, sd.getSubDirName());
			}
		}
		
	}
	
	
}
