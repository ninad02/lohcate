package lohcate;

import java.io.File;
import java.io.PrintStream;

import lohcate.clustering.Clustering;
import lohcate.clustering.ClusteringParams;
import lohcate.geneAnalysis.GeneEnrichment;
import lohcateEnums.SeqPlatform;

import com.martiansoftware.jsap.FlaggedOption;
import com.martiansoftware.jsap.JSAP;
import com.martiansoftware.jsap.JSAPResult;

import nutils.ArgumentParserUtils;
import nutils.CompareUtils;
import nutils.EnumMapSafe;
import nutils.IOUtils;

/** 
 * For now, this is just a wrapper class that serves as an entry point to LOHcate.
 * 
 * @author Ninad Dewal
 *
 */
public class LOHcate {

	public static final String SuffixInputFile = "germline.txt"; //".txt"; //.lohcateInput.txt";
	public static final boolean RunOld = true;
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		entryPoint(args);
	}
	
	// ========================================================================
	public static enum SubdirsDefault {
		VAFInputsNormalTumor(false, "naf_taf_inputs"),
		SitesClassified("sites_classified"),		
		
		Plots_VAF_2D        ("plots" + File.separator + "VAF_2D"),
		Plots_VAF_GenomeWide("plots" + File.separator + "VAF_Genomewide"),
		Plots_VAF_CopyNumber("plots" + File.separator + "VAF_CopyNum_Genomwide"),
		Plots_CopyNumber    ("plots" + File.separator + "copyNumber"),		
		Plots_Recurrence    ("plots" + File.separator + "recurrence"),
		
		Regions_GISTIC(       "regions" + File.separator + "gistic"),
		Regions_JISTIC(       "regions" + File.separator + "jistic"),
		Regions_BrowserTracks("regions" + File.separator + "browser_tracks"),
		
		GeneEnrichment("gene_enrichment"),
		Simulation("simulation");
		
		// --------------------------------------------------------------------
		String mDirName;
		boolean mIsOutputDir;

		// --------------------------------------------------------------------
		private SubdirsDefault(String dirName) {
			this(true, dirName);
		}

		// --------------------------------------------------------------------
		private SubdirsDefault(boolean isOutputDir, String dirName) {
			this.mIsOutputDir = isOutputDir;
			this.mDirName = dirName;			
		}
		
		// --------------------------------------------------------------------
		public String getSubDirName() { return mDirName; }
		
		// --------------------------------------------------------------------
		public boolean isOutputDir() { return mIsOutputDir; }
	}

	// ========================================================================
	// ========================================================================
	// ========================================================================
	public static class Subdirs {
		private String mDirRoot;
		private String mDirRootOutput;		
		private EnumMapSafe<SubdirsDefault, String> mSubDirNames;

		// ========================================================================
		public Subdirs(String dirRoot) {
			this(dirRoot, null);
		}
		
		// ========================================================================
		public Subdirs(String dirRoot, String outputDirRoot) {
			this.mDirRoot = dirRoot;
			this.mDirRootOutput = CompareUtils.isNull(outputDirRoot) ? mDirRoot : outputDirRoot;
			
			this.mSubDirNames = new EnumMapSafe<SubdirsDefault, String>(SubdirsDefault.class);
			initializeDefaultSubDirNames();
		}
		
		// ========================================================================
		public String getSubDirPath(SubdirsDefault sd) { return mSubDirNames.get(sd); }
		
		// ========================================================================
		public String getDirRoot() { return mDirRoot; }
		public String getOutputDirRoot() { return mDirRootOutput; }
		
		// ========================================================================
		public void setSubDirPath(SubdirsDefault sd, String customSubDirName) {
			String dirPathPrefix = sd.isOutputDir() ? mDirRootOutput : mDirRoot; 					
			mSubDirNames.put(sd, dirPathPrefix + File.separator + customSubDirName);
		}
		
		// ========================================================================
		private void initializeDefaultSubDirNames() {
			for (SubdirsDefault sd : SubdirsDefault.values()) {
				setSubDirPath(sd, sd.getSubDirName());
			}
		}
		
		// ========================================================================
		public void createDirectoryPathOnDisk(SubdirsDefault sd) {			
			IOUtils.createDirectoryPath(getSubDirPath(sd), false);					
		}
		
	}

	// ========================================================================
	// ========================================================================
	// ========================================================================
	// ========================================================================
	private static void entryPoint(String[] args) {
		String taskClustering = "clustering";
		String taskRegions = "regions";
		String taskGenes = "genes";
				
		String[] argsMajor = (args.length > 0) ? (new String[] { args[0], args[1] }) : (new String[] { " " });
		JSAP jsapTask = new JSAP();
		
		String lohcateTask = "TaskName";
		FlaggedOption task = new FlaggedOption(lohcateTask).setStringParser(JSAP.STRING_PARSER).setRequired(true)
				.setShortFlag('t').setLongFlag("task").setUsageName(lohcateTask);
		task.setHelp("Indicate which task LOHcate should perform: {" + taskClustering + ", " + taskRegions + ", " + taskGenes + "}");
		ArgumentParserUtils.registerJSAPParameter(jsapTask, task);
		
		String rootFolderPath = "RootFolderPath";
		FlaggedOption rootFolder = new FlaggedOption(rootFolderPath).setStringParser(JSAP.STRING_PARSER).setRequired(true)
				.setShortFlag(JSAP.NO_SHORTFLAG).setLongFlag("root").setUsageName(rootFolderPath);
		rootFolder.setHelp("Indicate the root directory of the data");
		ArgumentParserUtils.registerJSAPParameter(jsapTask, rootFolder);

		String rootFolderPathOutput = "RootFolderPathForOutput";
		FlaggedOption rootFolderOutput = new FlaggedOption(rootFolderPathOutput).setStringParser(JSAP.STRING_PARSER).setRequired(false)
				.setShortFlag(JSAP.NO_SHORTFLAG).setLongFlag("rootOutput").setUsageName(rootFolderPathOutput);
		rootFolderOutput.setHelp("Indicate the root directory of the output.  If not set, this defaults to the required root directory");
		ArgumentParserUtils.registerJSAPParameter(jsapTask, rootFolderOutput);
		
		JSAPResult jsapResult = ArgumentParserUtils.parseAndCheck(argsMajor, jsapTask, LOHcate.class.getName());		
		String taskName = jsapResult.getString(lohcateTask);		
		
		long sys_time_init = System.currentTimeMillis();			
		
		// ==================
		// CLUSTERING
		// ==================
		if (taskName.equals(taskClustering)) {
			
			String allelicBiasFile = "AllelicBiasFile";
			FlaggedOption allelicBias = new FlaggedOption(allelicBiasFile).setStringParser(JSAP.STRING_PARSER).setRequired(false)
					.setShortFlag(JSAP.NO_SHORTFLAG).setLongFlag("allelicBias").setUsageName(allelicBiasFile);
			allelicBias.setHelp("Specify a file that has allelic biases");
			ArgumentParserUtils.registerJSAPParameter(jsapTask, allelicBias);
			
			String simOutFileRootUsage = "SimulationOutputFilenameRoot";
			FlaggedOption simOutFileRoot = new FlaggedOption(simOutFileRootUsage).setStringParser(JSAP.STRING_PARSER).setRequired(false)
					.setShortFlag(JSAP.NO_SHORTFLAG).setLongFlag("simOutputFileRoot").setUsageName(simOutFileRootUsage);
			simOutFileRoot.setHelp("Specify the filename root for simulation/evaluation output");
			ArgumentParserUtils.registerJSAPParameter(jsapTask, simOutFileRoot);		
			
			String simParamsUsage = "SimulationParameters";
			FlaggedOption simParamsDef = new FlaggedOption(simParamsUsage).setStringParser(JSAP.STRING_PARSER).setRequired(false)
					.setShortFlag(JSAP.NO_SHORTFLAG).setLongFlag("simulation").setUsageName(simParamsUsage);
			simParamsDef.setHelp("Specifies whether to run internal testing by simulation");
			ArgumentParserUtils.registerJSAPParameter(jsapTask, simParamsDef);		
			
			ClusteringParams clusteringParams = ClusteringParams.GlobalClusteringParams;
			clusteringParams.registerClusteringParameters(jsapTask);
			
			jsapResult = ArgumentParserUtils.parseAndCheck(args, jsapTask, LOHcate.class.getName());			
			String rootFolderName       = jsapResult.getString(rootFolderPath);
			String rootOutputFolderName = jsapResult.getString(rootFolderPathOutput);
			String allelicBiasFilename = jsapResult.getString(allelicBiasFile);
			String simOutRootFilename  = jsapResult.getString(simOutFileRootUsage);
			String simParamsString     = jsapResult.getString(simParamsUsage);
			clusteringParams.configureParameters(jsapResult);
			
			CompareUtils.ensureTrue(rootFolderName != null, "ERROR: Must specify valid root folder name!");				
			LOHcate.Subdirs subdirs = new LOHcate.Subdirs(rootFolderName, rootOutputFolderName);
			
			LOHcateSimulator.LOHcateSimulatorParams simParams = null;
			PrintStream simOutputStream = System.out;			
			if (simParamsString != null) {				
				ClusteringParams.GlobalClusteringParams.setIsSimulation(true);
				simParams = new LOHcateSimulator.LOHcateSimulatorParams();
				simParams.setParamsWithString(simParamsString);
				simParams.printValues(System.out);		
			}			
			
			//SeqPlatform platform = SeqPlatform.getPlatform(Integer.parseInt(args[2]));
			Clustering.classifySites(subdirs,				
									 allelicBiasFilename,
					                 SeqPlatform.Illumina,
					                 simParams,
					                 simOutRootFilename); //args[2] --> 0::Illumina, 1::SOLiD
			IOUtils.closePrintStream(simOutputStream);		                 
			
			GeneEnrichment.getGeneEnrichment(subdirs.getSubDirPath(SubdirsDefault.SitesClassified), subdirs.getSubDirPath(SubdirsDefault.GeneEnrichment));
			
		} else if (taskName.equals(taskRegions)) {			
			jsapResult = ArgumentParserUtils.parseAndCheck(args, jsapTask, LOHcate.class.getName());			
			String rootFolderName = jsapResult.getString(rootFolderPath);
			LOHcate.Subdirs subdirs = new LOHcate.Subdirs(rootFolderName);
			
			System.out.println("Task: Segmentation of regions...");
			Regions.segmentRegionsAllFiles(subdirs.getSubDirPath(SubdirsDefault.SitesClassified),
								   subdirs.getSubDirPath(SubdirsDefault.Regions_GISTIC),
								   subdirs.getSubDirPath(SubdirsDefault.Regions_BrowserTracks));
								
		} else if (taskName.equals(taskGenes)) {
			jsapResult = ArgumentParserUtils.parseAndCheck(args, jsapTask, LOHcate.class.getName());			
			String rootFolderName = jsapResult.getString(rootFolderPath);
			LOHcate.Subdirs subdirs = new LOHcate.Subdirs(rootFolderName);
			
			GeneEnrichment.getGeneEnrichment(subdirs.getSubDirPath(SubdirsDefault.SitesClassified), subdirs.getSubDirPath(SubdirsDefault.GeneEnrichment));
		}
		
		/*
		String root = args[0]; //project directory
		switch (Integer.parseInt(args[1])) { //args[1] --> 'switchboard' parameter
			case 0:				
				break;
			case 1:
				
				break;
			case 2:
				
				break;
				
			// Everything below this point is Sidd's original code	
				/*
			case 5:
				for (int i = 0; i<Enrichment.cluster_names.length - 1; i++)
					Enrichment.getPathwayEnrichment(root + "/gene_enrichment.csv", root + "/kegg_pathways_roster.tsv", root + "/kegg/pathway_enrichment/" + Enrichment.cluster_names[i] + ".csv", i);
				break;
			case 6:
				Enrichment.annotatePathways(root + "/gene_enrichment.csv", root + "/kegg_pathways_roster.tsv", root + "/KEGG");
				break;
			case 7:
				Enrichment.getGOTermCounts(root + "/gene_enrichment.csv", root + "/gene_association.goa_human", root + "/GO/counts/go_term_counts.csv");
				for (int i = 0; i<Enrichment.cluster_names.length - 1; i++) {
					Enrichment.getGOTermCounts(root + "/gene_enrichment/gene_enrichment_top_" + Enrichment.cluster_names[i] + ".csv", root + "/gene_association.goa_human", root + "/GO/counts/go_term_counts_top_" + Enrichment.cluster_names[i] + ".csv");
					Enrichment.getGOTermEnrichment(root + "/GO/counts/go_term_counts.csv", root + "/GO/counts/go_term_counts_top_" + Enrichment.cluster_names[i] + ".csv", root + "/GO/enrichment/go_term_enrichment_top_" + Enrichment.cluster_names[i] + ".csv");
					
					Enrichment.getGOTermCounts(root + "/gene_enrichment/g1/gene_enrichment_top_" + Enrichment.cluster_names[i] + ".csv", root + "/gene_association.goa_human", root + "/GO/g1/counts/go_term_counts_top_" + Enrichment.cluster_names[i] + ".csv");
					Enrichment.getGOTermEnrichment(root + "/GO/g1/counts/go_term_counts.csv", root + "/GO/g1/counts/go_term_counts_top_" + Enrichment.cluster_names[i] + ".csv", root + "/GO/g1/enrichment/go_term_enrichment_top_" + Enrichment.cluster_names[i] + ".csv");
				}
				break;
				
		}*/
		
		System.out.println("Time elapsed: " + (System.currentTimeMillis()-sys_time_init)/1000 + " seconds");
	}
}
