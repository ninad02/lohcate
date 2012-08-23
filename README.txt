----- LOHcate ----- an automated analysis pipeline for LOH calling and visualization in cancer genomes

--------------------
D Wheeler & SG Reddy -- Human Genome Sequencing Center, Baylor College of Medicine
--------------------

Contact: siddharr@bcm.edu, sidreddy96@gmail.com


Analysis pipeline:

(1) scrape SNP calls from .../atlas-snp/filtered (Illumina) or .../snps (SOLiD)
(2) curate SNP calls (assign clusters, calculate allele fractions, &c.)
(3) generate mouse-ear plots
(4) define regions of event continuity
(5) score regions based on recurrence, exclusivity, event density, etc.
(6) parse CNV data into BED file. if n/a, then generate CNV calls from SNP array data.
(6) push region-score and CNV data into genome browser
(7) generate gene enrichment master table & multi-track histograms
(8) annotate KEGG Pathways with gene enrichment data
(9) generate GO term enrichment histograms


-----	-----


Input: SNP calls, CNV calls (or SNP array), DBSNP, KEGG Pathways, GO

Uses: 
 - Java, R, PennCNV (C++), Unix shell
 - DBSCAN density-based clustering (with slight modification), homebrew region definition & scoring algorithms

Output: 

 - Mouse ear plots of tumor VAF vs. normal AF in SNP calls
 - Genome browser (UCSC) visualizations of LOH/DUP events (and copy-number variation) across patient samples
 - Histograms detailing the distribution of LOH/DUP across genes
 - Histograms detailing the enrichment of Gene Ontology terms in genes with high LOH/DUP
 - Annotated pathways (from KEGG Pathways) that show the relative distribution of LOH/DUP across functional groups of genes

***ADVANTAGE: Can be applied to any cancer with SNP calls.
***DISADVANTAGE: Mouse-ear plot features might vary depending on cancer, sequencing platform (Illumina vs. SOLiD), and sequencing method (whole exome vs. wgs), and thus may not be amenable to our clustering (LOH-HET separation) methods.



**********SCRIPT DOCUMENTATION**********


(*) To scrape naf-taf-inputs (normal allele fraction vs. tumor allele fraction, &c.) from SNP calls, see .../bin/snp_scrape.sh
(*) To execute whole pipeline, see .../bin/master.sh

I've written a 'switchboard' into the main method of Script.java (which acts as a 'Runner' class for the pipeline):

	(*) "java Script [path_to_project_directory] 0 [0/1]" (0::Illumina, 1::SOLiD)

	Parses through naf-taf-inputs and, aside from trivial curation tasks like grabbing dbsnp allele frequencies from the DBSNP column and calculating 	variant allele fractions from coverage values (only for SOLiD), cluster mouse ear plots

	(*) "java Script [path_to_project_directory] 1"
	
	Defines regions of 'continuous' LOH, dup, &c. Assigns LOH/dup/&c. status to somatic variants that fall in the middle of these regions (but were 	not part of the LOH/dup side lobes or HET ball).
	
	(*) "java Script [path_to_project_directory] 2"
	
	Scores regions (defined above) by their recurrence across patient samples, density of LOH/dup/&c. variants, and exclusivity from HET variants.
	
	(*) "java Script [path_to_project_directory] 3"
	
	Parses curated SNP calls, regions, and region scores into BED files that can be uploaded to and viewed in the UCSC genome browser.
	
	(*) "java Script [path_to_project_directory] 4"
	
	Reorganizes the curated SNP call tables to show the distribution of LOH/dup/&c. variants across genes (as well as their recurrence across 		patients). Also generates lists of 'top' 20 genes by different variant counts.	
	
	(*) "java Script [path_to_project_directory] 5"
	
	Shows the enrichment of LOH/dup/&c. variants in all relevant pathways (those containing at least 1 mutated gene) in KEGG Pathways.
	
	(*) "java Script [path_to_project_directory] 6"
	
	Annotate pathways in KEGG Pathways by LOH/dup/&c. variant counts. (-) Light green :: relatively low variant count; (0) Gray :: variant count near 	mean; (+) Bright red :: relatively high variant count.
	
	(*) "java Script [path_to_project_directory] 7"
	
	Shows the enrichment of GO (Gene Ontology) terms in our 'top gene' lists, and thus, the different variant types (LOH/dup/&c.).

I've written a couple of R scripts to generate figures for some of the tables outputted by the scripts above:

(*) "R --slave --args [path_to_project_directory]/naf-taf-inputs [path_to_project_directory]/mouse_ears < genMouseEars.R"

Generates scatter plots ("mouse ear" plots) that show the distribution of a patient's variants' NAFs/TAFs. Colors variants by cluster (LOH/dup side lobes, HET ball, somatic tail, &c.)

(*) "./hist_gen.sh [path_to_project_directory]/gene_enrichment [path_to_project_directory]/gene_enrichment"

Generates histograms for gene enrichment tables.

(*) "R --slave --args [path_to_project_directory]/GO/enrichment [path_to_project_directory]/GO/enrichment < showGOTermEnrichment.R"

Generates histograms for GO term enrichment tables.


**********OUTPUT DOCUMENTATION**********

(*) .../naf-taf-inputs

contains 'raw' snp calls from atlas-snp / snp calling software used for SOLiD data --> tab-delimited (.txt) files -- for each patient -- with the following columns:

---Illumina: 
refName => chromosome #
coordinate => chromosomal position
variantBase-T => variant base in tumor sample (A/T/G/C)
TotCov-N => depth of coverage of position [coordinate] in normal sample
TotCov-T => depth of coverage of position [coordinate] in tumor sample
VarCov-N => # of reads in normal sample that map to position [coordinate] with a variant allele
VarCov-T => # of reads in turbo sample that map to position [coordinate] with a variant allele
VariantRatio_N => \frac{[VarCov-N]}{[TotCov-N]}
VariantRatio_T => \frac{[VarCov-T]}{[TotCov-T]}
dbsnp => dbsnp entry id (contains population allele frequencies for reference/variant bases)
MutationType => nonsynonymous / synonymous / intronic / intergenic / etc.
Hugo_Symbol => gene

---SOLiD:
Chrom <=> Illumina::refName
Coord <=> Illumina::coordinate
Ref => reference base for position [Coord]
TTotCov <=> Illumina::TotCov-T 
TVarCov <=> Illumina::VarCov-T
NTotCov <=> Illumina::TotCov-N
NVarCov <=> Illumina::VarCov-N
DBSNP <=> Illumina::dbsnp (since SOLiD naf-taf-inputs don't provide a 'variant base' column, there's little use for pop. allele frequencies)
ENTERZ_SYM <=> Illumina::Hugo_Symbol	
TYPE1 <=> Illumina::MutationType

(*) .../snps

contains 'curated' snp cals. comma-delimited (csv) files -- for each patient -- with the following columns:

chr => chromsome #
pos => chromosomal position
n_vaf => variant allele fraction in normal sample
t_vaf => variant allele fraction in tumor sample
allele_freq => population frequency of variant allele (from DBSNP)
gene => …well, gene
mutation_type => nonsynonymous / synonymous / intronic / intergenic / etc.
germ_som => "germline" or "somatic", depending on where the snp call is coming from (a .germline.txt naf-taf-input or .somatic.txt naf-taf-input)
cluster => allele-fraction plot cluster assignment (integer pointing to LOH sidelobe, DUP wedge, HET ball, etc.)

(*) .../mouse_ears

contains .png files -- one for each patient -- that plot tumor allele fraction (x-axis) against normal allele fraction (y-axis). these are the basis of our entire analysis; really the meat of the matter.

(*) .../regions

contains sub-directories for each chromosome --> 

(*) .../regions/chr[]
	
contains files for each allele-fraction (hereafter abbreviated AF) plot cluster, excluding the HET ball --> {".../regions/chr[]/dup.csv", ".../regions/chr[]/loh.csv", ".../regions/chr[]/roc-loh.csv"}

btw, "roc-loh" <=> right-of-center LOH

each csv file will have a series of rows --> column 1 :: patient barcode, successive columns :: chromosomal regions ("[start position]-[end position]") of supposedly contiguous LOH/DUP/whatever

(*) .../curated_snps

same as .../snps, except we've gone back and reassigned all somatic mutations that fall within regions of predicted LOH to "LOH" status in the 'cluster' column. this way, we don't lose any variant density in regions that might simply have more somatic mutations than germline. intended to bump up genes like JAK2 to the top of our gene enrichment lists in p. vera data sets.

(*) .../scored_regions

contains sub-directories for each chromosome -->

(*) .../scored_regions/chr[]

contains comma-delimited (csv) files -- for each AF plot cluster -- that contain the following columns:

region => chromosomal region ("[start]-[end]")
size => size of chromosomal region ([end] - [start])
recurrence => fraction of patients in which region is predicted
num_events => tot. # of LOH/DUP/&c. events occurring in region -- across patients
num_hets => tot. # of HET events occurring in region -- across patients
event_density => [num_events] / size
het_density => [num_hets] / size

(*) .../browser_tracks

contains all BED files for UCSC genome browser patient tracks -->

(*) .../browser_tracks/chr[]_[cluster:{dup, loh, &c.}]

contains BED files (in older iterations of the pipeline, they might be .csv files. but pay no heed to the file ext. 'tis misleading) -- for each patient -- that are formatted as such:

"browser position chr1:69511-249150116 <-- specifies overall chromosomal window (big enough to fit all LOH/DUP/&c. events, but could very well be leaving out HET variants that occur elsewhere on the chromosome)
browser hide all
track name="HEPBWP-021-T.csv" description=" " visibility=dense itemRgb="on"
chr1 13448184 13448199 row1 <-- arbitrary name attached to hash mark 0 + 13448184 13448199 205,201,201 <-- rob value that corresponds to mutation status (LOH/HET, nonsynonymous/synonymous)
chr1 89652071 89652094 row1 0 + 89652071 89652094 205,201,201
chr1 94473845 94473846 row2 0 + 94473845 94473846 205,201,201
chr1 204403659 204425028 row3 0 + 204403659 204425028 205,201,201
chr1 240077669 240371106 row4 0 + 240077669 240371106 205,201,201
…"

(*) .../browser_tracks/chr[]_[cluster].txt

a manifest file for patient tracks. you can copy-paste the contents of this file into the 'add custom tracks' field in the UCSC genome browser to batch-upload patient tracks.

(*) .../score_tracks

follows the same structure as .../browser_tracks, except .../browser_tracks/chr[]_[cluster] will contain files like this:

"browser position chr1:1115461-249150116
browser hide all
track name="recurrence" description=" " visibility=dense useScore=1
chr1 248084825 248084825 row1 200.0 <-- a score (can range from 200-900) that specifies the color of the chromosomal region (in shades of grey). higher score <=> darker. + 248084825 248084825
chr1 247006051 247006051 row1 900.0 + 247006051 247006051
..."

there are three files --> .../score_tracks/chr[]_[cluster]/recurrence.csv, .../score_tracks/chr[]_[cluster]/event_density.csv, and .../score_tracks/chr[]_[cluster]/het_density.csv. you will probably only be interested in the recurrence track, since the others are pretty clear from eyeballing the patient track. when you're dealing with >100 patients in the same genome browser window, a recurrence track that 'summarizes' the occurrence of a given region of predicted LOH/DUP/&c. across patients can be useful.

(*) .../gene_enrichment.csv

shows the distribution of variant/HET mutations across genes. looks like this:

"chr,range,gene,nonsynonymous,synonymous,germline,somatic,dup,loh,roc-loh,het,dup_recurrence,loh_recurrence,roc-loh_recurrence,het_recurrence
chr1,69270-69897,OR4F5,2.4849066497880004,2.833213344056216,3.367295829986474,0.0,0.0,0.0,0.0,1.9459101490553132,1,0,0,5
chr1,880091-892569,NOC2L,3.258096538021482,3.9318256327243257,4.68213122712422,0.0,0.0,0.0,0.0,3.258096538021482,0,0,0,10
chr1,897325-900505,KLHL17,0.0,3.4011973816621555,3.4965075614664802,0.0,0.0,1.0986122886681098,0.0,2.4849066497880004,0,1,1,8
..." 

looks prettier in JMP/Excel, but this lists how many nonsynonymous/synonymous/.../het mutations occur in a given gene -- across patients. '[cluster]_recurrence' columns give the tot. # patients in which at least 1 variant was found in a given gene.

(*) .../gene_enrichment

.../gene_enrichment/gene_enrichment_top_[cluster].csv pulls the top 20 genes in .../gene_enrichment.csv by [cluster] variant density (it also takes the natural log of mutation densities, so that when we start looking at histograms, high HET counts won't mask smaller values for LOH/DUP/&c.)

.../gene_enrichment/gene_enrichment/gene_enrichment_top_[cluster].png is a histogram version of its .csv counterpart, containing 4 different tracks that compare (1) LOH/DUP/&c. to HET, (2) nonsynonymous to synonymous, (3) germline to somatic, and (4) recurrence across patients.

.../gene_enrichment/g1 is the same as its parent, except that 'top 20' lists are generated with the stipulation that genes must have a 'region of variance' > 0 (in other words, we don't want to bias our lists to genes that have high density by virtue of only having 1 snp in them, and thus 'size' = 0).

(*) .../kegg/pathway_enrichment

contains comma-delimited (csv) files -- for each cluster -- that show the enrichment of different pathways in KEGG Pathways with LOH/DUP/&c. mutations. each file contains the following columns:

pathway => full name of pathway (string)
kegg_id => 5-digit # corresponding to pathway in KEGG Pathways
tau => Fisher's exact test --> \frac{\frac{|g_p|}{[|P|]}}{\frac{|g_x|}{|X|}}
|g_x| => # genes in .../gene_enrichment.csv that have [cluster] mutation density above the mean (across genes)
|X| => # genes in .../gene_enrichment.csv
|g_p| => # genes in pathway P that have [cluster] mutation density above the mean, as calculated in the .../gene_enrichment.csv pool of genes (not just the ones in pathway P)
|P| => # genes in pathway P that we have data for

(*) .../kegg/[cluster]

contains .png files for each pathway mentioned in .../kegg/pathway_enrichment/[cluster].csv. each gene-wiring diagram is annotated according to the following color scheme/gradient:

light green --> relatively low [cluster] mutation count
gray --> mutation count close to mean
bright red --> relatively high mutation count

(*) .../GO

contains Gene Ontology (GO) term enrichment data -->

.../GO/counts contains 'raw' GO term counts (as in, the # of times a given term occurs in a pool of genes -- non-trivial, since each gene will have multiple GO terms associated with it) -->

.../GO/counts/go_term_counts.csv --> contains term counts for the general pool of genes (every gene listed in .../gene_enrichment.csv)
.../GO/counts/go_term_counts_top_[cluster].csv --> contains term counts for the top 20 pool of genes, as specified in .../gene_enrichment/gene_enrichment_top_[cluster].csv

.../GO/enrichment contains csv (and .png files containing bar graphs) files for each cluster, listing the enrichment of all GO terms associated with genes in the top 20 gene list corresponding with the cluster. Enrichment is calculated using a fold change --> \frac{# occurrences of term in top 20 list}{# of occurrences in general pool of genes}

(*) …/GO/g1

much the same as .../GO, following the same scheme change as what happens between .../gene_enrichment and .../gene_enrichment/g1


**********NOTES ON USE**********

(*) Cluster analysis

 - NAF strip expander parameter --> # of std. deviations away from centre (y-coord, NAF) of HET ball that we should consider when looking for side lobes, wedge, &c.
 - Smaller eps and/or minPts will tighten up cluster capture (smaller HET ball / DUP wedge / whatever)
 - DUP wedge lasso parameter describes how much 'looser' you'd like the eps parameter to be when looking for DUP wedge variants (higher # <=> 'looser'). when you loosen up eps, you'll capture a slightly larger cluster (but with similar density, as long as you keep similar minPts parameter in both HET ball and DUP wedge detection). this results in a 
'lassoing' of DUP wedge points around HET ball.

(*) Region segmentation/scoring/parsing

(*) Enrichment analysis

See attached TeX document / slides.


*********MISC**********

Code comments will most likely not explain enough re cluster analysis, region segmentation, &c., so I've TeX'd a separate document / produced some slides that describe the algorithms I use (whether already existing, like DBSCAN and enrichment analysis, or home-brew, like segmentation/scoring/parsing).
