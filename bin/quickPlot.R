args=(commandArgs(TRUE))
library(ggplot2)

fin=args[1];#"/Users/siddharthreddy/S12/fracplots/scp/ovarian/naf-taf-inputs";
imgname=args[2];#"/Users/siddharthreddy/S12/fracplots/scp/ovarian/quick_plots";

pointSize=.5
alpha=99

files=list.files(path=fin, pattern="*txt", full.names=T)

if (!file.exists(imgname))
	dir.create(imgname);
for (i in 1:length(files)) {
	print(files[i])
	
	png(paste(paste(imgname, "/", sep=""), gsub(".txt", ".png", basename(files[i])), sep=""), width=500, height=400)
	  
  	a=read.table(files[i], sep="\t", header=T)
    
    #SOLiD
  	#tvan = a$TVarCov
  	#ttan = a$TTotCov
  	#nvan = a$NVarCov
  	#ntan = a$NTotCov
  	
	#taf = as.numeric(tvan) / as.numeric(ttan)
  	
  	#print(taf)
  	#print(naf)
  	
  	#Illumina
  	naf=a$VariantRatio_N
  	taf=a$VariantRatio_T
 
  	xlabel="Tumor Allele Fraction"
  	ylabel="Normal Allele Fraction"
  	title = gsub(".txt", "", files[i])
      
  	p1 <- ggplot(a, aes(x=taf, y=naf)) + geom_point()
  	plot(p1)
  	
  	dev.off()
}

