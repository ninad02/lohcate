args=(commandArgs(TRUE))
library(ggplot2)

fin=args[1];#"/Users/siddharthreddy/S12/fracplots/scp/ovarian/curated_snps";
imgname=args[2];#"/Users/siddharthreddy/S12/fracplots/scp/ovarian/";

pointSize=.5
alpha=99

files=list.files(path=fin, pattern="*csv", full.names=T)

if (!file.exists(imgname))
	dir.create(imgname);
for (i in 1:length(files)) {
	print(files[i])
	
	png(paste(paste(imgname, "/", sep=""), gsub(".csv", ".png", basename(files[i])), sep=""), width=500, height=400)
	  
  	a=read.table(files[i], sep=",", header=T)
   	    
  	naf=a$n_vaf
  	taf=a$t_vaf
  	cluster=a$cluster
 
  	xlabel="Tumor Allele Fraction"
  	ylabel="Normal Allele Fraction"
  	title = gsub(".csv", "", files[i])
      
  	p1 <- ggplot(a, aes(x=taf, y=naf, colour=factor(cluster))) + geom_point()
  	plot(p1)
  	
  	dev.off()
}

