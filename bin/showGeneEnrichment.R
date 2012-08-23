args=(commandArgs(TRUE))

library(ggplot2)
library(reshape2)

multiplot <- function(..., plotlist=NULL, cols) {
    require(grid)

    # Make a list from the ... arguments and plotlist
    plots <- c(list(...), plotlist)

    numPlots = length(plots)

    # Make the panel
    plotCols = cols                          # Number of columns of plots
    plotRows = ceiling(numPlots/plotCols) # Number of rows needed, calculated from # of cols

    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(plotRows, plotCols)))
    vplayout <- function(x, y)
        viewport(layout.pos.row = x, layout.pos.col = y)

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
        curRow = ceiling(i/plotCols)
        curCol = (i-1) %% plotCols + 1
        print(plots[[i]], vp = vplayout(curRow, curCol ))
    }

}

pointSize=1.5
alpha=99

indir = args[1] #.../gene_enrichment/...
outdir = args[2]

files=list.files(path=indir, pattern="*csv", full.names=T)
if (!file.exists(outdir))
	dir.create(outdir);
	
for (i in 1:length(files)) {
	
	png(paste(paste(outdir, "/", sep=""), gsub(".csv", ".png", basename(files[i])), sep=""), width=1400, height=500)
	  
	a = read.csv(file=files[i], head=TRUE, sep=",")
	
	nonsyn=a$nonsynonymous
	syn=a$synonymous
	chr=a$chr
	region=a$region
	gene=a$gene
	germ=a$germline
	som=a$somatic  
	het=a$het
	recurrence=a$recurrence
	event_count=a$event_count
	
	x <- data.frame(gene, chr, nonsyn, syn)
	mx <- melt(x, id.vars=1:2)
	p2 <- ggplot(mx, aes(x=gene, y=value, fill=variable)) + geom_bar(stat="identity") + facet_grid(~chr, scales = "free", space = "free")
	
	x <- data.frame(gene, chr, som, germ)
	mx <- melt(x, id.vars=1:2)
	p3 <- ggplot(mx, aes(x=gene, y=value, fill=variable)) + geom_bar(stat="identity") + facet_grid(~chr, scales = "free", space = "free")
	
	x <- data.frame(gene, chr, recurrence)
	p4 <- ggplot(x, aes(x=gene, y=recurrence)) + geom_bar(stat="identity") + facet_grid(~chr, scales = "free", space = "free")
	
	x <- data.frame(gene, chr, event_count, het)
	mx <- melt(x, id.vars=1:2)
	p1 <- ggplot(mx, aes(x=gene, y=value, fill=variable)) + geom_bar(stat="identity") + facet_grid(~chr, scales = "free", space = "free")
		
	multiplot(p1, p2, p3, p4, cols=1)  
	dev.off()

}
