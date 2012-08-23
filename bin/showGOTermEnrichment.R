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

#png("", width=500, height=750)

files=list.files(path=fin, pattern="*csv", full.names=T)
if (!file.exists(outdir))
	dir.create(outdir);
	
for (i in 1:length(files)) {  
	
	a= read.csv(file=files[i], head=TRUE, sep=",")
	enrichment=a$enrichment
	GO_term=a$term
	      
	x <- data.frame(GO_term, enrichment)
	p1 <- ggplot(x, aes(x=GO_term, y=enrichment)) + geom_histogram(width=.5, colour="lightblue", fill="lightblue") + opts(axis.text.y=theme_text(size=6)) + coord_flip()
	
	ggsave(paste(paste(outdir, "/", sep=""), gsub(".csv", ".pdf", basename(files[i])), sep=""), p1)  
	
	dev.off()

}
