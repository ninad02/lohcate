chmod +x snp_scrape.sh
./snp_scrape.sh $1 $2/naf-taf-inputs
java Script $1 0 0
R --slave --args $2/snps $2/mouse_ears < genMouseEars.R
java Script $2 1
java Script $2 2
java Script $2 3
java Script $2 4
chmod +x hist_gen.sh
./hist_gen.sh $2/gene_enrichment $2/gene_enrichment
java Script $2 5
java Script $2 6
java Script $2 7
R --slave --args $2/GO/enrichment $2/GO/enrichment < showGOTermEnrichment.R
