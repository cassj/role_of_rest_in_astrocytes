#!/usr/bin/Rscript
	 	
source("http://www.bioconductor.org/biocLite.R")

#default bioconductor packages 
biocLite()
	 	
biocLite("beadarray")
	 	
biocLite("IRanges")
	 	
biocLite("biomaRt")
	 	
biocLite("ChIPpeakAnno")

biocLite("BSgenome.Mmusculus.UCSC.mm9")

biocLite("GenomicRanges")`
