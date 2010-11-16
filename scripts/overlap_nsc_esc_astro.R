#!/usr/bin/Rscript

options(stringsAsFactors = FALSE);
args <- commandArgs(trailingOnly=TRUE)

library(IRanges)

#ns5.file<-"/mnt/johnson/NS5/PET/RangedData.R"
#esc.file<-"/mnt/johnson/ESC/PET/RangedData.R"
#astro.file<-"/mnt/rest_chip_astro/Macs/NA_peaks.RangedData.RData"
#output.prefix<-"/mnt/data/nsc_esc_astro_overlap"


ns5.file <- args[1]
esc.file <- args[2]
astro.file <- args[3]

output.prefix <- args[4]

esc.pet <- get(load(esc.file))
ns5.pet <- get(load(ns5.file))
astro.chipseq <- get(load(astro.file))


# get annotation for the positions
library(ChIPpeakAnno)
library(BSgenome.Mmusculus.UCSC.mm9)

library(biomaRt)
ensmart <- useMart("ensembl", dataset="mmusculus_gene_ensembl")

# fetch latest gene position info from ensembl.
#tss <- getAnnotation(ensmart, "TSS")
#exons <- getAnnotation(ensmart, "Exon")

save(exons, file="/mnt/work/exons.RData")
save(tss, file="/mnt/work/tss.RData")
load("/mnt/work/exons.RData")
load("/mnt/work/tss.RData")

# chromosome names are not the same. Note that this will break the NT_ stuff
# but we're not using them here so it doesn't matter
names(tss) <- paste("chr", names(tss), sep="")
names(exons) <- paste("chr", names(exons), sep="")


nearest.tss.es <- annotatePeakInBatch(esc.pet,
                                      AnnotationData=tss,
                                      PeakLocForDistance = "middle",
                                      FeatureLocForDistance = "start",
                                      output = "nearestStart"
                                      )

nearest.tss.ns <- annotatePeakInBatch(ns5.pet,
                                      AnnotationData=tss,
                                      PeakLocForDistance = "middle",
                                      FeatureLocForDistance = "start",
                                      output = "nearestStart"
                                      )

nearest.tss.astro <- annotatePeakInBatch(astro.chipseq,
                                      AnnotationData=tss,
                                      PeakLocForDistance = "middle",
                                      FeatureLocForDistance = "start",
                                      output = "nearestStart"
                                      )


#annotatePeakInBatch reorders stuff. Sort everything on space, start
esc.pet<-esc.pet[order(space(esc.pet),start(esc.pet)),]
ns5.pet<-ns5.pet[order(space(ns5.pet),start(ns5.pet)),]
astro.chipseq<-astro.chipseq[order(space(astro.chipseq),start(astro.chipseq)),]

nearest.tss.es<-nearest.tss.es[order(space(nearest.tss.es),start(nearest.tss.es)),]
nearest.tss.ns<-nearest.tss.ns[order(space(nearest.tss.ns),start(nearest.tss.ns)),]
nearest.tss.astro<-nearest.tss.astro[order(space(nearest.tss.astro),start(nearest.tss.astro)),]


##annotatePeakInBatch is really slow for overlaps, this is much faster
# findOverlaps doesn't reorder.
overlapping.exon.ns <- findOverlaps(ns5.pet, exons, type="any", select="first")
overlapping.exon.es <- findOverlaps(esc.pet, exons, type="any", select="first")
overlapping.exon.astro <- findOverlaps(astro.chipseq, exons, type="any", select="first")




esc.n <- dim(esc.pet)[1]
ns5.n <- dim(ns5.pet)[1]
astro.n <- dim(astro.chipseq)[1]

ns5.esc.ol <- findOverlaps(ns5.pet, esc.pet, maxgap=500, minoverlap=1)
astro.esc.ol <- findOverlaps(astro.chipseq, esc.pet, maxgap=500, minoverlap=1)
astro.ns5.ol <- findOverlaps(astro.chipseq, ns5.pet, maxgap=500, minoverlap=1)

#get overlap counts for venn
ns5.esc.n <- nrow(as.matrix(ns5.esc.ol))
astro.esc.n <- nrow(as.matrix(astro.esc.ol))
astro.ns5.n <- nrow(as.matrix(astro.ns5.ol))



#fetch gene symbol and description for all the ensembl ids. 
ens.annot <- getBM(mart=ensmart, attributes=c("ensembl_gene_id", "mgi_curated_gene_symbol","description" ))
rownames(ens.annot) <- ens.annot[,"ensembl_gene_id"]
ens.annot<-ens.annot[,-1]

#no overlap
esc <- do.call("rbind",lapply(names(ns5.esc.ol), function(n){

  n.e <- ns5.esc.ol[[n]]
  a.e <- astro.esc.ol[[n]]
  if(is.null(n.e) || is.null(a.e)) return ()

  n.e <- as.matrix(n.e)
  a.e <- as.matrix(a.e)

  this.esc <- as.data.frame(esc.pet[n])[,c("space",usecols)]
  es.nearest<-as.data.frame(nearest.tss.es[n])
  this.es <- cbind(this.esc, es.nearest, ens.annot[es.nearest[,"feature"], ])
    
  inds <- -1*(union(n.e[,"subject"], a.e[,"subject"]))
  es <- this.es[inds, ]
  
}))
              
ns5 <- do.call("rbind",lapply(names(ns5.esc.ol), function(n){

  n.e <- ns5.esc.ol[[n]]
  a.n <- astro.ns5.ol[[n]]
  if(is.null(n.e) || is.null(a.n)) return ()

  n.e <- as.matrix(n.e)
  a.n <- as.matrix(a.n)

  this.ns5 <- as.data.frame(ns5.pet[n])[,c("space",usecols)]
  ns.nearest<-as.data.frame(nearest.tss.ns[n])
  this.ns <- cbind(this.ns5, ns.nearest, ens.annot[ns.nearest[,"feature"],])
  
    
  inds <- -1*(union(n.e[,"query"], a.n[,"subject"]))
  ns <- this.ns[inds, ]
  
}))
         

astro <- do.call("rbind",lapply(names(ns5.esc.ol), function(n){

  a.e <- astro.esc.ol[[n]]
  a.n <- astro.ns5.ol[[n]]
  if(is.null(a.e) || is.null(a.n)) return ()

  a.e <- as.matrix(a.e)
  a.n <- as.matrix(a.n)

  this.astro <- as.data.frame(astro.chipseq[n])[,c("space",usecols)]
  astro.nearest<-as.data.frame(nearest.tss.astro[n])
  this.astro <- cbind(this.astro, astro.nearest, ens.annot[astro.nearest[,"feature"],])
    
  inds <- -1*(union(a.e[,"query"], a.n[,"query"]))
  astro <- this.astro[inds, ]
  
}))
         



#all three overlap (findOverlaps can't do this)
astro.both.ol <- lapply(names(ns5.esc.ol), function(x){
  n.e <- ns5.esc.ol[[x]]
  a.n <- astro.ns5.ol[[x]]
  if(is.null(n.e) || is.null(a.n)) return ()

  n.e <- as.matrix(n.e)
  a.n <- as.matrix(a.n)

  ol<-merge(a.n, n.e, by.x="subject", by.y="query")
  colnames(ol) <- c( "ns5", "astro","esc")

  return(ol)
  
} )
names(astro.both.ol) <- names(ns5.esc.ol)
astro.both.n <- sum(unlist(lapply( astro.both.ol,  function(x){nrow(x)} )))


#create a txt report file with these counts
filename <- paste(output.prefix,"report.txt", sep="_")
cat("ESC total:", dim(esc.pet), "\n" , file=filename)
cat("NS5:total:", dim(ns5.pet), "\n", file=filename, append=T)
cat("Astro:total:", dim(astro.chipseq), "\n", file=filename, append=T)
cat("NS/ES overlap:", ns5.esc.n, "\n", file=filename, append=T)
cat("Astro/ES overlap:",astro.esc.n,  "\n", file=filename, append=T)
cat("Astro/NS overlap:", astro.ns5.n, "\n", file=filename, append=T)
cat("NS/ES/Astro overlap:", astro.both.n, "\n", file=filename, append=T)


#get each venn subset:

nms<-names(ns5.esc.ol)
usecols <- c("start","end")


ns.es <- do.call("rbind",lapply(as.list(nms),function(n){
  inds <- as.matrix(ns5.esc.ol[[n]])
  ns <- as.data.frame(ns5.pet[n])[inds[,1], usecols]
  es <- as.data.frame(esc.pet[n])[inds[,2], usecols]
  
  ns.nearest<-as.data.frame(nearest.tss.ns[n])[inds[,1],]
  ns.nearest <- cbind(ns.nearest,  ens.annot[ns.nearest[,"feature"],]) 

  es.nearest<-as.data.frame(nearest.tss.es[n])[inds[,2],]
  es.nearest <- cbind(es.nearest,  ens.annot[es.nearest[,"feature"],]) 

  ns <- cbind(ns,ns.nearest)
  es <- cbind(es,es.nearest)
  
  colnames(ns) <- paste("ns5",colnames(ns), sep=".")
  colnames(es) <- paste("esc",colnames(es), sep=".")
  res <- cbind(space=n, ns,es)

}))


astro.es <- do.call("rbind",lapply(as.list(nms),function(n){
  if (is.null(astro.esc.ol[[n]])) return()
  inds <-  as.matrix(astro.esc.ol[[n]])
  astro <- as.data.frame(astro.chipseq[n])[inds[,1], usecols]
  es <- as.data.frame(esc.pet[n])[inds[,2], usecols]

  astro.nearest<-as.data.frame(nearest.tss.astro[n])[inds[,1],]
  astro.nearest <- cbind(astro.nearest,  ens.annot[astro.nearest[,"feature"],]) 

  es.nearest<-as.data.frame(nearest.tss.es[n])[inds[,2],]
  es.nearest <- cbind(es.nearest,  ens.annot[es.nearest[,"feature"],]) 
  
  astro <- cbind(astro,astro.nearest)
  es <- cbind(es,es.nearest)
  
  
  colnames(astro) <- paste("astro",colnames(astro), sep=".")
  colnames(es) <- paste("esc",colnames(es), sep=".")
  res <- cbind(space=n,astro,es)
}))


astro.ns <- do.call("rbind",lapply(as.list(nms),function(n){
  if (is.null(astro.ns5.ol[[n]])) return()
  inds <-  as.matrix(astro.ns5.ol[[n]])
  astro <- as.data.frame(astro.chipseq[n])[inds[,1], usecols]
  ns <- as.data.frame(ns5.pet[n])[inds[,2], usecols]

  astro.nearest<-as.data.frame(nearest.tss.astro[n])[inds[,1],]
  astro.nearest <- cbind(astro.nearest,  ens.annot[astro.nearest[,"feature"],]) 

  ns.nearest<-as.data.frame(nearest.tss.ns[n])[inds[,2],]
  ns.nearest <- cbind(ns.nearest,  ens.annot[ns.nearest[,"feature"],]) 

  astro <- cbind(astro,astro.nearest)
  ns <- cbind(ns,ns.nearest)

  colnames(astro) <- paste("astro",colnames(astro), sep=".")
  colnames(ns) <- paste("ns5",colnames(ns), sep=".")
  res <- cbind(space=n, astro,ns)
}))


astro.ns.es <-do.call("rbind", lapply(as.list(nms), function(n){
  this<-astro.both.ol[[n]]
  if(is.null(nrow(this))) return()

  astro <- as.data.frame(astro.chipseq[n])[this[,"astro"], usecols]
  ns <- as.data.frame(ns5.pet[n])[this[,"ns5"], usecols]
  es <- as.data.frame(esc.pet[n])[this[,"esc"], usecols]

  astro.nearest<-as.data.frame(nearest.tss.astro[n])[this[,"astro"],]
  astro.nearest <- cbind(astro.nearest,  ens.annot[astro.nearest[,"feature"],]) 

  ns.nearest<-as.data.frame(nearest.tss.ns[n])[this[,"ns5"],]
  ns.nearest <- cbind(ns.nearest,  ens.annot[ns.nearest[,"feature"],]) 

  es.nearest<-as.data.frame(nearest.tss.es[n])[this[,"esc"],]
  es.nearest <- cbind(es.nearest,  ens.annot[es.nearest[,"feature"],]) 

  astro <- cbind(astro,astro.nearest)
  ns <- cbind(ns, ns.nearest)
  es <- cbind(es, es.nearest)
  
  colnames(astro) <- paste("astro",colnames(astro), sep=".")
  colnames(ns) <- paste("ns5",colnames(ns), sep=".")
  colnames(es) <- paste("esc",colnames(es), sep=".")

  res <- cbind(space=n, astro,ns, es)
  
  
} ))



#and write the overlap info out in some useful format? 

write.csv(esc, file=paste(output.prefix, "esc.csv", sep="_"))
write.csv(ns5, file=paste(output.prefix, "ns5.csv", sep="_"))
write.csv(astro, file=paste(output.prefix, "astro.csv", sep="_"))
write.csv(ns.es, file=paste(output.prefix, "ns_es.csv", sep="_"))
write.csv(astro.es, file=paste(output.prefix, "astro_es.csv", sep="_"))
write.csv(astro.ns, file=paste(output.prefix, "astro_ns.csv", sep="_"))
write.csv(astro.ns.es, file=paste(output.prefix, "astro_ns_es.csv", sep="_"))
