#!/usr/bin/Rscript

options(stringsAsFactors = FALSE);
args <- commandArgs(trailingOnly=TRUE)

library(IRanges)

ns5.file<-"/mnt/johnson/NS5/PET/RangedData.R"
esc.file<-"/mnt/johnson/ESC/PET/RangedData.R"
astro.file<-"/mnt/rest_chip_astro/Macs/NA_peaks.RangedData.RData"
output.prefix<-"/mnt/data/nsc_esc_astro_overlap"


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

#save(exons, file="/mnt/work/exons.RData")
#save(tss, file="/mnt/work/tss.RData")
#load("/mnt/work/exons.RData")
#load("/mnt/work/tss.RData")

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

maxgap=100
ns5.esc.ol <- findOverlaps(ns5.pet, esc.pet, maxgap=maxgap, minoverlap=1)
astro.esc.ol <- findOverlaps(astro.chipseq, esc.pet, maxgap=maxgap, minoverlap=1)
astro.ns5.ol <- findOverlaps(astro.chipseq, ns5.pet, maxgap=maxgap, minoverlap=1)

#all three overlap (findOverlaps can't do this)
nms <- intersect(names(astro.chipseq), intersect(names(esc.pet), names(ns5.pet)))
astro.esc.ns5.ol <- lapply(nms, function(x){
  n.e <- ns5.esc.ol[[x]]
  a.n <- astro.ns5.ol[[x]]
  if(is.null(n.e) || is.null(a.n)) return ()

  n.e <- as.matrix(n.e)
  a.n <- as.matrix(a.n)

  ol<-merge(a.n, n.e, by.x="subject", by.y="query")
  colnames(ol) <- c( "ns5", "astro","esc")

  return(ol)
  
} )
names(astro.esc.ns5.ol) <- nms
astro.esc.ns5.n <- sum(unlist(lapply( astro.esc.ns5.ol,  function(x){nrow(x)} )))


#fetch gene symbol and description for all the ensembl ids. 
ens.annot <- getBM(mart=ensmart, attributes=c("ensembl_gene_id", "mgi_symbol","description" ))

# we need to merge the symbols and descriptions.
ids <- unique(ens.annot[,"ensembl_gene_id"])
ens.annot <- ens.annot[order(ens.annot[,"ensembl_gene_id"]),]
dups <- unique( ens.annot[which( duplicated(ens.annot[,"ensembl_gene_id"])),"ensembl_gene_id"] )
for(id in dups){
  inds <- which(ens.annot[,"ensembl_gene_id"]==id)
  ens.annot[inds,"mgi_symbol"] = paste(unique(ens.annot[inds,"mgi_symbol"]), collapse=",")
  ens.annot[inds,"description"] = paste( unique(ens.annot[inds,"description"]), collapse=",")
}

dups <- duplicated(ens.annot[,"ensembl_gene_id"])
if(sum(dups)>0) ens.annot <- ens.annot[!dups,]
rownames(ens.annot) <- ids
ens.annot<-ens.annot[,-1]
#save(ens.annot, file="/mnt/work/ensannot.RData")
#load("/mnt/work/ensannot.RData")



#no overlap
esc <- do.call("rbind",lapply(names(esc.pet), function(n){

  es.nearest<-as.data.frame(nearest.tss.es[n])
  this.es <- cbind(es.nearest, ens.annot[es.nearest[,"feature"], ])

  n.e <- ns5.esc.ol[[n]]
  a.e <- astro.esc.ol[[n]]

  rm.inds <-  NULL
  if(!is.null(n.e)) rm.inds <- c(rm.inds, as.matrix(n.e)[,"subject"])
  if(!is.null(a.e)) rm.inds <- c(rm.inds, as.matrix(a.e)[,"subject"])

  rm.inds <- -1*unique(rm.inds)
  if(length(rm.inds) > 0) this.es <- this.es[rm.inds,]
  return(this.es)
  
}))
              
ns5 <- do.call("rbind",lapply(names(ns5.pet), function(n){

  ns.nearest<-as.data.frame(nearest.tss.ns[n])
  this.ns <- cbind(ns.nearest, ens.annot[ns.nearest[,"feature"], ])

  n.e <- ns5.esc.ol[[n]]
  a.n <- astro.ns5.ol[[n]]

  rm.inds <-  NULL
  if(!is.null(n.e)) rm.inds <- c(rm.inds, as.matrix(n.e)[,"query"])
  if(!is.null(a.n)) rm.inds <- c(rm.inds, as.matrix(a.n)[,"subject"])

  rm.inds <- -1*unique(rm.inds)
  if(length(rm.inds)>0) this.ns <- this.ns[rm.inds,]
  return(this.ns)
  
}))
         

astro <- do.call("rbind",lapply(names(astro.chipseq), function(n){

  as.nearest<-as.data.frame(nearest.tss.astro[n])
  this.as <- cbind(as.nearest, ens.annot[as.nearest[,"feature"], ])

  a.e <- astro.esc.ol[[n]]
  a.n <- astro.ns5.ol[[n]]

  rm.inds <-  NULL
  if(!is.null(a.e)) rm.inds <- c(rm.inds, as.matrix(a.e)[,"query"])
  if(!is.null(a.n)) rm.inds <- c(rm.inds, as.matrix(a.n)[,"query"])

  rm.inds <- -1*unique(rm.inds)
  if(length(rm.inds)>0) this.as <- this.as[rm.inds,]
  return(this.as)

  
}))
         


         


usecols <- c("start","end")
ns.es <- do.call("rbind",lapply(intersect(names(ns5.pet), names(esc.pet)),function(n){

  #don't bother if we have no overlaps
  if (is.null(ns5.esc.ol[[n]])) return()

  inds <- as.data.frame(as.matrix(ns5.esc.ol[[n]]))

  # For the overlaps we *do* have, check they don't also overlap
  # with astro sites  
  all.ol <-  astro.esc.ns5.ol[[n]]
  if(!is.null(all.ol)){
    rm.inds <- -1*which(inds[,1] %in% all.ol[,"ns5"])
    if(length(rm.inds) > 0){
      inds <- inds[rm.inds,]
    }
  }
  if(nrow(inds)==0) return()
  
  ns <- as.data.frame(nearest.tss.ns[n])[inds[,1],]
  ns <- cbind(ns,  ens.annot[ns[,"feature"],]) 
  ns <- ns[,-1]
  
  es <- as.data.frame(nearest.tss.es[n])[inds[,2],]
  es <- cbind(es,  ens.annot[es[,"feature"],]) 
  es <- es[,-1]
  
  
  colnames(ns) <- paste("ns5",colnames(ns), sep=".")
  colnames(es) <- paste("esc",colnames(es), sep=".")
  res <- cbind(space=n, ns,es)

}))


astro.es <- do.call("rbind",lapply(intersect(names(astro.chipseq), names(esc.pet)),function(n){

  if (is.null(astro.esc.ol[[n]])) return()

  inds <- as.data.frame(as.matrix(astro.esc.ol[[n]]))
  all.ol <-  astro.esc.ns5.ol[[n]]
  if(!is.null(all.ol)){
    rm.inds <- -1*which(inds[,1] %in% all.ol[,"astro"])
    if(length(rm.inds) > 0){
      inds <- inds[rm.inds,]
    }
  }
  if(nrow(inds)==0) return()
  
  astro <- as.data.frame(nearest.tss.astro[n])[inds[,1],]
  astro <- cbind(astro,  ens.annot[astro[,"feature"],]) 
  astro <- astro[,-1]
  
  es <- as.data.frame(nearest.tss.es[n])[inds[,2],]
  es <- cbind(es,  ens.annot[es[,"feature"],]) 
  es <- es[,-1]
  
  colnames(astro) <- paste("astro",colnames(astro), sep=".")
  colnames(es) <- paste("esc",colnames(es), sep=".")
  res <- cbind(space=n,astro,es)
}))




astro.ns <- do.call("rbind",lapply(intersect(names(astro.chipseq), names(ns5.pet)),function(n){

  if (is.null(astro.ns5.ol[[n]])) return()

  inds <- as.data.frame(as.matrix(astro.ns5.ol[[n]]))
  all.ol <-  astro.esc.ns5.ol[[n]]
  if(!is.null(all.ol)){
    rm.inds <- -1*which(inds[,1] %in% all.ol[,"astro"])
    if(length(rm.inds) > 0){
      inds <- inds[rm.inds,]
    }
  }
  if(nrow(inds)==0) return()
  
  astro <- as.data.frame(nearest.tss.astro[n])[inds[,1],]
  astro <- cbind(astro,  ens.annot[astro[,"feature"],]) 
  astro <- astro[,-1]
  
  ns <- as.data.frame(nearest.tss.ns[n])[inds[,2],]
  ns <- cbind(ns,  ens.annot[ns[,"feature"],]) 
  ns <- ns[,-1]
  
  colnames(astro) <- paste("astro",colnames(astro), sep=".")
  colnames(ns) <- paste("ns5",colnames(ns), sep=".")
  res <- cbind(space=n, astro,ns)
}))


astro.ns.es <- do.call("rbind",lapply(names(astro.esc.ns5.ol), function(n){
  if(is.null(astro.esc.ns5.ol[[n]])){return()}
  
  inds <-  astro.esc.ns5.ol[[n]]
  
  ast <- as.data.frame(nearest.tss.astro[n]) [inds[,"astro"],]
  ast <- cbind(ast, ens.annot[ast[,"feature"],]) 
  ast <- ast[,-1]
                                      
  ns <- as.data.frame(nearest.tss.ns[n])[inds[,"ns5"],]
  ns <- cbind(ns,  ens.annot[ns[,"feature"],]) 
  ns <- ns[,-1]

  es <- as.data.frame(nearest.tss.es[n])[inds[,"esc"],]
  es <- cbind(es,  ens.annot[es[,"feature"],])
  es <- es[,-1]
  
  colnames(ast) <- paste("astro", colnames(ast),sep=".")
  colnames(ns) <- paste("ns5",colnames(ns), sep=".")
  colnames(es) <- paste("esc",colnames(es), sep=".")
  
  res <- cbind(space=n, ast,ns,es)
  
}))

#and write the overlap info out in some useful format? 

write.csv(esc, file=paste(output.prefix, "esc.csv", sep="_"))
write.csv(ns5, file=paste(output.prefix, "ns5.csv", sep="_"))
write.csv(astro, file=paste(output.prefix, "astro.csv", sep="_"))
write.csv(ns.es, file=paste(output.prefix, "ns_es.csv", sep="_"))
write.csv(astro.es, file=paste(output.prefix, "astro_es.csv", sep="_"))
write.csv(astro.ns, file=paste(output.prefix, "astro_ns.csv", sep="_"))
write.csv(astro.ns.es, file=paste(output.prefix, "astro_ns_es.csv", sep="_"))


# total counts
esc.n <- dim(esc.pet)[1]
ns5.n <- dim(ns5.pet)[1]
astro.n <- dim(astro.chipseq)[1]

# venn counts
esc.uniq.n <- dim(esc)[1]
ns5.uniq.n <- dim(ns5)[1]
astro.uniq.n <- dim(astro)[1]

ns.es.n <- dim(ns.es)[1]
astro.es.n <- dim(astro.es)[1]
astro.ns.n <- dim(astro.ns)[1]

astro.es.ns.n <- dim(astro.ns.es)[1]


#create a txt report file with these counts
filename <- paste(output.prefix,"report.txt", sep="_")

cat("ESC Total:", esc.n, "\n", file=filename)
cat("NS5 Total:", ns5.n, "\n", file=filename, append=T)
cat("Astro Total:", astro.n, "\n", file=filename, append=T)

cat("ESC Unique:", esc.uniq.n, "\n", file=filename, append=T)
cat("NS5 Unique:", ns5.uniq.n, "\n", file=filename, append=T)
cat("Astro Unique:", astro.uniq.n, "\n", file=filename, append=T)

cat("NS5/ESC Overlap:", ns.es.n, "\n", file=filename, append=T)
cat("Astro/ESC Overlap:", astro.es.n, "\n", file=filename, append=T)
cat("Astro/NS5 Overlap:", astro.ns.n, "\n", file=filename, append=T)

cat("Astro/ESC/NS5 Overlap:", astro.es.ns.n, "\n", file=filename, append=T)
