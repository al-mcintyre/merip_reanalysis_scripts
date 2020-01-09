load.fi <- function(fi, anno, annotation.order,bsg){
  peaks <- deq:::import.peaks(fi,anno,annotation.order)
  peaks$peaks$seqs <- BSgenome::getSeq(bsg,peaks$peaks)
  peaks$peaks$length <- GenomicRanges::width(peaks$peaks)
  peaks$peaks$dracn <- count.matches("[AGT][AG]AC",peaks$peaks$seqs) #GGACT also did not reveal differences among peak callers
  return(peaks)
}
