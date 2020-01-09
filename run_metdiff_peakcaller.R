library(doParallel)
library(foreach)
library(MeTDiff)

#MetDiff function copy without the TREATED BAMS
metdiff.copy <- function(
  IP_BAM,INPUT_BAM,
  GENOME = NA,
  UCSC_TABLE_NAME = "knownGene",
  GENE_ANNO_GTF = NA,
  TXDB = NA,
  OUTPUT_DIR=NA,
  EXPERIMENT_NAME="MeTDiff_output",
  WINDOW_WIDTH=50,
  SLIDING_STEP=50,
  FRAGMENT_LENGTH=100,
  READ_LENGTH=NA,
  MINIMAL_PEAK_LENGTH=FRAGMENT_LENGTH/2,
  PEAK_CUTOFF_PVALUE=NA,
  PEAK_CUTOFF_FDR=0.05,
  FOLD_ENRICHMENT=1,
  DIFF_PEAK_CUTOFF_PVALUE=NA,
  DIFF_PEAK_CUTOFF_FDR=0.05,
  DIFF_PEAK_ABS_FOLD_CHANGE=1,
  MINIMAL_MAPQ=30,
  REMOVE_LOCAL_TAG_ANOMALITIES=TRUE,
  POISSON_MEAN_RATIO=1,
  TESTING_MODE=NA,
  SAVE_INTERMEDIATE=TRUE
){
  
  # Wrap parameters ##################################################
  PARAMETERS=list();
  PARAMETERS$GENE_ANNO_GTF=GENE_ANNO_GTF
  PARAMETERS$IP_BAM=IP_BAM
  PARAMETERS$INPUT_BAM=INPUT_BAM
  PARAMETERS$GENOME = GENOME
  PARAMETERS$UCSC_TABLE_NAME=UCSC_TABLE_NAME
  UCSC_TABLE_NAME = UCSC_TABLE_NAME
  PARAMETERS$FRAGMENT_LENGTH=FRAGMENT_LENGTH
  PARAMETERS$READ_LENGTH=READ_LENGTH
  PARAMETERS$WINDOW_WIDTH=WINDOW_WIDTH
  PARAMETERS$SLIDING_STEP=SLIDING_STEP
  PARAMETERS$MINIMAL_PEAK_LENGTH=MINIMAL_PEAK_LENGTH
  PARAMETERS$PEAK_CUTOFF_PVALUE=PEAK_CUTOFF_PVALUE
  PARAMETERS$PEAK_CUTOFF_FDR=PEAK_CUTOFF_FDR
  PARAMETERS$FOLD_ENRICHMENT=FOLD_ENRICHMENT
  PARAMETERS$MINIMAL_MAPQ=MINIMAL_MAPQ
  PARAMETERS$REMOVE_LOCAL_TAG_ANOMALITIES=REMOVE_LOCAL_TAG_ANOMALITIES
  PARAMETERS$POISSON_MEAN_RATIO=POISSON_MEAN_RATIO
  PARAMETERS$DIFF_PEAK_CUTOFF_PVALUE=DIFF_PEAK_CUTOFF_PVALUE
  PARAMETERS$DIFF_PEAK_CUTOFF_FDR=DIFF_PEAK_CUTOFF_FDR
  PARAMETERS$DIFF_PEAK_ABS_FOLD_CHANGE=DIFF_PEAK_ABS_FOLD_CHANGE
  PARAMETERS$OUTPUT_DIR=OUTPUT_DIR
  PARAMETERS$EXPERIMENT_NAME=EXPERIMENT_NAME
  PARAMETERS$TESTING_MODE=TESTING_MODE
  PARAMETERS$SAVE_INTERMEDIATE=SAVE_INTERMEDIATE
  PARAMETERS$TXDB=TXDB
  
  # check annotation
  if (is.na(PARAMETERS$GENOME) & is.na(PARAMETERS$GENE_ANNO_GTF)) { 
    stop("must specify the genome assembly or provide a gene gtf file for exomePeak to work!", 
         call. = TRUE, domain = NULL)}
  
  # dependent variables
  if (is.na(PARAMETERS$READ_LENGTH)) {PARAMETERS$READ_LENGTH=MeTDiff:::.get.bam.read.length(PARAMETERS$IP_BAM[1])}
  if (is.na(PARAMETERS$MINIMAL_PEAK_LENGTH)) {PARAMETERS$MINIMAL_PEAK_LENGTH=PARAMETERS$FRAGMENT_LENGTH/2}
  if (is.na(PARAMETERS$PEAK_CUTOFF_PVALUE)) {PARAMETERS$PEAK_CUTOFF_TYPE="FDR"} else  {PARAMETERS$PEAK_CUTOFF_TYPE="PVALUE"}
  if (is.na(PARAMETERS$DIFF_PEAK_CUTOFF_PVALUE)) {PARAMETERS$DIFF_PEAK_CUTOFF_TYPE="FDR"} else  {PARAMETERS$DIFF_PEAK_CUTOFF_TYPE="PVALUE"} 
  if (is.na(PARAMETERS$OUTPUT_DIR)) {PARAMETERS$OUTPUT_DIR=getwd()}
  
  # algrithm ##################################################
  
  # read gene annotation
  ANNOTATION = MeTDiff:::.read.gtf(PARAMETERS)
  ANNOTATION_BATCH_ID = MeTDiff:::.divide.anno.into.batches(ANNOTATION)
  
  # index bam files
  # get bam file
  bam=c(PARAMETERS$IP_BAM,PARAMETERS$INPUT_BAM) #removed TREATED BAMS
  no_bam_files=length(bam)
  for (ibam in 1:no_bam_files) {file=bam[ibam]
                                writeLines(file)
                                if (! file.exists(paste(file,'.bai',sep=""))) {
                                  print(paste("Stage: index bam file", file))
                                  indexBam(file)
                                }}
  
  # get reads count report
  SAMPLE_ID = MeTDiff:::.get.sample.id(PARAMETERS)
  BAM_CHRS = MeTDiff:::.get.bam.chrs(PARAMETERS$IP_BAM[1])
  
  # get reads count
  print("Get Reads Count ...")
  print("This step may take a few hours ...")
  if (is.na(PARAMETERS$TESTING_MODE)) {no_batch_to_run=max(ANNOTATION_BATCH_ID)} else {no_batch_to_run=PARAMETERS$TESTING_MODE}
  noGroups=ceiling(no_batch_to_run/170);
  group_bar=round(seq(from = 1, to = no_batch_to_run +1,length.out=noGroups+1))
  # get space
  READS_COUNT=MeTDiff:::.get.reads.count(1,PARAMETERS,ANNOTATION,ANNOTATION_BATCH_ID,BAM_CHRS,no_bam_files,bam)
  READS_COUNT=READS_COUNT[integer(0),]
  READS_COUNT_FORMAT=READS_COUNT
  # groups
  for (igroup in 1:noGroups){
    print( paste(as.character(signif(igroup*100/noGroups, digits = 3)),"%") )
    temp=list()
    batch_batch=group_bar[igroup]:(group_bar[igroup+1]-1)
    reads_count_group=READS_COUNT_FORMAT
    no_batch=length(batch_batch)
    for (ibatch in 1:no_batch){
      temp[[ibatch]]=MeTDiff:::.get.reads.count(batch_batch[ibatch],PARAMETERS,ANNOTATION,ANNOTATION_BATCH_ID,BAM_CHRS,no_bam_files,bam)
    }
    for (ibatch in 1:no_batch){
      reads_count_group=rbind(reads_count_group,temp[[ibatch]])
    }
  READS_COUNT=rbind(READS_COUNT,reads_count_group)}
  READS_COUNT=MeTDiff:::.help.minorm(data.frame(READS_COUNT),SAMPLE_ID)

  # peak calling
  PEAK = MeTDiff:::.peak.call.module(READS_COUNT,SAMPLE_ID,PARAMETERS)
  
  # store the result
  dir.create(paste(PARAMETERS$OUTPUT_DIR,PARAMETERS$EXPERIMENT_NAME,sep='/'),recursive =TRUE,showWarnings = FALSE)
  dir=paste(PARAMETERS$OUTPUT_DIR,PARAMETERS$EXPERIMENT_NAME,sep='/')

  # peak, no diff
  TOTAL_PEAK_RESULT = NA
  TOTAL_PEAK_RESULT=MeTDiff:::.get.table.peak.result(PEAK,ANNOTATION,READS_COUNT,SAMPLE_ID,PARAMETERS,ANNOTATION_BATCH_ID,PEAK$loci2peak_merged)
  
  # save result
  if (PARAMETERS$SAVE_INTERMEDIATE==TRUE) {
    tmp_rs =list(ANNOTATION=ANNOTATION,
                 ANNOTATION_BATCH_ID=ANNOTATION_BATCH_ID,
                 PARAMETERS=PARAMETERS,
                 READS_COUNT=READS_COUNT,
                 SAMPLE_ID=SAMPLE_ID,BAM_CHRS=BAM_CHRS)
     save("tmp_rs", file=paste(dir,'metdiff.Rdata',sep='/'))
  }

  # from .report.diff.peak.based on result (removing differential expression information)
  write.table(TOTAL_PEAK_RESULT,file=paste(dir,"peak.xls",sep="/"), sep="\t",row.names =FALSE,quote = FALSE)
  temp = TOTAL_PEAK_RESULT
  names(temp)[1]=paste("#",names(temp)[1])
  write.table(temp,file=paste(dir,"peak.bed",sep="/"), sep="\t",col.names=FALSE,row.names =FALSE,quote = FALSE)

  return(PEAK)
}

#run the MeTDiff peak caller if it hasn't been run
run.metdiff.peakcaller.notreat <- function(dir,exp,species,fraglen,HG38.GTF,MM10.GTF){
  writeLines(as.character(exp))
  writeLines(as.character(fraglen))
  #if (! file.exists(paste0(dir,"/",exp,"_MeTDiff_output/MeTDiff_output/peak.bed"))){
  if (species == "hg38"){
      gtf <- HG38.GTF } else {
          gtf <- MM10.GTF }
  input.bams <- list.files(paste0(dir,'/alignments'),pattern=paste0("^",exp,"_input_[1-9].star.sorted.bam$"),full.names=TRUE)
  ip.bams <- list.files(paste0(dir,'/alignments'),pattern=paste0("^",exp,"_IP_[1-9].star.sorted.bam$"),full.names=TRUE) 
  metdiff.peaks <- metdiff.copy(GENE_ANNO_GTF = gtf, IP_BAM = ip.bams,INPUT_BAM = input.bams, FRAGMENT_LENGTH = fraglen,
                           OUTPUT_DIR = paste0(dir,"/",exp,"_MeTDiff_output")) 
}

#run MeTDiff if it hasn't been run
run.metdiff <- function(dir,exp,t.exp,species,fraglen,HG38.GTF,MM10.GTF){
  write(paste(dir,exp,t.exp,species,fraglen),file='metdiff_starting.txt',append=TRUE)
  if (! file.exists(paste0(dir,"/",exp,"_",t.exp,"_MeTDiff_output/MeTDiff_output/peak.bed"))){
    if (species == "hg38"){
        gtf <- HG38.GTF} else {
            gtf <- MM10.GTF}
    input.bams <- list.files(paste0(dir,'/alignments'),pattern=paste0("^",exp,"_input_[1-9].star.sorted.bam$"),full.names=TRUE)
    ip.bams <- list.files(paste0(dir,'/alignments'),pattern=paste0("^",exp,"_IP_[1-9].star.sorted.bam$"),full.names=TRUE)
    t.input.bams <- list.files(paste0(dir,'/alignments'),pattern=paste0("^",t.exp,"_input_[1-9].star.sorted.bam$"),full.names=TRUE) 
    t.ip.bams <- list.files(paste0(dir,'/alignments'),pattern=paste0("^",t.exp,"_IP_[1-9].star.sorted.bam$"),full.names=TRUE)
    if (length(input.bams)>1 & length(t.input.bams)>1){
        metdiff.peaks <- MeTDiff::metdiff(GENE_ANNO_GTF = gtf, IP_BAM = ip.bams,INPUT_BAM = input.bams, 
                            TREATED_IP_BAM = t.ip.bams, TREATED_INPUT_BAM = t.input.bams, FRAGMENT_LENGTH = fraglen,
                             OUTPUT_DIR = paste0(dir,"/",exp,"_",t.exp,"_MeTDiff_output"))
    }
  }
  write(paste(dir,exp,t.exp,species,fraglen),file='metdiff_done.txt',append=TRUE)
}


#command line parameters to specify: summary file (eg. exp_summary.txt), number of cores, gtf file location for hg38, gtf file location for mm10
args <- commandArgs(TRUE)
summary <- args[1]
ncores <- as.integer(args[2])
HG38.GTF <- args[3]
MM10.GTF <- args[4]

#set number of cores to use
cl <- makeCluster(ncores)

#run through all experiments for Figures 1,2, and 3 in the exp_summary.txt file provided
exp.summary <- read.csv(summary,sep=' ',row.names=NULL)
fig1.summary <- exp.summary[which(exp.summary$fig == "1a" | exp.summary$fig == "1d"),]
fig23.summary <- exp.summary[which(exp.summary$fig == 2 | exp.summary$fig == 3),]

#run MeTDiff on multiple cores to identify peaks (very slow) - saves results to files
registerDoParallel(cl)
n <- length(fig1.summary$study)
if (n > 0){
    writeLines(paste(summary, HG38.GTF,MM10.GTF))
foreach(dir=fig1.summary$study,exp=fig1.summary$control,species=fig1.summary$species,
        fraglen=fig1.summary$fraglen,.combine='c',.packages = 'MeTDiff') %dopar% {
        run.metdiff.peakcaller.notreat(dir,exp,species,fraglen,HG38.GTF,MM10.GTF)
    }
}

n <- length(fig23.summary$study)
if (n > 0){
sub.summary <- sub.summary[which(sub.summary$control.type == "negative"),]
        writeLines(paste(as.character(sub.summary$study[1]),sub.summary$control[1],sub.summary$species[1]))
        foreach(dir=sub.summary$study,exp=sub.summary$control,species=sub.summary$species,
            fraglen=sub.summary$fraglen,.combine='c',.packages = 'MeTDiff') %dopar% {
            run.metdiff.peakcaller.notreat(dir,exp,species,fraglen)}
fig23.summary <- fig23.summary[which(fig23.summary$control.type != "negative"),]
foreach(dir=fig23.summary$study,exp=fig23.summary$control,t.exp=fig23.summary$stimulus,
            species=fig23.summary$species,fraglen=fig23.summary$fraglen,
            .combine='c',.packages = 'MeTDiff') %dopar% {
        run.metdiff(dir,exp,t.exp,species,fraglen,HG38.GTF,MM10.GTF)
    }   
}

stopCluster(cl)
