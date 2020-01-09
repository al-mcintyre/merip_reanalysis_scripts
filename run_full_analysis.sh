#!/bin/bash

NUM_THREADS=16 #set number of threads to run analyses
#must provide locations for hg38 and mm10 gtf reference files
HG38_GTF=''
MM10_GTF=''

#select steps and analyses to run
RUN_STAR=false #change to true to start analysis from star alignments (requires you have downloaded and named files in accordance with exp_summary.txt)
RUN_MACS2=false #run peak calling using MACS2 (requires alignments in bam files)
RUN_METDIFF=false #run peak calling using MeTDiff (requires alignments in bam files)
RUN_KALLISTO=false #run kallisto to estimate gene expression for Figure 1 (requires fastqs named appropriately)
PLOT_FIG1=true #plot Figure 1, requires results from RUN_STAR, RUN_KALLISTO, and ≥ 1 peak caller
PLOT_FIG2=true #plot Figure 2, requires results from RUN_STAR, RUN_KALLISTO, and ≥ 1 peak caller
PLOT_FIG3=true #plot Figure 3, requires results from RUN_STAR and ≥ 1 peak caller
PLOT_FIG4=true #plot Figure 4A, requires results from RUN_STAR and ≥ 1 peak caller
PLOT_FIG5=true #plot Figure 5, required files all present in folder "fig5"
#to produce most figures from provided peak and deq summary files: set RUN_MACS2, RUN_METDIFF, RUN_KALLISTO, and RUN_DEQ to false, all plots to true
#see specific sections below for more details

## Step 1: align data, produce bams

#run STAR alignments for paired end and single end data
#this requires the summary file exp_summary.txt, with directory names in the first column and species in the third column
#also requires Trimmomatic jar and adapter fastas - need to set paths in run_star.1.sh script 
#each directory must contain the raw fastq files for a different data set within a folder named "fastqs"
#(fastq names  should be formatted as ${condition}_${IP or input}_${rep #}_R${1 or 2}.fastq.gz, eg. stress_IP_1_R1.fastq.gz)
#this step will generate bam files within a folder called "alignments"

if $RUN_STAR; then 
#set locations of star indexes below for human and mouse
SI_HG38=''
SI_MM10=''

while read p; do
    DIR=$(echo $p | cut -d' ' -f 1)
    if [ ! $DIR == "study"  ]; then
        REPS=$(echo $p | cut -d' ' -f 11)
        GENOME=$(echo $p | cut -d' ' -f 3)
        if [ "$GENOME" == "mm10" ]; then
            STAR_INDEX=$SI_MM10
        elif [ "$GENOME" == "hg38" ]; then 
            STAR_INDEX=$SI_HG38
        fi
        echo $DIR $REPS $STAR_INDEX $NUM_THREADS
        bash run_star.1.sh $DIR $REPS $STAR_INDEX $NUM_THREADS
    fi
done < exp_summary.txt
fi

## Step 2: call peaks 
#run MACS2 and/or MeTDiff to call peaks based on IP and input bam files
#this requires the successful completion of STAR, i.e.  bam files for each data set within a folder named "alignments"
#this step will generate bed files in folders macs2_results and/or {$condition(s)}_MeTDiff_output/MeTDiff_output
#also produces the peak vs read count plots for SFig2 and the DRAC motif enrichment plot for SFig1g

if [ ! -d fig1 ]; then mkdir fig1; fi
if [ ! -d fig2 ]; then mkdir fig2; fi

PEAKCALLERS=()
if $RUN_MACS2 || $RUN_METDIFF; then
if $RUN_METDIFF; then
    Rscript run_metdiff_peakcaller.R exp_summary.txt $NUM_THREADS $HG38_GTF $MM10_GTF
    echo 'done MeTDiff peak calling'
    PEAKCALLERS+=('metdiff')
fi
if $RUN_MACS2; then PEAKCALLERS+=('macs2'); fi
while read p; do
    DIR=$(echo $p | cut -d' ' -f1)
    if [ "$DIR" != "study" ]; then
    LABEL=$(echo $p | cut -d' ' -f2)
    CELLLINE=$(echo $p | cut -d' ' -f4)
    COND=$(echo $p | cut -d' ' -f6)
    YEAR=$(echo $p | cut -d' ' -f13)
    FRAGLEN=$(echo $p | cut -d' ' -f9)
    REPS=$(echo $p | cut -d' ' -f 11)
    FIG=$(echo $p | cut -d' ' -f 14)

    if $RUN_MACS2; then
        #bash run_macs2.2.sh $DIR $FRAGLEN $REPS $NUM_THREADS
        if [ "$FIG" == "2" ]; then 
            continue #bash plot_sfig2.sh $DIR $CELLLINE $COND $YEAR $REPS $LABEL macs2
        fi
    fi
    if $RUN_METDIFF; then
        for FI in $DIR/*_MeTDiff_output/MeTDiff_output/peak.bed; do
            finame="${FI%.*}"
            bedtools sort -i $FI > $finame.sorted.bed
        done
        if [ "$FIG" == "2" ]; then
            bash plot_sfig2.sh $DIR $CELLLINE $COND $YEAR $REPS $LABEL metdiff
        fi
    fi
    fi
done < exp_summary.txt

#plot peaks vs reads comparison and peak DRAC motif enrichment
for PEAKCALLER in "${PEAKCALLERS[@]}"; do
    #Rscript plot_sfig2.R $PEAKCALLER 'IP'
    #Rscript plot_sfig2.R $PEAKCALLER 'input'
    if [ "$PEAKCALLER" == "macs2" ]; then
        Rscript plot_sfig1g.R $MM10_GTF $HG38_GTF $PEAKCALLER exp_summary.txt
    fi
done
fi

## Step 3: estimate gene expression for Figures 1 and 2 using kallisto 
#this requires kallisto indexes for reference transcriptomes, fastq files, and estimated fragment lengths
#output is stored in subfolders labelled kallisto_results

#set locations of kallisto indexes for human and mouse
IDX_HG38=''
IDX_MM10=''

if $RUN_KALLISTO; then
while read p; do

    FIG=$(echo $p | cut -d' ' -f14)
    if [ "$FIG" != "fig" ] && ([ "$FIG" == "1" ] || [ "$FIG" == "4" ] || [ "$FIG" == "xiao" ] || [ "$FIG" == "2" ]); then
        DIR=$(echo $p | cut -d' ' -f1)
        COND=$(echo $p | cut -d' ' -f6)
        FRAGLEN=$(echo $p | cut -d' ' -f9)
        REPS=$(echo $p | cut -d' ' -f11)
        GENOME=$(echo $p | cut -d' ' -f3)
        if [ "$GENOME" == "mm10" ]; then
            IDX=$IDX_MM10
        elif [ "$GENOME" == "hg38" ]; then
            IDX=$IDX_HG38
        fi
        bash run_kallisto.sh $DIR $COND $FRAGLEN $IDX $REPS $NUM_THREADS $FIG
    fi
done < exp_summary.txt
fi

##Step 4: plot most subfigures for Figure 1 and prep files for Figure 2

#for Figure 2, option to modify the threshold for mean gene expression (recommend 10-50) and the peak caller (metdiff or macs2)
PEAKCALLER=macs2
THRESH=10

if $PLOT_FIG1; then
    while read p; do
    FIG=$(echo $p | cut -d' ' -f14)
    if [ "$FIG" = "1" ] || [ "$FIG" = "2" ] || [ "$FIG" = "xiao" ]; then
        DIR=$(echo $p | cut -d' ' -f1)
        LABEL=$(echo $p | cut -d' ' -f2)
        READLEN=$(echo $p | cut -d' ' -f8)
        FRAGLEN=$(echo $p | cut -d' ' -f9)
        PE=$(echo $p | cut -d' ' -f10)
        REPS=$(echo $p | cut -d' ' -f11)
        GENOME=$(echo $p | cut -d' ' -f3)
        if [ "$GENOME" == "mm10" ]; then
            IDX=$IDX_MM10
            GTF=$MM10_GTF
        elif [ "$GENOME" == "hg38" ]; then
            IDX=$IDX_HG38
            GTF=$HG38_GTF
        fi
        SAMPLE=$(echo $p | cut -d' ' -f 6)
        if [ "$FIG" = "1" ] || [ "$FIG" = "2" ]; then 
            bash plot_fig1a.sh $DIR $SAMPLE $READLEN $FRAGLEN $PE true $GTF $THRESH $LABEL
        else
            bash plot_fig1a.sh $DIR $SAMPLE $READLEN $FRAGLEN $PE false $GTF $THRESH $LABEL
        fi
        if [ "$FIG" = "1" ]; then 
            bash plot_fig1bc.sh $DIR $SAMPLE
            Rscript plot_fig1b.R $DIR $SAMPLE
        fi
    fi
    done < exp_summary.txt
    python plot_fig1a_summary.py #plots Supp Fig 1c
    #Rscript plot_sfig1a.R $MM10_GTF $HG38_GTF exp_summary.txt #need to install MeTPeak and exomePeak
fi

## Step 5: plot subfigures for Figure 2 (heatmaps)
if $PLOT_FIG2; then
    bash plot_fig2.sh exp_summary.txt $PEAKCALLER $THRESH $MM10_GTF $HG38_GTF #must run plot_fig1a.sh first
    Rscript plot_fig2.R $PEAKCALLER $THRESH exp_summary.txt 2
    if [ "$PEAKCALLER" != "metdiff" ]; then 
        Rscript plot_fig2.R $PEAKCALLER $THRESH exp_summary.txt xiao
    fi
fi

## Step 6: plot subfigures for Figure 3
RUN_DEQ=FALSE #change to TRUE to rerun DEQ (requires bam files)
if $PLOT_FIG3; then
    if [ ! -d fig3 ]; then mkdir fig3; fi
    Rscript plot_fig3.R $PEAKCALLER exp_summary.txt $MM10_GTF $HG38_GTF $RUN_DEQ
fi

## Step 7: plot subfigures for Figure 4
RUN_DEQ=FALSE #change to TRUE to rerun DEQ (requires bam files)
if $PLOT_FIG4; then
    if [ ! -d fig4 ]; then mkdir fig4; fi
    bash plot_fig4prep.sh exp_summary.txt
    Rscript plot_fig4.R $PEAKCALLER $THRESH exp_summary.txt $MM10_GTF $HG38_GTF $RUN_DEQ
    Rscript plot_fig4_mergedpeaks.R merged $THRESH exp_summary_mergedpeaks.txt $MM10_GTF $HG38_GTF $RUN_DEQ
fi

## Step 8: plot subfigures for Figure 5
if $PLOT_FIG5; then 
    if [ ! -d fig5 ]; then mkdir fig5; fi
    Rscript plot_fig5.R
fi
