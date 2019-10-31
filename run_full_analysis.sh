#!/bin/bash

NUM_THREADS=16 #set number of threads to run analyses
HG38_GTF='/athena/masonlab/scratch/users/abm237/REFS/gencode_beds_hg38/gencode.v25.annotation.gtf'
MM10_GTF='x'

#select steps and analyses to run
RUN_STAR=false #change to true to start analysis from star alignments (requires you have downloaded and named files in accordance with exp_summary_test.txt)
RUN_MACS2=true #run peak calling using MACS2 (requires alignments in bam files)
RUN_METDIFF=false #run peak calling using MeTDiff (requires alignments in bam files)
RUN_KALLISTO=false #run kallisto to estimate gene expression for Figure 1 (requires fastqs named appropriately)
PLOT_FIG1=false #plot Figure 1, requires results from RUN_STAR, RUN_KALLISTO, and ≥ 1 peak caller
PLOT_FIG2=true #plot Figure 2, requires results from RUN_STAR and ≥ 1 peak caller
PLOT_FIG3=false #plot Figure 3, requires results from RUN_STAR and ≥ 1 peak caller
PLOT_FIG4=false #plot Figure 4, required files all present in folder "fig4"
#see specific sections below for more details

#todo: plot_sfig1d.R + plot_sfig1a.R (requires fig2 neg control to test), plot_fig1cd, plot_fig2.R, plot_fig3.R

## Step 1: align data, produce bams

#run STAR alignments for paired end and single end data
#this requires the summary file exp_summary_test.txt, with directory names in the first column and species in the third column
#also requires Trimmomatic jar and adapter fastas - need to set paths in run_star.1.sh script 
#each directory must contain the raw fastq files for a different data set within a folder named "fastqs"
#(fastq names  should be formatted as ${condition}_${IP or input}_${rep #}_R${1 or 2}.fastq.gz, eg. stress_IP_1_R1.fastq.gz)
#this step will generate bam files within a folder called "alignments"

if $RUN_STAR; then 
#set locations of star indexes below for human and mouse
SI_HG38='/athena/masonlab/scratch/databases/hg38/star/'
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
done < exp_summary_test.txt
fi

## Step 2: call peaks 

#run MACS2 and/or MeTDiff to call peaks based on IP and input bam files
#this requires the successful completion of STAR above and bam files for each data set within a folder named "alignments"

PEAKCALLERS=()
if $RUN_MACS2 || $RUN_METDIFF; then
if $RUN_METDIFF; then
    #Rscript run_metdiff_peakcaller.R exp_summary_test.txt $NUM_THREADS $HG38_GTF $MM10_GTF
    PEAKCALLERS+=('metdiff')
fi
while read p; do
    DIR=$(echo $p | cut -d' ' -f1)
    if [ "$DIR" != "study" ]; then
    LABEL=$(echo $p | cut -d' ' -f2)
    CELLLINE=$(echo $p | cut -d' ' -f4)
    COND=$(echo $p | cut -d' ' -f6)
    YEAR=$(echo $p | cut -d' ' -f13)
    FRAGLEN=$(echo $p | cut -d' ' -f9)
    REPS=$(echo $p | cut -d' ' -f 11)

    if $RUN_MACS2; then
        PEAKCALLERS+=('macs2')
        #bash run_macs2.2.sh $DIR $FRAGLEN $REPS $NUM_THREADS
        bash plot_sfig1e.sh $DIR $CELLLINE $COND $YEAR $REPS $LABEL macs2
    fi
    if $RUN_METDIFF; then
        for FI in $DIR/*_MeTDiff_output/MeTDiff_output/peak.bed; do
            finame="${FI%.*}"
            bedtools sort -i $FI > $finame.sorted.bed
            echo $finame $finame.sorted.bed 
        done
        #bash plot_sfig1e.sh $DIR $CELLLINE $COND $YEAR $REPS $LABEL metdiff
    fi
    fi
done < exp_summary_test.txt
#plot peaks vs reads comparison and peak DRAC motif enrichment
for PEAKCALLER in "${PEAKCALLERS[@]}"; do
    Rscript plot_sfig1e.R $PEAKCALLER 'IP'
    Rscript plot_sfig1e.R $PEAKCALLER 'input'
    if [ "$PEAKCALLER" == "macs2" ]; then
        Rscript plot_sfig1d.R $MM10_GTF $HG38_GTF $PEAKCALLER exp_summary_test.txt
    fi
done
fi

## Step 3: estimate gene expression for Figure 1  using kallisto 

#set locations of kallisto indexes for human and mouse
IDX_HG38='/athena/masonlab/scratch/databases/hg38/kallisto/Homo_sapiens.GRCh38.rel79.cdna.all.kallidx'
IDX_MM10=''

if $RUN_KALLISTO; then
while read p; do

    FIG=$(echo $p | cut -d' ' -f14)
    if [ "$FIG" != "fig" ] && ([ "$FIG" == "1a" ] || [ "$FIG" == "3" ] || [ "$FIG" == "xiao" ] || [ "$FIG" == "1d" ]); then
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
done < exp_summary_test.txt
fi

##Step 4: plot most subfigures for Figure 1

#for Figure 1c-d, option to modify the threshold for mean gene expression (recommend 10-50) and the peak caller (metdiff or macs2)
PEAKCALLER=macs2
THRESH=10

if $PLOT_FIG1; then
while read p; do

    FIG=$(echo $p | cut -d' ' -f14)
    if [ "$FIG" = "1a" ] || [ "$FIG" = "1d" ] || [ "$FIG" = "xiao" ]; then
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
        #bash plot_fig1a.sh $DIR $SAMPLE $READLEN $FRAGLEN $PE true $GTF $THRESH $LABEL
        #bash plot_fig1b.sh $DIR $SAMPLE
        #Rscript plot_fig1b.R $DIR $SAMPLE
    fi
done < exp_summary_test.txt
#python plot_fig1a_summary.py
Rscript plot_sfig1a.R $MM10_GTF $HG38_GTF

#bash plot_fig1cd.sh $PEAKCALLER $THRESH $MM10_GTF $HG38_GTF
#Rscript plot_fig1cd.R
fi

## Step 5: plot subfigures for Figure 2

if $PLOT_FIG2; then
while IFS= read -r p; do
    DIR=$(echo $p | cut -d' ' -f1)
for PEAKCALLER in "${PEAKCALLERS[@]}"; do
    if [ "$PEAKCALLER" == "macs2" ]; then
    outdir=$DIR/macs2_results
    if [ -d $outdir ] ; then
        echo $outdir
        for fi in $outdir/*.narrowPeak; do awk -v OFS="\t" '($4>=1){print $1,$2,$3}' $fi > $fi.peaks.bed; done
    fi
    fi
done
done < exp_summary_test.txt
Rscript plot_fig2.R $PEAKCALLER $THRESH
fi

## Step 6: plot subfigures for Figure 3

if $PLOT_FIG3; then
bash plot_fig3prep.sh
Rscript plot_fig3.R $PEAKCALLER $THRESH 
fi

## Step 7: plot subfigures for Figure 4
if $PLOT_FIG4; then 
Rscript plot_fig4.R
fi
