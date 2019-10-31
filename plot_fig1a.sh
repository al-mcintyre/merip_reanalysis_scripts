DIR=$1 
COND=$2 
READLEN=$3
FRAGLEN=$4
PE=$5
PLOT=$6
GTF=$7
THRESH=$8
LABEL=$9

min=1
m6a_bed=$DIR/macs2_results/${COND}_min$min.bed

awk -v OFS='\t' '{name[FNR]=$1;len[FNR]=$2;eff[FNR]=$3;a[FNR]+=$4;b[FNR]+=$5;c[FNR]++;}END{for(i=1;i<=FNR;i++)print name[i],len[i],eff[i],a[i]/c[i],b[i]/c[i];}' $DIR/kallisto_results/${COND}*/abundance.tsv >  $DIR/kallisto_results/${COND}_kallisto_mean_abundance.tsv

python get_most_highly_expressed.py $DIR/kallisto_results/${COND}_kallisto_mean_abundance.tsv $READLEN True $GTF $THRESH
set -- $DIR/macs2_results/${COND}*_2.macs2_peaks.narrowPeak
if [ $min -gt 1 ] && [ -f "$DIR/macs2_results/${COND}*_${min}.macs2_peaks.narrowPeak" ] ; then
    echo '>=2 reps'
    bedtools multiinter -i $DIR/macs2_results/${COND}*_*.macs2_peaks.narrowPeak | awk -v OFS="\t" -v thresh=$min '($4>=thresh){print $1,$2,$3}' | bedtools merge -i - > $m6a_bed
else 
    echo '1 rep'
    cp $DIR/macs2_results/${COND}*_1.macs2_peaks.narrowPeak $m6a_bed
fi
wc -l $m6a_bed

if $PLOT ; then
    bedtools intersect -a $DIR/kallisto_results/${COND}_highest_expression_transcripts_all.gtf -b $m6a_bed -wa -c > $DIR/kallisto_results/${COND}_highest_expression_transcripts_counts.gtf 
    python plot_fig1a_counts_vs_expression.py $DIR/kallisto_results/${COND}_highest_expression_transcripts_counts.gtf $READLEN $COND $PE $FRAGLEN $LABEL
fi 
