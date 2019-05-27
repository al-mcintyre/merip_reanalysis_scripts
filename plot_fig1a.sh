dir=$1 #neuron_engel #neg_control_mocks #neuron_engel 
cells=$2 #basal #mock #basal 
genome=$3 #mm10 #hg38 #mm10
readlen=$4 #200 #50 #200 #50 #200 #read length for Engel paper ~100 bp per end, but sometimes completely overlapped
plot=$5

min=1
m6a_bed=$dir/macs2_results/${cells}_min$min.bed

bash run_kallisto.sh $dir $cells $genome

awk -v OFS='\t' '{name[FNR]=$1;len[FNR]=$2;eff[FNR]=$3;a[FNR]+=$4;b[FNR]+=$5;c[FNR]++;}END{for(i=1;i<=FNR;i++)print name[i],len[i],eff[i],a[i]/c[i],b[i]/c[i];}' $dir/kallisto_results/${cells}*/abundance.tsv >  $dir/kallisto_results/${cells}_kallisto_mean_abundance.tsv

python get_most_highly_expressed.py $dir/kallisto_results/${cells}_kallisto_mean_abundance.tsv $genome $readlen True
set -- $dir/macs2_results/${cells}*_2.macs2_peaks.narrowPeak
if [ -f "$1" ] ; then # "$dir/macs2_results/${cells}*_2.macs2_peaks.narrowPeak" ]; then 
    echo '>=2 reps'
    bedtools multiinter -i $dir/macs2_results/${cells}*_*.macs2_peaks.narrowPeak | awk -v OFS="\t" -v thresh=$min '($4>=thresh){print $1,$2,$3}' | bedtools merge -i - > $m6a_bed
else 
    echo '1 rep'
    cp $dir/macs2_results/${cells}*_1.macs2_peaks.narrowPeak $m6a_bed
fi
wc -l $m6a_bed

if $plot ; then
    bedtools intersect -a $dir/kallisto_results/${cells}_highest_expression_transcripts_all.gtf -b $m6a_bed -wa -c > $dir/kallisto_results/${cells}_highest_expression_transcripts_counts.gtf 
    python plot_counts_vs_expression.py $dir/kallisto_results/${cells}_highest_expression_transcripts_counts.gtf $readlen
fi 
