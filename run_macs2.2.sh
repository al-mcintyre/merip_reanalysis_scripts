summary=fig3_experiment_summary.txt

while read p; do
#dir=$1
#FRAGLEN=$2 8th col
dir=$(echo $p | cut -d' ' -f1) #2 for fig1
FRAGLEN=$(echo $p | cut -d' ' -f8)
echo $dir $FRAGLEN
if [ "$dir" != "directory" ]; then
    if [ "$FRAGLEN" -ne 100 ]; then
cd $dir
mkdir macs2_results2
for fi in fastqs/basal*.gz; do
    COND=$(basename $fi)
    COND=${COND%_IP*}
    COND=${COND%_input*}
    echo $fi $COND
    for REP in 1 2 3 4 5 6 7; do
        echo alignments/${COND}_IP_${REP}.star.sorted.bam
        if [ -f alignments/${COND}_input_${REP}.star.sorted.bam ]; then
            #transcriptome size estimate based on ucsc table browser summary for gencode v 26, hg38 
            if [ ! -f macs2_results2/${COND}_${REP}.macs2_peaks.narrowPeak ]; then
                echo macs2_results2/${COND}_${REP}.macs2
                macs2 callpeak -t alignments/${COND}_IP_${REP}.star.sorted.bam -c alignments/${COND}_input_${REP}.star.sorted.bam --nomodel --extsize $FRAGLEN -g 100e6 -n macs2_results2/${COND}_${REP}.macs2 -f BAM --verbose 3
            fi
        fi
    done
done
cd ..
fi
fi

done <$summary
