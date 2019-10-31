#uses gnu parallel (sem) to run multiple processes for each directory

DIR=$1
FRAGLEN=$2
REPS=$3
NUM_THREADS=$4

cd $DIR
if [ ! -d macs2_results ]; then
    mkdir macs2_results
fi
for fi in fastqs/*IP*_1*trim*.gz; do
    COND=$(basename $fi)
    COND=${COND%_IP*}
    for REP in $(seq 1 $REPS); do
        if [ -f alignments/${COND}_input_${REP}.star.sorted.bam ]; then
            #transcriptome size estimate based on ucsc table browser summary for gencode v26 hg38 and mm10
            if [ ! -f macs2_results/${COND}_${REP}.macs2_peaks.narrowPeak ]; then
                sem -j$NUM_THREADS macs2 callpeak -t alignments/${COND}_IP_${REP}.star.sorted.bam -c alignments/${COND}_input_${REP}.star.sorted.bam --nomodel --extsize $FRAGLEN -g 100e6 -n macs2_results/${COND}_${REP}.macs2 -f BAM --verbose 3
                echo 'Calling peaks for ' $COND ' rep ' $REP
            fi
        fi
    done
done
sem --wait
cd ..
