macsdir=/pbtech_mounts/cmlab_store006/abm237/bin/anaconda2-4.1.1/lib/python2.7/site-packages/MACS2

dir=$1
cd $dir
mkdir macs2_results
for fi in fastqs/*.gz; do
    COND=$(basename $fi)
    COND=${COND%_IP*}
    COND=${COND%_input*}
    echo $fi $COND
    for REP in 1 2 3 4 5 6; do
        echo alignments/${COND}_IP_${REP}.star.sorted.bam
        if [ -f alignments/${COND}_input_${REP}.star.sorted.bam ]; then
            #transcriptome size estimate based on ucsc table browser summary for gencode v 26, hg38 
            if [ ! -f macs2_results/${COND}_${REP}.macs2_peaks.narrowPeak ]; then
                echo macs2_results/${COND}_${REP}.macs2
                macs2 callpeak -t alignments/${COND}_IP_${REP}.star.sorted.bam -c alignments/${COND}_input_${REP}.star.sorted.bam --nomodel --extsize 100 -g 100e6 -n macs2_results/${COND}_${REP}.macs2 -f BAM --verbose 3
            fi
        fi
    done
done
