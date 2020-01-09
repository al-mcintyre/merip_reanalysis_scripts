DIR=$1
CELLLINE=$2
COND=$3
YEAR=$4
REPS=$5
LABEL=$6
PEAKCALLER=$7

for FRACTION in IP input; do
    OUTFI=fig2/reads_vs_${PEAKCALLER}_peaks_${FRACTION}.txt
    echo 'here'
    echo $OUTFI
    if [ ! -f "$OUTFI" ]; then
        echo "label year cell_line rep reads peaks" > $OUTFI
    fi
    echo $LABEL
    for REP in $(seq 1 $REPS); do
        bam=${DIR}/alignments/${COND}_${FRACTION}_${REP}.star.sorted.bam
        if [ -f $bam  ]; then
            num_reads=$(samtools flagstat $bam | grep mapped | cut -d' ' -f 1 | head -1)
            if [ "$PEAKCALLER" == "macs2" ]; then
                num_peaks=$(wc -l $DIR/macs2_results/${COND}_${REP}.macs2_summits.bed | cut -d' ' -f 1)
            else
                num_peaks=$(wc -l $DIR/${COND}_MeTDiff_output/MeTDiff_output/peak.bed | cut -d' ' -f 1)
            fi
            echo ${CELLLINE}_${LABEL} $YEAR $COND $REP $num_reads $num_peaks >> $OUTFI
        fi
    done
done
