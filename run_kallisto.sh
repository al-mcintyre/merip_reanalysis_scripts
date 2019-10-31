#for SE reads, assuming illumina adapters add ~120 bases to bioanalyzer trace - 
DIR=$1
COND=$2
FRAGLEN=$3
IDX=$4
REPS=$5
NUM_THREADS=$6
FIG=$7

cd $DIR
echo $DIR
if [ ! -d kallisto_results ]; then 
    mkdir kallisto_results
fi

if [ "$FIG" == "xiao" ]; then
    COND=$(echo $COND | cut -d'_' -f1)
    REP=$(echo $COND | cut -d'_' -f2)
fi

for fi in fastqs/$COND*input*trim*.gz; do
    echo $fi $IDX
    for REP in $(seq 1 $REPS); do
        if [ -f fastqs/${COND}_input_${REP}.trim.fastq.gz ]; then
            SAMPLE=${COND}_${REP}
            if [ ! -d kallisto_results/$SAMPLE.kallisto ]; then
                echo 'running kallisto for ' fastqs/${COND}_input_${REP}.fastq.gz
                kallisto quant -t $NUM_THREADS -i $IDX -o kallisto_results/$SAMPLE.kallisto --single -l $FRAGLEN -s 20 fastqs/${COND}_input_${REP}.trim.fastq.gz 
            fi
        fi
        if [ -f fastqs/${COND}_input_${REP}_R1.trim.fastq.gz ]; then
            SAMPLE=${COND}_${REP}
            if [ ! -d kallisto_results/$SAMPLE.kallisto ]; then
                echo 'running kallisto for ' fastqs/${COND}_input_${REP}_R1.trim.fastq.gz
                kallisto quant -t $NUM_THREADS -i $IDX -o kallisto_results/$SAMPLE.kallisto fastqs/${COND}_input_${REP}_R1.trim.fastq.gz fastqs/${COND}_input_${REP}_R2.trim.fastq.gz 
            fi
        fi 
    done
done
cd ..
