#for SE reads, assuming illumina adapters add ~120 bases to bioanalyzer trace - guessing here total of 180 would be approximately correct (50 bp reads)
indir=$1
label=$2
genome=$3

if [ "$genome" == "mm10" ]; then
   kallindex=/athena/masonlab/scratch/users/abm237/REFS/gencode_mm10/mm10.kallisto.idx
elif [ "$genome" == "hg38" ]; then
   kallindex=/athena/masonlab/scratch/databases/hg38/kallisto/Homo_sapiens.GRCh38.rel79.cdna.all.kallidx
   #kallindex=/athena/masonlab/scratch/users/abm237/REFS/gencode_hg38/hg38.kmer29.kallisto.idx
fi
cd $indir
echo $indir
if [ ! -d kallisto_results ]; then 
    mkdir kallisto_results
fi
for fi in fastqs/$label*input*trim*.gz; do
    COND=$(basename $fi)
    COND=${COND%_IP*}
    COND=${COND%_input*}
    echo $fi $COND
    for REP in 1 2 3 4 5 6 7 ; do
        if [ -f fastqs/${COND}_input_${REP}.fastq.gz ]; then
            sample=${COND}_${REP}
            if [ ! -d kallisto_results/$sample.kallisto ]; then
                echo fastqs/${COND}_input_${REP}.fastq.gz
                kallisto quant -i $kallindex -o kallisto_results/$sample.kallisto --single -l 180 -s 20 fastqs/${COND}_input_${REP}.trim.fastq.gz
            fi
        fi
        if [ -f fastqs/${COND}_input_${REP}_R1.fastq.gz ]; then
            sample=${COND}_${REP}
            if [ ! -d kallisto_results/$sample.kallisto ]; then
                echo fastqs/${COND}_input_${REP}_R1.trim.fastq.gz
                kallisto quant -i $kallindex -o kallisto_results/$sample.kallisto fastqs/${COND}_input_${REP}_R1.trim.fastq.gz fastqs/${COND}_input_${REP}_R2.trim.fastq.gz 
            fi
        fi 
    done
done
cd ..
