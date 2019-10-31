#requires trimmomatic.jar and adapter files in the same directory (or replace locations)

READDIR=$1
REPS=$2
STAR_INDEX=$3
NUM_THREADS=$4

cd $READDIR

if [ ! -d alignments ]; then
    mkdir alignments
fi
for fi in fastqs/*IP_1*.fastq.gz; do
    COND=$(basename $fi)
    COND=${COND%_IP*}
    
    for EXP in input IP; do   
        for REP in $(seq 1 $REPS); do
            f=${COND}_${EXP}_${REP}
            echo $f

            PRETRIM=fastqs/${f}.fastq.gz
            PRETRIM1=fastqs/${f}_R1.fastq.gz
            PRETRIM2=fastqs/${f}_R2.fastq.gz
            FASTQ1=fastqs/${f}_R1.trim.fastq.gz
            FASTQ2=fastqs/${f}_R2.trim.fastq.gz
            
            if [ -f $PRETRIM1 ] || [ -f $FASTQ1 ] ; then
                if [ ! -f $FASTQ1 ]; then
                    java -jar trimmomatic-0.36.jar PE -phred33 $PRETRIM1 $PRETRIM2 $FASTQ1 fastqs/${f}_R1_unpaired.fq.gz $FASTQ2 fastqs/${f}_R2_unpaired.fq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:15
                fi
                if [ -f fastqs/${f}_R1.trim.fastq.gz ]; then
                    if [ ! -f alignments/${f}.star.sorted.bam.bai ] && [ ! -f ${f}.star.sorted.bam.bai ]; then
                        STAR --genomeDir $STAR_INDEX --runThreadN $NUM_THREADS --readFilesIn fastqs/${f}_R1.fastq.gz fastqs/${f}_R2.fastq.gz --outFilterMultimapNmax 1 --readFilesCommand zcat --outFileNamePrefix ${f}.star. --genomeLoad LoadAndKeep  --outSAMstrandField intronMotif
                    fi
                fi
            fi
            if [ -f $PRETRIM ] || [ -f fastqs/$f.trim.fastq.gz ]; then
                if [ ! -f fastqs/$f.trim.fastq.gz ]; then
                    java -jar trimmomatic-0.36.jar SE -phred33 -threads $NUM_THREADS fastqs/$f.fastq.gz fastqs/$f.trim.fastq ILLUMINACLIP:TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:15
                    gzip fastqs/$f.trim.fastq
                fi
                if [ -f fastqs/$f.trim.fastq.gz ]; then
                    if [ ! -f alignments/${f}.star.sorted.bam.bai ] && [ ! -f ${f}.star.sorted.bam.bai ]; then
                        STAR --genomeDir $STAR_INDEX --runThreadN $NUM_THREADS --readFilesIn fastqs/${f}.trim.fastq.gz --outFilterMultimapNmax 1 --readFilesCommand zcat --outFileNamePrefix ${f}.star. --genomeLoad LoadAndKeep  --outSAMstrandField intronMotif
                    fi
                fi
            fi
            if [ -f ${COND}_${EXP}_${REP}.star.Aligned.out.sam ]; then 
                samtools view -b ${COND}_${EXP}_${REP}.star.Aligned.out.sam > ${COND}_${EXP}_${REP}.star.bam
                rm ${COND}_${EXP}_${REP}.star.Aligned.out.sam
                samtools sort -o ${COND}_${EXP}_${REP}.star.sorted.bam ${COND}_${EXP}_${REP}.star.bam
                rm ${COND}_${EXP}_${REP}.star.bam
                samtools index ${COND}_${EXP}_${REP}.star.sorted.bam
            fi
            
        done
    done
done

mv *.bam* alignments/
mv *.star* alignments/
STAR --genomeDir ${STAR_INDEX} --genomeLoad Remove --outFileNamePrefix star_genome_remove
rm -r star_genome_remove*
