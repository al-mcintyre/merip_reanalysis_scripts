#script to run STAR for single end reads
READDIR='' #specify experiment directory (reads should be stored within a directory inside called "fastqs" and should be formatted as ${condition}_${IP or input}_${rep #}.fastq.gz, eg. stress_IP_1.fastq.gz -- or just modify the below) 
STAR_INDEX='' #specify location of star index for the appropriate genome
cd $READDIR
NUM_THREADS=16 #as you like

mkdir alignments
for fi in fastqs/*IP_[1-9].fastq.gz; do
    COND=$(basename $fi)
    COND=${COND%_IP*}
    COND=${COND%_input*}
    for EXP in input IP; do   
        for REP in 1 2 3 4 5 6 7; do
            f=${COND}_${EXP}_${REP}
            ls fastqs/${f}.fastq.gz 
            if [ -f fastqs/${f}.fastq.gz ]; then
                if [ ! -f fastqs/$f.trim.fastq.gz ]; then
                    java -jar ~/bin/Trimmomatic-0.36/trimmomatic-0.36.jar SE -phred33 -threads $NUM_THREADS fastqs/$f.fastq.gz fastqs/$f.trim.fastq ILLUMINACLIP:TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:15
                    gzip fastqs/$f.trim.fastq
                fi
                if [ -f fastqs/$f.trim.fastq.gz ]; then
                    if [ ! -f alignments/${f}.star.sorted.bam.bai ]; then
                      if [ ! -f ${f}.star.sorted.bam.bai ]; then
                        STAR --genomeDir $STAR_INDEX --runThreadN $NUM_THREADS --readFilesIn fastqs/${f}.trim.fastq.gz --outFilterMultimapNmax 1 --readFilesCommand zcat --outFileNamePrefix ${f}.star. --genomeLoad LoadAndKeep  --outSAMstrandField intronMotif

                        samtools view -b ${COND}_${EXP}_${REP}.star.Aligned.out.sam > ${COND}_${EXP}_${REP}.star.bam 
                        rm ${COND}_${EXP}_${REP}.star.Aligned.out.sam
                        samtools sort -o ${COND}_${EXP}_${REP}.star.sorted.bam ${COND}_${EXP}_${REP}.star.bam
                        rm ${COND}_${EXP}_${REP}.star.bam
                        samtools index ${COND}_${EXP}_${REP}.star.sorted.bam
                      fi 
                    fi
                fi
            fi
        done
    done
done

mv *.bam* alignments/
mv *.star* alignments/
STAR --genomeDir ${STAR_INDEX} --genomeLoad Remove --outFileNamePrefix star_genome_remove
