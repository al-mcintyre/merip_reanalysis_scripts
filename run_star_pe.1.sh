#script to run STAR for paired end reads
READDIR='' #specify experiment directory (reads should be stored within a directory inside called "fastqs" and should be formatted as ${condition}_${IP or input}_${rep #}_R${1 or 2}.fastq.gz, eg. stress_IP_1_R1.fastq.gz -- or just modify the below) 
STAR_INDEX='' #specify location of star index for the appropriate genome
cd $READDIR
NUM_THREADS=16 #as you like

mkdir alignments
for fi in fastqs/*IP_1*R1.fastq.gz; do
    COND=$(basename $fi)
    COND=${COND%_IP*}
    COND=${COND%_input*}
    for EXP in input IP; do   
        for REP in 1 2 3 4 5 6 7; do
            f=${COND}_${EXP}_${REP}
            PRETRIM1=fastqs/${f}_R1.fastq.gz
            PRETRIM2=fastqs/${f}_R2.fastq.gz
            FASTQ1=fastqs/${f}_R1.trim.fastq.gz
            FASTQ2=fastqs/${f}_R2.trim.fastq.gz
            ls $PRETRIM1 
            if [ -f $PRETRIM1 ]; then
                if [ ! -f $FASTQ1 ]; then
                    java -jar ~/bin/Trimmomatic-0.36/trimmomatic-0.36.jar PE -phred33 $PRETRIM1 $PRETRIM2 $FASTQ1 fastqs/${f}_R1_unpaired.fq.gz $FASTQ2 fastqs/${f}_R2_unpaired.fq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:15
                fi
                if [ -f fastqs/${f}_R1.trim.fastq.gz ]; then
                    if [ ! -f alignments/${f}.star.sorted.bam.bai ]; then
                      if [ ! -f ${f}.star.sorted.bam.bai ]; then
                        STAR --genomeDir $STAR_INDEX --runThreadN $NUM_THREADS --readFilesIn fastqs/${f}_R1.fastq.gz fastqs/${f}_R2.fastq.gz --outFilterMultimapNmax 1 --readFilesCommand zcat --outFileNamePrefix ${f}.star. --genomeLoad LoadAndKeep  --outSAMstrandField intronMotif
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
