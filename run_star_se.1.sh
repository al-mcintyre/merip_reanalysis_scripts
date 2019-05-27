#!/bin/bash -l
#$ -N star
#$ -j y
#$ -cwd
#$ -pe smp 8
#$ -l h_vmem=5G
#$ -l h_rt=24:00:0
#$ -l athena=True

#unset LD_LIBRARY_PATH
#spack load star@2.6.1a%gcc@6.3.0
spack load samtools@1.9%gcc@6.3.0

READDIR=huang_histonemeth_2019 #zhang_alkbh5_2017 #rubio_hcmv_2018 #li_ythdf2_2018 #liu_endometrialtumours_2018 #su_r2hg_2018 #lichinchi_hiv_2016 #cerchietti_heatshock_unpublished #anders_stress_2017 #lichinchi_hiv_2016 #tan_kshv_2018 #hess_ftoko_2013 #barbieri_leukemia_2017 #geula_development_2015 #li_leukemia_2017 #fustin_circadian_2013 #tan_kshv_2018 #geula_development_2015 #tirumuru_hiv_2016 #geula_development_2015 #hesser_kshv_2018 #neg_control_mocks #zhou_heatshock_2015 #neg_control_mocks #mettl3_yu #stresses_dominissini #original_meyer #neg_control_mocks #mirna_chen #mettl314_wang #mettl3_batista #zikv_lichinchi2 #translation_lin #mettl14_huang #original_meyer #neg_control_mocks #$1
GENOME=hg38 #mm10 #hg38 #mm10 #hg38 #$2
PREFIX='' #brain #hek293T
cp TruSeq3-SE.fa $READDIR
cd $READDIR
if [ "$GENOME" == "mm10" ]; then
    STAR_INDEX=/athena/masonlab/scratch/users/abm237/REFS/mm10_star/
elif [ "$GENOME" == "hg38" ]; then 
    STAR_INDEX=/athena/masonlab/scratch/databases/hg38/star/
fi
NUM_THREADS=16
echo $STAR_INDEX 

mkdir alignments
for fi in fastqs/${PREFIX}*IP_[123].fastq.gz; do
    COND=$(basename $fi)
    COND=${COND%_IP*}
    COND=${COND%_input*}
    for EXP in input IP; do   
        for REP in 1 2 3; do
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
