DIR=$1
COND=$2
min=1
m6a_bed=$DIR/macs2_results/${COND}_min$min.bed

outdir=$DIR/macs2_results

#get set of all peaks from all reps
for fi in $outdir/*.narrowPeak; do awk -v OFS="\t" '($4>=1){print $1,$2,$3}' $fi > $fi.peaks.bed; done
bedtools multiinter -i $outdir/${COND}*.macs2_peaks.narrowPeak.peaks.bed | awk -v OFS="\t" '($4>=1){print $1,$2,$3}' | bedtools merge -i - > $outdir/${COND}_replicate_intersect.bed
cp $outdir/${COND}_replicate_intersect.bed $outdir/${COND}_single_rep.bed.tmp

#intersect each rep with the full set and save intersect file 
for fi in $outdir/${COND}*narrowPeak.peaks.bed; do echo $fi; bedtools intersect -a $outdir/${COND}_single_rep.bed.tmp -b $fi -wa -c | cut -f 4 > intersect.tmp; paste $outdir/${COND}_replicate_intersect.bed intersect.tmp > $outdir/${COND}_replicate_intersect.tmp; mv $outdir/${COND}_replicate_intersect.tmp $outdir/${COND}_replicate_intersect.bed; done #head $outdir/${COND}_replicate_intersect.bed; done
rm *.tmp
python plot_fig1c_replicate_overlap.py $DIR ${COND} 

num_reps=$(head -n 1 $outdir/${COND}_replicate_intersect.bed | awk '{print NF-3}')
echo $num_reps
for n in $(seq $num_reps); do 
    echo $n
    awk -v nreps=$n -v OFS="\t" '{for(i=4; i<=NF; i++) if($i >=1) j+=1; if (j == nreps) {print $1,$2,$3}; j=0; i=4}' $outdir/${COND}_replicate_intersect.bed > $outdir/${COND}_replicate_intersect_${n}.bed
done
