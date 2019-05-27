dir=engel_neuron_2018 #neg_control_mocks  #neuron_engel 
cells=basal #mock #basal 
min=1
m6a_bed=$dir/macs2_results/${cells}_min$min.bed

outdir=$dir/macs2_results

#get set of all peaks from all reps
for fi in $outdir/*.narrowPeak; do awk -v OFS="\t" '($4>=1){print $1,$2,$3}' $fi > $fi.peaks.bed; done
bedtools multiinter -i $outdir/${cells}*.macs2_peaks.narrowPeak.peaks.bed | awk -v OFS="\t" '($4>=1){print $1,$2,$3}' | bedtools merge -i - > $outdir/${cells}_replicate_intersect.bed
cp $outdir/${cells}_replicate_intersect.bed $outdir/${cells}_single_rep.bed.tmp

#intersect each rep with the full set and save intersect file 
for fi in $outdir/${cells}*narrowPeak.peaks.bed; do echo $fi; bedtools intersect -a $outdir/${cells}_single_rep.bed.tmp -b $fi -wa -c | cut -f 4 > intersect.tmp; paste $outdir/${cells}_replicate_intersect.bed intersect.tmp > $outdir/${cells}_replicate_intersect.tmp; mv $outdir/${cells}_replicate_intersect.tmp $outdir/${cells}_replicate_intersect.bed; head $outdir/${cells}_replicate_intersect.bed; done
rm *.tmp
python plot_replicate_overlap_fig1b.py $dir ${cells} 

num_reps=$(head -n 1 $outdir/${cells}_replicate_intersect.bed | awk '{print NF-3}')
echo $num_reps
for n in $(seq $num_reps); do 
    echo $n
    awk -v nreps=$n -v OFS="\t" '{for(i=4; i<=NF; i++) if($i >=1) j+=1; if (j == nreps) {print $1,$2,$3}; j=0; i=4}' $outdir/${cells}_replicate_intersect.bed > $outdir/${cells}_replicate_intersect_${n}.bed
done
