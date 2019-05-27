if [ ! -d peak_overlaps ]; then 
    mkdir peak_overlaps
fi

infi=fig1cd_exp_summary.txt
summary=peak_overlaps/peak_overlaps_summary.txt
overwrite=false
min_reps=1
if [ ! -f "$summary" ]; then
    echo "cellline1 cellline2 exp1 exp2 gene_overlap min_reps peak_overlap total_peaks" > $summary
fi

while IFS= read -r line1; do
    echo $line1
    cellline1=$(echo $line1 | cut -d' ' -f1)
    dir1=$(echo $line1 | cut -d' ' -f2 )
    exp1=$(echo $line1 | cut -d' ' -f3)
    species1=$(echo $line1 | cut -d' ' -f4)
    sample1=$(echo $line1 | cut -d' ' -f5)
    label1=${cellline1}_${exp1}
    readlen1=$(echo $line1 | cut -d' ' -f6)

    while read -r line2; do
        cellline2=$(echo $line2 | cut -d' ' -f1)
        dir2=$(echo $line2 | cut -d' ' -f2)
        exp2=$(echo $line2 | cut -d' ' -f3)
        species2=$(echo $line2 | cut -d' ' -f4)
        sample2=$(echo $line2 | cut -d' ' -f5)
        label2=${cellline2}_${exp2}
        readlen2=$(echo $line2 | cut -d' ' -f6)

        if [ "$label1" != "$label2" ] && ( [ "$species1" == "$species2" ] && ( $overwrite || ( [ ! -f peak_overlaps/${label1}_${label2}_1.bed ] && [ ! -f peak_overlaps/${label2}_${label1}_1.bed ] ) ) ) ; then

            echo $label1 $label2

            bash plot_fig1a.sh $dir1 $sample1 $species1 $readlen1 false
            bash plot_fig1a.sh $dir2 $sample2 $species2 $readlen2 false

            comm <(sort $dir1/kallisto_results/${sample1}_transcripts_over50.gtf) <(sort $dir2/kallisto_results/${sample2}_transcripts_over50.gtf) -12 > peak_overlaps/${label1}_${label2}.gtf
            ngenes=$(cut -d'"' -f 2 peak_overlaps/${label1}_${label2}.gtf | sort | uniq | wc -l)
            bedtools intersect -a $dir1/macs2_results/${sample1}_min${min_reps}.bed -b peak_overlaps/${label1}_${label2}.gtf -wa | bedtools merge -i - > peak_overlaps/${label1}_${label2}_1.bed
            bedtools intersect -a $dir2/macs2_results/${sample2}_min${min_reps}.bed -b peak_overlaps/${label1}_${label2}.gtf -wa | bedtools merge -i - > peak_overlaps/${label1}_${label2}_2.bed
            peaks1=$(wc -l peak_overlaps/${label1}_${label2}_1.bed | cut -d' ' -f1 )
            peaks2=$(wc -l peak_overlaps/${label1}_${label2}_2.bed | cut -d' ' -f1 )
            overlap1=$(bedtools intersect -a peak_overlaps/${label1}_${label2}_1.bed -b peak_overlaps/${label1}_${label2}_2.bed -wa | bedtools merge -i - | wc -l) #or use sort | uniq since should all be wa entries
            overlap2=$(bedtools intersect -a peak_overlaps/${label1}_${label2}_2.bed -b peak_overlaps/${label1}_${label2}_1.bed -wa | bedtools merge -i - | wc -l)
            echo ${cellline1} ${cellline2} ${exp1} ${exp2} $ngenes ${min_reps} $overlap1 $peaks1 >> $summary
            echo ${cellline2} ${cellline1} ${exp2} ${exp1} $ngenes ${min_reps} $overlap2 $peaks2 >> $summary

        fi
    done < $infi
done < $infi
