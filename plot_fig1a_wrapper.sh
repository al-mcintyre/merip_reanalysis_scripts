if [ ! -d peak_overlaps ]; then 
    mkdir peak_overlaps
fi

infi=fig1a_exp_summary.txt
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

    bash plot_fig1a.sh $dir1 $sample1 $species1 $readlen1 true
done < $infi
