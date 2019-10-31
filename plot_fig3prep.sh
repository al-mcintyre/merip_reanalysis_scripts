min=1

while IFS= read -r line; do
    dir=$(echo $line | cut -d' ' -f 1)

    for ind in 4 5 7; do

        cells=$(echo $line | cut -d' ' -f $ind)
        if [ -d $dir ] ; then 
            echo $dir $cells

            m6a_bed=$dir/macs2_results/${cells}_min$min.bed
            set -- $dir/macs2_results/${cells}*_2.macs2_peaks.narrowPeak
            if [ -f "$1" ] ; then # "$dir/macs2_results/${cells}*_2.macs2_peaks.narrowPeak" ]; then 
                echo '>=2 reps'
                bedtools multiinter -i $dir/macs2_results/${cells}*_*.macs2_peaks.narrowPeak | awk -v OFS="\t" -v thresh=$min '($4>=thresh){print $1,$2,$3}' | bedtools merge -i - > $m6a_bed
            else
                echo '1 rep'
                cp $dir/macs2_results/${cells}*_1.macs2_peaks.narrowPeak $m6a_bed
            fi

            for fi in $dir/macs2_results/*.narrowPeak; do awk -v OFS="\t" '($4>=1){print $1,$2,$3}' $fi > $fi.peaks.bed; done
        fi
    done
done < fig3_experiment_summary.txt #control_summary.txt
