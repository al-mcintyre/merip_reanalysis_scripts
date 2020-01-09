min=1

while IFS= read -r line; do
    dir=$(echo $line | cut -d' ' -f 1)

    fig=$(echo $line | cut -d' ' -f 14)
    if [ "$fig" == "4" ] || [ "$fig" == "kshv" ] || [ "$fig" == "hiv" ] || [ "$fig" == "hcmv" ]; then
        for ind in 6 7 ; do #take sixth and seventh columns from experiment summary (baseline and stimulus conditions)

        cond=$(echo $line | cut -d' ' -f $ind)
        if [ -d $dir ] ; then 
            #echo $dir $cond
            #m6a_bed=$dir/macs2_results/${cond}_min$min.bed
            #set -- $dir/macs2_results/${cond}*_2.macs2_peaks.narrowPeak
            #bedtools multiinter -i $dir/macs2_results/*${cond}_*.macs2_peaks.narrowPeak | awk -v OFS="\t" -v thresh=$min '($4>=thresh){print $1,$2,$3}' | bedtools merge -i - > $m6a_bed
            for fi in $dir/macs2_results/*.narrowPeak; do awk -v OFS="\t" '($4>=1){print $1,$2,$3}' $fi > $fi.peaks.bed; done
        fi
      done
  fi
done < $1
