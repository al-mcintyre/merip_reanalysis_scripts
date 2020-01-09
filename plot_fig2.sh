INFI=$1
PEAKCALLER=$2
THRESH=$3
MM10_GTF=$4
HG38_GTF=$5
outdir=fig2/peak_overlaps_${PEAKCALLER}_$THRESH
if [ ! -d $outdir ]; then
    mkdir $outdir
fi
summary=fig2/peak_overlaps_summary_${PEAKCALLER}_${THRESH}.txt
overwrite=false
min_reps=1
if [ ! -f "$summary" ]; then
    echo "cellline1 cellline2 exp1 exp2 gene_overlap min_reps peak_overlap total_peaks" > $summary
fi

while IFS= read -r line1; do
    echo $line1
    cellline1=$(echo $line1 | cut -d' ' -f4)
    dir1=$(echo $line1 | cut -d' ' -f1)
    label1=$(echo $line1 | cut -d' ' -f2)
    species1=$(echo $line1 | cut -d' ' -f3)
    sample1=$(echo $line1 | cut -d' ' -f6)
    readlen1=$(echo $line1 | cut -d' ' -f8)
    fraglen1=$(echo $line1 | cut -d' ' -f9)
    fig1=$(echo $line1 | cut -d' ' -f14)
    barcode1=${label1}_${sample1}

    if [ "$fig1" == "2" ] || [ "$fig1" == "xiao" ] ; then
      while read -r line2; do
        cellline2=$(echo $line2 | cut -d' ' -f4)
        dir2=$(echo $line2 | cut -d' ' -f1)
        label2=$(echo $line2 | cut -d' ' -f2)
        species2=$(echo $line2 | cut -d' ' -f3)
        sample2=$(echo $line2 | cut -d' ' -f6)
        readlen2=$(echo $line2 | cut -d' ' -f8)
        fraglen2=$(echo $line2 | cut -d' ' -f9)
        fig2=$(echo $line2 | cut -d' ' -f14)
        barcode2=${label2}_${sample2}

        if [ "$barcode1" != "$barcode2" ] && ( [ "$fig1" == "$fig2" ] && [ "$species1" == "$species2" ] && ( $overwrite || ( [ ! -f $outdir/${barcode1}_${barcode2}_1.bed ] && [ ! -f $outdir/${barcode2}_${barcode1}_1.bed ] ) ) ) ; then

            if [ "$species1" == "hg38" ]; then
                GTF=$HG38_GTF
            else
                GTF=$MM10_GTF
            fi

            #bash plot_fig1a.sh $dir1 $sample1 $readlen1 $fraglen1 'x' false $GTF $THRESH $label1
            #bash plot_fig1a.sh $dir2 $sample2 $readlen2 $fraglen2 'x' false $GTF $THRESH $label2

            comm <(sort $dir1/kallisto_results/${sample1}_transcripts_over${THRESH}.gtf) <(sort $dir2/kallisto_results/${sample2}_transcripts_over${THRESH}.gtf) -12 > $outdir/${barcode1}_${barcode2}.gtf
            ngenes=$(cut -d'"' -f 2 $outdir/${barcode1}_${barcode2}.gtf | sort | uniq | wc -l)

            if [ "$PEAKCALLER" == "metdiff" ]; then
                bedtools intersect -a $dir1/${sample1}_MeTDiff_output/MeTDiff_output/peak.sorted.bed -b $outdir/${barcode1}_${barcode2}.gtf -wa | bedtools merge -i - > $outdir/${barcode1}_${barcode2}_1.bed
                bedtools intersect -a $dir2/${sample2}_MeTDiff_output/MeTDiff_output/peak.sorted.bed -b $outdir/${barcode1}_${barcode2}.gtf -wa | bedtools merge -i - > $outdir/${barcode1}_${barcode2}_2.bed
            else
                bedtools intersect -a $dir1/macs2_results/${sample1}_min${min_reps}.bed -b $outdir/${barcode1}_${barcode2}.gtf -wa | bedtools merge -i - > $outdir/${barcode1}_${barcode2}_1.bed
                bedtools intersect -a $dir2/macs2_results/${sample2}_min${min_reps}.bed -b $outdir/${barcode1}_${barcode2}.gtf -wa | bedtools merge -i - > $outdir/${barcode1}_${barcode2}_2.bed
            fi
            peaks1=$(wc -l $outdir/${barcode1}_${barcode2}_1.bed | cut -d' ' -f1 )
            peaks2=$(wc -l $outdir/${barcode1}_${barcode2}_2.bed | cut -d' ' -f1 )
            overlap1=$(bedtools intersect -a $outdir/${barcode1}_${barcode2}_1.bed -b $outdir/${barcode1}_${barcode2}_2.bed -wa | bedtools merge -i - | wc -l) #or use sort | uniq since should all be wa entries
            overlap2=$(bedtools intersect -a $outdir/${barcode1}_${barcode2}_2.bed -b $outdir/${barcode1}_${barcode2}_1.bed -wa | bedtools merge -i - | wc -l)
            echo ${cellline1} ${cellline2} ${label1} ${label2} $ngenes ${min_reps} $overlap1 $peaks1 >> $summary
            echo ${cellline2} ${cellline1} ${label2} ${label1} $ngenes ${min_reps} $overlap2 $peaks2 >> $summary
        fi
    done < $INFI
  fi
done < $INFI
