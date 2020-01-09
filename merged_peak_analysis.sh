#cat tirumuru_hiv_2016/macs2_results/jurkat_[hm]*.narrowPeak.peaks.bed lichinchi_hiv_2016/macs2_results/*.narrowPeak.peaks.bed | sort -k1,1 -k2,2n | bedtools merge -i - > merged_hiv_peaks.bed

cat tan_kshv_2018/macs2_results/human_iSLK_latent*.narrowPeak.peaks.bed tan_kshv_2018/macs2_results/human_iSLK_kshv_48h*.narrowPeak.peaks.bed hesser_kshv_2018/macs2_results/mock*.narrowPeak.peaks.bed hesser_kshv_2018/macs2_results/kshv*.narrowPeak.peaks.bed | sort -k1,1 -k2,2n | bedtools merge -i - > merged_kshv_peaks.bed

cat winkler_hcmv_2018/macs2_results/hcmv_6h*.narrowPeak.peaks.bed winkler_hcmv_2018/macs2_results/hcmv_72h*.narrowPeak.peaks.bed rubio_hcmv_2018/macs2_results/hcmv*.narrowPeak.peaks.bed rubio_hcmv_2018/macs2_results/mock*.narrowPeak.peaks.bed | sort -k1,1 -k2,2n | bedtools merge -i - > merged_dsDNA_peaks.bed


