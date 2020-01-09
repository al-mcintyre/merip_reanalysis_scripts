## Scripts for paper on MeRIP-seq reproducibility and detection of changes in m6A

BioRxiv preprint [here](https://www.biorxiv.org/content/10.1101/657130v1). 

This repository is for scripts to generate the figures. Bam files are too large to host here, but it is possible to download the fastqs from GEO/SRA and rerun alignment using STAR.
Alternatively, run_full_analysis.sh provides a wrapper script to generate plots from saved results (requires all R packages and gtf files for mouse and human). 

Please see R package [DEQ](https://www.github.com/al-mcintyre/deq) to run DESeq2, edgeR, and QNB from input and IP bam files as described in our paper.
