# lrsim: ultrafast and authentic long-read simulator

## Introduction
lrsim is a very authentic long-read simulator, it can simulate error-prone long-read sequencing data with real read length distribution (Currently only provides ONT Ultra Long dist file, but you can generate your own dist file).
Furthermore, lrsim is ultrafast, can generate 30x whole-genome long reads within minutes utilizing 32 CPU cores (comparing to days that badread demands, even with the multi-processing optimized VISOR version).

## Compiling and installation
```shell script
make
make install
```

## Usage
```shell script
#unzip the ldist file
gunzip distdata/HG002_ONT_UL.ldist.gz
lrsim -t 32 -f distdata/HG002_ONT_UL.ldist -b 30G haplotype0.fa haplotype1.fa others.fa > reads.fq
```

## Get your own distribution file
You can generate your own distribution file with our script:
```shell script
samtools stats -@7 real_seq_bam.bam | grep ^RL | cut -f 2- | python3 extractDistData.py > real_seq_bam.ldist
```