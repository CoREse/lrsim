# lrsim: ultrafast and authentic long-read simulator

## Introduction
lrsim is an authentic long-read simulator, it can simulate error-prone long-read sequencing data with real read length distribution and insertion:deletion ratio (Currently only provides ONT Ultra Long model file, but you can generate your own model file), it also provides a lot of custom parameters for you to tweak with.
Furthermore, lrsim is ultrafast, can generate 30x whole-genome long reads within minutes (usually limited by your IO speed) utilizing 32 CPU cores (comparing to days that badread demands, even with the multi-processing optimized VISOR version).

## Compiling and installation
You should first clone or download the source, and then enter the source directory.

Subsequently:
```shell script
git submodule update --init --recursive
make
#Install to system default location
make install
#Or if you want to install it to location you specify (will install lrsim to /PATH/TO/YOUR/LOCAL/bin):
PREFIX=/PATH/TO/YOUR/LOCAL make install
```

## Usage
```shell script
#unzip the model file
gunzip -c models/HG002_ONT_UL.lrsm.gz > models/HG002_ONT_UL.lrsm
#or simply run
bash unzipmodels.sh
#run simulation
lrsim -t 32 -m models/HG002_ONT_UL.lrsm -b 30G haplotype0.fa haplotype1.fa others.fa > reads.fq
```

## Get your own model file
You can generate your own model file to imitate your real long-read data with our script:
```shell script
samtools stats -@7 real_seq_bam.bam | ./extractModel.py > models/real_seq_bam.lrsm
```