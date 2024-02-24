#!/bin/bash
source activate amromics
if [ ! -f GCF_000240185.1_ASM24018v2_genomic.fna.gz ]; then
    wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/240/185/GCF_000240185.1_ASM24018v2/GCF_000240185.1_ASM24018v2_genomic.fna.gz
fi
if [ ! -f GCF_000364385.2_ASM36438v2_genomic.fna.gz ]; then
    wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/364/385/GCF_000364385.2_ASM36438v2/GCF_000364385.2_ASM36438v2_genomic.fna.gz
fi

rm -f sra*txt
awk -F"," '/Illumina/{print $1}' acc_list.csv| parallel -j 16 'fasterq-dump --progress --split-3 {} && gzip {}_1.fastq && gzip {}_2.fastq'
awk -F"," '/Pacbio|ONT/{print $1}' acc_list.csv| parallel -j 16 'fasterq-dump --progress --split-3 {} && gzip {}.fastq'
cp SRR24335776.fastq.gz SRR24335776_1.fastq.gz
cp SRR24335778.fastq.gz SRR24335778_1.fastq.gz
