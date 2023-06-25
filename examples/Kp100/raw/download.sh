#!/bin/bash
source activate amromics
if [ ! -f GCF_000240185.1_ASM24018v2_genomic.fna.gz ]; then
    wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/240/185/GCF_000240185.1_ASM24018v2/GCF_000240185.1_ASM24018v2_genomic.fna.gz
fi
if [ ! -f GCF_000364385.2_ASM36438v2_genomic.fna.gz ]; then
    wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/364/385/GCF_000364385.2_ASM36438v2/GCF_000364385.2_ASM36438v2_genomic.fna.gz
fi

### Faster downloading:
# Get FTP link from EBI
rm -f sra*txt
awk -F"," '/Illumina/{print $1}' acc_list.csv| parallel -j 16 'curl -X GET "https://www.ebi.ac.uk/ena/portal/api/filereport?accession={}&fields=fastq_ftp&result=read_run" 2>/dev/null | sed -n "2 p" >> sra_short.ftp.txt'
awk -F"," '/Pacbio|ONT/{print $1}' acc_list.csv| parallel -j 16 'curl -X GET "https://www.ebi.ac.uk/ena/portal/api/filereport?accession={}&fields=fastq_ftp&result=read_run" 2>/dev/null | sed -n "2 p" >> sra_long.ftp.txt'

# FTP download is slow too, convert to fasq link to download by ascp
sed 's/;/\t/g' sra_short.ftp.txt | sed 's/ftp.sra.ebi.ac.uk\//era-fasp@fasp.sra.ebi.ac.uk:/g' > sra_short.ascp.txt
totN=$(cat sra_short.ascp.txt|wc -l)
curN=0
KEYPATH="${CONDA_PREFIX}"/opt/aspera/connect/etc/asperaweb_id_dsa.openssh
while read s r1 r2
do
	curN=$((curN + 1))
	echo "Downloading $s...($curN/$totN) short-reads data:"
	[ -f ${s}_1.fastq.gz ] && echo "${s}_1.fastq.gz exist...skip" || ascp -QT -l 300m -P33001 -i "${KEYPATH}" "${r1}" .	
	[ -f ${s}_2.fastq.gz ] && echo "${s}_2.fastq.gz exist...skip" || ascp -QT -l 300m -P33001 -i "${KEYPATH}" "${r2}" .	
	echo "Finish download $s."
done < sra_short.ascp.txt


sed 's/;/\t/g' sra_long.ftp.txt | sed 's/ftp.sra.ebi.ac.uk\//era-fasp@fasp.sra.ebi.ac.uk:/g' > sra_long.ascp.txt
totN=$(cat sra_long.ascp.txt|wc -l)
curN=0
while read s r
do
	curN=$((curN + 1))
	echo "Downloading $s...($curN/$totN) long-reads data:"
	compgen -G ${s}*fastq.gz > /dev/null && echo "${s} exist...skip" || ascp -QT -l 300m -P33001 -i "${KEYPATH}" "${r}" .	
	echo "Finish download $s."
done < sra_long.ascp.txt

