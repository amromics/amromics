#!/bin/bash
source activate amromics
if [ ! -f GCF_000240185.1_ASM24018v2_genomic.fna.gz ]; then
    wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/240/185/GCF_000240185.1_ASM24018v2/GCF_000240185.1_ASM24018v2_genomic.fna.gz
fi

if [ ! -f GCF_000364385.2_ASM36438v2_genomic.fna.gz ]; then
    wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/364/385/GCF_000364385.2_ASM36438v2/GCF_000364385.2_ASM36438v2_genomic.fna.gz
fi

### Download by fasterq-dump too slow
#for ac in SRR8607448 SRR8607464 SRR8607467 SRR8607459 SRR8607460 SRR8607456 SRR8607451 SRR8607450 SRR8607452 SRR8607471 SRR8607469 SRR8607462 SRR8607457 SRR8607455 SRR8607461 SRR8607458 SRR8607470 SRR8607449 SRR8607465 SRR8607454 SRR8607468 SRR8607463 SRR8607466 SRR8607453;do
#    if [ ! -f ${ac}_2.fastq.gz ];then
#        echo " Downloading Accession ${ac} ..."
#        fasterq-dump --progress --split-3 ${ac} && gzip ${ac}_1.fastq && gzip ${ac}_2.fastq
#    else
#        echo " Accession ${ac} has been downloaded!"
#    fi
#done

### Faster downloading:
# Get FTP link from EBI
:>sra.ftp.txt
parallel -j 32 'curl -X GET "https://www.ebi.ac.uk/ena/portal/api/filereport?accession={}&fields=fastq_ftp&result=read_run" 2>/dev/null | tail -n1 >> sra.ftp.txt' ::: SRR8607448 SRR8607464 SRR8607467 SRR8607459 SRR8607460 SRR8607456 SRR8607451 SRR8607450 SRR8607452 SRR8607471 SRR8607469 SRR8607462 SRR8607457 SRR8607455 SRR8607461 SRR8607458 SRR8607470 SRR8607449 SRR8607465 SRR8607454 SRR8607468 SRR8607463 SRR8607466 SRR8607453

# FTP download is slow too, convert to fasq link to download by ascp
sed 's/;/\t/g' sra.ftp.txt | sed 's/ftp.sra.ebi.ac.uk\//era-fasp@fasp.sra.ebi.ac.uk:/g' > sra.ascp.txt
totN=$(cat sra.ascp.txt|wc -l)
curN=0
KEYPATH="${CONDA_PREFIX}"/opt/aspera/connect/etc/asperaweb_id_dsa.openssh
while read s r1 r2
do
	curN=$((curN + 1))
	echo "Downloading $s...($curN/$totN):"
	[ -f ${s}_1.fastq.gz ] && echo "${s}_1.fastq.gz exist...skip" || ascp -QT -l 300m -P33001 -i "${KEYPATH}" "${r1}" .	
	[ -f ${s}_2.fastq.gz ] && echo "${s}_2.fastq.gz exist...skip" || ascp -QT -l 300m -P33001 -i "${KEYPATH}" "${r2}" .	
	echo "Finish download $s."
done < sra.ascp.txt

