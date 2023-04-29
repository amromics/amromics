import os
import shutil
import csv
import logging
import gzip
import glob
import re
import multiprocessing
from amromics.libs.bioseq import read_sequence_file
from Bio import SeqIO
from Bio.Seq import Seq
from amromics.utils.command import run_command
from amromics.utils.utils import get_open_func
logger = logging.getLogger(__name__)
NUM_CORES_DEFAULT = multiprocessing.cpu_count()
def annotate_prokka(prefix_name,assembly,genus=None,species=None, strain=None,gram=None, base_dir='.',  overwrite=False,timing_log=None, threads=0):
    """
        Run annotation process using prokka
        :param read_data: result holder
        :param base_dir: working directory
        :return: path to output file in result holder
    """
    if threads == 0:
        threads = NUM_CORES_DEFAULT

    path_out = os.path.join(base_dir, prefix_name+'_prokka' )
    if not os.path.exists(path_out):
        os.makedirs(path_out)

    annotation_gbk=  os.path.join(path_out, prefix_name + '.gbk.gz')
    annotation_gff= path_out+'/'+str(prefix_name)+'.gff.gz'
    annotation_faa = path_out+'/'+str(prefix_name)+'.faa.gz'
    annotation_ffn = path_out+'/'+str(prefix_name)+'.ffn.gz'
    annotation_fna = path_out+'/'+str(prefix_name)+'.fna.gz'
    if os.path.isfile(annotation_gff) and os.path.isfile(annotation_gbk) and (not overwrite):
        # Dont run again if gff/gbk file exists
        logger.info('GFF and GBK files found, skip annotating')
        return annotation_gff,annotation_faa,annotation_ffn,annotation_fna,annotation_gbk
    gunzip_fasta=assembly
    if assembly.endswith('.gz'):
        gunzip_fasta = os.path.join(path_out, prefix_name + '.fin')
        cmd = 'gunzip -c {} > {}'.format(assembly, gunzip_fasta)
        run_command(cmd)
    cmd = 'prokka --force --cpus {threads} --addgenes --mincontiglen 200'.format(threads=threads)
    cmd += ' --prefix {sample_id} --locus {sample_id} --outdir {path} '.format(sample_id=prefix_name, path=path_out)
    if not genus ==None and genus:
        cmd += ' --genus ' +genus
    if not species  ==None and species:
        species = species.replace(' ','_')
        cmd += ' --species ' + species
    if not strain  ==None and strain:
        cmd += ' --strain ' + strain
    if not gram  ==None and gram:
        cmd += ' --gram ' + gram
    cmd += ' ' + gunzip_fasta
    cmd = "bash -c '{}'".format(cmd)
    ret = run_command(cmd, timing_log)
    if ret != 0:
        raise Exception('Command {} returns non-zero ()!'.format(cmd, ret))

    for file_name in glob.glob(os.path.join(path_out, '*')):
        ext = file_name[-3:]
        if ext in ['gff', 'gbk', 'ffn','faa','fna']: # fna?
            run_command('gzip {}'.format(file_name))
        else:
            os.remove(file_name)

    return annotation_gff,annotation_faa,annotation_ffn,annotation_fna,annotation_gbk
def parseGFF(sample_id,gff_file_in,base_dir):

    path_out = os.path.join(base_dir, sample_id+'_prokka' )

    if not os.path.exists(path_out):
        os.makedirs(path_out)
    annotation_gff= path_out+'/'+str(sample_id)+'.gff.gz'
    annotation_ffn = path_out+'/'+str(sample_id)+'.ffn.gz'
    annotation_fna=path_out+'/'+str(sample_id)+'.fna'
    found_fasta = False
    open_func = get_open_func(gff_file_in)
    fasta_file=open(annotation_fna,'wt')
    with open_func(gff_file_in,'rt') as in_fh,gzip.open(annotation_gff, 'wt') as gff_re:
        last_cds=""
        id_number=0


        bed_records = []
        gene_index = 0
        seq_id = None
        suffix = 1
        for line in in_fh:

            if found_fasta == True:
                gff_re.write(line.replace('/','_'))
                fasta_file.write(line)
                continue
            if re.match(r"^##FASTA", line) != None:
                found_fasta = True
                gff_re.write(line)
                continue
            if re.match(r"^#", line) != None:
                gff_re.write(line)
                continue
            line = line.rstrip('\n')
            cells = line.split('\t')
            if cells[2] != 'CDS':
                continue


            ##print(line)
            #print(gff_file_in)
            strand = cells[6]
            start = int(cells[3])
            end = int(cells[4])
            length = end - start + 1
            if length % 3 != 0:
                continue
            cells[0] = cells[0].replace('-','_')
            if seq_id != cells[0]:
                seq_id = cells[0]
                gene_index = 0
            tags=cells[8].split(';')
            gene_id=None
            gene_name = ''
            gene_product = ''
            for tag in tags:
                gene = re.match(r"^gene=(.+)", tag)
                if gene != None:
                    gene_name = gene.group(1)
                    gene_name = re.sub(r'\W', '_', gene_name)
                    continue
                product = re.match(r"^product=(.+)", tag)
                if product != None:
                    gene_product = product.group(1)
            gene_id = sample_id + '_' +f'{suffix:05d}'
            #print(gene_id)
            suffix += 1
            newdes="ID="+gene_id+";"
            for i in range(2,len(tags)):
                newdes=newdes+tags[i]+";"
            newline=""
            for i in range(len(cells)-1):
                newline=newline+cells[i]+"\t"
            newline=newline+newdes+"\n"
            gff_re.write(newline)
            bed_records.append((cells[0].strip(), start - 1, end, gene_id, strand,gene_product))
            gene_index += 1
        fasta_file.close()
        seqs = {}
        count_contigs=0
        for seq in read_sequence_file(annotation_fna):
            seqs[seq.name] = seq
            count_contigs=count_contigs+1
        print(count_contigs)
        print(len(seqs.keys()))
        gene_seqs = []
        for bed_record in bed_records:
            (seq_id, start, end, gene_id, strand,gene_product) = bed_record
            if seq_id not in seqs.keys():
                print(seq_id+" not in fna")
                continue
            seq = seqs[seq_id]
            subseq = Seq(str(seq.sequence)[start:end])
            gene_seq = SeqIO.SeqRecord(seq=subseq,id= gene_id,description=gene_product)
            gene_seqs.append(gene_seq)
        with gzip.open(annotation_ffn, 'wt') as ffn_fh:
            SeqIO.write(gene_seqs, ffn_fh, 'fasta')
    return annotation_gff,annotation_ffn
