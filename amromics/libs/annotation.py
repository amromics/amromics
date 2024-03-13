import os
import logging
import re
import multiprocessing
from amromics.libs.bioseq import read_sequence_file
from Bio import SeqIO
from Bio.Seq import Seq
import shutil
from amromics.utils.command import run_command
from amromics.utils.utils import get_open_func
from BCBio import GFF
logger = logging.getLogger(__name__)
NUM_CORES_DEFAULT = multiprocessing.cpu_count()

def annotate_prokka(sample, base_dir='.', timing_log=None, threads=4):
    """
        Run annotation process using prokka
        :param read_data: result holder
        :param base_dir: working directory
        :return: path to output file in result holder
    """
    sample_id = sample['id']

    path_out = os.path.join(base_dir, sample_id+'_prokka' )
    if not os.path.exists(path_out):
        os.makedirs(path_out)

    annotation_gff = os.path.join(path_out,  sample_id +'.gff')
    annotation_faa = os.path.join(path_out,  sample_id +'.faa')
    annotation_ffn = os.path.join(path_out,  sample_id +'.ffn')
    annotation_gbk = os.path.join(path_out,  sample_id +'.gbk')
    cmd = 'prokka --force --cpus {threads} --addgenes --mincontiglen 200'.format(threads=threads)
    cmd += ' --prefix {sample_id} --locus {sample_id} --outdir {path} '.format(sample_id=sample_id, path=path_out)
    if sample['genus']:
        cmd += ' --genus ' + sample['genus']

    if sample['species']:
        species = sample['species'].replace(' ','_')
        cmd += ' --species ' + species

    if sample['strain']:
        cmd += ' --strain ' + sample['strain']
    if sample['gram']:
        cmd += ' --gram ' + sample['gram']
    cmd += ' ' + sample['assembly']

    cmd = "bash -c '{}'".format(cmd)
    ret = run_command(cmd, timing_log)
    if ret != 0:
        raise Exception('Command {} returns non-zero ()!'.format(cmd, ret))

    shutil.move(annotation_gff, sample['annotation_gff'])
    shutil.move(annotation_faa, sample['annotation_faa'])
    shutil.move(annotation_ffn, sample['annotation_ffn'])
    shutil.move(annotation_gbk, sample['annotation_gbk'])
    shutil.rmtree(path_out)

    #return annotation_gff,annotation_faa,annotation_ffn,annotation_fna,annotation_gbk


def parseGFF(sample, gff_file_in):

    sample_id = sample['id']
    annotation_gff, annotation_faa, annotation_ffn,annotation_gbk = sample['annotation_gff'], sample['annotation_faa'], sample['annotation_ffn'], sample['annotation_gbk']
    assemply_file = sample['assembly']

    found_fasta = False
    open_func = get_open_func(gff_file_in)
    #fasta_file=open(path_out+'/'+str(sample_id)+'.fna','w')

    bed_records = []
    with open_func(gff_file_in,'rt') as in_fh, open(annotation_gff, 'w') as gff_re, open(assemply_file,'w') as fna_fh:
        last_cds=""
        id_number=0
        gene_index = 0
        seq_id = None
        suffix = 1
        for line in in_fh:
            if found_fasta == True:
                gff_re.write(line.replace('/','_')) #TODO: why this replacement???
                fna_fh.write(line)
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
    #fasta_file.close()
    seqs = {}
    #count_contigs=0
    for seq in read_sequence_file(sample['assembly']):
        seqs[seq.name] = seq
        #    count_contigs=count_contigs+1
        #print(count_contigs)
        #print(len(seqs.keys()))
    gene_seqs = []
    protein_seqs = []
    for bed_record in bed_records:
        (seq_id, start, end, gene_id, strand,gene_product) = bed_record
        if seq_id not in seqs.keys():
            print(seq_id+" not in fna")
            continue
        seq = seqs[seq_id]
        subseq = Seq(str(seq.sequence)[start:end])
        gene_seq = SeqIO.SeqRecord(seq=subseq,id= gene_id,description=gene_product)

        if strand=='+' :
            protein_seq = gene_seq.translate(table=11, stop_symbol='')
        else:
            protein_seq = gene_seq.reverse_complement().translate(table=11, stop_symbol='')
            gene_seq.seq=subseq.reverse_complement()
        gene_seqs.append(gene_seq)
        protein_seq.id = gene_id
        protein_seq.description = gene_product
        protein_seqs.append(protein_seq)

    with open(annotation_ffn, 'w') as ffn_fh:
        SeqIO.write(gene_seqs, ffn_fh, 'fasta')

    with open(annotation_faa, 'w') as faa_fh:
        SeqIO.write(protein_seqs, faa_fh, 'fasta')
    #convert gff to gbk
    convertGFF2Genbank(annotation_gbk,annotation_gff,assemply_file)
    # for file_name in glob.glob(os.path.join(path_out, '*')):
    #     ext = file_name[-3:]
    #     if ext in ['fna']: # fna?
    #         run_command('gzip {}'.format(file_name))
    #return annotation_gff,annotation_faa,annotation_ffn,annotation_fna
def convertGFF2Genbank(genbank_file, gff_file, fasta_file=None, molecule_type="DNA"):

    if fasta_file:
        fasta_input = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))
    else:
        fasta_input = {}
    gff_iter = GFF.parse(gff_file, fasta_input)
    SeqIO.write(_check_gff(_fix_ncbi_id(gff_iter), molecule_type), genbank_file, "genbank")
def _fix_ncbi_id(fasta_iter):
    """GenBank identifiers can only be 16 characters; try to shorten NCBI.
    """
    for rec in fasta_iter:
        if len(rec.name) > 16 and rec.name.find("|") > 0:
            new_id = [x for x in rec.name.split("|") if x][-1]
            print("Warning: shortening NCBI name %s to %s" % (rec.id, new_id))
            rec.id = new_id
            rec.name = new_id
        yield rec


def _check_gff(gff_iterator, molecule_type):
    """Check GFF files before feeding to SeqIO to be sure they have sequences.
    """
    for rec in gff_iterator:
        if "molecule_type" not in rec.annotations:
            rec.annotations["molecule_type"] = molecule_type
        yield _flatten_features(rec)


def _flatten_features(rec):
    """Make sub_features in an input rec flat for output.

    GenBank does not handle nested features, so we want to make
    everything top level.
    """
    out = []
    for f in rec.features:
        cur = [f]
        while len(cur) > 0:
            nextf = []
            for curf in cur:
                out.append(curf)
                if len(curf.sub_features) > 0:
                    nextf.extend(curf.sub_features)
            cur = nextf
    rec.features = out
    return rec
