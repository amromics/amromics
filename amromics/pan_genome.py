import os
import shutil
import re
import logging
from Bio import SeqIO
from Bio.Seq import Seq
from amromics.utils import run_command

logger = logging.getLogger(__name__)

def get_input(handle, report):
    """
    Get and extract ggf files and assembly files to temp folder

    Parameters
    -------
    -------
    """
    handle['samples'] = []
    for sample in report['samples']:
        sample_id = sample['id']
        gffgz_file = os.path.join(sample['annotation'], sample_id + '.gff.gz')
        gff_file = os.path.join(handle['temp_folder'], sample_id + '.gff')
        if run_command('gunzip -c {} > {}'.format(gffgz_file, gff_file)) != 0:
            raise Exception('Cannot get {}'.format(gffgz_file))

        assembly_gz_file = report['assembly']
        assembly_file = os.path.join(handle['temp_folder'], sample_id + '.fna')
        if run_command('gunzip -c {} > {}'.format(assembly_gz_file, assembly_file)) != 0:
            raise Exception('Cannot get {}'.format(assembly_gz_file))
        
        sample_handle['id'] = sample_id
        sample_handle['gff_file'] = gff_file
        sample_handle['assembly_file'] = assembly_file
        handle['samples'].append(sample_handle)

    return handle


def create_bed_file_from_gff(ggf_file, bed_file):
    """
    Take in gff file and return bed file. Filter out small genes(<18 base)

    Parameters
    -------
    -------
    """
    in_handle = open(ggf_file,'r')
    out_handle = open(bed_file, 'w')
    found_fasta = 0
    for line in in_handle:
        if found_fasta == 1:
            continue
        if re.match(r"^##FASTA", line) != None:
            found_fasta = 1
            continue
        if re.match(r"^##", line) != None:
            continue
        line = line.rstrip('\n')
        cells = line.split('\t')
        if cells[2] != 'CDS':
            continue
            
        #filter out small genes (<18 base)
        start = int(cells[3])
        end = int(cells[4])
        if end - start < 18:
            continue

        gene_id = re.findall(r"ID=(.+?);",cells[8])
        gene_id = gene_id[0]
        trand = cells[6]
        seq_id = cells[0]
        out_handle.write(seq_id+'\t'+str(start -1)+'\t'+ str(end)+'\t'+gene_id+'\t'+trand+'\n')
    in_handle.close()
    out_handle.close()


def extract_proteins(sample, out_dir, timing_log=None):
    """
    Take in GFF file and fasta file and create protein sequences in FASTA format

    Parameters
    -------
    -------
    """
    fna_file = sample['assembly_file']
    ggf_file = sample['gff_file']
    sample_id = sample['id']
    # extract nucleotide region
    extracted_fna_file = os.path.join(out_dir, sample_id +'.extracted.fna')
    bed_file = os.path.join(out_dir, sample_id +'.bed')
    create_bed_file_from_gff(ggf_file, bed_file)
    cmd = f"bedtools getfasta -s -fi {fna_file} -bed {bed_file} -fo {extracted_fna_file} -name > /dev/null 2>&1"
    ret = run_command(cmd, timing_log)
    if ret != 0:
        raise Exception('Error running bedtools')
    # translate nucleotide to protein
    faa_file = os.path.join(out_dir, sample_id +'.faa')
    faa_hd = open(faa_file, 'w')
    for seq_record in SeqIO.parse(extracted_fna_file, "fasta"):
        headers = seq_record.id.split(':')
        seq_record.id = headers[0]
        seq_record.name = ''
        seq_record.description = ''
        seq_record.seq = seq_record.seq.translate(table=11)
        SeqIO.write(seq_record, faa_hd, 'fasta')
    faa_hd.close()

    sample['bed'] = bed_file
    sample['extracted_fna_file'] = extracted_fna_file
    sample['faa_file'] = extracted_fna_file
    return sample
    

#def combine_proteins():


#def run_cdhit_iterative():


#def all_against_all_blast():


#def cluster_mcl():


#def reinflate_clusters():


#def split_paralogs():


#def label_group():


#def annotate_group():


def run_pan_genome_analysis(report, collection_dir='.', threads=8, overwrite=False, timing_log=None):
    # Check if any sample has been updated
    for sample in report['samples']:
        overwrite = overwrite or sample['updated']

    pan_genome_folder = os.path.join(collection_dir, 'pan_genome')
    temp_folder = os.path.join(collection_dir, 'temp')
    pan_genome_output = os.path.join(pan_genome_folder,'')
    report['pan_genome'] = pan_genome_folder
    
    # Check if pan-genome has run
    if os.path.isfile(pan_genome_output) and (not overwrite):
        logger.info('pan-genome has run and the input has not changed, skip')
        return report

    if not os.path.exists(pan_genome_folder):
        os.makedirs(pan_genome_folder)
    if not os.path.exists(temp_folder):
        os.makedirs(temp_folder)
    
    handle = {}
    handle['main_folder'] = pan_genome_folder
    handle['temp_folder'] = temp_folder

    handle = get_input(handle=handle,report=report)

    for sample in handle['samples']:
        sample_temp_dir = os.path.join(temp_folder, sample['id'])
        sample = extract_proteins(sample, out_dir = handle['temp_folder'])



sample = {
    'assembly_file' : '/home/ted/ubuntu/AMR/amromics-vis/dev/BMP071.fna',
    'gff_file' : '/home/ted/ubuntu/AMR/amromics-vis/dev/BMP071.gff',
    'id':'BMP071'
}

extract_proteins(sample, out_dir = '/home/ted/ubuntu/AMR/amromics-vis/dev')



## TO-DO ##
# consistence between handle and report
# check protein translate
# filter protein sequence