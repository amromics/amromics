import os
import shutil
import re
import json
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
        gff_file = os.path.join(handle['temp_dir'], sample_id + '.gff')
        if run_command('gunzip -c {} > {}'.format(gffgz_file, gff_file)) != 0:
            raise Exception('Cannot get {}'.format(gffgz_file))

        assembly_gz_file = report['assembly']
        assembly_file = os.path.join(handle['temp_dir'], sample_id + '.fna')
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


def extract_proteins(handle, timing_log=None):
    """
    Take in GFF file and fasta file and create protein sequences in FASTA format

    Parameters
    -------
    -------
    """
    temp_dir = handle['temp_dir']
    for sample in handle['samples']:
        fna_file = sample['assembly_file']
        ggf_file = sample['gff_file']
        sample_id = sample['id']
        # extract nucleotide region
        extracted_fna_file = os.path.join(temp_dir, sample_id +'.extracted.fna')
        bed_file = os.path.join(temp_dir, sample_id +'.bed')
        create_bed_file_from_gff(ggf_file, bed_file)
        cmd = f"bedtools getfasta -s -fi {fna_file} -bed {bed_file} -fo {extracted_fna_file} -name > /dev/null 2>&1"
        ret = run_command(cmd, timing_log)
        if ret != 0:
            raise Exception('Error running bedtools')
        # translate nucleotide to protein
        faa_file = os.path.join(temp_dir, sample_id +'.faa')
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
        sample['faa_file'] = faa_file
    return handle
    

def combine_proteins(handle, timing_log=None):
    """
    Take in multiple FASTA sequences containing proteomes and concat them together 
    and output a FASTA file

    Parameters
    -------
    -------
    """
    temp_dir = handle['temp_dir']
    combined_faa_file = os.path.join(temp_dir, 'combined.faa')
    faa_file_list =[]
    for sample in handle['samples']:
        faa_file = sample['faa_file']
        faa_file_list.append(faa_file)

    cmd = "cat {} > {}".format(" ".join(faa_file_list),combined_faa_file)
    ret = run_command(cmd, timing_log)
    if ret != 0:
        raise Exception('Error combining protein sequences')

    handle['combined_faa_file'] = combined_faa_file
    return handle


def parse_cluster_file(cd_hit_cluster_file):
    cd_hit_cluster_fh = open(cd_hit_cluster_file, 'r')
    clusters = {}
    for line in cd_hit_cluster_fh:
        if re.match(r"^>", line) != None:
            cluster_name = re.findall(r'^>(.+)$', line)
            cluster_name = cluster_name[0]
            clusters[cluster_name] = {}
            clusters[cluster_name]['gene_names'] = []
        else:
            result = re.findall(r'[\d]+\t[\w]+, >(.+)\.\.\. (.+)$', line)
            if len(result) == 1:
                gene_name = result[0][0]
                identity = result[0][1]
                if identity == '*':
                    clusters[cluster_name]['representative'] = gene_name
                else:
                    clusters[cluster_name]['gene_names'].append(gene_name)

    return clusters


def run_cd_hit_iterative(handle, threads, timing_log):
    temp_dir = handle['temp_dir']
    combined_faa_file = handle['combined_faa_file']

    cd_hit_cluster_output = os.path.join(temp_dir, 'cluster')
    cd_hit_cluster_file = os.path.join(temp_dir, 'cluster.clstr')

    removed_cluster_file = os.path.join(temp_dir, 'removed_cluster')
    removed_cluster_fh = open(removed_cluster_file, 'a')

    percent_match = 1
    greater_than_or_equal = True
    number_of_samples = len(handle['samples'])
    
    while percent_match >= 0.98:
        cmd = f'cd-hit -i {combined_faa_file} -o {cd_hit_cluster_output} -s {percent_match} -c {percent_match} -T {threads} -M 0 -g 1 -d 256 '
        ret = run_command(cmd, timing_log)
        if ret != 0:
            raise Exception('Error running cd-hit')
        
        clusters = parse_cluster_file(cd_hit_cluster_file)

        if percent_match != 1:
            greater_than_or_equal = False
        full_cluster_gene_names =[]
        for cluster_name in clusters.keys():
            representative_gene = clusters[cluster_name]['representative']
            other_genes = clusters[cluster_name]['gene_names']
            if greater_than_or_equal == True and len(other_genes) >= number_of_samples -1:
                full_cluster_gene_names.append(representative_gene)
                full_cluster_gene_names.extend(other_genes)
                removed_cluster_fh.write(representative_gene+'\t'+'\t'.join(other_genes)+'\n')
            if greater_than_or_equal == False and len(other_genes) == number_of_samples -1:
                full_cluster_gene_names.append(representative_gene)
                full_cluster_gene_names.extend(other_genes)
                removed_cluster_fh.write(representative_gene+'\t'+'\t'.join(other_genes)+'\n')
        

        cluster_filtered_faa_file = combined_faa_file + '.filtered'
        cluster_filtered_faa_fh = open(cluster_filtered_faa_file, 'w')
        for seq_record in SeqIO.parse(combined_faa_file, "fasta"):
            if seq_record.id in full_cluster_gene_names:
                SeqIO.write(seq_record, cluster_filtered_faa_fh, 'fasta')
        cluster_filtered_faa_fh.close()

        shutil.move(cluster_filtered_faa_file, combined_faa_file)
        percent_match -=0.005

    removed_cluster_fh.close()
    
    handle['cd_hit_cluster_output'] = cd_hit_cluster_output
    handle['cd_hit_cluster_file'] = cd_hit_cluster_file
    handle['removed_cluster_file'] = removed_cluster_file
    return handle

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
    temp_dir = os.path.join(collection_dir, 'temp')
    pan_genome_output = os.path.join(pan_genome_folder,'')
    report['pan_genome'] = pan_genome_folder
    
    # Check if pan-genome has run
    if os.path.isfile(pan_genome_output) and (not overwrite):
        logger.info('pan-genome has run and the input has not changed, skip')
        return report

    if not os.path.exists(pan_genome_folder):
        os.makedirs(pan_genome_folder)
    if not os.path.exists(temp_dir):
        os.makedirs(temp_dir)
    
    handle = {}
    handle['main_dir'] = pan_genome_folder
    handle['temp_dir'] = temp_dir

    handle = get_input(handle=handle,report=report)
    handle = extract_proteins(handle, timing_log=timing_log)
    handle = combine_proteins(handle, timing_log=timing_log)
    handle = run_cd_hit_iterative(handle, threads=threads, timing_log=timing_log)



if __name__ == '__main__':
    
    handle = {
        'main_dir' : '/home/ted/ubuntu/AMR/amromics-vis/dev',
        'temp_dir' : '/home/ted/ubuntu/AMR/amromics-vis/dev/temp',
        'samples':[
            {
                'assembly_file' : '/home/ted/ubuntu/AMR/amromics-vis/dev/BMP071.fna',
                'gff_file' : '/home/ted/ubuntu/AMR/amromics-vis/dev/BMP071.gff',
                'id':'BMP071'
            },
            {
                'assembly_file' : '/home/ted/ubuntu/AMR/amromics-vis/dev/BMP077.fna',
                'gff_file' : '/home/ted/ubuntu/AMR/amromics-vis/dev/BMP077.gff',
                'id':'BMP077'
            }
        ]
    }
    threads = 4
    timing_log = '/home/ted/ubuntu/AMR/amromics-vis/dev/time.log'
    handle = extract_proteins(handle, timing_log=timing_log)
    handle = combine_proteins(handle, timing_log=timing_log)
    handle = run_cd_hit_iterative(handle, threads=threads, timing_log=timing_log)
    print(json.dumps(handle, indent=4, sort_keys=True))
    



## TO-DO ##
# handle or report ???
# check protein translate (compare with perl script)
    # filtered protein sequence produced by roary did not have * character in the middle of sequence => check unfiltered protein sequence 
    # protein sequence with the same id are different from each other => run roary perl scrip individually
# filter protein sequence
# warming of CD-HIT: Discarding invalid sequence or sequence without identifier and description!