import os
import shutil
import re
import json
import logging
from Bio import SeqIO
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
    """
    Parse cd-hit .clstr file

    Parameters
    -------
    -------
    """
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
    # convert to a simple dictionary
    clusters_new ={}
    for cluster_name in clusters.keys():
        clusters_new[clusters[cluster_name]['representative']] = clusters[cluster_name]['gene_names']

    return clusters_new


def run_cd_hit_iterative(handle, threads, timing_log):
    """
    Run CD-HIT iteratively

    Parameters
    -------
    -------
    """
    temp_dir = handle['temp_dir']
    combined_faa_file = handle['combined_faa_file']

    cd_hit_cluster_represent = os.path.join(temp_dir, 'cluster')
    cd_hit_cluster_file = os.path.join(temp_dir, 'cluster.clstr')

    excluded_cluster_file = os.path.join(temp_dir, 'excluded_cluster')
    if os.path.isfile(excluded_cluster_file):
        os.remove(excluded_cluster_file)
    excluded_cluster_fh = open(excluded_cluster_file, 'a')

    percent_match = 1
    greater_than_or_equal = True
    number_of_samples = len(handle['samples'])
    
    while percent_match >= 0.98:
        cmd = f'cd-hit -i {combined_faa_file} -o {cd_hit_cluster_represent} -s {percent_match} -c {percent_match} -T {threads} -M 0 -g 1 -d 256 > /dev/null'
        ret = run_command(cmd, timing_log)
        if ret != 0:
            raise Exception('Error running cd-hit')
        
        clusters = parse_cluster_file(cd_hit_cluster_file)

        if percent_match != 1:
            greater_than_or_equal = False
        full_cluster_gene_names =[]
        for cluster_represent in clusters.keys():
            other_genes = clusters[cluster_represent]
            if greater_than_or_equal == True and len(other_genes) >= number_of_samples -1:
                full_cluster_gene_names.append(cluster_represent)
                full_cluster_gene_names.extend(other_genes)
                excluded_cluster_fh.write(cluster_represent+'\t'+'\t'.join(other_genes)+'\n')
            if greater_than_or_equal == False and len(other_genes) == number_of_samples -1:
                full_cluster_gene_names.append(cluster_represent)
                full_cluster_gene_names.extend(other_genes)
                excluded_cluster_fh.write(cluster_represent+'\t'+'\t'.join(other_genes)+'\n')

        cluster_filtered_faa_file = combined_faa_file + '.filtered'
        cluster_filtered_faa_fh = open(cluster_filtered_faa_file, 'w')
        for seq_record in SeqIO.parse(combined_faa_file, "fasta"):
            if seq_record.id in full_cluster_gene_names:
                SeqIO.write(seq_record, cluster_filtered_faa_fh, 'fasta')
        cluster_filtered_faa_fh.close()

        shutil.move(cluster_filtered_faa_file, combined_faa_file)
        percent_match -=0.005

    excluded_cluster_fh.close()
    
    cluster_json = os.path.join(temp_dir, 'cluster.json')
    with open(cluster_json, 'w') as cluster_json_hd:
        json.dump(clusters, cluster_json_hd, indent=4, sort_keys=True)


    handle['cd_hit_cluster_represent'] = cd_hit_cluster_represent
    handle['cd_hit_cluster_file'] = cd_hit_cluster_file
    handle['excluded_cluster_file'] = excluded_cluster_file
    handle['cluster_json'] = cluster_json
    return handle


def chunk_fasta_file(fasta_file, out_dir):
    """
    Take a fasta file and chunk it up into smaller files with the length of 200000

    Parameters
    -------
    -------
    """
    if os.path.exists(out_dir):
        shutil.rmtree(out_dir)
        os.makedirs(out_dir)
    else:
        os.makedirs(out_dir)
    
    chunked_file_list = []
    chunk_number = 0
    current_chunk_length = 0
    chunked_file = os.path.join(out_dir, str(chunk_number) + '.seq')
    chunked_fh = open(chunked_file, 'a')
    chunked_file_list.append(chunked_file)
    for seq_record in SeqIO.parse(fasta_file, "fasta"):
        if current_chunk_length > 200000:
            chunked_fh.close()
            chunk_number += 1
            current_chunk_length = 0
            chunked_file = os.path.join(out_dir, str(chunk_number) + '.seq')
            chunked_file_list.append(chunked_file)
            chunked_fh = open(chunked_file, 'a')
            SeqIO.write(seq_record, chunked_fh, 'fasta')
        chunked_file = os.path.join(out_dir, str(chunk_number) + '.seq')
        SeqIO.write(seq_record, chunked_fh, 'fasta')
        current_chunk_length += len(seq_record.seq)
    chunked_fh.close()
    return chunked_file_list


def all_against_all_blast(handle, threads, timing_log):
    """
    Run all against all blast in parallel

    Parameters
    -------
    -------
    """
    out_dir = os.path.join(handle['temp_dir'], 'blast')
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    # make blast database
    fasta_file = handle['cd_hit_cluster_represent']
    blast_db = os.path.join(out_dir, 'output_contigs')
    cmd = f"makeblastdb -in {fasta_file} -dbtype prot -out {blast_db} -logfile /dev/null"
    ret = run_command(cmd, timing_log)
    if ret != 0:
        raise Exception('Error running makeblastdb')
    
    # chunk fasta file
    chunk_dir = os.path.join(out_dir, 'chunk_files')
    chunked_file_list = chunk_fasta_file(fasta_file, chunk_dir)

    # run parallel all-against-all blast
    blast_cmds_file = os.path.join(out_dir,"blast_cmds.txt")
    blast_cmds_fh = open(blast_cmds_file,'w')
    blast_output_file_list = []
    for chunked_file in chunked_file_list:
        blast_output_file = os.path.splitext(chunked_file)[0] + '.out'
        blast_output_file_list.append(blast_output_file)
        cmd = f"blastp -query {chunked_file} -db {blast_db} -evalue 1E-6 -num_threads {threads} -outfmt 6 -max_target_seqs 2000 " + "| awk '{ if ($3 > 98) print $0;}' 2> /dev/null 1> " + blast_output_file
        blast_cmds_fh.write(cmd + '\n')
    blast_cmds_fh.close()  
    cmd = f"parallel -a {blast_cmds_file}"
    ret = run_command(cmd, timing_log)
    if ret != 0:
        raise Exception('Error running parallel all-against-all blast')

    # combining blast results
    blast_result_file = os.path.join(handle['temp_dir'], 'blast_results')
    if os.path.isfile(blast_result_file):
        os.remove(blast_result_file)
    for blast_output_file in blast_output_file_list:
        os.system(f'cat {blast_output_file} >> {blast_result_file}')

    handle['blast_result_file'] = blast_result_file
    return handle


def cluster_with_mcl(handle, threads, timing_log):
    """
    Take blast results and outputs clustered results

    Parameters
    -------
    -------
    """
    out_dir = handle['temp_dir']
    blast_results = handle['blast_result_file']
    output_mcl_file = os.path.join(out_dir, 'uninflated_mcl_groups')
    cmd = f"mcxdeblast -m9 --score r --line-mode=abc {blast_results} 2> /dev/null | mcl - --abc -I 1.5 -o {output_mcl_file} > /dev/null 2>&1"
    ret = run_command(cmd, timing_log)
    if ret != 0:
        raise Exception('Error running mcl')
    handle['uninflated_mcl_groups'] = output_mcl_file
    return handle


def reinflate_clusters(handle):
    """
    Take the clusters file from cd-hit and use it to inflate the output of MCL

    Parameters
    -------
    -------
    """
    with open(handle['cluster_json'],'r') as fh:
        clusters = json.load(fh)
    mcl_file = handle['uninflated_mcl_groups']
    mcl_fh = open(mcl_file, 'r')
    output_file = os.path.join(handle['temp_dir'], 'inflated_unsplit_mcl_groups')
    output_fh = open(output_file, 'w')
    # Inflate genes from cdhit which were sent to mcl
    for line in mcl_fh:
        inflated_genes = []
        line = line.rstrip('\n')
        genes = line.split('\t')
        for gene in genes:
            inflated_genes.append(gene)
            if gene in clusters:
                inflated_genes.extend(clusters[gene])
                del clusters[gene]
        output_fh.write('\t'.join(inflated_genes)+'\n')
    mcl_fh.close()

    # Inflate any clusters that were in the clusters file but not sent to mcl
    for gene in clusters.keys():
        output_fh.write(gene + '\t' + '\t'.join(cluster[gene])+ '\n')

    # Add clusters which were excluded
    excluded_cluster_file = handle['excluded_cluster_file']
    excluded_cluster_fh = open(excluded_cluster_file, 'r')
    for line in excluded_cluster_fh:
        output_fh.write(line)
    excluded_cluster_fh.close()
    output_fh.close()

    handle['inflated_unsplit_mcl_groups'] = output_file
    return handle

#def split_paralogs():


def label_group(handle):
    """
    Add labels to the groups

    Parameters
    -------
    -------
    """
    input_file = handle['inflated_unsplit_mcl_groups']
    ouput_file = os.path.join(handle['temp_dir'], 'labeled_mcl_groups')
    input_fh = open(input_file, 'r')
    output_fh = open(ouput_file, 'w')
    counter = 1
    for line in input_fh:
        output_fh.write('groups_' + str(counter) + '\t' + line )
        counter += 1
    input_fh.close()
    output_fh.close()
    handle['labeled_mcl_groups'] = ouput_file
    return handle


def annotate_group():
    """
    Take in a group file and associated GFF files for the isolates 
    and update the group name to the gene name

    Parameters
    -------
    -------
    """


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
    handle = all_against_all_blast(handle, threads=threads, timing_log=timing_log)
    handle = cluster_with_mcl(handle, threads=threads, timing_log=timing_log)
    handle = reinflate_clusters(handle)
    handle = label_group(handle)


if __name__ == '__main__':
    handle = {
        "blast_result_file": "/home/ted/ubuntu/AMR/amromics-vis/dev/temp/blast_results",
        "cd_hit_cluster_file": "/home/ted/ubuntu/AMR/amromics-vis/dev/temp/cluster.clstr",
        "cd_hit_cluster_represent": "/home/ted/ubuntu/AMR/amromics-vis/dev/temp/cluster",
        "cluster_json": "/home/ted/ubuntu/AMR/amromics-vis/dev/temp/cluster.json",
        "combined_faa_file": "/home/ted/ubuntu/AMR/amromics-vis/dev/temp/combined.faa",
        "excluded_cluster_file": "/home/ted/ubuntu/AMR/amromics-vis/dev/temp/excluded_cluster",
        "inflated_unsplit_mcl_groups": "/home/ted/ubuntu/AMR/amromics-vis/dev/temp/inflated_unsplit_mcl_groups",
        "labeled_mcl_groups": "/home/ted/ubuntu/AMR/amromics-vis/dev/temp/labeled_mcl_groups",
        "main_dir": "/home/ted/ubuntu/AMR/amromics-vis/dev",
        "samples": [
            {
                "assembly_file": "/home/ted/ubuntu/AMR/amromics-vis/dev/BMP071.fna",
                "bed": "/home/ted/ubuntu/AMR/amromics-vis/dev/temp/BMP071.bed",
                "extracted_fna_file": "/home/ted/ubuntu/AMR/amromics-vis/dev/temp/BMP071.extracted.fna",
                "faa_file": "/home/ted/ubuntu/AMR/amromics-vis/dev/temp/BMP071.faa",
                "gff_file": "/home/ted/ubuntu/AMR/amromics-vis/dev/BMP071.gff",
                "id": "BMP071"
            },
            {
                "assembly_file": "/home/ted/ubuntu/AMR/amromics-vis/dev/BMP077.fna",
                "bed": "/home/ted/ubuntu/AMR/amromics-vis/dev/temp/BMP077.bed",
                "extracted_fna_file": "/home/ted/ubuntu/AMR/amromics-vis/dev/temp/BMP077.extracted.fna",
                "faa_file": "/home/ted/ubuntu/AMR/amromics-vis/dev/temp/BMP077.faa",
                "gff_file": "/home/ted/ubuntu/AMR/amromics-vis/dev/BMP077.gff",
                "id": "BMP077"
            }
        ],
        "temp_dir": "/home/ted/ubuntu/AMR/amromics-vis/dev/temp",
        "uninflated_mcl_groups": "/home/ted/ubuntu/AMR/amromics-vis/dev/temp/uninflated_mcl_groups"
    }
    threads = 4
    timing_log = '/home/ted/ubuntu/AMR/amromics-vis/dev/time.log'
    handle = extract_proteins(handle, timing_log=timing_log)
    handle = combine_proteins(handle, timing_log=timing_log)
    handle = run_cd_hit_iterative(handle, threads=threads, timing_log=timing_log)
    handle = all_against_all_blast(handle, threads=threads, timing_log=timing_log)
    handle = cluster_with_mcl(handle, threads=threads, timing_log=timing_log)
    handle = reinflate_clusters(handle)
    handle = label_group(handle)
    print(json.dumps(handle, indent=4, sort_keys=True))
    



## TO-DO ##
# handle or report ???
# check protein translate (compare with perl script)
    # filtered protein sequence produced by roary did not have * character in the middle of sequence => check unfiltered protein sequence 
    # protein sequence with the same id are different from each other => run roary perl scrip individually
# filter protein sequence
# warming of CD-HIT: Discarding invalid sequence or sequence without identifier and description!