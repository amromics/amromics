import os
import shutil
import re
import json
#from datetime import datetime
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


def parse_gff_file(ggf_file, bed_out_file, sample_id, dictionary):
    """
    Parse gff file, create a bed file and return a dictionary contain gen id 
    and gene annotation. Filter out small genes(<18 base)

    Parameters
    -------
    -------
    """
    in_handle = open(ggf_file,'r')
    out_handle = open(bed_out_file, 'w')
    found_fasta = 0
    last_seq_id = None
    count = 1
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

        seq_id = cells[0]
        if seq_id != last_seq_id:
            count = 1
            last_seq_id = seq_id
        #gene_id = re.findall(r"ID=(.+?);",cells[8])
        #gene_id = gene_id[0]
        gene_id = sample_id + '_' + seq_id + '_' + str(count)
        count = count + 1
        trand = cells[6]
        out_handle.write(seq_id+'\t'+str(start -1)+'\t'+ str(end)+'\t'+gene_id+'\t'+trand+'\n')\

        # create annotation dictionary
        #if gene_id in dictionary:
        #    gene_id = gene_id + datetime.now().strftime("%M%S%f")
        dictionary[gene_id] = {}
        dictionary[gene_id]['sample_id'] = sample_id
        gene_name = re.findall(r"Name=(.+?);",cells[8])
        if len(gene_name) != 0:
            gene_name = gene_name[0]
            dictionary[gene_id]['name'] = gene_name
        gene_product = re.findall(r"product=(.+?)$",cells[8])
        if len(gene_product) != 0:
            gene_product = gene_product[0]
            dictionary[gene_id]['product'] = gene_product

    in_handle.close()
    out_handle.close()


def extract_proteins(handle, timing_log=None):
    """
    Take in GFF file and fasta file and create protein sequences in FASTA format.
    Create json file contain annotation information of genes

    Parameters
    -------
    -------
    """
    temp_dir = handle['temp_dir']
    gene_annotation = {}
    for sample in handle['samples']:
        fna_file = sample['assembly_file']
        ggf_file = sample['gff_file']
        sample_id = sample['id']
        # extract nucleotide region
        extracted_fna_file = os.path.join(temp_dir, sample_id +'.extracted.fna')
        bed_file = os.path.join(temp_dir, sample_id +'.bed')
        parse_gff_file(ggf_file, bed_file, sample_id, gene_annotation)

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
    
    handle['gene_annotation'] = gene_annotation
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
    for cluster_name in clusters:
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

    excluded_cluster = []

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
        for cluster_represent in clusters:
            other_genes = clusters[cluster_represent]
            this_cluster = []
            if greater_than_or_equal == True and len(other_genes) >= number_of_samples -1:
                this_cluster.append(cluster_represent)
                this_cluster.extend(other_genes)
            if greater_than_or_equal == False and len(other_genes) == number_of_samples -1:
                this_cluster.append(cluster_represent)
                this_cluster.extend(other_genes)
            if len(this_cluster) != 0:
                full_cluster_gene_names.extend(this_cluster)
                excluded_cluster.append(this_cluster)

        cluster_filtered_faa_file = combined_faa_file + '.filtered'
        cluster_filtered_faa_fh = open(cluster_filtered_faa_file, 'w')
        for seq_record in SeqIO.parse(combined_faa_file, "fasta"):
            if seq_record.id in full_cluster_gene_names:
                SeqIO.write(seq_record, cluster_filtered_faa_fh, 'fasta')
        cluster_filtered_faa_fh.close()

        shutil.move(cluster_filtered_faa_file, combined_faa_file)
        percent_match -=0.005

    handle['cd_hit_cluster_represent'] = cd_hit_cluster_represent
    handle['cd_hit_cluster_file'] = cd_hit_cluster_file
    handle['excluded_cluster'] = excluded_cluster
    handle['cd_hit_cluster'] = clusters
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
    output_mcl_file = os.path.join(out_dir, 'uninflated_mcl_clusters')
    cmd = f"mcxdeblast -m9 --score r --line-mode=abc {blast_results} 2> /dev/null | mcl - --abc -I 1.5 -o {output_mcl_file} > /dev/null 2>&1"
    ret = run_command(cmd, timing_log)
    if ret != 0:
        raise Exception('Error running mcl')
    handle['uninflated_mcl_clusters'] = output_mcl_file
    return handle


def reinflate_clusters(handle):
    """
    Take the clusters file from cd-hit and use it to reinflate the output of MCL

    Parameters
    -------
    -------
    """
    inflated_clusters = []
    clusters = handle['cd_hit_cluster']

    # Inflate genes from cdhit which were sent to mcl
    mcl_file = handle['uninflated_mcl_clusters']
    mcl_fh = open(mcl_file, 'r')
    for line in mcl_fh:
        inflated_genes = []
        line = line.rstrip('\n')
        genes = line.split('\t')
        for gene in genes:
            inflated_genes.append(gene)
            if gene in clusters:
                inflated_genes.extend(clusters[gene])
                del clusters[gene]
        inflated_clusters.append(inflated_genes)
    mcl_fh.close()

    # Inflate any clusters that were in the clusters file but not sent to mcl
    for gene in clusters:
        inflated_genes = []
        inflated_genes.append(gene)
        inflated_genes.extend(cluster[gene])
        inflated_clusters.append(inflated_genes)

    # Add clusters which were excluded
    excluded_cluster = handle['excluded_cluster']
    for cluster in excluded_cluster:
        inflated_clusters.append(cluster)

    handle['inflated_unsplit_clusters'] = inflated_clusters
    del handle['cd_hit_cluster']
    del handle['excluded_cluster']
    return handle

def find_paralogs(cluster, gene_annotation):
    samples = {}
    for gene_id in cluster:
        sample_id = gene_annotation[gene_id]['sample_id']
        if sample_id not in samples:
            samples[sample_id] = []
        samples[sample_id].append(sample_id)
    
    # pick paralogs with the smallest number of genes
    smallest_number = 1000000
    paralog_genes = None
    for sample_id in samples:
        genes = samples[sample_id]
        count = len(genes)
        if count > 1 and count < smallest_number:
            paralog_genes = genes
            smallest_number = count

    return paralog_genes


def create_orthologs(cluster, paralog_genes, gene_annotation, gene_to_cluster_index):
    # get neighbour gene
    neighbour_gene_dictionary = {}
    for gene_id in cluster:
        match = re.search(r'^(.+_)(\d+)$', gene_id)
        gene_id_prefix = match.group(1)
        gene_position = int(match.group(2))
        before_position = gene_position - 5
        after_position = gene_position + 5
        neighbour_genes = []
        for i in range(before_position, after_position+1):
            neighbour_gene_id = gene_id_prefix + str(i)
            if neighbour_gene_id in gene_annotation and neighbour_gene_id != gene_id:
                neighbour_genes.append(neighbour_gene_id)
        neighbour_gene_dictionary[gene_id] = neighbour_genes

    # find cluster indices of all the neighbour genes of each paralog gene
    cluster_indices_around_paralogs = []
    for p in paralog_genes:
        neighbours_of_p = neighbour_gene_dictionary[p]
        cluster_indices_around_p = []
        for neighbour_gene in neighbours_of_p:
            cluster_index = gene_to_cluster_index[neighbour_gene]
            cluster_indices_around_p.append(cluster_index)
        cluster_indices_around_paralogs.append(cluster_indices_around_p)
    
    # create data structure to hold new clusters
    new_clusters = []
    for p in paralog_genes:
        new_clusters.append([p])
    new_clusters.append([]) # extra "leftovers" list to gather genes that don't share CGN with any paralog gene

    # add other members of the cluster to their closest match
    for g in cluster:
        if g in paralog_genes:
            continue

        neighbour_genes_of_g = neighbour_gene_dictionary[g]
        if len(neighbour_genes_of_g) == 0:
            new_clusters[-1].append(g)

        # find paralog gene which is closest match with g
        best_score = 0
        best_score_index = -1  # -1 is the index of "leftovers" list
        for p in range(len(paralog_genes)):
            cluster_indices_around_p = cluster_indices_around_paralogs[p]
            score_of_p = 0
            for neighbour_gene in neighbour_genes_of_g:
                cluster_index = gene_to_cluster_index[neighbour_gene]
                if cluster_index in cluster_indices_around_p:
                    score_of_p += 1
            score_of_p = score_of_p / len(neighbour_genes_of_g)
            if score_of_p > best_score:
                best_score = score_of_p
                best_score_index = p

            new_clusters[best_score_index].append(g)

    # check for "leftovers", remove if absent
    if len(new_clusters[-1]) == 0:
        del new_clusters[-1]
    
    return new_clusters


def split_paralogs(handle):
    gene_annotation = handle['gene_annotation']
    unsplit_clusters = handle['inflated_unsplit_clusters']

    clusters_not_paralogs = []

    # run iteratively
    out_clusters = unsplit_clusters
    for i in range(5):
        in_clusters = out_clusters
        out_clusters = []
        any_paralogs = 0
        for cluster in in_clusters:
            if len(cluster) == 1:
                out_clusters.append(cluster)
                continue
            first_gene = cluster[0]
            if first_gene in clusters_not_paralogs:
                out_clusters.append(cluster)
                continue

            # check paralogs
            paralog_genes = find_paralogs(cluster, gene_annotation)

            if paralog_genes == None:
                clusters_not_paralogs.append(first_gene)
                out_clusters.append(cluster)
                continue
            
            # convert in_clusters so we can find the cluster index of gene
            gene_to_cluster_index = {}
            for index in range(len(in_clusters)):
                genes = in_clusters[index]
                for gene in genes:
                    gene_to_cluster_index[gene] = index

            # split paralogs
            orthologs_clusters = create_orthologs(cluster, paralog_genes, gene_annotation, gene_to_cluster_index)
            out_clusters.extend(orthologs_clusters)
            any_paralogs = 1

        # check if next iteration is required
        if any_paralogs == 0:
            break

    split_clusters = out_clusters
    handle['split_clusters'] = split_clusters
    del handle['inflated_unsplit_clusters']
    return handle

def label_cluster(handle):
    """
    Add labels to the cluster

    Parameters
    -------
    -------
    """
    unlabeled_clusters = handle['split_clusters']

    # Add labels to the clusters
    labeled_clusters = {}
    counter = 1
    for cluster in unlabeled_clusters:
        labeled_clusters['groups_' + str(counter)] = cluster
        counter += 1

    handle['labeled_clusters'] = labeled_clusters
    del handle['split_clusters']
    return handle


def annotate_cluster(handle):
    """
    Update the cluster name to the gene name

    Parameters
    -------
    -------
    """
    clusters = handle['labeled_clusters']
    gene_annotation = handle['gene_annotation']
    annotated_cluster = {}
    for cluster_name in clusters:
        cluster_new_name = cluster_name
        cluster_product = ''
        gene_name_count = {}
        max_number = 0
        gene_id_list = clusters[cluster_name]
        for gene_id in gene_id_list:
            if 'name' in gene_annotation[gene_id]:
                gene_name = gene_annotation[gene_id]['name']
                gene_name_count[gene_name] = gene_name_count.get(gene_name, 1) + 1
                if gene_name_count[gene_name] > max_number:
                    cluster_new_name = gene_name
                    max_number = gene_name_count[gene_name]
                    if 'product' in gene_annotation[gene_id]:
                        cluster_product = gene_annotation[gene_id]['product']
        if cluster_product == '':
            cluster_product =[]
            for gene_id in gene_id_list:
                if 'product' in gene_annotation[gene_id]:
                    gene_product = gene_annotation[gene_id]['product']
                    if gene_product not in cluster_product:
                        cluster_product.append(gene_product)
            if len(cluster_product) > 0:
                cluster_product = ', '.join(cluster_product)
            else:
                cluster_product = 'hypothetical protein'
        annotated_cluster[cluster_new_name] = {'gene_id':gene_id_list, 'product':cluster_product}
    del handle['labeled_clusters']
    del handle['gene_annotation']
    handle['annotated_cluster'] = annotated_cluster
    return handle

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
    handle = split_paralogs(handle)
    handle = label_cluster(handle)
    handle = annotate_cluster(handle)
    


if __name__ == '__main__':
    handle = {
        "blast_result_file": "/home/ted/ubuntu/AMR/amromics-vis/dev/temp/blast_results",
        "cd_hit_cluster_file": "/home/ted/ubuntu/AMR/amromics-vis/dev/temp/cluster.clstr",
        "cd_hit_cluster_represent": "/home/ted/ubuntu/AMR/amromics-vis/dev/temp/cluster",
        "combined_faa_file": "/home/ted/ubuntu/AMR/amromics-vis/dev/temp/combined.faa",
        "gene_annotation": "/home/ted/ubuntu/AMR/amromics-vis/dev/temp/annotation.json",
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
        "uninflated_mcl_clusters": "/home/ted/ubuntu/AMR/amromics-vis/dev/temp/uninflated_mcl_clusters",
    }
    threads = 4
    timing_log = '/home/ted/ubuntu/AMR/amromics-vis/dev/time.log'
    handle = extract_proteins(handle, timing_log=timing_log)
    json.dump(handle['gene_annotation'], open('gene_annotation', 'w'), indent=4, sort_keys=True)
    handle = combine_proteins(handle, timing_log=timing_log)
    handle = run_cd_hit_iterative(handle, threads=threads, timing_log=timing_log)
    handle = all_against_all_blast(handle, threads=threads, timing_log=timing_log)
    handle = cluster_with_mcl(handle, threads=threads, timing_log=timing_log)
    json.dump(handle['excluded_cluster'], open('excluded_cluster', 'w'), indent=4, sort_keys=True)
    json.dump(handle['cd_hit_cluster'], open('cd_hit_cluster', 'w'), indent=4, sort_keys=True)
    handle = reinflate_clusters(handle)
    json.dump(handle['inflated_unsplit_clusters'], open('inflated_unsplit_clusters', 'w'), indent=4, sort_keys=True)
    handle = split_paralogs(handle)
    json.dump(handle['split_clusters'], open('split_clusters', 'w'), indent=4, sort_keys=True)
    handle = label_cluster(handle)
    json.dump(handle['labeled_clusters'], open('labeled_clusters', 'w'), indent=4, sort_keys=True)
    handle = annotate_cluster(handle)
    json.dump(handle['annotated_cluster'], open('annotated_cluster', 'w'), indent=4, sort_keys=True)
    
    #print(json.dump(handle['annotated_cluster'], 'annotated_cluster', indent=4, sort_keys=True))
    #print(json.dump(handle['annotated_cluster'], 'annotated_cluster' indent=4, sort_keys=True))



## TO-DO ##
# handle or report ???
# input from gene bank
# check protein translate (compare with perl script)
    # filtered protein sequence produced by roary did not have * character in the middle of sequence => check unfiltered protein sequence 
    # protein sequence with the same id are different from each other => run roary perl scrip individually
# filter protein sequence
# warming of CD-HIT: Discarding invalid sequence or sequence without identifier and description!