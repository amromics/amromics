import os
import shutil
import re
import json
import gzip
import csv
from datetime import datetime
import logging
from Bio import SeqIO
import pandas as pd
from amromics.utils import run_command

logger = logging.getLogger(__name__)


def parse_gff_file(ggf_file, bed_out_file, fasta_out_file, sample_id, dictionary):
    """
    Parse gff file, create a bed file and return a dictionary contain gen id 
    and gene annotation. Filter out small genes(<18 base)

    Parameters
    -------
    -------
    """
    found_fasta = 0
    last_seq_id = None
    count = 1
    with gzip.open(ggf_file,'rt') as in_fh, open(bed_out_file, 'w') as bed_fh, open(fasta_out_file, 'w') as fasta_fh:
        for line in in_fh:
            if found_fasta == 1:
                fasta_fh.write(line)
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
            length = end - start
            if length < 18:
                continue

            seq_id = cells[0]
            
            #if seq_id != last_seq_id:
            #    count = 1
            #    last_seq_id = seq_id
            #gene_id = sample_id + '_' + seq_id + '_' + str(count)
            #count = count + 1
            
            gene_id = re.findall(r"ID=(.+?);",cells[8])
            gene_id = gene_id[0]
            #if gene_id in dictionary:
            #    gene_id = gene_id + datetime.now().strftime("%M%S%f")
            
            trand = cells[6]
            row = [seq_id, str(start -1), str(end), gene_id, '1', trand]
            bed_fh.write('\t'.join(row)+ '\n')

            # create annotation dictionary
            dictionary[gene_id] = {}
            dictionary[gene_id]['sample_id'] = sample_id
            dictionary[gene_id]['length'] = length
            gene_name = re.findall(r"Name=(.+?);",cells[8])
            if len(gene_name) != 0:
                gene_name = gene_name[0]
                dictionary[gene_id]['name'] = gene_name
            gene_product = re.findall(r"product=(.+?)$",cells[8])
            if len(gene_product) != 0:
                gene_product = gene_product[0]
                dictionary[gene_id]['product'] = gene_product
            dictionary[gene_id]['contig'] = seq_id


def extract_proteins(report, timing_log=None):
    """
    Take in GFF file and create protein sequences in FASTA format.
    Create json file contain annotation information of each gene_id

    Parameters
    -------
    -------
    """
    temp_dir = os.path.join(report['temp_dir'], 'samples')
    if not os.path.exists(temp_dir):
        os.makedirs(temp_dir)
    gene_annotation = {}
    for sample in report['samples']:
        sample_id = sample['id']
        sample_dir = os.path.join(temp_dir, sample_id)
        if not os.path.exists(sample_dir):
            os.makedirs(sample_dir)
        # parse gff file
        ggf_file = os.path.join(sample['annotation'], sample_id + '.gff.gz')
        fna_file = os.path.join(sample_dir, sample_id + '.fna')
        bed_file = os.path.join(sample_dir, sample_id + '.bed')
        parse_gff_file(ggf_file, bed_file, fna_file, sample_id, gene_annotation)
        # extract nucleotide region
        extracted_fna_file = os.path.join(sample_dir, sample_id +'.extracted.fna')
        cmd = f"bedtools getfasta -s -fi {fna_file} -bed {bed_file} -fo {extracted_fna_file} -name > /dev/null 2>&1"
        ret = run_command(cmd, timing_log)
        if ret != 0:
            raise Exception('Error running bedtools')
        # translate nucleotide to protein
        faa_file = os.path.join(sample_dir, sample_id +'.faa')
        with open(faa_file, 'w') as faa_hd:
            for seq_record in SeqIO.parse(extracted_fna_file, "fasta"):
                headers = seq_record.id.split(':')
                seq_record.id = headers[0]
                seq_record.name = ''
                seq_record.description = ''
                seq_record.seq = seq_record.seq.translate(table=11)
                SeqIO.write(seq_record, faa_hd, 'fasta')
        sample['bed'] = bed_file
        sample['extracted_fna_file'] = extracted_fna_file
        sample['faa_file'] = faa_file
    
    report['gene_annotation'] = gene_annotation
    return report


def combine_proteins(report, timing_log=None):
    """
    Take in multiple FASTA sequences containing proteomes and concat them together 
    and output a FASTA file

    Parameters
    -------
    -------
    """
    temp_dir = report['temp_dir']
    combined_faa_file = os.path.join(temp_dir, 'combined.faa')
    faa_file_list =[]
    for sample in report['samples']:
        faa_file = sample['faa_file']
        faa_file_list.append(faa_file)

    cmd = "cat {} > {}".format(" ".join(faa_file_list),combined_faa_file)
    ret = run_command(cmd, timing_log)
    if ret != 0:
        raise Exception('Error combining protein sequences')
    
    report['combined_faa_file'] = combined_faa_file
    return report


def parse_cluster_file(cd_hit_cluster_file):
    """
    Parse cd-hit .clstr file

    Parameters
    -------
    -------
    """
    
    clusters = {}
    with open(cd_hit_cluster_file, 'r') as fh:
        for line in fh:
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


def run_cd_hit_iterative(report, threads, timing_log):
    """
    Run CD-HIT iteratively

    Parameters
    -------
    -------
    """
    temp_dir = report['temp_dir']
    combined_faa_file = report['combined_faa_file']

    remain_faa_file = os.path.join(temp_dir, 'remain.faa')
    shutil.copyfile(combined_faa_file, remain_faa_file)
    report['remain_faa_file'] = remain_faa_file

    cd_hit_cluster_fasta = os.path.join(temp_dir, 'cluster')
    cd_hit_cluster_file = os.path.join(temp_dir, 'cluster.clstr')

    excluded_cluster = []

    greater_than_or_equal = True
    number_of_samples = len(report['samples'])
    
    lower = 0.98
    step = 0.005
    percent_match = 1
    while percent_match >= lower:
        cmd = f'cd-hit -i {remain_faa_file} -o {cd_hit_cluster_fasta} -s {percent_match} -c {percent_match} -T {threads} -M 0 -g 1 -d 256 > /dev/null'
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

        cluster_filtered_faa_file = remain_faa_file + '.filtered'
        with open(cluster_filtered_faa_file, 'w') as fh:
            for seq_record in SeqIO.parse(remain_faa_file, "fasta"):
                if seq_record.id not in full_cluster_gene_names:
                    SeqIO.write(seq_record, fh, 'fasta')
        shutil.move(cluster_filtered_faa_file, remain_faa_file)
        percent_match -= step

    cmd = f'cd-hit -i {remain_faa_file} -o {cd_hit_cluster_fasta} -s {lower} -c {lower} -T {threads} -M 0 -g 1 -d 256 > /dev/null'
    ret = run_command(cmd, timing_log)
    if ret != 0:
        raise Exception('Error running cd-hit')
    clusters = parse_cluster_file(cd_hit_cluster_file)

    report['cd_hit_cluster_fasta'] = cd_hit_cluster_fasta
    report['cd_hit_cluster_file'] = cd_hit_cluster_file
    report['excluded_cluster'] = excluded_cluster
    report['cd_hit_cluster'] = clusters
    return report


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
    chunked_fh = open(chunked_file, 'w')
    chunked_file_list.append(chunked_file)
    for seq_record in SeqIO.parse(fasta_file, "fasta"):
        if current_chunk_length > 200000:
            chunked_fh.close()
            chunk_number += 1
            current_chunk_length = 0
            chunked_file = os.path.join(out_dir, str(chunk_number) + '.seq')
            chunked_file_list.append(chunked_file)
            chunked_fh = open(chunked_file, 'w')
            SeqIO.write(seq_record, chunked_fh, 'fasta')
        chunked_file = os.path.join(out_dir, str(chunk_number) + '.seq')
        SeqIO.write(seq_record, chunked_fh, 'fasta')
        current_chunk_length += len(seq_record.seq)
    chunked_fh.close()
    return chunked_file_list


def all_against_all_blast(report, threads, timing_log):
    """
    Run all against all blast in parallel

    Parameters
    -------
    -------
    """
    out_dir = os.path.join(report['temp_dir'], 'blast')
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    # make blast database
    fasta_file = report['cd_hit_cluster_fasta']
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
    blast_output_file_list = []
    with open(blast_cmds_file,'w') as fh:
        for chunked_file in chunked_file_list:
            blast_output_file = os.path.splitext(chunked_file)[0] + '.out'
            blast_output_file_list.append(blast_output_file)
            cmd = f"blastp -query {chunked_file} -db {blast_db} -evalue 1E-6 -num_threads {threads} -outfmt 6 -max_target_seqs 2000 " + "| awk '{ if ($3 >= 98) print $0;}' 2> /dev/null 1> " + blast_output_file
            fh.write(cmd + '\n')
    cmd = f"parallel -a {blast_cmds_file}"
    ret = run_command(cmd, timing_log)
    if ret != 0:
        raise Exception('Error running parallel all-against-all blast')

    # combining blast results
    blast_result_file = os.path.join(report['temp_dir'], 'blast_results')
    if os.path.isfile(blast_result_file):
        os.remove(blast_result_file)
    for blast_output_file in blast_output_file_list:
        os.system(f'cat {blast_output_file} >> {blast_result_file}')

    report['blast_result_file'] = blast_result_file
    return report


def cluster_with_mcl(report, threads, timing_log):
    """
    Take blast results and outputs clustered results

    Parameters
    -------
    -------
    """
    out_dir = report['temp_dir']
    blast_results = report['blast_result_file']
    output_mcl_file = os.path.join(out_dir, 'uninflated_mcl_clusters')
    cmd = f"mcxdeblast -m9 --score r --line-mode=abc {blast_results} 2> /dev/null | mcl - --abc -I 1.5 -o {output_mcl_file} > /dev/null 2>&1"
    ret = run_command(cmd, timing_log)
    if ret != 0:
        raise Exception('Error running mcl')
    report['uninflated_mcl_clusters'] = output_mcl_file
    return report


def reinflate_clusters(report):
    """
    Take the clusters file from cd-hit and use it to reinflate the output of MCL

    Parameters
    -------
    -------
    """
    inflated_clusters = []
    clusters = report['cd_hit_cluster']

    # Inflate genes from cdhit which were sent to mcl
    mcl_file = report['uninflated_mcl_clusters']
    with open(mcl_file, 'r') as fh:
        for line in fh:
            inflated_genes = []
            line = line.rstrip('\n')
            genes = line.split('\t')
            for gene in genes:
                inflated_genes.append(gene)
                if gene in clusters:
                    inflated_genes.extend(clusters[gene])
                    del clusters[gene]
            inflated_clusters.append(inflated_genes)

    # Inflate any clusters that were in the clusters file but not sent to mcl
    for gene in clusters:
        inflated_genes = []
        inflated_genes.append(gene)
        inflated_genes.extend(clusters[gene])
        inflated_clusters.append(inflated_genes)

    # Add clusters which were excluded
    excluded_cluster = report['excluded_cluster']
    for cluster in excluded_cluster:
        inflated_clusters.append(cluster)

    report['inflated_unsplit_clusters'] = inflated_clusters
    return report

def find_paralogs(cluster, gene_annotation):
    samples = {}
    for gene_id in cluster:
        sample_id = gene_annotation[gene_id]['sample_id']
        if sample_id not in samples:
            samples[sample_id] = []
        samples[sample_id].append(gene_id)
    
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
        contig = gene_annotation[gene_id]['contig']
        neighbour_genes = []
        for i in range(before_position, after_position+1):
            neighbour_gene_id = gene_id_prefix + '{:0>5}'.format(i)
            if neighbour_gene_id in gene_annotation and neighbour_gene_id != gene_id and gene_annotation[neighbour_gene_id]['contig'] == contig:
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
            continue

        # find paralog gene which is closest match with g
        best_score = 0
        best_score_index = -1  # -1 is the index of "leftovers" list
        for p,_ in enumerate(paralog_genes):
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


def split_paralogs(report):
    gene_annotation = report['gene_annotation']
    unsplit_clusters = report['inflated_unsplit_clusters']

    clusters_not_paralogs = []

    # run iteratively
    out_clusters = unsplit_clusters
    for i in range(50):
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
            for index, genes in enumerate(in_clusters):
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
    report['split_clusters'] = split_clusters
    return report

def label_cluster(report):
    """
    Add labels to the cluster

    Parameters
    -------
    -------
    """
    unlabeled_clusters = report['split_clusters']

    # Add labels to the clusters
    labeled_clusters = {}
    counter = 1
    for cluster in unlabeled_clusters:
        labeled_clusters['groups_' + str(counter)] = cluster
        counter += 1

    report['labeled_clusters'] = labeled_clusters
    return report


def annotate_cluster(report):
    """
    Update the cluster name to the gene name

    Parameters
    -------
    -------
    """
    clusters = report['labeled_clusters']
    gene_annotation = report['gene_annotation']
    annotated_clusters = {}
    for cluster_name in clusters:
        cluster_new_name = cluster_name
        cluster_product = None
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
        if cluster_product == None:
            cluster_product =[]
            for gene_id in gene_id_list:
                if 'product' in gene_annotation[gene_id]:
                    gene_product = gene_annotation[gene_id]['product']
                    if gene_product not in cluster_product:
                        cluster_product.append(gene_product)
            if len(cluster_product) > 0:
                cluster_product = ', '.join(cluster_product)
            else:
                cluster_product = 'unknown'
        # check if cluster_new_name is already exist
        if cluster_new_name in annotated_clusters:
            cluster_new_name += '_' + datetime.now().strftime("%M%S%f")
        annotated_clusters[cluster_new_name] = {'gene_id':gene_id_list, 'product':cluster_product}
    report['annotated_clusters'] = annotated_clusters
    return report


def create_spreadsheet(report):
    spreadsheet_file = os.path.join(report['pan_genome'], 'gene_presence_absence.csv.gz')
    annotated_clusters = report['annotated_clusters']
    gene_annotation = report['gene_annotation']
    with gzip.open(spreadsheet_file, 'wt') as fh:
        writer = csv.writer(fh, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)

        # write header
        header = ['Gene', 'Non-unique Gene name', 'Annotation', 'No. isolates', 'No. sequences', 'Avg sequences per isolate', 'Genome Fragment','Order within Fragment', 'Accessory Fragment','Accessory Order with Fragment', 'QC','Min group size nuc', 'Max group size nuc', 'Avg group size nuc' ]
        for sample in report['samples']:
            header.append(sample['id'])
        writer.writerow(header)

        # write row
        for cluster in annotated_clusters:
            sample_dict = {}
            for gene in annotated_clusters[cluster]['gene_id']:
                sample_id = gene_annotation[gene]['sample_id']
                if sample_id not in sample_dict:
                    sample_dict[sample_id] = []
                sample_dict[sample_id].append(gene)
            
            row = []
            # Gene
            row.append(cluster)
            # Non-unique Gene name
            row.append("")
            # Annotation
            row.append(annotated_clusters[cluster]['product'])
            # No. isolates
            row.append(len(sample_dict))
            # No. sequences
            row.append(len(annotated_clusters[cluster]['gene_id']))
            # Avg sequences per isolate
            row.append("")
            # Genome Fragment
            row.append("")
            # Order within Fragment
            row.append("")
            # Accessory Fragment
            row.append("")
            # Accessory Order with Fragment
            row.append("")
            # QC
            row.append("")
            # Min group size nuc
            row.append("")
            # Max group size nuc
            row.append("")
            # Avg group size nuc
            row.append("")
            # sample columns
            for sample in report['samples']:
                sample_id = sample['id']
                if sample_id in sample_dict:
                    gene_list = sample_dict[sample_id]
                    row.append('\t'.join(gene_list))
                else:
                    row.append('')
            writer.writerow(row)
    report['spreadsheet'] = spreadsheet_file
    return report


def create_rtab(report):
    rtab_file = os.path.join(report['pan_genome'], 'gene_presence_absence.Rtab')
    annotated_clusters = report['annotated_clusters']
    gene_annotation = report['gene_annotation']
    with open(rtab_file, 'w') as fh:
        writer = csv.writer(fh, delimiter='\t')

        # write header
        header = ['Gene']
        for sample in report['samples']:
            header.append(sample['id'])
        writer.writerow(header)

        # write row
        for cluster in annotated_clusters:
            row = []
            # Gene
            row.append(cluster)
            # Samples
            sample_dict = {}
            for gene in annotated_clusters[cluster]['gene_id']:
                sample_id = gene_annotation[gene]['sample_id']
                if sample_id not in sample_dict:
                    sample_dict[sample_id] = []
                sample_dict[sample_id].append(gene)
            for sample in report['samples']:
                sample_id = sample['id']
                gene_list = sample_dict.get(sample_id, [])
                row.append(len(gene_list))
            writer.writerow(row)
    report['rtab'] = rtab_file
    return report


def create_summary(report):
    rtab_file = report['rtab']
    cluster_df = pd.read_csv(rtab_file, sep='\t', index_col='Gene')

    num_core = 0
    num_soft_core = 0
    num_shell = 0
    num_cloud = 0
    num_sample = len(report['samples'])
    for cluster, row in cluster_df.iterrows():
        absent = len(row[row == 0])
        percent = (num_sample - absent) / num_sample
        if percent >= 0.99:
            num_core += 1
        elif percent >= 0.95:
            num_soft_core += 1
        elif percent >= 0.15:
            num_shell += 1
        else:
            num_cloud += 1
    total = num_core + num_soft_core + num_shell + num_cloud

    summary_file = os.path.join(report['pan_genome'], 'summary_statistics.txt')
    with open(summary_file, 'w') as fh:
        fh.write('Core genes' + '\t' + '(99% <= strains <= 100%)' + '\t'+ str(num_core) + '\n')
        fh.write('Soft core genes' + '\t' + '(95% <= strains < 99%)' + '\t'+ str(num_soft_core) + '\n')
        fh.write('Shell genes' + '\t' + '(15% <= strains < 95%)' + '\t' + str(num_shell) + '\n')
        fh.write('Cloud genes' + '\t' + '(0% <= strains < 15%)' + '\t'+ str(num_cloud) + '\n')
        fh.write('Total genes' + '\t' + '(0% <= strains <= 100%)' + '\t'+ str(total))
    report['summary'] = summary_file
    return report


def create_representative_fasta(report):
    unsplit_clusters = report['inflated_unsplit_clusters']
    gene_annotation = report['gene_annotation']
    combined_fasta = report['combined_faa_file']
    representative_fasta = os.path.join(report['pan_genome'], 'representative.fasta')

    representative_list = []
    for cluster in unsplit_clusters:
        length_max = 0
        representative = None
        for gene_id in cluster:
            length = gene_annotation[gene_id]['length']
            if length > length_max:
                representative = gene_id
                length_max = length
        representative_list.append(representative)
    
    with open(representative_fasta, 'w') as fh:
        for seq_record in SeqIO.parse(combined_fasta, 'fasta'):
            if seq_record.id in representative_list:
                SeqIO.write(seq_record, fh, 'fasta')
    
    report['representative_fasta'] = representative_fasta
    return report


def run_pan_genome_analysis(report, collection_dir='.', threads=8, overwrite=False, timing_log=None):
    # Check if any sample has been updated
    for sample in report['samples']:
        overwrite = overwrite or sample['updated']

    pan_genome_folder = os.path.join(collection_dir, 'pan_genome')
    temp_dir = os.path.join(collection_dir, 'temp')
    pan_genome_output = os.path.join(pan_genome_folder,'summary_statistics.txt')
    report['pan_genome'] = pan_genome_folder
    report['temp_dir'] = temp_dir
    # Check if pan-genome has run
    if os.path.isfile(pan_genome_output) and (not overwrite):
        logger.info('Pan-genome has run and the input has not changed, skip')
        return report

    if not os.path.exists(pan_genome_folder):
        os.makedirs(pan_genome_folder)
    if not os.path.exists(temp_dir):
        os.makedirs(temp_dir)

    report = extract_proteins(report, timing_log=timing_log)
    report = combine_proteins(report, timing_log=timing_log)
    report = run_cd_hit_iterative(report, threads=threads, timing_log=timing_log)
    report = all_against_all_blast(report, threads=threads, timing_log=timing_log)
    report = cluster_with_mcl(report, threads=threads, timing_log=timing_log)
    report = reinflate_clusters(report)
    report = split_paralogs(report)
    report = label_cluster(report)
    report = annotate_cluster(report)
    report = create_spreadsheet(report)
    report = create_rtab(report)
    report = create_summary(report)
    report = create_representative_fasta(report)

    # clean
    #del report['cd_hit_cluster']
    #del report['excluded_cluster']
    #del report['inflated_unsplit_clusters']
    #del report['split_clusters']
    #del report['labeled_clusters']
    #del report['gene_annotation']
    #del report['annotated_clusters']
    #del report['blast_result_file']
    #del report['cd_hit_cluster_fasta']
    #del report['cd_hit_cluster_file']
    #del report['combined_faa_file']
    #del report['remain_faa_file']
    #del report['uninflated_mcl_clusters']
    #for sample in report['samples']:
    #    del sample['bed']
    #    del sample['extracted_fna_file']
    #    del sample['faa_file']
    #shutil.rmtree(temp_dir)

    return report

    
    
if __name__ == '__main__':
    report = json.load(open('dev/report1.json', 'r'))
    report = run_pan_genome_analysis(
        report, 
        collection_dir='dev/pan1', 
        threads=4, 
        overwrite=True
    )
    temp_dir = report['temp_dir']
    json.dump(report['gene_annotation'], open(temp_dir+'/gene_annotation', 'w'), indent=4, sort_keys=True)
    json.dump(report['cd_hit_cluster'], open(temp_dir+'/cd_hit_cluster', 'w'), indent=4, sort_keys=True)
    json.dump(report['excluded_cluster'], open(temp_dir+'/excluded_cluster', 'w'), indent=4, sort_keys=True)
    json.dump(report['inflated_unsplit_clusters'], open(temp_dir+'/inflated_unsplit_clusters', 'w'), indent=4, sort_keys=True)
    json.dump(report['split_clusters'], open(temp_dir+'/split_clusters', 'w'), indent=4, sort_keys=True)
    json.dump(report['labeled_clusters'], open(temp_dir+'/labeled_clusters', 'w'), indent=4, sort_keys=True)
    json.dump(report['annotated_clusters'], open(temp_dir+'/annotated_clusters', 'w'), indent=4, sort_keys=True)
    json.dump(report, open('dev/pan1/report1_output.json', 'w'), indent=4, sort_keys=True)

## TODO ##
# filter protein sequence?
# improve split paralogs