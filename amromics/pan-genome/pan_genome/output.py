import os
import csv
import sys
import gzip
import logging
import shutil
from datetime import datetime

from Bio import SeqIO

logger = logging.getLogger(__name__)

def classify_cluster(num_sample, total, count):
    percent = num_sample / total
    if percent >= 0.99:
        count[0] += 1 # core
    elif percent >= 0.95:
        count[1] += 1 # softcore
    elif percent >= 0.15:
        count[2] += 1 # shell
    else:
        count[3] += 1 # cloud

def write_summary(count_list, out_dir):
    total = sum(count_list)
    summary_file = os.path.join(out_dir, 'summary_statistics.txt')
    with open(summary_file, 'w') as fh:
        fh.write('Core genes' + '\t' + '(99% <= strains <= 100%)' 
                 + '\t'+ str(count_list[0]) + '\n')
        fh.write('Soft core genes' + '\t' + '(95% <= strains < 99%)' 
                 + '\t'+ str(count_list[1]) + '\n')
        fh.write('Shell genes' + '\t' + '(15% <= strains < 95%)' 
                 + '\t' + str(count_list[2]) + '\n')
        fh.write('Cloud genes' + '\t' + '(0% <= strains < 15%)' 
                 + '\t'+ str(count_list[3]) + '\n')
        fh.write('Total genes' + '\t' + '(0% <= strains <= 100%)' 
                 + '\t'+ str(total))

def create_output(clusters, clusters_annotation, 
                  gene_dictionary, samples, out_dir):
    starttime = datetime.now()

    cluster_info_file = os.path.join(out_dir, 'cluster_info.csv')
    presence_absence_file = os.path.join(
        out_dir, 'gene_presence_absence.csv.gz')
    
    # 4 value corespond to count of core, softcore, shell, cloud gene
    count_list = [0] * 4 

    with open(cluster_info_file, 'w') as fh_1, \
         gzip.open(presence_absence_file, 'wt') as fh_2:
        writer_1 = csv.writer(
            fh_1, delimiter=',', quotechar='"', quoting=csv.QUOTE_ALL)
        writer_2 = csv.writer(
            fh_2, delimiter=',', quotechar='"', quoting=csv.QUOTE_ALL)

        # write header
        header_1 = ['ID','Gene', 'Annotation', 'No. isolates', 'No. sequences', 
                    'Avg sequences per isolate', 'Min group size nuc', 
                    'Max group size nuc', 'Avg group size nuc' ]
        writer_1.writerow(header_1)
        
        header_2 = ['ID']
        for sample in samples:
            header_2.append(sample['id'])
        writer_2.writerow(header_2)
        total_sample = len(samples)

        # write row
        cluster_id = 0
        for cluster, cluster_annotation in zip(clusters, clusters_annotation):
            row_1 = []
            sample_dict = {}
            length_list = []
            for gene in cluster:
                sample_id = gene_dictionary[gene][0]
                sample_dict.setdefault(sample_id, []).append(gene)
                length = gene_dictionary[gene][2]
                length_list.append(length)
            
            # ID
            row_1.append(cluster_id)
            # Gene
            row_1.append(cluster_annotation[0])
            # Annotation
            row_1.append(cluster_annotation[1])
            # No. isolates
            num_sample = len(sample_dict)
            row_1.append(num_sample)
            classify_cluster(num_sample, total_sample, count_list)
            # No. sequences
            row_1.append(len(cluster))
            # Avg sequences per isolate
            avg_seq = len(cluster) / len(sample_dict)
            row_1.append(round(avg_seq,2))
            # Min group size nuc
            row_1.append(min(length_list))
            # Max group size nuc
            row_1.append(max(length_list))
            # Avg group size nuc
            nuc_size = sum(length_list) / len(length_list)
            row_1.append(round(nuc_size,0))

            writer_1.writerow(row_1)

            # samples
            row_2 = []
            row_2.append(cluster_id)
            for sample in samples:
                sample_id = sample['id']
                if sample_id in sample_dict:
                    gene_list = sample_dict[sample_id]
                    row_2.append('\t'.join(gene_list))
                else:
                    row_2.append('')
            writer_2.writerow(row_2)

            cluster_id += 1
    
    write_summary(count_list, out_dir)

    elapsed = datetime.now() - starttime
    logging.info(f'Create output -- time taken {str(elapsed)}')


def update_output(
        old_clusters, new_clusters, new_clusters_annotation, 
        gene_dictionary, new_samples, temp_dir, collection_dir):
    starttime = datetime.now()
    new_cluster_info_file = os.path.join(
        temp_dir, 'cluster_info.csv')
    new_presence_absence_file = os.path.join(
        temp_dir, 'gene_presence_absence.csv.gz')
    old_cluster_info_file = os.path.join(
        collection_dir, 'cluster_info.csv')
    old_presence_absence_file = os.path.join(
        collection_dir, 'gene_presence_absence.csv.gz')
    
    # 4 value corespond to count of core, softcore, shell, cloud gene
    count_list = [0] * 4 

    with open(old_cluster_info_file, 'r') as in_fh_1, \
         gzip.open(old_presence_absence_file, 'rt') as in_fh_2:
        with open(new_cluster_info_file, 'w') as out_fh_1, \
             gzip.open(new_presence_absence_file, 'wt') as out_fh_2:
            csv.field_size_limit(sys.maxsize)
            writer_1 = csv.writer(
                out_fh_1, delimiter=',', quotechar='"', quoting=csv.QUOTE_ALL)
            writer_2 = csv.writer(
                out_fh_2, delimiter=',', quotechar='"', quoting=csv.QUOTE_ALL)
            reader_1 = csv.reader(in_fh_1, delimiter=',')
            reader_2 = csv.reader(in_fh_2, delimiter=',')
            
            # update header
            header_1 = next(reader_1)
            header_2 = next(reader_2)
            num_old_samples = len(header_2[1:]) # dont count id column
            for sample in new_samples:
                header_2.append(sample['id'])
            writer_1.writerow(header_1)
            writer_2.writerow(header_2)
            total_sample = len(header_2) - 1 # dont count id column
            
            # update old row
            cluster_id = 0
            for row_1, row_2, cluster in zip(reader_1, reader_2, old_clusters):
                num_seq = int(row_1[4])
                num_iso = int(row_1[3])
                min_len = int(row_1[6])
                max_len = int(row_1[7])
                new_sample_dict = {}
                new_length_list = []
                for new_gene in cluster:
                    sample_id = gene_dictionary[new_gene][0]
                    new_sample_dict.setdefault(sample_id, []).append(new_gene)
                    length = gene_dictionary[new_gene][2]
                    if length > max_len:
                        max_len = length
                    if length < min_len:
                        min_len = length
                    new_length_list.append(length)
                
                # No. isolates
                num_sample = num_iso + len(new_sample_dict)
                row_1[3] = num_sample
                classify_cluster(num_sample, total_sample, count_list)
                # No. sequences
                row_1[4] = num_seq + len(cluster)
                # Avg sequences per isolate
                avg_seq = row_1[4] / row_1[3]
                row_1[5] = round(avg_seq,2)
                # Min group size nuc
                row_1[6] = min_len
                # Max group size nuc
                row_1[7] = max_len
                # Avg group size nuc
                nuc_size = (
                    (float(row_1[8]) * num_seq + len(new_length_list)) 
                    / row_1[4])
                row_1[8] = round(nuc_size,0)
                
                writer_1.writerow(row_1)

                for sample in new_samples:
                    sample_id = sample['id']
                    if sample_id in new_sample_dict:
                        gene_list = new_sample_dict[sample_id]
                        row_2.append('\t'.join(gene_list))
                    else:
                        row_2.append('')
                
                writer_2.writerow(row_2)
                cluster_id += 1


            # write new row
            for cluster, cluster_annotation in zip(new_clusters, 
                                                   new_clusters_annotation):
                row_1 = []
                sample_dict = {}
                length_list = []
                for gene in cluster:
                    sample_id = gene_dictionary[gene][0]
                    sample_dict.setdefault(sample_id, []).append(gene)
                    length = gene_dictionary[gene][2]
                    length_list.append(length)
                
                # ID
                row_1.append(cluster_id)
                cluster_id += 1
                # Gene
                row_1.append(cluster_annotation[0])
                # Annotation
                row_1.append(cluster_annotation[1])
                # No. isolates
                num_sample = len(sample_dict)
                row_1.append(num_sample)
                classify_cluster(num_sample, total_sample, count_list)
                # No. sequences
                row_1.append(len(cluster))
                # Avg sequences per isolate
                avg_seq = len(cluster) / len(sample_dict)
                row_1.append(round(avg_seq,2))
                # Min group size nuc
                row_1.append(min(length_list))
                # Max group size nuc
                row_1.append(max(length_list))
                # Avg group size nuc
                nuc_size = sum(length_list) / len(length_list)
                row_1.append(round(nuc_size,0))
                
                writer_1.writerow(row_1)
                
                # sample
                row_2 = []
                row_2.append(cluster_id)
                for i in range(0, num_old_samples):
                    row_2.append('')
                for sample in new_samples:
                    sample_id = sample['id']
                    if sample_id in sample_dict:
                        gene_list = sample_dict[sample_id]
                        row_2.append('\t'.join(gene_list))
                    else:
                        row_2.append('')
                
                writer_2.writerow(row_2)
    
    shutil.move(new_cluster_info_file, old_cluster_info_file)
    shutil.move(new_presence_absence_file, old_presence_absence_file)
    
    write_summary(count_list, collection_dir)

    elapsed = datetime.now() - starttime
    logging.info(f'Update output -- time taken {str(elapsed)}')


def create_rtab(clusters, gene_dictionary, samples, out_dir):
    starttime = datetime.now()
    rtab_file = os.path.join(out_dir, 'gene_presence_absence.Rtab.gz')
    with gzip.open(rtab_file, 'wt') as fh:
        writer = csv.writer(fh, delimiter='\t')

        # write header
        header = ['ID']
        for sample in samples:
            header.append(sample['id'])
        writer.writerow(header)

        # write row
        cluster_id = 0
        for cluster in clusters:
            row = []
            # ID
            row.append(cluster_id)
            cluster_id += 1
            # Samples
            sample_dict = {}
            for gene in cluster:
                sample_id = gene_dictionary[gene][0]
                sample_dict.setdefault(sample_id, []).append(gene)
            for sample in samples:
                sample_id = sample['id']
                gene_list = sample_dict.get(sample_id, [])
                row.append(len(gene_list))
            writer.writerow(row)
    elapsed = datetime.now() - starttime
    logging.info(f'Create Rtab -- time taken {str(elapsed)}')
    return rtab_file

def update_rtab(old_file, old_clusters, new_clusters, 
                gene_dictionary, new_samples, temp_dir):
    starttime = datetime.now()
    new_file = os.path.join(temp_dir, 'gene_presence_absence.Rtab.gz')
    with gzip.open(new_file, 'wt') as out_fh, \
         gzip.open(old_file, 'rt') as in_fh:
        csv.field_size_limit(sys.maxsize)
        writer = csv.writer(out_fh, delimiter='\t')
        reader = csv.reader(in_fh, delimiter='\t')
        
        # update header
        header = next(reader)
        num_old_samples = len(header[1:])
        for sample in new_samples:
            header.append(sample['id'])
        writer.writerow(header)
        
        # update old row
        cluster_id = 0
        for row, cluster in zip(reader, old_clusters):
            cluster_id += 1
            new_sample_dict = {}
            for gene in cluster:
                sample_id = gene_dictionary[gene][0]
                new_sample_dict.setdefault(sample_id, []).append(gene)
    
            for sample in new_samples:
                sample_id = sample['id']
                gene_list = new_sample_dict.get(sample_id, [])
                row.append(len(gene_list))

            writer.writerow(row)

        # write new row
        for cluster in new_clusters:
            row = []
            # ID
            row.append(cluster_id)
            cluster_id += 1
            # Samples
            sample_dict = {}
            for gene in cluster:
                sample_id = gene_dictionary[gene][0]
                sample_dict.setdefault(sample_id, []).append(gene)
            for i in range(0, num_old_samples):
                row.append(0)
            for sample in new_samples:
                sample_id = sample['id']
                gene_list = sample_dict.get(sample_id, [])
                row.append(len(gene_list))
            writer.writerow(row)

    shutil.move(new_file, old_file)
    elapsed = datetime.now() - starttime
    logging.info(f'Update Rtab -- time taken {str(elapsed)}')
    return old_file


def create_summary(rtab_file, out_dir):
    starttime = datetime.now()
    num_core = 0
    num_soft_core = 0
    num_shell = 0
    num_cloud = 0
    with gzip.open(rtab_file, 'rt') as fh:
        for line in fh:
            line = line.rstrip()
            cells = line.split('\t')
            if cells[0] == 'Gene':
                num_sample = len(cells) - 1
                continue
            num_zero = 0
            for cell in cells:
                if cell == '0':
                    num_zero += 1
            percent = (num_sample - num_zero) / num_sample
            if percent >= 0.99:
                num_core += 1
            elif percent >= 0.95:
                num_soft_core += 1
            elif percent >= 0.15:
                num_shell += 1
            else:
                num_cloud += 1
    total = num_core + num_soft_core + num_shell + num_cloud

    summary_file = os.path.join(out_dir, 'summary_statistics.txt')
    with open(summary_file, 'w') as fh:
        fh.write('Core genes' + '\t' + '(99% <= strains <= 100%)' 
                 + '\t'+ str(num_core) + '\n')
        fh.write('Soft core genes' + '\t' + '(95% <= strains < 99%)' 
                 + '\t'+ str(num_soft_core) + '\n')
        fh.write('Shell genes' + '\t' + '(15% <= strains < 95%)' 
                 + '\t' + str(num_shell) + '\n')
        fh.write('Cloud genes' + '\t' + '(0% <= strains < 15%)' 
                 + '\t'+ str(num_cloud) + '\n')
        fh.write('Total genes' + '\t' + '(0% <= strains <= 100%)' 
                 + '\t'+ str(total))
    elapsed = datetime.now() - starttime
    logging.info(f'Create summary -- time taken {str(elapsed)}')
    return summary_file


def create_representative_fasta(
        gene_to_cluster, clusters_annotation, 
        source_fasta, out_fasta, mode='w'):
    for seq_record in SeqIO.parse(open(source_fasta), "fasta"):
        gene_id = seq_record.id
        if gene_id in gene_to_cluster:
            index = gene_to_cluster[gene_id]
            annotation = clusters_annotation[index]
            seq_record.description = "~~~".join(annotation)
            clusters_annotation[index].append(seq_record)

    with open(out_fasta, mode) as out_fh:
        for cluster in clusters_annotation:
            seq_record=cluster[2]
            SeqIO.write(seq_record, out_fh, 'fasta')

# def create_representative_fasta(clusters, gene_dictionary,
#                                 fasta_list, out_dir):
#     starttime = datetime.now()
#     representative_fasta = os.path.join(out_dir, 'representative.fasta')
#     representative_list = set()
#     for cluster in clusters:
#         length_max = 0
#         representative = None
#         for gene_id in cluster:
#             length = gene_dictionary[gene_id][2]
#             if length > length_max:
#                 representative = gene_id
#                 length_max = length
#         representative_list.add(representative)
#     utils.create_fasta_include(
#         fasta_file_list=fasta_list, 
#         include_list=representative_list, 
#         output_file=representative_fasta
#         )
#     elapsed = datetime.now() - starttime
#     logging.info(f'Create representative fasta -- time taken {str(elapsed)}')
    
#     return representative_fasta


def write_gene_dictionary(gene_dictionary, out_dir, mode='w'):
    # starttime = datetime.now()
    
    with open(os.path.join(out_dir, 'gene_dictionary.tsv'),mode) as fh:
        writer = csv.writer(fh, delimiter='\t')
        for gene in gene_dictionary:
            row = []
            row.append(gene)
            row.extend(gene_dictionary[gene])
            writer.writerow(row)
    
    # elapsed = datetime.now() - starttime
    # logging.info(f'Export gene annotation -- time taken {str(elapsed)}')

def read_gene_dictionary(annotation_file):
    # starttime = datetime.now()
    
    gene_dictionary = {}
    with open(annotation_file,'r') as fh:
        csv_reader = csv.reader(fh, delimiter='\t')
        for row in csv_reader:
            row[3] = int(row[3])
            gene_dictionary[row[0]] = tuple(row[1:])
    
    # elapsed = datetime.now() - starttime
    # logging.info(f'Import gene annotation -- time taken {str(elapsed)}')
    return gene_dictionary

def write_gene_position(gene_position, out_dir, mode='w'):
    # starttime = datetime.now()
    
    with open(os.path.join(out_dir, 'gene_position.tsv'),mode) as fh:
        writer = csv.writer(fh, delimiter='\t')
        for sample in gene_position:
            for seq in gene_position[sample]:
                row = []
                row.append(sample)
                row.append(seq)
                row.append(';'.join(gene_position[sample][seq]))
                writer.writerow(row)
    
    # elapsed = datetime.now() - starttime
    # logging.info(f'Export gene annotation -- time taken {str(elapsed)}')

def read_gene_position(position_file):
    # starttime = datetime.now()
    
    gene_position = {}
    with open(position_file,'r') as fh:
        csv_reader = csv.reader(fh, delimiter='\t')
        for row in csv_reader:
            ls = row[1].split(';')
            gene_position[row[0]] = ls
    
    # elapsed = datetime.now() - starttime
    # logging.info(f'Import gene position -- time taken {str(elapsed)}')
    return gene_position