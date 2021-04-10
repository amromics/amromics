import os
import shutil
import re
import json
import gzip
import csv
import sqlite3
from datetime import datetime
import logging
from Bio import SeqIO
import pandas as pd
from amromics.utils import run_command

logger = logging.getLogger(__name__)



def parse_gff_file(ggf_file, bed_out_file, fasta_out_file, sample_name):
    """
    Parse gff file, filter out small genes(<18 base). Create a bed file, 
    fasta file and a dictionary contain gene annotation. 

    Parameters
    -------
    -------
    """
    found_fasta = 0
    last_contig_name = None
    pos = 0
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

            contig_name = cells[0]
            if contig_name != last_contig_name:
                pos = 0
                last_contig_name = contig_name
            pos += 1
            
            cds_name = re.findall(r"ID=(.+?);",cells[8])
            cds_name = cds_name[0]
            
            trand = cells[6]
            row = [contig_name, str(start -1), str(end), cds_name, '1', trand]
            bed_fh.write('\t'.join(row)+ '\n')

            # Contig table
            cur.execute('SELECT id FROM Contig WHERE name = ? LIMIT 1', (contig_name,))
            try:
                contig_id = cur.fetchone()[0]
            except:
                cur.execute('INSERT OR IGNORE INTO Contig (name) VALUES (?)', (contig_name,))
                conn.commit()
                if cur.rowcount != 1 :
                    print('Error inserting contig: ',contig_name)
                    continue
                contig_id = cur.lastrowid
            
            # Cds table  
            cur.execute('SELECT id FROM Sample WHERE name = ? LIMIT 1', (sample_name,))
            sample_id = cur.fetchone()[0]

            cur.execute('SELECT id FROM Cds WHERE name = ? LIMIT 1', (cds_name,))
            try:
                cds_id = cur.fetchone()[0]
                cds_name += '_' + datetime.now().strftime("%M%S%f")
            except:
                pass
            cur.execute('''INSERT OR IGNORE INTO Cds 
                            (name, sample_id, contig_id, position, trand, start, end, length) 
                            VALUES (?,?,?,?,?,?,?,?)''', 
                            (cds_name, sample_id, contig_id, pos, trand, start, end, length)
                        )
            conn.commit()
            if cur.rowcount != 1 :
                print('Error inserting cds: ',cds_name)
                continue
            cds_id = cur.lastrowid

            # Gene table
            gene_name = re.findall(r"Name=(.+?);",cells[8])
            if len(gene_name) != 0:
                gene_name = gene_name[0]
                cur.execute('SELECT id FROM Gene WHERE name = ? LIMIT 1', (gene_name,))
                try:
                    gene_id = cur.fetchone()[0]
                except:
                    cur.execute('INSERT OR IGNORE INTO Gene (name) VALUES (?)', (gene_name,))
                    conn.commit()
                    if cur.rowcount != 1 :
                        print('Error inserting gene: ',gene_name)
                        continue
                    gene_id = cur.lastrowid
                cur.execute('UPDATE Cds SET gene_id=? WHERE id=?', (gene_id, cds_id))
                conn.commit()
            
            # Product table
            product_name = re.findall(r"product=(.+?)$",cells[8])
            if len(product_name) != 0:
                product_name = product_name[0]
                cur.execute('SELECT id FROM Product WHERE name = ? LIMIT 1', (product_name,))
                try:
                    product_id = cur.fetchone()[0]
                except:
                    cur.execute('INSERT OR IGNORE INTO Product (name) VALUES (?)', (product_name,))
                    conn.commit()
                    if cur.rowcount != 1 :
                        print('Error inserting product: ',product_name)
                        continue
                    product_id = cur.lastrowid
                cur.execute('UPDATE Cds SET product_id=? WHERE id=?', (product_id, cds_id))
                conn.commit()


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
    for sample in report['samples']:
        sample_id = sample['id']
        sample_dir = os.path.join(temp_dir, sample_id)
        if not os.path.exists(sample_dir):
            os.makedirs(sample_dir)
        # parse gff file
        ggf_file = os.path.join(sample['annotation'], sample_id + '.gff.gz')
        fna_file = os.path.join(sample_dir, sample_id + '.fna')
        bed_file = os.path.join(sample_dir, sample_id + '.bed')
        parse_gff_file(ggf_file, bed_file, fna_file, sample_id)
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
                dna = str(seq_record.seq)
                seq_record.seq = seq_record.seq.translate(table=11)
                protein = str(seq_record.seq)
                cur.execute('UPDATE Cds SET dna =?, protein=? WHERE name=?', (dna, protein, headers[0]))
                SeqIO.write(seq_record, faa_hd, 'fasta')
        cur.commit()
        sample['bed'] = bed_file
        sample['extracted_fna_file'] = extracted_fna_file
        sample['faa_file'] = faa_file
    
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

    cd_hit_cluster_fasta = os.path.join(temp_dir, 'cluster')
    cd_hit_cluster_file = os.path.join(temp_dir, 'cluster.clstr')

    excluded_cluster = []

    greater_than_or_equal = True
    number_of_samples = len(report['samples'])
    
    for percent_match in [1, 0.99, 0.985, 0.98, 0.98]:
        cmd = f'cd-hit -i {combined_faa_file} -o {cd_hit_cluster_fasta} -s {percent_match} -c {percent_match} -T {threads} -M 0 -g 1 -d 256 > /dev/null'
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
        with open(cluster_filtered_faa_file, 'w') as fh:
            for seq_record in SeqIO.parse(combined_faa_file, "fasta"):
                if seq_record.id not in full_cluster_gene_names:
                    SeqIO.write(seq_record, fh, 'fasta')
        shutil.move(cluster_filtered_faa_file, combined_faa_file)

    # insert into database
    for cluster in excluded_cluster:
        cur.execute('INSERT INTO Cluster DEFAULT VALUES')
        conn.commit()
        cluster_id = cur.lastrowid
        for cds_name in cluster:
            cur.execute('UPDATE Cds SET cluster_id=? WHERE name=?', (cluster_id, cds_name))
    
    for repre_cds_name in clusters:
        cur.execute('INSERT INTO Cluster DEFAULT VALUES')
        conn.commit()
        cluster_id = cur.lastrowid
        cur.execute('UPDATE Cds SET cluster_id=? WHERE name=?', (cluster_id, repre_cds_name))
        for cds_name in clusters[repre_cds_name]:
            cur.execute('UPDATE Cds SET cluster_id=? WHERE name=?', (cluster_id, cds_name))

    report['cd_hit_cluster_fasta'] = cd_hit_cluster_fasta
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
            cmd = f"blastp -query {chunked_file} -db {blast_db} -evalue 1E-6 -num_threads {threads} -outfmt 6 -max_target_seqs 2000 " + "| awk '{ if ($3 > 95) print $0;}' 2> /dev/null 1> " + blast_output_file
            fh.write(cmd + '\n')
    cmd = f"parallel --bar -j {threads} -a {blast_cmds_file}"
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

    mcl_file = report['uninflated_mcl_clusters']
    with open(mcl_file, 'r') as fh:
        for line in fh:
            line = line.rstrip('\n')
            cur.execute('INSERT INTO Cluster DEFAULT VALUES')
            conn.commit()
            new_cluster_id = cur.lastrowid
            cds_names = line.split('\t')
            for cds_name in cds_names:
                cur.execute('SELECT cluster_id FROM Cds WHERE name = ? LIMIT 1', (cds_name,))
                old_cluster_id = cur.fetchone()[0]
                cur.execute('UPDATE Cds SET cluster_id=? WHERE cluster_id=?', (new_cluster_id, old_cluster_id))
                cur.execute('DELETE FROM Cluster WHERE id=?', (old_cluster_id,))

    return report

def find_paralogs(cds_ids):
    samples = {}
    for (cds_id,) in cds_ids:
        cur.execute('SELECT sample_id FROM Cds WHERE id = ? LIMIT 1', (cds_id,))
        sample_id = cur.fetchone()[0]
        if sample_id not in samples:
            samples[sample_id] = []
        samples[sample_id].append(cds_id)

    # pick paralogs with the smallest number of genes
    smallest_number = 1000000 
    paralog_cds_ids = None
    for sample_id in samples:
        sample_cds_ids = samples[sample_id]
        count = len(sample_cds_ids)
        if count > 1 and count < smallest_number:
            paralog_cds_ids = sample_cds_ids
            smallest_number = count

    return paralog_cds_ids

def find_neighbour(cds_id):
    cur.execute('SELECT contig_id, position FROM Cds WHERE id = ? LIMIT 1', (cds_id,))
    (contig_id, position) = cur.fetchone()
    pos_list = range(position-5, position)
    pos_list.extend(range(position+1, position+6))
    cur.execute('SELECT id FROM Cds WHERE contig_id = ? AND position IN ({})'.format(
        ','.join(['?'] * len(pos_list)), 
        (contig_id, pos_list)))
    return cur.fetchall()


def create_orthologs(cds_ids, paralog_cds_ids, cluster_column):
    # create new cluster
    for cds_id in paralog_cds_ids:
        cur.execute('INSERT Splitcluster (paralogs) VALUES (1)')
        conn.commit()
        splitcluster_id = cur.lastrowid
        cur.execute('UPDATE Cds SET splitcluster_id=? WHERE id=?', (splitcluster_id, cds_id))
    
    # create left over cluster
    cur.execute('INSERT Splitcluster (paralogs) VALUES (1)')
    conn.commit()
    splitcluster_id = cur.lastrowid

    # add other members of the cluster to their closest match
    for cds_id in cds_ids:
        if cds_id in paralog_cds_ids:
            continue

        neighbour_of_cds = find_neighbour(cds_id)
        if len(neighbour_of_cds) == 0:
            cur.execute('UPDATE Cds SET splitcluster_id=? WHERE id=?', (splitcluster_id, cds_id))
            continue

        # find paralog gene which is closest match with g
        best_score = 0
        belong = None
        for paralog_cds_id in paralog_cds_ids:
            neighbour_of_paralog = find_neighbour(paralog_cds_id)
            cur.execute('SELECT {} FROM Cds WHERE id IN ({})'.format(
                cluster_column,
                ','.join(['?'] * len(neighbour_of_paralog)), 
                (neighbour_of_paralog,)))
            clusters_around_p = cur.fetchall()
            
            score_of_p = 0
            for (neighbour_cds_id,) in neighbour_of_cds:
                cur.execute(f'SELECT {cluster_column} FROM Cds WHERE id = ? LIMIT 1', (neighbour_cds_id,))
                a = cur.fetchone()
                if a in clusters_around_p:
                    score_of_p += 1
            score_of_p = score_of_p / len(neighbour_of_cds)
            if score_of_p > best_score:
                best_score = score_of_p
                belong = paralog_cds_id

        if belong == None:
            cur.execute('UPDATE Cds SET splitcluster_id=? WHERE id=?', (splitcluster_id, cds_id))
        else:
            cur.execute('''UPDATE Cds SET splitcluster_id=(
                SELECT splitcluster_id FROM Cds WHERE id = ?
            ) WHERE id=?''', (belong, cds_id,))

    # check for "leftovers", remove if absent
    cur.excecute('SELECT id FROM Cds WHERE splitcluster_id = ?', (splitcluster_id,))
    try:
        cur.fetchone()[0]
    except:
        cur.execute('DELETE FROM Splitcluster WHERE id=?', (splitcluster_id,))

def split_paralogs(report):
    cur.execute('SELECT id FROM Cluster WHERE paralogs=1')
    for (cluster_id,) in cur:
        cur.execute('SELECT id FROM Cds WHERE cluster_id = ?', (cluster_id,))
        cds_ids = cur.fetchall()
        if len(cds_ids) == 1:
            cds_id = cds_ids[0][0]
            cur.execute('INSERT Splitcluster DEFAULT VALUES')
            conn.commit()
            splitcluster_id = cur.lastrowid
            cur.execute('UPDATE Cds SET splitcluster_id=? WHERE id=?', (splitcluster_id, cds_id))
            cur.execute('UPDATE Cluster SET paralogs=0 WHERE id=?', (cluster_id,))
            continue
        # check paralogs
        paralog_cds_ids = find_paralogs(cds_ids)
        if paralog_cds_ids == None:
            cur.execute('INSERT Splitcluster DEFAULT VALUES')
            conn.commit()
            splitcluster_id = cur.lastrowid
            cur.execute('UPDATE Cluster SET paralogs=0 WHERE id=?', (cluster_id,))
            for (cds_id,) in cds_ids:
                cur.execute('UPDATE Cds SET splitcluster_id=? WHERE id=?', (splitcluster_id, cds_id))
            continue
        
        # split paralogs
        create_orthologs(cds_ids, paralog_cds_ids, 'cluster_id')

    # run iteratively
    for i in range(50):
        cur.execute('SELECT id FROM Splitcluster WHERE paralogs=1')
        for (old_splitcluster_id,) in cur:
            cur.execute('SELECT id FROM Cds WHERE splitcluster_id = ?', (old_splitcluster_id,))
            cds_ids = cur.fetchall()
            if len(cds_ids) == 1:
                cds_id = cds_ids[0][0]
                cur.execute('UPDATE Splitcluster SET paralogs=0 WHERE id=?', (old_splitcluster_id,))
                continue
            # check paralogs
            paralog_cds_ids = find_paralogs(cds_ids)
            if paralog_cds_ids == None:
                cur.execute('UPDATE Splitcluster SET paralogs=0 WHERE id=?', (old_splitcluster_id,))
                continue

            # delete old cluster
            cur.execute('DELETE FROM Splitcluster WHERE id=?', (old_splitcluster_id,))
            
            # split paralogs
            create_orthologs(cds_ids, paralog_cds_ids, 'splitcluster_id')


def annotate_cluster(report):
    """
    Update the cluster name to the gene name

    Parameters
    -------
    -------
    """
    cur.execute('SELECT id FROM Splitcluster')
    count = 1
    for (splitcluster_id,) in cur:
        # update cluster name
        cur.execute('''
        SELECT 
            gene_id,
            COUNT(gene_id) AS num_gene 
        FROM Cds 
        WHERE splitcluster_id=? 
        GROUP BY gene_id
        ORDER BY num_gene DESC
        LIMIT 1
        ''',(splitcluster_id,))
        try:
            gene_id = cur.fetchone()[0]
            cur.execute('''
            UPDATE Splitcluster
            SET name=(SELECT name FROM Gene WHERE id=?)
            WHERE id=?  
            ''', (gene_id, splitcluster_id))
        except:
            cluster_name = 'group_' + str(count)
            count += 1
            cur.execute('''
            UPDATE Splitcluster
            SET name=?
            WHERE id=?  
            ''', (cluster_name, splitcluster_id))

        # update cluster product
        cur.execute('''
        SELECT 
            product_id,
            COUNT(product_id) AS num_product 
        FROM Cds 
        WHERE splitcluster_id=? 
        GROUP BY product_id
        ORDER BY num_product DESC
        LIMIT 1
        ''',(splitcluster_id,))
        try:
            product_id = cur.fetchone()[0]
            cur.execute('''
            UPDATE Splitcluster
            SET product_id=?
            WHERE id=?  
            ''', (product_id, splitcluster_id))
        except:
            pass


def create_spreadsheet(report):
    spreadsheet_file = os.path.join(report['pan_genome'], 'gene_presence_absence.csv.gz')
    with gzip.open(spreadsheet_file, 'wt') as fh:
        writer = csv.writer(fh, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)

        # write header
        header = ['Gene', 'Non-unique Gene name', 'Annotation', 'No. isolates', 'No. sequences', 'Avg sequences per isolate', 'Genome Fragment','Order within Fragment', 'Accessory Fragment','Accessory Order with Fragment', 'QC','Min group size nuc', 'Max group size nuc', 'Avg group size nuc' ]
        cur.execute('SELECT name FROM Sample')
        for (sample_name,) in cur:
            header.append(sample_name)
        writer.writerow(header)

        # write row
        cur.execute('SELECT id, name, product FROM Splitcluster')
        for (splitcluster_id, splitcluster_name, splitcluster_product) in cur:
            row = []
            # Gene
            row.append(splitcluster_name)
            # Non-unique Gene name
            row.append("")
            # Annotation
            if splitcluster_product = None:
                row.append("")
            else:
                row.append(splitcluster_product)
            # No. isolates
            cur.execute('SELECT sample_id FROM Cds WHERE splitcluster_id=? GROUP BY sample_id',(splitcluster_id,))
            row.append(len(cur.fetchone()))
            # No. sequences
            cur.execute('SELECT count(id) FROM Cds WHERE splitcluster_id=?',(splitcluster_id,))
            row.append(cur.fetchone()[0])
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
            cur.execute('SELECT id, name FROM Sample')
            for (sample_id, sample_name) in cur:
                cur.execute('SELECT name FROM Cds WHERE splitcluster_id=? AND sample_id=?', (splitcluster_id, sample_id))
                result = cur.fetchall()
                if len(result) == 0:
                     row.append('')
                else:
                    cds_names = []
                    for (cds_name,) in result:
                        cds_names.append(cds_name)
                    row.append('\t'.join(gene_list))
            writer.writerow(row)
    report['spreadsheet'] = spreadsheet_file
    return report


def create_rtab(report):
    rtab_file = os.path.join(report['pan_genome'], 'gene_presence_absence.Rtab')
    with open(rtab_file, 'w') as fh:
        writer = csv.writer(fh, delimiter='\t')

        # write header
        header = ['Gene']
        cur.execute('SELECT name FROM Sample')
        for (sample_name,) in cur:
            header.append(sample_name)
        writer.writerow(header)

        # write row
        cur.execute('SELECT id, name FROM Splitcluster')
        for (splitcluster_id, splitcluster_name) in cur:
            row = []
            # Gene
            row.append(splitcluster_name)
            # Samples
            cur.execute('SELECT id, name FROM Sample')
            for (sample_id, sample_name) in cur:
                cur.execute('SELECT id FROM Cds WHERE splitcluster_id=? AND sample_id=?', (splitcluster_id, sample_id))
                row.append(len(cur.fetchall()))
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
        report['summary'] = os.path.join(pan_genome_folder,'summary_statistics.txt')
        report['rtab'] = os.path.join(report['pan_genome'], 'gene_presence_absence.Rtab')
        report['spreadsheet'] = os.path.join(report['pan_genome'], 'gene_presence_absence.csv.gz')
        return report

    if not os.path.exists(pan_genome_folder):
        os.makedirs(pan_genome_folder)
    if not os.path.exists(temp_dir):
        os.makedirs(temp_dir)

    # create database
    conn = sqlite3.connect(os.path.join(pan_genome_folder, 'pan_genome.sqlite'))
    cur = conn.cursor()
    cur.executescript('''
        CREATE TABLE IF NOT EXISTS Cds (
            id              INTEGER NOT NULL PRIMARY KEY AUTOINCREMENT UNIQUE,
            name            TEXT UNIQUE,
            sample_id       INTEGER,
            contig_id       INTEGER,
            position        INTEGER,
            trand           TEXT,
            start           INTEGER,
            end             INTEGER,
            length          INTEGER,
            gene_id         INTEGER,
            product_id      INTEGER,
            dna             TEXT,
            protein         TEXT,
            cluster_id      INTEGER,
            splitcluster_id      INTEGER
        );
        CREATE TABLE IF NOT EXISTS Cluster (
            id              INTEGER NOT NULL PRIMARY KEY AUTOINCREMENT UNIQUE,
            name            TEXT UNIQUE,
            product_id      INTEGER,
            paralogs        INTEGER DEFAULT 1
        );
        CREATE TABLE IF NOT EXISTS Spitcluster (
            id              INTEGER NOT NULL PRIMARY KEY AUTOINCREMENT UNIQUE,
            name            TEXT UNIQUE,
            product_id      INTEGER,
            paralogs        INTEGER DEFAULT 0
        );        
        CREATE TABLE IF NOT EXISTS Gene (
            id              INTEGER NOT NULL PRIMARY KEY AUTOINCREMENT UNIQUE,
            name            TEXT UNIQUE
        );
        CREATE TABLE IF NOT EXISTS Sample (
            id              INTEGER NOT NULL PRIMARY KEY AUTOINCREMENT UNIQUE,
            name            TEXT UNIQUE
        );
        CREATE TABLE IF NOT EXISTS Contig (
            id              INTEGER NOT NULL PRIMARY KEY AUTOINCREMENT UNIQUE,
            name            TEXT UNIQUE
        );
        CREATE TABLE IF NOT EXISTS Product (
            id              INTEGER NOT NULL PRIMARY KEY AUTOINCREMENT UNIQUE,
            name            TEXT UNIQUE
        )
    ''')
    conn.commit()

    for sample in report['samples']:
        sample_name = sample['id']
        cur.execute('INSERT OR IGNORE INTO Sample (name) VALUES (?)', (sample_name,))
    conn.commit()

    report = extract_proteins(report, timing_log=timing_log)
    report = combine_proteins(report, timing_log=timing_log)
    report = run_cd_hit_iterative(report, threads=threads, timing_log=timing_log)
    report = all_against_all_blast(report, threads=threads, timing_log=timing_log)
    report = cluster_with_mcl(report, threads=threads, timing_log=timing_log)
    report = reinflate_clusters(report)
    report = split_paralogs(report)
    report = annotate_cluster(report)
    report = create_spreadsheet(report)
    report = create_rtab(report)
    report = create_summary(report)

    # clean
    del report['blast_result_file']
    del report['cd_hit_cluster_fasta']
    del report['combined_faa_file']
    del report['uninflated_mcl_clusters']
    for sample in report['samples']:
        del sample['bed']
        del sample['extracted_fna_file']
        del sample['faa_file']
    #shutil.rmtree(temp_dir)

    return report

    
if __name__ == '__main__':
    report = json.load(open('dev/report', 'r'))
    report = run_pan_genome_analysis(
        report, 
        collection_dir='dev', 
        threads=8, 
        overwrite=True, 
        timing_log='dev/time.log'
    )
    json.dump(report, open('dev/report_output.json', 'w'), indent=4, sort_keys=True)

## TODO ##
# filter protein sequence ?
# why run cd-hit iteratively that way
# improve split paralogs