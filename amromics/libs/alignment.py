import os,sys,shutil
import logging
import multiprocessing
import gzip
from Bio import SeqIO
import csv
import pandas as pd
import json
import re
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import amromics.libs.bioseq as bioseq
import amromics.libs.msa2vcf as msa2vcf
from amromics.libs.bioseq import read_sequence_file, write_fasta
from amromics.utils.command import run_command
from datetime import datetime
logger = logging.getLogger(__name__)
NUM_CORES_DEFAULT = multiprocessing.cpu_count()
def run_alignment_by_parsnp(pan_folder,ffn_dir,base_dir, overwrite=False,  timing_log=None,threads=0):
    """
        Run aligment process to create both multi-alignment  and phylogeny tree for each gene in gene clusters
        :param read_data: result holder
        :param ffn_dir: path to folder of .ffn (output of prokka)
        :param base_dir: working directory
        :return:
    """
    gene_cluster_file=pan_folder+'/gene_presence_absence.csv'
    dict_cds={}
    for root, dirs, files in os.walk(ffn_dir):
        for _file in files:
            if _file.endswith('.ffn.gz'):
                with gzip.open(os.path.join(root, _file), 'rt') as fn:
                    for seq in SeqIO.parse(fn, 'fasta'):
                        dict_cds[seq.id] = seq


    #make folder contains sequences for each gene
    alignment_dir=os.path.join(base_dir,'alignments')
    if (not overwrite) and os.path.exists(alignment_dir):
        return alignment_dir
    if not os.path.exists(alignment_dir):
        os.makedirs(alignment_dir)

    gene_df = pd.read_csv(gene_cluster_file, dtype=str)
    gene_df.fillna('', inplace=True)

    sample_columns = list(gene_df.columns)[14:]
    for _, row in gene_df.iterrows():
        gene_id = row['Gene']
        gene_list = []
        for sample_column in sample_columns:
            if row[sample_column]:
                # roary can pool together genes from the same sample and tab-separate them
                for sample_gene in row[sample_column].split('\t'):
                    gene_list.append(sample_gene)
                    # TODO: make sure all samples in this gene have not updated

        gene_list = sorted(gene_list)
        # Only analyse if there are more than 3 genes
        if len(gene_list) < 3:
            logger.info('There are too few genes for {} skipping'.format(gene_id))
            continue

        gene_dir = os.path.join(alignment_dir, gene_id)
        # Check if done before
        gene_list_json = os.path.join(gene_dir, 'gene_list.json')
        # if os.path.isfile(os.path.join(gene_dir, 'parsnp.tree')) and (not overwrite):
        if os.path.isfile(gene_list_json):
            with open(gene_list_json) as fn:
                existing_gene_list = json.load(fn)
                if gene_list == existing_gene_list:
                    logger.info('Phylogeny for gene {} done, skipping'.format(gene_id))
                    continue  # for _, row

        gene_file_dir = os.path.join(gene_dir, 'files')
        if not os.path.exists(gene_file_dir):
            os.makedirs(gene_file_dir)

        gene_files = []
        for sample_gene in gene_list:
            gene_file = os.path.join(gene_file_dir, sample_gene + '.fasta')
            SeqIO.write(dict_cds[sample_gene], gene_file, 'fasta')
            gene_files.append(gene_file)

        # Use the first gene as the reference
        cmd = 'parsnp -d {} -r {} -o {} -p {}'.format(
            ' '.join(gene_files[1:]), gene_files[0], gene_dir, threads)
        ret = run_command(cmd)
        # if ret != 0:
        #     raise Exception('error')

        with open(gene_list_json, 'w') as fn:
            json.dump(gene_list, fn)
        #run_command('gzip {}'.format(os.path.join(gene_dir, 'parsnp.xmfa')))
        #run_command('gzip {}'.format(os.path.join(gene_dir, 'parsnp.ggr')))

        if os.path.exists(gene_file_dir):
            shutil.rmtree(gene_file_dir)
        #clean up
        run_command('rm -f ' + os.path.join(gene_dir, '*.ini ') + os.path.join(gene_dir, '*block* '))
        shutil.rmtree(os.path.join(gene_dir, 'blocks'), True)
        shutil.rmtree(os.path.join(gene_dir, 'tmp'), True)


    return alignment_dir
def run_protein_alignment(pan_folder, collection_dir, threads=8, overwrite=False, timing_log=None, rate_coverage=0.15):
    """
    Align protein sequence by mafft

    Parameters
    ----------
    report: object
        A report object
    collection_dir: str
        working directory of the collection
    threads: int
        number of threads to use
    overwrite: bool
        whether to overwrite existing result even if input did not change
    timing_log: str
        file to log timing
    min_cover: float
        Only genes present in at least min_cover%  of samples will run alignment
    Returns
        report object
    -------
    """
    alignment_dir = os.path.join(collection_dir, 'alignments')

    gene_cluster_file =pan_folder + '/gene_presence_absence.Rtab'
    gene_df = pd.read_csv(gene_cluster_file, sep='\t', index_col='Gene')
    gene_df.fillna('', inplace=True)

    min_number=int(rate_coverage*len(gene_df.columns))
    cmds_file = os.path.join(alignment_dir,"align_cmds")
    with open(cmds_file,'w') as cmds:
        for gene_id, row in gene_df.iterrows():
            # Only align if there are at least 2 sequences
            if row.sum() < 2:
                continue
            # Only align if the number  of sequences > min_number
            if row.sum() < min_number:
                continue
            gene_id = re.sub(r'\W+', '', gene_id)
            gene_dir = os.path.join(alignment_dir, gene_id)

            # check if done before
            gene_aln_file = os.path.join(gene_dir, gene_id + '.faa.aln')
            if (not overwrite) and os.path.isfile(gene_aln_file):
                continue

            gene_seq_file = os.path.join(gene_dir, gene_id + '.faa')
            if not os.path.isfile(gene_seq_file):
                #TODO logger.info('{} does not exist'.format(gene_aln_file))
                continue

            cmd = f"mafft --auto --quiet --thread 1 {gene_seq_file} > {gene_aln_file}"
            cmds.write(cmd + '\n')

    cmd = f"parallel --bar -j {threads} -a {cmds_file}"
    ret = run_command(cmd, timing_log)
    #report['alignments'] = alignment_dir
    return alignment_dir


def create_nucleotide_alignment(pan_folder, collection_dir, threads=8, overwrite=False, timing_log=None,min_cover=0.15):
    """
    Create nucleotide alignment according to protein alignment

    Parameters
    ----------
    report: object
        A report object
    collection_dir: str
        working directory of the collection
    threads: int
        number of threads to use
    overwrite: bool
        whether to overwrite existing result even if input did not change
    timing_log: str
        file to log timing
    Returns
        report object
    -------
    """
    logger.info('Creating nucleotide alignment')
    alignment_dir = os.path.join(collection_dir, 'alignments')

    gene_cluster_file = pan_folder + '/gene_presence_absence.Rtab'
    gene_df = pd.read_csv(gene_cluster_file, sep='\t', index_col='Gene')
    gene_df.fillna('', inplace=True)
    min_number=int(min_cover*len(gene_df.columns))
    for gene_id, row in gene_df.iterrows():
        # Only run if there are at least 2 sequences
        if row.sum() < 2:
            continue
        # Only align if the number  of sequences > min_number
        if row.sum() < min_number:
            continue
        gene_id = re.sub(r'\W+', '', gene_id)
        gene_dir = os.path.join(alignment_dir, gene_id)

        # check if done before
        nucleotide_aln_file = os.path.join(gene_dir, gene_id + '.fna.aln')
        if (not overwrite) and os.path.isfile(nucleotide_aln_file):
            continue

        protein_aln_file = os.path.join(gene_dir, gene_id + '.faa.aln')
        if not os.path.isfile(protein_aln_file):
            #TODO: logger.info('{} does not exist'.format(protein_aln_file))
            continue
        protein_dict = {}
        with open(protein_aln_file) as fh:
            for seq_record in SeqIO.parse(fh, 'fasta'):
                protein_dict[seq_record.id] = str(seq_record.seq)

        nucleotide_seq_file = os.path.join(gene_dir, gene_id + '.fna')
        nucleotide_dict = {}
        for seq_record in SeqIO.parse(nucleotide_seq_file, 'fasta'):
            nucleotide_dict[seq_record.id] = str(seq_record.seq)

        with open(nucleotide_aln_file, 'wt') as fh:
            for seq_id in protein_dict.keys():
                protein = protein_dict[seq_id]
                nucleotide = nucleotide_dict[seq_id]
                result = ''
                codon_pos = 0
                for c in protein:
                    if c == '-':
                        result += '---'
                    else:
                        result += nucleotide[codon_pos * 3: codon_pos * 3 + 3]
                        codon_pos += 1
                new_record = SeqRecord(Seq(result), id = seq_id, description = '')
                SeqIO.write(new_record, fh, 'fasta')

        # remove input files
        os.remove(nucleotide_seq_file)
        protein_seq_file = os.path.join(gene_dir, gene_id + '.faa')
        os.remove(protein_seq_file)

    return alignment_dir
def create_core_gene_alignment(pan_folder, collection_dir, threads=8, overwrite=False, timing_log=None):
    """
    Concatenate all the nucleotide alignment of core genes to create core gene alignment

    Parameters
    ----------
    report: object
        A report object
    collection_dir: str
        working directory of the collection
    threads: int
        number of threads to use
    overwrite: bool
        whether to overwrite existing result even if input did not change
    timing_log: str
        file to log timing
    Returns
        report object
    -------
    """
    logger.info('Creating core gene alignment')
    alignment_dir = os.path.join(collection_dir, 'alignments')
    phylogeny_folder = os.path.join(collection_dir, 'phylogeny')
    #report['phylogeny'] = phylogeny_folder
    if not os.path.exists(phylogeny_folder):
        os.makedirs(phylogeny_folder)
    # check if done before
    core_gene_aln_file = os.path.join(phylogeny_folder, 'core_gene_alignment.aln')
    if os.path.isfile(core_gene_aln_file) and (not overwrite):
        logger.info('Core gene alignment exists and input has not changed, skipping')
        return phylogeny_folder

    gene_cluster_file = pan_folder + '/gene_presence_absence.Rtab'
    gene_df = pd.read_csv(gene_cluster_file, sep='\t', index_col='Gene')
    gene_df.fillna('', inplace=True)

    seq_dict = {}
    samples=gene_df.columns.tolist()
    for sample in samples:
        seq_dict[sample]= ''
    sample_list = seq_dict.keys()
    for gene_id, row in gene_df.iterrows():
        # Only run if it is core gene
        if len(row[row == 0]) != 0:
            continue
        gene_id = re.sub(r'\W+', '', gene_id)
        gene_dir = os.path.join(alignment_dir, gene_id)

        nucleotide_aln_file = os.path.join(gene_dir, gene_id + '.fna.aln')
        if not os.path.isfile(nucleotide_aln_file):
            #TODO logger.info('{} does not exist'.format(nucleotide_aln_file))
            continue
        cluster_dict = {}
        with open(nucleotide_aln_file, 'rt') as fh:
            for seq_record in SeqIO.parse(fh, 'fasta'):
                if '-' in seq_record.id:
                    gene_sample_id=seq_record.id[seq_record.id.rfind('-')+1:]
                else:
                    gene_sample_id=seq_record.id
                sample_name = re.findall(r'^(.+)_', gene_sample_id)
                sample_name = sample_name[0]
                if sample_name not in sample_list:
                    print(sample_list)
                    raise Exception(f'Error concatenating gene alignment: {sample_name} is not a sample id')
                cluster_dict[sample_name] = str(seq_record.seq)

        for sample_name in cluster_dict:
            seq_dict[sample_name] += cluster_dict[sample_name]

    with open(core_gene_aln_file, 'w') as fh:
        for sample in sample_list:
            new_record = SeqRecord(Seq(seq_dict[sample]), id = sample, description = '')
            SeqIO.write(new_record, fh, 'fasta')

    return core_gene_aln_file
def writeTempSeqFile(temp_seqs_dir, gene_id,seq):

    if not os.path.exists(temp_seqs_dir):
        os.makedirs(temp_seqs_dir)
    temp_file=os.path.join(temp_seqs_dir,gene_id+".fasta")
    f=open(temp_file,"w")
    f.write(">"+str(gene_id)+"\n")
    f.write(seq)
    f.close()
    return temp_file
def get_gene_sequences(samples, pan_folder,sample_col, collection_dir):
    """
    Create protein sequences and nucleotide sequences for each gene cluster

    Parameters
    ----------
    pan_folder: string
        pan_folder
    collection_dir: str
        working directory of the collection
    threads: int
        number of threads to use
    overwrite: bool
        whether to overwrite existing result even if input did not change
    timing_log: str
        file to log timing
    Returns
        report object
    -------
    """
    logger.info('Getting sequences of gene clusters')
    gene_cluster_file = pan_folder + '/gene_presence_absence.csv'
    dict_nucleotide = {}

    ffn_files = [sample['annotation_ffn'] for sample in samples]
    for ffn_file_path in ffn_files:
        ffn_binfile_path=ffn_file_path.replace('.ffn','.bin')
        if os.path.exists(ffn_binfile_path):
            os.remove(ffn_binfile_path)
        with open(ffn_binfile_path, 'wb') as binary_fh:
            with open(ffn_file_path, 'r') as fasta_fh:
                for seq_record in SeqIO.parse(fasta_fh, 'fasta'):
                    start_position = binary_fh.tell()
                    binary_fh.write(str(seq_record.seq).encode('utf-8'))
                    end_position = binary_fh.tell()
                    #seq_record.seq = seq_record.seq
                    #seq_record = SeqRecord(seq_record.seq, id=seq_record.id, description = '')
                    #dict_nucleotide[seq_record.id]=writeTempSeqFile(temp_dna_seq_folder,seq_record.id,str(seq_record.seq))
                    dict_nucleotide[seq_record.id] =(ffn_binfile_path,start_position, end_position)
                    #dict_nucleotide[seq_record.id] = seq_record
    dict_prot = {}
    faa_files = [sample['annotation_faa'] for sample in samples]
    for faa_file_path in faa_files:
        faa_binfile_path=faa_file_path.replace('.faa','.bin')
        if os.path.exists(faa_binfile_path):
            os.remove(faa_binfile_path)
        with open(faa_binfile_path, 'wb') as binary_fh:
            with open(faa_file_path, 'r') as fasta_fh:
                for seq_record in SeqIO.parse(fasta_fh, 'fasta'):
                    start_position = binary_fh.tell()
                    binary_fh.write(str(seq_record.seq).encode('utf-8'))
                    end_position = binary_fh.tell()
                    #seq_record.seq = seq_record.seq
                    #seq_record = SeqRecord(seq_record.seq, id=seq_record.id, description = '')
                    #dict_prot[seq_record.id] = writeTempSeqFile(temp_prot_seq_folder,seq_record.id,str(seq_record.seq))
                    #dict_prot[seq_record.id] = seq_record
                    dict_prot[seq_record.id] = (faa_binfile_path,start_position, end_position)
    #print (dict_nucleotide)
    # make folder contains sequences for each gene
    alignment_dir = os.path.join(collection_dir, 'alignments')
    if not os.path.exists(alignment_dir):
        os.makedirs(alignment_dir)
    #report['alignments'] = alignment_dir

    gene_df = pd.read_csv(gene_cluster_file, dtype=str)
    gene_df.fillna('', inplace=True)
    sample_columns = list(gene_df.columns)[sample_col:]
    for _, row in gene_df.iterrows():
        gene_id = row['Gene']
        gene_id = re.sub(r'\W+', '', gene_id)
        gene_dir = os.path.join(alignment_dir, gene_id)
        if not os.path.exists(gene_dir):
            os.makedirs(gene_dir)

        # check if done before
        protein_seq_file = os.path.join(gene_dir, gene_id + '.faa')
        nucleotide_seq_file = os.path.join(gene_dir, gene_id + '.fna')

        gene_list = []
        for sample_column in sample_columns:
            if row[sample_column]:
                # roary can pool together genes from the same sample and tab-separate them
                for sample_gene in row[sample_column].split('\t'):
                    if '-' in sample_gene:
                        sample_gene=sample_gene.split('-')[-1]
                    gene_list.append(sample_gene)
        gene_list = sorted(gene_list)
        #print(gene_list)
        with open(protein_seq_file, 'w') as prot_fh, open(nucleotide_seq_file, 'w') as nucl_fh:
            for sample_gene in gene_list:
                #nu_seq_record = dict_nucleotide[sample_gene]
                #pro_seq_record=dict_prot[sample_gene]
                #pro_seq = nu_seq_record.seq.translate(table=11)
                #pro_seq, isReversed = translateDNA2Prot(nu_seq_record)
                #pro_seq_record = SeqRecord(pro_seq, id = nu_seq_record.id, description = '')
                #if isReversed:
                #    nu_seq_record.seq=nu_seq_record.seq.reverse_complement()
                #for record in SeqIO.parse(nu_seq_record,"fasta"):
                #    SeqIO.write(record, nucl_fh, 'fasta')
                #for record in SeqIO.parse(pro_seq_record,"fasta"):
                #    SeqIO.write(record, prot_fh, 'fasta')
                #record_nu_dict = SeqIO.index(nu_seq_record, "fasta")
                #SeqIO.write(record_nu_dict[sample_gene], nucl_fh, 'fasta')

                #record_prot_dict = SeqIO.index(pro_seq_record, "fasta")
                #SeqIO.write(record_prot_dict[sample_gene], prot_fh, 'fasta')
                # for record in SeqIO.parse(nu_seq_record,"fasta"):
                #     if record.id==sample_gene:
                #         SeqIO.write(record, nucl_fh, 'fasta')
                #         break
                # for record in SeqIO.parse(pro_seq_record,"fasta"):
                #     if record.id==sample_gene:
                #         SeqIO.write(record, prot_fh, 'fasta')
                #         break
                nu_seq_file,start_nu_seq,end_nu_seq = dict_nucleotide[sample_gene]
                nu_record=SeqRecord(Seq(binary_to_record(nu_seq_file,start_nu_seq,end_nu_seq)), id = sample_gene, description = '')
                #print(nu_record.seq)
                SeqIO.write(nu_record, nucl_fh, 'fasta')
                pro_seq_file,start_pro_seq,end_pro_seq = dict_prot[sample_gene]
                pro_record=SeqRecord(Seq(binary_to_record(pro_seq_file,start_pro_seq,end_pro_seq)), id = sample_gene, description = '')
                SeqIO.write(pro_record, prot_fh, 'fasta')
    # if os.path.exists(temp_dna_seq_folder):
    #     shutil.rmtree(temp_dna_seq_folder)
    # if os.path.exists(temp_prot_seq_folder):
    #     shutil.rmtree(temp_prot_seq_folder)
    return alignment_dir
def binary_to_record(binary_file, start_position, end_position):
    record=None
    with open(binary_file, 'rb') as binary_fh:

        binary_fh.seek(start_position)


        data = binary_fh.read(end_position - start_position)


        record =data.decode('utf-8')


    return record
def translateDNA2Prot(sr):
    #print(type(sr))
    prot_d1=sr.seq.translate(table=11)
    count_d1=prot_d1.count('*')
    prot_d2=sr.seq.reverse_complement().translate(table=11)
    count_d2=prot_d2.count('*')
    if prot_d2.endswith('*'):
        #print("seq.reverse_complement() called, d1="+str(count_d1)+",count_d2="+str(count_d2))
        return prot_d2, True
    else:
        return prot_d1, False

def runGeneAlignment(samples, pan_folder,sample_col,collection_dir, threads=8, overwrite=False, timing_log=None,rate_coverage=0.15):
    stime = datetime.now()
    alignment_dir=get_gene_sequences(samples, pan_folder,sample_col,collection_dir=collection_dir)
    elapsed = datetime.now() - stime
    logger.info(f'Get gene sequences -- time taken {str(elapsed)}')
    stime = datetime.now()
    alignment_dir=run_protein_alignment(pan_folder, collection_dir=collection_dir, overwrite=overwrite,threads=threads,timing_log=timing_log,rate_coverage=rate_coverage)
    elapsed = datetime.now() - stime
    logger.info(f'Protein alignment -- time taken {str(elapsed)}')
    stime = datetime.now()
    alignment_dir=create_nucleotide_alignment(pan_folder, collection_dir=collection_dir, overwrite=overwrite,threads=threads,timing_log=timing_log)
    elapsed = datetime.now() - stime
    logger.info(f'Nucliotide alignment -- time taken {str(elapsed)}')
    return alignment_dir
def copyClustersPanta(collection_dir,panta_folder):
    alignment_dir = os.path.join(collection_dir, 'alignments')
    if os.path.exists(alignment_dir):
        shutil.rmtree(alignment_dir)
    panta_cluster=os.path.join(panta_folder,"clusters")
    if not os.path.exists(panta_cluster):
        return None
    shutil.copytree(panta_cluster, alignment_dir)
    return alignment_dir
def runVCFCallingFromGeneAlignment(samples,pangenome_folder, collection_dir, threads=8, overwrite=False, timing_log=None):
    """
    Call VCFs files between representative seq and other seq

    Parameters
    ----------
    pangenome folder: path
        Path to panta output
    collection_dir: str
        working directory of the collection
    threads: int
        number of threads to use
    overwrite: bool
        whether to overwrite existing result even if input did not change
    timing_log: str
        file to log timing
    Returns
        alignment folder
    -------
    """
    alignment_dir = os.path.join(collection_dir, 'alignments')
    vcf_dir = os.path.join(collection_dir, 'VCFs')
    try:

        if not os.path.exists(vcf_dir):
            os.makedirs(vcf_dir)
        gene_cluster_file =pangenome_folder + '/gene_presence_absence.Rtab'
        gene_df = pd.read_csv(gene_cluster_file, sep='\t', index_col='Gene')
        gene_df.fillna('', inplace=True)
        map_sample_vcf={}
        ref_pan=[]
        map_sample_prot_vcf={}
        for gene_id, row in gene_df.iterrows():
            # Only align if there are at least 2 sequences
            if row.sum() < 2:
                continue
            gene_id = re.sub(r'\W+', '', gene_id)
            gene_dir = os.path.join(alignment_dir, gene_id)
            # check if done before

            gene_aln_file = os.path.join(gene_dir, gene_id + '.fna.aln')
            if not os.path.isfile(gene_aln_file):
                #continue if gene is not aligned
                continue

            gene_prot_aln_file = os.path.join(gene_dir, gene_id + '.faa.aln')
            if not os.path.isfile(gene_prot_aln_file):
                #continue if gene is not aligned
                continue
            rep_name=None
            rep_seq=None
            for sample in samples:
                isFound=False
                for seq in SeqIO.parse(gene_aln_file, 'fasta'):

                    if '-' in seq.name:
                        gene_sample_id= seq.name[seq.name.rfind('-')+1:]
                    else:
                        gene_sample_id=seq.name
                    s=gene_sample_id[:gene_sample_id.rfind('_')]
                    if sample['id']==s:
                        rep_seq=seq
                        isFound=True
                        break
                if isFound:
                    break
            rep_seq.seq=rep_seq.seq.replace('-','')
            rep_name=rep_seq.name
            ref_pan.append(rep_seq)
            # for seq in read_sequence_file(gene_aln_file):
            #     rep_name=seq.name
            #     ref_pan.append(seq)
            #     break
            if rep_name==None:
                continue
            map_gene_vcf=msa2vcf.go(gene_aln_file,rep_name,gene_dir)
            map_sample_vcf=createSampleVcfDict(map_gene_vcf,map_sample_vcf,vcf_dir)

            map_gene_prot_vcf=msa2vcf.go(gene_prot_aln_file,rep_name,gene_dir,isProt=True)
            map_sample_prot_vcf=createSampleVcfDict(map_gene_prot_vcf,map_sample_prot_vcf,vcf_dir)
        #print(map_sample_vcf)
        #print("size of dict map sample and vcf : "+str(sys.getsizeof(map_sample_vcf)))
        for s in map_sample_vcf.keys():
            #print("generate vcf for sample "+s)
            vcf_sample_dir=os.path.join(vcf_dir,s)
            if not os.path.exists(vcf_sample_dir):
                continue
            vcf_file = os.path.join(vcf_sample_dir,s+".vcf")

            generateSampleVcfFile(map_sample_vcf[s],vcf_file,s)
        for s in map_sample_prot_vcf.keys():
            vcf_sample_dir=os.path.join(vcf_dir,s)
            if not os.path.exists(vcf_sample_dir):
                continue
            vcf_file = os.path.join(vcf_sample_dir,s+".prot.vcf")

            generateSampleVcfFile(map_sample_prot_vcf[s],vcf_file,s)
        print("write "+str(len(ref_pan))+" sequences to "+os.path.join(vcf_dir,"pangenome_reference.fasta"))
        with open(os.path.join(vcf_dir,"pangenome_reference.fasta"),'w') as o:
            SeqIO.write(ref_pan,o,"fasta")
        #write_fasta(os.path.join(vcf_dir,"pangenome_reference.fasta"),ref_pan)
    except Exception as error:
        logger.error("Error create vcf files:"+str(error))

    #cmd = f"parallel --bar -j {threads} -a {cmds_file}"
    #ret = run_command(cmd, timing_log)
    #report['alignments'] = alignment_dir
    return vcf_dir
def createSampleVcfDict(map_vcf,dict,vcf_dir):
    for g in map_vcf.keys():
        if '-' in g:
            gene_sample_id= g[g.rfind('-')+1:]
        else:
            gene_sample_id=g
        s=gene_sample_id[:gene_sample_id.rfind('_')]
        id=gene_sample_id[gene_sample_id.rfind('_')+1:]
        vcf_sample_dir=os.path.join(vcf_dir,s)
        if not os.path.exists(vcf_sample_dir):
            os.makedirs(vcf_sample_dir)
        # with open(vcf_sample_dir+"/"+id +".vcf", 'w') as f:
        #     for line in map_gene_vcf[g]:
        #         f.write(line+"\n")

        if not s in dict.keys():
             dict[s]=[]
        #map_sample_vcf[s].append(os.path.join(vcf_sample_dir,id+".vcf"))
        vcf_obj={"gene":g,"vcf":map_vcf[g]}
        dict[s].append(vcf_obj)
    return dict
def generateSampleVcfFile(list_vcf,vcf_file,sample_name):
    with open(vcf_file,'w') as vcf:
        vcf.write("##fileformat=VCFv4.2\n")
        vcf.write("##source="+os.path.basename(sys.argv[0])+"\n")
        vcf.write("##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">\n")
        vcf.write("##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">\n")
        vcf.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
        #vcf.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"+str_list_query+"\n")
        vcf.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\t"+sample_name+"\n")
        num_gene=len(list_vcf)
        for i in range(num_gene):
            f = open(list_vcf[i]["vcf"], "r")
            for line in f:
                if line.startswith("#"):
                    continue
                #vcf.write(line+"\t"+genPresentMark(num_gene,i)+"\n")
                vcf.write(line)
            f.close()
    run_command('gzip -f {}'.format(vcf_file))
def genPresentMark(total, index):
    str_i=""
    for i in range(0,total):
        if not i == index:
            str_i=str_i+"\t"+"."
        else:
            str_i=str_i+"\t1"
    return str_i
