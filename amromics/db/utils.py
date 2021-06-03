# -*- coding: utf-8 -*-
"""
    Utilities


Revision history:
----------------
2019-08-16: MDC Created

"""
from __future__ import division, print_function, absolute_import

import logging
import os, glob,shutil
import amromics.libs.bioseq as bioseq
from collections import OrderedDict
from amromics.db.amr_gene import AMRGene, AMRGeneDB
import xml.etree.ElementTree as ET
import wget
import urllib
import gzip
import shutil
from zipfile import ZipFile
import re
import amromics.utils.google_cloud as gcp
logger = logging.getLogger(__name__)


def check_arg_annot(gene_file, protein_file):
    """
    Function to check consistency between the gene file and the protein file

    :param gene_file:
    :param protein_file:
    :return:
    """
    annot_protein_db = {}
    annot_gene_db = {}

    for seq in bioseq.read_sequence_file(protein_file):
        annot_protein_db[seq.name] = seq

    for seq in bioseq.read_sequence_file(gene_file):
        annot_gene_db[seq.name] = seq
    mismatch = 0

    for key in annot_gene_db:
        gene_seq = annot_gene_db[key]

        trans, ORF = bioseq.translate(gene_seq.sequence)
        protein_seq = annot_protein_db[key]

        if not bioseq.identical_protein(trans, protein_seq.sequence):
            mismatch += 1
            print(key, ORF, protein_seq.sequence)
            print(key, ORF, trans)
        if ORF > 0:
            print(key, ORF, gene_seq.length())

    print(mismatch / len(annot_protein_db))


def check_genes(gene_list):
    problems = 0
    big_problems = 0
    not_complete = 0
    for seq in gene_list:
        trans, ORF = bioseq.translate(seq.sequence, reverse=False)
        if trans[-1] != bioseq.STOP_AA:
            not_complete += 1

        cds_len = len(trans) * 3
        if cds_len != len(seq.sequence):
            if seq.length() - cds_len >= 3:
                big_problems += 1
                print(seq.name, ORF, cds_len, seq.sequence, trans)
            else:
                problems += 1
    print(len(gene_list), 'total')
    print(problems, 'problems')
    print(big_problems, 'big problems')
    print(not_complete, 'not_complete')


def extract_arg_annot(db_path='submodules/arg-annot/'):
    """
    Extract genes from arg_annot database
    :param db_path:
    :return:
    """
    current_id = 0

    gene_file = os.path.join(db_path, 'genes.fna')
    #protein_file = os.path.join(db_path, 'proteins.faa')
    #check_arg_annot(gene_file, protein_file)

    # Read the protein file
    amr_db = AMRGeneDB()
    problems = 0
    not_complete = 0
    for seq in bioseq.read_sequence_file(gene_file):
        trans, ORF = bioseq.translate(seq.sequence, reverse=False)
        if trans[-1] != bioseq.STOP_AA:
            print(seq.name)
            print(seq.sequence)
            print(trans)
            print('---------------------------------------')
            not_complete += 1


        toks = seq.name.split(':')
        cds_len = len(trans) * 3
        gene = AMRGene(id='AMG' + str(current_id), cds=seq.sequence[ORF:ORF + cds_len])
        current_id += 1
        gene.external_ids = 'arg_annot_id={}:{}-{}'.format(toks[0], ORF + 1, cds_len)
        if cds_len != len(seq.sequence):
            print(seq.name, ORF, cds_len, len(seq.sequence))
            problems += 1
        gene.external_ids += ';genbank_id='+toks[1]

        toks = toks[2].split('-')
        factor = 1
        if toks[0].startswith('c'):
            start_pos = int(toks[0][1:])
            factor = -1
        else:
            start_pos = int(toks[0])

        gene.external_ids += str(start_pos + ORF)+'-'+str(start_pos + ORF + cds_len * factor)
        amr_db.add_gene(gene)
        #seq.set_desc('genbank_id='+toks[1] + ':' + toks[2]+':+')  # TODO make sure it is + strand
        #annot_db.append(seq)
    print(problems, ' problems')
    print(not_complete, 'not_complete')
    return amr_db




def extract_arpcard(path='submodules/card/data/'):
    for nuc_file in glob.glob(path + 'nucleotide_fasta_protein_*.fasta'):
        protein_file = nuc_file.replace('nucleotide', 'protein')

        nuc_reader = bioseq.read_sequence_file(nuc_file)
        protein_reader = bioseq.read_sequence_file(protein_file)
        aprcard = []
        for nuc_seq in nuc_reader:
            protein_seq = protein_reader.__next__()
            nuc_id = 'None'
            toks = nuc_seq.name.split('|')
            for tok in toks:
                if tok.startswith('ARO:'):
                    nuc_id = tok
            protein_id = 'None'
            for tok in protein_seq.name.split('|'):
                if tok.startswith('ARO:'):
                    protein_id = tok

            if nuc_id != protein_id:
                print('Not expected ', nuc_id, '!=', protein_id)
            translate_protein, _ = bioseq.translate(nuc_seq.sequence)

            if not bioseq.identical_protein(translate_protein, protein_seq.sequence):
                #print(translate_protein, protein_seq.sequence)
                print(nuc_seq.name, nuc_seq.length() / 3,  len(translate_protein), ' vs ', protein_seq.name,  protein_seq.length())
                print(translate_protein)
                print(protein_seq.sequence)

            pos = toks[3].split('-')
            start = int(pos[0]) + 1
            desc = 'genbank_id=' + toks[1] + ':' + str(start) + '-' + pos[1] + ':' + toks[2] + ';CARD=' + toks[
                4] + ';' + nuc_seq.desc
            nuc_seq.set_desc(desc)

            aprcard.append(nuc_seq)
    return aprcard



def extract_resfinder(db_path='submodules/resfinder_db/'):
    res_finder_genes = []
    config_file = os.path.join(db_path, 'config')
    resfinder_gene_fn = os.path.join(db_path, 'genes.fna')
    with open(config_file, 'r') as fn, open(resfinder_gene_fn, 'w') as resfn:
        for line in fn:
            if line.startswith('#'):
                continue
            toks = line.split('\t')
            if len(toks) < 3:
                continue

            prefix = toks[0]
            amr_class = toks[1]
            gene_file = os.path.join(db_path, prefix + '.fsa')

            for seq in bioseq.read_sequence_file(gene_file):
                name = seq.name
                toks = name.split('_')
                for i in range(len(toks)):
                    if toks[i].isdigit():
                        break
                desc = 'genbank_id=' + ('_'.join(toks[i+1:]))
                seq.set_desc(desc+';'+amr_class)
                res_finder_genes.append( seq)
                resfn.write(seq.format_fasta())
    return res_finder_genes




def extract_mega_res(db_path='submodules/megares/megares_v1.01/'):
    """TODO: to fix this"""
    gene_file = os.path.join(db_path, 'megares_database_v1.01.fasta')
    megares_genes = OrderedDict()

    for seq in bioseq.read_sequence_file(gene_file):
        megares_genes[seq.name] = seq

    megares_links = {}
    link_file = os.path.join(db_path, 'megares_to_external_header_mappings_v1.01.tsv')
    with open(link_file,'r') as fn:
        ## ignore the first line
        line = fn.readline()
        for line in fn.readlines():
            toks = line.strip().split('\t')
            if len(toks) != 3:
                logger.warning('Expect 3 tokens, found only {}: {}'.format(len(toks), line.strip()))
                exit(1)
            key = toks[1]
            if key not in megares_links.keys():
                megares_links[key] = {}

            megares_links[key][toks[0]] = toks[2]
    return megares_genes, megares_links

def downloadfile(url,savefile):
    #save the bak file before download new file if file exists
    if os.path.exists(savefile):
        os.rename(savefile,savefile+'.bak')
    try:
        print('\nDownloading file {}'.format(savefile))
        wget.download(url,savefile)
        if os.path.exists(savefile+'.bak'):
            os.remove(savefile+'.bak')
        return savefile
    except BaseException as error:
        print(str(error))
        if os.path.exists(savefile+'.bak'):
            os.rename(savefile+'.bak',savefile)
        return None
def get_mlst_db():
    '''
        Download PubMLST from ftp folder, read dbases.xml to get url for downloading
        combine sequences into one .fa file and make blast db
    '''
    url = 'https://pubmlst.org/data/dbases.xml'
    webdata = urllib.request.urlopen(url)
    tree = ET.parse(webdata)
    root = tree.getroot()
    mlst_dir='db/mlst'
    pubmlst_dir=os.path.join(mlst_dir,'pubmlst')
    if not os.path.exists(pubmlst_dir):
        os.makedirs(pubmlst_dir)

    for species in root:
        scheme_dir=''
        for profiles in species.iter('profiles'):
            for child in profiles:
                if child.tag=='url':
                    #scheme=os.path.basename(child.text).strip().replace('.txt','')
                    toks = child.text.split("/")
                    scheme_name = toks[toks.index("db")+1].split("_")[1]
                    scheme_num = toks[toks.index("schemes")+1]
                    scheme = scheme_name + (("_"+scheme_num) if scheme_num!="1" else "")
                    scheme_dir=os.path.join(pubmlst_dir,scheme)
                    if not os.path.exists(scheme_dir):
                        os.makedirs(scheme_dir)
                    #print(os.path.basename(child.text).strip())
                    downloadfile(child.text,os.path.join(scheme_dir,scheme+".txt"))
        for loci in species.iter('loci'):
            for locus in loci:
                for child in locus:
                    if child.tag=='url':
                        print(os.path.join(scheme_dir,locus.text.strip()))
                        downloadfile(child.text,os.path.join(scheme_dir,locus.text.strip()+".tfa"))
    ignore_species={'afumigatus','blastocystis','calbicans','cglabrata','ckrusei','ctropicalis','csinensis','kseptempunctata','sparasitica','tvaginalis'}
    for rm_species in ignore_species:
        rm_dir=os.path.join(pubmlst_dir,rm_species)
        if os.path.exists(rm_dir):
            shutil.rmtree(rm_dir)
    #make blastdb:
    blast_dir=os.path.join(mlst_dir,'blast')
    if not os.path.exists(blast_dir):
        os.makedirs(blast_dir)
    mlst_file=os.path.join(blast_dir,'mlst.fa')
    try:
        os.remove(mlst_file)
    except OSError:
        pass
    for root, dirs, files in os.walk(pubmlst_dir):
        for _file in files:
            if _file.endswith('.tfa'):
                infile=str(root)+'/'+_file
                scheme = os.path.basename(str(root))
                catcmd='cat {tfa} | sed -e \'s/>/>{scheme}./g\' >>{mlst}'.format(
                    scheme=scheme,
                    tfa=infile,
                    mlst=mlst_file

                )
                #print(catcmd)
                os.system(catcmd)
    cmd="makeblastdb -hash_index -in {blastfile} -dbtype nucl -title PubMLST -parse_seqids".format(
        blastfile=mlst_file
    )
    print (cmd)
    os.system(cmd)
def get_virulome_db():
    '''
        Download vfdb from mgc.ac.cn, rewrite sequence'name
    '''
    vfdb_dir='db/vfdb'
    if not os.path.exists(vfdb_dir):
        os.makedirs(vfdb_dir)
    vfdb_file_zip=os.path.join(vfdb_dir,'sequences.gz')
    vfdb_file=os.path.join(vfdb_dir,'sequences.fa')
    vfdb_file_out=os.path.join(vfdb_dir,'sequences')
    vfdb_file_zip=downloadfile('http://www.mgc.ac.cn/VFs/Down/VFDB_setA_nt.fas.gz',vfdb_file_zip)
    if not vfdb_file_zip==None:
        with gzip.open(vfdb_file_zip, 'rb') as f_in:
            with open(vfdb_file, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        with open(vfdb_file_out,'w') as f:
            for seq in bioseq.read_sequence_file(vfdb_file):
                #print(seq.name+';'+seq.desc)
                accession=''
                z=re.match(r"^(\w+)\(\w+\|(\w+)(\.\d+)?\)", seq.name)
                if z:
                    accession=z[2]
                else:
                    accession=seq.name
                z=re.match(r"^\((.*?)\)", seq.desc)
                gene=''
                if z:
                    gene=z[1]
                else:
                    gene=accession
                seq.set_name('vfdb~~~'+gene+'~~~'+accession)
                f.write(seq.format_fasta(max_line=100))
    #clean up
    if os.path.exists(vfdb_file_zip):
        os.remove(vfdb_file_zip)
    if os.path.exists(vfdb_file):
        os.remove(vfdb_file)
    #make blast db
    cmd="makeblastdb -in {blastfile} -dbtype nucl -title vfdb".format(
        blastfile=vfdb_file_out
    )
    print (cmd)
    os.system(cmd)
def get_plasmidfinder():
    plasmidfinder_dir='db/plasmidfinder'
    if not os.path.exists(plasmidfinder_dir):
        os.makedirs(plasmidfinder_dir)
    url='https://bitbucket.org/genomicepidemiology/plasmidfinder_db/get/HEAD.zip'
    plasmidfinder_zip=os.path.join(plasmidfinder_dir,'plasmidfinder.zip')
    plasmidfinder_unzip=os.path.join(plasmidfinder_dir,'temp')
    plasmidfinder_out=os.path.join(plasmidfinder_dir,'sequences')
    plasmidfinder_zip=downloadfile(url,plasmidfinder_zip)
    if not plasmidfinder_zip==None:

        with ZipFile(plasmidfinder_zip, 'r') as zipObj:
            zipObj.extractall(plasmidfinder_unzip)
        with open(plasmidfinder_out,'w') as f:
            for root, dirs, files in os.walk(plasmidfinder_unzip):
                for _file in files:
                    if _file.endswith(('.fsa')):
                        #print(str(root)+'/'+_file)
                        for seq in bioseq.read_sequence_file(str(root)+'/'+_file):
                            accession=''
                            gene=''
                            #parse string to get gene and acc
                            z=re.match(r"^(.*)_(([A-Z]+|NC_)\d+(\.\d+)?)", seq.name)
                            if z:
                                accession=z[2]
                                gene=z[1]
                            else:
                                accession=seq.name

                            seq.set_desc(seq.name)
                            seq.set_name('plasmidfinder~~~'+gene+'~~~'+accession+'~~~')
                            f.write(seq.format_fasta(max_line=100))
    #cleanup
    if os.path.exists(plasmidfinder_zip):
        os.remove(plasmidfinder_zip)
    if os.path.exists(plasmidfinder_unzip):
        shutil.rmtree(plasmidfinder_unzip)
    cmd="makeblastdb -in {blastfile} -dbtype nucl -title plasmidfinder".format(
        blastfile=plasmidfinder_out
    )
    print (cmd)
    os.system(cmd)
def get_pmlst():
    '''
    Get db from ge.cbs.dtu.dk, format db like https://github.com/tseemann/mlst
    '''
    pmlst_dir='db/pmlst'
    if not os.path.exists(pmlst_dir):
        os.makedirs(pmlst_dir)
    blast_dir=os.path.join(pmlst_dir,'blast')
    pubmlst_dir=os.path.join(pmlst_dir,'pubmlst')
    if not os.path.exists(blast_dir):
        os.makedirs(blast_dir)
    url='https://bitbucket.org/genomicepidemiology/pmlst_db/get/HEAD.zip'
    pmlst_zip=os.path.join(pmlst_dir,'pmlst.zip')
    pmlst_unzip=os.path.join(pmlst_dir,'temp')

    pmlst_zip=downloadfile(url,pmlst_zip)
    pmlst_file=os.path.join(blast_dir,'pmlst.fa')
    try:
        os.remove(pmlst_file)
    except OSError:
        pass
    if not pmlst_zip==None:
        with ZipFile(pmlst_zip, 'r') as zipObj:
            zipObj.extractall(pmlst_unzip)
        for root, dirs, files in os.walk(pmlst_unzip):
            for _file in files:
                if _file.endswith(('.fsa')):
                    #print(str(root)+'/'+_file)
                    infile=str(root)+'/'+_file
                    scheme = os.path.basename(_file).strip().replace('.fsa','')
                    catcmd='cat {fsa} | sed -e \'s/>/>{scheme}./g\' >>{pmlst}'.format(
                        scheme=scheme,
                        fsa=infile,
                        pmlst=pmlst_file

                    )
                    #print(catcmd)
                    os.system(catcmd)
                if _file.endswith(('.txt.clean')):
                    scheme=os.path.basename(_file).strip().replace('.txt.clean','')
                    scheme_dir=os.path.join(pubmlst_dir,scheme)
                    if not os.path.exists(scheme_dir):
                        os.makedirs(scheme_dir)
                    #move profile .txt to scheme folder
                    shutil.move(os.path.join(root,_file), os.path.join(scheme_dir,scheme+'.txt'))
    #cleanup
    if os.path.exists(pmlst_zip):
        os.remove(pmlst_zip)
    if os.path.exists(pmlst_unzip):
        shutil.rmtree(pmlst_unzip)
    cmd="makeblastdb -hash_index -in {blastfile} -dbtype nucl -title pMLST -parse_seqids".format(
        blastfile=pmlst_file
    )
    print (cmd)
    os.system(cmd)
    #https://bitbucket.org/genomicepidemiology/intfinder_db
def get_integron():
    integron_dir='db/integron'
    if not os.path.exists(integron_dir):
        os.makedirs(integron_dir)
    url='https://bitbucket.org/genomicepidemiology/intfinder_db/get/HEAD.zip'
    integron_zip=os.path.join(integron_dir,'integron.zip')
    integron_unzip=os.path.join(integron_dir,'temp')
    integron_out=os.path.join(integron_dir,'sequences')
    integron_zip=downloadfile(url,integron_zip)
    if not integron_zip==None:

        with ZipFile(integron_zip, 'r') as zipObj:
            zipObj.extractall(integron_unzip)
        with open(integron_out,'w') as f:
            for root, dirs, files in os.walk(integron_unzip):
                for _file in files:
                    if _file.endswith(('.fsa')):
                        #print(str(root)+'/'+_file)
                        for seq in bioseq.read_sequence_file(str(root)+'/'+_file):
                            accession=''
                            gene=''
                            #parse string to get gene and acc
                            z=re.match(r"^(.*)_(([A-Z]+|NC_)\d+(\.\d+)?)", seq.name)
                            if z:
                                accession=z[2]
                                gene=z[1]
                            else:
                                accession=seq.name

                            seq.set_desc(seq.name)
                            seq.set_name('integron~~~'+gene+'~~~'+accession+'~~~')
                            f.write(seq.format_fasta(max_line=100))
    #cleanup
    if os.path.exists(integron_zip):
        os.remove(integron_zip)
    if os.path.exists(integron_unzip):
        shutil.rmtree(integron_unzip)
    cmd="makeblastdb -in {blastfile} -dbtype nucl -title integron".format(
        blastfile=integron_out
    )
    print (cmd)
    os.system(cmd)
def get_kraken2():
    '''
        Get kraken db
    '''
    #url='https://genome-idx.s3.amazonaws.com/kraken/k2_standard_8gb_20200919.tar.gz'
    #for better download speed
    url='ftp://ftp.ccb.jhu.edu/pub/data/kraken2_dbs/old/minikraken2_v1_8GB_201904.tgz'
    kraken2_db='db/kraken2/'
    if not os.path.exists(kraken2_db):
        os.makedirs(kraken2_db)
    k2std_zip=os.path.join(kraken2_db,'k2std.tar.gz')
    k2std_unzip=os.path.join(kraken2_db,'k2std')
    k2std_zip=downloadfile(url,k2std_zip)
    if not k2std_zip==None:
        shutil.unpack_archive(k2std_zip, k2std_unzip)
        for (root,dirs,files) in os.walk(k2std_unzip, topdown=True):
            for f in files:
                filepath=os.path.join(root, f)
                shutil.move(filepath, os.path.join(k2std_unzip,f))
    if os.path.exists(k2std_zip):
        os.remove(k2std_zip)
def get_prophage():
    #url='http://phaster.ca/downloads/prophage_virus.db'
    url='https://gdurl.com/y2eu'
    prophage_dir='db/prophage'
    if not os.path.exists(prophage_dir):
        os.makedirs(prophage_dir)
    prophage_file=os.path.join(prophage_dir,'sequences')
    prophage_file=downloadfile(url,prophage_file)
    with open(integron_out,'w') as f:
        for seq in bioseq.read_sequence_file(prophage_file):
            accession=''
            gene=''
            #parse string to get gene and acc
            z=re.match(r"^PROPHAGE_(.*)\|(.*)\|(.*)\|(([A-Z]+|NP_)\d+(\.\d+)?)", seq.name)
            if z:
                accession=z[4]
                gene=z[1]
            else:
                accession=seq.name

                seq.set_desc(seq.name)
            seq.set_name('integron~~~'+gene+'~~~'+accession+'~~~')
            f.write(seq.format_fasta(max_line=100))
    cmd="makeblastdb -in {blastfile} -dbtype prot -title prophage".format(
        blastfile=prophage_file
    )
    print (cmd)
    os.system(cmd)
def update_amrfinderplus():
    '''
        Get amrfinderplus db and save to custom folder
    '''
    amrfinderplus_dir='db/amrfinderplus/data/'
    if not os.path.exists(amrfinderplus_dir):
        os.makedirs(amrfinderplus_dir)

    cmd="amrfinder_update -d {amrfinderplusdb}".format(amrfinderplusdb=amrfinderplus_dir)
    os.system(cmd)

def get_integrall_db():
    integrall_dir='db/integrall'
    if not os.path.exists(integrall_dir):
        os.makedirs(integrall_dir)
    url='http://integrall.bio.ua.pt/?getFastaAll'
    integrall_file=os.path.join(integrall_dir,'integrall.fasta')
    integrall_file_out=os.path.join(integrall_dir,'sequences')
    integrall_file=downloadfile(url,integrall_file)
    if not integrall_file==None:
        with open(integrall_file_out,'w') as f:
            for seq in bioseq.read_sequence_file(integrall_file):
                print(seq.name)
                if '?' in seq.name:
                    continue
                z= seq.name.split('|')
                y=seq.desc.split('|')
                if 'Δ' in y[1]:
                    y[1]=y[1].replace('Δ','beta-')
                seq.set_name('integrall~~~'+y[1]+'~~~'+z[1])
                seq.set_desc(y[0])
                f.write(seq.format_fasta(max_line=100))
    cmd="makeblastdb -in {blastfile} -dbtype nucl -title integrall".format(
        blastfile=integrall_file_out
    )
    print (cmd)
    os.system(cmd)
def setup_db():
    get_kraken2()
    get_mlst_db()
    get_virulome_db()
    get_plasmidfinder()
    update_amrfinderplus()
    get_pmlst()
    get_integrall_db()
    get_integron()
    #get_prophage()



    #A dictionary mapping forein ID into mine
    #my_dict = {}

    # arg_annot_path = '/misc/workspace/amromics/bifx/submodules/arg-annot/'
    # arg_annot_genes = extract_arg_annot(db_path=arg_annot_path)
    # logger.info('Extract {} genes from arg_annot'.format(len(arg_annot_genes)))
    #
    # card_path = '/misc/workspace/amromics/bifx/submodules/card/data/'
    # card_genes = extract_arpcard(card_path)
    # logger.info('Extract {} genes from card'.format(len(card_genes)))
    #
    # resfinder_path = '/misc/workspace/amromics/bifx/submodules/resfinder_db'
    # resfinder_genes = extract_resfinder(resfinder_path)
    # logger.info('Extract {} genes from resfinder'.format(len(resfinder_genes)))
    #
    # megares_path = '/misc/workspace/amromics/bifx/submodules/megares/megares_v1.01/'
    # megares_genes, megares_links = extract_mega_res(db_path=megares_path)
    # logger.info('Extract {} genes from megares'.format(len(megares_genes)))


    #resfinder_genes = get_resfinder('/home/minhduc/Projects/AMRBrowser/amrbpipeline/submodules/resfinder_db')
    #arg_annot = get_arg_annot()
    #db_path = '/home/minhduc/Projects/AMRBrowser/amrbpipeline/submodules/arg-annot/'
    #gene_file = os.path.join(db_path, 'genes.fna')
    #protein_file = os.path.join(db_path, 'proteins.faa')
    #check_arg_annot(gene_file, protein_file)

    #get_arpcard('/home/minhduc/Projects/AMRBrowser/amrbpipeline/submodules/card/data/')

    #for mega_name in megares_links.keys():
    #    link_dict = megares_links[mega_name]
    #    mega_seq = megares_genes[mega_name]
    #    for db_name in link_dict.keys():
    #        if db_name == 'Resfinder':
    #            resfinder_name = link_dict[db_name]
    #            print(mega_seq.sequence)
    #            resfinder_gene = resfinder_genes[resfinder_name]
                #print (mega_seq.sequence)
#                print(resfinder_gene.sequence == mega_seq.sequence)
    #            break

    return None


def setup_db_from_bucket(google_json_file, bucket_name='amr-pipeline-db', download_dir='db'):
    # first, if db folder exists, move it to another folder
    if os.path.exists(download_dir):
        shutil.move(download_dir, download_dir + '_old')
    gcp.gcp_copy_bucket(google_json_file, bucket_name, download_dir)
    # create symlink for download_dir + '/amrfinderplus/data/latest'
    amrfinder_latest = download_dir + '/amrfinderplus/data/latest'
    amr_datadir = download_dir + '/amrfinderplus/data'
    if not os.path.exists(amrfinder_latest):
        # create this symlink poitng to the first dir in the parent dir
        amrdatasets = [ name for name in os.listdir(amr_datadir) if os.path.isdir(os.path.join(amr_datadir, name)) ]
        os.symlink(amrdatasets[0], amrfinder_latest)

if __name__ == "__main__":
    setup_db()

    # logging.basicConfig(level=logging.DEBUG,
    #                     format='%(asctime)s [%(name)s] %(levelname)s : %(message)s')
    # logger = logging.getLogger(__name__)
    # for gene_file in [#'/home/minhduc/Projects/amromics_bifx/submodules/arg-annot/genes.fna',
    #                   '/misc/Backup/Dropbox/Projects/AMRBrowser/amrbpipeline/submodules/resfinder_db/genes.fna',
    #                   # '/home/minhduc/Projects/amromics_bifx/submodules/megares/megares_v1.01/megares_database_v1.01.fasta',
    #                   # '/misc/Backup/Dropbox/Projects/AMRBrowser/amrbpipeline/submodules/card/data/nucleotide_fasta_protein_homolog_model.fasta',
    #                   # '/misc/Backup/Dropbox/Projects/AMRBrowser/amrbpipeline/submodules/card/data/nucleotide_fasta_protein_knockout_model.fasta',
    #                   # '/misc/Backup/Dropbox/Projects/AMRBrowser/amrbpipeline/submodules/card/data/nucleotide_fasta_protein_overexpression_model.fasta',
    #                   # '/misc/Backup/Dropbox/Projects/AMRBrowser/amrbpipeline/submodules/card/data/nucleotide_fasta_protein_variant_model.fasta'
    #                   ]:
    #     print(gene_file)
    #     gene_list = []
    #     for seq in bioseq.read_sequence_file(gene_file):
    #         gene_list.append(seq)
    #     check_genes(gene_list)


    print('Done')





###################################################################################

def check_gene_protein(cds_list, protein_list):
    assert(len(cds_list) == len(protein_list))
    for i, gene in enumerate(cds_list):
        trans, ORF = bioseq.translate(gene.sequence)
        if not bioseq.identical_protein(trans, protein_list[i].sequence):
            print(gene.name)


def update_resfinder_db(res_finder_path):
    pass





def cluster_genes_cdhit(protein_file, gene_file):
    """
    Clustering genes
    :param protein_file:
    :param gene_file:
    :return:
    """
    cdhit_out = os.path.join(os.path.dirname(protein_file), 'cdhit_proteins')
    cmd = 'cdhit -i {protein_file} -o {cdhit_out}'.format(
        protein_file=protein_file,
        cdhit_out=cdhit_out
    )
    logger.info('Running ' + cmd)
    os.system(cmd)

    cdhit_out = os.path.join(os.path.dirname(gene_file), 'cdhit_genes')
    cmd = 'cdhit-est -i {gene_file} -o {cdhit_out}'.format(
        gene_file=gene_file,
        cdhit_out=cdhit_out
    )
    logger.info('Running ' + cmd)
    os.system(cmd)

    #protein_cltr = cdhit_out + '.clstr'
    #with open(protein_cltr,'r') as cltr_fn:



#VRSKNFSWRYSLAATVLLLSPFDLLASLGMDMYLPAVPFMPNALGTTASTVQLTLATYLVMIGAGQLLFGPLSDRLGRRPVLLGGGLAYVVASMGLAFTSLAEVFLGLRILQACGASACLVSTFATVRDIYAGREESNVIYGILGSMLAMVPAVGPLLGALVDMWLGWRAIFAFLGLGMIAASAAAWRFWPETRVQRVTGLQWSQLLLPVKCLNFWLYTLCYAAGMGSFFVFFSIAPGLIMGRQGVSQLGFSLLFATVAIAMVFTARFMGRVIPKWGSPSVLRMGMGCLIAGAVLLAITEIWASQSVLGFIAPMWLVGIGVATAVSVAPNGALQGFDHVAGTVTAVYFCLGGVLLGSIGTLIISLLPRNTAWPVVVYCLTLATVVLGLSCVSRAKGSRGQGEHDVVALQSAESTSNPNR
#VRSKNFSWRYSLAATVLLLSPFDLLASLGMDMYLPAVPFMPNALGTTASTVQLTLATYLVMIGAGQLLFGPLSDRLGRRPVLLGGGLAYVVASMGLAFTSLAEVFLGLRILQACGASACLVSTFATVRDIYAGREESNVIYGILGSMLAMVPAVGPLLGALVDMWLGWRAIFAFLGLGMIAASAAAWRFWPETRVQRVTGLQWSQLLLPVKCLNFWLYTLCYAAGMGSFFVFFSIAPGLIMGRQGVSQLGFSLLFATVAIAMVFTARFMGRVIPKWGSPSVLRMGMGCLIAGAVLLAITEIWASQSVLGFIAPMWLVGIGVATAVSVAPNGALQGFDHVAGTVTAVYFCLGGVLLGSIGTLIISLLPRNTAWPVVVYCLTLATVVLGLSCVSRAKGSRGQGEHDVVALQSAESTSNPNR


def get_arg_annot(db_path='submodules/arg-annot/'):
    """
    Not completed
    """
    gene_file = os.path.join(db_path, 'genes.fna')
    protein_file = os.path.join(db_path, 'proteins.faa')
    annot_protein_db = {}
    annot_gene_db = {}

    # Read the protein file
    for seq in bioseq.read_sequence_file(protein_file):
        annot_protein_db[seq.name] = seq

    for seq in bioseq.read_sequence_file(gene_file):
        annot_gene_db[seq.name] = seq

    # Sort the protein sequences by length -- cdhit would do that anyway
    sorted_db = sorted(annot_protein_db.values(), key=lambda s: s.length(), reverse=True)
    sorted_db_file = os.path.join(db_path, 'protein_sorted.faa')
    name_dict = {}

    with open(sorted_db_file, 'w') as ofn:
        for i in range(len(sorted_db)):
            seq = sorted_db[i]
            seq.set_desc('ARG-ANNOT=' + seq.name)
            sname = 'P' + str(i)
            name_dict[sname] = seq.name
            seq.set_name(sname)
            ofn.write(seq.format_fasta())

    #Run cdhit
    cdhit_out = os.path.join(db_path, 'protein_sorted_cdhit')
    cmd = 'cdhit -i {protein_file} -o {cdhit_out}'.format(
        protein_file=sorted_db_file,
        cdhit_out=cdhit_out
    )
    logger.info('Running ' + cmd)
    ret = os.system(cmd)
    if ret != 0:
        exit(ret)

    cluster_dict = OrderedDict()
    cluster_prefix = '>Cluster '
    cluster_id = ''
    with open(cdhit_out + '.clstr', 'r') as cdhit_fn:
        for line in cdhit_fn:
            line = line.strip()
            if line.startswith(cluster_prefix):
                cluster_id = line[len(cluster_prefix):]
                cluster_id = 'ARP' + ('0' * (6 - len(cluster_id))) + cluster_id
                cluster_dict[cluster_id] = []
            else:
                toks = line.split()
                name = toks[2][1:-3]
                #print(toks[2], name, cluster_id)
                cluster_dict[cluster_id].append(name_dict[name])

    fam_protein_file = 'database/fam_protein.faa'
    fam_gene_file = 'database/fam_gene.fna'

    tot_protein_file = 'database/tot_protein.faa'
    tot_gene_file = 'database/tot_gene.faa'

    with open(fam_gene_file, 'w') as fam_gene, \
        open(fam_protein_file, 'w') as fam_protein, \
        open(tot_gene_file, 'w') as tot_gene, \
        open(tot_protein_file, 'w') as tot_protein:
        for cluster_id in cluster_dict.keys():
            cluster = cluster_dict[cluster_id]
            for idx, name in enumerate(cluster):
                protein_seq = annot_protein_db[name]
                protein_seq.set_name(cluster_id + '.' + str(idx))
                tot_protein.write(protein_seq.format_fasta())

                gene_seq = annot_gene_db[name]
                gene_seq.set_desc(protein_seq.get_desc())
                gene_seq.set_name(cluster_id + '.' + str(idx))
                tot_gene.write(gene_seq.format_fasta())

                if idx == 0:
                    protein_seq.set_name(cluster_id)
                    protein_seq.set_desc('')
                    fam_protein.write(protein_seq.format_fasta())
                    gene_seq.set_name(cluster_id)
                    gene_seq.set_desc('')
                    fam_gene.write(gene_seq.format_fasta())

                    translated = bioseq.translate(gene_seq.sequence)
                    if translated[-1:] == '_':
                        translated = translated[:-1]
                    print(cluster_id,translated.count('_'), translated == protein_seq.sequence, translated, protein_seq.sequence)

                # end if
            # end for
    return annot_protein_db
# blastn -task blastn -subject submodules/arg-annot/gene_sorted.fna -query submodules/arg-annot/gene_sorted.fna -outfmt "6 qseqid qlen qstart qend sseqid slen sstart send length pident nident evalue score bitscore qseq sseq" > submodules/arg-annot/gene_sorted.blastn
# blastp -matrix BLOSUM80 -subject proteins.faa -query proteins.faa -outfmt "6 qseqid  qlen qstart qend sseqid slen sstart send length pident nident evalue score bitscore qseq sseq"
