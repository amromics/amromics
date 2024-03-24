import os
import shutil
import csv
import logging
import multiprocessing
import pandas as pd
import amromics.libs.bioseq as bioseq
import amromics.libs.element_finder as element_finder
from amromics.utils.command import run_command
import amromics.libs.mlst as mlst
logger = logging.getLogger(__name__)
NUM_CORES_DEFAULT = multiprocessing.cpu_count()

def detect_amr_abricate(sample, threads=8, timing_log=None):
    """
    Run abricate to identify resistant genes

    Parameters
    ----------
    sample:
        a dictionary-like object containing various attributes for a sample
    sample_dir: str
        the directory of the sample
    threads: int
        number of threads to use
    overwrite:bool
        whether to overwrite the existing result
    timing_log: str
        log file
    Returns
    -------
        path to resistant gene file
    """

    assembly = sample['assembly']
    resistome_fn = sample['resistome']
    prefix = resistome_fn[:-4]
    dbs=['ncbi','megares','ecoh','argannot','card','resfinder']

    numError=0
    outputfiles=[]
    for db in dbs:
        outfile= prefix + '_'+db+'.tsv'
        cmd = 'abricate --quiet --threads {threads} --nopath --db {db} {infile} > {outfile}'.format(
            threads=threads,
            db=db,
            infile=assembly,
            outfile=outfile)
        if run_command(cmd, timing_log) != 0:
            numError=numError+1
        else:
            outputfiles.append(outfile)

    if numError==len(dbs):
        raise Exception('Error running amr')
    
    combined_tsv = pd.concat([pd.read_csv(f,sep='\t') for f in outputfiles ])
    combined_tsv.sort_values(['SEQUENCE','START'],ascending=[True, True],inplace=True)
    combined_tsv.to_csv(resistome_fn, index=False,sep='\t', encoding='utf-8-sig')
    

def detect_amr_amrfinder(prefix_name,faa_file,fna_file,gff_file,genus=None,species=None,  base_dir='.', db='db/amrfinderplus/data/latest', timing_log=None, threads=0):
    """
        Run AMR analysis, using AMRfinderPlus for searching amr genes, virulome genes and point mutaions.
        :param read_data: result holder
        :param base_dir: working directory
        :return: path to output file in result holder
    """
    if threads == 0:
        threads = NUM_CORES_DEFAULT
    path_out = os.path.join(base_dir,   prefix_name+'_amrfinder')
    if not os.path.exists(path_out):
        os.makedirs(path_out)
    #AMR profiling with CARD. TODO: replace by consensus db later
    ret_out = os.path.join(path_out, prefix_name + '_amr.tsv')
    amr_out = os.path.join(path_out, prefix_name+ '_resistome.tsv')
    virulen_out = os.path.join(path_out, prefix_name + '_virulome.tsv')
    point_out = os.path.join(path_out, prefix_name + '_point.tsv')
    #using abricate
    # cmd = 'abricate --quiet --threads {threads} --nopath --db card {infile} > {outfile}'.format(threads=threads,infile=read_data['assembly'],outfile=amr_out)
    # cmd = "bash -c '{}'".format(cmd)
    # if run_command(cmd) != 0:
    #     return None
    #using build-in function
    #element_finder.search_amr(sample=read_data['assembly'],output=amr_out,threads=threads)
    #using AMRFinderPlus

    #process files in prokka folder, prepare for using amrfinder
    #move files from prokka to temp folder
    temp_dir = os.path.join(base_dir, 'amr_temp')
    if not os.path.exists(temp_dir):
        os.makedirs(temp_dir)

    temp_gff_file=os.path.join(temp_dir,prefix_name+'.gff')
    source_gff_file=None

    source_gff_file=gff_file
    #add Name property to column 9 of gff file (AMRfinder need it!) and remove #fasta section
    if not source_gff_file==None:
        destination= open(temp_gff_file, "w" )
        #source= open( source_gff_file, "r" )
        with open(source_gff_file) as source:
            for line in source:
                if line.startswith('##FASTA'):
                    break
                # no longer to replace ID with Name, Amrfinder official support prokka format
                destination.write( line )

        #source.close()
        destination.close()    
    cmd = 'amrfinder -d {database} -p {faa_file}  -n {fna_file} -g {gff_file} --plus --threads {threads} -o {outfile} -a prokka'\
    .format(
        database=db,
        faa_file=faa_file,
        fna_file=fna_file,
        gff_file=temp_gff_file,
        threads=threads,
        outfile=ret_out
    )
    #full option if has --Genus
    if not genus==None:
        organism = genus.capitalize()
        if not species==None and not species=='':
            organism = species.replace(' ','_')
        organisms = ['Campylobacter', 'Enterococcus_faecalis', 'Enterococcus_faecium', 'Escherichia', 'Klebsiella', 'Salmonella', 'Staphylococcus_aureus', 'Staphylococcus_pseudintermedius', 'Vibrio_cholerae']
        if organism in organisms:
            cmd = 'amrfinder -d {database} -p {faa_file} -O {organism}  -n {fna_file} -g {gff_file} --plus --threads {threads} -o {outfile} -a prokka'\
            .format(
                database=db,
                faa_file=faa_file,
                organism=organism,
                fna_file=fna_file,
                gff_file=temp_gff_file,
                threads=threads,
                outfile=ret_out
            )
        else:
            cmd = 'amrfinder -d {database} -p {faa_file} -n {fna_file} -g {gff_file} --plus --threads {threads} -o {outfile} -a prokka'\
            .format(
                database=db,
                faa_file=faa_file,
                fna_file=fna_file,
                gff_file=temp_gff_file,
                threads=threads,
                outfile=ret_out
            )

    cmd = "bash -c '{}'".format(cmd)
    if run_command(cmd,timing_log) != 0:
        return None

    #clean up:
    if os.path.exists(temp_dir):
        shutil.rmtree(temp_dir)
    #proccess output files:
    virulen=[]
    amr=[]
    point=[]
    header=[]
    with open(ret_out) as tsvfile:
        reader = csv.DictReader(tsvfile, dialect='excel-tab')
        for row in reader:
            header=row.keys()
            if row['Element type']=='VIRULENCE':
                virulen.append(row)
            elif row['Element subtype']=='POINT':
                point.append(row)
            else:
                amr.append(row)
    with open(amr_out, 'w', newline='') as csvfile:

        writer = csv.DictWriter(csvfile, fieldnames=header,delimiter='\t')
        writer.writeheader()
        for row in amr:
            writer.writerow(row)

    with open(point_out, 'w', newline='') as csvfile:

        writer = csv.DictWriter(csvfile, fieldnames=header,delimiter='\t')
        writer.writeheader()
        for row in point:
            writer.writerow(row)

    with open(virulen_out, 'w', newline='') as csvfile:

        writer = csv.DictWriter(csvfile, fieldnames=header,delimiter='\t')
        writer.writeheader()
        for row in virulen:
            writer.writerow(row)

    if os.path.exists(ret_out):
        os.remove(ret_out)

    return amr_out,point_out,virulen_out

###Virulome profiling using abricate with VFDB
def detect_virulome(sample, threads=4):
    """
    Run in-house script to identify virulent genes using VFDB    
    """
    element_finder.search_virulome(sample['assembly'], sample['virulome'],threads=threads)
    

###Find plasmid's origin of replication using abricate with plasmidfinder db
def detect_plasmid(sample, threads=1):
    element_finder.search_plasmid(sample['assembly'], sample['plasmid'],threads=threads)
    

def detect_pmlst(prefix_name,assembly,  base_dir='.', threads=0):
    if threads == 0:
        threads = NUM_CORES_DEFAULT

    path_out = os.path.join(base_dir,prefix_name+'_pmlst' )
    if not os.path.exists(path_out):
        os.makedirs(path_out)

    #Plasmid finder
    pmlst_out = os.path.join(path_out, prefix_name + '_pmlst.tsv')
    if os.path.isfile(pmlst_out):
        return pmlst_out

    # cmd = 'abricate --quiet --threads {threads} --nopath --db plasmidfinder {infile} > {outfile}'.format(threads=threads,infile=read_data['assembly'],outfile=oriREP_out)
    # cmd = "bash -c '{}'".format(cmd)
    # if run_command(cmd) != 0:
    #     return None    

    m=mlst.find_mlst(query_file=assembly,blastdb='db/pmlst/blast/pmlst.fa',mlstdb='db/pmlst/pubmlst',num_threads=threads)
    with open(pmlst_out, 'w') as f:
        f.write("%s\t%s\t%s"%(m['file'],m['scheme'],m['st']))
        for gene in m['profile']:
            f.write("\t%s"%gene)
        f.write("\n")        
    
    return pmlst_out

def detect_insertion_sequence(prefix_name,assembly,  base_dir='.', threads=0):
    """
        Run isescan for searching IS
        :param read_data: result holder
        :param base_dir: working directory
        :return: path to output file in result holder
    """
    # TODO: include overwrite mode
    if threads == 0:
        threads = NUM_CORES_DEFAULT

    path_out = os.path.join(base_dir, prefix_name+"_isescan")
    if not os.path.exists(path_out):
        os.makedirs(path_out)

    isescan_out = os.path.join(path_out, prefix_name + '_is.tsv')
    if os.path.isfile(isescan_out):
        return isescan_out 

    cmd = 'isescan.py --nthread {threads} --seqfile {asm} --output {output}  '.format(
        threads=threads,
        asm=assembly,
        output=path_out
    )
    if run_command(cmd) != 0:
        return None
    #if os.path.exists(path_out+'/prediction'):
    #    shutil.rmtree(path_out+'/prediction')
    #shutil.copytree('prediction', path_out+'/prediction')
    #read_data['is'] = path_out+'/prediction'
    isout=None
    for root, dirs, files in os.walk(path_out, topdown=False):
        for name in files:
            if name.endswith('.raw'):
                isout=os.path.join(root, name)
    #if os.path.exists('prediction'):
    #    shutil.rmtree('prediction')    
    return isout


def detect_integron(prefix_name,assembly, base_dir='.', timing_log=None,threads=0):
    # TODO: include overwrite mode
    if threads == 0:
        threads = NUM_CORES_DEFAULT
    #path_out = os.path.join(base_dir, 'integron_finder_' + read_data['sample_id'])
    path_out = os.path.join(base_dir,  prefix_name+'_integrall' )
    if not os.path.exists(path_out):
        os.makedirs(path_out)
    #Plasmid finder
    integron_out = os.path.join(path_out, prefix_name + '_integron.tsv')

    if os.path.isfile(integron_out):
        return integron_out
    
    element_finder.search_integrall(sample=assembly, output=integron_out,threads=threads)    
    return integron_out


def detect_prophage(prefix_name,faa_file, base_dir='.', timing_log=None,threads=0):
    #TODO: include overwrite mode
    if threads == 0:
        threads = NUM_CORES_DEFAULT

    path_out = os.path.join(base_dir, 'element_finder_' +prefix_name)
    if not os.path.exists(path_out):
        os.makedirs(path_out)

    #Plasmid finder
    prophage_out = os.path.join(path_out,prefix_name + '_prophage.tsv')
    if os.path.isfile(prophage_out):
        return prophage_out   
        
    element_finder.search_prophage(sample=faa_file,output=prophage_out,threads=threads)    
    return prophage_out
