# AMRomics   


## Welcome to AMRomics

**AMRomics** is a software package to analyze microbial genomics data collections.
It provides a pipeline that bundles the current best practice for 
multiple aspects of AMR analyses. The pipeline analysis results can be 
represented and visualized via a web application. The web application also 
provides efficient data management.
 
AMRomics is written in python, it includes the followings dependencies:
 * blast (known to work with 2.10.1+)
 * samtools (1.11)
 * trimmomatic (0.39)
 * spades (3.15.2)
 * shovill (1.1.0)
 * flye (2.9)
 * prokka (1.14.6)
 * mlst (2.19.6)
 * abricate (1.0.1 | Database: vfdb ecoli_vf ecoh card megares resfinder argannot ncbi plasmidfinder)
 * roary (3.13.0) 
 * iqtree (2.1.2)
 * fasttree

## Installation

The simplest method is installed via conda:

0. Make sure a conda version is installed in the computer. If not done already, download and install the appropriate conda, such as anaconda from 
   https://repo.anaconda.com/archive/
   
1. Create a conda environment with all the necessary dependencies: From the repository directory run
```bash
git clone --recursive https://github.com/amromics/amromics.git
cd amromics
conda create -y -c conda-forge -c defaults --name amromics python=3.10 mamba
source activate amromics
mamba install -y -c conda-forge -c bioconda -c anaconda -c etetoolkit -c rpetit3 -c defaults --file requirements.txt
```

2. Recommended: install amromics library into conda environment
```bash
source activate amromics  #If not have already activated
pip install .
```
3. Install panta submodule:
```bash
cd submodules/panta
mamba install -y -c conda-forge -c bioconda -c anaconda -c defaults  --file requirements.txt
pip install .
cd ../..
```   
4. Setup MLST database: AMRomics requires a copy of pubMLST database set up on the folder that AMRomics pipeline is run from. We make available the database in the accompanied file `db.tar.gz` updated on Feb 1, 2024. Just unzip the tarball 
```bash
tar zxvf db.tar.gz
```
Alternatively, you can update the latest database by running the following command line. Note the running make take some time depending on the bandwidth network and the responsiveness of pubMLST server
```
./amr-analysis.py download_db
```

## Usage

### Input preparation

AMRomics takes in as input a collection of bacterial samples in various data types: sequencing reads (fastq format), assembly (fasta format) or assembly with annotations (gff3 where annotations are followed by assembly in fasta). The list of samples in the collection is summarized in a tsv file with the following column headers:
- `sample_id`: An unique ID for each sample. sample_id has to be one word i.e., no space and contains only alphanumerical, _, and - characters
- `sample_name`: A description of the sample
- `input_type`: can take in one of `gff` (for annotations in gff3 format, expecting the assembly at the end of the gff file as in the input from prokka), `asm`, `assembly` (for assembly in fasta format), `Illumina` (for Illumina sequencing reads in fastq format), `pacbio-raw`, `pacbio-hifi`, `pacbio-corr` (for Pacbio sequencing reads in fastq format), `nano-raw`, `nano-hq`, `nano-corr` (for Nanopore sequencing reads in fastq format).
- `files`: path to the input data file(s). For Illumina input type, two fastq files can be supplied if paired-end sequencing is used, and they are separated by a semicolon. For all other input type, only one input file is expected.
- `trim`: if the value is `True` or `yes`, the read data will be trimmed - only applied for Illumina sequencing reads
- `genus`: the Genus of the sample, eg `Escherichia`
- `species`: the species of the sample, eg. `coli`
- `gsize`: an estimate of the genome size, in number, eg 5000000, needed only for subsample Illumina sequencing data to 100x
- `metadata`: any relevant information for the sample, in the format `key1:value1; key2:value2`. For example: `Geographic Location:Houston,USA;Insert Date:8/8/2017;Host Name:Human, Homo sapiens;ampicillin:Resistant;aztreonam:Resistant;ciprofloxacin:Resistant;gentamicin:Susceptible;tetracycline:Susceptible`

We provide several examples in the folder `examples` - See below.


AMRomics pipeline can be invoked with one command line `./amr-analysis.py`. Its usage is as follows:

```bash
./amr-analysis.py  pg -h
usage: amromics pg [-h] [-t THREADS] [-m MEMORY] -c COLLECTION_ID [-n COLLECTION_NAME] -i INPUT [--work-dir WORK_DIR] [--time-log TIME_LOG]
                   [--method METHOD] [--genetree GENETREE] [--progressive PROGRESSIVE] [--tree TREE] [--overwrite OVERWRITE] [--initdb {True,False}]

Pan-genome analysis of a collection

options:
  -h, --help            show this help message and exit
  -t THREADS, --threads THREADS
                        Number of threads to use, 0 for all (default: 0)
  -m MEMORY, --memory MEMORY
                        Amount of memory in Gb to use (default: 30)
  -c COLLECTION_ID, --collection-id COLLECTION_ID
                        Collection ID (default: None)
  -n COLLECTION_NAME, --collection-name COLLECTION_NAME
                        Collection name (default: )
  -i INPUT, --input INPUT
                        Input file (default: None)
  --work-dir WORK_DIR   Working directory (default: data/work)
  --time-log TIME_LOG   Time log file (default: None)
  --method METHOD       panta or roary (default: panta)
  --genetree GENETREE   Run phylogenty for each gene cluster or not (default: False)
  --progressive PROGRESSIVE
                        Run pangenome in progressive mode (default: False)
  --tree TREE           fasttree or iqtree (default: fasttree)
  --overwrite OVERWRITE
                        Force overwrite exist results (default: False)
  --initdb {True,False}
                        Init full database (default: False)

```

#### Examples

We prepare several collections from public data. To download the raw data for these collections

To download the miniature dataset that consists of 5 samples: 1 in sequence assembly, 2 with Illumina sequencing reads, 1 with Nanopore and 1 with Pacbio sequencing reads.
```bash
cd examples/Kp24/raw
./download_kp4.sh
cd ../../
```

To download the 24 sample collection
```bash
cd examples/Kp24/raw
./download_kp24.sh
cd ../../
```

To download the 90 sample collection
```bash
cd examples/Kp24/raw
./download_kp24.sh
cd ../../
```

The following command will run that 24 samples through the pipeline, and import the results
to the web-app for visualization:

```bash
./amr-analysis.py pg --time-log k24_time.log  -t 7 -m 25 -c KpClinicalGRBZ -i examples/Kp24/Kp24.tsv --work-dir data/work  -n "Collection of 24 clinical isolates from Greek and Brazil"
```
#### Run with progressive mode:
```bash
./amr-analysis.py pg --time-log k24_time.log  -t 7 -m 25 -c KpClinicalGRBZ --progressive True -i examples/Kp89/Kp89.tsv --work-dir data/work  -n "Collection of 24+89 clinical isolates from Greek and Brazil"
```
