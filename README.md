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

### Examples

#### Miniature dataset

#### Case study

We prepare a collection of dataset for public data base. To download the raw data,
```bash
cd examples/Kp24/raw
./download.sh
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

<!--

#### Prepare input file
- Data file inputted for analysis needs to be in *.tsv* format 
((To-do: Check if .tsv format is required)) and follows specific requirements. 
Please check the sample input file *data/samples/set1.tsv* for an example.
- Note:
  + Column names need to be as follow:
    - sample_id	
    - sample_name	
    - input_type	
    - files	
    - genus	
    - species	
    - strain	
    - gram	
    - metadata
  + *gram* column should be empty. ((To-do: Delete gram column?))
  + *metadata* is empty or in the format: key1:value1;key2:value2;...  
  For example: Geographic Location:Houston,USA;Insert Date:8/8/2017;Host Name:Human, Homo sapiens;ampicillin:Resistant;aztreonam:Resistant;ciprofloxacin:Resistant;gentamicin:Susceptible;tetracycline:Susceptible
#### Run pipeline and export visualization data to web application

-->
