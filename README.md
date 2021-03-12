# Pastrami
Pastrami is a novel, scalable computational algorithm for rapid human ancestry estimation at population-, subcontinental- and continental-levels.  Pastrami works on two key methodologies: exact haplotype matching and non-negative least square (NNLS) optimization.

Codebase stage: development
Developers and maintainers: Andrew Conley, Lavanya Rishishwar

## Installation
We are working on improving the user experience with Pastrami.  This section will be updated in near future with easier installation procedure for pastrami and it's dependencies.

### Dependencies
Pastrami requires the following OS/programs/modules to run:
* Linux/Mac OS - tested on RedHat OS, expected to run on Mac environments, will possibly work on WSL though Mac and WSL haven't been tested so far
* Plink2 (recommended over Plink, but both will work) - https://www.cog-genomics.org/plink/2.0/; should be accessibly through $PATH
* Python version 3.8 or higher
* Python libraries = Pathos, numpy, scipy, pandas

#### Plink2 installation
Plink2 installation is quite straightforward.  Select the appropriate binary for your processor architecture from this page:https://www.cog-genomics.org/plink/2.0/.  An example installation is shown below for Intel 64-bit architecture:
```
# Download the zip file
wget http://s3.amazonaws.com/plink2-assets/plink2_linux_x86_64_20210203.zip

# Unzip the zip file
unzip plink2_linux_x86_64_20210203.zip

# Optionally, create a local bin and add it to your PATH (ignore these steps if you know what you are doing)
mkdir -p ~/bin
echo "PATH=$HOME/bin:$PATH" >> ~/.bashrc
source ~/.bashrc

# Move the file to your local bin
mv plink2 ~/bin

# Make sure there is an executable by the name plink
ln -s ~/bin/plink2 ~/bin/plink
```

#### Native pip-based installation
Assuming Plink/Plink2 is available system-wide, the following commands will install the python modules and clone the git repo
```
# Install python modules
pip install pathos numpy pandas

# Download Pastrami
git clone https://github.com/healthdisparities/pastrami

# Give it executable permissions
cd pastrami
chmod +x pastrami.py

# Run Pastrami
./pastrami.py
```

#### Conda based installation
Assuming Plink/Plink2 is available system-wide, the following commands will install the python modules and clone the git repo
```
# Create a conda environment with Python version = 3.8
conda create -n pastrami python=3.8

# Activate the newly created conda environment
conda activate pastrami

# Install core packages
conda install -c conda-forge pathos numpy pandas

# Download Pastrami
git clone https://github.com/healthdisparities/pastrami

# Give it executable permissions
cd pastrami
chmod +x pastrami.py

# Run Pastrami
./pastrami.py

# Deactivate environment, when not using pastrami
conda deactivate pastrami
```
## Quickstart guide
This section will be populated in near future with small example datasets and commands for analyzing them.

## Basic usage
### General program structure overview
Pastrami is a five step process that can be run as a single command or using a set of sequential commands.  Each of the five step within Pastrami can be accessed through the following subcommands:
1. *hapmake*: Create haplotype map file
2. *build*: Create a database from reference genotype file (in TPED/TFAM format)
3. *query*: Query the database using a query genotype file (in TPED/TFAM format)
4. *coanc*: Create copying fraction matrix by performing pairwise individual comparisons
5. *aggregate*: Estimate ancestral fractions and aggregate them into fine-grain and sub-continental ancestry

```
usage: pastrami.py [--help] {all,hapmake,build,query,coanc,aggregate} ...

pastrami.py - population scale haplotype copying

optional arguments:
  --help, -h, --h

pastrami.py commands:
  {all,hapmake,build,query,coanc,aggregate}
    all                                      Perform full analysis
    hapmake                                  Create haplotypes
    build                                    Build reference set
    query                                    Query reference set
    coanc                                    Individual v. individual co-ancestry
    aggregate                                Aggregate results and estimate ancestries
```

These steps are all grouped within the ***all*** subcommand.

### Input file format
The genotype and individual information is expected in TPED/TFAM format.  These formats are described at length within Plink's documentation and will not be reiterated here for brevity's sake.
* TPED format description: https://www.cog-genomics.org/plink/1.9/formats#tped
* TFAM format description: https://www.cog-genomics.org/plink/1.9/formats#tfam

If your files are in VCF format, you can convert them into TPED/TFAM format using the following plink command:
```
plink --vcf input.vcf.gz --recode transpose --out output_prefix
```

### A typical run
The ***all*** subcommand is the one that you will use if you have a set of reference individuals and a set of query individuals (presumably with some admixture) and you are interested in performing a one-off ancestry estimation analysis:

```
# Consolidated workflow:
./pastrami.py all --reference-prefix reference --query-prefix query --pop-group pop2group.txt --haplotype chrom.hap --out-prefix full_run -v --threads 20
```
The full description of each commandline option is provided below and can be access from within the script by running the "all" subcommand without any arguments:

```
usage: pastrami.py all [-h] [--reference-prefix <PREFIX>] [--query-prefix <TSV>] [--haplotypes <TSV>]
                       [--pop-group pop2group.txt] [--out-prefix out-prefix] [--map-dir maps/] [--min-snps MinSNPs]
                       [--max-snps MaxSNPs] [--max-rate MaxRate] [--per-individual] [--log-file run.log] [--threads N]
                       [--verbose]

optional arguments:
  -h, --help                                   show this help message and exit

Required Input options:
  --reference-prefix <PREFIX>                  Prefix for the reference TPED/TFAM input files
  --query-prefix <TSV>                         Prefix for the query TPED/TFAM input files
  --haplotypes <TSV>                           File of haplotype positions
  --pop-group pop2group.txt                    File containing population to group (e.g., tribes to region) mapping

Output options:
  --out-prefix out-prefix                      Output prefix (default: pastrami) for creating following sets of files.
                                               <prefix>.pickle, <prefix>_query.tsv, <prefix>.tsv, <prefix>.fam,
                                               <prefix>_fractions.Q, <prefix>_paintings.Q, <prefix>_estimates.Q,
                                               <prefix>_fine_grain_estimates.Q

Optional arguments for hapmake stage:
  --map-dir maps/                              Directory containing genetic maps: chr1.map, chr2.map, etc
  --min-snps MinSNPs                           Minimum number of SNPs in a haplotype block (default: 7)
  --max-snps MaxSNPs                           Maximum number of SNPs in a haplotype block (default: 20)
  --max-rate MaxRate                           Maximum recombination rate (default: 0.3)

Runtime options:
  --per-individual                             Generate per-individual copying rather than per-population copying
  --log-file run.log, -l run.log, --l run.log  File containing log information (default: run.log)
  --threads N, -t N                            Number of concurrent threads (default: 4)
  --verbose, -v                                Print program progress information on screen
```

Consider an example scenario where you want to assess the European ancestry estimates of a group of query individuals.  If your European reference files are named ***euro_ref.tped, euro_ref.tfam*** and your query individual files are named ***my_query.tped, my_query.tfam***, and your population to group mapping is stored within the pop2group.txt file (described more below), your ***all*** subcommand will look as follows:
```
./pastrami.py all --reference-prefix euro_ref --query-prefix my_query --pop-group pop2group.txt --haplotype chrom.hap --out-prefix full_run -v --threads 20
```
I prefer to run pastrami the the verbose flag and my hardware can support 20 threads.  Feel free to turn off the verbose flag if you so desire, and change the number of threads in accordance with your hardware capability.  This command will create the following outputs (listed in order of importance):
1. full_run_estimates.Q: ADMIXTURE-style ancestry fraction estimates at *group level*. This is the output results that you are interested in.
2. full_run_fine_grain_estimates.Q: ADMIXTURE-style ancestry fraction estimates at *population level*. Population level estimates are generally harder to tease apart and are likely to be less reliable.  These results are aggregated in full_run_estimates.Q file.
3. full_run.pickle: this is a Python pickle object that stores the processed reference file (consider this as the "database"). Save this file if you'd like to continue query against the same database in future.
4. full_run_query.tsv: this is your processed query file that can be used.  *This file is not currently being used in any process and will be removed in future.*
5. full_run.tsv: this is the copying fractions file that contains all the reference and query individuals. *This file is not currently being used in any process and will be removed in future.*
6. full_run.fam: these are the individuals that were processed as part of this run
7. full_run_fractions.Q: these are the ancestry fractions (not estimates) from which the final estimates are calculated from. *This file is not currently being used in any process and will be removed in future.*
8. full_run_paintings.Q: these are the ancestry paintings which are used alongside ancestry fractions to calculate ancestry estimates. *This file is not currently being used in any process and will be removed in future.*


### Stepwise execution
Pastrami is designed in a way that the reference files can be processed once and utilized many times in future.  For the following use-case, consider that your reference individuals are African in origin and are stored in ***afro_ref.tped*** and ***afro_ref.tfam*** files.  Your query individuals are still stored in ***my_query.tped*** and ***my_query.tfam*** files.  Your population-to-group mapping is stored in a file named pop2group.txt.  The following series of commands will allow you the pipeline in steps and save the intermediate outputs.

```
# step 1: build a haplotype file
./pastrami.py hapmake --map-dir map_data --haplotypes african.hap -v
# This step shouldn't take more than a few seconds, the remaining step will be time consuming in nature

# step 2: build database (pickle file)
./pastrami.py build --reference-pickle african.pickle --reference-prefix afro_ref --haplotype african.hap -v --threads 20
# african.pickle is the output file

# step 3: querying the database
./pastrami.py query --reference-pickle african.pickle --query-prefix my_query --combined-out african.tsv  --query-out african_query.tsv -v --threads 20 

# step 4: aggregate the results
## combine the individual list for aggregate subcommand
cat afro_ref.tfam my_query.tfam > input.tfam
## run the subcommand
./pastrami.py aggregate --pop-group pop2group.txt --pastrami-output african.tsv --pastrami-fam input.tfam --out-prefix aggregate_output -v --threads 20 
```


### Pop-group mapping
Based on general expectation and our experience in this area, deconvoluting population-level ancestry is very hard (and at times not possible) when populations are genetically extremely close.  In these cases, we strongly recommend grouping closely related populations into groups.  E.g., grouping African tribes into subcontinental African groups.  This information can be supplied to Pastrami in a simple tab-separated file (pop2group.txt file in the previous examples).  The structure of the file looks like this:

```
#Population   Group
GWD           Western
MSL           Western
Yacouba       Western
Ahizi         Western
Fon           Nigerian
Bariba        Nigerian
Yoruba        Nigerian
```
The whitespace separating the column is tab (shown as spaces here for formatting reasons).  Any lines starting with # are ignored.


