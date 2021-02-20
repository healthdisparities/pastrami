# Notes on running pastrami
# Dec 2020

conda create -n pastrami python=3.8
conda activate pastrami
conda install -c conda-forge pathos numpy pandas

# Haplotype
# run makeHaplotypeBlocks.R

# "database building"
./pastrami.py build --reference-pickle african.pickle --reference-prefix maskedAfricanChrAllReference --threads 20 --haplotype africanHaplotypes.tsv

# european.pickle is the output pickle object

# "querying the database"
./pastrami.py query --reference-pickle african.pickle --query-prefix maskedAfricanChrAllQuery --combined-out african.tsv  --query-out african.query.tsv --threads 20


# Post Pastrami
./postPastramiAfrican.R african.tsv

##################################################
# Goal: Desired commandline call
./pastrami.py all --reference-prefix maskedAfricanChrAllReference --query-prefix maskedAfricanChrAllQuery --out out_calls


##################################################

# Steps to follow with updated pastrami
# Feb 19, 2021

# step 1: build a haplotype file
./pastrami.py -v hapmake --map-dir map_data --hap-file african.hap

# step 2: build database (pickle file)
./pastrami.py -v --threads 20 build --reference-pickle african.pickle --reference-prefix maskedAfricanChrAllReference --haplotype african.hap

# step 3: querying the database
./pastrami.py --threads 20 query --reference-pickle african.pickle --query-prefix maskedAfricanChrAllQuery --combined-out african.tsv  --query-out african.query.tsv

# step 4: aggregate the results
 cat maskedAfricanChrAllQuery.tfam maskedAfricanChrAllReference.tfam > input.tfam
./pastrami.py -v --threads 20 aggregate --pop-group pop2group.txt --pastrami-output african.tsv --pastrami-fam input.tfam
