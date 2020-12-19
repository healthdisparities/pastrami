#!/usr/bin/env python3
"""Pastrami - Population scale haplotype copying script"""

__author__ = "Andrew Conley, Lavanya Rishishwar"
__copyright__ = "Copyright 2020, Andrew Conley, Lavanya Rishishwar"
__credits__ = ["Andrew Conely", "Lavanya Rishishwar"]
__license__ = "GPL"
__version__ = "0.1"
__maintainer__ = "Andrew Conley, Lavanya Rishishwar"
__email__ = "aconley@ihrc.com; lrishishwar@ihrc.com"
__status__ = "Development"

from argparse import ArgumentParser, HelpFormatter
import datetime
import logging
import math
import os.path
import pandas as pd
import pathos.multiprocessing as mp
import pickle
import subprocess
import sys

VERSION = 0.1
PROGRAM_NAME = "pastrami.py"


class Colors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'


class HaplotypeMaker:
    pass


class Analysis:

    def __init__(self, opts, unknown_args):

        self.chromosomes = list(range(1, 23))
        self.min_haplotype_occurences = 0

        self.sub_command = opts.sub_command

        self.reference_tped_file = None
        self.reference_tfam_file = None
        self.query_tped_file = None
        self.query_tfam_file = None
        self.reference_tfam = None
        self.query_tfam = None
        self.haplotypes = None
        self.reference_individual_populations = None
        self.reference_population_counts = None
        self.query_copying_fractions = None
        self.combined_copying_fractions = None

        # Build sub-command options
        if self.sub_command == 'build':
            self.reference_pickle_output_file = opts.reference_pickle_out
            self.reference_prefix = opts.reference_prefix
            self.haplotype_file = opts.haplotypes
            self.reference_output_file = opts.reference_out
            self.reference_tpeds = pd.Series([], dtype=pd.StringDtype())
            self.reference_background = pd.Series([], dtype=pd.StringDtype())

        # Query options
        if self.sub_command == 'query':
            self.reference_pickle_file = opts.reference_pickle
            self.query_prefix = opts.query_prefix
            self.query_output_file = opts.query_out
            self.combined_output_file = opts.combined_out
            self.query_tpeds = pd.Series([], dtype=pd.StringDtype())
            self.query_individuals = None

        # Co-ancestry options
        if self.sub_command == 'coanc':
            self.haplotype_file = opts.haplotypes
            self.reference_prefix = opts.reference_prefix
            self.query_prefix = opts.query_prefix

            self.refernce_output_file = opts.reference_out
            self.query_output_file = opts.query_out
            self.query_combined_file = opts.combined_out

        # Either values
        self.reference_individuals = None
        self.reference_populations = None
        self.reference_haplotype_counts = None
        self.reference_haplotype_fractions = None
        self.reference_copying_fractions = None

        # Any errors we encounter
        self.errors = []

        # The actual queue for the analysis
        self.analysis = []

        self.threads = opts.threads
        self.verbosity = 1
        self.fake_run = False

        # Verbosity levels and colors
        self.error_color = Colors.FAIL
        self.main_process_verbosity = 1
        self.warning_color = Colors.WARNING
        self.warning_verbosity = 1
        self.main_process_color = Colors.OKGREEN
        self.sub_process_verbosity = 2
        self.sub_process_color = Colors.OKBLUE
        self.command_verbosity = 3

    @staticmethod
    def error_out(message=None):
        if message is not None:
            sys.exit(f"Error: {message}")
        else:
            sys.exit("The program encountered an error and has to exit.")

    @staticmethod
    def validate_file(the_file):
        return os.path.isfile(the_file)

    def validate_file_size(self, the_file):
        if not self.fake_run:
            return os.stat(the_file).st_size > 0
        else:
            return True

    def print_and_log(self, text, verbosity, color=Colors.ENDC):
        if verbosity <= self.verbosity:
            time_string = '[' + str(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")) + ']'
            print(time_string + ' ' + color + text + Colors.ENDC)

    def validate_file_and_size_or_error(self, the_file, error_prefix='The file',
                                        presence_suffix='doesn\'t exist',
                                        size_suffix='is size 0'):
        if not self.validate_file(the_file) and not self.fake_run:
            self.print_and_log(' '.join([error_prefix, the_file, presence_suffix]), 0, Colors.FAIL)
            self.error_out()

        if not self.validate_file_size(the_file) and not self.fake_run:
            self.print_and_log(' '.join([error_prefix, the_file, size_suffix]), 0, Colors.FAIL)
            self.error_out()

    def validate_haplotypes(self):
        if self.haplotype_file is None:
            self.errors += [self.sub_command + ' requires --haplotypes']
            return

        self.validate_file_and_size_or_error(self.haplotype_file, 'Haplotype file')

    def validate_reference_prefix(self):
        if self.reference_prefix is None:
            self.errors += [self.sub_command + ' requires --reference-prefix']
            return

        self.reference_tped_file = self.reference_prefix + '.tped'
        self.reference_tfam_file = self.reference_prefix + '.tfam'

        [self.validate_file_and_size_or_error(the_file=i) for i in [self.reference_tped_file, self.reference_tfam_file]]

    def validate_query_prefix(self):
        if self.query_prefix is None:
            self.errors += [self.sub_command + ' requires --query-prefix']

        self.query_tped_file = self.query_prefix + '.tped'
        self.query_tfam_file = self.query_prefix + '.tfam'

        [self.validate_file_and_size_or_error(the_file=i) for i in [self.query_tped_file, self.query_tfam_file]]

    def validate_reference_pickle(self):
        if self.reference_pickle_file is None:
            self.errors += [self.sub_command + ' requires --query-prefix']

        self.validate_file_and_size_or_error(self.reference_pickle_file, 'Reference pickle')

    def validate_options(self):
        if self.sub_command == 'build':
            self.validate_reference_prefix()
            self.validate_haplotypes()
            self.analysis += ['build_reference_set']

        if self.sub_command == 'query':
            self.validate_reference_pickle()
            self.validate_query_prefix()
            self.analysis += ['query_reference_set']

        if self.sub_command == 'coanc':
            self.validate_reference_prefix()
            if self.query_prefix is not None:
                self.validate_query_prefix()
            self.validate_haplotypes()
            self.analysis += ['build_coanc']

        if len(self.analysis) == 0:
            self.errors = self.errors + ['Nothing to do!']

    def go(self):

        print(self.analysis)
        while True:
            step = self.analysis[0]
            self.analysis = self.analysis[1:]
            function = getattr(self, step)
            function()

            if len(self.analysis) == 0:
                break

    def load_reference_pickle(self):
        print('Loading reference pickle ' + self.reference_pickle_file)
        pickle_file_handle = open(self.reference_pickle_file, 'rb')
        old_pickle = pickle.load(pickle_file_handle)
        self.reference_tfam = old_pickle.reference_tfam
        # self.reference_haplotype_counts = old_pickle.reference_haplotype_counts
        self.reference_haplotype_fractions = old_pickle.reference_haplotype_fractions
        self.reference_populations = old_pickle.reference_populations
        self.reference_background = old_pickle.reference_background
        self.haplotypes = old_pickle.haplotypes
        self.reference_copying_fractions = old_pickle.reference_copying_fractions
        print('Reference pickle file loaded!')
        del old_pickle

    def load_reference_tfam(self):
        self.reference_tfam = pd.read_table(self.reference_tfam_file, index_col=None, header=None, sep=' ')
        self.reference_tfam.index = self.reference_tfam.iloc[:, 1]
        self.reference_tfam.columns = ['population', 'id', 'x1', 'x2', 'x3', 'x4']
        self.reference_individuals = pd.Series(self.reference_tfam.index.values)
        print('Found', self.reference_tfam.shape[0], 'reference individuals')

        # Figure out the reference populations
        self.reference_populations = self.reference_tfam.iloc[:, 0].unique()
        print('Found', len(self.reference_populations), 'refeerence populations')
        self.reference_individual_populations = self.reference_tfam.iloc[:, 0]

        self.reference_population_counts = self.reference_tfam.iloc[:, 0].value_counts()
        self.reference_population_counts = self.reference_population_counts[self.reference_populations]

    def load_haplotypes(self):
        self.haplotypes = pd.read_table(self.haplotype_file, index_col=None, header=None)
        self.haplotypes.columns = ['chromosome', 'left', 'right', 'cM']
        self.haplotypes = self.haplotypes.loc[self.haplotypes['chromosome'].apply(lambda x: x in self.chromosomes), :]
        print('Found', self.haplotypes.shape[0], 'haplotypes\n')

    def load_query_tfam(self):
        self.query_tfam = pd.read_table(self.query_tfam_file, index_col=None, header=None, sep=' ')
        self.query_tfam.index = self.query_tfam.iloc[:, 1]
        self.query_tfam.columns = ['population', 'id', 'x1', 'x2', 'x3', 'x4']
        self.query_individuals = self.query_tfam.index.values
        print('Found', self.query_tfam.shape[0], 'query individuals')

    @staticmethod
    def pull_chromosome_tped(prefix, chromosome):

        # TODO: Create the files in a temporary directory
        # Pull the genotypes for a given chromosome from a given prefix
        chromosome_file_prefix = os.path.join('.', str(chromosome) + '.tmp')
        command = ' '.join(['plink',
                            '--tfile', prefix,
                            '--recode transpose',
                            '--chr', str(chromosome),
                            '--out', chromosome_file_prefix,
                            '--keep-allele-order'])
        # print(command)
        subprocess.call(command, shell=True)
        chromosome_tped_file = chromosome_file_prefix + '.tped'

        # Read the thing and remove the temporary files
        chromosome_tped = pd.read_table(chromosome_tped_file, index_col=None, header=None, sep=' ').iloc[:, 4::2]
        command = ' '.join(
            ['rm', ' '.join([chromosome_file_prefix + '.' + i for i in ['tped', 'tfam', 'nosex', 'log']])])
        # print(command)
        subprocess.call(command, shell=True)

        return chromosome_tped

    def load_reference_tpeds(self):

        # Read in all of the reference TPEDs
        # print(self.chromosomes)
        pool = mp.Pool(processes=self.threads)
        results = pool.map(lambda x: self.pull_chromosome_tped(self.reference_prefix, x), self.chromosomes)
        for i in range(len(self.chromosomes)):
            self.reference_tpeds[str(self.chromosomes[i])] = results[i]
            self.reference_tpeds[str(self.chromosomes[i])].columns = self.reference_individuals
        print('Loaded', len(self.reference_tpeds), 'refernece TPEDs')

    def load_query_tpeds(self):

        # Read in all of the reference TPEDs
        print(self.chromosomes)
        pool = mp.Pool(processes=self.threads)
        results = pool.map(lambda x: self.pull_chromosome_tped(self.query_prefix, x), self.chromosomes)
        for i in range(len(self.chromosomes)):
            self.query_tpeds[str(self.chromosomes[i])] = results[i]
            self.query_tpeds[str(self.chromosomes[i])].columns = self.query_individuals
        print('Loaded', len(self.query_tpeds), 'query TPEDs')

    def build_reference_set(self):

        # Load up the reference data
        self.load_haplotypes()
        self.load_reference_tfam()
        self.load_reference_tpeds()

        #        self.haplotypes = self.haplotypes.iloc[self.haplotypes.iloc[:,0].values in self.chromosomes, :]

        # Parallelize across equal-sized chunks of haplotypes for good speed
        chunk_size = min([30, math.ceil(self.haplotypes.shape[0] / self.threads)])
        print('Splitting haplotypes into chunks of', chunk_size)
        chunk_indices = list(range(0, math.ceil(self.haplotypes.shape[0] / chunk_size) * chunk_size + 1, chunk_size))
        # chunks = len(chunk_indices) - 1
        chunk_indices[-1] = self.haplotypes.shape[0]
        chunk_indices = [chunk_indices[i:i + 2] for i in range(len(chunk_indices) - 1)]

        # Find the haplotypes and corresponding genotypes for each chunk
        chunk_haplotypes = []
        chunk_genotypes = []
        for this_chunk in chunk_indices:

            chunk_haplotypes += [self.haplotypes.iloc[this_chunk[0]:this_chunk[1], :]]
            this_chunk_genotypes = []

            for this_haplotype in range(this_chunk[0], this_chunk[1]):
                this_chromosome, left_snp, right_snp, cM = self.haplotypes.iloc[this_haplotype, :]
                this_chromosome, left_snp, right_snp = [int(i) for i in [this_chromosome, left_snp, right_snp]]
                this_chunk_genotypes += [self.reference_tpeds[str(this_chromosome)].iloc[left_snp:right_snp, :]]
            chunk_genotypes += [this_chunk_genotypes]

        chunk_data = [[chunk_indices[i]] + [chunk_haplotypes[i]] + [chunk_genotypes[i]] + [i] for i in
                      range(len(chunk_indices))]

        del self.reference_tpeds

        pool = mp.Pool(processes=self.threads)

        # Gives us back a list of lists
        # Dim 1 - Chunk
        # Dim 2 - Results of a chunk
        results = pool.map(self.chunk_build_reference_copying_fractions, chunk_data)
        pool.close()

        # Find the population counts for each reference haplotype
        self.reference_haplotype_counts = results[0][0]
        for i in range(1, len(results)):
            self.reference_haplotype_counts = pd.concat([self.reference_haplotype_counts, results[i][0]], axis=0)
        print(self.reference_haplotype_counts.shape)

        # Find the population fractions for each reference haplotype
        self.reference_haplotype_fractions = results[0][1]
        for i in range(1, len(results)):
            self.reference_haplotype_fractions = pd.concat([self.reference_haplotype_fractions, results[i][1]], axis=0)

        # Find the reference indvidiual-refernce population copying rates
        # Here the copying fraction across all of the haplotypes but scaled 0-1
        self.reference_copying_fractions = results[0][2]
        for i in range(1, len(results)):
            self.reference_copying_fractions += results[i][2]
        self.reference_copying_fractions = self.reference_copying_fractions.T.apply(lambda x: x / x.sum()).T

        # Find the background individual-reference copying rates.
        # We define the background as the lowest copying rate seen for any reference individual for any reference population.
        # Here the min of reach columin in the copying fractions
        self.reference_background = self.reference_copying_fractions.min(axis=0)

        # Adjust the reference individual v. population copying fractions by the backgroud and save
        self.reference_copying_fractions -= self.reference_background
        self.reference_copying_fractions = self.reference_copying_fractions.T.apply(lambda x: x / x.sum()).T
        self.reference_copying_fractions.to_csv(self.reference_output_file, sep='\t')

        # Pickle this reference set
        pickle_file = open(self.reference_pickle_output_file, 'wb')
        pickle.dump(self, pickle_file)
        pickle_file.close()

    def chunk_build_reference_copying_fractions(self, chunk_data):
        chunk_indices, chunk_haplotypes, chunk_genotypes, chunk_number = chunk_data

        print('Starting chunk ' + str(chunk_number))
        chunk_haplotype_counts = pd.Series([], dtype=pd.StringDtype())
        chunk_haplotype_fractions = pd.Series([], dtype=pd.StringDtype())

        chunk_reference_individual_counts = pd.DataFrame(0,
                                                         index=self.reference_individuals,
                                                         columns=self.reference_populations)

        # print('Found', chunk_haplotypes.shape[0], 'haplotypes for this chunk')
        # print('Found', len(chunk_haplotypes), 'corresponding genotype sets')

        for haplotype in range(len(chunk_haplotypes)):

            chromosome, left_snp, right_snp, cM = chunk_haplotypes.iloc[haplotype, :]
            chromosome, left_snp, right_snp = [int(i) for i in [chromosome, left_snp, right_snp]]
            haplotype_index = str(chromosome) + ':' + str(left_snp) + "-" + str(right_snp)
            # print(haplotype_index)

            haplotype_genotypes = chunk_genotypes[haplotype]

            # Find haplotypes from those individuals with no missingness
            have_missing = haplotype_genotypes.apply(lambda x: (x != 0).all())
            haplotype_genotypes = haplotype_genotypes.iloc[:, have_missing.values]
            haplotype_strings = haplotype_genotypes.apply(lambda x: ''.join(x.astype(str)))

            # Count the number of times each haplotype is seen in each population
            haplotype_counts = pd.DataFrame(
                {'population': self.reference_tfam.loc[haplotype_strings.index, 'population'],
                 'haplotype': haplotype_strings})
            haplotype_counts = haplotype_counts.groupby(['population', 'haplotype']).size().reset_index()
            haplotype_counts = haplotype_counts.pivot(index='haplotype', columns='population').fillna(0)
            haplotype_counts.columns = haplotype_counts.columns.droplevel()

            # Find the counts of each haplotype in each population
            distinct_haplotypes = haplotype_strings.unique()
            haplotype_counts_with_zeroes = pd.DataFrame(0, index=distinct_haplotypes,
                                                        columns=self.reference_populations)
            haplotype_counts_with_zeroes.loc[haplotype_counts.index, haplotype_counts.columns] += haplotype_counts
            haplotype_counts = haplotype_counts_with_zeroes
            # print(haplotype_counts.shape)

            # Keep track of the counts for each haplotype on this chromosome
            chunk_haplotype_counts[haplotype_index] = haplotype_counts
            chunk_haplotype_fractions[haplotype_index] = haplotype_counts.T.apply(lambda x: x / x.sum()).T * cM

            # Look for where the reference individuals fall and tally up their copying counts
            for individual in self.reference_individuals:

                if individual not in haplotype_strings.index:
                    continue

                individual_haplotype = haplotype_strings[individual]
                individual_population = self.reference_tfam.loc[individual, 'population']

                # print(individual)

                # Figure out if we need to remove this asshole individual from the reference set,
                # i.e., don't allow a haplotype to copy from itself
                if individual_haplotype in haplotype_counts.index:

                    # Don't try to count this if the reference individual was the only instance of the haplotype
                    if haplotype_counts.loc[individual_haplotype, :].sum() == 1:
                        continue

                    # Otherwise, use this haplotype, but subtract 1 for this reference individuals count
                    haplotype_counts_to_use = haplotype_counts.copy().loc[individual_haplotype, :]
                    haplotype_counts_to_use[individual_population] -= 1
                    # print(haplotype_counts_to_use)

                    # And normalize to 1 with the remaining haplotypes and update the individual
                    haplotype_counts_to_use = haplotype_counts_to_use / haplotype_counts_to_use.sum() * cM
                    # print(haplotype_counts_to_use)
                    chunk_reference_individual_counts.loc[individual, :] += haplotype_counts_to_use

        return [chunk_haplotype_counts, chunk_haplotype_fractions, chunk_reference_individual_counts, chunk_number]

    def query_reference_set(self):

        self.load_query_tfam()
        self.load_query_tpeds()
        self.load_reference_pickle()

        # Parallelize across equal-sized chunks of haplotypes for good speed
        chunk_size = min([30, math.ceil(self.haplotypes.shape[0] / self.threads)])
        # print('Splitting haplotypes into chunks of', chunk_size)
        chunk_indices = list(range(0, math.ceil(self.haplotypes.shape[0] / chunk_size) * chunk_size + 1, chunk_size))
        # chunks = len(chunk_indices) - 1
        chunk_indices[-1] = self.haplotypes.shape[0]
        chunk_indices = [chunk_indices[i:i + 2] for i in range(len(chunk_indices) - 1)]

        # Find the haplotypes and corresponding genotypes for each chunk
        chunk_haplotypes = []
        chunk_reference_fractions = []
        chunk_genotypes = []
        for this_chunk in chunk_indices:

            chunk_haplotypes += [self.haplotypes.iloc[this_chunk[0]:this_chunk[1], :]]
            chunk_reference_fractions += [self.reference_haplotype_fractions[this_chunk[0]:this_chunk[1]]]

            this_chunk_genotypes = []
            for this_haplotype in range(this_chunk[0], this_chunk[1]):
                this_chromosome, left_snp, right_snp, cM = self.haplotypes.iloc[this_haplotype, :]
                this_chromosome, left_snp, right_snp = [int(i) for i in [this_chromosome, left_snp, right_snp]]
                this_chunk_genotypes += [self.query_tpeds[str(this_chromosome)].iloc[left_snp:right_snp, :]]
            chunk_genotypes += [this_chunk_genotypes]

        chunk_data = [
            [chunk_indices[i]] + [chunk_haplotypes[i]] + [chunk_genotypes[i]] + [chunk_reference_fractions[i]] + [i] for
            i in range(len(chunk_indices))]

        # Fuck the GIL
        del self.reference_haplotype_fractions
        del self.query_tpeds

        pool = mp.Pool(processes=self.threads)

        results = pool.map_async(self.chunk_query_reference_copying_fractions, chunk_data)
        results.wait()
        results = results.get()
        pool.close()
        pool.join()

        # Find the query indvidiual-refernce population copying rates
        self.query_copying_fractions = results[0]
        for i in range(1, len(results)):
            self.query_copying_fractions += results[i]
        self.query_copying_fractions = self.query_copying_fractions.T.apply(lambda x: x / x.sum()).T

        # Remove the background seen in the reference set, making sure they don't go below 0
        self.query_copying_fractions -= self.reference_background
        self.query_copying_fractions[self.query_copying_fractions < 0] = 0
        self.query_copying_fractions = self.query_copying_fractions.T.apply(lambda x: x / x.sum()).T
        self.query_copying_fractions.to_csv(self.query_output_file, sep='\t')

        self.combined_copying_fractions = pd.concat([self.reference_copying_fractions, self.query_copying_fractions],
                                                    axis=0)
        self.combined_copying_fractions.to_csv(self.combined_output_file, sep='\t')

    def chunk_query_reference_copying_fractions(self, chunk_data):

        chunk_indices, chunk_haplotypes, chunk_genotypes, chunk_reference_fractions, chunk_number = chunk_data

        # print('Starting chunk ' + str(chunk_number))
        # chunk_haplotype_counts = pd.Series()
        # chunk_haplotype_fractions = pd.Series()
        chunk_query_individual_fractions = pd.DataFrame(0,
                                                        index=self.query_individuals,
                                                        columns=self.reference_populations)

        # print('Found', chunk_haplotypes.shape[0], 'haplotypes for this chunk')
        # print('Found', len(chunk_haplotypes), 'corresponding genotype sets')
        # print('Found', len(chunk_reference_fractions), 'corresponding fraction sets')

        # For each haplotype, get all of the counts for each individual
        for haplotype in range(len(chunk_haplotypes)):

            chromosome, left_snp, right_snp, cM = chunk_haplotypes.iloc[haplotype, :]
            chromosome, left_snp, right_snp = [int(i) for i in [chromosome, left_snp, right_snp]]

            haplotype_index = str(chromosome) + ':' + str(left_snp) + '-' + str(right_snp)

            # print(haplotype_index)
            haplotype_fractions = chunk_reference_fractions[haplotype_index]
            haplotype_genotypes = chunk_genotypes[haplotype]

            # Find haplotypes from those individuals with no missingness
            haplotype_strings = haplotype_genotypes.apply(lambda x: ''.join(x.astype(str)))

            # Find the copying fraction for each individual
            for individual in haplotype_strings.index.values:
                individual_haplotype = haplotype_strings[individual]

                if individual_haplotype in haplotype_fractions.index:
                    chunk_query_individual_fractions.loc[individual, :] += haplotype_fractions.loc[individual_haplotype,
                                                                           :]

        return chunk_query_individual_fractions

    def build_coanc(self):

        # Load up the reference data
        self.load_haplotypes()
        self.load_reference_tfam()
        self.load_query_tfam()

        # Keep track of all of the results
        # reference_results = pd.Series()
        # query_results = pd.Series()

        pool = mp.Pool(processes=self.threads)

        # For memory GIL BS, run each chromosome on its own
        copying_counts = pd.DataFrame(0, index=self.reference_individuals, columns=self.reference_individuals)
        haplotype_counts = pd.Series([0] * len(self.reference_individuals), index=self.reference_individuals)

        for chromosome in self.chromosomes:

            # Just the haps for this chromosome
            chromosome_haplotypes = self.haplotypes.loc[self.haplotypes.iloc[:, 0] == chromosome, :]

            # Just the TPEDs for this chromosome
            reference_tped = self.pull_chromosome_tped(self.reference_prefix, chromosome)
            reference_tped.columns = self.reference_individuals
            query_tped = self.pull_chromosome_tped(self.query_prefix, chromosome)

            # Parallelize across equal-sized chunks of haplotypes for good speed
            chunk_size = min([30, math.ceil(chromosome_haplotypes.shape[0] / self.threads)])
            # print('Splitting haplotypes into chunks of', chunk_size)
            chunk_indices = list(
                range(0, math.ceil(chromosome_haplotypes.shape[0] / chunk_size) * chunk_size + 1, chunk_size))
            # chunks = len(chunk_indices) - 1
            chunk_indices[-1] = chromosome_haplotypes.shape[0]
            chunk_indices = [chunk_indices[i:i + 2] for i in range(len(chunk_indices) - 1)]

            # Find the haplotypes and corresponding genotypes for each chunk
            chunk_haplotypes = []
            chunk_reference_genotypes = []
            chunk_query_genotypes = []
            for this_chunk in chunk_indices:

                chunk_haplotypes += [chromosome_haplotypes.iloc[this_chunk[0]:this_chunk[1], :]]

                this_chunk_reference_genotypes = []
                this_chunk_query_genotypes = []

                for this_haplotype in range(this_chunk[0], this_chunk[1]):
                    this_chromosome, left_snp, right_snp = chromosome_haplotypes.iloc[this_haplotype, :]
                    this_chunk_reference_genotypes += [reference_tped.iloc[left_snp:right_snp, :]]
                    this_chunk_query_genotypes += [query_tped.iloc[left_snp:right_snp, :]]

                chunk_reference_genotypes += [this_chunk_reference_genotypes]
                chunk_query_genotypes += [this_chunk_query_genotypes]

            chunk_data = [[chunk_indices[i]] +
                          [chunk_haplotypes[i]] +
                          [chunk_reference_genotypes[i]] +
                          [chunk_query_genotypes[i]] + [i] for i in range(len(chunk_indices))]

            del reference_tped
            del query_tped

            results = pool.map_async(self.chunk_build_coanc, chunk_data)
            results.wait()
            results = results.get()
            for i in range(len(results)):
                copying_counts += results[i][0]
                haplotype_counts += results[i][1]

        pool.close()
        pool.join()
        copying_counts = copying_counts.T.apply(lambda x: x / haplotype_counts).T
        copying_counts.to_csv('aVA.tsv', sep='\t')

    def chunk_build_coanc(self, chunk_data):

        chunk_indices, chunk_haplotypes, chunk_reference_genotypes, chunk_query_genotypes, chunk_number = chunk_data
        # print('Coanc chunk', chunk_number)

        # chunk_haplotype_fractions = pd.DataFrame(0, index = self.reference_individuals, columns = self.reference_individuals)

        chunk_haplotype_fractions = pd.Series([], dtype=pd.StringDtype())
        for individual in self.reference_individuals:
            chunk_haplotype_fractions[individual] = pd.Series([0] * len(self.reference_individuals),
                                                              index=self.reference_individuals)

        chunk_individual_counts = pd.Series([0] * len(self.reference_individuals), index=self.reference_individuals)

        for haplotype in range(len(chunk_haplotypes)):

            chromosome, left_snp, right_snp = chunk_haplotypes.iloc[haplotype, :]
            haplotype_index = str(chromosome) + ':' + str(left_snp) + "-" + str(right_snp)
            # print(haplotype_index)

            haplotype_genotypes = chunk_reference_genotypes[haplotype]

            # Find haplotypes from those individuals with no missingness
            have_missing = haplotype_genotypes.apply(lambda x: (x != 0).all())
            haplotype_reference_individuals = self.reference_individuals.loc[have_missing.values]
            haplotype_genotypes = haplotype_genotypes.iloc[:, have_missing.values]
            haplotype_strings = haplotype_genotypes.apply(lambda x: ''.join(x.astype(str)))

            # Count the number of times each haplotype is seen in each individual
            haplotype_counts = pd.DataFrame({'individual': haplotype_reference_individuals,
                                             'haplotype': haplotype_strings.values})
            haplotype_counts = haplotype_counts.groupby(['individual', 'haplotype']).size().reset_index()
            # print(haplotype_counts.head())
            haplotype_counts = haplotype_counts.pivot(index='haplotype', columns='individual').fillna(0)
            haplotype_counts.columns = haplotype_counts.columns.droplevel()

            # Find the counts of each haplotype across all of the individuals, adding in 0's
            distinct_haplotypes = haplotype_strings.unique()
            haplotype_counts_with_zeroes = pd.DataFrame(0, index=distinct_haplotypes,
                                                        columns=self.reference_individuals)
            haplotype_counts_with_zeroes.loc[haplotype_counts.index, haplotype_counts.columns] += haplotype_counts
            haplotype_counts = haplotype_counts_with_zeroes

            # print(haplotype_counts.head())

            # Look for where the reference individuals fall and tally up their copying counts
            for individual in haplotype_reference_individuals:

                if individual not in haplotype_strings.index:
                    continue
                individual_haplotype = haplotype_strings[individual]

                if individual_haplotype in haplotype_counts.index:
                    chunk_haplotype_fractions[individual] = chunk_haplotype_fractions[individual] + \
                                                            haplotype_counts.loc[individual_haplotype, :]

                    chunk_individual_counts[individual] += 1
        return [pd.concat(chunk_haplotype_fractions.to_dict(), axis=1).T, chunk_individual_counts]


# TODO: Implement checks for dependency progams
def check_dependencies():
    pass


# TODO: Implement safe way of executing commands
def run_command(command):
    pass


# TODO: Move all the screen print statements to logger
def init_logger(log_file=None):
    if log_file is None:
        logging.basicConfig(level=logging.INFO, format='[%(asctime)s] %(message)s',
                            datefmt='%m/%d/%Y %I:%M:%S %p')
    else:
        logging.basicConfig(filename=log_file, filemode='w', level=logging.DEBUG, format='[%(asctime)s] %(message)s',
                            datefmt='%m/%d/%Y %I:%M:%S %p')


if __name__ == '__main__':
    parser = ArgumentParser(prog=PROGRAM_NAME,
                            add_help=False,
                            description=f'''
                            {PROGRAM_NAME} - population scale haplotype copying
                            ''',
                            formatter_class=lambda prog: HelpFormatter(prog, width=120, max_help_position=120))

    parser.add_argument('--help', '-h', '--h', action='store_true', default=False)

    subparsers = parser.add_subparsers(title=f"{PROGRAM_NAME} commands")

    # Make the build sub-command parser
    build_parser = subparsers.add_parser('build', help='Build reference set',
                                         formatter_class=lambda prog: HelpFormatter(prog, width=120,
                                                                                    max_help_position=120))
    build_parser.set_defaults(sub_command='build')

    # Add build input arguments
    build_input_group = build_parser.add_argument_group('Input options')
    build_input_group.add_argument('--haplotypes', required=False, default=None, metavar='<TSV>',
                                   help='File of haplotype positions')
    build_input_group.add_argument('--reference-prefix', required=False, default=None, metavar='<PREFIX>',
                                   help='Prefix for the reference TPED/TFAM input files')

    # Add build output arguments
    build_output_group = build_parser.add_argument_group('Output options')
    build_output_group.add_argument('--reference-pickle-out', required=False, default=None, metavar='<OUTPUT>',
                                    help='The reference pickle file output')
    build_output_group.add_argument('--reference-out', required=False, default=None, metavar='<OUTPUT>',
                                    help='The reference copying matrix output')

    build_running_group = build_parser.add_argument_group('Running options')
    build_running_group.add_argument('--threads', required=False, type=int, default=1, metavar='<INT>',
                                     help='How many parallel threads to run')
    build_running_group.add_argument('--per-individual', required=False, default=False, action='store_true',
                                     help='Generate per-individual copying rather than per-population copying')

    # Make the query sub-command
    query_parser = subparsers.add_parser('query', help='Query reference set',
                                         formatter_class=lambda prog: HelpFormatter(prog, width=120,
                                                                                    max_help_position=120))
    query_parser.set_defaults(sub_command='query')

    query_input_group = query_parser.add_argument_group('Input options')
    query_input_group.add_argument('--reference-pickle', required=False, default=None, metavar='<PICKLE>',
                                   help='A pre-made reference pickle')
    query_input_group.add_argument('--query-prefix', required=False, default=None, metavar='<PREFIX>',
                                   help='Prefix for the query TPED/TFAM input files')

    query_output_group = query_parser.add_argument_group('Output options')
    query_output_group.add_argument('--query-out', required=False, default=None, metavar='<OUTPUT>',
                                    help='The query copying matrix output')
    query_output_group.add_argument('--combined-out', required=False, default=None, metavar='<OUTPUT>',
                                    help='The combined reference/query copying matrix output')

    query_running_group = query_parser.add_argument_group('Running options')
    query_running_group.add_argument('--threads', required=False, type=int, default=1, metavar='<INT>',
                                     help='How many parallel threads to run')

    # Make the co-ancestry sub-command
    coanc_parser = subparsers.add_parser('coanc', help='Individual v. individual co-ancestry',
                                         formatter_class=lambda prog: HelpFormatter(prog, width=120,
                                                                                    max_help_position=120))
    coanc_parser.set_defaults(sub_command='coanc')

    coanc_input_group = coanc_parser.add_argument_group('Input options')
    coanc_input_group.add_argument('--haplotypes', required=False, default=None, metavar='<TSV>',
                                   help='File of haplotype positions')
    coanc_input_group.add_argument('--reference-prefix', required=False, default=None, metavar='<PICKLE>',
                                   help='Prefix for the reference TPED/TFAM input files')
    coanc_input_group.add_argument('--query-prefix', required=False, default=None, metavar='<PREFIX>',
                                   help='Prefix for the query TPED/TFAM input files')

    coanc_output_group = coanc_parser.add_argument_group('Output options')
    coanc_output_group.add_argument('--reference-out', required=False, default=None, metavar='<OUTPUT>',
                                    help='The reference v. reference copying matrix output')
    coanc_output_group.add_argument('--query-out', required=False, default=None, metavar='<OUTPUT>',
                                    help='The query v. reference copying matrix output')
    coanc_output_group.add_argument('--combined-out', required=False, default=None, metavar='<OUTPUT>',
                                    help='The all v. reference copying matrix output')

    coanc_running_group = coanc_parser.add_argument_group('Running options')
    coanc_running_group.add_argument('--threads', required=False, type=int, default=1, metavar='<INT>',
                                     help='How many parallel threads to run')

    # Make the haplotype maker command
    coanc_parser = subparsers.add_parser('hapmake', help='Create haplotypes',
                                         formatter_class=lambda prog: HelpFormatter(prog, width=120,
                                                                                    max_help_position=120))
    coanc_parser.set_defaults(sub_command='hapmake')

    coanc_input_group = coanc_parser.add_argument_group('Input options')
    coanc_input_group.add_argument('--map-file', required=False, default=None, metavar='',
                                   help='File of haplotype positions')

    options, unknown_arguments = parser.parse_known_args()

    if options.help or "sub_command" not in options:
        print(Colors.HEADER)
        parser.print_help()
        print(Colors.ENDC)
        sys.exit(0)

    if options.sub_command == 'build' and options.help:
        print(Colors.HEADER)
        build_parser.print_help()
        print(Colors.ENDC)
        sys.exit(0)

    if options.sub_command == 'query' and options.help:
        print(Colors.HEADER)
        query_parser.print_help()
        print(Colors.ENDC)
        sys.exit(0)

    if options.sub_command == 'coanc' and options.help:
        print(Colors.HEADER)
        coanc_parser.print_help()
        print(Colors.ENDC)
        sys.exit(0)

    # Start the analysis
    analysis = Analysis(options, unknown_arguments)

    analysis.validate_options()

    if len(analysis.errors) > 0:
        print(Colors.HEADER)
        if options.sub_command == 'build':
            build_parser.print_help()
        elif options.sub_command == 'query':
            query_parser.print_help()
        elif options.sub_command == 'coanc':
            coanc_parser.print_help()
        print(Colors.ENDC)
        print(Colors.FAIL + '\n\nErrors:')
        [print(i) for i in analysis.errors]
        print(Colors.ENDC)
        sys.exit()

    # If we're still good, start the actual analysis
    analysis.go()
