#!/usr/bin/env python3
"""Pastrami - Population scale haplotype copying script"""

__author__ = "Andrew Conley, Lavanya Rishishwar"
__copyright__ = "Copyright 2021, Andrew Conley, Lavanya Rishishwar"
__credits__ = ["Andrew Conely", "Lavanya Rishishwar", "Shivam Sharma", "Emily Norris", "Sonali Gupta"]
__license__ = "GPL"
__version__ = "0.3"
__maintainer__ = "Andrew Conley, Lavanya Rishishwar"
__email__ = "aconley@ihrc.com; lrishishwar@ihrc.com"
__status__ = "Development"
__title__ = "pastrami.py"

# Standard modules
import logging
import math
import os.path
import pickle
import random
import re
import shlex
import shutil
import statistics
import string
import subprocess
import sys
from argparse import ArgumentParser, HelpFormatter

py_version = sys.version_info
if py_version[0] < 3 or py_version[1] < 4:
    sys.exit(f"Error: {__title__} requires Python version 3.4+ to work. Please install a newer version of Python.")

# Additional installs
try:
    import numpy as np
except ModuleNotFoundError as err:
    sys.exit(f"Error: Numpy not found. Please install numpy prior to running {__title__}")
try:
    from scipy.optimize import minimize
except ModuleNotFoundError as err:
    sys.exit(f"Error: Scipy not found. Please install scipy prior to running {__title__}")
try:
    import pandas as pd
except ModuleNotFoundError as err:
    sys.exit(f"Error: Pandas not found. Please install pandas prior to running {__title__}")
try:
    import pathos.multiprocessing as mp
except ModuleNotFoundError as err:
    sys.exit(f"Error: Pathos not found. Please install pathos prior to running {__title__}")

VERSION = __version__
PROGRAM_NAME = __title__


class Colors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'


class Support:
    @staticmethod
    def error_out(message: str = None):
        if message is not None:
            sys.exit(Colors.FAIL + f"Error: {message}" + Colors.ENDC)
        else:
            sys.exit(Colors.FAIL + "The program encountered an error and has to exit." + Colors.ENDC)

    @staticmethod
    def validate_file(the_file: str):
        return os.path.isfile(the_file)

    @staticmethod
    def validate_file_size(the_file: str, fake_run: str = False):
        if fake_run:
            return True
        else:
            return os.stat(the_file).st_size > 0

    @staticmethod
    def validate_file_and_size_or_error(the_file: str, error_prefix: str = 'The file',
                                        presence_suffix: str = 'doesn\'t exist',
                                        size_suffix: str = 'is size 0', fake_run: bool = False):
        if not Support.validate_file(the_file=the_file) and not fake_run:
            error_message = ' '.join([error_prefix, the_file, presence_suffix])
            Support.error_out(message=error_message)

        if not Support.validate_file_size(the_file=the_file) and not fake_run:
            error_message = ' '.join([error_prefix, the_file, size_suffix])
            Support.error_out(message=error_message)

    @staticmethod
    def validate_dir(the_dir: str):
        return os.path.isdir(the_dir)

    # TODO: This is a hastily written method, needs error fixing
    @staticmethod
    def validate_dir_or_error(the_dir: str, error_prefix: str = "The dir", presence_suffix: str = "doesn't exist",
                              fake_run: bool = False):
        if not Support.validate_dir(the_dir=the_dir) and not fake_run:
            error_message = ' '.join([error_prefix, the_dir, presence_suffix])
            Support.error_out(message=error_message)

    # TODO: Implement checks for dependency progams
    @staticmethod
    def check_dependencies(program_list: list = None) -> list:
        errors = []
        for program in program_list:
            if shutil.which(program) is None:
                errors.append(program)
        return errors

    @staticmethod
    def find_plink_binary():
        if shutil.which("plink") is not None:
            return "plink"
        elif shutil.which("plink2") is not None:
            return "plink2"
        else:
            return None

    @staticmethod
    def run_command(command_str: str = None, command_list: list = None, shell=False):
        if command_str is None and command_list is None:
            raise ValueError("Support.run_command() was called without any command to execute.")
        try:
            if command_str is not None:
                logging.info(f"Attempting to run: {command_str}")
                output = subprocess.check_output(shlex.split(command_str), encoding="utf-8", shell=shell)

            else:
                logging.info(f"Attempting to run: " + " ".join([str(x) for x in command_list]))
                output = subprocess.check_output(command_list, encoding="utf-8", shell=shell)
        except subprocess.CalledProcessError as e:
            logging.error(f"Encountered an error executing the command: ")
            if command_str is not None:
                logging.error(command_str)
            else:
                logging.error(command_list)
            logging.error(f"Error details:")
            logging.error(f"Exit code={e.returncode}")
            logging.error(f"Error message={e.output}")
            sys.exit(1)
        # logging.info(f"Command output = {output}")
        logging.info("Command executed without raising any exceptions")
        return output

    @staticmethod
    def validate_filename(filename: str):
        if re.match(r"^[a-zA-Z0-9_.-]+$", filename):
            return True
        else:
            return False

    @staticmethod
    def validate_output_prefix(out_prefix: str):
        parent, prefix = os.path.split(out_prefix)
        if parent != "":
            if not Support.validate_dir(parent):
                Support.safe_dir_create(parent)
        return Support.validate_filename(prefix)

    @staticmethod
    def safe_dir_create(this_dir: str):
        try:
            os.makedirs(this_dir)
        except IOError:
            print(f"I don't seem to have access to output prefix directory. Are the permissions correct?")
            sys.exit(1)

    @staticmethod
    def safe_dir_rm(this_dir: str):
        try:
            os.rmdir(this_dir)
        except IOError:
            print(f"I don't seem to have access to output prefix directory. Are the permissions correct?")
            sys.exit(1)

    @staticmethod
    def merge_fam_files(infile1: str, infile2: str, outputfile: str):
        """Merged two TFAM files into a single one (for aggregate function)

        Parameters
        ----------
        infile1 : str
            Input TFAM file #1 (e.g., reference TFAM file)

        infile2: str
            Input TFAM file #2 (e.g., query TFAM file)

        outputfile: str
            Output TFAM file

        Returns
        -------
        None
        """
        with open(outputfile, "w") as out_handle:
            with open(infile1, "r") as infile1_handle:
                for line in infile1_handle:
                    out_handle.write(line)
            with open(infile2, "r") as infile2_handle:
                for line in infile2_handle:
                    out_handle.write(line)

    @staticmethod
    def create_pop_group_from_tfam(tfam_in_file: str, tsv_out_file: str):
        """Takes unique population names from the input file and make them as group

        Parameters
        ----------
        tfam_in_file : str
            Input TFAM file

        tsv_out_file: str
            Output TSV file

        Returns
        -------
        None
        """
        unique_populations = {}
        with open(tfam_in_file, "r") as f_in:
            for line in f_in:
                pop = line.strip().split()[0]
                unique_populations[pop] = True

        unique_populations = sorted(unique_populations.keys())
        with open(tsv_out_file, "r") as f_out:
            f_out.write("#Population\tGroup\n")
            for this_pop in unique_populations:
                f_out.write(f"{this_pop}\t{this_pop}\n")

    @staticmethod
    def init_logger(log_file, verbosity):
        """Configures the logging for printing
        Returns
        -------
        None
            Logger behavior is set based on the Inputs variable
        """
        try:
            logging.basicConfig(filename=log_file, filemode="w", level=logging.DEBUG,
                                format=f"[%(asctime)s] %(message)s",
                                datefmt="%m-%d-%Y %I:%M:%S %p")
        except FileNotFoundError:
            print(f"The supplied location for the log file '{log_file}'" +
                  f"doesn't exist. Please check if the location exists.")
            sys.exit(1)
        except IOError:
            print(f"I don't seem to have access to make the log file." +
                  f"Are the permissions correct or is there a directory with the same name?")
            sys.exit(1)

        if verbosity:
            console = logging.StreamHandler()
            console.setLevel(logging.INFO)
            formatter = logging.Formatter(fmt=f"[%(asctime)s] %(message)s", datefmt="%m-%d-%Y %I:%M:%S %p")
            console.setFormatter(formatter)
            logging.getLogger().addHandler(console)


class Analysis:
    chromosomes = list(range(1, 23))
    fake_run = False
    debug = True
    min_haplotype_occurences = 0
    optim_step_size = 0.0001
    error_threshold = 1e-8
    optim_iterations = 10
    tolerance = 1e-8
    ancestry_fraction_postfix = "_fractions.Q"
    ancestry_painting_postfix = "_paintings.Q"
    pop_estimates_postfix = "_estimates.Q"
    outsourced_optimizer_pop_estimates_postfix = "_outsourcedOptimizer_estimates.Q"

    finegrain_estimates_postfix = "_fine_grain_estimates.Q"
    program_list = {'required': ["plink"]}

    def __init__(self, opts):
        # General attributes
        self.threads = opts.threads
        self.log_file = opts.log_file
        self.verbosity = opts.verbosity
        self.pool = None

        # Any errors we encounter
        self.errors = []
        # The actual queue for the analysis
        self.analysis = []

        # Verbosity levels and colors
        self.error_color = Colors.FAIL
        self.main_process_verbosity = 1
        self.warning_color = Colors.WARNING
        self.warning_verbosity = 1
        self.main_process_color = Colors.OKGREEN
        self.sub_process_verbosity = 2
        self.sub_process_color = Colors.OKBLUE
        self.command_verbosity = 3

        # Sub-commands
        self.sub_command = opts.sub_command
        self.plink_command = None

        self.ancestry_infile = None
        self.combined_copying_fractions = None
        self.combined_output_file = None
        self.fam_infile = None
        self.haplotype_file = None
        self.haplotypes = None
        self.map_dir = None
        self.max_rate = None
        self.max_snps = None
        self.min_snps = None
        self.out_prefix = None
        self.pop_group_file = None
        self.query_combined_file = None
        self.query_copying_fractions = None
        self.query_output_file = None
        self.query_prefix = None
        self.query_tfam = None
        self.query_tfam_file = None
        self.query_tped_file = None
        self.reference_copying_fractions = None
        self.reference_haplotype_counts = None
        self.reference_haplotype_fractions = None
        self.reference_individual_populations = None
        self.reference_individuals = None
        self.reference_output_file = None
        self.reference_pickle_output_file = None
        self.reference_population_counts = None
        self.reference_populations = None
        self.reference_prefix = None
        self.reference_tfam = None
        self.reference_tfam_file = None
        self.reference_tped_file = None

        # All sub-command options
        if self.sub_command == 'all':
            # TODO: Test all subcommand
            # Required options
            self.reference_prefix = opts.reference_prefix
            self.query_prefix = opts.query_prefix
            self.out_prefix = opts.out_prefix
            self.map_dir = opts.map_dir
            self.haplotype_file = opts.haplotypes
            self.pop_group_file = opts.pop_group_file
            # Outputs to be made
            self.reference_pickle_output_file = None
            self.reference_output_file = None
            self.reference_pickle_file = None
            self.query_output_file = None
            self.combined_output_file = None
            self.ancestry_infile = None
            self.fam_infile = None
            # Hapmake
            self.min_snps = opts.min_snps
            self.max_snps = opts.max_snps
            self.max_rate = opts.max_rate
            # Build options
            self.reference_tpeds = pd.Series([], dtype=pd.StringDtype())
            self.reference_background = pd.Series([], dtype=pd.StringDtype())
            # Query options
            self.query_tpeds = pd.Series([], dtype=pd.StringDtype())
            self.query_individuals = None
            # aggregate options
            self.ref_pop_group_map = {}
            self.ind_pop_map = {}
            self.pop_ind_map = {}  # reverse map of ind_pop_map, sacrificing memory for speed later on
            self.reference_individual_dict = {}
            self.reference_pops = {}
            self.ancestry_fractions = {}
            self.af_header = []
            self.painting_vectors = {}
            self.painting_vectors_keys = []
            self.fine_grain_estimates = {}

        # Hapmake sub-command options
        if self.sub_command == 'hapmake':
            self.min_snps = opts.min_snps
            self.max_snps = opts.max_snps
            self.max_rate = opts.max_rate
            self.map_dir = opts.map_dir
            self.haplotype_file = opts.haplotypes

        # Build sub-command options
        if self.sub_command == 'build':
            self.reference_pickle_output_file = opts.reference_pickle_out
            self.reference_prefix = opts.reference_prefix
            self.haplotype_file = opts.haplotypes
            self.reference_output_file = opts.reference_out
            self.reference_tpeds = pd.Series([], dtype=pd.StringDtype())
            self.reference_background = pd.Series([], dtype=pd.StringDtype())

        # Query sub-command options
        if self.sub_command == 'query':
            self.reference_pickle_file = opts.reference_pickle
            self.query_prefix = opts.query_prefix
            self.query_output_file = opts.query_out
            self.combined_output_file = opts.combined_out
            self.query_tpeds = pd.Series([], dtype=pd.StringDtype())
            self.query_individuals = None

        # Co-ancestry sub-command options
        if self.sub_command == 'coanc':
            self.haplotype_file = opts.haplotypes
            self.reference_prefix = opts.reference_prefix
            self.query_prefix = opts.query_prefix
            self.reference_output_file = opts.reference_out
            self.query_output_file = opts.query_out
            self.query_combined_file = opts.combined_out

        # Aggregate sub-command options
        if self.sub_command == 'aggregate':
            self.ancestry_infile = opts.ancestry_infile
            self.pop_group_file = opts.pop_group_file
            self.out_prefix = opts.out_prefix
            self.fam_infile = opts.fam_infile
            self.ref_pop_group_map = {}
            self.ind_pop_map = {}
            self.pop_ind_map = {}  # reverse map of ind_pop_map, sacrificing memory for speed later on
            self.reference_individual_dict = {}
            self.reference_pops = {}
            self.ancestry_fractions = {}
            self.af_header = []
            self.painting_vectors = {}
            self.painting_vectors_keys = []
            self.fine_grain_estimates = {}

        Support.init_logger(log_file=self.log_file, verbosity=self.verbosity)

    """ 
        [Class section] Run control
    """

    def validate_options(self):

        plink_command = Support.find_plink_binary()
        if plink_command is None:
            self.errors += ["Can't find plink or plink2. Please make sure the binary exists as one of those two names"]
        else:
            self.plink_command = plink_command

        if self.sub_command == "all":
            self.validate_and_set_all_subcommand_options()
            # self.analysis += ['build_reference_set', 'query_reference_set', 'build_coanc', 'post_pastrami']
            self.analysis += ['build_reference_set', 'query_reference_set', 'post_pastrami']

        if self.sub_command == 'hapmake':
            self.validate_map_dir()
            self.analysis += ['make_haplotypes']

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

        if self.sub_command == 'aggregate':
            self.validate_ancestry_infile()
            self.validate_pop_group_file()
            self.validate_fam_infile()
            self.analysis += ['post_pastrami']

        if len(self.analysis) == 0:
            self.errors = self.errors + ['Nothing to do!']

    def __str__(self):
        long_string = f"""
        Class constants: 
                chromosomes                   =  {Analysis.chromosomes}
                fake_run                      =  {Analysis.fake_run}
                debug                         =  {Analysis.debug}
                min_haplotype_occurences      =  {Analysis.min_haplotype_occurences}
                optim_step_size               =  {Analysis.optim_step_size}
                error_threshold               =  {Analysis.error_threshold}
                optim_iterations              =  {Analysis.optim_iterations}
                tolerance                     =  {Analysis.tolerance}
                ancestry_fraction_postfix     =  {Analysis.ancestry_fraction_postfix}
                ancestry_painting_postfix     =  {Analysis.ancestry_painting_postfix}
                pop_estimates_postfix         =  {Analysis.pop_estimates_postfix}
                finegrain_estimates_postfix   =  {Analysis.finegrain_estimates_postfix}
                
        Instance variables:
                * General program parameters
                log_file                      =  {self.log_file}
                threads                       =  {self.threads}
                verbosity                     =  {self.verbosity}
                * Verbosity options
                command_verbosity             =  {self.command_verbosity}
                main_process_verbosity        =  {self.main_process_verbosity}
                sub_process_verbosity         =  {self.sub_process_verbosity}
                warning_verbosity             =  {self.warning_verbosity}
                * Subcommand to be executed
                sub_command                   =  {self.sub_command}
                * Hapmake-specific parameter options
                max_rate                      =  {self.max_rate}
                max_snps                      =  {self.max_snps}
                min_snps                      =  {self.min_snps}
                * Hapmake-specific input/output options
                out_prefix                    =  {self.out_prefix}
                map_dir                       =  {self.map_dir}
                haplotype_file                =  {self.haplotype_file}
                * Query files input/output options
                query_prefix                  =  {self.query_prefix}
                query_tfam_file               =  {self.query_tfam_file}
                query_tped_file               =  {self.query_tped_file}
                query_output_file             =  {self.query_output_file}
                query_combined_file           =  {self.query_combined_file}
                * Reference files input/output options
                reference_prefix              =  {self.reference_prefix}
                reference_tfam_file           =  {self.reference_tfam_file}
                reference_tped_file           =  {self.reference_tped_file}
                reference_output_file         =  {self.reference_output_file}
                reference_pickle_output_file  =  {self.reference_pickle_output_file}
                * Combined query-reference file location
                combined_output_file          =  {self.combined_output_file}
                * Aggregate-specific options
                pop_group_file                =  {self.pop_group_file}
                ancestry_infile               =  {self.ancestry_infile}
                fam_infile                    =  {self.fam_infile}
        """
        return long_string

    # TODO: Print a summary of what parameters were provided, what needs to be performed
    def summarize_run(self):
        logging.info(self.main_process_color + str(self) + Colors.ENDC)
        logging.info(self.main_process_color + f"Analysis to perform: " + ",".join(self.analysis) + Colors.ENDC)

    def go(self):
        self.summarize_run()
        # self.pool = mp.Pool(processes=self.threads)
        self.pool = mp.ProcessingPool(nodes=self.threads)
        while True:
            step = self.analysis[0]
            self.analysis = self.analysis[1:]
            function = getattr(self, step)
            function()
            if len(self.analysis) == 0:
                break
        # self.pool.terminate()

    """ 
        [Class section] Functions for validating file
    """

    def validate_and_set_all_subcommand_options(self):
        self.validate_reference_prefix()
        self.validate_query_prefix()

        Support.validate_output_prefix(self.out_prefix)

        if self.haplotype_file is None:
            if self.map_dir is None:
                self.errors += [self.sub_command + ' requires --haplotypes or --map-dir']
                return
            else:
                self.validate_map_dir()
                self.haplotype_file = self.out_prefix + ".hap"

        if self.pop_group_file is None:
            self.pop_group_file = self.out_prefix + ".pop_group.tsv"
            Support.create_pop_group_from_tfam(tfam_in_file=self.reference_tfam_file, tsv_out_file=self.pop_group_file)
        else:
            self.validate_pop_group_file()

        self.reference_pickle_output_file = self.out_prefix + ".pickle"
        self.reference_output_file = self.out_prefix + ".hap"
        self.reference_pickle_file = self.reference_pickle_output_file
        self.query_output_file = self.out_prefix + "_query.tsv"
        self.combined_output_file = self.out_prefix + ".tsv"
        self.ancestry_infile = self.combined_output_file
        self.fam_infile = self.out_prefix + ".fam"
        Support.merge_fam_files(infile1=self.query_tfam_file,
                                infile2=self.reference_tfam_file,
                                outputfile=self.fam_infile)

    def validate_haplotypes(self):
        if self.haplotype_file is None:
            self.errors += [self.sub_command + ' requires --haplotypes']
            return

        Support.validate_file_and_size_or_error(the_file=self.haplotype_file, error_prefix='Haplotype file',
                                                fake_run=self.fake_run)

    def validate_reference_prefix(self):
        if self.reference_prefix is None:
            self.errors += [self.sub_command + ' requires --reference-prefix']
            return

        self.reference_tped_file = self.reference_prefix + '.tped'
        self.reference_tfam_file = self.reference_prefix + '.tfam'

        for i in [self.reference_tped_file, self.reference_tfam_file]:
            Support.validate_file_and_size_or_error(the_file=i, fake_run=self.fake_run)

    def validate_query_prefix(self):
        if self.query_prefix is None:
            self.errors += [self.sub_command + ' requires --query-prefix']
            return

        self.query_tped_file = self.query_prefix + '.tped'
        self.query_tfam_file = self.query_prefix + '.tfam'

        for i in [self.query_tped_file, self.query_tfam_file]:
            Support.validate_file_and_size_or_error(the_file=i, fake_run=self.fake_run)

    def validate_reference_pickle(self):
        if self.reference_pickle_file is None:
            self.errors += [self.sub_command + ' requires --query-prefix']
            return

        Support.validate_file_and_size_or_error(the_file=self.reference_pickle_file,
                                                error_prefix='Reference pickle', fake_run=self.fake_run)

    def validate_map_dir(self):
        if self.map_dir is None:
            self.errors += [self.sub_command + ' requires --map-dir']
            return

        Support.validate_dir_or_error(the_dir=self.map_dir, error_prefix='Map directory', fake_run=self.fake_run)

    def validate_ancestry_infile(self):
        if self.ancestry_infile is None:
            self.errors += [self.sub_command + ' requires --pastrami-output']
            return

        Support.validate_file_and_size_or_error(the_file=self.ancestry_infile,
                                                error_prefix='Pastrami\' query output',
                                                fake_run=self.fake_run)

    # TODO: If user doesn't supply pop-group file, create one based on the TFAM file
    def validate_pop_group_file(self):
        if self.pop_group_file is None:
            self.errors += [self.sub_command + ' requires --pop-group']
            return

        Support.validate_file_and_size_or_error(the_file=self.pop_group_file,
                                                error_prefix='Population group mapping file',
                                                fake_run=self.fake_run)

    def validate_fam_infile(self):
        if self.fam_infile is None:
            self.errors += [self.sub_command + ' requires --pastrami-fam']
            return

        Support.validate_file_and_size_or_error(the_file=self.fam_infile,
                                                error_prefix='FAM input',
                                                fake_run=self.fake_run)

    """ 
        [Class section] Haplotype maker
    """

    def process_hapmap_file(self, chrom: int):
        logging.info(f"[hapmake|chr{chrom}] Started processing")
        haplotypes = ""
        map_data = []
        with open(os.path.join(self.map_dir, f"chr{chrom}.map"), "r") as f:
            for line in f:
                (position, cmorgan, snp) = line.rstrip().split("\t")
                map_data.append([int(position), float(cmorgan), snp])
        logging.info(f"[hapmake|chr{chrom}] File read")

        left_snp = 0
        right_snp = 0
        snps = False

        for row in range(len(map_data)):
            right_snp += 1
            if right_snp >= len(map_data):
                break

            # If the two SNPs have a recombination rate great than the max rate
            if map_data[right_snp][1] - map_data[left_snp][1] >= self.max_rate:
                if right_snp - left_snp >= self.min_snps:
                    snps = True
                else:
                    left_snp = right_snp
                    right_snp += 1

            # If the haplotype is long enough
            if right_snp - left_snp >= self.max_snps:
                snps = True

            # If snps isn't False, then save the range of the window
            if snps is True:
                haplotypes += f"{chrom}\t{left_snp}\t{right_snp}\t{map_data[right_snp - 2][1] - map_data[left_snp][1]}\n"
                snps = False
                left_snp = right_snp
        logging.info(f"[hapmake|chr{chrom}] All haplotypes discovered")
        return haplotypes

    def make_haplotypes(self):
        logging.info(f"[hapmake|chrAll] Starting pool for processing")
        # pool = mp.Pool(processes=self.threads)
        # self.pool.restart()
        results = self.pool.map(self.process_hapmap_file, range(1, 23))
        # results.wait()
        # results = results.get()
        with open(self.haplotype_file, "w") as f:
            f.write("".join(results) + "\n")
        logging.info(f"[hapmake|chrAll] Files written")
        # self.pool.close()
        # self.pool.join()
        # self.pool.terminate()

    """ 
        [Class section] Core Pastrami code - to be fragmented further in future
    """

    def load_reference_pickle(self):
        logging.info('Loading reference pickle ' + self.reference_pickle_file)
        pickle_file_handle = open(self.reference_pickle_file, 'rb')
        old_pickle = pickle.load(pickle_file_handle)
        self.reference_tfam = old_pickle.reference_tfam
        # self.reference_haplotype_counts = old_pickle.reference_haplotype_counts
        self.reference_haplotype_fractions = old_pickle.reference_haplotype_fractions
        self.reference_populations = old_pickle.reference_populations
        self.reference_background = old_pickle.reference_background
        self.haplotypes = old_pickle.haplotypes
        self.reference_copying_fractions = old_pickle.reference_copying_fractions
        del old_pickle
        logging.info('Reference pickle file loaded!')

    def load_reference_tfam(self):
        self.reference_tfam = pd.read_table(self.reference_tfam_file, index_col=None, header=None, sep=' ')
        self.reference_tfam.index = self.reference_tfam.iloc[:, 1]
        self.reference_tfam.columns = ['population', 'id', 'x1', 'x2', 'x3', 'x4']
        self.reference_individuals = pd.Series(self.reference_tfam.index.values)
        logging.info(f"Found {self.reference_tfam.shape[0]} reference individuals")

        # Figure out the reference populations
        self.reference_populations = self.reference_tfam.iloc[:, 0].unique()
        logging.info(f"Found {len(self.reference_populations)} refeerence populations")
        self.reference_individual_populations = self.reference_tfam.iloc[:, 0]

        self.reference_population_counts = self.reference_tfam.iloc[:, 0].value_counts()
        self.reference_population_counts = self.reference_population_counts[self.reference_populations]

    def load_haplotypes(self):
        self.haplotypes = pd.read_table(self.haplotype_file, index_col=None, header=None)
        self.haplotypes.columns = ['chromosome', 'left', 'right', 'cM']
        self.haplotypes = self.haplotypes.loc[self.haplotypes['chromosome'].apply(lambda x: x in self.chromosomes), :]
        logging.info(f"Found {self.haplotypes.shape[0]} haplotypes")

    def load_query_tfam(self):
        self.query_tfam = pd.read_table(self.query_tfam_file, index_col=None, header=None, sep=' ')
        self.query_tfam.index = self.query_tfam.iloc[:, 1]
        self.query_tfam.columns = ['population', 'id', 'x1', 'x2', 'x3', 'x4']
        self.query_individuals = self.query_tfam.index.values
        logging.info(f"Found {self.query_tfam.shape[0]} query individuals")

    def pull_chromosome_tped(self, prefix, chromosome, temp_dir):
        # Pull the genotypes for a given chromosome from a given prefix
        chromosome_file_prefix = os.path.join(temp_dir, str(chromosome) + '.tmp')
        logging.info(f"Loading SNPs from chr{chromosome}")
        command = [self.plink_command, '--tfile', prefix,
                   '--recode transpose',
                   '--chr', str(chromosome),
                   '--out', chromosome_file_prefix,
                   '--keep-allele-order']
        # TODO: Move the command to Support.run_command
        # Support.run_command(command_str=" ".join(command), shell=True)
        # Support.run_command(command_list=command)
        subprocess.check_output(" ".join(command), shell=True)

        if not os.path.isfile(chromosome_file_prefix + ".tped"):
            print(f"Failed to create {chromosome_file_prefix}.tped")
            sys.exit(1)

        # Read the thing and remove the temporary files
        chromosome_tped_file = chromosome_file_prefix + '.tped'
        chromosome_tped = pd.read_table(chromosome_tped_file, index_col=None, header=None, sep=' ').iloc[:, 4::2]

        # TODO: Move this to Python's rm command?
        command = ['rm'] + [chromosome_file_prefix + '.' + i for i in ['tped', 'tfam', 'nosex', 'log']]

        # TODO: Move the command to Support.run_command
        # Support.run_command(command_list=command)
        subprocess.check_output(" ".join(command), shell=True)
        logging.info(f"Finished loading SNPs from chr{chromosome}")
        return chromosome_tped

    def load_reference_tpeds(self):
        # Read in all of the reference TPEDs
        # print(self.chromosomes)
        temp_name = "pastrami_" + ''.join(random.choice(string.ascii_letters) for _ in range(10))
        Support.safe_dir_create(temp_name)
        logging.info(f"Loading reference TPEDs with {self.threads} threads")
        # pool = mp.Pool(processes=self.threads)
        # self.pool.restart()
        results = self.pool.map(
            lambda x: self.pull_chromosome_tped(prefix=self.reference_prefix, chromosome=x, temp_dir=temp_name),
            self.chromosomes)
        # results.wait()
        # results = results.get()
        # self.pool.close()
        # self.pool.join()
        # self.pool.terminate()
        Support.safe_dir_rm(temp_name)
        for i in range(len(self.chromosomes)):
            self.reference_tpeds[str(self.chromosomes[i])] = results[i]
            self.reference_tpeds[str(self.chromosomes[i])].columns = self.reference_individuals
        logging.info(f"Loaded {len(self.reference_tpeds)} referenece TPEDs")

    def load_query_tpeds(self):
        # Read in all of the reference TPEDs
        logging.info("Chromosomes to process: ")
        logging.info(self.chromosomes)
        temp_name = "pastrami_" + ''.join(random.choice(string.ascii_letters) for _ in range(10))
        Support.safe_dir_create(temp_name)
        logging.info(f"Loading query TPEDs with {self.threads} threads")
        # pool = mp.Pool(processes=self.threads)
        # self.pool.restart()
        results = self.pool.map(
            lambda x: self.pull_chromosome_tped(prefix=self.query_prefix, chromosome=x, temp_dir=temp_name),
            self.chromosomes)
        # results.wait()
        # results = results.get()
        # self.pool.close()
        # self.pool.join()
        # self.pool.terminate()
        Support.safe_dir_rm(temp_name)
        for i in range(len(self.chromosomes)):
            self.query_tpeds[str(self.chromosomes[i])] = results[i]
            self.query_tpeds[str(self.chromosomes[i])].columns = self.query_individuals
        logging.info(f"Loaded {len(self.query_tpeds)} query TPEDs")

    def build_reference_set(self):
        # Load up the reference data
        self.load_haplotypes()
        self.load_reference_tfam()
        self.load_reference_tpeds()

        #        self.haplotypes = self.haplotypes.iloc[self.haplotypes.iloc[:,0].values in self.chromosomes, :]

        # Parallelize across equal-sized chunks of haplotypes for good speed
        chunk_size = min([30, math.ceil(self.haplotypes.shape[0] / self.threads)])
        logging.info(f"Splitting haplotypes into chunks of {chunk_size}")
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

        # noinspection PyTypeChecker
        chunk_data = [[chunk_indices[i]] + [chunk_haplotypes[i]] + [chunk_genotypes[i]] + [i] for i in
                      range(len(chunk_indices))]

        del self.reference_tpeds

        # Gives us back a list of lists
        # Dim 1 - Chunk
        # Dim 2 - Results of a chunk
        # self.pool.restart()
        results = self.pool.map(self.chunk_build_reference_copying_fractions, chunk_data)
        # results.wait()
        # results = results.get()
        # self.pool.close()
        # self.pool.join()
        # self.pool.terminate()

        # Find the population counts for each reference haplotype
        self.reference_haplotype_counts = results[0][0]
        for i in range(1, len(results)):
            self.reference_haplotype_counts = pd.concat([self.reference_haplotype_counts, results[i][0]], axis=0)
        logging.info(self.reference_haplotype_counts.shape)

        # Find the population fractions for each reference haplotype
        self.reference_haplotype_fractions = results[0][1]
        for i in range(1, len(results)):
            self.reference_haplotype_fractions = pd.concat([self.reference_haplotype_fractions, results[i][1]], axis=0)

        # Find the reference indvidiual-reference population copying rates
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

        logging.info('Starting chunk ' + str(chunk_number))
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
        logging.info(f"Splitting haplotypes into chunks of {chunk_size}")
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

        # noinspection PyTypeChecker
        chunk_data = [
            [chunk_indices[i]] + [chunk_haplotypes[i]] + [chunk_genotypes[i]] + [chunk_reference_fractions[i]] + [i] for
            i in range(len(chunk_indices))]

        # Fuck the GIL
        del self.reference_haplotype_fractions
        del self.query_tpeds

        logging.info("Calculating reference copying fractions...")
        # pool = mp.Pool(processes=self.threads)
        # self.pool.restart()
        results = self.pool.map(self.chunk_query_reference_copying_fractions, chunk_data)
        # results.wait()
        # results = results.get()
        # self.pool.close()
        # self.pool.join()
        # self.pool.terminate()
        logging.info("Reference copying fractions calculated")
        logging.info("Finding the query indvidiual-reference population copying rates...")
        # Find the query indvidiual-reference population copying rates
        self.query_copying_fractions = results[0]
        for i in range(1, len(results)):
            self.query_copying_fractions += results[i]
        self.query_copying_fractions = self.query_copying_fractions.T.apply(lambda x: x / x.sum()).T

        # Remove the background seen in the reference set, making sure they don't go below 0
        self.query_copying_fractions -= self.reference_background
        self.query_copying_fractions[self.query_copying_fractions < 0] = 0
        self.query_copying_fractions = self.query_copying_fractions.T.apply(lambda x: x / x.sum()).T
        logging.info(f"Writing to files {self.query_output_file} and {self.combined_output_file}")
        self.query_copying_fractions.to_csv(self.query_output_file, sep='\t')

        self.combined_copying_fractions = pd.concat([self.reference_copying_fractions, self.query_copying_fractions],
                                                    axis=0)
        self.combined_copying_fractions.to_csv(self.combined_output_file, sep='\t')

    def chunk_query_reference_copying_fractions(self, chunk_data):

        chunk_indices, chunk_haplotypes, chunk_genotypes, chunk_reference_fractions, chunk_number = chunk_data

        chunk_query_individual_fractions = pd.DataFrame(0,
                                                        index=self.query_individuals,
                                                        columns=self.reference_populations)

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
                    chunk_query_individual_fractions.loc[individual, :] += \
                        haplotype_fractions.loc[individual_haplotype, :]

        return chunk_query_individual_fractions

    def build_coanc(self):

        # Load up the reference data
        self.load_haplotypes()
        self.load_reference_tfam()
        self.load_query_tfam()

        # pool = mp.Pool(processes=self.threads)
        # self.pool.restart()
        # For memory GIL BS, run each chromosome on its own
        copying_counts = pd.DataFrame(0, index=self.reference_individuals, columns=self.reference_individuals)
        haplotype_counts = pd.Series([0] * len(self.reference_individuals), index=self.reference_individuals)

        # TODO: This can potentially be parallelized
        for chromosome in self.chromosomes:

            # Just the haps for this chromosome
            chromosome_haplotypes = self.haplotypes.loc[self.haplotypes.iloc[:, 0] == chromosome, :]
            temp_name = "pastrami_" + ''.join(random.choice(string.ascii_letters) for _ in range(10))
            Support.safe_dir_create(temp_name)
            logging.info(f"Extracting reference TPEDs with {self.threads} threads for coanc")
            # Just the TPEDs for this chromosome
            reference_tped = self.pull_chromosome_tped(prefix=self.reference_prefix, chromosome=chromosome,
                                                       temp_dir=temp_name)
            reference_tped.columns = self.reference_individuals
            logging.info(f"Extracting query TPEDs with {self.threads} threads for coanc")
            query_tped = self.pull_chromosome_tped(prefix=self.query_prefix, chromosome=chromosome,
                                                   temp_dir=temp_name)
            Support.safe_dir_rm(temp_name)
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

            # noinspection PyTypeChecker
            chunk_data = [[chunk_indices[i]] +
                          [chunk_haplotypes[i]] +
                          [chunk_reference_genotypes[i]] +
                          [chunk_query_genotypes[i]] + [i] for i in range(len(chunk_indices))]

            del reference_tped
            del query_tped

            results = self.pool.map(self.chunk_build_coanc, chunk_data)
            # results.wait()
            # results = results.get()
            for i in range(len(results)):
                copying_counts += results[i][0]
                haplotype_counts += results[i][1]

        # self.pool.close()
        # self.pool.join()
        # self.pool.terminate()
        copying_counts = copying_counts.T.apply(lambda x: x / haplotype_counts).T
        # TODO: Figure out what this aVA.tsv file is
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

            # chromosome, left_snp, right_snp = chunk_haplotypes.iloc[haplotype, :]
            # haplotype_index = str(chromosome) + ':' + str(left_snp) + "-" + str(right_snp)
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

    """ 
        [Class section] Aggregate functions that are run post main Pastrami code
    """

    def post_pastrami(self):
        # Read files
        self.set_pop_group_map()
        self.process_fam_file()
        # Calculate and print ancestry fractions
        self.process_ancestry_fractions()
        self.print_ancestry_fractions()
        # Calculate and print painting vectors
        self.painting_vector_median()
        self.print_ancestry_paintings()
        # Calculate and print fine_grain and population_level estimates
        self.get_fine_grain_estimates()
        self.print_fine_grain_estimates()
        self.print_population_level_estimates()

    def set_pop_group_map(self):
        """Reads in the pop x group file of the format:
        #Population	Group
        GWD	Western
        MSL	Western
        Yacouba	Western
        Ahizi	Western
        Fon	Nigerian

        the two columns are tab separated. The file is read and stored in self.pop_group_map dictionary.

        Returns
        -------
        None
        Sets self.pop_group_map values and moves on
        """
        logging.info("Reading in population grouping file")
        with open(self.pop_group_file, "r") as f:
            for line in f:
                # ignore comment line
                if line.startswith("#"):
                    continue
                else:
                    # TODO: implement check for at least two columns in the input file
                    pop, group = line.rstrip().split("\t")
                    self.ref_pop_group_map[pop] = group
        logging.info("Population grouping file read")

    def process_fam_file(self):
        """Read in and save the FAM file in the self.ind_pop_map dictionary

        Returns
        -------
        None
        Sets self.ind_pop_map values and moves on
        """
        logging.info("Processing FAM file")
        # variables used for logging purposes only
        ref_ind_count = 0
        line_number = 0
        self.reference_individual_dict = {}
        with open(self.fam_infile, "r") as f:
            for line in f:
                # ignore comment lines
                if line.startswith("#"):
                    continue
                else:
                    # TODO: ensure that the input file has at least 3 columns
                    col_split = line.rstrip().split()
                    pop, ind = col_split[0], col_split[1]
                    ind = re.sub(r"\.[12]$", "", ind)
                    self.ind_pop_map[ind] = pop
                    if pop in self.ref_pop_group_map:
                        if ind not in self.reference_individual_dict:
                            self.reference_individual_dict[ind] = True
                            ref_ind_count += 1
                        else:
                            continue
                    line_number += 1
        logging.info(f"FAM file read, read {line_number} lines and found {ref_ind_count} reference individuals")

    def process_ancestry_fractions(self) -> None:
        """Does a bunch of things:
        1. Reads in the ancestry file and stores it in a dictionary of lists (key: individual -> list:ancestry frac.)
        2. Filters out populations (columns) not present in pop_group_map dict
        3. Filters out individuals (rows) not present in reference_individuals dict
        4. Average out fractions for each chromosome
        5. Scales the fractions by dividing each value by column minimum
        6. Second scaling by row - each value divided by sum(row)
        7. Populates the reverse map of population to individual

        Returns
        -------
        None
        Sets object values and moves on
        """
        logging.info("Reading in the Pastrami ancestry fractions")
        self.ancestry_fractions = {}
        with open(self.ancestry_infile, "r") as f:
            # Process the header containing the population names
            full_header = f.readline().strip().split("\t")
            indices = []  # will hold the columns to keep
            self.af_header = []  # will hold the ancestry_fraction header

            # Column filter to retain only reference populations
            for pop_index in range(len(full_header)):
                if full_header[pop_index] in self.ref_pop_group_map:
                    self.af_header.append(full_header[pop_index])
                    indices.append(pop_index + 1)  # +1 for individual id
            column_mins = [1] * len(self.af_header)
            logging.info(f"Found {len(self.af_header)} reference populations in the ancestry fraction file")

            # Iterate through the rest of the lines
            for line in f:
                columns = line.rstrip().split("\t")
                ind_id = re.sub(r"\.[12]", "", columns[0])

                # Column filter to retain only reference populations
                fractions = [float(columns[x]) for x in indices]

                # Populate the pop_ind_map dictionary of lists
                if ind_id not in self.ind_pop_map:
                    logging.error(f"Error: Individual id {ind_id} was not found in input file {self.ancestry_infile}")
                    sys.exit(1)

                this_pop = self.ind_pop_map[ind_id]
                if this_pop not in self.pop_ind_map:
                    self.pop_ind_map[this_pop] = []
                    self.pop_ind_map[this_pop].append(ind_id)
                    self.ancestry_fractions[ind_id] = fractions
                elif ind_id not in self.ancestry_fractions:
                    self.pop_ind_map[this_pop].append(ind_id)
                    self.ancestry_fractions[ind_id] = fractions
                else:
                    # Average the ancestry fractions and identify column minimums
                    # LR: I am assuming that we will only see an individual
                    # twice -> dictated by the regex and if condition above
                    for i in range(len(fractions)):
                        self.ancestry_fractions[ind_id][i] = (self.ancestry_fractions[ind_id][i] +
                                                              fractions[i]) / 2
                        if ind_id in self.reference_individual_dict and \
                                column_mins[i] > self.ancestry_fractions[ind_id][i]:
                            # noinspection PyTypeChecker
                            column_mins[i] = self.ancestry_fractions[ind_id][i]

        logging.info(f"Found {len(self.ancestry_fractions.keys())} reference individuals in the ancestry fraction file")

        # Scale factors by column minimum, followed by row sum
        logging.info("File read, scaling the fractions")
        for ind_id in self.ancestry_fractions:
            row_sum = 0
            # Column minimum scaling
            for i in range(len(self.af_header)):
                self.ancestry_fractions[ind_id][i] = self.ancestry_fractions[ind_id][i] - column_mins[i]
                if self.ancestry_fractions[ind_id][i] < 0:
                    self.ancestry_fractions[ind_id][i] = 0
                row_sum += self.ancestry_fractions[ind_id][i]
            # Row sum scaling
            # TODO: Check why row_sum is 0 in some cases
            if row_sum == 0:
                row_sum = 1
            for i in range(len(self.af_header)):
                self.ancestry_fractions[ind_id][i] = self.ancestry_fractions[ind_id][i] / row_sum
        logging.info("Ancestry fractions processed")

    def shrinkage(self, i_was_in_the_pool=0.1) -> None:
        for ind_id in self.ancestry_fractions:
            row_len = len(self.af_header)
            for i in range(row_len):
                self.ancestry_fractions[ind_id][i] += ((1 / row_len) - self.ancestry_fractions[ind_id][
                    i]) * i_was_in_the_pool

    def rescale(self) -> None:
        """Rescales ancestry_fractions so that summation of all values equals 1

        Returns
        -------
        None
        """
        for ind_id in self.ancestry_fractions:
            row_sum = 0
            for i in range(len(self.af_header)):
                if self.ancestry_fractions[ind_id][i] < 0:
                    self.ancestry_fractions[ind_id][i] = 0
                else:
                    row_sum += self.ancestry_fractions[ind_id][i]
            for i in range(len(self.af_header)):
                self.ancestry_fractions[ind_id][i] /= row_sum

    def painting_vector_median(self) -> None:
        """Calculates paint vector median for each population. Following steps are performed:
        1. For each reference population, number of individuals are calculated
        2. Reference populations with number of individuals less than 4 are discarded
        3. For the remaining populations, for each column, population mean and sd are calculated
        4. Population mean and sd are used to calculate population Z-score for each column
        5. Individuals with anamolous fractions, defined as fraction outside -5 <= Z <= 5, are discarded
        6. Population-scale per column medians are calculated and stored as painting vectors

        Returns
        -------
        None. Sets the value of painting vectors.
        """
        logging.info(f"Beginning calculation of median painting vectors")
        populations = sorted(self.ref_pop_group_map.keys())
        logging.info(f"{len(populations)} populations to process")
        self.painting_vectors = {}

        for this_pop in populations:
            # Throw away populations with 3 or less individuals
            if len(self.pop_ind_map[this_pop]) <= 3:
                logging.info(f"Population {this_pop} contains less than 4 individuals, discarding.")
                continue

            # Extract the individuals of interest for easy access
            pop_inds = self.pop_ind_map[this_pop]

            # Create 0 list for mean and sd calculation; mean and sd will be used for Z-score calculation
            this_mean = [0] * len(self.af_header)
            this_sd = [0] * len(this_mean)

            # Calculate population mean
            for ind_id in pop_inds:
                for col_id in range(len(self.af_header)):
                    this_mean[col_id] += self.ancestry_fractions[ind_id][col_id]
            this_mean = [i / len(pop_inds) for i in this_mean]

            # Calculate population SD
            for ind_id in pop_inds:
                for col_id in range(len(self.af_header)):
                    this_sd[col_id] += (self.ancestry_fractions[ind_id][col_id] - this_mean[col_id]) ** 2
            this_sd = [(i / (len(pop_inds) - 1)) ** 0.5 for i in this_sd]

            # Identify outliers; defined as individuals with ancestry fraction with Z-score > +/- 5
            keep_set = []
            for ind_id in pop_inds:
                drop_this_ind = False
                for col_id in range(len(self.af_header)):
                    this_z = 1.0
                    if this_sd[col_id] == 0:
                        #Standard deviation
                        logging.error(f"Error: Individual id {ind_id} has a standard deviation of {this_sd[col_id]}. This might be the case if your TFAM file did not have an entry for each haplotype of an indvidual. Each individual should have 2 enteries in the TFAM file.")
                        sys.exit(1)
                    else:
                        this_z = (self.ancestry_fractions[ind_id][col_id] - this_mean[col_id]) / this_sd[col_id]
                    # print(f"{ind_id} and {col_id}: {this_z}")
                    if abs(this_z) > 5:
                        drop_this_ind = True
                        logging.debug(f"Individual id {ind_id} contains outlier value " +
                                      f"(value for {self.af_header[col_id]}=" +
                                      f"{round(self.ancestry_fractions[ind_id][col_id], 3)}). Population mean = " +
                                      f"{round(this_mean[col_id], 3)} and sd = {round(this_sd[col_id], 3)}")
                        break
                if not drop_this_ind:
                    keep_set.append(ind_id)

            dropped = len(pop_inds) - len(keep_set)
            if dropped > 0:
                logging.info(f"For reference population {this_pop}, dropped {dropped} outlier individuals")

            # Calculate median ancestry fraction vectors from non-outlier individuals
            self.painting_vectors[this_pop] = []
            for col_id in range(len(self.af_header)):
                this_fraction = []
                for ind_id in keep_set:
                    this_fraction.append(self.ancestry_fractions[ind_id][col_id])
                self.painting_vectors[this_pop].append(statistics.median(this_fraction))

        self.painting_vectors_keys = list(self.painting_vectors.keys())
        logging.info(
            f"Finished calculating painting vectors for {len(self.painting_vectors)} reference populations")

    def regularized_rss(self, par: np.array, data: np.array) -> float:
        """Calculates regularized residual sum of square error and adds a non-negative penalty to it

        Parameters
        ----------
        par : np.array
            1-D numpy array containing parameters from which error needs to be computed

        data: np.array
            1-D numpy array containing the ancestry fractions of an individual

        Returns
        -------
        RSS error + penalty
        """
        error = 0
        non_negative_count = 0

        for value in par:
            if value > 0:
                non_negative_count += 1
        # non-negative count penalty
        penalty_multiplier = (non_negative_count - 1) * 0.25

        for col_id in range(len(self.af_header)):
            predicted = 0
            for row_id in range(len(self.painting_vectors_keys)):
                predicted += self.painting_vectors[self.painting_vectors_keys[row_id]][col_id] * par[row_id]

            error += (predicted - data[col_id]) ** 2
        # error calculated as residual sum of square (RSS) of predicted and actual
        penalty = error * penalty_multiplier
        # final unit that is to be minimized
        error_penalty = error + penalty
        # logging.debug(f"Error: {error}; Penalty: {penalty}; sum: {error_penalty}")
        return error_penalty

    def regularize_this_row(self, row: str) -> list:
        """Runs the minimize function on a single row from the ancestry fractions

        Parameters
        ----------
        row : str
            Ancestry fraction row name to process.

        Returns
        -------
        List
            A list of optimized estimates. The optimization can get stuck at local minima and
            may produce sub-optimal results.
        """
        logging.info(f"Optimizing {row}")
        data = np.array(self.ancestry_fractions[row])

        # default set to be same as data
        par = data
        this_estimate = minimize(fun=self.regularized_rss, args=data, x0=par, method='L-BFGS-B',
                                 bounds=[(0, 1) for _ in range(len(par))],
                                 options={'eps': Analysis.optim_step_size, 'maxiter': Analysis.optim_iterations},
                                 tol=Analysis.tolerance)
        logging.info(f"Optimized penalty for {row} (def case): " + str(this_estimate["fun"]))

        if this_estimate["fun"] <= Analysis.error_threshold:
            return this_estimate["x"]

        # case 2: equally weighted
        par = np.array([1 / len(self.painting_vectors)] * len(self.painting_vectors))
        new_estimate = minimize(fun=self.regularized_rss, args=data, x0=par, method='L-BFGS-B',
                                bounds=[(0, 1) for _ in range(len(par))],
                                options={'eps': Analysis.optim_step_size, 'maxiter': Analysis.optim_iterations},
                                tol=Analysis.tolerance)
        if new_estimate["fun"] <= Analysis.error_threshold:
            return new_estimate["x"]
        
        if new_estimate["fun"] < this_estimate["fun"]:
            this_estimate = new_estimate
        logging.info(f"Optimized penalty for {row} (case 2): " + str(this_estimate["fun"]))

        # case 3: try some additional random runs
        for iteration in range(len(self.painting_vectors)):
            par = np.random.dirichlet(np.ones(len(self.painting_vectors)),size=1)[0]
            new_estimate = minimize(fun=self.regularized_rss, args=data, x0=par, method='L-BFGS-B',
                                    bounds=[(0, 1) for _ in range(len(par))],
                                    options={'eps': Analysis.optim_step_size, 'maxiter': Analysis.optim_iterations},
                                    tol=Analysis.tolerance)
            
            if new_estimate["fun"] <= Analysis.error_threshold:
                return new_estimate["x"]
            
            if new_estimate["fun"] < this_estimate["fun"]:
                this_estimate = new_estimate
            logging.info(f"Optimized penalty for {row} (case 3, iteration {iteration}): " + str(new_estimate["fun"]))

        return this_estimate["x"]

    def get_fine_grain_estimates(self):
        """Calculates fine grain estimates of individuals listed within ancestry_fractions dictionary

        Returns
        -------
        None.
            Sets self.fine_grain_estimates to the fine grain estimates calculated
        """
        logging.info(f"Beginning calculation of NNLS estimates")
        # pool = mp.Pool(processes=self.threads)
        # self.pool.restart()
        
        ### The following line initiates the original Pastrami optimizer function. 
        #results = self.pool.map(self.regularize_this_row, self.ancestry_fractions.keys())
        
        # results.wait()
        # results = results.get()
        # self.pool.close()
        # self.pool.join()
        # self.pool.terminate()
        logging.info(f"Finished all NNLS estimate calculations")
        
        #Outosurcing the same information to an R optimizer. 
        #R optimizer gives us more accurate optimization results.
        logging.info(f"Optimizing the same estimates using an outsourced R optimizer.")
        outsourcing_command_list = ["Rscript", "outsourcedOptimizer.R", "-p", self.out_prefix + Analysis.ancestry_painting_postfix,
                               "-f",  self.out_prefix + Analysis.ancestry_fraction_postfix,
                               "-o", self.out_prefix + Analysis.outsourced_optimizer_pop_estimates_postfix]
        
        try:
            logging.info(f"Attempting to run: {' '.join(outsourcing_command_list)}")
            output = subprocess.check_output(outsourcing_command_list)

        except subprocess.CalledProcessError as e:
            #Outsourcing failed.
            logging.error(f"Encountered an error executing the command: {outsourcing_command_list}")
            logging.error(f"Exit code={e.returncode}")
            logging.error(f"Error message={e.output}")
        
        logging.info(f"Finished all outsourced optimizations.")
        #R outsourcing optimizer worked!        
        
        #Now let's incorporate the outsourced optimzer numbers back to the results variable.
        with open(self.out_prefix + Analysis.outsourced_optimizer_pop_estimates_postfix) as f:
            outsourced_optimizer_file_content = f.read()
        
        '''
        Convert the outsourced optimizer matrix into a list of np array 
        This should be similar to the results list varibale declared above in this function.   
        The TSV file created by outsourced optimizer is transformed to the following variable ([sample1, sample2, ...] for 4 populations):
            [array([0.      , 0.432669, 0.      , 0.567331]), array([0.      , 0.335317, 0.      , 0.664683]), ... ]
        '''
        outsourced_optimizer_results = [np.array([0 if (j == "NA" or not j.replace(".","").isdigit()) else float(j) for j in i.split("\t")[1:]], dtype='float') for i in outsourced_optimizer_file_content.split("\n")[1:] if i != ""]
        
        #Put the results back to the original results variable of this function.
        results = outsourced_optimizer_results
        
        #Delete the outsourced optimizer variables.
        del outsourced_optimizer_file_content, outsourced_optimizer_results
        
        #Carry on with the standard Pastrami code.
        ind_ids = list(self.ancestry_fractions.keys())
        for i in range(len(ind_ids)):
            this_sum = sum(results[i])
            # TODO: Investigate why sum is sometimes 0
            if this_sum == 0:
                this_sum = 1
            for j in range(len(results[i])):
                results[i][j] = results[i][j] / this_sum
            self.fine_grain_estimates[ind_ids[i]] = results[i]

    def print_ancestry_fractions(self):
        """Prints ancestry fractions to a file

        Returns
        -------
        None.
            Prints the ancestry fractions to {out_prefix}_fractions.Q file
        """
        outfile = self.out_prefix + Analysis.ancestry_fraction_postfix
        logging.info(f"Writing ancestry fractions to {outfile}")
        with open(outfile, "w") as f:
            f.write("Ind\t" + "\t".join(self.af_header) + "\n")
            for ind_id in self.ancestry_fractions:
                f.write(f"{ind_id}\t" + "\t".join([str(x) for x in self.ancestry_fractions[ind_id]]) + "\n")
        logging.info(f"{outfile} successfully created")

    def print_ancestry_paintings(self):
        """Prints ancestry paintings to a file

        Returns
        -------
        None.
            Prints the ancestry fractions to {out_prefix}_paintings.Q file
        """
        outfile = self.out_prefix + Analysis.ancestry_painting_postfix
        logging.info(f"Writing ancestry paintings to {outfile}")
        with open(outfile, "w") as f:
            f.write("Pop\t" + "\t".join(self.af_header) + "\n")
            for pop in self.painting_vectors:
                f.write(pop + "\t" + "\t".join([str(i) for i in self.painting_vectors[pop]]) + "\n")
        logging.info(f"{outfile} successfully created")

    def print_fine_grain_estimates(self):
        """Prints fine grain estimates to a file

        Returns
        -------
        None.
            Prints the fine grain estimates into {out_prefix}_fine_grain_estimates.Q file
        """
        outfile = self.out_prefix + Analysis.finegrain_estimates_postfix
        logging.info(f"Writing fine grain estimates to {outfile}")
        with open(outfile, "w") as f:
            f.write("Id\t" + "\t".join(self.painting_vectors.keys()) + "\n")
            for ind_id in self.fine_grain_estimates:
                f.write(f"{ind_id}\t" + "\t".join([str(x) for x in self.fine_grain_estimates[ind_id]]) + "\n")
        logging.info(f"{outfile} successfully created")

    def print_population_level_estimates(self):
        """Calculates population-level estimates and prints them to a file

        Returns
        -------
        None.
            Prints the population-level estimates to {out_prefix}_fractions.Q file
        """
        pop_level_estimates = {}

        column_name = list(self.painting_vectors.keys())
        for pop_index in range(len(column_name)):
            pop = column_name[pop_index]
            if pop not in self.ref_pop_group_map:
                raise ValueError(f"Can't seem to find {pop} in the {self.pop_group_file}.")
            pop_groups = self.ref_pop_group_map[pop]
            if pop_groups not in pop_level_estimates:
                pop_level_estimates[pop_groups] = []
            pop_level_estimates[pop_groups].append(pop_index)

        pop_groups = sorted(pop_level_estimates.keys())

        outfile = self.out_prefix + Analysis.pop_estimates_postfix
        logging.info(f"Writing population-level estimates to {outfile}")
        with open(outfile, "w") as f:
            f.write("Ind\t" + "\t".join(pop_groups) + "\n")
            for ind_id in self.fine_grain_estimates:
                f.write(ind_id)
                for this_pop_group in pop_groups:
                    this_group_sum = 0
                    for pop_index in pop_level_estimates[this_pop_group]:
                        this_group_sum += self.fine_grain_estimates[ind_id][pop_index]
                    f.write(f"\t{this_group_sum}")
                f.write("\n")
        logging.info(f"{outfile} created successfully")

    """ 
        End of class
    """


if __name__ == '__main__':
    parser = ArgumentParser(prog=PROGRAM_NAME,
                            add_help=False,
                            description=f'''
                            {PROGRAM_NAME} - population scale haplotype copying
                            ''',
                            formatter_class=lambda prog: HelpFormatter(prog, width=120, max_help_position=120))

    # TODO: Make these options work in way that they can be used after the subcommand
    parser.add_argument('--help', '-h', '--h', action='store_true', default=False)

    subparsers = parser.add_subparsers(title=f"{PROGRAM_NAME} commands")

    # Make the "all" sub-command parser
    all_parser = subparsers.add_parser('all', help='Perform full analysis',
                                       formatter_class=lambda prog: HelpFormatter(prog, width=120,
                                                                                  max_help_position=120))
    all_parser.set_defaults(sub_command='all')

    # Add build input arguments
    all_input_group = all_parser.add_argument_group('Required Input options')
    all_input_group.add_argument('--reference-prefix', required=False, default=None, metavar='<PREFIX>',
                                 help='Prefix for the reference TPED/TFAM input files')
    all_input_group.add_argument('--query-prefix', required=False, default=None, metavar='<TSV>',
                                 help='Prefix for the query TPED/TFAM input files')
    all_input_group.add_argument('--haplotypes', required=False, default=None, metavar='<TSV>',
                                 help='File of haplotype positions')
    all_input_group.add_argument('--pop-group', required=False, default=None, metavar='pop2group.txt', type=str,
                                 help='File containing population to group (e.g., tribes to region) mapping',
                                 dest="pop_group_file")

    # Add build output arguments
    all_output_group = all_parser.add_argument_group('Output options')
    all_output_group.add_argument('--out-prefix', required=False, default="pastrami", metavar='out-prefix',
                                  type=str, dest="out_prefix",
                                  help="""Output prefix (default: %(default)s) for creating following sets of files.\n 
                                         <prefix>.pickle, \n
                                         <prefix>_query.tsv, \n
                                         <prefix>.tsv, \n
                                         <prefix>.fam, \n
                                         <prefix>_fractions.Q, \n
                                         <prefix>_paintings.Q, \n
                                         <prefix>_estimates.Q, \n
                                         <prefix>_fine_grain_estimates.Q\n""")

    all_hapmake_group = all_parser.add_argument_group('Optional arguments for hapmake stage')
    all_hapmake_group.add_argument('--map-dir', required=False, default=None, metavar='maps/', type=str,
                                   help='Directory containing genetic maps: chr1.map, chr2.map, etc',
                                   dest="map_dir")
    all_hapmake_group.add_argument('--min-snps', required=False, default=7, metavar='MinSNPs', type=int,
                                   help='Minimum number of SNPs in a haplotype block (default: %(default)s)',
                                   dest="min_snps")
    all_hapmake_group.add_argument('--max-snps', required=False, default=20, metavar='MaxSNPs', type=int,
                                   help='Maximum number of SNPs in a haplotype block (default: %(default)s)',
                                   dest="max_snps")
    all_hapmake_group.add_argument('--max-rate', required=False, default=0.3, metavar='MaxRate', type=float,
                                   help='Maximum recombination rate (default: %(default)s)',
                                   dest="max_rate")

    all_runtime_group = all_parser.add_argument_group('Runtime options')
    all_runtime_group.add_argument('--per-individual', required=False, default=False, action='store_true',
                                   help='Generate per-individual copying rather than per-population copying')
    all_runtime_group.add_argument('--log-file', "-l", "--l", required=False, default="run.log", metavar='run.log',
                                   type=str,
                                   help='File containing log information (default: %(default)s)', dest="log_file")
    all_runtime_group.add_argument('--threads', "-t", required=False, default=4, metavar='N', type=int,
                                   help='Number of concurrent threads (default: %(default)s)', dest="threads")
    all_runtime_group.add_argument('--verbose', "-v", required=False, action="store_true",
                                   help='Print program progress information on screen', dest="verbosity")

    # Make the "haplotype" maker command
    hapmake_parser = subparsers.add_parser('hapmake', help='Create haplotypes',
                                           formatter_class=lambda prog: HelpFormatter(prog, width=120,
                                                                                      max_help_position=120))
    hapmake_parser.set_defaults(sub_command='hapmake')

    hapmake_input_group = hapmake_parser.add_argument_group('Input options')
    hapmake_input_group.add_argument('--map-dir', required=False, default=None, metavar='maps/', type=str,
                                     help='Directory containing genetic maps: chr1.map, chr2.map, etc',
                                     dest="map_dir")
    hapmake_input_group.add_argument('--min-snps', required=False, default=7, metavar='MinSNPs', type=int,
                                     help='Minimum number of SNPs in a haplotype block (default: %(default)s)',
                                     dest="min_snps")
    hapmake_input_group.add_argument('--max-snps', required=False, default=20, metavar='MaxSNPs', type=int,
                                     help='Maximum number of SNPs in a haplotype block (default: %(default)s)',
                                     dest="max_snps")
    hapmake_input_group.add_argument('--max-rate', required=False, default=0.3, metavar='MaxRate', type=float,
                                     help='Maximum recombination rate (default: %(default)s)',
                                     dest="max_rate")

    hapmake_output_group = hapmake_parser.add_argument_group('Output options')
    hapmake_output_group.add_argument('--haplotypes', required=False, default="out.haplotypes",
                                      metavar='out.haplotypes',
                                      type=str,
                                      help='Output file containing haplotypes (default: %(default)s)')

    hapmake_runtime_group = hapmake_parser.add_argument_group('Runtime options')
    hapmake_runtime_group.add_argument('--log-file', "-l", "--l", required=False, default="run.log", metavar='run.log',
                                       type=str,
                                       help='File containing log information (default: %(default)s)', dest="log_file")
    hapmake_runtime_group.add_argument('--threads', "-t", required=False, default=4, metavar='N', type=int,
                                       help='Number of concurrent threads (default: %(default)s)', dest="threads")
    hapmake_runtime_group.add_argument('--verbose', "-v", required=False, action="store_true",
                                       help='Print program progress information on screen', dest="verbosity")

    # Make the "build" sub-command parser
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

    build_runtime_group = build_parser.add_argument_group('Runtime options')
    build_runtime_group.add_argument('--per-individual', required=False, default=False, action='store_true',
                                     help='Generate per-individual copying rather than per-population copying')
    build_runtime_group.add_argument('--log-file', "-l", "--l", required=False, default="run.log", metavar='run.log',
                                     type=str,
                                     help='File containing log information (default: %(default)s)', dest="log_file")
    build_runtime_group.add_argument('--threads', "-t", required=False, default=4, metavar='N', type=int,
                                     help='Number of concurrent threads (default: %(default)s)', dest="threads")
    build_runtime_group.add_argument('--verbose', "-v", required=False, action="store_true",
                                     help='Print program progress information on screen', dest="verbosity")

    # Make the "query" sub-command
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

    query_runtime_group = query_parser.add_argument_group('Runtime options')
    query_runtime_group.add_argument('--log-file', "-l", "--l", required=False, default="run.log", metavar='run.log',
                                     type=str,
                                     help='File containing log information (default: %(default)s)', dest="log_file")
    query_runtime_group.add_argument('--threads', "-t", required=False, default=4, metavar='N', type=int,
                                     help='Number of concurrent threads (default: %(default)s)', dest="threads")
    query_runtime_group.add_argument('--verbose', "-v", required=False, action="store_true",
                                     help='Print program progress information on screen', dest="verbosity")

    # Make the "co-ancestry" sub-command
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

    coanc_runtime_group = coanc_parser.add_argument_group('Runtime options')
    coanc_runtime_group.add_argument('--log-file', "-l", "--l", required=False, default="run.log", metavar='run.log',
                                     type=str,
                                     help='File containing log information (default: %(default)s)', dest="log_file")
    coanc_runtime_group.add_argument('--threads', "-t", required=False, default=4, metavar='N', type=int,
                                     help='Number of concurrent threads (default: %(default)s)', dest="threads")
    coanc_runtime_group.add_argument('--verbose', "-v", required=False, action="store_true",
                                     help='Print program progress information on screen', dest="verbosity")

    # Make the "aggregate" maker command
    aggregate_parser = subparsers.add_parser('aggregate', help='Aggregate results and estimate ancestries',
                                             formatter_class=lambda prog: HelpFormatter(prog, width=120,
                                                                                        max_help_position=120))
    aggregate_parser.set_defaults(sub_command='aggregate')

    aggregate_input_group = aggregate_parser.add_argument_group('Input options')
    aggregate_input_group.add_argument('--pastrami-output', required=False, default=None, metavar='ancestry.tsv',
                                       type=str, dest="ancestry_infile",
                                       help="Output file generated from Pastrami's query subcommand")
    aggregate_input_group.add_argument('--pop-group', required=False, default=None, metavar='pop2group.txt', type=str,
                                       help='File containing population to group (e.g., tribes to region) mapping',
                                       dest="pop_group_file")
    aggregate_input_group.add_argument('--pastrami-fam', required=False, default=None, metavar='african.fam', type=str,
                                       help='File containing individual and population mapping in FAM format',
                                       dest="fam_infile")

    aggregate_output_group = aggregate_parser.add_argument_group('Output options')
    aggregate_output_group.add_argument('--out-prefix', required=False, default="pastrami", metavar='out-prefix',
                                        type=str, dest="out_prefix",
                                        help="""Output prefix for ancestry estimates files (default: %(default)s). 
                                         Four files are created: 
                                         <prefix>_fractions.Q, 
                                         <prefix>_paintings.Q, 
                                         <prefix>_estimates.Q, 
                                         <prefix>_fine_grain_estimates.Q\n""")

    aggregate_runtime_group = aggregate_parser.add_argument_group('Runtime options')
    aggregate_runtime_group.add_argument('--log-file', "-l", "--l", required=False, default="run.log",
                                         metavar='run.log', type=str, dest="log_file",
                                         help='File containing log information (default: %(default)s)')
    aggregate_runtime_group.add_argument('--threads', "-t", required=False, default=4, metavar='N', type=int,
                                         help='Number of concurrent threads (default: %(default)s)', dest="threads")
    aggregate_runtime_group.add_argument('--verbose', "-v", required=False, action="store_true",
                                         help='Print program progress information on screen', dest="verbosity")

    options, unknown_arguments = parser.parse_known_args()

    if len(unknown_arguments) > 0:
        print(Colors.HEADER)
        print("User specificed unknown arguments: " + " ".join([str(x) for x in unknown_arguments]))
        print("Please see the correct usage below:")
        parser.print_help()
        print(Colors.ENDC)
        sys.exit(0)

    if options.help or "sub_command" not in options:
        print(Colors.HEADER)
        parser.print_help()
        print(Colors.ENDC)
        sys.exit(0)

    if options.sub_command == 'all' and options.help:
        print(Colors.HEADER)
        all_parser.print_help()
        print(Colors.ENDC)
        sys.exit(0)

    if options.sub_command == 'hapmake' and options.help:
        print(Colors.HEADER)
        hapmake_parser.print_help()
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

    if options.sub_command == 'aggregate' and options.help:
        print(Colors.HEADER)
        aggregate_parser.print_help()
        print(Colors.ENDC)
        sys.exit(0)

    # Start the analysis
    analysis = Analysis(options)

    analysis.validate_options()

    if len(analysis.errors) > 0:
        print(Colors.HEADER)
        if options.sub_command == 'all':
            all_parser.print_help()
        elif options.sub_command == 'build':
            build_parser.print_help()
        elif options.sub_command == 'query':
            query_parser.print_help()
        elif options.sub_command == 'coanc':
            coanc_parser.print_help()
        elif options.sub_command == "hapmake":
            hapmake_parser.print_help()
        elif options.sub_command == "aggregate":
            aggregate_parser.print_help()
        print(Colors.ENDC)
        print(Colors.FAIL + '\n\nErrors:')
        print("\n".join(analysis.errors))
        print(Colors.ENDC)
        sys.exit()
    
    # If we're still good, start the actual analysis
    analysis.go()
