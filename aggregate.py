#!/usr/bin/env python3

"""Script for post-processing pastrami results"""

__author__ = "Andrew Conley, Lavanya Rishishwar"
__copyright__ = "Copyright 2020, Andrew Conley, Lavanya Rishishwar"
__credits__ = ["Andrew Conely", "Lavanya Rishishwar"]
__license__ = "GPL"
__version__ = "0.1"
__maintainer__ = "Andrew Conley, Lavanya Rishishwar"
__email__ = "aconley@ihrc.com; lrishishwar@ihrc.com"
__status__ = "Development"

from argparse import ArgumentParser, HelpFormatter
import sys
import statistics
import pathos.multiprocessing as mp
import re
import pandas as pd
import numpy as np
import logging
from scipy.optimize import minimize

PROGRAM_NAME = "aggregate.py"


class Colors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'


class Aggregate:
    debug = True
    optim_step_size = 0.0001

    def __init__(self, opts):
        self.ancestry_infile = opts.ancestry_infile
        self.pop_group_file = opts.pop_group_file
        self.out_prefix = opts.out_prefix
        self.threads = opts.threads
        self.fam_infile = opts.fam_infile
        self.verbosity = opts.verbosity
        self.log_file = opts.log_file

        self.pop_group_map = {}
        self.ind_pop_map = {}
        self.pop_ind_map = {}  # reverse map of ind_pop_map, sacrificing memory for speed later on
        self.reference_individuals = {}
        self.reference_pops = {}
        self.ancestry_fractions = {}
        self.af_header = []
        self.painting_vectors = {}
        self.painting_vectors_keys = []
        self.fine_grain_estimates = {}

        self.errors = []

        self.init_logger()

    def validate_options(self):
        pass

    def go(self):
        function = getattr(self, "post_pastrami")
        function()

    def post_pastrami(self):
        self.set_pop_group_map()
        self.process_fam_file()
        self.process_ancestry_fractions()
        self.painting_vector_median()
        self.get_fine_grain_estimates()
        self.print_ancestry_fractions()
        self.print_ancestry_paintings()
        self.print_fine_grain_estimates()
        self.print_population_level_estimates()

    def set_pop_group_map(self):
        logging.info("Reading in population grouping file")
        with open(self.pop_group_file, "r") as f:
            for line in f:
                if line.startswith("#"):
                    continue
                else:
                    pop, group = line.rstrip().split("\t")
                    self.pop_group_map[pop] = group
        logging.info("Population grouping file read")

    def process_fam_file(self):
        logging.info("Processing FAM file")
        ref_ind_count = 0
        line_number = 0
        with open(self.fam_infile, "r") as f:
            for line in f:
                if line.startswith("#"):
                    continue
                else:
                    pop, ind, extra = line.rstrip().split("\t", 2)
                    ind = re.sub(r"\.[12]", "", ind)
                    self.ind_pop_map[ind] = pop
                    if pop in self.pop_group_map:
                        self.reference_individuals[ind] = True
                        ref_ind_count += 1
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
                if full_header[pop_index] in self.pop_group_map:
                    self.af_header.append(full_header[pop_index])
                    indices.append(pop_index + 1)  # +1 for individual id
            column_mins = [1] * len(self.af_header)
            logging.info(f"Found {len(self.af_header)} reference populations in the ancestry fraction file")

            # Iterate through the rest of the lines
            ind_seen = {}
            for line in f:
                columns = line.rstrip().split("\t")
                ind_id = re.sub(r"\.[12]", "", columns[0])

                # Column filter to retain only reference populations
                fractions = [float(columns[x]) for x in indices]

                # Row filter to retain only reference individuals
                if ind_id not in self.reference_individuals:
                    continue

                # Populate the pop_ind_map dictionary of lists
                this_pop = self.ind_pop_map[ind_id]
                if this_pop not in self.pop_ind_map:
                    self.pop_ind_map[this_pop] = []
                if ind_id not in ind_seen:
                    self.pop_ind_map[this_pop].append(ind_id)
                    ind_seen[ind_id] = 1

                # Average the ancestry fractions and identify column minimums
                if ind_id in self.ancestry_fractions:
                    # @LR: I am assuming that we will only see an individual
                    # twice -> dictated by the regex and if condition above
                    for i in range(len(fractions)):
                        self.ancestry_fractions[ind_id][i] = (self.ancestry_fractions[ind_id][i] +
                                                              fractions[i]) / 2
                        if column_mins[i] > self.ancestry_fractions[ind_id][i]:
                            column_mins[i] = self.ancestry_fractions[ind_id][i]
                else:
                    self.ancestry_fractions[ind_id] = fractions

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
            for i in range(len(self.af_header)):
                self.ancestry_fractions[ind_id][i] = self.ancestry_fractions[ind_id][i] / row_sum
        logging.info("Ancestry fractions processed")

    def process_ancestry_fractions_pdway(self):
        logging.info("Reading in the Pastrami ancestry fractions")
        self.ancestry_fractions = pd.read_table(self.ancestry_infile, index_col=0, header=0, sep="\t")
        self.ancestry_fractions.aggregate()
        logging.info("Ancestry fraction read")

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
        populations = sorted(self.pop_group_map.keys())
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
                    this_z = (self.ancestry_fractions[ind_id][col_id] - this_mean[col_id]) / this_sd[col_id]
                    # print(f"{ind_id} and {col_id}: {this_z}")
                    if abs(this_z) > 5:
                        drop_this_ind = True
                        logging.debug(
                            f"Individual id {ind_id} contains outlier value (value for {self.af_header[col_id]}={round(self.ancestry_fractions[ind_id][col_id], 3)}). Population mean = {round(this_mean[col_id], 3)} and sd = {round(this_sd[col_id], 3)}")
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
        logging.info(f"Finished calculating painting vectors for {len(self.painting_vectors)} reference populations")

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
        logging.debug(f"Error: {error}; Penalty: {penalty}; sum: {error_penalty}")
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
        par = np.array([1 / len(self.painting_vectors)] * len(self.painting_vectors))
        data = np.array(self.ancestry_fractions[row])
        this_estimate = minimize(fun=self.regularized_rss, args=data, x0=par, method='L-BFGS-B',
                                 bounds=[(0, 1) for _ in range(len(par))], options={'eps': Aggregate.optim_step_size})
        logging.info(f"Optimized penalty for {row}: " + str(this_estimate["fun"]))
        return this_estimate["x"]

    def print_ancestry_fractions(self):
        """Prints ancestry fractions to a file

        Returns
        -------
        None.
            Prints the ancestry fractions to {out_prefix}_fractions.Q file
        """
        outfile = self.out_prefix + "_fractions.Q"
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
        outfile = self.out_prefix + "_paintings.Q"
        logging.info(f"Writing ancestry paintings to {outfile}")
        with open(outfile, "w") as f:
            f.write("Pop\t" + "\t".join(self.af_header) + "\n")
            for pop in self.painting_vectors:
                f.write(pop + "\t" + "\t".join([str(i) for i in self.painting_vectors[pop]]) + "\n")
        logging.info(f"{outfile} successfully created")

    def get_fine_grain_estimates(self):
        """Calculates fine grain estimates of individuals listed within ancestry_fractions dictionary

        Returns
        -------
        None.
            Sets self.fine_grain_estimates to the fine grain estimates calculated
        """
        logging.info(f"Beginning calculation of NNLS estimates")
        pool = mp.Pool(processes=self.threads)
        results = pool.map(self.regularize_this_row, self.ancestry_fractions.keys())
        logging.info(f"Finished all NNLS estimate calculations")

        ind_ids = list(self.ancestry_fractions.keys())
        for i in range(len(ind_ids)):
            self.fine_grain_estimates[ind_ids[i]] = results[i]

    def print_fine_grain_estimates(self):
        """Prints fine grain estimates to a file

        Returns
        -------
        None.
            Prints the fine grain estimates into {out_prefix}_fine_grain_estimates.Q file
        """
        outfile = self.out_prefix + "_fine_grain_estimates.Q"
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
        pops = list(self.painting_vectors.keys())
        for pop_index in range(len(pops)):
            pop = pops[pop_index]
            if pop not in self.pop_group_map:
                raise ValueError(f"Can't seem to find {pop} in the {self.pop_group_file}.")
            pop_groups = self.pop_group_map[pop]
            if pop_groups not in pop_level_estimates:
                pop_level_estimates[pop_groups] = []
            pop_level_estimates[pop_groups].append(pop_index)

        pop_groups = sorted(pop_level_estimates.keys())

        outfile = self.out_prefix + "_estimates.Q"
        logging.info(f"Writing population-level estimates to {outfile}")
        with open(outfile, "w") as f:
            f.write("Ind\t" + "\t".join(pop_groups) + "\n")
            for ind_id in self.ancestry_fractions:
                f.write(ind_id)
                for this_pop in pop_groups:
                    this_group_sum = 0
                    for pop_index in pop_level_estimates[this_pop]:
                        this_group_sum += self.ancestry_fractions[ind_id][pop_index]
                    f.write(f"\t{this_group_sum}")
                f.write("\n")
        logging.info(f"{outfile} created successfully")

    def init_logger(self):
        """Configures the logging for printing
        Returns
        -------
        None
            Logger behavior is set based on the Inputs variable
        """
        try:
            logging.basicConfig(filename=self.log_file, filemode="w", level=logging.DEBUG,
                                format=f"[%(asctime)s] %(message)s",
                                datefmt="%m-%d-%Y %I:%M:%S %p")
        except FileNotFoundError:
            print(
                f"The supplied location for the log file '{self.log_file}' doesn't exist. Please check if the location exists.")
            sys.exit(1)
        except IOError:
            print(
                f"I don't seem to have access to make the log file. Are the permissions correct or is there a directory with the same name?")
            sys.exit(1)

        if self.verbosity:
            console = logging.StreamHandler()
            console.setLevel(logging.INFO)
            formatter = logging.Formatter(fmt=f"[%(asctime)s] %(message)s", datefmt="%m-%d-%Y %I:%M:%S %p")
            console.setFormatter(formatter)
            logging.getLogger().addHandler(console)


if __name__ == '__main__':
    parser = ArgumentParser(prog=PROGRAM_NAME,
                            add_help=False,
                            description=f'''
                            {PROGRAM_NAME} - population scale haplotype copying
                            ''',
                            formatter_class=lambda prog: HelpFormatter(prog, width=120, max_help_position=120))

    parser.add_argument('--help', '-h', '--h', action='store_true', default=False)
    subparsers = parser.add_subparsers(title=f"{PROGRAM_NAME} commands")

    # Make the haplotype maker command
    aggregate_parser = subparsers.add_parser('aggregate', help='Create haplotypes',
                                             formatter_class=lambda prog: HelpFormatter(prog, width=120,
                                                                                        max_help_position=120))
    aggregate_parser.set_defaults(sub_command='aggregate')

    aggregate_input_group = aggregate_parser.add_argument_group('Input options')
    aggregate_input_group.add_argument('--pastrami-output', required=True, default=None, metavar='ancestry.tsv',
                                       type=str,
                                       help="Output file generated from Pastrami's query subcommand",
                                       dest="ancestry_infile")
    aggregate_input_group.add_argument('--pop-group', required=True, default=None, metavar='pop2group.txt', type=str,
                                       help='File containing population to group (e.g., tribes to region) mapping',
                                       dest="pop_group_file")
    aggregate_input_group.add_argument('--pastrami-fam', required=True, default=None, metavar='african.fam', type=str,
                                       help='File containing individual and population mapping in FAM format',
                                       dest="fam_infile")
    aggregate_input_group.add_argument('--threads', required=False, default=4, metavar='Threads', type=int,
                                       help='Number of concurrent threads (default: %(default)s)', dest="threads")
    aggregate_input_group.add_argument('--verbose', "-v", required=False, action="store_true",
                                       help='Print program progress information on screen', dest="verbosity")
    aggregate_output_group = aggregate_parser.add_argument_group('Output options')
    aggregate_output_group.add_argument('--out-prefix', required=False, default="pastrami", metavar='out-prefix',
                                        type=str, dest="out_prefix",
                                        help="""Output prefix for ancestry estimates files (default: %(default)s). 
                                     Four files are created:\n 
                                     * <prefix>_fractions.Q, \n
                                     * <prefix>_paintings.Q,\n
                                     * <prefix>_estimates.Q, \n
                                     * <prefix>_fine_grain_estimates.Q\n""")
    aggregate_output_group.add_argument('--log-file', required=False, default="run.log", metavar='run.log', type=str,
                                        help='File containing log information (default: %(default)s)', dest="log_file")

    options, unknown_arguments = parser.parse_known_args()
    # options, unknown_arguments = parser.parse_known_args(
    #     "aggregate --verbose --pastrami-output african.tsv --pop-group pop2group.txt --pastrami-fam maskedAfricanChrAllAboveMinimum.fam --threads 1".split(
    #         " "))

    if options.help or "sub_command" not in options:
        print(Colors.HEADER)
        parser.print_help()
        print(Colors.ENDC)
        sys.exit(0)

    if options.sub_command == 'aggregate' and options.help:
        print(Colors.HEADER)
        aggregate_parser.print_help()
        print(Colors.ENDC)
        sys.exit(0)

    # Start the analysis
    analysis = Aggregate(options)

    analysis.validate_options()

    if len(analysis.errors) > 0:
        print(Colors.HEADER)
        if options.sub_command == 'aggregate':
            aggregate_parser.print_help()
        print(Colors.ENDC)
        print(Colors.FAIL + '\n\nErrors:')
        [print(i) for i in analysis.errors]
        print(Colors.ENDC)
        sys.exit()

    # If we're still good, start the actual analysis
    analysis.go()
