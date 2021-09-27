#!/usr/bin/env python

"""Describe the script here"""

__author__ = "Lavanya Rishishwar"
__copyright__ = "Copyright 2020, Lavanya Rishishwar"
__credits__ = ["Lavanya Rishishwar"]
__license__ = "GPL"
__version__ = "0.1"
__maintainer__ = "Lavanya Rishishwar"
__email__ = "lavanyarishishwar@gmail.com"
__status__ = "Development"

import logging
import random
import re
import sys
import os
import pandas as pd
import subprocess
from argparse import ArgumentParser, HelpFormatter

VERSION = 0.1
PROGRAM_NAME = "simulate_individuals.py"


class Colors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'


class Simulation:
    def __init__(self, opts):
        # General attributes
        self.ref_pops = opts.ref_pops
        self.ref_percs = opts.ref_percs
        self.log_file = opts.log_file
        self.verbose = opts.verbose
        self.hap_file = opts.hap_file
        self.dir = opts.dir
        self.out_prefix = opts.out_prefix
        self.num_inds = opts.num
        self.errors = []
        self.haplotypes = []
        self.ref_tfam = {}
        self.ref_tped = {}

        self.simulated_tfam = []
        self.simulated_tped = pd.DataFrame()
        self.simulated_fractions = {}
        self.ref_for_sim_tfam = {}
        self.trace = []  # will store tracing of how we made individuals
        self.ref_sub_tfam = None
        self.ref_sub_tped = None
        os.environ['NUMEXPR_MAX_THREADS'] = str(opts.threads)

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
            print(f"The supplied location for the log file '{self.log_file}'" +
                  f"doesn't exist. Please check if the location exists.")
            sys.exit(1)
        except IOError:
            print(f"I don't seem to have access to make the log file." +
                  f"Are the permissions correct or is there a directory with the same name?")
            sys.exit(1)

        if self.verbose:
            console = logging.StreamHandler()
            console.setLevel(logging.INFO)
            formatter = logging.Formatter(fmt=f"[%(asctime)s] %(message)s", datefmt="%m-%d-%Y %I:%M:%S %p")
            console.setFormatter(formatter)
            logging.getLogger().addHandler(console)

    def validate_options(self):
        # Are reference populations supplied?
        if self.ref_pops is None:
            self.errors += ['Reference populations (-r or --reference-pops) is required']
        else:
            self.ref_pops = self.ref_pops.split(",")

        # Are reference proportions supplied?
        if self.ref_percs is None:
            self.errors += ['Reference percentages (-p or --reference-percentages) is required']
        else:
            # Are they numbers?
            conversion_error = False
            total_sum = 0
            self.ref_percs = self.ref_percs.split(",")
            for i in range(len(self.ref_percs)):
                perc = self.ref_percs[i]
                if not re.match(r"\d+\.?\d*", perc):
                    self.errors += [f"{perc} is not a number!"]
                    conversion_error = True
                else:
                    self.ref_percs[i] = float(self.ref_percs[i])
                    total_sum += self.ref_percs[i]

            # Scale them between 0 to 1
            if not conversion_error:
                for i in range(len(self.ref_percs)):
                    self.ref_percs[i] /= total_sum

        # Is haplotype file defined?
        if self.hap_file is None:
            self.errors += ['Haplotype file (-f or --haplotype-file) is required']
        # Does it exist?
        elif not os.path.isfile(self.hap_file):
            self.errors += ["Haplotype file (-f or --haplotype-file) doesn't exist"]

        # Does the TPED directory exist?
        if not os.path.isdir(self.dir):
            self.errors += ["TPED/TFAM directory doesn't exists"]
        else:
            # Do the populations defined exists?
            if self.ref_pops is not None:
                for i in range(len(self.ref_pops)):
                    prefix = os.path.join(self.dir, self.ref_pops[i])
                    if not os.path.isfile(prefix + ".tped"):
                        self.errors += [f"TPED file for {prefix} was not found"]
                    if not os.path.isfile(prefix + ".tfam"):
                        self.errors += [f"TFAM file for {prefix} was not found"]

        # If ref_pops and ref_percs are provided, ensure same amount of values are provided
        if self.ref_pops is not None and self.ref_percs is not None:
            if len(self.ref_pops) != len(self.ref_percs):
                self.errors += [f"Number of reference percentages do not match number of reference populations"]

    def go(self):
        self.read_haplotype_file()
        self.read_ref_pops()
        self.simulate_ind()
        self.write_sim_ind()
        self.write_genomic_fractions()

    def read_haplotype_file(self):
        logging.info("Opening haplotype file")
        last_chrom = None
        last_stop = 0
        with open(self.hap_file, "r") as f:
            for line in f:
                line = line.strip()
                if line == "":
                    continue
                chrom, start, stop, cm = line.split()
                chrom = int(chrom)
                start = int(start)
                stop = int(stop)
                if last_chrom is None:
                    last_chrom = chrom
                if last_chrom != chrom:
                    self.haplotypes.append([last_chrom, last_stop, -1])
                    last_chrom = chrom
                    last_stop = 0
                if start - last_stop != 0:
                    self.haplotypes.append([last_chrom, last_stop, start])
                    logging.debug(f"[HAP] Had to add {last_chrom}:{last_stop}-{start}")
                last_stop = stop
                self.haplotypes.append([last_chrom, start, last_stop])

            self.haplotypes.append([last_chrom, last_stop, -1])
        logging.info("Read haplotype file")

    def read_ref_pops(self):
        logging.info("Reading reference files...")
        for this_pop in self.ref_pops:
            file = os.path.join(self.dir, this_pop)
            self.ref_tfam[this_pop] = []
            this_col_names = []
            # Read TFAM file
            logging.info(f"Reading {file}.tfam")
            with open(file + ".tfam", "r") as f:
                for line in f:
                    family_id, ind_id, father, mother, sex, pheno = line.rstrip().split()
                    self.ref_tfam[this_pop].append(ind_id)
                    logging.debug(f"For pop {this_pop}, added {ind_id}")
                    this_col_names += [ind_id + "_A", ind_id + "_B"]

            logging.info(f"Reading {file}.tped")
            self.ref_tped[this_pop] = pd.read_table(file + ".tped", index_col=None, header=None, sep=' ')
            self.ref_tped[this_pop].columns = ["chr", "snp", "sex", "position"] + this_col_names
            logging.info(f"References read")

    def simulate_ind(self):
        # TODO: Make it generalizable. For now, I am making the simplistic assumption that everyone has the same set of variants
        logging.info(f"Starting simulation")
        self.ref_for_sim_tfam = {}
        self.simulated_tped = self.ref_tped[self.ref_pops[0]].loc[:, ['chr', 'snp', 'sex', 'position']]

        # TODO: parallelize this loop
        for sim_i in range(self.num_inds):
            logging.info(f"Simulating individual number {sim_i}")
            this_trace1 = []
            this_trace2 = []
            genome_avg = {x: 0 for x in self.ref_pops}
            genome_avg["total"] = 0
            pop_to_ind = {}

            # 1 and 2 for each chrom haplotype
            this_ind_tped1 = pd.Series([], dtype=pd.Int8Dtype())
            this_ind_tped2 = pd.Series([], dtype=pd.Int8Dtype())
            for i in range(len(self.haplotypes)):
                random_inds = []  # Individuals that will be used as "parents" for this individual
                # Go through each reference population
                for this_pop in self.ref_pops:
                    # Pull a random individual from the reference population
                    random_ind = random.choice(self.ref_tfam[this_pop])
                    random_inds.append(random_ind)
                    # if this_pop not in self.ref_for_sim_tfam:
                    #     self.ref_for_sim_tfam[this_pop] = {}
                    # self.ref_for_sim_tfam[this_pop][random_ind] = True
                    pop_to_ind[random_ind] = this_pop

                source_ind1 = random.choices(random_inds, weights=self.ref_percs, k=1)[0]
                source_ind2 = random.choices(random_inds, weights=self.ref_percs, k=1)[0]

                logging.debug(f"Individuals pulled for haplotype {i + 1}: {source_ind1} and {source_ind2}")

                source_pop1 = pop_to_ind[source_ind1]
                source_pop2 = pop_to_ind[source_ind2]

                genome_avg[source_pop1] += 1
                genome_avg[source_pop2] += 1
                genome_avg["total"] += 2

                # TODO: randomize chromosome selection (.1_A or .2_A)
                source_ind1 = re.sub(r"\.[12]", ".1_A", source_ind1)
                source_ind2 = re.sub(r"\.[12]", ".2_A", source_ind2)

                # TODO: decide if it's even worth tracking this or if genome average is sufficient
                this_trace1.append(f"{source_pop1}/{source_ind1}")
                this_trace2.append(f"{source_pop2}/{source_ind2}")

                chrom = self.haplotypes[i][0]
                start = self.haplotypes[i][1]
                stop = self.haplotypes[i][2]

                if stop != -1:
                    this_hap1 = self.ref_tped[source_pop1][source_ind1][
                                    self.ref_tped[source_pop1]["chr"] == chrom].iloc[
                                start:stop]
                    this_hap2 = self.ref_tped[source_pop2][source_ind2][
                                    self.ref_tped[source_pop2]["chr"] == chrom].iloc[
                                start:stop]
                else:
                    this_hap1 = self.ref_tped[source_pop1][source_ind1][
                                    self.ref_tped[source_pop1]["chr"] == chrom].iloc[
                                start:]
                    this_hap2 = self.ref_tped[source_pop2][source_ind2][
                                    self.ref_tped[source_pop2]["chr"] == chrom].iloc[
                                start:]
                logging.debug(f"For {chrom}:{start}-{stop}, I pulled {this_hap1.shape[0]}")

                this_ind_tped1 = pd.concat([this_ind_tped1, this_hap1])
                this_ind_tped2 = pd.concat([this_ind_tped2, this_hap2])

            print(this_ind_tped1.shape[0])
            this_ind_tped1.reset_index(drop=True, inplace=True)
            this_ind_tped2.reset_index(drop=True, inplace=True)

            self.simulated_fractions[f"SIM{sim_i}"] = {x: genome_avg[x] * 100 / genome_avg["total"] for x in
                                                       self.ref_pops}

            self.simulated_tfam.append(["SIM", f"SIM{sim_i}.1"])
            self.simulated_tfam.append(["SIM", f"SIM{sim_i}.2"])
            logging.info(f"Adding SIM{sim_i}.1_A")
            self.simulated_tped[f"SIM{sim_i}.1_A"] = this_ind_tped1
            logging.info(f"Adding SIM{sim_i}.1_B")
            self.simulated_tped[f"SIM{sim_i}.1_B"] = this_ind_tped1
            logging.info(f"Adding SIM{sim_i}.2_A")
            self.simulated_tped[f"SIM{sim_i}.2_A"] = this_ind_tped2
            logging.info(f"Adding SIM{sim_i}.2_B")
            self.simulated_tped[f"SIM{sim_i}.2_B"] = this_ind_tped2
            self.trace.append(this_trace1)
            self.trace.append(this_trace2)
        logging.info("Simulation completed")

    def write_sim_ind(self):
        logging.info(f"Writing simulated reference file ref_{self.out_prefix}.tfam")
        # simulated_ref = self.ref_tped[self.ref_pops[0]].loc[:, ['chr', 'snp', 'sex', 'position']]
        # with open(f"ref_{self.out_prefix}.tfam", "w") as ref:
        #     for this_pop in self.ref_for_sim_tfam:
        #         for this_ind in self.ref_for_sim_tfam[this_pop]:
        #             ref.write(f"{this_pop} {this_ind} 0 0 0 -9\n")
        #             this_ind = re.sub(pattern=r"\.[12]", repl="", string=this_ind)
        #             simulated_ref[this_ind + ".1_A"] = self.ref_tped[this_pop][this_ind + ".1_A"]
        #             simulated_ref[this_ind + ".1_B"] = self.ref_tped[this_pop][this_ind + ".1_B"]
        #             simulated_ref[this_ind + ".2_A"] = self.ref_tped[this_pop][this_ind + ".2_A"]
        #             simulated_ref[this_ind + ".2_B"] = self.ref_tped[this_pop][this_ind + ".2_B"]
        # logging.info(f"Writing simulated reference file ref_{self.out_prefix}.tped")
        # simulated_ref.to_csv(f"ref_{self.out_prefix}.tped", header=None, index=None, sep=' ', mode='w')
        tfam_cmd = f"bash -c \" cat {os.path.join(self.dir, self.ref_pops[0])}.tfam "
        tped_cmd = f"bash -c \" paste -d' ' {os.path.join(self.dir, self.ref_pops[0])}.tped "
        for this_pop in self.ref_pops[1:]:
            file = os.path.join(self.dir, this_pop)
            tped_cmd += f" <(cut -f5- -d' ' {file}.tped) "
            tfam_cmd += f" {file}.tfam "
        tped_cmd += f"> ref_{self.out_prefix}.tped \""
        tfam_cmd += f"> ref_{self.out_prefix}.tfam \""

        logging.debug(f"Running: \n{tped_cmd}\n")
        try:
            proc = subprocess.Popen(tped_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                    shell=True, universal_newlines=True)
            proc.communicate()
        except subprocess.CalledProcessError as e:
            logging.error(f"Encountered an error executing the command: ")
            logging.error(f"Error details:")
            logging.error(f"Exit code={e.returncode}")
            logging.error(f"Error message={e.output}")

        logging.debug(f"Running: \n{tfam_cmd}\n")
        try:
            proc = subprocess.Popen(tfam_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                    shell=True, universal_newlines=True)
            proc.communicate()
        except subprocess.CalledProcessError as e:
            logging.error(f"Encountered an error executing the command: ")
            logging.error(f"Error details:")
            logging.error(f"Exit code={e.returncode}")
            logging.error(f"Error message={e.output}")

        logging.info(f"Writing simulated query file query_{self.out_prefix}.tfam")
        with open(f"query_{self.out_prefix}.tfam", "w") as query:
            for entry in self.simulated_tfam:
                query.write(f"{entry[0]} {entry[1]} 0 0 0 -9\n")
        logging.info(f"Writing simulated query file query_{self.out_prefix}.tped")
        self.simulated_tped.to_csv(f"query_{self.out_prefix}.tped", header=False, index=False, sep=' ', mode='w')

    def write_genomic_fractions(self):
        logging.info(f"Writing simulated ancestry fractions to query_{self.out_prefix}.simfrac.tsv")
        with open(f"query_{self.out_prefix}.simfrac.tsv", "w") as f:
            f.write("Ind\t" + "\t".join(self.ref_pops) + "\n")
            for ind in self.simulated_fractions:
                f.write(f"{ind}\t" + "\t".join([str(self.simulated_fractions[ind][x]) for x in self.ref_pops]) + "\n")


if __name__ == '__main__':
    parser = ArgumentParser(prog=PROGRAM_NAME,
                            add_help=False,
                            description=f'''
                            {PROGRAM_NAME} - simulate admixed individuals
                            ''',
                            formatter_class=lambda prog: HelpFormatter(prog, width=120, max_help_position=120))

    parser.add_argument('--help', '-h', '--h', action='store_true', default=False)

    # Add build input arguments
    all_input_group = parser.add_argument_group('Required Input options')
    all_input_group.add_argument('--reference-pops', '-r', required=False, default=None, metavar='POP1,POP2,POP3...',
                                 help='''Comma separated list of prefix for the reference TPED/TFAM input files.
                                 Filenames are expected to be same as population names 
                                 (POP1.tped, POP1.tfam, POP2.tped, POP2.tfam,...)''',
                                 dest="ref_pops", type=str)
    all_input_group.add_argument('--reference-percentages', '-p', required=False, default=None, metavar='60,20,10...',
                                 help='Prefix for the query TPED/TFAM input files (will be normalized to 0 to 1)',
                                 dest="ref_percs", type=str)
    all_input_group.add_argument('--haplotype-file', '-f', required=False, default=None, metavar='<haplotype-file.txt>',
                                 help='Haplotype file containing chr\thap_start\thap_stop',
                                 dest="hap_file", type=str)

    # Add runtime options
    all_runtime_group = parser.add_argument_group('Runtime options')
    all_runtime_group.add_argument('--directory', '-d', required=False, default="./", metavar='<directory path>',
                                   help='Directory containing the reference TPED/TFAM files',
                                   dest="dir", type=str)
    all_runtime_group.add_argument('--num', '-n', required=False, default=100, metavar='<number of sim individuals>',
                                   help='Number of individuals to simulate',
                                   dest="num", type=int)
    all_runtime_group.add_argument('--verbose', '-v', required=False,
                                   help='Verbosity levels', dest="verbose", action='store_true')
    all_runtime_group.add_argument('--threads', '-t', required=False, default=8, metavar='<threads>',
                                   help='Number of threads to use',
                                   dest="threads", type=int)

    # Add build output arguments
    all_output_group = parser.add_argument_group('Output options')
    all_output_group.add_argument('--out-prefix', '-o', required=False, default="out", metavar='<out-prefix>',
                                  type=str, dest="out_prefix",
                                  help="Output prefix for creating TPED/TFAM file")
    all_output_group.add_argument('--log', '-l', required=False, default="run.log", metavar='<log-file>',
                                  type=str, dest="log_file",
                                  help="Log file")

    options, unknown_arguments = parser.parse_known_args()

    if len(unknown_arguments) > 0:
        print(Colors.HEADER)
        print("User specificed unknown arguments: " + " ".join([str(x) for x in unknown_arguments]))
        print("Please see the correct usage below:")
        parser.print_help()
        print(Colors.ENDC)
        sys.exit(1)

    simulate = Simulation(options)
    simulate.validate_options()
    if len(simulate.errors) > 0:
        print(Colors.FAIL)
        print("Error(s) encountered, cannot proceed:")
        print("\n".join(simulate.errors))
        print("\n\nSee usage:")
        print(Colors.ENDC)
        print(Colors.HEADER)
        parser.print_help()
        print(Colors.ENDC)
        sys.exit(1)
    simulate.init_logger()
    simulate.go()
