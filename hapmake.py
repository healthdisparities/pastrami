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

from argparse import ArgumentParser, HelpFormatter
import sys
import os
import pandas as pd
import pathos.multiprocessing as mp

PROGRAM_NAME = "hapmake.py"


class Colors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'


class HapMake:
    def __init__(self, opts):
        self.min_snps = opts.min_snps
        self.max_snps = opts.max_snps
        self.max_rate = opts.min_snps
        self.map_dir = opts.map_dir
        self.hap_file = opts.hap_file
        self.threads = opts.threads
        self.errors = []

    def validate_options(self):
        pass

    def go(self):
        function = getattr(self, "make_haplotypes")
        function()

    def process_hapmap_file(self, chrom: int):
        haplotypes = []
        map_data = []
        with open(os.path.join(self.map_dir, f"chr{chrom}.map"), "r") as f:
            for line in f:
                (position, cmorgan, snp) = line.rstrip().split("\t")
                map_data.append([int(position), float(cmorgan), snp])

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
                haplotypes.append([str(chrom), str(left_snp - 1), str(right_snp - 1),
                                   str(map_data[right_snp - 2][1] - map_data[left_snp][1])])
        return haplotypes

    def make_haplotypes(self):
        pool = mp.Pool(processes=self.threads)
        results = pool.map(self.process_hapmap_file, range(1, 23))
        with open(self.hap_file, "w") as f:
            for chrom_hap in results:
                for row in chrom_hap:
                    f.write("\t".join(row) + "\n")


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
    hapmake_parser = subparsers.add_parser('hapmake', help='Create haplotypes',
                                           formatter_class=lambda prog: HelpFormatter(prog, width=120,
                                                                                      max_help_position=120))
    hapmake_parser.set_defaults(sub_command='hapmake')

    hapmake_input_group = hapmake_parser.add_argument_group('Input options')
    hapmake_input_group.add_argument('--map-dir', required=True, default=None, metavar='maps/', type=str,
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
    hapmake_input_group.add_argument('--threads', required=False, default=4, metavar='Threads', type=int,
                                     help='Number of concurrent threads (default: %(default)s)', dest="threads")
    hapmake_input_group = hapmake_parser.add_argument_group('Output options')
    hapmake_input_group.add_argument('--hap-file', required=False, default=None, metavar='out.haplotypes', type=str,
                                     help='Output file containing haplotypes (default: %(default)s)', dest="hap_file")

    options, unknown_arguments = parser.parse_known_args()

    if options.help or "sub_command" not in options:
        print(Colors.HEADER)
        parser.print_help()
        print(Colors.ENDC)
        sys.exit(0)

    if options.sub_command == 'hapmake' and options.help:
        print(Colors.HEADER)
        hapmake_parser.print_help()
        print(Colors.ENDC)
        sys.exit(0)

    # Start the analysis
    analysis = HapMake(options)

    analysis.validate_options()

    if len(analysis.errors) > 0:
        print(Colors.HEADER)
        if options.sub_command == 'hapmake':
            hapmake_parser.print_help()
        print(Colors.ENDC)
        print(Colors.FAIL + '\n\nErrors:')
        [print(i) for i in analysis.errors]
        print(Colors.ENDC)
        sys.exit()

    # If we're still good, start the actual analysis
    analysis.go()
