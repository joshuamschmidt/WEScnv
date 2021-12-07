#!/usr/bin/env python
import argparse
import os
import sys
import textwrap
import mmap
import numpy as np
import pandas as pd
from difflib import SequenceMatcher

'''required and optional argument parser'''

parser = argparse.ArgumentParser(prog='assignBioSex',
                        formatter_class=argparse.RawDescriptionHelpFormatter,
                        description=textwrap.dedent('''\
                        Per Sample MosDepth summaries to bioSex
                        ------------------------------------------------
                        Generates a table of sample + bioSex assignments

                        Column headers are sample names derived from
                        the list of input files and the common file
                        {suffix} e.g. SAMPLE_1{.mosdepth.summary.txt}.
                        Works just as well when data is mean coverage.

                        !!!!!Python needs to be >=3.9!!!!!!
                        '''),
                        add_help=False,
                        epilog="Questions, bugs etc?\njoshmschmidt1@gmail.com\ngithub.com/joshuamschmidt")
parser._action_groups.pop()
required = parser.add_argument_group('required arguments')
optional = parser.add_argument_group('optional arguments')

# Add back help
optional.add_argument(
    '-h',
    '--help',
    action='help',
    default=argparse.SUPPRESS,
    help='show this help message and exit'
)

parser.add_argument('files', metavar='FS', type=str, nargs='+',
                    help='files to convert')

required.add_argument('--suffix', type=str, dest='file_suffix',
                    help='common file suffix')
parser.set_defaults(file_suffix=".mosdepth.summary.txt")

optional.add_argument('--threads', type=int, dest='threads',
                       help='multithreaded mode with int threads:NOT IMPLEMENTED')

optional.add_argument('--out', type=argparse.FileType('w'),
                       dest='out_file',
                       help='output file (or stdout if not set)',
                       default=sys.stdout)

# '''class for counts/coverage files'''
class coverageSummaries():
    def __init__(self, files, file_suffix: str):
        self.files = files
        self.file_suffix = file_suffix
        self.n_files = len(self.files)
        self.summaries = []
        self.make_summaries()
    def make_summaries(self):
        samples_info = []
        for i, file in enumerate(self.files):
            sample_name = file.removesuffix(self.file_suffix)
            coverageDF = pd.read_csv(file, sep='\t')
            target_mean = coverageDF.query('chrom == "total_region"').iloc[0]['mean']
            autosomes = np.asarray(coverageDF[coverageDF.chrom.str.contains('_region')][:22]['mean'])
            autosome_mean = np.mean(np.asarray(coverageDF[coverageDF.chrom.str.contains('_region')][:22]['mean']))
            autosome_sd = np.std(np.asarray(coverageDF[coverageDF.chrom.str.contains('_region')][:22]['mean']))
            x_mean = coverageDF.query('chrom == "chrX_region"').iloc[0]['mean']
            mean_xToA = np.mean(x_mean / autosomes)
            sd_xToA = np.std(x_mean / autosomes)
            y_mean = coverageDF.query('chrom == "chrY_region"').iloc[0]['mean']
            mean_yToA = np.mean(y_mean / autosomes)
            sd_yToA = np.std(y_mean / autosomes)
            yToX = y_mean / x_mean
            bioSex='UNASSIGNED'
            # Assign sex
            if mean_xToA >= 0.45 and mean_xToA <= 0.65:
                if mean_yToA > 0.2:
                    bioSex="MALE:XY"
                if mean_yToA < 0.02:
                    bioSex="FEMALE:XO"
            elif mean_xToA >= 0.85 and mean_xToA <=1.25:
                if mean_yToA <0.02:
                    bioSex="FEMALE:XX"
                elif mean_yToA > 0.2:
                     bioSex="MALE:XXY"
            elif mean_xToA > 1.25:
                if mean_yToA <0.02:
                    bioSex="FEMALE:XXX"
            summary = [sample_name, bioSex, autosome_mean, autosome_sd, x_mean, mean_xToA, sd_xToA, y_mean, mean_yToA, sd_yToA]
            samples_info.append(summary)
        df = pd.DataFrame(samples_info)
        df.columns = ['Sample_ID','bioSex','m_autosome','sd_autosome','X','m_X/A','sd_X/A','Y','m_Y/A','sd_Y/A']
        self.summaries = df

def main():
    args = parser.parse_args()
    assert all(args.file_suffix in file for file in args.files), "files have different suffixes"
    samples = coverageSummaries(args.files, file_suffix=args.file_suffix)
    samples.summaries.to_csv(path_or_buf=args.out_file, sep='\t', encoding='utf-8',index=False, float_format='%3f',header=True)

if __name__ == '__main__':
    main()

