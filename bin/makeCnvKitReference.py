#!/usr/bin/env python
import argparse
import os
import sys
import numpy as np
import pandas as pd


'''required and optional argument parser'''

parser = argparse.ArgumentParser(prog='makeCnvKitReference',
                        formatter_class=argparse.RawDescriptionHelpFormatter,
                        description=textwrap.dedent('''\
                        Make CNVKIT reference file set per sample
                        --------------------------------------------
                        Generates a CNVKIT reference for a sample.

                        Requires sample_id, file of matched reference
                        ids and globs of targetcoverage and antitargetcoverage
                        files

                        Runs the CNVKIT reference script
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

parser.add_argument('--sample_id', type=str, dest='sample_id',
                    help='target sample id')

required.add_argument('--matched_ref', type=str, dest='matched_ref',
                    help='file of matched ref samples')

required.add_argument('--refFasta', type=str, dest='refFasta',
                    help='reference genome fasta file')

parser.add_argument('--target_files', type=str, nargs='+', dest='target_files',
                    help='target coverage files')

parser.add_argument('--anti_target_files', type=str, nargs='+', dest='anti_target_files',
                    help='antittarget coverage files')


class coverageReference():
    def __init__(self, sample_id, matched_ref_file, target_files, anti_target_files):
        self.sample = sample_id
        self.matched_ref_file = matched_ref_file
        self.target_files = target_files
        self.anti_target_files = anti_target_files
    
    def make_matched_ref_names(self):
    	with open(self.matched_ref_file) as file:
    		ref_lines = [line.rstrip() for line in file]
        ref_lines = ref_lines[1:]
        ref_lines.sort()
        self.matched_target_files =  [ref + '.targetcoverage.cnn' for ref in ref_lines]
        self.matched_anti_target_files = [ref + '.antitargetcoverage.cnn' for ref in ref_lines]
        self.matched_all_files = self.matched_target_files + self.matched_anti_target_files
        self.matched_all_files.sort()
        self.matched_all_files_string= (' ').joinelf.matched_all_files
       
       # actually files should be staged by nextflow so just suffix is sufficient 
       #self.filtered_target_files = [file for file in self.target_files if os.path.basename(file) in self.matched_target_files]
       #self.filtered_anti_target_files = [file if os.path.basename(file) in self.matched_anti_target_files for file in self.anti_target_files]


# '''main''' 
def main():
    args = parser.parse_args()
    refObject=coverageReference(sample_id=args.sample_id, matched_ref_file=args.ref_panel, target_files=args.target_files, anti_target_files=args.anti_target_files)
    outfile=args.sample_id + '.Reference.cnn'
    command = 'cnvkit.py reference ' + refObject.matched_all_files_string + ' -f ' + args.refFasta + ' -o ' + outfile
    os.system(command)
if __name__ == '__main__':
    main()
