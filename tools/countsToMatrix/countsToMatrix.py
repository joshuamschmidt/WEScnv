#!/usr/bin/env python
import argparse
import os
import sys
import textwrap
import mmap
import numpy as np
import pandas as pd
from difflib import SequenceMatcher
#from multiprocessing.dummy import Pool as ThreadPool

'''required and optional argument parser'''

parser = argparse.ArgumentParser(prog='countsToMatrix',
                        formatter_class=argparse.RawDescriptionHelpFormatter,
                        description=textwrap.dedent('''\
                        Per Sample Target BED Read Counts to Matrix
                        --------------------------------------------
                        Generates a Matrix of read counts with dims 
                        n_targets * n_samples to stdout (or outfile 
                        if given).
                        
                        Column headers are sample names derived from
                        the list of input files and the common file
                        {suffix} e.g. SAMPLE_1{-counts.txt}.
                        Works just as well when data is mean coverage.
                        
                        Optional: 
                        A n-col BED file --> first n-cols matrix.
                        Convert counts to FKPM (requires BED).
                        Multithread.
                        
                        Warnings:
                        * The BED file can have a header, so long
                        as col1 is one of chr, Chr, or CHR.
                        
                        !!!!! Aside from a basic comparison of the first 
                        line of the first count file to the first
                        line of the bed file, order consistency is NOT
                        checked. It is up to the user to correctly sort
                        all inputs !!!!!!
                        
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

optional.add_argument('--bed', type=str, dest='bed_file',
                       help='bed file with n cols, to form first n cols of output matrix',
                       default=None)

optional.add_argument('--merge-bed', dest='merge_bed',
                       help='combine bed file and counts',
                       action='store_true')
parser.set_defaults(merge_bed=False)

optional.add_argument('--fpkm', dest='fpkm',
                       help='convert counts to FPKM',
                       action='store_true')
parser.set_defaults(fpkm=False)

optional.add_argument('--column', type=int, dest='count_column',
                       help='column number of the count data. Default is 5')
parser.set_defaults(count_column=5)

optional.add_argument('--threads', type=int, dest='threads',
                       help='multithreaded mode with int threads:NOT IMPLEMENTED')

optional.add_argument('--out', type=argparse.FileType('w'),
                       dest='out_file',
                       help='output file (or stdout if not set)',
                       default=sys.stdout)
                       
# '''class for counts/coverage files'''
class Counts():
    def __init__(self, files, count_column: int, file_suffix: str):
        self.files = files
        self.file_suffix = file_suffix
        self.count_column = count_column
        self.n_files = len(self.files)
        self.sample_names = []
        self.n_rows = None
        self.first_key=None
        self.get_n_rows()
        self.count_array = np.empty([self.n_rows,self.n_files], dtype=np.int64)
        self.fill_count_array()
        #self.sample_names = [file.removesuffix(file_suffix) for file in self.files]
    
    def get_n_rows(self):
        file=self.files[0]
        first_fileDF=pd.read_csv(file, sep='\t',header=None)
        lines = len(first_fileDF.index)
        self.n_rows = lines
    
    def fill_count_array(self):
        count_index = self.count_column - 1
        for i, file in enumerate(self.files):
            countDF = pd.read_csv(file, sep='\t',header=None)
            count_array=countDF.iloc[:,count_index].values
            #check dtype of count file. change type of self.count_array
            if i==0:
                self.count_array = self.count_array.astype(count_array.dtype)
                self.first_key = str(countDF.iloc[0,0]) + "_" + str(countDF.iloc[0,1])
            assert self.n_rows == np.size(count_array), f"File {file} has an unexpected number of rows. Expected {self.n_rows}. Observed {np.size(count_array)}"
            self.count_array[:,i] = count_array
            self.sample_names.append(file.removesuffix(self.file_suffix))
    
    def convert_to_fkpm(self,bedObject):
        #with open(self.files[0]) as f:
        #    first_line = f.readline().strip()
        #    first_array = first_line.split("\t")
        #    self.first_key = str(first_array[0]) + "_" + str(first_array[1])
        assert self.first_key == bedObject.first_key, f"count files are not sorted in the same order as the bed file. File key: {self.first_key}, Bed key: {bedObject.first_key}"
        assert self.n_rows == bedObject.n_rows, f"count files and bed file differ in nummber of features. File: {self.n_rows}, Bed: {bedObject.n_rows}" 
        perM_scaling_factors = self.count_array.sum(axis=0) / 1e6
        fpm_scaled_count_array = self.count_array / perM_scaling_factors
        fpkm_array = np.array(fpm_scaled_count_array / bedObject.kb_lengths[:,None])
        self.count_array = fpkm_array
    
    def array_to_df(self):
        self.count_array = pd.DataFrame(self.count_array)
        self.count_array.columns = self.sample_names
    
    def append_bed(self, bed):
       assert isinstance(self.count_array, pd.DataFrame), "you need to convert to DF before this! (array_to_df)"
       self.count_array = pd.concat([bed.bed_data,self.count_array],axis=1)


# '''class for Bedfiles'''
class BedFile():
    def __init__(self, bedfile):
        self.bed_file = bedfile
        self.read_bed()
        self.n_rows = len(self.bed_data.index)
        self.first_key = str(self.bed_data.iloc[0,0]) + "_" + str(self.bed_data.iloc[0,1])
        self.kb_lengths = (self.bed_data.iloc[:,2].values - self.bed_data.iloc[:,1].values)/1000
    
    def read_bed(self):
        first_line = pd.read_csv(self.bed_file, sep='\t', header=None, nrows=1)
        first_value = first_line.iloc[0,0]
        if first_value=='chr' or first_value=='Chr' or first_value=='CHR':
             self.bed_data = pd.read_csv(self.bed_file, sep='\t',dtype=object)
             n1, n2 = list(self.bed_data.columns)[1:3]
             self.bed_data = self.bed_data.astype({n1: int, n2: int})
        else:
            self.bed_data = pd.read_csv(self.bed_file, sep='\t',header=None,dtype=object)
            col_names = list(self.bed_data.columns)
            col_names = ["user_"+str(name+1) for name in  col_names]
            col_names[0:3] = 'chr', 'start', 'end'
            self.bed_data.columns = col_names

# '''main''' 
def main():
    args = parser.parse_args()
    match_string = SequenceMatcher(None, args.files[0], args.files[1]).find_longest_match(0, len(args.files[0]), 0, len(args.files[1]))
    file_suffix = args.files[0][match_string.a: match_string.a + match_string.size]
    assert all(file_suffix in file for file in args.files), "files have different suffixes"
    counts=Counts(files=args.files,count_column=args.count_column, file_suffix=file_suffix)
    if args.bed_file is not None:
        bed = BedFile(args.bed_file)
    
    if args.fpkm:
        assert args.bed_file is not None, "Conversion to FKPM requires a BED file"
        counts.convert_to_fkpm(bed)
    
    counts.array_to_df()
    if args.merge_bed:
        counts.append_bed(bed)
    # '%4g' should leave ints as ints, but print floats to four decimals.....
    counts.count_array.to_csv(path_or_buf=args.out_file, sep='\t', encoding='utf-8',index=False, float_format='%4g',header=True)

if __name__ == '__main__':
    main()
