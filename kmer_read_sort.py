#! /usr/bin/python

import numpy
import subprocess
import argparse
from multiprocessing import Pool

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--read_file_1')
    #parser.add_argument('--read_file_2')
    #parser.add_argument('--genome_file')
    parser.add_argument('--threads')
    return parser.parse_args()

def kmer_len_gen(s):
    seq_length = len(s)
    kmer_len = seq_length / 10
    return kmer_len


def main():
    args = parse_args()
# GET READ LEN TO EXTRAPOLATE KMER LENGTH
    line_count = 0
    header = ''
    seq = ''
    seperate = ''
    qual = ''
    read_line_count = 0
    read_list = []
    qual_list = []
    with open(args.read_file_1) as read1:
        for line in read1:
            # read in read file and count the number of reads in the entire file
            # this information will be used for the location of the read within the file
            #
            if len(line) > 0:
                line_count+=1
                read_line_count+=1

                if line_count == 1:
                    header = line
                elif line_count == 2:
                    seq = line
                    read_list.append(seq)
                elif line_count == 3:
                    seperate = line
                elif line_count == 4:
                    qual = line
                    qual_list.append(qual)
                    line_count = 0

                    # subprocess.call(['./kmer_mapper.py', '--read_seq' , seq])

    p = Pool(args.threads)
    print(p.map(kmer_len_gen, read_list))








if __name__ == '__main__':
    main()
