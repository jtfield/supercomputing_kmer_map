#! /usr/bin/python

import numpy
import sys
import argparse

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--read_file_1')
    #parser.add_argument('--read_file_2')
    #parser.add_argument('--genome_file')
    return parser.parse_args()


def main():
    args = parse_args()
# GET READ LEN TO EXTRAPOLATE KMER LENGTH
    line_count = 0
    header = ''
    seq = ''
    seperate = ''
    qual = ''
    read_line_count = 0
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
                elif line_count == 3:
                    seperate = line
                elif line_count == 4:
                    qual = line
                    line_count = 0
                    # seq_length = len(seq)
                    #print(seq_length)
                    # kmer_len = seq_length / 10
                    # print(kmer_len)







if __name__ == '__main__':
    main()
