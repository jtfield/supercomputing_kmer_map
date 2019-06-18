#!/usr/bin/env python


import argparse


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--read_file_1')
    # parser.add_argument('--genome_file')
    # parser.add_argument('--threads')
    # parser.add_argument('--read_file_2')
    return parser.parse_args()


# get length of reads and find appropriate kmer lengths
def kmer_len_gen(s):
    kmer_lengths = []
    for single_read in s:
        seq_length = len(single_read)
        kmer_len = seq_length / 10
        kmer_lengths.append(kmer_len)
    return kmer_lengths



def main():
    args = parse_args()
    # GET READ LEN TO EXTRAPOLATE KMER LENGTH
    line_count = 0
    header = ''
    seq = ''
    seperate = ''
    qual = ''
    read_list = []
    qual_list = []
    with open(args.read_file_1) as read1:
        for line in read1:
            # read in read file and count the number of reads in the entire file
            # this information will be used for the location of the read within the file
            #
            if len(line) > 0:
                line_count+=1

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


#PRODUCES LIST OF LENGTHS OF VARIOUS KMERS IN ORDER FROM FASTQ FILE
    kmer_sizes = kmer_len_gen(read_list)
    print(kmer_sizes)
    max_kmer_size = max(kmer_sizes)
    print(max_kmer_size)

if __name__ == '__main__':
    main()
