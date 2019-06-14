#! /usr/bin/python

import numpy
import subprocess
import argparse
from multiprocessing import Pool
import parmap

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

def kmer_hash_gen(read, k_size):
    seq_len = len(read)
    half_read = seq_len / 2
    kmer_start_pos = half_read - (k_size / 2)
    kmer_end_pos = half_read + (k_size / 2)
    return read[kmer_start_pos:kmer_end_pos]

def kmer_convert_to_bit(kmer_list):
    # A = 00
    # C = 01
    # G = 10
    # T = 11
    hash_kmers = []
    for letter in kmer_list:
        if letter == 'A':
            hash_kmers.append(00)
        elif letter == 'C':
            hash_kmers.append(01)
        elif letter == 'G':
            hash_kmers.append(10)
        elif letter == 'T':
            hash_kmers.append(11)
    return hash_kmers


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

                    # subprocess.call(['./kmer_mapper.py', '--read_seq' , seq])

    p = Pool(int(args.threads))


    #print(p.map(kmer_len_gen, read_list))
    kmer_sizes = (p.map(kmer_len_gen, read_list))
    # print(kmer_sizes)

    # the pm_pbar=True requres a pip install tqdm
    kmer_seqs = parmap.starmap(kmer_hash_gen, zip(read_list , kmer_sizes), pm_pbar=True)
    print(kmer_seqs)

    kmer_convert = (p.map(kmer_convert_to_bit, kmer_seqs))
    print(kmer_convert)








if __name__ == '__main__':
    main()
