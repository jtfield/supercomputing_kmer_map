#!/usr/bin/env python

import sys
import argparse

def parse_args():
    parser = argparse.ArgumentParser()
    # parser.add_argument('--read_file_1')
    parser.add_argument('--genome_file')
    # parser.add_argument('--threads')
    # parser.add_argument('--read_file_2')
    return parser.parse_args()


# get length of reads and find appropriate kmer lengths
def kmer_len_gen(s):
    kmer_lengths = 4
    return kmer_lengths

# def kmer_hash_gen(read, k_size):
# #################################################################
# # code for getting kmer from the beginning of read
#     kmer_list = []
#     for sing_read in read:
#         seq_len = len(sing_read)
#         # half_read = seq_len / 2
#         # kmer_start_pos = half_read - (size / 2)
#         # kmer_end_pos = half_read + (size / 2)
#         kmer = sing_read[:int(k_size)]
#         kmer_list.append(kmer)
#     return(kmer_list)

def genome_reader(genome_file):
    gen_seq = ""
    with open(genome_file) as seq:
        for line in seq:
            #line = line.strip()
            line = line.strip('\n')
            #line = line.strip(' ')
            gen_seq = gen_seq + line
    return gen_seq

# def kmer_convert_to_bit(full_kmer_list):
#     # A = 00
#     # C = 01
#     # G = 10
#     # T = 11
#     hash_kmers = []
#     for seq in full_kmer_list:
#         #print("converting kmers to hash")
#         kmer_seq = []
#         for letter in seq:
#             if letter == 'A':
#                 kmer_seq.append('00')
#             elif letter == 'C':
#                 kmer_seq.append('01')
#             elif letter == 'G':
#                 kmer_seq.append('10')
#             elif letter == 'T':
#                 kmer_seq.append('11')
#         hash_kmers.append(kmer_seq)
#     return hash_kmers

# # convert genome or msa into hash form
# def gen_convert_to_bit(genome_file):
#     hash_gen = []
#     full_genome = genome_file
#     with open(full_genome) as sequence:
#         for line in sequence:
#             line = line.upper()
#             if ">" not in line:
#                 for letter in line:
#                     if letter == 'A':
#                         hash_gen.append('00')
#                     elif letter == 'C':
#                         hash_gen.append('01')
#                     elif letter == 'G':
#                         hash_gen.append('10')
#                     elif letter == 'T':
#                         hash_gen.append('11')
#     return hash_gen

# use the max kmer size to chunk up the genome into kmer sized chunks
# these chunks should move up the genome 1 nucleotide hash at a time
# def split_genome(genome_hash_list, max_kmer_len):
#     genome_hash_chunk_array = {}
#     nuc_pos = 0
#     # with open(genome_hash_list) as sequence:
#     for line in genome_hash_list:
#         if ">" not in line:
#             while max_kmer_len < len(line):
#                 print(line)
#                 # window = line[nuc_pos:max_kmer_len]
#                 # window = window.upper()
#                 # window_hash = hash(window)
#                 # nuc_pos+=1
#                 max_kmer_len+=1
#                 #
#                 # if not window_hash in genome_hash_chunk_array:
#                 #     # add kmer to the hash table and add
#                 #     genome_hash_chunk_array[window_hash] = [nuc_pos]
#                 #
#                 # elif window_hash in genome_hash_chunk_array:
#                 #     genome_hash_chunk_array[window_hash].append(nuc_pos)
#             return genome_hash_chunk_array


def split_genome(genome_hash_list, max_kmer_len):
    genome_hash_chunk_array = {}
    nuc_pos = 0
    while nuc_pos < len(genome_hash_list):
        window = genome_hash_list[nuc_pos:max_kmer_len+nuc_pos]
        #window = line[nuc_pos:max_kmer_len]
        window = window.upper()
        window_hash = hash(window)
        nuc_pos+=1
        #
        if not window_hash in genome_hash_chunk_array:
            # add kmer to the hash table and add
            genome_hash_chunk_array[window_hash] = [nuc_pos]

        elif window_hash in genome_hash_chunk_array:
            genome_hash_chunk_array[window_hash].append(nuc_pos)
    return genome_hash_chunk_array



def main():
    args = parse_args()

    # # GET READ LEN TO EXTRAPOLATE KMER LENGTH
    # line_count = 0
    # header = ''
    # seq = ''
    # seperate = ''
    # qual = ''
    # read_list = []
    # qual_list = []
    # with open(args.read_file_1) as read1:
    #     for line in read1:
    #         # read in read file and count the number of reads in the entire file
    #         # this information will be used for the location of the read within the file
    #         #
    #         if len(line) > 0:
    #             line_count+=1
    #
    #             if line_count == 1:
    #                 header = line
    #             elif line_count == 2:
    #                 seq = line
    #                 read_list.append(seq)
    #             elif line_count == 3:
    #                 seperate = line
    #             elif line_count == 4:
    #                 qual = line
    #                 qual_list.append(qual)
    #                 line_count = 0


#PRODUCES LIST OF LENGTHS OF VARIOUS KMERS IN ORDER FROM FASTQ FILE
    # kmer_sizes1 = kmer_len_gen(read_list)
    #print(kmer_sizes1)
    # max_kmer_size1 = max(kmer_sizes1)
    # print(max_kmer_size1)

    # PRODUCES SEQUENCES FOR KMERS OF APPROPRIATE LENGTHS
    # kmer_seqs1 = kmer_hash_gen(read_list , kmer_sizes1)
    # #print(kmer_seqs1)
    #
    #
    # # CONVERTS KMERS INTO HASH VALUES (A = 00, C = 01, G = 10, T = 11)
    # kmer_convert1 = kmer_convert_to_bit(kmer_seqs1)
    #print(kmer_convert1)

    # CONVERT GENOME OR MSA TO HASH
    #gen_hasher1 = gen_convert_to_bit(args.genome_file)
    #print(gen_hasher1)
    kmer_sizes1 = 9

    genome = genome_reader(args.genome_file)
    #print(genome)

    # genome = "ACGTACTTTAAACTACTAGGGATACTAG"


    # CONVERT HASH GENOME TO LIST OF KMERS FOR MAPPING
    hash_gen_chunker1 = split_genome(genome, kmer_sizes1)
    print(hash_gen_chunker1)
    print(sys.getsizeof(hash_gen_chunker1))

    # # MAP ALL READ KMERS TO ALL GENOME KMERS AND CONSTRUCT DICT
    # # DICT CONTAINS INFO ON BEST MAP LOCATION WITHIN THE GENOME KMERS
    # kmer_search1 = kmer_matcher(kmer_convert1, hash_gen_chunker1)
    # # print(kmer_search1)
    #
    # # CONVERT READS TO BYTES
    # read_convert = kmer_convert_to_bit(read_list)
    # # print(read_convert)
    #
    # #
    # kmer_and_read_matcher = read_kmer_matcher(kmer_search1, read_convert, gen_hasher1)
    #print(kmer_and_read_matcher)

if __name__ == '__main__':
    main()

# import hashlib
# import sys

# # hash for integer unchanged
# print('Hash for 181 is:', hash(181))
#
# # hash for decimal
# print('Hash for 181.23 is:',hash(181.23))
#
# # hash for string
# print('Hash for Python is:', hash('Python'))
#
# print('Hash for A is:', hash('A'))
# print(sys.getsizeof(hash('ACGT')))
#
# x = 2
# print(sys.getsizeof(x))

# print(sys.getsizeof('A'))
# print(sys.getsizeof('00'))
#
# hash_object = hashlib.md5(b'Hello World')
# print(sys.getsizeof(hash_object.hexdigest()))




# s = 'ACGTACGT'
#
#
# def string_hash(seq):
#
#     for item in seq:
#         x = hash(item)
#         return x
#
# spot = string_hash(s)
# print(spot)
# print(sys.getsizeof(spot))
