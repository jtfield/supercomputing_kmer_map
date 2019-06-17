#! /usr/bin/python

import numpy
import subprocess
import argparse
from multiprocessing import Pool
import parmap

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--read_file_1')
    parser.add_argument('--genome_file')
    parser.add_argument('--threads')
    # parser.add_argument('--read_file_2')
    return parser.parse_args()

# get length of reads and find appropriate kmer lengths
def kmer_len_gen(s):
    seq_length = len(s)
    kmer_len = seq_length / 10
    return kmer_len


# get kmer sequences in list form
def kmer_hash_gen(read, k_size):
    seq_len = len(read)
    half_read = seq_len / 2
    kmer_start_pos = half_read - (k_size / 2)
    kmer_end_pos = half_read + (k_size / 2)
    return read[kmer_start_pos:kmer_end_pos]

# convert kmers into hash forms
def kmer_convert_to_bit(kmer_list):
    # A = 00
    # C = 01
    # G = 10
    # T = 11
    hash_kmers = []
    for letter in kmer_list:
        if letter == 'A':
            hash_kmers.append('00')
        elif letter == 'C':
            hash_kmers.append('01')
        elif letter == 'G':
            hash_kmers.append('10')
        elif letter == 'T':
            hash_kmers.append('11')
    return hash_kmers


# convert genome or msa into hash form
def gen_convert_to_bit(genome_file):
    hash_gen = []
    full_genome = genome_file
    with open(full_genome) as sequence:
        for line in sequence:
            line = line.upper()
            if ">" not in line:
                for letter in line:
                    if letter == 'A':
                        hash_gen.append('00')
                    elif letter == 'C':
                        hash_gen.append('01')
                    elif letter == 'G':
                        hash_gen.append('10')
                    elif letter == 'T':
                        hash_gen.append('11')
    return hash_gen

# use the max kmer size to chunk up the genome into kmer sized chunks
# these chunks should move up the genome 1 nucleotide hash at a time
def split_genome(genome_hash_list, max_kmer_len):
    genome_hash_chunk_array = []
    nuc_pos = 0
    while max_kmer_len != len(genome_hash_list):
        window = genome_hash_list[nuc_pos:max_kmer_len]
        genome_hash_chunk_array.append(window)
        nuc_pos+=1
        max_kmer_len+=1
    return genome_hash_chunk_array


# convert list of lists to dictionary of lists with positions as keys
def list_to_dict(list_of_lists):
    list_count = 0
    list_dict = {}
    for item in list_of_lists:
        list_count+=1
        list_dict[list_count] = item
    return list_dict


# compare the kmer-bit string from the reads to each kmer-bit string from the genome
# def kmer_matcher(read_kmer_dict, genome_kmer_dict):
#     matching_kmers = cmp(read_kmer_dict, genome_kmer_dict)
#     return matching_kmers


# create new dictionaries containing a single kmer from the read_kmers and the entire dictionary of genome kmers
# this will be inefficient but i don't know how to do it any other way yet
def dict_combiner(read_kmer_dict, genome_kmer_dict):
    new_combined_dict = {}
    new_combined_list = []
    combined_list_count = 0
    for key1, kmer1 in read_kmer_dict.items():
        combined_list_count+=1
        temp_combined_dict = {}
        temp_read_kmer_dict = {}
        #temp_genome_kmer_dict = {}
        read_kmer_key = "kmer_num_" + str(key1)
        temp_read_kmer_dict[read_kmer_key] = kmer1

        genome_kmer_counter = "genome_kmer_" + str(combined_list_count)
        read_kmer_counter = "read_kmer_" + str(combined_list_count)
        temp_combined_dict[read_kmer_counter] = temp_read_kmer_dict
        temp_combined_dict[genome_kmer_counter] = genome_kmer_dict

        # new_combined_list.append(temp_read_kmer_dict)
        # new_combined_list.append(genome_kmer_dict)
        # combined_list_count+=1
        #new_combined_dict[combined_list_count] = new_combined_list
        combined_read_genome_dict = "combined_read_genome_dict_" + str(combined_list_count)
        new_combined_dict[combined_read_genome_dict] = temp_combined_dict
    return new_combined_dict

# def kmer_search(kmer_dicts):
#     list_counter = 0
#     read_kmer = []
#     for item in kmer_dicts:
#         list_counter+=1
        # print(item)
        # print(list_counter)
        # print('\n')
        # if list_counter == 1:
        #     print(item)
        #     read_kmer.append(item)
        #     if list_counter > 1:
        #         return item



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


    #PRODUCES LIST OF LENGTHS OF VARIOUS KMERS IN ORDER FROM FASTQ FILE
    kmer_sizes = (p.map(kmer_len_gen, read_list))
    max_kmer_size = max(kmer_sizes)


    # PRODUCES SEQUENCES FOR KMERS OF APPROPRIATE LENGTHS
    # the pm_pbar=True requres a pip install tqdm
    kmer_seqs = parmap.starmap(kmer_hash_gen, zip(read_list , kmer_sizes), pm_pbar=True)
    #print(kmer_seqs)

    # CONVERTS KMERS INTO HASH VALUES (A = 00, C = 01, G = 10, T = 11)
    kmer_convert = (p.map(kmer_convert_to_bit, kmer_seqs))
    #print(kmer_convert)

    # CONVERT LIST OF KMERS TO A dictionary
    convert_list_to_dict = list_to_dict(kmer_convert)
    # print(convert_list_to_dict)

    # CONVERT GENOME OR MSA TO HASH
    gen_hasher = gen_convert_to_bit(args.genome_file)
    # print(gen_hasher)

    # CONVERT HASH GENOME TO LIST OF KMERS FOR MAPPING
    hash_gen_chunker = split_genome(gen_hasher, max_kmer_size)


    # CONVERT GENOME LIST OF KMERS TO DICT
    convert_genome_list_to_dict = list_to_dict(hash_gen_chunker)
    # print(convert_genome_list_to_dict)


    # ATTEMPTING TO COMPARE EACH READ KMER TO EACH GENOME KMER IN A PARALLEL PROCESS
    combined_dicts = dict_combiner(convert_list_to_dict, convert_genome_list_to_dict)
    print(combined_dicts)

    # kmer_match_search = kmer_search(combined_dicts[1])
    # print(kmer_match_search)

    # kmer_match_search = kmer_search(combined_dicts[2])
    # print(kmer_match_search)


if __name__ == '__main__':
    main()
