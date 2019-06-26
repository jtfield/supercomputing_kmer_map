#!/usr/bin/env python

import sys
import argparse


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--read_file_1')
    parser.add_argument('--genome_file')
    parser.add_argument('--kmer_size')
    # parser.add_argument('--read_file_2')
    return parser.parse_args()


# get kmer sequences in list form
def kmer_hash_gen(reads, kmer_size):
#################################################################
# code for getting kmer from the beginning of read
    #size = args.kmer_size
    kmer_table = {}
    read_count = 0
    kmer_read_position = 0
    for sing_read in reads:
        read_count+=1
        # seq_len = len(sing_read)
        # half_read = seq_len / 2
        # kmer_start_pos = half_read - (size / 2)
        # kmer_end_pos = half_read + (size / 2)
        sing_read = sing_read.upper()
        kmer = sing_read[:int(kmer_size)]
        kmer_hash = hash(kmer)
        kmer_table[kmer_hash] = read_count
    return(kmer_table)


# chunks each complete read into kmers, hashes them and adds them to a hash table
# def full_read_kmer_hash_gen(reads, kmer_size):
#     kmer_table = {}
#     read_count = 0
#     for sing_read in reads:
#         read_count+=1
#         # """Yield successive n-sized chunks from l."""
#         parts = [sing_read[i:i+int(kmer_size)] for i in range(0, len(sing_read), int(kmer_size))]
#         kmer_position = 0
#         read_kmer_chunks_table = {}
#         for kmer in parts:
#             kmer = kmer.strip('\n')
#             kmer = hash(kmer)
#             kmer_position+=1
#             read_kmer_chunks_table[kmer_position] = kmer
#         kmer_table[read_count] = read_kmer_chunks_table
#     return kmer_table

def full_read_kmer_hash_gen(reads, kmer_size):
    kmer_table = {}
    read_count = 0
    for sing_read in reads:
        read_count+=1
        # """Yield successive n-sized chunks from l."""
        parts = [sing_read[i:i+int(kmer_size)] for i in range(0, len(sing_read), int(kmer_size))]
        kmer_position = 0
        read_kmer_chunks_table = {}
        for kmer in parts:
            kmer = kmer.strip('\n')
            kmer = hash(kmer)
            kmer_position+=1
            read_kmer_chunks_table[kmer] = kmer_position
        kmer_table[read_count] = read_kmer_chunks_table
    return kmer_table



# read in genome and strip out newlines to put all sequences on the same line
def genome_reader(genome_file):
    gen_seq = ""
    with open(genome_file) as seq:
        for line in seq:
            #line = line.strip()
            line = line.strip('\n')
            #line = line.strip(' ')
            gen_seq = gen_seq + line
    return gen_seq



# hash genome and add it to a hash table with kmers corresponding to the positions of the first nucleotide
# of the kmer
def split_genome(genome_hash_list, max_kmer_len):
    genome_hash_chunk_array = {}
    nuc_pos = 0
    while nuc_pos < len(genome_hash_list):
        window = genome_hash_list[int(nuc_pos):int(max_kmer_len)+int(nuc_pos)]
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



# verbose genome hash table gen
def verbose_split_genome(genome_hash_list, max_kmer_len):
    genome_hash_chunk_array = {}
    nuc_pos = 0
    while nuc_pos < len(genome_hash_list):
        window = genome_hash_list[int(nuc_pos):int(max_kmer_len)+int(nuc_pos)]
        #window = line[nuc_pos:max_kmer_len]
        window = window.upper()
        window_hash = hash(window)
        nuc_pos+=1
        #
        if not window_hash in genome_hash_chunk_array:
            # add kmer to the hash table and add
            genome_hash_chunk_array[window_hash] = [window]

        elif window_hash in genome_hash_chunk_array:
            genome_hash_chunk_array[window_hash].append(window)
    return genome_hash_chunk_array


# loop through both hash tables and identify shared sequences
# uses single kmers per read table
def hash_table_kmer_matcher(genome_kmer_hash_table, read_hash_table):

    kmer_match_locations = {}
    for read_seq, read_pos in read_hash_table.items():
        if read_seq in genome_kmer_hash_table:
            gen_pos = genome_kmer_hash_table[read_seq]
            kmer_match_locations[read_pos] = gen_pos
    return kmer_match_locations


# loop through both hash tables and identify shared sequences
# uses multiple kmers per read
def hash_table_multi_kmer_matcher(genome_kmer_hash_table, read_hash_table):
    total_matchs_table = {}
    for read_number, kmer_dict in read_hash_table.items():
        kmer_match_locations = {}
        for read_seq, read_pos in kmer_dict.items():
            if not read_pos in kmer_match_locations:
                if read_seq in genome_kmer_hash_table:
                    gen_pos = genome_kmer_hash_table[read_seq]
                    kmer_match_locations[read_pos] = [gen_pos]
            elif read_pos in kmer_match_locations:
                if read_seq in genome_kmer_hash_table:
                    gen_pos = genome_kmer_hash_table[read_seq]
                    kmer_match_locations[read_pos].append(gen_pos)
        if len(kmer_match_locations) != 0:

            total_matchs_table[read_number] = kmer_match_locations

    return total_matchs_table


# splits genome into kmers in reverse order and adds them to the hash table
def split_genome_reverse(genome_hash_list, max_kmer_len):
    reverse_genome_hash_chunk_array = {}
    nuc_pos = 0
    while nuc_pos < len(genome_hash_list):
        #window = genome_hash_list[int(nuc_pos):int(max_kmer_len)+int(nuc_pos)]
        window = genome_hash_list[(len(genome_hash_list) - int(max_kmer_len)) - int(nuc_pos):len(genome_hash_list) - int(nuc_pos)]
        #window = line[nuc_pos:max_kmer_len]
        window = window.upper()
        window_hash = hash(window)
        nuc_pos+=1
        #
        if not window_hash in reverse_genome_hash_chunk_array:
            # add kmer to the hash table and add
            reverse_genome_hash_chunk_array[window_hash] = [nuc_pos]

        elif window_hash in reverse_genome_hash_chunk_array:
            reverse_genome_hash_chunk_array[window_hash].append(nuc_pos)
    return reverse_genome_hash_chunk_array






def main():
    args = parse_args()

    size = args.kmer_size
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
                    # qual_list.append(qual)
                    line_count = 0


#PRODUCES HASH TABLE OF KMERS FROM THE NEW READ FILE
    # read_kmer_hash_table = kmer_hash_gen(read_list, size)
    # print(read_kmer_hash_table)

    full_read_kmer_hash_table = full_read_kmer_hash_gen(read_list, size)
    # print(full_read_kmer_hash_table)
    # itemlist = list(full_read_kmer_hash_table.items())
    # for item in itemlist:
    #     print(item)

    del read_list

#STRIPS NEW LINES FROM GENOME FILE
    genome = genome_reader(args.genome_file)
    # print(genome)


#PRODUCES HASH TABLE OF KMERS FROM THE GENOME OR LOCI
    genome_kmer_hash_table = split_genome(genome, size)
    # print(genome_kmer_hash_table)




    # genome_kmer_hash_table = verbose_split_genome(genome, size)
    # itemlist = list(genome_kmer_hash_table.items())
    # itemlist = sorted(itemlist, key=lambda x: -len(x[1]))
    #     #print(list(map(lambda x: (x[0], x[1][0], len(x[1])), itemlist)))
    # lenlist = list(map(lambda x: (x[0], len(x[1])), itemlist))
    # print("\n".join(list(map(str, lenlist))))


#MATCH READ KMERS TO GENOME KMERS AND COLLECT LOCATION INFO IN NEW TABLE
    # kmer_matching = hash_table_kmer_matcher(genome_kmer_hash_table, read_kmer_hash_table)
    # print(kmer_matching)
    #
    multi_kmer_match = hash_table_multi_kmer_matcher(genome_kmer_hash_table, full_read_kmer_hash_table)
    print(multi_kmer_match)


    # itemlist = list(multi_kmer_match.items())
    # for item in itemlist:
    #     for key, value in item[1].items():
    #         for value_1 in value:
    #             print(value_1)

#PERFORM GENOME KMER SPLIT FOR REVERSE COMPLIMENT KMERS
    # reverse_genome_kmer_hash_table = split_genome_reverse(genome, size)
    # print(reverse_genome_kmer_hash_table)



if __name__ == '__main__':
    main()
