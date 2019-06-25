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


# loop through both hash tables and identify shared sequences
def hash_table_kmer_matcher(genome_kmer_hash_table, read_hash_table):

    kmer_match_locations = {}
    for read_seq, read_pos in read_hash_table.items():
        if read_seq in genome_kmer_hash_table:
            gen_pos = genome_kmer_hash_table[read_seq]
            kmer_match_locations[read_pos] = gen_pos
    return kmer_match_locations


# get the length of the reads, read in the reads to a new hash table
# def full_read_hash_table(read_list):
#     read_table = {}
#     read_count = 0
#     for sing_read in read_list:
#         read_count+=1
#         sing_read = sing_read.upper()
#         read_hash = hash(sing_read)
#         read_table[read_hash] = read_count
#     return(read_table)



# # begin mapping the read kmers to the genome kmers
# def kmer_matcher(read_kmer_list, genome_kmer_list):
#     # initialize list to collect all information about kmer mappings
#     # initialize a counter to keep track of which kmer is which
#     read_kmer_counter = 0
#     master_kmer_stats_dict = []
#
#     # loop through the read kmers and prepare to map each read kmer to each genome kmer
#     for read_kmer in read_kmer_list:
#         read_kmer_counter+=1
#
#         # initialize a list to hold the stats for each read kmer
#         # this will be added to the master kmer info list
#         read_kmer_stats = {}
#         read_kmer_stats['read_number'] = read_kmer_counter
#
#         # get the length of the kmer to identify mismatch threshold
#         mismatch_threshold = 0
#         kmer_len = len(read_kmer)
#         four_percent_of_kmer = 0.04 * kmer_len
#         if four_percent_of_kmer < 1:
#             mismatch_threshold = 1
#         elif four_percent_of_kmer >= 1:
#             mismatch_threshold = round(four_percent_of_kmer)
#         read_kmer_stats['mismatch_threshold'] = mismatch_threshold
#         read_kmer_stats['kmer_length'] = kmer_len
#         # read_kmer_stats['read_kmer_seq'] = read_kmer
#
#         # THIS IS UNNECESSARY TO THE CODE. IT JUST HALTS MAPPING AFTER ONE READ KMER
#         # FOR TESTING PURPOSES ONLY
#         ##########################################################################
#         # if read_kmer_counter > 1:
#         #     break
#         # else:
#         ##########################################################################
#
#
#             # keep track of the position of the genome kmer being mapped to
#         genome_kmer_count = 0
#         genome_kmer_position_tracker = {}
#
#         best_match_loc = {}
#         best_match_score = {}
#         best_match_counter = 0
#         best_match_position = 0
#
#         for genome_kmer in genome_kmer_list:
#             genome_kmer_count+=1
#
#
#             mismatch_counter = 0
#             correct_match_counter = 0
#             for read_kmer_letter, genome_kmer_letter in zip(read_kmer, genome_kmer):
#
#                 if read_kmer_letter == genome_kmer_letter:
#                     correct_match_counter+=1
#                 elif read_kmer_letter == genome_kmer_letter:
#                     mismatch_counter+=1
#
#             # if mismatch_counter >= mismatch_threshold:
#             #     print("skip")
#             #
#             # elif mismatch_counter > mismatch_threshold:
#             genome_kmer_position_tracker[genome_kmer_count] = correct_match_counter
#
#
#             if correct_match_counter > best_match_counter:
#                 best_match_counter = correct_match_counter
#                 best_match_position = genome_kmer_count
#
#             elif correct_match_counter <= best_match_counter:
#                 best_match_counter = best_match_counter
#                 best_match_position = best_match_position
#
#
#             #read_kmer_stats['genome_kmer_mapping_stats'] = genome_kmer_position_tracker
#         # best_match_loc[best_match_position] = best_match_counter
#         # best_match_loc['best_match_loc'] = best_match_position
#         # best_match_score['best_match_score'] = best_match_counter
#         read_kmer_stats['best_match_score'] = best_match_counter
#         read_kmer_stats['best_match_loc'] = best_match_position
#
#
#
#         master_kmer_stats_dict.append(read_kmer_stats)
#     return master_kmer_stats_dict
#     # for key in master_kmer_stats_dict:
#     #     find_max = key['genome_kmer_mapping_stats']
#     #
#     #     inverse = [(value, key) for key, value in find_max.items()]
#     #     print max(inverse)[1]
#
#
#
#     #return master_kmer_stats_dict
#         # print("read_kmer_number")
#         # print(read_kmer_counter)
#         # print("read matches to loci")
#         # print(correct_match_counter)
#         # print("read mismatches loci")
#         # print(mismatch_counter)
#
# # match kmers to reads and begin mapping reads to the locations associated with the kmers
# def read_kmer_matcher(kmer_list_of_dicts, read_list, byte_genome):
#     # print(byte_genome)
#     for kmer_dict in kmer_list_of_dicts:
#         read_count = 0
#         kmer_id = kmer_dict['read_number']
#         pos_id = kmer_dict['best_match_loc']
#
#         genome_pos = pos_id - 1
#
#         read_count = 0
#         combined_list = []
#         for read in read_list:
#             read_len = len(read)
#             read_count+=1
#             segment_end = (pos_id + read_len) - 1
#             # combined_list = []
#             if kmer_id == read_count:
#                 # print(genome_pos)
#                 # print(segment_end)
#                 genome_segment = byte_genome[genome_pos:segment_end]
#                 combined_list.append(genome_segment)
#                 combined_list.append(read)
#                 # print(read)
#         print(combined_list)
#                 # print("start over")
#                 # genome_segment = byte_genome[genome_pos:read_len]
#                 # print(read)
#                 # print(genome_segment)
#                 # print("start over")
#
#
#
#         # for read_seq in read_list:
#         #     read_len = len(read_seq)
#         #     read_count+=1
#         #     kmer_loc = 0
#         #     if kmer_id == read_count:
#         #         kmer_loc = kmer_dict['best_match_loc']
#         #
#         #         genome_letter_position = 0
#         #         read_len_counter = 0
#         #         mapping_comparison_segment = []
#         #         genome_segment = byte_genome[int(kmer_loc):int(read_len)]
#         #         # mapping_comparison_segment.append(read_seq)
#         #         # genome_segment = []
#         #         # for letter in byte_genome:
#         #         #     genome_letter_position+=1
#         #         #     if genome_letter_position >= kmer_loc:
#         #         #         read_len_counter+=1
#         #         #         if read_len_counter <= read_len:
#         #                     # genome_segment.append(letter)
#         #
#         #         # print(read_seq)
#         #         # print(genome_segment)
#         #     mapping_comparison_segment.append(read_seq)
#         #     mapping_comparison_segment.append(genome_segment)
#             #print(mapping_comparison_segment)




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
                    qual_list.append(qual)
                    line_count = 0


#PRODUCES HASH TABLE OF KMERS FROM THE NEW READ FILE
    read_kmer_hash_table = kmer_hash_gen(read_list, size)
    #print(read_kmer_hash_table)

#STRIPS NEW LINES FROM GENOME FILE
    genome = genome_reader(args.genome_file)
    # print(genome)


#PRODUCES HASH TABLE OF KMERS FROM THE GENOME OR LOCI
    genome_kmer_hash_table = split_genome(genome, size)
    #print(genome_kmer_hash_table)

#MATCH READ KMERS TO GENOME KMERS AND COLLECT LOCATION INFO IN NEW TABLE
    kmer_matching = hash_table_kmer_matcher(genome_kmer_hash_table, read_kmer_hash_table)
    print(kmer_matching)


    del read_kmer_hash_table
    del genome_kmer_hash_table

#ADD FULL READS TO A HASH TABLE AFTER HASHING
    #read_hashing = full_read_hash_table(read_list)
    # print(read_hashing)

#ADD READ SIZED GENOME WINDOW CHUNKS TO HASH TABLE
    # genome_hashing =


    # # PRODUCES SEQUENCES FOR KMERS OF APPROPRIATE LENGTHS
    # kmer_seqs1 = kmer_hash_gen(read_list , kmer_sizes1)
    # # print(kmer_seqs1)
    # # counters = 0
    #
    # # for line in read_list:
    # #     print(line)
    # #
    # # with open(args.genome_file) as sequence:
    # #     for line in sequence:
    # #         print(line.upper())
    #
    # # CONVERTS KMERS INTO HASH VALUES (A = 00, C = 01, G = 10, T = 11)
    # kmer_convert1 = kmer_convert_to_bit(kmer_seqs1)
    # # counters = 0
    # # for i in kmer_convert:
    # #     counters+=1
    # # print(counters)
    #
    # # CONVERT GENOME OR MSA TO HASH
    # gen_hasher1 = gen_convert_to_bit(args.genome_file)
    # # print(gen_hasher1)
    #
    # # CONVERT HASH GENOME TO LIST OF KMERS FOR MAPPING
    # hash_gen_chunker1 = split_genome(gen_hasher1, max_kmer_size1)
    # # hash_gen_counter = 0
    # # for i in hash_gen_chunker1:
    # #     hash_gen_counter+=1
    # #     print("genome kmer")
    # #     print(hash_gen_counter)
    # #     print(i)
    #
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
    # #print(kmer_and_read_matcher)

if __name__ == '__main__':
    main()
