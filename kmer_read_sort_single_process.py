#!/usr/bin/env python


import argparse


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--read_file_1')
    parser.add_argument('--genome_file')
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

# get kmer sequences in list form
def kmer_hash_gen(read, k_size):
#################################################################
# code for getting kmer from the beginning of read
    kmer_list = []
    for sing_read, size in zip(read, k_size):
        seq_len = len(sing_read)
        # half_read = seq_len / 2
        # kmer_start_pos = half_read - (size / 2)
        # kmer_end_pos = half_read + (size / 2)
        kmer = sing_read[:int(size)]
        kmer_list.append(kmer)
    return(kmer_list)

##############################################################
# code for getting kmer from the middle of read
    # kmer_list = []
    # for sing_read, size in zip(read, k_size):
    #     seq_len = len(sing_read)
    #     half_read = seq_len / 2
    #     kmer_start_pos = half_read - (size / 2)
    #     kmer_end_pos = half_read + (size / 2)
    #     kmer = sing_read[kmer_start_pos:kmer_end_pos]
    #     kmer_list.append(kmer)
    # return(kmer_list)

# convert kmers into hash forms
def kmer_convert_to_bit(full_kmer_list):
    # A = 00
    # C = 01
    # G = 10
    # T = 11
    hash_kmers = []
    for seq in full_kmer_list:
        #print("converting kmers to hash")
        kmer_seq = []
        for letter in seq:
            if letter == 'A':
                kmer_seq.append('00')
            elif letter == 'C':
                kmer_seq.append('01')
            elif letter == 'G':
                kmer_seq.append('10')
            elif letter == 'T':
                kmer_seq.append('11')
        hash_kmers.append(kmer_seq)
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


# begin mapping the read kmers to the genome kmers
def kmer_matcher(read_kmer_list, genome_kmer_list):
    # initialize list to collect all information about kmer mappings
    # initialize a counter to keep track of which kmer is which
    read_kmer_counter = 0
    master_kmer_stats_dict = []

    # loop through the read kmers and prepare to map each read kmer to each genome kmer
    for read_kmer in read_kmer_list:
        read_kmer_counter+=1

        # initialize a list to hold the stats for each read kmer
        # this will be added to the master kmer info list
        read_kmer_stats = {}
        read_kmer_stats['read_number'] = read_kmer_counter

        # get the length of the kmer to identify mismatch threshold
        mismatch_threshold = 0
        kmer_len = len(read_kmer)
        four_percent_of_kmer = 0.04 * kmer_len
        if four_percent_of_kmer < 1:
            mismatch_threshold = 1
        elif four_percent_of_kmer >= 1:
            mismatch_threshold = round(four_percent_of_kmer)
        read_kmer_stats['mismatch_threshold'] = mismatch_threshold
        read_kmer_stats['kmer_length'] = kmer_len
        # read_kmer_stats['read_kmer_seq'] = read_kmer

        # THIS IS UNNECESSARY TO THE CODE. IT JUST HALTS MAPPING AFTER ONE READ KMER
        # FOR TESTING PURPOSES ONLY
        ##########################################################################
        # if read_kmer_counter > 1:
        #     break
        # else:
        ##########################################################################


            # keep track of the position of the genome kmer being mapped to
        genome_kmer_count = 0
        genome_kmer_position_tracker = {}

        best_match_loc = {}
        best_match_score = {}
        best_match_counter = 0
        best_match_position = 0

        for genome_kmer in genome_kmer_list:
            genome_kmer_count+=1


            mismatch_counter = 0
            correct_match_counter = 0
            for read_kmer_letter, genome_kmer_letter in zip(read_kmer, genome_kmer):

                if read_kmer_letter == genome_kmer_letter:
                    correct_match_counter+=1
                elif read_kmer_letter == genome_kmer_letter:
                    mismatch_counter+=1

            # if mismatch_counter >= mismatch_threshold:
            #     print("skip")
            #
            # elif mismatch_counter > mismatch_threshold:
            genome_kmer_position_tracker[genome_kmer_count] = correct_match_counter


            if correct_match_counter > best_match_counter:
                best_match_counter = correct_match_counter
                best_match_position = genome_kmer_count

            elif correct_match_counter <= best_match_counter:
                best_match_counter = best_match_counter
                best_match_position = best_match_position


            #read_kmer_stats['genome_kmer_mapping_stats'] = genome_kmer_position_tracker
        # best_match_loc[best_match_position] = best_match_counter
        # best_match_loc['best_match_loc'] = best_match_position
        # best_match_score['best_match_score'] = best_match_counter
        read_kmer_stats['best_match_score'] = best_match_counter
        read_kmer_stats['best_match_loc'] = best_match_position



        master_kmer_stats_dict.append(read_kmer_stats)
    return master_kmer_stats_dict
    # for key in master_kmer_stats_dict:
    #     find_max = key['genome_kmer_mapping_stats']
    #
    #     inverse = [(value, key) for key, value in find_max.items()]
    #     print max(inverse)[1]



    #return master_kmer_stats_dict
        # print("read_kmer_number")
        # print(read_kmer_counter)
        # print("read matches to loci")
        # print(correct_match_counter)
        # print("read mismatches loci")
        # print(mismatch_counter)

# match kmers to reads and begin mapping reads to the locations associated with the kmers
def read_kmer_matcher(kmer_list_of_dicts, read_list, byte_genome):
    # print(byte_genome)
    for kmer_dict in kmer_list_of_dicts:
        read_count = 0
        kmer_id = kmer_dict['read_number']
        pos_id = kmer_dict['best_match_loc']

        genome_pos = pos_id - 1

        read_count = 0
        for read in read_list:
            read_len = len(read)
            read_count+=1
            segment_end = (pos_id + read_len) - 1
            combined_list = []
            if kmer_id == read_count:
                print(genome_pos)
                print(segment_end)
                print(byte_genome[genome_pos:segment_end])
                print(read)
                print("start over")
                # genome_segment = byte_genome[genome_pos:read_len]
                # print(read)
                # print(genome_segment)
                # print("start over")



        # for read_seq in read_list:
        #     read_len = len(read_seq)
        #     read_count+=1
        #     kmer_loc = 0
        #     if kmer_id == read_count:
        #         kmer_loc = kmer_dict['best_match_loc']
        #
        #         genome_letter_position = 0
        #         read_len_counter = 0
        #         mapping_comparison_segment = []
        #         genome_segment = byte_genome[int(kmer_loc):int(read_len)]
        #         # mapping_comparison_segment.append(read_seq)
        #         # genome_segment = []
        #         # for letter in byte_genome:
        #         #     genome_letter_position+=1
        #         #     if genome_letter_position >= kmer_loc:
        #         #         read_len_counter+=1
        #         #         if read_len_counter <= read_len:
        #                     # genome_segment.append(letter)
        #
        #         # print(read_seq)
        #         # print(genome_segment)
        #     mapping_comparison_segment.append(read_seq)
        #     mapping_comparison_segment.append(genome_segment)
            #print(mapping_comparison_segment)




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
    kmer_sizes1 = kmer_len_gen(read_list)
    # print(kmer_sizes1)
    max_kmer_size1 = max(kmer_sizes1)
    # print(max_kmer_size1)

    # PRODUCES SEQUENCES FOR KMERS OF APPROPRIATE LENGTHS
    kmer_seqs1 = kmer_hash_gen(read_list , kmer_sizes1)
    # print(kmer_seqs1)
    # counters = 0

    # for line in read_list:
    #     print(line)
    #
    # with open(args.genome_file) as sequence:
    #     for line in sequence:
    #         print(line.upper())

    # CONVERTS KMERS INTO HASH VALUES (A = 00, C = 01, G = 10, T = 11)
    kmer_convert1 = kmer_convert_to_bit(kmer_seqs1)
    # counters = 0
    # for i in kmer_convert:
    #     counters+=1
    # print(counters)

    # CONVERT GENOME OR MSA TO HASH
    gen_hasher1 = gen_convert_to_bit(args.genome_file)
    # print(gen_hasher1)

    # CONVERT HASH GENOME TO LIST OF KMERS FOR MAPPING
    hash_gen_chunker1 = split_genome(gen_hasher1, max_kmer_size1)
    # hash_gen_counter = 0
    # for i in hash_gen_chunker1:
    #     hash_gen_counter+=1
    #     print("genome kmer")
    #     print(hash_gen_counter)
    #     print(i)

    # MAP ALL READ KMERS TO ALL GENOME KMERS AND CONSTRUCT DICT
    # DICT CONTAINS INFO ON BEST MAP LOCATION WITHIN THE GENOME KMERS
    kmer_search1 = kmer_matcher(kmer_convert1, hash_gen_chunker1)
    # print(kmer_search1)

    # CONVERT READS TO BYTES
    read_convert = kmer_convert_to_bit(read_list)
    # print(read_convert)

    #
    kmer_and_read_matcher = read_kmer_matcher(kmer_search1, read_convert, gen_hasher1)
    #print(kmer_and_read_matcher)

if __name__ == '__main__':
    main()
