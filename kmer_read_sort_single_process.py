#!/usr/bin/env python

import sys
import argparse
import statistics


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--read_file_1')
    parser.add_argument('--genome_file')
    parser.add_argument('--kmer_size')
    parser.add_argument('--read_file_2')
    parser.add_argument('--read_size')
    return parser.parse_args()



def read_file_reader(read_file):
    # size = args.kmer_size
    # GET READ LEN TO EXTRAPOLATE KMER LENGTH
    line_count = 0
    header = ''
    seq = ''
    seperate = ''
    qual = ''
    read_list = []
    qual_list = []
    with open(read_file) as read1:
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
    return read_list


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
            if '>' in line:
                continue
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
                    if not len(gen_pos) > 5:
                        kmer_match_locations[read_pos] = list(gen_pos)
                    # print(kmer_match_locations[read_pos])
                    # kmer_match_locations[read_pos] = []
                    # kmer_match_locations.setdefault(read_pos, []).append(gen_pos)
            elif read_pos in kmer_match_locations:
                if read_seq in genome_kmer_hash_table:
                    gen_pos = genome_kmer_hash_table[read_seq]
                    if not len(gen_pos) > 5:
                        kmer_match_locations[read_pos].append(gen_pos)
                    # kmer_match_locations.setdefault(read_pos, []).append(gen_pos)
        if len(kmer_match_locations) != 0:

            total_matchs_table[read_number] = kmer_match_locations

    return total_matchs_table


# reverse complement genome
# def reverse_complement_genome(genome_seq):


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
        # window_hash = window
        nuc_pos+=1
        #
        if not window_hash in reverse_genome_hash_chunk_array:
            # add kmer to the hash table and add
            reverse_genome_hash_chunk_array[window_hash] = [nuc_pos]

        elif window_hash in reverse_genome_hash_chunk_array:
            reverse_genome_hash_chunk_array[window_hash].append(nuc_pos)
    reverse_genome_hash_chunk_array.pop('', None)
    return reverse_genome_hash_chunk_array


# get location information from the output of kmer matching and use this to map full reads to the location in the genome
def full_read_location_grabber(kmer_match_hash_table, read_len, kmer_len):
    read_map_locations = {}
    for read_num, dict in kmer_match_hash_table.items():
        if len(dict) >= 3:
            read_map_locations[read_num] = []
            locs_dict = {}
            locs_list = []

            # this is just here becuase i had to add a nested list to the dictionary in the kmer_match_hash_table step
            for kmer_num, match_loc in dict.items():
                # print(len(match_loc[0]))
                if len(match_loc) <= 1:

                    locs_dict[kmer_num] = match_loc
                    locs_list.append(match_loc[0])
                else:
                    locs_dict[kmer_num] = match_loc
                    locs_list.append(match_loc[0])
                # else:
                #     for item in match_loc[0]:
                #         # locs_dict[kmer_num].append(item)
                #         locs_list.append(item)

                #DON'T KNOW HOW TO LOGICALLY HANDLE KMERS THAT MAP TO MULTIPLE POSITIONS

            kmer_range = max(locs_list) - min(locs_list)

            # kmer_range = int(kmer_range)
            read_len = int(read_len)
            if kmer_range <= read_len:
                if min(locs_dict.keys()) == 1:
                    first_kmer = min(locs_dict.keys())
                    first_kmer_genome_match = locs_dict[first_kmer][0]
                    read_map_locations[read_num] = first_kmer_genome_match
                elif min(locs_dict.keys()) != 1:
                    first_kmer = min(locs_dict.keys())
                    missing_nucs = int(first_kmer) * int(kmer_len)
                    read_start = locs_dict[first_kmer][0] - missing_nucs
                    read_map_locations[read_num] = read_start


            elif kmer_range > read_len:

                while kmer_range > read_len and len(locs_list) >= 3:
                    locs_list.remove(max(locs_list))
                    locs_list.remove(min(locs_list))
                    # print(locs_list)
                    kmer_range = max(locs_list) - min(locs_list)

            if min(locs_dict.keys()) == 1:
                first_kmer = min(locs_dict.keys())
                first_kmer_genome_match = locs_dict[first_kmer][0]
                read_map_locations[read_num] = first_kmer_genome_match
            elif min(locs_dict.keys()) != 1:
                first_kmer = min(locs_dict.keys())
                missing_nucs = int(first_kmer) * int(kmer_len)
                read_start = locs_dict[first_kmer][0] - missing_nucs
                read_map_locations[read_num] = read_start
            # first_kmer = min(locs_dict.keys())
            # first_kmer_genome_match = locs_dict[first_kmer][0]
            # read_map_locations[read_num] = first_kmer_genome_match

    return read_map_locations

# def forward_map_reads(read_locations_table, read_list, genome_file):
#     for read_num, read_loc in read_locations_table.items():
#         read_counter = 0
#
#         for read in read_list:
#             read_length = len(read)
#
#             read_counter+=1
#             if read_num == read_counter:
#                 genome_seq = genome_file[int(read_loc) - 1 :(int(read_loc) + int(read_length)) - 1]
#                 genome_seq = genome_seq.upper()
#                 print(read)
#                 print(genome_seq)
#                 print('new_seqs')
#                 # for nuc1 in read:
#                 #     for nuc2 in genome_seq:
#                 #         if nuc1 == nuc2:


def forward_map_reads(read_locations_table, read_file, genome_file):
    # print(read_locations_table)
    with open(read_file) as rf:
        read_counter = 0
        while True:
            read_counter+=1
            rf.readline()
            seq = rf.readline()
            # print(seq)
            rf.readline()
            qual = rf.readline()
            if read_counter in read_locations_table:
                # continue
                # print(read_locations_table.get(read_counter))
                read_loc = read_locations_table.get(read_counter)
                read_length = len(seq)
                genome_seq = genome_file[int(read_loc) - 1 :(int(read_loc) + int(read_length)) - 1]
                genome_seq = genome_seq.upper()
                # print(seq)
                # print(genome_seq)
                # print('\n')
                # print('\n')
                # print('new_seqs')
            if len(seq) == 0:
                break




    # for read_num, read_loc in read_locations_table.items():
    #     read_counter = 0
    #
    #     for read in read_list:
    #         read_length = len(read)
    #
    #         read_counter+=1
    #         if read_num == read_counter:
    #             genome_seq = genome_file[int(read_loc) - 1 :(int(read_loc) + int(read_length)) - 1]
    #             genome_seq = genome_seq.upper()
    #             print(read)
    #             print(genome_seq)
    #             print('new_seqs')
    #             # for nuc1 in read:
    #             #     for nuc2 in genome_seq:
    #             #         if nuc1 == nuc2:


def main():
    args = parse_args()




    size = args.kmer_size
    input_file = open(args.read_file_1, 'rb')

    #STRIPS NEW LINES FROM GENOME FILE
    genome = genome_reader(args.genome_file)

    genome_kmer_hash_table = split_genome(genome, size)
    # print(genome_kmer_hash_table)

    while(True):
        read_read_file = []
        quals = []
        line_counter = 0
        read_counter = 0
        for lines in range(1000000):
            line_counter+=1
            if line_counter == 1:
                input_file.readline()
            if line_counter == 2:

                seq = input_file.readline()
                # print(seq)
                read_counter+=1
            if line_counter == 3:

                input_file.readline()
            if line_counter == 4:

                input_file.readline()
                line_counter = 0
                read_read_file.append(seq)

            #PRODUCES HASH TABLE OF KMERS FROM THE NEW READ FILE
            # read_kmer_hash_table = kmer_hash_gen(read_list, size)
            # print(read_kmer_hash_table)

            # read_read_file = read_file_reader(args.read_file_1)
            # # print(read_read_file)
            # print(sys.getsizeof(read_read_file))

            full_read_kmer_hash_table = full_read_kmer_hash_gen(read_read_file, size)
            # print(full_read_kmer_hash_table)
            # itemlist = list(full_read_kmer_hash_table.items())
            # for item in itemlist:
            #     print(item)
            print(sys.getsizeof(full_read_kmer_hash_gen))


            # #STRIPS NEW LINES FROM GENOME FILE
            #     genome = genome_reader(args.genome_file)
            # print(genome)
            # print(sys.getsizeof(genome))

            #PRODUCES HASH TABLE OF KMERS FROM THE GENOME OR LOCI
            genome_kmer_hash_table = split_genome(genome, size)
            # print(genome_kmer_hash_table)
            print(sys.getsizeof(genome_kmer_hash_table))



            # genome_kmer_hash_table = verbose_split_genome(genome, size)
            # itemlist = list(genome_kmer_hash_table.items())
            # itemlist = sorted(itemlist, key=lambda x: -len(x[1]))
            #     #print(list(map(lambda x: (x[0], x[1][0], len(x[1])), itemlist)))
            # lenlist = list(map(lambda x: (x[0], len(x[1])), itemlist))
            # print("\n".join(list(map(str, lenlist))))


            #MATCH READ KMERS TO GENOME KMERS AND COLLECT LOCATION INFO IN NEW TABLE
            # kmer_matching = hash_table_kmer_matcher(genome_kmer_hash_table, read_kmer_hash_table)
            # print(kmer_matching)


            # del read_read_file
            multi_kmer_match = hash_table_multi_kmer_matcher(genome_kmer_hash_table, full_read_kmer_hash_table)
            # print(multi_kmer_match)
            print(sys.getsizeof(multi_kmer_match))


            # itemlist = list(multi_kmer_match.items())
            # for item in itemlist:
            #     for key, value in item[1].items():
            #         for value_1 in value:
            #             print(value_1)

            #PERFORM GENOME KMER SPLIT FOR REVERSE COMPLIMENT KMERS
        #     reverse_genome_kmer_hash_table = split_genome_reverse(genome, size)
        #     # print(reverse_genome_kmer_hash_table)
        #
            #READ IN SECOND SET OF READS
        #     read_read_file_2 = read_file_reader(args.read_file_2)
        #     # print(read_read_file_2)
        #
            #HASH READS FROM FILE TWO AND ADD THEM TO HASH TABLE
        #     second_full_read_kmer_hash_table = full_read_kmer_hash_gen(read_read_file_2, size)
        #     # print(second_full_read_kmer_hash_table)
        #
            #RUN HASH MATCH ON KMERS FROM SECOND SET OF READS AGAINST REVERSE KMER HASHED GENOME
        #     second_multi_kmer_match = hash_table_multi_kmer_matcher(reverse_genome_kmer_hash_table, second_full_read_kmer_hash_table)
        #     # print(second_multi_kmer_match)


            del full_read_kmer_hash_table
            del genome_kmer_hash_table



        #GRAB MAPPING LOCATIONS AND CHECK IF MULTIPLE KMERS FROM A SINGLE READ FALL INTO THE SAME AREA OF THE GENOME
            match_loc = full_read_location_grabber(multi_kmer_match, args.read_size, size)
            # print(match_loc)
            print(sys.getsizeof(match_loc))
            del multi_kmer_match

        #TESTING READ MAPPING TO THE GENOME
            read_mapping = forward_map_reads(match_loc, args.read_file_1, genome)
            # print(read_mapping)
            print(sys.getsizeof(read_mapping))

        print('#########################################################################################################')
        if len(seq) == 0:
            break



if __name__ == '__main__':
    main()
