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
    kmer_list = []
    for sing_read, size in zip(read, k_size):
        seq_len = len(sing_read)
        half_read = seq_len / 2
        kmer_start_pos = half_read - (size / 2)
        kmer_end_pos = half_read + (size / 2)
        kmer = sing_read[kmer_start_pos:kmer_end_pos]
        kmer_list.append(kmer)
    return(kmer_list)

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

def kmer_matcher(read_kmer_list, genome_kmer_list):
    read_kmer_counter = 0
    for read_kmer in read_kmer_list:
        kmer_len = len(read_kmer)
        four_percent_of_kmer = 0.04 * kmer_len
        print(four_percent_of_kmer)
        read_kmer_counter+=1
        if read_kmer_counter > 1:
            break
        else:
            genome_kmer_count = 0
            for genome_kmer in genome_kmer_list:
                genome_kmer_count+=1
                # print("new_genome_kmer")
                # print("read_kmer_len {}".format(kmer_len))
                mismatch_counter = 0
                correct_match_counter = 0
                for read_kmer_letter, genome_kmer_letter in zip(read_kmer, genome_kmer):
                    print("set")
                    #print("read_kmer_len {}".format(kmer_len))
                    # print(read_kmer_letter)
                    # print(genome_kmer_letter)
        #         if read_kmer_letter == genome_kmer_letter:
        #             correct_match_counter+=1
        #         elif read_kmer_letter == genome_kmer_letter:
        #             mismatch_counter+=1
        # print("read_kmer_number")
        # print(read_kmer_counter)
        # print("read matches to loci")
        # print(correct_match_counter)
        # print("read mismatches loci")
        # print(mismatch_counter)



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
    # print(kmer_sizes)
    max_kmer_size = max(kmer_sizes)
    # print(max_kmer_size)

    # PRODUCES SEQUENCES FOR KMERS OF APPROPRIATE LENGTHS
    kmer_seqs = kmer_hash_gen(read_list , kmer_sizes)
    # print(kmer_seqs)

    # CONVERTS KMERS INTO HASH VALUES (A = 00, C = 01, G = 10, T = 11)
    kmer_convert = kmer_convert_to_bit(kmer_seqs)
    # print(kmer_convert)

    # CONVERT GENOME OR MSA TO HASH
    gen_hasher = gen_convert_to_bit(args.genome_file)
    # print(gen_hasher)

    # CONVERT HASH GENOME TO LIST OF KMERS FOR MAPPING
    hash_gen_chunker = split_genome(gen_hasher, max_kmer_size)
    # print(hash_gen_chunker)

    kmer_search = kmer_matcher(kmer_convert, hash_gen_chunker)
    print(kmer_search)

if __name__ == '__main__':
    main()
