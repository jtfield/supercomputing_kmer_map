#! /usr/bin/python

import numpy
import sys
import argparse

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--read_seq')
    return parser.parse_args()


def main():
    args = parse_args()

    read_seq = args.read_seq

    seq_length = len(read_seq)
    #print(seq_length)
    kmer_len = seq_length / 10
    print(kmer_len)
    print("waffle")


if __name__ == '__main__':
    main()
