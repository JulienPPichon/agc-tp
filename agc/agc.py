#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""OTU clustering"""

import argparse
import sys
import os
import gzip
import statistics
from collections import Counter
# https://github.com/briney/nwalign3
# ftp://ftp.ncbi.nih.gov/blast/matrices/
import nwalign3 as nw

__author__ = "Your Name"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["Your Name"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Your Name"
__email__ = "your@email.fr"
__status__ = "Developpement"


def isfile(path):
    """Check if path is an existing file.
      :Parameters:
          path: Path to the file
    """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def get_arguments():
    """Retrieves the arguments of the program.
      Returns: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', '-amplicon_file', dest='amplicon_file', type=isfile, required=True, 
                        help="Amplicon is a compressed fasta file (.fasta.gz)")
    parser.add_argument('-s', '-minseqlen', dest='minseqlen', type=int, default = 400,
                        help="Minimum sequence length for dereplication")
    parser.add_argument('-m', '-mincount', dest='mincount', type=int, default = 10,
                        help="Minimum count for dereplication")
    parser.add_argument('-c', '-chunk_size', dest='chunk_size', type=int, default = 100,
                        help="Chunk size for dereplication")
    parser.add_argument('-k', '-kmer_size', dest='kmer_size', type=int, default = 8,
                        help="kmer size for dereplication")
    parser.add_argument('-o', '-output_file', dest='output_file', type=str,
                        default="OTU.fasta", help="Output file")
    return parser.parse_args()
    
    
def read_fasta(amplicon_file, minseqlen):
    """Take a fasta file and a minimum sequence length as argument.
        Return a generator of sequence with the length filtered.
    """"
    with gzip.open(amplicon_file, "r") as filin:
        for line in filin:
            if len(line) >= minseqlen:
                yield line.strip()


def dereplication_fulllength(amplicon_file, minseqlen, mincount):
    """Take a generator of sequence, a minimum sequence length and a minimum number of
        occurence as arguments.
        Return a generator of sequence with at least 40 occurences
    """
    occurence_dict = {}
    for seq in read_fasta(amplicon_file, minseqlen):
        if seq not in occurence_dict:
            occurence_dict[seq] = 1
        else:
            occurence_dict[seq] += 1
    for seq, nb_seq in sorted(occurence_dict.items(), key = lambda item: item[1], reverse = True):
        if nb_seq < 40:
        	break
        else:
        	yield [seq, nb_seq]
        	
        	
def get_chunks(sequence, chunk_size):
    """Take a sequence and a size as arguments.
        Cut the sequence in the correct number of chunks if there is at least
        four chunks returned.
    """
    chunks = []
    if len(sequence) < 4 * chunk_size:
        pass
    else:
        if len(sequence) % chunk_size == 0:
            nb_chunks = len(sequence) / chunk_size
            for chunk in range(nb_chunks - 1):
                chunks.append[sequence[chunk * chunk_size:(chunk + 1) * chunk_size]]
        else:
            nb_chunks = len(sequence) / chunk_size
            for chunk in range(nb_chunks - 2):
                chunks.append[sequence[chunk * chunk_size:(chunk + 1) * chunk_size]]
    return chunks


def cut_kmer(sequence, kmer_size):
    """Take a sequence and a siez of kmer as arguments.
        Return all the possible kmer as generator.
    """
    for i in range(len(sequence) - kmer_size + 1):
        yield sequence[i:i+kmer_size]
        
        
def get_unique_kmer(kmer_dict, sequence, id_seq, kmer_size):
    """Take a dictionnary with a kmer as index and a list of id as values.
    """
    for kmer in cut_kmer(sequence, kmer_size):
        if kmer not in kmer_dict:
            kmer_dict[kmer] = []
        kmer_dict[kmer].append(id_seq)
    return kmer_dict


def search_mates(kmer_dict, sequence, kmer_size):
    return [i[0] for i in Counter([ids for kmer in cut_kmer(sequence, km_size) 
        if kmer in kmer_dict for ids in kmer_dict[kmer]]).most_common(8)]


#==============================================================
# Main program
#==============================================================
def main():
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()
    get_chunks(args.amplicon_file, args.minseqlen, args.mincount, args.chunk_size)
    


if __name__ == '__main__':
    main()
