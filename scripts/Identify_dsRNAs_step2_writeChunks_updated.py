###################################################################
#Script Name    : Identify_dsRNAs_step2_writeChunks_updated.py
#Description    : Go through 50 bp windows and get those with at least 3 editing sites
#Args_in        : All ES sites and a minimum number of people (This is five people)
#Args_out   :
#Author     : Mudra Choudhury
#Email      : mudrachoudhury3@gmail.com
#Date       :
#Last Edited    : 4/24/19
###################################################################

import os
import collections
import argparse
from collections import defaultdict
from itertools import islice
import sys

parser = argparse.ArgumentParser(description='Input Files')
parser.add_argument('-i', metavar='input_file', type=str, help='input file with all editing sites')
parser.add_argument('-n', metavar='min_es', type=str, help='minimum number of editing sites per window')
parser.add_argument('-s', metavar='window_size', type=str, help='size of sliding window for each editing site')
parser.add_argument('-c', metavar='chromosome', type=str, help='Chromosome that the script will be analyzing, ex: chr12')


args = parser.parse_args()
input_fn = args.i
min_es_in_window = args.n
window_size = args.s
chromosome_argument = args.c

split_input = input_fn.split("_")

output_fn = "".join(("../data/All_EES_",window_size,"_window_regions_w_min_",min_es_in_window,"_ES_per_window_using_REDI_2020_",chromosome_argument,"_updated.bed"))


def slidingWindow(sequence,winSize,step=1):
    """Returns a generator that will iterate through
    the defined chunks of input sequence.  Input sequence
    must be iterable."""
    # Verify the inputs
    try: it = iter(sequence)
    except TypeError:
        raise Exception("**ERROR** sequence must be iterable.")
    if not ((type(winSize) == type(0)) and (type(step) == type(0))):
        raise Exception("**ERROR** type(winSize) and type(step) must be int.")
    if step > winSize:
        raise Exception("**ERROR** step must not be larger than winSize.")
    if winSize > len(sequence):
        raise Exception("**ERROR** winSize must not be larger than sequence length.")
    # Pre-compute number of chunks to emit
    numOfChunks = ((len(sequence)-winSize)/step)+1
    # Do the work
    for i in range(0,numOfChunks*step,step):
        yield sequence[i:i+winSize]


def get_EES_windows(chromosome, positions, es_set, window_size, extra_info): #get editing enriched windows around an editing site
    ees_regions = []
    min_es = int(min_es_in_window)
    step = 1
    windows = slidingWindow(positions, window_size, step) #retrieve all windows of the window size within
    for window in windows:
        #print window
        intersection = es_set.intersection(window)
        #print len(intersection), min_es
        if len(intersection) >= min_es:
            #print len(intersection)
            start, end = (window[0]-1), window[-1]
            ees_line = [chromosome, str(start), str(end), "NA", "NA"] + extra_info
            ees_regions.append(ees_line)
    return ees_regions


chromosomes = [chromosome_argument]
#chromosomes = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", 'chrX', "chrY"]
#chromosomes = ["chr20", "chr21", "chr22", 'chrX', "chrY"]
#chromosomes = ["chr21","chr22"]

# Go through every editing site (50bp windows) and get number of editing sites in all the possible windows for each editing site
#for each chromosome, aggregate all editing sites
#All_EES_regions = []
start_print = "".join(("Starting script...\n"))
sys.stderr.write(start_print)

for chrom in chromosomes:
    All_EES_regions = []
    print chrom
    chrom_positions_es_set = set()
    n = 0
    input_file = open(input_fn, "r")
    for line in input_file:
        if n == 0:
            n += 1
            continue
        else:
            split_line = line.split("\t")
            chrom_positions_es_set.add(int(split_line[1]))
    #print chrom_positions_es_set
    #Go through the 50 bp windows and check if they have at least 3 editing sites
    i = 0
    input_file = open(input_fn, "r")
    for es_line in input_file:
        if i == 0:
            i += 1
            continue
        split_es_line = ((es_line).rstrip()).split("\t")
        #if split_es_line[0] == chrom:
        #print(chrom)
        position = int(split_es_line[1])
        extra_info = [split_es_line[2], split_es_line[3], split_es_line[4]]
        start_window, end_window = (position - int(window_size)), (position + int(window_size))
        window_seq = range(start_window,end_window)
        #print window_seq
        list_of_ees_regions = get_EES_windows(chrom, window_seq, chrom_positions_es_set, int(window_size), extra_info)
        #if list_of_ees_regions != []:
        #    print(list_of_ees_regions)
        All_EES_regions = All_EES_regions + list_of_ees_regions
    #if chrom == chromosomes[0]:
        if i % 1000 == 0:
            line_no = str(i)
            log_print = "".join(("processed line: ",line_no,"\n"))
            sys.stderr.write(log_print)
            if i == 1000:
                with open(output_fn, 'w') as file:
                    file.writelines('\t'.join(j) + '\n' for j in All_EES_regions)
                    All_EES_regions = []
            else:
                with open(output_fn, 'a') as file:
                    file.writelines('\t'.join(j) + '\n' for j in All_EES_regions)
                    All_EES_regions = []
        i += 1
        #with open(output_fn, 'a') as file:
        #    file.writelines('\t'.join(i) + '\n' for i in All_EES_regions)

with open(output_fn, 'a') as file:
	file.writelines('\t'.join(i) + '\n' for i in All_EES_regions)

print "python script completed"
