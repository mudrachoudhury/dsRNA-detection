###################################################################
#Script Name    : Identify_dsRNAs_step2_writeChunks_updated.py
#Description    : Go through 50 bp windows and get those with at least 3 editing sites
#Args_in        : All ES sites and a minimum number of people
#Args_out   : Output file contains windows with at least 3 editing sites
#Author     : Mudra Choudhury
#Email      : mudrachoudhury3@gmail.com
#Date       :
#Last Edited    : 4/24/19
###################################################################

# Import necessary libraries
import os
import collections
import argparse
from collections import defaultdict
from itertools import islice
import sys

# Set up argument parser for command line inputs
parser = argparse.ArgumentParser(description='Input Files')
parser.add_argument('-i', metavar='input_file', type=str, help='input file with all editing sites')
parser.add_argument('-n', metavar='min_es', type=str, help='minimum number of editing sites per window')
parser.add_argument('-s', metavar='window_size', type=str, help='size of sliding window for each editing site')
parser.add_argument('-c', metavar='chromosome', type=str, help='Chromosome that the script will be analyzing, ex: chr12')

# Parse arguments from command line
args = parser.parse_args()
input_fn = args.i
min_es_in_window = args.n
window_size = args.s
chromosome_argument = args.c

# Process the input filename
split_input = input_fn.split("_")

# Construct the output filename
output_fn = "".join(("../data/All_EES_",window_size,"_window_regions_w_min_",min_es_in_window,"_ES_per_window_using_REDI_2020_",chromosome_argument,"_updated.bed"))

# Define a function for sliding window operation
def slidingWindow(sequence,winSize,step=1):
    """Returns a generator that will iterate through the defined chunks of input sequence.
    Input sequence must be iterable."""
    # Verify inputs and throw exceptions if invalid
    try: it = iter(sequence)
    except TypeError:
        raise Exception("**ERROR** sequence must be iterable.")
    if not ((type(winSize) == type(0)) and (type(step) == type(0))):
        raise Exception("**ERROR** type(winSize) and type(step) must be int.")
    if step > winSize:
        raise Exception("**ERROR** step must not be larger than winSize.")
    if winSize > len(sequence):
        raise Exception("**ERROR** winSize must not be larger than sequence length.")
    
    # Calculate number of chunks and yield them
    numOfChunks = ((len(sequence)-winSize)/step)+1
    for i in range(0,numOfChunks*step,step):
        yield sequence[i:i+winSize]

# Function to get Editing Enriched Sites (EES) windows
def get_EES_windows(chromosome, positions, es_set, window_size, extra_info):
    """Extracts windows with a minimum number of editing sites."""
    ees_regions = []
    min_es = int(min_es_in_window)
    step = 1
    windows = slidingWindow(positions, window_size, step)
    for window in windows:
        intersection = es_set.intersection(window)
        if len(intersection) >= min_es:
            start, end = (window[0]-1), window[-1]
            ees_line = [chromosome, str(start), str(end), "NA", "NA"] + extra_info
            ees_regions.append(ees_line)
    return ees_regions

# Define chromosomes to process
chromosomes = [chromosome_argument]

# Initialize processing
start_print = "Starting script...\n"
sys.stderr.write(start_print)

# Main processing loop for each chromosome
for chrom in chromosomes:
    All_EES_regions = []
    print(chrom)
    chrom_positions_es_set = set()
    n = 0
    input_file = open(input_fn, "r")
    for line in input_file:
        if n == 0:  # Skip header
            n += 1
            continue
        else:
            split_line = line.split("\t")
            chrom_positions_es_set.add(int(split_line[1]))
    
    i = 0
    input_file = open(input_fn, "r")
    for es_line in input_file:
        if i == 0:  # Skip header
            i += 1
            continue
        split_es_line = ((es_line).rstrip()).split("\t")
        position = int(split_es_line[1])
        extra_info = [split_es_line[2], split_es_line[3], split_es_line[4]]
        start_window, end_window = (position - int(window_size)), (position + int(window_size))
        window_seq = range(start_window,end_window)
        list_of_ees_regions = get_EES_windows(chrom, window_seq, chrom_positions_es_set, int(window_size), extra_info)
        All_EES_regions += list_of_ees_regions
        if i % 1000 == 0:
            log_print = "processed line: " + str(i) + "\n"
            sys.stderr.write(log_print)
            with open(output_fn, 'a' if i != 1000 else 'w') as file:
                file.writelines


print "python script completed"
