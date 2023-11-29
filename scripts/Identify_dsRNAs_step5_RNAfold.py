###################################################################
# Script Name: Identify_dsRNAs_step5_RNAfold.py
# Description: This script extracts the sequence and saves the percent basepaired 
#              and the free energy for determining high confidence dsRNAs.
# Args_in: All ES (Editing Site) sites and a minimum number of people (This is five people).
# Args_out: Outputs a file of EERs (Editing Enriched Regions) and their RNA fold metrics.
# Author: Mudra Choudhury
# Email: 
# Date: 02/2019
# Last Edited: 
###################################################################

# Import necessary libraries
from collections import defaultdict
import sys
import pickle
import random
import csv
import os
import argparse
import re
from itertools import groupby
from Bio.Seq import Seq

# Functions defined in the script:
# 1. get_longest_dsRNA: Finds the longest dsRNA with less than a specified percentage of mismatches.
# 2. get_longest_stem: Identifies the longest stem structure with less than a specified percentage of mismatches.
# 3. features_of_annotations: Analyzes the RNA structure and calculates various features like percentage of bulge, stem, loop, etc.
# 4. updateMinorSpace: Updates the regions of the RNA that are not part of the stem-loop structure.
# 5. annotateFold: Annotates the RNA fold structure, identifying different structural elements like loops, stems, bulges, etc.
# 6. findOccurrences: Helper function to find occurrences of a character in a string.

def get_longest_dsRNA(grouped_annot, max_mismatch_perc):
    #Get longest dsRNA with less than max% mismatches:
    new_dsRNA_length = 0
    new_dsRNA_bulge_length = 0
    new_dsRNA_mismatch_perc = float(100)
    longest_dsRNA_length = 0
    longest_dsRNA_mistmatch_perc = float(0)
    i = 0
    while i < len(grouped_annot):
        #print grouped_annot[i][0], grouped_annot[i][1]
        if i == 0:
            if grouped_annot[i][0] == "ds":
                #print "In test1"
                new_dsRNA_length += grouped_annot[i][1]
                i += 1
            else:
                #print "In test2"
                i += 1
        else:
            if grouped_annot[i-1][0] == "ds" and grouped_annot[i][0] == "bulge" and grouped_annot[i+1][0] == "ds":
                #print "In test3"
                new_dsRNA_bulge_length += grouped_annot[i][1]
                new_dsRNA_length += grouped_annot[i + 1][1]
                #print new_dsRNA_bulge_length
                #print new_dsRNA_length
                i += 2
            elif grouped_annot[i-1][0] == "ds" and grouped_annot[i][0] == "bulge" and grouped_annot[i+1][0] != "ds":
                #print "In test4"
                new_dsRNA_total_length = new_dsRNA_length + new_dsRNA_bulge_length
                new_dsRNA_mismatch_perc = float(new_dsRNA_bulge_length)/float(new_dsRNA_total_length)
                #print new_dsRNA_mismatch_perc
                if new_dsRNA_total_length >= longest_dsRNA_length and new_dsRNA_mismatch_perc <= max_mismatch_perc:
                    longest_dsRNA_length = new_dsRNA_total_length
                    longest_dsRNA_mistmatch_perc = new_dsRNA_mismatch_perc
                new_dsRNA_bulge_length, new_dsRNA_length = 0, 0
                i += 2
            elif grouped_annot[i-1][0] == "ds" and grouped_annot[i][0] != "bulge":
                #print "In test5"
                new_dsRNA_total_length = new_dsRNA_length + new_dsRNA_bulge_length
                new_dsRNA_mismatch_perc = float(new_dsRNA_bulge_length)/float(new_dsRNA_total_length)
                #print new_dsRNA_mismatch_perc
                if new_dsRNA_total_length >= longest_dsRNA_length and new_dsRNA_mismatch_perc <= max_mismatch_perc:
                    longest_dsRNA_length = new_dsRNA_total_length
                    longest_dsRNA_mistmatch_perc = new_dsRNA_mismatch_perc
                new_dsRNA_bulge_length, new_dsRNA_length = 0, 0
                i += 1
            elif grouped_annot[i-1][0] != "ds" and grouped_annot[i][0] == "ds":
                #print "In test6"
                new_dsRNA_length += grouped_annot[i][1]
                i += 1
            elif grouped_annot[i-1][0] != "ds" and grouped_annot[i][0] != "ds":
                #print "In test7"
                i += 1
    #Just in case the last segment is a singular dsRNA with no mismatches, calculate the length and mismatches of the last one
    if new_dsRNA_length != 0:
        new_dsRNA_total_length = new_dsRNA_length + new_dsRNA_bulge_length
        new_dsRNA_mismatch_perc = float(new_dsRNA_bulge_length)/float(new_dsRNA_total_length)
        if new_dsRNA_total_length >= longest_dsRNA_length and new_dsRNA_mismatch_perc <= max_mismatch_perc:
            longest_dsRNA_length = new_dsRNA_total_length
            longest_dsRNA_mistmatch_perc = new_dsRNA_mismatch_perc
    longest_dsRNA_info = [longest_dsRNA_length, round(longest_dsRNA_mistmatch_perc,2)]
    #print longest_dsRNA_info
    return longest_dsRNA_info

def get_longest_stem(grouped_annot, max_mismatch_perc):
    #Get longest stem with less than max% mismatches:
    new_stem_length = 0
    new_stem_bulge_length = 0
    new_stem_mismatch_perc = float(100)
    longest_stem_length = 0
    longest_stem_mistmatch_perc = float(0)
    i = 0
    while i < len(grouped_annot):
        #print i
        #print grouped_annot[i][0], grouped_annot[i][1]
        if i == 0:
            if grouped_annot[i][0] == "stem":
                #print "In stem test1"
                new_stem_length += grouped_annot[i][1]
                i += 1
            else:
                #print "In stem test2"
                i += 1
        else:
            if grouped_annot[i-1][0] == "stem" and grouped_annot[i][0] == "bulge" and grouped_annot[i+1][0] == "stem":
                #print "In test3"
                new_stem_bulge_length += grouped_annot[i][1]
                new_stem_length += grouped_annot[i + 1][1]
                #print new_stem_bulge_length
                #print new_stem_length
                i += 2
            elif grouped_annot[i-1][0] == "stem" and grouped_annot[i][0] == "bulge" and grouped_annot[i+1][0] != "stem":
                #print "In test4"
                new_stem_total_length = new_stem_length + new_stem_bulge_length
                new_stem_mismatch_perc = float(new_stem_bulge_length)/float(new_stem_total_length)
                #print new_stem_mismatch_perc
                if new_stem_total_length >= longest_stem_length and new_stem_mismatch_perc <= max_mismatch_perc:
                    longest_stem_length = new_stem_total_length
                    longest_stem_mistmatch_perc = new_stem_mismatch_perc
                new_stem_bulge_length, new_stem_length = 0, 0
                i += 2
            elif grouped_annot[i-1][0] == "stem" and grouped_annot[i][0] != "bulge":
                #print "In test5"
                new_stem_total_length = new_stem_length + new_stem_bulge_length
                new_stem_mismatch_perc = float(new_stem_bulge_length)/float(new_stem_total_length)
                #print new_stem_mismatch_perc
                if new_stem_total_length >= longest_stem_length and new_stem_mismatch_perc <= max_mismatch_perc:
                    longest_stem_length = new_stem_total_length
                    longest_stem_mistmatch_perc = new_stem_mismatch_perc
                new_stem_bulge_length, new_stem_length = 0, 0
                i += 1
            elif grouped_annot[i-1][0] != "stem" and grouped_annot[i][0] == "stem":
                #print "In test6"
                new_stem_length += grouped_annot[i][1]
                i += 1
            elif grouped_annot[i-1][0] != "stem" and grouped_annot[i][0] != "stem":
                #print "In test7"
                i += 1
     #Just in case the last segment is a singular stem with no mismatches, calculate the length and mismatches of the last one
    if new_stem_length != 0:
        new_stem_total_length = new_stem_length + new_stem_bulge_length
        new_stem_mismatch_perc = float(new_stem_bulge_length)/float(new_stem_total_length)
        if new_stem_total_length >= longest_stem_length and new_stem_mismatch_perc <= max_mismatch_perc:
            longest_stem_length = new_stem_total_length
            longest_stem_mistmatch_perc = new_stem_mismatch_perc
    longest_stem_info = [longest_stem_length, round(longest_stem_mistmatch_perc,2)]
    #print longest_stem_info
    return longest_stem_info

def features_of_annotations(annotated_string, max_mismatch_perc):
    total_length = float(len(annotated_string))
    ### Get percentage of each type
    bulge_perc = round((float(annotated_string.count("bulge"))/total_length), 2)
    stem_perc = round((float(annotated_string.count("stem"))/total_length), 2)
    ss_perc = round((float(annotated_string.count("ss"))/total_length), 2)
    ds_perc = round((float(annotated_string.count("ds"))/total_length), 2)
    loop_perc = round((float(annotated_string.count("loop"))/total_length), 2)

    #print "DS: ", ds_perc, "Stem: ", stem_perc, "Bulge: ", bulge_perc, "Loop: ", loop_perc, "ss: ", ss_perc
    ### Get length of each type:
    grouped_annot = [(k, sum(1 for i in g)) for k,g in groupby(annotated_string)]

    #bulge_to_ds_string = ['ds' if word == 'bulge' for word in words]
    #print grouped_annot
    #print bulge_to_ds_string

    #Get total number of regions for each type. Think about this a little bit... how can we weight these numbers by length of the dsRNA?
    total_loop_segments, total_stem_segments, total_ds_segments, total_ss_segments, total_bulge_segments = 0, 0, 0, 0, 0
    for segment in grouped_annot:
        if segment[0] == "loop":
            total_loop_segments += 1
        if segment[0] == "stem":
            total_stem_segments += 1
        if segment[0] == "ds":
            total_ds_segments += 1
        if segment[0] == "ss":
            total_ss_segments += 1
        if segment[0] == "bulge":
            total_bulge_segments += 1
    #print "Total length: ", total_length, "Total loop segments: ", total_loop_segments, loop_perc, "Total stem segments: ", total_stem_segments, stem_perc, "Total ds segments: ", total_ds_segments, ds_perc, "Total ss segments: ", total_ss_segments, ss_perc, "Total bulge segments: ", total_bulge_segments, bulge_perc

    #Get longest dsRNA with less than 10% mismatches:
    longest_dsRNA_info = get_longest_dsRNA(grouped_annot, max_mismatch_perc)
    longest_stem_info = get_longest_stem(grouped_annot, max_mismatch_perc)
    #print "Longest dsRNA and percent mismatches: ", longest_dsRNA_info
    #print "Longest stem and percent mismatches: ", longest_stem_info
    all_info = [str(longest_dsRNA_info[0]), str(longest_dsRNA_info[1]), str(longest_stem_info[0]), str(longest_stem_info[1]), str(total_ds_segments), str(ds_perc), str(total_stem_segments), str(stem_perc), str(total_bulge_segments), str(bulge_perc), str(total_loop_segments), str(loop_perc), str(total_ss_segments), str(ss_perc)]
    #print all_info
    return all_info

def updateMinorSpace(stems, position):
        mSpace = []
        notSpace = []

        for stem in stems:
                for i in range(stem[0], stem[1] + 1):
                        notSpace.append(i)

        for i in range(0, position):
                if (i not in notSpace):
                        mSpace.append(i)

        return mSpace

def annotateFold(foldStructure):
    '''
    print foldStructure
    for i,j in enumerate(foldStructure):
            print i, j
    '''

    #create binding dict via stack
    bindingDict = {}

    for i, j in enumerate(foldStructure):
        bindingDict[i] = None

    stack = []
    for i, j in enumerate(foldStructure):
        if j == '(':
            #print "stacked"
            stack.append(i)
        elif j == ')':
            #print "popped"
            bp1 = stack.pop()
            bp2 = i
            bindingDict[bp1] = bp2
            bindingDict[bp2] = bp1

    #find all loops/have to make sure flanking nt are auto-binding
    loops = []
    p = re.compile('[(][.]{1,100}[)]')
    loopsIter = p.finditer(foldStructure)
    for m in loopsIter:
            loops.append(m.span())
    #print loops

    minorSpace = []
    stemLoops = []

    for loop in loops:
        i = loop[1]-1
        while True:
            if (i<len(foldStructure)) and (foldStructure[i]==')') and not((i - bindingDict[i] < 0) or (bindingDict[i] in minorSpace)):
                #keep track of last valid index(not None or '.')
                lastvalid = i
                #   print 'lastvalid'
                #   print lastvalid
            if i > len(foldStructure) - 1:
                #create boundary
                stemLoops.append([bindingDict[lastvalid], lastvalid])
                #print stemLoops, i
                minorSpace = updateMinorSpace(stemLoops, i) #Minor space are those regions not within a loop
                #print stemLoops
                #print minorSpace
                break
            if bindingDict[i]==None:
                #skip internal "."
                i += 1
                continue
            if (i - bindingDict[i] < 0) or (bindingDict[i] in minorSpace):
                #create boundary and designate stem loops from minorspace
                stemLoops.append([bindingDict[lastvalid], lastvalid])
                minorSpace = updateMinorSpace(stemLoops, i)
                #print minorSpace
                break
            else:
                i += 1

    #capture last nt into minor structure
    for x in range(i, len(foldStructure)):
        if x not in minorSpace:
            minorSpace.append(x)
            #print x

    #Search for dotted/non base pair regions
    strand = []
    a = re.compile('[.]{1,1000}')
    dotIter = a.finditer(foldStructure)
    for m in dotIter:
            strand.append(m.span())
    #print strand

    singlestrand = []
    bulge = []
    for dot in strand:
        left = dot[0]
        right = dot[1]
        #check boundaries for single strand
        if (left==0) or (right>len(foldStructure)-1):
            singlestrand.append((left,right-1))
            continue
        if (foldStructure[left-1]==')' and foldStructure[right]=='('):
            singlestrand.append((left,right-1))
            continue
        if ((foldStructure[left-1]==')' and foldStructure[right]==')') or (foldStructure[left-1]=='(' and foldStructure[right]=='(')):
            #check for bulges
            bulge.append((left,right-1))

    #remove single strands from minor space
    for ss in singlestrand:
        for i in range(ss[0], ss[1]+1):
            if i in minorSpace:
                minorSpace.remove(i)
    #remove bulges from minor space
    for bu in bulge:
        for j in range(bu[0], bu[1]+1):
            if j in minorSpace:
                minorSpace.remove(j)

    #whatever is left in the minor space is db...
    doubleSpace = []
    for i in minorSpace:
            doubleSpace.append(i)

    #Annotate the nucleotides
    aDict = {}

    loopVals = []
    for loop in loops:
        l = loop[0] + 1
        ll = loop[1] - 1
        loopVals.extend(range(l, ll))

    bVals = []
    for b in bulge:
        for i in range(b[0], b[1] + 1):
                bVals.append(i)

    for stemLoop in stemLoops:
        for i, j in enumerate(range(stemLoop[0], stemLoop[1] + 1)):
            #aDict.setdefault(j, []).append('stemLoop')
            if j in loopVals and j not in bVals:
                aDict.setdefault(j, []).append('loop')
            elif j not in bVals:
                aDict.setdefault(j, []).append('stem')

    for b in bulge:
        for i in range(b[0], b[1] + 1):
            aDict.setdefault(i, []).append('bulge')

    for s in singlestrand:
        for i in range(s[0], s[1] + 1):
            aDict.setdefault(i, []).append('ss')

    for ds in doubleSpace:
        aDict.setdefault(ds, []).append('ds')

    annotated_string = []
    for num in aDict:
        #annotated_string.append('%s:%s' % (num, ','.join(aDict[num])))
        annotated_string.append(",".join(aDict[num]))


    return annotated_string



# Main program execution starts here
if __name__ == "__main__":
    # Set up argument parser for command line inputs
    parser = argparse.ArgumentParser(description='Input Files')
    parser.add_argument('--i', metavar='input_bed_file', type=str, help='Bed file with EES regions of certain length')
    parser.add_argument('--o', metavar='output_bed_file', type=str, help='Output Bed file with EES regions and RNAfold info attached')
    args = parser.parse_args()

    # Assign arguments to variables
    input_fn = args.i
    out_fn = args.o
    max_mismatch_perc = 0.2 # Set maximum mismatch percentage

    # Load the whole hg19 genome from a pickle file
    genome_seq_fn = "/home/stran/hyperediting_pipeline_EX/hg19/genome.pickle"
    print "Loading the pickle of whole hg19 genome..."
    seq_pickle_open = open(genome_seq_fn, 'r')
    print "Read pickle"
    genome_pickle = pickle.load(seq_pickle_open)
    print "Done loading genome."


    def findOccurrences(s, ch):
        return [i for i, letter in enumerate(s) if letter == ch]


    random_int = random.randint(1, 100000)
    tmp_fn="".join(("temp_seq_",str(random_int),".fa"))

    All_lines = [["chromosome", "start", "end", "strand", "length", "MFE", "AMFE", "MFEI", "MEA", "Perc_bp_paired", "ENSMBL_MFE", "ENSMBL_D","longest_dsRNA_length", "longest_dsRNA_prop_mismatch", "longest_stem_length", "longest_stem_perc_mismatch", "total_ds_segments", "ds_prop", "total_stem_segments", "stem_prop", "total_bulge_segments", "bulge_prop", "total_loop_segments", "loop_prop", "total_ss_segments", "ss_prop"]]
    i = 0
    regions_fn=open(input_fn,"r")
    #regions_fn=(open(input_fn,"r")).readlines()
    #regions_fn = regions_fn[8500:9000] #for debugging
    for line in regions_fn:
        i=i+1
        if i % 1000 == 0:
            print "Iteration: ", i, "_________________"
        #if i > 5: #for debugging
        #    break
        split_line=line.split("\t")
        chrom, start, end = split_line[0], int(split_line[1]), int(split_line[2]) #EES region
        strand = (split_line[3]).rstrip()
        region_seq = genome_pickle[chrom][start:end]

        if strand == "-":
            seq = Seq(region_seq)
            rc_region_seq=str(seq.reverse_complement())
            region_seq = rc_region_seq

        with open(tmp_fn,'w') as out:
            out.write(region_seq)

        ##### Implement RNAfold
        RNAfold_string="".join(('~/software/ViennaRNA/bin/RNAfold --MEA <',tmp_fn))
        result=os.popen(RNAfold_string).read()
        #print result
        result_split = result.split("\n")

        ### Get maximum free energy and total length
        try:
            bp_result = result_split[1]
        except IndexError:
            os.remove(tmp_fn)
            continue
        bp_result_split = bp_result.split(" ")
        total_length = float(len(bp_result_split[0]))
        pairing_string = bp_result_split[0]
        free_energy_temp = bp_result_split[1]
        #print free_energy_temp
        #print free_energy_temp[1:-1]
        free_energy_temp2 = re.sub("[^0-9^.^-]","",free_energy_temp)
        if free_energy_temp2 == "":
            free_energy = float("nan")
        else:
            free_energy = float(free_energy_temp2)
        #free_energy = float(free_energy_temp[1:-1])
        #print free_energy

        ### Get ensembl free energy
        line3_result = result_split[2]
        line3_result_split = line3_result.split(" ")
        ensembl_free_energy_temp = line3_result_split[1]
        ensembl_free_energy = ensembl_free_energy_temp[1:-1]
        #print ensembl_free_energy

        ### Get distance to ensembl
        line4_result = result_split[3]
        line4_result_split = line4_result.split(" ")
        dist_to_ensembl_temp = (line4_result_split[2]).split("=")
        if len(dist_to_ensembl_temp) > 1:
            dist_to_ensembl_temp2 = dist_to_ensembl_temp[1]
            dist_to_ensembl = dist_to_ensembl_temp2[0:-1]
        else:
            dist_to_ensembl = float("nan")
        #print dist_to_ensembl

        ### Get MEA: maximum error accuracy
        line5_result = result_split[4]
        line5_result_split = line5_result.split(" ")
        MEA_temp = (line5_result_split[2]).split("=")
        if len(MEA_temp) > 1:
            MEA_temp2 = MEA_temp[1]
            MEA = MEA_temp2[0:-1]
        else:
            MEA = float("nan")
        #print MEA

        if free_energy != float("nan"):
            AMFE = free_energy/total_length
        else:
            AMFE = float("nan")


        #unpaired_length = findOccurrences(bp_result_split[0], '[.]')
        unpaired_length = pairing_string.count('.')
        unpaired_length = float(unpaired_length)
        paired_length = total_length - unpaired_length
        paired_perc = (total_length - unpaired_length) / total_length

        #### Get GC content and MFEI
        GC_count = float(result_split[0].count("G") + result_split[0].count("g") + result_split[0].count("C") + result_split[0].count("c"))
        MFEI = (100 * free_energy)/total_length/GC_count
        #print MFEI

        ########################### This section of the script annotates the features of the sequence structure from mfold ##############################
        annotated_string = annotateFold(pairing_string)
        #print annotated_string
        annotation_info = features_of_annotations(annotated_string, max_mismatch_perc)

        os.remove(tmp_fn)

        new_line = [chrom, str(start), str(end), str(strand), str(total_length), str(free_energy), str(AMFE), str(MFEI), str(MEA), str(paired_perc), str(ensembl_free_energy), str(dist_to_ensembl)] + annotation_info
        #print new_line
        All_lines.append(new_line)

    #print All_lines
    with open(out_fn, 'w') as file:
        file.writelines('\t'.join(i) + '\n' for i in All_lines)
    print "Step 5 job completed"
