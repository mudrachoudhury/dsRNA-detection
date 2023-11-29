###################################################################
# Script Name: Identify_dsRNAs_step4_length.R
# Description: This script plots dsRNA lengths according to merge distances 
#              selected to help identify the proper cutoffs for dsRNA lengths.
# Args_in: input bed file of merged EES regions.
# Args_out: plots of lengths distribution.
# Author: Mudra Choudhury
# Email: [email]
# Date: 02/24/2019
# Last Edited: 
###################################################################

# Load required libraries
library("optparse")
library(data.table)
library(ggplot2)
library(ggrepel)
library(dplyr)

# Define command line options
option_list = list(
  make_option(c("-i", "--input_fn"), type="character", default=NULL, 
              help="input bed file of dsRNA regions", metavar="character"),
  make_option(c("-l", "--length"), type="character", default=NULL, 
              help="length cutoff", metavar="character"),  
  make_option(c("-o", "--output_fn"), type="character", 
              help="output bed file of dsRNA regions within length", metavar="character"))

# Create option parser and parse arguments
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

# Assign arguments to variables
in_fn = opt$input_fn
length = as.numeric(opt$length)

# Define merge distances and other parameters for the analysis
merge_dist_list = c(50, 100, 150, 200, 400, 600, 800, 1000, 1500, 2000)
indv_list = c(5)
min_length = as.numeric(300)
max_length = as.numeric(4000)

# Prepare output PDF file name
out_pdf = paste0("../plots/Pick_length_cutoff_", as.character(min_length), "-", as.character(max_length), "_large_merge_distances.pdf")

# Begin PDF output
pdf(out_pdf, width = 6, height = 5)
for (indv in indv_list){
    all_region_lengths = data_frame()
    final_region_lengths = data_frame()
    
    for (distance in merge_dist_list){
        # Generate input filename based on distance
        in_fn= paste0("../data/Merged_", distance, "_dist_All_EES_50_window_regions_w_min_3_ES_per_window_using_REDI_2020_combined.bed")

        # Process input bed file
        bed_file = read.table(in_fn, header=FALSE, sep = "\t")
        lengths = bed_file$V3 - bed_file$V2
        all_region_lengths = rbind(all_region_lengths, lengths)
        
        # Filter regions based on length and write to file
        keep_rows = which(lengths < max_length & lengths > min_length)
        output_bed_file = bed_file[keep_rows,]
        new_lengths = output_bed_file$V3 - output_bed_file$V2
        final_region_lengths = rbind(final_region_lengths, new_lengths)
        out_fn = paste0("../data/", as.character(distance), "_length_cutoff.bed")
        write.table(output_bed_file, file = out_fn, sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
    }
    
    # Plot all region lengths
    region_lengths_df = as.data.frame(t(all_region_lengths))
    colnames(region_lengths_df) <- sapply(merge_dist_list, as.character)
    w.plot <- melt(region_lengths_df) 
    plot = ggplot(aes(x=value, y=stat(count), colour=variable), data=w.plot) +
            geom_density() + 
            ggtitle("All region lengths at diff merge distances") + 
            labs(color = "Merge distance cutoff", x = "Region Length") + 
            theme_classic()
    print(plot)
    
    # Plot final region lengths
    final_region_lengths_df = as.data.frame(t(final_region_lengths))
    colnames(final_region_lengths_df) <- sapply(merge_dist_list, as.character)
    n.plot <- melt(final_region_lengths_df) 
    plot2 = ggplot(aes(x=value, y=stat(count), colour=variable), data=n.plot) +
            geom_density() + 
            ggtitle("Final region lengths for at diff merge distances") + 
            labs(color = "Merge distance cutoff", x = "Region Length") + 
            theme_classic()
    print(plot2)
}

# End PDF output
dev.off()

