###################################################################
#Script Name	: Identify_dsRNAs_step4_length.R
#Description	: 
#Args_in	: input bed file of merged EES regions
#Args_out	:
#Author		: Mudra Choudhury
#Email		: mudrachoudhury3@gmail.com
#Date		: 02/24/2019
#Last Edited    :
###################################################################
args = commandArgs(trailingOnly = TRUE)
library("optparse")
library(data.table)
library(ggplot2)
library(ggrepel)
library(dplyr)

########Take in arguments###############################
option_list = list(
  make_option(c("-i", "--input_fn"), type="character", default=NULL, 
              help="input bed file of dsRNA regions", metavar="character"),
  make_option(c("-l", "--length"), type="character", default=NULL, 
              help="length cutoff", metavar="character"),  
  make_option(c("-o", "--output_fn"), type="character", 
              help="output bed file of dsRNA regions within length", metavar="character"))

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
########################################################

args = commandArgs(trailingOnly = TRUE)
in_fn = opt$input_fn
length = as.numeric(opt$length)

merge_dist_list = c(50, 100, 150, 200, 400,600,800,1000,1500,2000)
indv_list = c(5)
min_length = as.numeric(300)
max_length = as.numeric(4000)
min_paste = as.character(min_length)
max_paste = as.character(max_length)
out_pdf = paste0("../plots/Pick_length_cutoff_", min_paste, "-",max_paste,"_large_merge_distances.pdf")


pdf(out_pdf, width = 6, height = 5)
for (indv in indv_list){
  all_region_lengths = data_frame()
  final_region_lengths = data_frame()
  for (distance in merge_dist_list){
    #print(indv)
    #print(distance)
    ############################ Input files ###############################
    #If you want to take in arguments then remove these lines and the for loop above.
    in_fn= paste0("../data/Merged_",distance,"_dist_All_EES_50_window_regions_w_min_3_ES_per_window_using_REDI_2020_combined.bed")
    #########################################################################

    input_split = unlist(strsplit(in_fn,"[/]"))
    out_prefix_temp1 = input_split[3]
    out_prefix_temp2 = unlist(strsplit(out_prefix_temp1,"[.]"))
    output_prefix=paste0(input_split[1],"/",input_split[2],"/",out_prefix_temp2[1])
    out_fn= paste0(output_prefix,"_",as.character(min_length),"-",as.character(max_length),"_length_cutoff.bed")
    print(out_fn)

    # If you want the length to be within a certain range
    #bp_range = (length*.10)
    #max_length = length + bp_range
    #min_length = length - bp_range


    bed_file=read.table(in_fn, header=FALSE, sep = "\t")
    lengths = bed_file$V3 - bed_file$V2
    all_region_lengths = rbind(all_region_lengths, lengths)
    
    #keep_rows = which(lengths > length)
    keep_rows = which(lengths < max_length & lengths > min_length)
    output_bed_file = bed_file[keep_rows,]
    new_lengths = output_bed_file$V3 - output_bed_file$V2
    final_region_lengths = rbind(final_region_lengths, new_lengths)

    write.table(output_bed_file, file = out_fn, sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
    
    #print("Done")

    }
    
    region_lengths_df=as.data.frame(t(all_region_lengths))
    colnames(region_lengths_df) <- sapply(merge_dist_list, as.character)
    library(ggplot2)
    library(reshape2)
    w.plot <- melt(region_lengths_df) 
    
    p <- ggplot(aes(x=value, stat(count), colour=variable), data=w.plot)
    plot = p + geom_density() + ggtitle(paste0("All region lengths at diff merge distances")) + labs(color = "Merge distance cutoff", x = "Region Length") + theme_classic()
    print(plot)
    
    final_region_lengths_df=as.data.frame(t(final_region_lengths))
    colnames(final_region_lengths_df) <- sapply(merge_dist_list, as.character)
    library(ggplot2)
    library(reshape2)
    n.plot <- melt(final_region_lengths_df) 
    
    p2 <- ggplot(aes(x=value, stat(count), colour=variable), data=n.plot)
    plot2 = p2 + geom_density() + ggtitle(paste0("Final region lengths for at diff merge distances")) + labs(color = "Merge distance cutoff", x = "Region Length") + theme_classic()
    print(plot2)
    

}
dev.off()
