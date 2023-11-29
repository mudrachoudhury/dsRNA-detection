###################################################################
#Script Name	: Identify_dsRNAs_step6_RNAFold_cutoffs.R
#Description	: #plot the MFE vs MEA and apply cutoff to RNA_fold output regions. plot MFEI value
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
library(cowplot)

input_fns = args[1] #file containing paths to all input files
MFE_cutoff = as.numeric(args[2])
MEA_cutoff = as.numeric(args[3])
#MFE_cutoff = -200
#MEA_cutoff = 300
stem_ds_length_cutoff = 200
stem_ds_mistmatch_cutoff = 0.20
#AMFE_cutoff = -0.4 #for right now we don't need to apply an AMFE cutoff. The script can be changed to apply it at the end if needed.
max_length = 5000
output_pdf = "../plots/Pick_RNAfold_thresholds_stem_length_and_mm.pdf"
axis_title_size = 10
axis_size = 10
title_size = 12

#input_files = c("../data/Merged_100/Merged_100_combined_output_w_RNA_fold.bed",
#                "../data/Merged_200/Merged_200_combined_output_w_RNA_fold.bed",
#                "../data/Merged_400/Merged_400_combined_output_w_RNA_fold.bed",
#                "../data/Merged_600/Merged_600_combined_output_w_RNA_fold.bed",
#		"../data/Merged_800/Merged_800_combined_output_w_RNA_fold.bed",
#		"../data/Merged_1500/Merged_1500_combined_output_w_RNA_fold.bed",
#		"../data/Merged_2000/Merged_2000_combined_output_w_RNA_fold.bed",
#                "../data/Merged_1000/Merged_1000_combined_output_w_RNA_fold.bed")

input_files = c("../data/Merged_1000_EERs/Merged_1000_combined_output_w_RNA_fold.bed")
############ Plot the frequency of the different dsRNA statstics (like AMFE, MFE, stem length) to pick a cutoff

pdf(output_pdf, width=5, height=3)

#input_files = read.table(input_fns, sep = "\n", header = FALSE) #so far we are not reading in the files, I have hard coded them
for (file_name in input_files){
  print(file_name)
  split_fn = unlist(strsplit(file_name, "/"))
  merge_dist = unlist(strsplit(split_fn[3], "_"))[2]
  out_prefix=unlist(strsplit(split_fn[4],"[.]"))[1]
  #out_bed_fn=paste0("../",num_indv,"/",out_prefix,"_DsStemLength_",as.character(stem_ds_length_cutoff),"_mmProp_",as.character(stem_ds_mistmatch_cutoff),"_rank_AMFE_",as.character(AMFE_cutoff),"_cutoff.bed") #uncomment if applying an AFE cutoff
  out_bed_fn=paste0("../data/Merged_", merge_dist,"_EERs/",out_prefix,"_DsStemLength_",as.character(stem_ds_length_cutoff),"_mmProp_",as.character(stem_ds_mistmatch_cutoff),"_rank_AMFE.bed") #comment if applying AMFE cutoff
  file = as.data.frame(read.table(file_name, sep = "\t", header = TRUE))
  title_prefix = paste0("Merged ", merge_dist," distance")
  file$MEA = as.numeric(as.character(file$MEA))
  file$MFE = as.numeric(as.character(file$MFE))
  file$MFEI = as.numeric(as.character(file$MFEI))
  file$longest_dsRNA_length = as.numeric(as.character(file$longest_dsRNA_length))
  file$longest_dsRNA_prop_mismatch = as.numeric(as.character(file$longest_dsRNA_prop_mismatch))
  file$longest_stem_length = as.numeric(as.character(file$longest_stem_length))
  file$longest_stem_perc_mismatch = as.numeric(as.character(file$longest_stem_perc_mismatch))
  file$length = as.numeric(as.character(file$length))
  file$AMFE = as.numeric(as.character(file$AMFE))
  
  
  #####Get indeces with longer stem and indeces with longer ds regions for more accurate correlation plotting
  ds_indeces = which(file$longest_dsRNA_length > stem_ds_length_cutoff & file$longest_dsRNA_prop_mismatch < stem_ds_mistmatch_cutoff)
  stem_indeces = which(file$longest_stem_length >= stem_ds_length_cutoff & file$longest_stem_perc_mismatch < stem_ds_mistmatch_cutoff)
  new_file = file[sort(unique(union(ds_indeces, stem_indeces))),]
  
  ###### This is a file before applying the stem/ds length cutoff. The before and after is for plotting and looking at the number of dsRNAs we lose by applying this cutoff.
  before_cutoff_longer_ds_index = which(file$longest_dsRNA_length > file$longest_stem_length)
  before_cutoff_longer_stem_index = which(file$longest_stem_length > file$longest_dsRNA_length)
  before_cutoff_longer_ds_RNAs = file[before_cutoff_longer_ds_index,]
  before_cutoff_longer_stem_RNAs = file[before_cutoff_longer_stem_index,]
  before_cutoff_longer_stem_RNAs$longest_ds_or_stem_length = before_cutoff_longer_stem_RNAs$longest_stem_length
  before_cutoff_longer_stem_RNAs$longest_ds_or_stem_mm_prop = before_cutoff_longer_stem_RNAs$longest_stem_perc_mismatch
  before_cutoff_longer_ds_RNAs$longest_ds_or_stem_length = before_cutoff_longer_ds_RNAs$longest_dsRNA_length
  before_cutoff_longer_ds_RNAs$longest_ds_or_stem_mm_prop = before_cutoff_longer_ds_RNAs$longest_dsRNA_prop_mismatch
  file_wo_cutoff = rbind(before_cutoff_longer_stem_RNAs, before_cutoff_longer_ds_RNAs)
  
  
  
  ###### After stem/ds length cutoff
  longer_ds_index = which(new_file$longest_dsRNA_length > new_file$longest_stem_length)
  longer_stem_index = which(new_file$longest_stem_length > new_file$longest_dsRNA_length)
  longer_ds_RNAs = new_file[longer_ds_index,]
  longer_stem_RNAs = new_file[longer_stem_index,]
  longer_stem_RNAs$longest_ds_or_stem_length = longer_stem_RNAs$longest_stem_length
  longer_stem_RNAs$longest_ds_or_stem_mm_prop = longer_stem_RNAs$longest_stem_perc_mismatch
  longer_ds_RNAs$longest_ds_or_stem_length = longer_ds_RNAs$longest_dsRNA_length
  longer_ds_RNAs$longest_ds_or_stem_mm_prop = longer_ds_RNAs$longest_dsRNA_prop_mismatch
  new_file_w_cutoff = rbind(longer_stem_RNAs,longer_ds_RNAs)
  
  
  ### correlation plots
  
  ### MEA vs MFE
  corr_plot <- ggplot(file,aes(x=MEA,y=MFE)) + geom_point(alpha = 0.2) + ggtitle(paste0(title_prefix," wo cutoff")) +
    labs(x="Maximum Expected Accuracy", y="Minimum Free Energy") + 
    theme(axis.text=element_text(size=axis_size), axis.title=element_text(size=axis_title_size,face="bold"), title = element_text(size=title_size)) + geom_smooth(method='lm')
  
  ### length vs AMFE
  corr_plot2 <- ggplot(file,aes(x=length,y=AMFE)) + geom_point(alpha = 0.2) +
    labs(x="dsRNA length wo cutoff", y="Adjusted Minimum Free Energy (AMFE)", subtitle = paste0(title_prefix)) + 
    theme(axis.text=element_text(size=axis_size), axis.title=element_text(size=axis_title_size, face = "bold"), title = element_text(size=title_size)) + stat_smooth(method = "loess", formula = y ~ x, size = 1)
  
  ### length vs MFE
  corr_plot3 <- ggplot(file,aes(x=length,y=MFE)) + geom_point(alpha = 0.2) +
    labs(x="dsRNA length wo cutoff", y="Minimum Free Energy (MFE)") + 
    theme(axis.text=element_text(size=axis_size), axis.title=element_text(size=axis_title_size,face="bold"), title = element_text(size=title_size)) + stat_smooth(method = "loess", formula = y ~ x, size = 1)
  
  
  #another plot of AMFE old vs AMFE new
  corr_plot2.5 <- ggplot(new_file_w_cutoff,aes(x=length,y=AMFE)) + geom_point(alpha = 0.2) +
    labs(x="total dsRNA length of those w stem/ds cutoff", y="Adjusted Minimum Free Energy (AMFE)") + 
    theme(axis.text=element_text(size=axis_size), axis.title=element_text(size=axis_title_size, face = "bold"), title = element_text(size=title_size)) + stat_smooth(method = "loess", formula = y ~ x, size = 1)
  
  
    ## original longest ds or stem region vs proportion of mismatches
  corr_plot4 <- ggplot(file_wo_cutoff,aes(x=longest_ds_or_stem_length,y=longest_ds_or_stem_mm_prop)) + geom_point(alpha = 0.2) + 
    labs(x="longest ds/stem length wo cutoff", y="mismatch proportion", subtitle = paste0(title_prefix)) + 
    theme(axis.text=element_text(size=axis_size), axis.title=element_text(size=title_size,face="bold"),title = element_text(size=title_size))
  
  ## longest ds or stem region after applying bp and mm cutoffs
  corr_plot4.5 <- ggplot(new_file_w_cutoff, aes(x=longest_ds_or_stem_length,y=longest_ds_or_stem_mm_prop)) + geom_point(alpha = 0.2) +
    labs(x="longest ds/stem length w cutoff", y="mismatch proportion") + 
    theme(axis.text=element_text(size=axis_size), axis.title=element_text(size=axis_title_size,face="bold"),title = element_text(size=title_size))

  ## longest length vs AMFE wo cutoff
  corr_plot5 <- ggplot(file_wo_cutoff,aes(x=longest_ds_or_stem_length,y=AMFE)) + geom_point(alpha = 0.2) +
    labs(x="longest ds/stem length wo cuoff", y="AMFE", subtitle = paste0(title_prefix)) + 
    theme(axis.text=element_text(size=axis_size), axis.title=element_text(size=axis_title_size,face="bold"), title = element_text(size=title_size)) + stat_smooth(method = "lm")
  
  ## longest length vs AMFE w cutoff
  corr_plot5.5 <- ggplot(new_file_w_cutoff,aes(x=longest_ds_or_stem_length,y=AMFE)) + geom_point(alpha = 0.2) +
    labs(x="longest ds/stem length w cutoff", y="AMFE") + 
    theme(axis.text=element_text(size=axis_size), axis.title=element_text(size=axis_title_size,face="bold"), title = element_text(size=title_size)) + stat_smooth(method = "lm")
  
  
  ##### histograms
  MFE_hist <- hist(file$MFE, cex.axis=0.5, cex.main=0.5, cex.lab=0.5, col="blue", breaks=20, main=paste0("Freq of MFE: ", title_prefix), xlab = "Minimum Free Energy")
  MEA_hist <- hist(file$MEA, cex.axis=0.5, cex.main=0.5, cex.lab=0.5, col="red", breaks=20, xlim = c(0,1500), main=paste0("Freq of MEA: ", title_prefix), xlab = "Maximum Expected Accuracy")
  AMFE_hist <- hist(file$AMFE, cex.axis=0.5, cex.main=0.5, cex.lab=0.5, col="purple", breaks=20, main=paste0("Freq of AMFE: ", title_prefix), xlab = "AMFE")
  MFEI_hist <- hist(file$MFEI, cex.axis=0.5, cex.main=0.5, cex.lab=0.5, col="purple", breaks=20, main=paste0("Freq of MFEI: ", title_prefix), xlab = "MFEI")
  
  # print(corr_plot)
  # print(plot_grid(corr_plot3, corr_plot2))
  # print(plot_grid(corr_plot2, corr_plot2.5))
  # print(plot_grid(corr_plot4, corr_plot4.5))
  print(plot_grid(corr_plot5, corr_plot5.5))
  # print(MFE_hist)
  # print(AMFE_hist)
  # print(MEA_hist)
  
  #### rank the dsRNAs according to AMFE:
  #new_file_w_AMFE_cutoff = new_file_w_cutoff[which(new_file_w_cutoff$AMFE <= AMFE_cutoff),] #if you would like to apply an AMFE cutoff then use this line
  new_file_w_AMFE_cutoff <- new_file_w_cutoff #if you would like to apply the AMFE cutoff, comment this line.
  new_file_w_AMFE_length_cutoff = new_file_w_AMFE_cutoff[which(new_file_w_AMFE_cutoff$length<=max_length),]
  ranked_dsRNAs = new_file_w_AMFE_length_cutoff[order(new_file_w_AMFE_length_cutoff$AMFE),]
  ranked_dsRNAs$rank = c(1:nrow(ranked_dsRNAs))

    
  write.table(ranked_dsRNAs, file=out_bed_fn, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

}
dev.off()


