library(ggplot2)
library(dplyr)
library(magrittr)
library(ggvenn)
library(Biostrings)
library(seqinr)

setwd("/Users/calin/MyDrive/faculty/projects/TCS/sequences/theBIGone/oligo_design/fasta_input/")

TM2end4o = readAAStringSet("HK_2TM_1EC_HAMP_4oligo_wpos_TM2end.fasta")
TM2end4o <- data.frame(seq_name=names(TM2end4o), sequence=paste(TM2end4o))

TM2start4o = readAAStringSet("HK_2TM_1EC_HAMP_4oligo_wpos_TM2start.fasta")
TM2start4o <- data.frame(seq_name=names(TM2start4o), sequence=paste(TM2start4o))

TM2end5o = readAAStringSet("HK_2TM_1EC_HAMP_5oligo_wpos_TM2end.fasta")
TM2end5o <- data.frame(seq_name=names(TM2end5o), sequence=paste(TM2end5o))

TM2start5o = readAAStringSet("HK_2TM_1EC_HAMP_5oligo_wpos_TM2start.fasta")
TM2start5o <- data.frame(seq_name=names(TM2start5o), sequence=paste(TM2start5o))

TM2end = rbind(TM2end5o, TM2end4o) %>%
  mutate(LengthAASeq=nchar(sequence),
         Length.nt = LengthAASeq*3)

TM2start = rbind(TM2start5o, TM2start4o) %>%
  mutate(LengthAASeq=nchar(sequence),
         Length.nt = LengthAASeq*3)

#Sanity check to make sure no bad AA
TM2end_bad_AA <- TM2end %>%
  filter(grepl("B", sequence, fixed=TRUE) | grepl("J", sequence, fixed=TRUE) | grepl("O", sequence, fixed=TRUE) | grepl("U", sequence, fixed=TRUE) | grepl("X", sequence, fixed=TRUE) | grepl("Z", sequence, fixed=TRUE))

TM2start_bad_AA <- TM2start %>%
  filter(grepl("B", sequence, fixed=TRUE) | grepl("J", sequence, fixed=TRUE) | grepl("O", sequence, fixed=TRUE) | grepl("U", sequence, fixed=TRUE) | grepl("X", sequence, fixed=TRUE) | grepl("Z", sequence, fixed=TRUE))

rm(TM2end_bad_AA,TM2start_bad_AA)


#Sanity check to make sure starts with Met
TM2end_bad_M <- TM2end %>%
  filter(!startsWith(sequence, 'M'))
TM2start_bad_M <- TM2start %>%
  filter(!startsWith(sequence, 'M'))

#fix M for A0A6V1PG63

TM2start$sequence[TM2start$seq_name=="A0A6V1PG63"]<-paste("M",TM2start$sequence[TM2start$seq_name=="A0A6V1PG63"],sep="")
TM2end$sequence[TM2end$seq_name=="A0A6V1PG63"]<-paste("M",TM2end$sequence[TM2end$seq_name=="A0A6V1PG63"],sep="")

TM2end_bad_M <- TM2end %>%
  filter(!startsWith(sequence, 'M'))
TM2start_bad_M <- TM2start %>%
  filter(!startsWith(sequence, 'M'))

rm(TM2end_bad_M, TM2start_bad_M)

#300mers
dropsynth_bins = c(258,391,601,811,1021,1231,1441,1651)
dropsynth_bins_lower = c(1,258,391,601,811,1021,1231,1441)
############################################################
# length filtered after here
###########################################################
TM2endds <- TM2end %>%
  select(-LengthAASeq) %>%
  filter(Length.nt < max(dropsynth_bins)) %>%
  mutate(DSOligos = 0) %>%
    mutate(DSOligos = 
             case_when(
      Length.nt >= dropsynth_bins_lower[1] & Length.nt < dropsynth_bins[1] ~ 1,
      Length.nt >= dropsynth_bins_lower[2] & Length.nt < dropsynth_bins[2] ~ 2,
      Length.nt >= dropsynth_bins_lower[3] & Length.nt < dropsynth_bins[3] ~ 3,
      Length.nt >= dropsynth_bins_lower[4] & Length.nt < dropsynth_bins[4] ~ 4,
      Length.nt >= dropsynth_bins_lower[5] & Length.nt < dropsynth_bins[5] ~ 5,
      Length.nt >= dropsynth_bins_lower[6] & Length.nt < dropsynth_bins[6] ~ 6,
      Length.nt >= dropsynth_bins_lower[7] & Length.nt < dropsynth_bins[7] ~ 7,
      Length.nt >= dropsynth_bins_lower[8] & Length.nt < dropsynth_bins[8] ~ 8)) %>%
  mutate(NumOligos = as.character(DSOligos))


#how many proteins in each bind size, no filtering
TM2endds %>% group_by(DSOligos) %>%
  summarize(n = n()) 

names(TM2endds)
##################################
TM2startds <- TM2start %>%
  select(-LengthAASeq) %>%
  filter(Length.nt < max(dropsynth_bins)) %>%
  mutate(DSOligos = 0) %>%
  mutate(DSOligos = 
           case_when(
             Length.nt >= dropsynth_bins_lower[1] & Length.nt < dropsynth_bins[1] ~ 1,
             Length.nt >= dropsynth_bins_lower[2] & Length.nt < dropsynth_bins[2] ~ 2,
             Length.nt >= dropsynth_bins_lower[3] & Length.nt < dropsynth_bins[3] ~ 3,
             Length.nt >= dropsynth_bins_lower[4] & Length.nt < dropsynth_bins[4] ~ 4,
             Length.nt >= dropsynth_bins_lower[5] & Length.nt < dropsynth_bins[5] ~ 5,
             Length.nt >= dropsynth_bins_lower[6] & Length.nt < dropsynth_bins[6] ~ 6,
             Length.nt >= dropsynth_bins_lower[7] & Length.nt < dropsynth_bins[7] ~ 7,
             Length.nt >= dropsynth_bins_lower[8] & Length.nt < dropsynth_bins[8] ~ 8)) %>%
  mutate(NumOligos = as.character(DSOligos))


#how many proteins in each bind size, no filtering
TM2startds %>% group_by(DSOligos) %>%
  summarize(n = n()) 

names(TM2startds)

#TM2end
# 3 alone
TM2endds_temp <- TM2endds %>%
  filter(DSOligos==3)
fileout_name = paste("sort_seqsR_HK_2TM_1EC_HAMP_3only_oligo_wpos_TM2end.fasta",sep="")
write.fasta(sequences=as.list(TM2endds_temp$sequence), 
            names=as.list(TM2endds_temp$seq_name),
            file.out=fileout_name, 
            open = "w", nbchar = 60, as.string = TRUE)

# 4 alone
TM2endds_temp <- TM2endds %>%
  filter(DSOligos==4)
fileout_name = paste("sort_seqsR_HK_2TM_1EC_HAMP_4only_oligo_wpos_TM2end.fasta",sep="")
write.fasta(sequences=as.list(TM2endds_temp$sequence), 
            names=as.list(TM2endds_temp$seq_name),
            file.out=fileout_name, 
            open = "w", nbchar = 60, as.string = TRUE)

#have 5 oligo alone
TM2endds_temp <- TM2endds %>%
  filter(DSOligos==5)
fileout_name = paste("sort_seqsR_HK_2TM_1EC_HAMP_5only_oligo_wpos_TM2end.fasta",sep="")
write.fasta(sequences=as.list(TM2endds_temp$sequence), 
            names=as.list(TM2endds_temp$seq_name),
            file.out=fileout_name, 
            open = "w", nbchar = 60, as.string = TRUE)

#TM2start
#
#
#do 3 alone
TM2startds_temp <- TM2startds %>%
  filter(DSOligos==3)
fileout_name = paste("sort_seqsR_HK_2TM_1EC_HAMP_3only_oligo_wpos_TM2start.fasta",sep="")
write.fasta(sequences=as.list(TM2startds_temp$sequence), 
            names=as.list(TM2startds_temp$seq_name),
            file.out=fileout_name, 
            open = "w", nbchar = 60, as.string = TRUE)

#do 4 alone
TM2startds_temp <- TM2startds %>%
  filter(DSOligos==4)
fileout_name = paste("sort_seqsR_HK_2TM_1EC_HAMP_4only_oligo_wpos_TM2start.fasta",sep="")
write.fasta(sequences=as.list(TM2startds_temp$sequence), 
            names=as.list(TM2startds_temp$seq_name),
            file.out=fileout_name, 
            open = "w", nbchar = 60, as.string = TRUE)
