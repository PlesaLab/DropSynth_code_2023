#!/bin/bash
#protocol.sh
# redirect stdout/stderr to a file
exec > >(tee -i log.txt)
exec 2>&1

echo "Hi, $USER. This is the main DropSynth oligo generation script v3."
echo "Please read the README.md before running"

echo "Screening assembly and amplification primers."
#screen primers
#python primer_screen.py


############################
## For 4 oligo HAMP end:
echo "Generating genes and splitting into oligo payload. He4"
python split_genes_for_oligos.py -c HK_He_4_oligos.conf > logs/HK_He_4_oligos.log

############################
## For 5 oligo HAMP end:
cat fasta_input/HK_2TM_1EC_HAMP_5oligo_wpos_HAMPend.fasta db_oligo/HK4_He_out.extra_split.proteins > fasta_input/HK_2TM_1EC_HAMP_5oligo_wpos_HAMPend_w_extra.fasta
echo "Generating genes and splitting into oligo payload. He5"
python split_genes_for_oligos.py -c HK_He_5_oligos.conf > logs/HK_He_5_oligos.log

#python divide_into_libs.py -c HK_He_5_oligos.conf > logs/HK_He_5_oligos_divide.log

###################
#./sort_libs.sh
echo "Adding amplification primers (BEFORE BLAT), microbead barcodes, and nicking sites. He4"
python oligobuffergen.py -c HK_He_4_oligos.conf > logs/HK_He_4_oligos_buffergen.log
echo "Adding amplification primers (BEFORE BLAT), microbead barcodes, and nicking sites. He5"
python oligobuffergen.py -c HK_He_5_oligos.conf > logs/HK_He_5_oligos_buffergen.log
#./sort_libs.sh

###################

echo "Screening amp and asm primers for off-target hybridization."
# screen primers for off-target hybridization
./BLAT_primers.sh
python BLAT_hits_parse.py

echo "NOTE: Depending on primer hits above you may have to rerun oligobuffergen for all libs."
#depending on hits you may have to rerun oligobuffergen

echo "Final oligo check of design rules."
#check as many design rules as possible
#do not resuse code in this script to avoid propagating errors
python final_oligo_check.py

echo "Making CSV file with oligos."
#make a csv file with all oligo seqs to order microarray
python FASTA_to_csv_for_Chip.py

echo "Checking oligos for homopolymers."
#check homopolymers
python homopolymer_stats.py

echo "Done."
