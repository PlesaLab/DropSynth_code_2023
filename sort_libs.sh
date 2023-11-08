#!/bin/bash
#protocol.sh
# redirect stdout/stderr to a file
exec > >(tee -i log.txt)
exec 2>&1

echo "Starting sort."

input_folder=db_libs
output_folder=db_oligo_sorted
output_lib=db_libs_sorted


declare -a libnum=("3" "4")
declare -a libname=("01" "02" "03" "04" "05" "06" "07" "08")
COUNTER=0
for i in "${libnum[@]}"
do
	cp "$input_folder"/HPseqs_"$i"_oligos_split_Lib1.oligos "$output_folder"/Lib"${libname[$COUNTER]}"_"$i"oli_cod1.oligos
	cp "$input_folder"/HPseqs_"$i"_oligos_split_Lib1.proteins "$output_folder"/Lib"${libname[$COUNTER]}"_"$i"oli_cod1.proteins
	cp "$input_folder"/HPseqs_"$i"_oligos_split_Lib1.proteins "$output_lib"/Lib"${libname[$COUNTER]}"_"$i"oli_cod1.proteins
	
	cp "$input_folder"/HPseqs_"$i"_oligos_split_Lib1-finaloligos.fasta "$output_folder"/Lib"${libname[$COUNTER]}"_"$i"oli_cod1-finaloligos.fasta
	cp "$input_folder"/HPseqs_"$i"_oligos_split_Lib1-finaloligos.fasta "$output_lib"/Lib"${libname[$COUNTER]}"_"$i"oli_cod1-finaloligos.fasta
	
	cp "$input_folder"/HPseqs_"$i"_oligos_split_Lib1.full_wRE_wPrim.genes "$output_folder"/Lib"${libname[$COUNTER]}"_"$i"oli_cod1.full_wRE_wPrim.genes
	cp "$input_folder"/HPseqs_"$i"_oligos_split_Lib1.full_wRE_noPrim.genes "$output_folder"/Lib"${libname[$COUNTER]}"_"$i"oli_cod1.full_wRE_noPrim.genes
	cp "$input_folder"/HPseqs_"$i"_oligos_split_Lib1.full_nRE_nPrim.genes "$output_folder"/Lib"${libname[$COUNTER]}"_"$i"oli_cod1.full_nRE_nPrim.genes
	let COUNTER=COUNTER+1
	#################################
	cp "$input_folder"/HPseqs_"$i"_oligos_split_Lib1.codon2.oligos "$output_folder"/Lib"${libname[$COUNTER]}"_"$i"oli_cod2.oligos
	cp "$input_folder"/HPseqs_"$i"_oligos_split_Lib1.codon2.proteins "$output_folder"/Lib"${libname[$COUNTER]}"_"$i"oli_cod2.proteins
	cp "$input_folder"/HPseqs_"$i"_oligos_split_Lib1.codon2.proteins "$output_lib"/Lib"${libname[$COUNTER]}"_"$i"oli_cod2.proteins

	cp "$input_folder"/HPseqs_"$i"_oligos_split_Lib1-finaloligos.codon2.fasta "$output_folder"/Lib"${libname[$COUNTER]}"_"$i"oli_cod2-finaloligos.fasta
	cp "$input_folder"/HPseqs_"$i"_oligos_split_Lib1-finaloligos.codon2.fasta "$output_lib"/Lib"${libname[$COUNTER]}"_"$i"oli_cod2-finaloligos.fasta
	
	cp "$input_folder"/HPseqs_"$i"_oligos_split_Lib1.codon2.full_wRE_wPrim.genes "$output_folder"/Lib"${libname[$COUNTER]}"_"$i"oli_cod2.full_wRE_wPrim.genes
	cp "$input_folder"/HPseqs_"$i"_oligos_split_Lib1.codon2.full_wRE_noPrim.genes "$output_folder"/Lib"${libname[$COUNTER]}"_"$i"oli_cod2.full_wRE_noPrim.genes
	cp "$input_folder"/HPseqs_"$i"_oligos_split_Lib1.codon2.full_nRE_nPrim.genes "$output_folder"/Lib"${libname[$COUNTER]}"_"$i"oli_cod2.full_nRE_nPrim.genes
	let COUNTER=COUNTER+1
	#################################
	cp "$input_folder"/HPseqs_"$i"_oligos_split_Lib2.oligos "$output_folder"/Lib"${libname[$COUNTER]}"_"$i"oli_cod1.oligos
	cp "$input_folder"/HPseqs_"$i"_oligos_split_Lib2.proteins "$output_folder"/Lib"${libname[$COUNTER]}"_"$i"oli_cod1.proteins
	cp "$input_folder"/HPseqs_"$i"_oligos_split_Lib2.proteins "$output_lib"/Lib"${libname[$COUNTER]}"_"$i"oli_cod1.proteins

	cp "$input_folder"/HPseqs_"$i"_oligos_split_Lib2-finaloligos.fasta "$output_folder"/Lib"${libname[$COUNTER]}"_"$i"oli_cod1-finaloligos.fasta
	cp "$input_folder"/HPseqs_"$i"_oligos_split_Lib2-finaloligos.fasta "$output_lib"/Lib"${libname[$COUNTER]}"_"$i"oli_cod1-finaloligos.fasta
	
	cp "$input_folder"/HPseqs_"$i"_oligos_split_Lib2.full_wRE_wPrim.genes "$output_folder"/Lib"${libname[$COUNTER]}"_"$i"oli_cod1.full_wRE_wPrim.genes
	cp "$input_folder"/HPseqs_"$i"_oligos_split_Lib2.full_wRE_noPrim.genes "$output_folder"/Lib"${libname[$COUNTER]}"_"$i"oli_cod1.full_wRE_noPrim.genes
	cp "$input_folder"/HPseqs_"$i"_oligos_split_Lib2.full_nRE_nPrim.genes "$output_folder"/Lib"${libname[$COUNTER]}"_"$i"oli_cod1.full_nRE_nPrim.genes
	let COUNTER=COUNTER+1
	#################################
	cp "$input_folder"/HPseqs_"$i"_oligos_split_Lib2.codon2.oligos "$output_folder"/Lib"${libname[$COUNTER]}"_"$i"oli_cod2.oligos
	cp "$input_folder"/HPseqs_"$i"_oligos_split_Lib2.codon2.proteins "$output_folder"/Lib"${libname[$COUNTER]}"_"$i"oli_cod2.proteins
	cp "$input_folder"/HPseqs_"$i"_oligos_split_Lib2.codon2.proteins "$output_lib"/Lib"${libname[$COUNTER]}"_"$i"oli_cod2.proteins

	cp "$input_folder"/HPseqs_"$i"_oligos_split_Lib2-finaloligos.codon2.fasta "$output_folder"/Lib"${libname[$COUNTER]}"_"$i"oli_cod2-finaloligos.fasta
	cp "$input_folder"/HPseqs_"$i"_oligos_split_Lib2-finaloligos.codon2.fasta "$output_lib"/Lib"${libname[$COUNTER]}"_"$i"oli_cod2-finaloligos.fasta
	
	cp "$input_folder"/HPseqs_"$i"_oligos_split_Lib2.codon2.full_wRE_wPrim.genes "$output_folder"/Lib"${libname[$COUNTER]}"_"$i"oli_cod2.full_wRE_wPrim.genes
	cp "$input_folder"/HPseqs_"$i"_oligos_split_Lib2.codon2.full_wRE_noPrim.genes "$output_folder"/Lib"${libname[$COUNTER]}"_"$i"oli_cod2.full_wRE_noPrim.genes
	cp "$input_folder"/HPseqs_"$i"_oligos_split_Lib2.codon2.full_nRE_nPrim.genes "$output_folder"/Lib"${libname[$COUNTER]}"_"$i"oli_cod2.full_nRE_nPrim.genes
	let COUNTER=COUNTER+1
done

echo "Done."
