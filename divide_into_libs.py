#!/usr/bin/python

#this file will take output from split_genes_for_oligos.py
#and split them into separate libraries of a chosen size

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import math
import sys, getopt
import configparser
import csv

##################################
#FUNCTIONS:
def getBuildOligos(filename,num_oligos):
    handle = open(filename,'r')
    buildoligos = []
    countOligos = 0
    localCount = num_oligos
    for line in handle:
        if line.startswith('>'):
            #make sure the previous sequence had the correct number of oligos
            if localCount != num_oligos:
                print("Wrong num oligos. Before: "+str(buildoligos[-1]))
                print("Total oligos: " +str(localCount)+". Expected: "+str(num_oligos))
            buildoligos.append([line.strip()])
            localCount = 0
        elif line != "":
            buildoligos[-1].append(Seq(line.strip()))
            countOligos += 1
            localCount += 1
        
    handle.close()
    return buildoligos

def getBuildOligos_separateReq(filename,num_oligos,req_headers_list):
    handle = open(filename,'r')
    buildoligos = []
    buildoligos_req = []
    countOligos = 0
    localCount = num_oligos
    req_flag = 0
    for line in handle:
        if line.startswith('>'):
            #make sure the previous sequence had the correct number of oligos
            if localCount != num_oligos:
                print("Wrong num oligos. Before: "+str(buildoligos[-1]))
                print("Total oligos: " +str(localCount)+". Expected: "+str(num_oligos))
            local_id = line[1:-1]
            if local_id in req_headers_list:
                buildoligos_req.append([line.strip()])
                req_flag = 1
            else:
                buildoligos.append([line.strip()])
                req_flag = 0
            localCount = 0
        elif line != "":
            if req_flag == 1:
                buildoligos_req[-1].append(Seq(line.strip()))
            else:
                buildoligos[-1].append(Seq(line.strip()))
            countOligos += 1
            localCount += 1
        
    handle.close()
    return buildoligos, buildoligos_req

def obtainFastaSequences(filename):
    handle = open(filename)
    records = []
    for seqrecord in SeqIO.parse(handle, "fasta"):
        records.append(seqrecord)
        #print(seqrecord.id)
        #print(len(seqrecord))
    #print(len(records))
    return records


def obtainFastaSequences_separateReq(filename,req_headers_list):
    handle = open(filename)
    records = []
    records_req = []
    for seqrecord in SeqIO.parse(handle, "fasta"):
        if seqrecord.id in req_headers_list:
            records_req.append(seqrecord)
        else:
            records.append(seqrecord)
        #print(seqrecord.id)
        #print(len(seqrecord))
    #print(len(records))
    return records, records_req

def outputOligoLibrary_notFasta(constructs,filename,append_or_write,index):
    handle = open(filename,append_or_write)
    if index == 0:
        codon_str = "Codon1"
    else:
        codon_str = "Codon2"
    for construct in constructs:
        basename = construct[0]
        handle.write(basename+'\n')
        for index in range(1,len(construct)):
            handle.write(str(construct[index])+'\n')
    handle.close()

def splitLibs(filename, location, split_lib_size, num_libs, num_oligos, output_location, rectype, codon2on, seqs_to_require_file):
    
    #sequences which must be included in each split
    require_filetype = seqs_to_require_file.split(".")[-1]
    req_headers_list = []
    if require_filetype == "fasta":
        seqs_to_require = obtainFastaSequences(seqs_to_require_file)
        for ihd in seqs_to_require:
            req_headers_list.append(ihd.id)

    elif require_filetype == "csv":
        with open(seqs_to_require_file, newline='') as f:
            reader = csv.reader(f)
            req_headers_list = list(reader)

    else:
        print("Something went wrong with the seqs_to_require_file")
    

    full_file_name = location+"/"+filename
    buildoligos, buildoligos_req = getBuildOligos_separateReq(full_file_name+".oligos",num_oligos,req_headers_list)
    total_rec_length = len(buildoligos) #req seq removed from count
    #the total number of records which must be in each library
    print("Found "+str(len(buildoligos_req))+" required buildoligos")

    num_seq_required = len(req_headers_list)

    calc_lib_optimal_size = (total_rec_length)/num_libs+num_seq_required

    if (calc_lib_optimal_size != split_lib_size):
        print("Warning Calculated Optimal Library Size (",str(calc_lib_optimal_size),") not equal to Lib Size in Config File!")

    #number of records remaining in the last lib
    rec_remainder_counter = (total_rec_length+num_seq_required*num_libs) % split_lib_size

    #compensante for required seqs
    max_num_libs_counter = int(math.ceil((total_rec_length-num_seq_required+num_seq_required*num_libs)/split_lib_size))
    if rec_remainder_counter != 0:
        print("total number of records: "+str(total_rec_length)+" (not including req. seqs.), to make "+ str(max_num_libs_counter-1)+" full libs of "+str(split_lib_size)+" recs and one partial lib of "+str(rec_remainder_counter)+" records.")
    else:
        print("total number of records: "+str(total_rec_length)+" (not including req. seqs.), to make "+ str(max_num_libs_counter)+" full libs")
    print(str(num_seq_required)+" sequences are required in each sublibrary")

    #this help indexing in last lib
    if rec_remainder_counter == 0: #perfect fit
        rec_remainder_counter = split_lib_size

    #obtainFastaSequences(filename)
    gene_nt_filename = location+'/'+filename+".full_wRE_wPrim.genes"
    gene_nt, gene_nt_req = obtainFastaSequences_separateReq(gene_nt_filename,req_headers_list)

    gene_nt_no_primer_or_RE_filename = location+'/'+filename+".full_nRE_nPrim.genes"
    gene_nt_no_primer_or_RE, gene_nt_no_primer_or_RE_req = obtainFastaSequences_separateReq(gene_nt_no_primer_or_RE_filename,req_headers_list)
    gene_nt_no_primer_filename = location+'/'+filename+".full_wRE_noPrim.genes"
    gene_nt_no_primer, gene_nt_no_primer_req = obtainFastaSequences_separateReq(gene_nt_no_primer_filename,req_headers_list)
    if rectype == 'AA':
        protein_aa_filename = location+'/'+filename+".proteins"
        protein_aa, protein_aa_req = obtainFastaSequences_separateReq(protein_aa_filename,req_headers_list)
    if rectype == 'DNA':
        dna_filename = location+'/'+filename+".dna"
        dna_recs, dna_recs_req = obtainFastaSequences_separateReq(dna_filename,req_headers_list)

    if codon2on == "yes":
        buildoligos2, buildoligos2_req = getBuildOligos_separateReq(full_file_name+".codon2.oligos",num_oligos,req_headers_list)

        #obtainFastaSequences(filename)
        gene_nt_filename2 = location+'/'+filename+".codon2.full_wRE_wPrim.genes"
        gene_nt2, gene_nt2_req = obtainFastaSequences_separateReq(gene_nt_filename2,req_headers_list)

        gene_nt_no_primer_or_RE_filename2 = location+'/'+filename+".codon2.full_nRE_nPrim.genes"
        gene_nt_no_primer_or_RE2, gene_nt_no_primer_or_RE2_req = obtainFastaSequences_separateReq(gene_nt_no_primer_or_RE_filename2,req_headers_list)
        gene_nt_no_primer_filename2 = location+'/'+filename+".codon2.full_wRE_noPrim.genes"
        gene_nt_no_primer2, gene_nt_no_primer2_req = obtainFastaSequences_separateReq(gene_nt_no_primer_filename2,req_headers_list)
        if rectype == 'AA':
            protein_aa_filename2 = location+'/'+filename+".codon2.proteins"
            protein_aa2, protein_aa2_req = obtainFastaSequences_separateReq(protein_aa_filename2,req_headers_list)
        if rectype == 'DNA':
            dna_filename2 = location+'/'+filename+".codon2.dna"
            dna_recs2, dna_recs2_req = obtainFastaSequences_separateReq(dna_filename2,req_headers_list)

    libcounter = 1
    if num_libs <= max_num_libs_counter:#make sure you dont go over
        while libcounter < num_libs+1:
            #always add controls first, then remaining records
            ###############################################
            filename_split_lib = output_location+"/"+filename+"_split_Lib"+str(libcounter)+".oligos"
            #first is 0 to 383 - Num_req_recs
            #second is 384- Num_req_recs to 767- Num_req_recs
            if num_seq_required > 0:
                
                if libcounter == max_num_libs_counter: #last lib may be partial
                    start_index = (split_lib_size-num_seq_required)*(libcounter-1)
                    end_index = (split_lib_size-num_seq_required)*(libcounter-1)+rec_remainder_counter
                    print("split_lib_size: "+str(split_lib_size)
                        +" num_seq_required: "+str(num_seq_required)
                        +" libcounter: "+str(libcounter)
                        +" rec_remainder_counter: "+str(rec_remainder_counter))
                    print("last lib wreq start: "+str(start_index)+" end: "+str(end_index))
                elif libcounter == 1:
                    start_index = 0
                    end_index = (split_lib_size-num_seq_required)*libcounter
                    print("first lib wreq start: "+str(start_index)+" end: "+str(end_index))
                else:
                    start_index = (split_lib_size-num_seq_required)*(libcounter-1)
                    end_index = (split_lib_size-num_seq_required)*libcounter
                    print("middle lib wreq start: "+str(start_index)+" end: "+str(end_index))
            else:

                if libcounter == max_num_libs_counter: #last lib may be partial
                    start_index = split_lib_size*(libcounter-1)
                    end_index = split_lib_size*(libcounter-1)+rec_remainder_counter
                    print("last lib noreq start: "+str(start_index)+" end: "+str(end_index))
                elif libcounter == 1:
                    start_index = 0
                    end_index = split_lib_size*libcounter
                    print("first lib noreq start: "+str(start_index)+" end: "+str(end_index))
                else:
                    start_index = split_lib_size*(libcounter-1)
                    end_index = split_lib_size*libcounter
                    print("middle lib noreq start: "+str(start_index)+" end: "+str(end_index))

            #concatenate required recs onto this portion of the lib:
            temp_build_oligos = buildoligos_req + buildoligos[start_index:end_index]
            outputOligoLibrary_notFasta(temp_build_oligos,filename_split_lib,'w',0)
            print("Wrote "+str(len(temp_build_oligos))+" records to "+filename_split_lib)
            
            ###############################################
            temp_filename = output_location+"/"+filename+"_split_Lib"+str(libcounter)+".full_wRE_wPrim.genes"
            gene_nt_fileout = open(temp_filename,'w')

            temp_gene_nt = gene_nt_req+gene_nt[start_index:end_index]
            
            SeqIO.write(temp_gene_nt, gene_nt_fileout, "fasta")
            gene_nt_fileout.close()
            print("Wrote "+str(len(temp_gene_nt))+" records to "+temp_filename)

            ###############################################
            temp_filename = output_location+"/"+filename+"_split_Lib"+str(libcounter)+".full_nRE_nPrim.genes"
            gene_nt_no_primer_or_RE_fileout = open(temp_filename,'w')

            temp_gene_nt_no_primer_or_RE = gene_nt_no_primer_or_RE_req + gene_nt_no_primer_or_RE[start_index:end_index]

            SeqIO.write(temp_gene_nt_no_primer_or_RE, gene_nt_no_primer_or_RE_fileout, "fasta")
            gene_nt_no_primer_or_RE_fileout.close()
            print("Wrote "+str(len(temp_gene_nt_no_primer_or_RE))+" records to "+temp_filename)

            ###############################################
            temp_filename = output_location+"/"+filename+"_split_Lib"+str(libcounter)+".full_wRE_noPrim.genes"
            gene_nt_no_primer_fileout = open(temp_filename,'w')

            temp_gene_nt_no_primer = gene_nt_no_primer_req + gene_nt_no_primer[start_index:end_index]

            SeqIO.write(temp_gene_nt_no_primer, gene_nt_no_primer_fileout, "fasta")
            gene_nt_no_primer_fileout.close()
            print("Wrote "+str(len(temp_gene_nt_no_primer))+" records to "+temp_filename)

            ###############################################
            if rectype == 'AA':
                temp_filename = output_location+"/"+filename+"_split_Lib"+str(libcounter)+".proteins"
                protein_aa_fileout = open(temp_filename,'w')

                temp_protein_aa = protein_aa_req + protein_aa[start_index:end_index]

                SeqIO.write(temp_protein_aa, protein_aa_fileout, "fasta")
                protein_aa_fileout.close()
                print("Wrote "+str(len(temp_protein_aa))+" records to "+temp_filename)

            ###############################################
            if rectype == 'DNA':
                temp_filename = output_location+"/"+filename+"_split_Lib"+str(libcounter)+".dna"
                dna_fileout = open(temp_filename,'w')
                temp_dna_recs = dna_recs_req + dna_recs[start_index:end_index]
                SeqIO.write(temp_dna_recs, dna_fileout, "fasta")
                dna_fileout.close()
                print("Wrote "+str(len(temp_dna_recs))+" records to "+temp_filename)

            ###############################################
            if codon2on == "yes":

                ###############################################
                filename_split_lib2 = output_location+"/"+filename+"_split_Lib"+str(libcounter)+".codon2.oligos"
                #first is 0 to 383
                #second is 384 to 767
                temp_build_oligos2 = buildoligos2_req+buildoligos2[start_index:end_index]

                outputOligoLibrary_notFasta(temp_build_oligos2,filename_split_lib2,'w',0)
                print("Wrote "+str(len(temp_build_oligos2))+" records to "+filename_split_lib2)

                ###############################################
                temp_filename = output_location+"/"+filename+"_split_Lib"+str(libcounter)+".codon2.full_wRE_wPrim.genes"
                gene_nt_fileout2 = open(temp_filename,'w')

                temp_gene_nt2 = gene_nt2_req + gene_nt2[start_index:end_index]

                SeqIO.write(temp_gene_nt2, gene_nt_fileout2, "fasta")
                gene_nt_fileout2.close()
                print("Wrote "+str(len(temp_gene_nt2))+" records to "+temp_filename)

                ###############################################
                temp_filename = output_location+"/"+filename+"_split_Lib"+str(libcounter)+".codon2.full_nRE_nPrim.genes"
                gene_nt_no_primer_or_RE_fileout2 = open(temp_filename,'w')

                temp_gene_nt_no_primer_or_RE2 = gene_nt_no_primer_or_RE2_req + gene_nt_no_primer_or_RE2[start_index:end_index]

                SeqIO.write(temp_gene_nt_no_primer_or_RE2, gene_nt_no_primer_or_RE_fileout2, "fasta")
                gene_nt_no_primer_or_RE_fileout2.close()
                print("Wrote "+str(len(temp_gene_nt_no_primer_or_RE2))+" records to "+temp_filename)

                ###############################################
                temp_filename = output_location+"/"+filename+"_split_Lib"+str(libcounter)+".codon2.full_wRE_noPrim.genes"
                gene_nt_no_primer_fileout2 = open(temp_filename,'w')

                temp_gene_nt_no_primer2 = gene_nt_no_primer2_req + gene_nt_no_primer2[start_index:end_index]

                SeqIO.write(temp_gene_nt_no_primer2, gene_nt_no_primer_fileout2, "fasta")
                gene_nt_no_primer_fileout2.close()
                print("Wrote "+str(len(temp_gene_nt_no_primer2))+" records to "+temp_filename)

                ###############################################
                if rectype == 'AA':
                    temp_filename = output_location+"/"+filename+"_split_Lib"+str(libcounter)+".codon2.proteins"
                    protein_aa_fileout2 = open(temp_filename,'w')

                    temp_protein_aa2 = protein_aa2_req + protein_aa2[start_index:end_index]

                    SeqIO.write(temp_protein_aa2, protein_aa_fileout2, "fasta")
                    protein_aa_fileout2.close()
                    print("Wrote "+str(len(temp_protein_aa2))+" records to "+temp_filename)

                ###############################################
                if rectype == 'DNA':
                    temp_filename = output_location+"/"+filename+"_split_Lib"+str(libcounter)+".codon2.dna"
                    dna_fileout2 = open(temp_filename,'w')

                    temp_dna_recs2 = dna_recs2_req + dna_recs2[start_index:end_index]

                    SeqIO.write(temp_gene_nt2, dna_fileout2, "fasta")
                    dna_fileout2.close()
                    print("Wrote "+str(len(temp_gene_nt2))+" records to "+temp_filename)

                ###############################################


            libcounter += 1

def main(argv):
   configfile = ''
   try:
      opts, args = getopt.getopt(argv,"hc:",["configfile="])
   except getopt.GetoptError:
      print('divide_into_libs.py -c <configfile>')
      sys.exit(2)
   for opt, arg in opts:
      if opt == '-h':
         print('divide_into_libs.py -c <configfile>')
         sys.exit()
      elif opt in ("-c", "--config"):
         inputfile = arg
   print('Config file is: ', inputfile)
   return inputfile

if __name__ == '__main__':

    config_path = main(sys.argv[1:])
    config = configparser.ConfigParser()
    config.read_file(open("lib_config/"+config_path))

    #the base name of this library
    filename = config.get('lib_specific', 'output_file_names').replace('.fasta','')

    #where did split_genes_for_oligos.py put its output files?
    location = config.get('lib_specific', 'output_folder')

    #how many genes per lib 384, 1536, etc...
    lib_size = int(config.get('Lib_splitting', 'div_lib_size'))

    number_of_libs = int(config.get('Lib_splitting', 'number_of_libs'))

    # number of oligos for each split gene
    num_oligos_input = int(config.get('lib_specific', 'num_oligos'))

    #are these proteins or DNA sequences
    rectype_in = config.get('lib_specific', 'input_type')

    #two different codon optimizations for each protein?
    if config.get('second_lib', 'gen_two_codon') == "True":
        codon2on = "yes"
    else:
        codon2on = "no"

    output_location = config.get('Lib_splitting', 'output_location')


    seqs_to_require_file = config.get('Lib_splitting', 'seqs_to_require_file')

    needs_split = config.get('Lib_splitting', 'needs_split') == "True"


    #filename, location, size, number, num oligos, rectype, codon2on
    #input_array = [["HPseqs_3_oligos", "db_oligo", "384" ,"5", "3", "AA", "yes"],
    #                ["HPseqs_4_oligos", "db_oligo", "384", "2", "4", "AA", "yes"],
    #                ["HPseqs_5_oligos", "db_oligo", "384", "2", "5", "AA", "yes"],
    #                ["HPseqs_6_oligos", "db_oligo", "384", "2", "6", "AA", "yes"],
    #                ["HPseqs_7_oligos", "db_oligo", "384", "2", "7", "AA", "yes"],
    #                ["HPseqs_8_oligos", "db_oligo", "384", "2", "8", "AA", "yes"]]

    if needs_split:
        input_array = [[filename, location, lib_size, number_of_libs, num_oligos_input, rectype_in, codon2on, seqs_to_require_file]]

        for i in input_array:
            print(i)
            splitLibs(i[0],i[1],int(i[2]), int(i[3]), int(i[4]),output_location,i[5],i[6],i[7])
    else:
        print("This library does not need to be split! Move along...")

