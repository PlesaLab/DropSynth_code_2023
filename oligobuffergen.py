#!/usr/bin/python

#this script does the following:
#for each library
#get oligos
#add typeIIs cutters (BtsI)
#check that no new restriction sites added
#pad the oligos to same length
#add barcode with nicking sites (Nt.BspQI)
#check that no new restriction sites added
#add amplification primers
#check that no new restriction sites added
#output

#python 3.7
#updated 09.27.2021 Calin Plesa

#to run use
#python oligobuffergen.py -c <configfile>
#python oligobuffergen.py -c HP_Lib_1.conf

#to profile
#python -m cProfile -o profile_08.txt oligobuffergen.py -c <configfile>
#echo 'stats' | python3 -m pstats profile_08.txt > profile_08_humanreadable.txt


from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Restriction import *
#import pylab #for the plots
import sys, getopt
import configparser
import random
import re
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

def obtainFastaSequences(filename):
    handle = open(filename)
    records = []
    for seqrecord in SeqIO.parse(handle, "fasta"):
        records.append(seqrecord)
        #print(seqrecord.id)
        #print(len(seqrecord))
    #print(len(records))
    return records


def occurrencesoverlap(string, sub):
    count = start = 0
    while True:
        start = string.find(sub, start) + 1
        if start > 0:
            count+=1
        else:
            return count

def checkRestrictionSites(seqrecord, REbatch, verbosemode=True):
    temp_str = str(seqrecord.seq)
    fail_flag = 0

    key_array = []
    key_count_allowed = dict()
    key_overlap = dict()
    for key in REbatch:
        siteseq=REbatch[key][1]
        key_array.append(siteseq) #add RE site to list
        key_overlap[siteseq] = REbatch[key][2]
        key_count_allowed[siteseq] = REbatch[key][3]
        #print("key "+key)
        #print("overlap "+REbatch[key][2])
        #print("count "+str(REbatch[key][3]))


    for key in key_array:
        #loop over all matches, if any
        #this site is not allowed at all:
        if key_count_allowed[key] == 0:
            for m in re.finditer(key, str(seqrecord.seq)):
                fail_flag = 1
        #these sites are allowed once
        elif key_count_allowed[key] > 0:
            if key_overlap[key] == "Yes":
                #count occurances with overlap
                count_sites = occurrencesoverlap(str(seqrecord.seq),key)
                if count_sites != key_count_allowed[key]:
                    if verbosemode:
                        print("Found "+str(count_sites)+" sites of "+key+" instead of "+str(key_count_allowed[key]))
                        print(temp_str)
                    fail_flag = 1
            elif key_overlap[key] == "No":
                count_sites = str(seqrecord.seq).count(key)
                #print(key)
                #print(count_sites)
                #print(key_count_allowed[key])
                if count_sites != key_count_allowed[key]:
                    if verbosemode:
                        print("Found "+str(count_sites)+" sites of "+key+" instead of "+str(key_count_allowed[key]))
                        print(temp_str)
                    fail_flag = 1

    return fail_flag

# Reversing a list using reversed()
def Reverse(lst):
    return [ele for ele in reversed(lst)]

def addBtsI(constructs, REbatch):
    #add BtsI sites and check that only two exist

    num_oligosin = len(constructs)
    fail_flag = dict()
    for i in range(0,num_oligosin):
        fail_flag[i] = 0

    newconstructs = []
    for construct in constructs:
        newconstruct = []
        newconstruct.append(construct[0])
        for i in range(1,len(construct)):#length includes header
            newoligo = Seq("GCAGTG") + construct[i] + Seq("CACTGC")
            
            fail_flag[i] = checkRestrictionSites(SeqRecord(newoligo, id=""), REbatch)
            #rb = RestrictionBatch([BtsI, BspQI, BsaI])
            #seqsearch = rb.search(newoligo)
            if fail_flag[i] == 1:
                if i==1:
                    #first oligo
                    print("\naddBtsI: Bad number of restriction sites in first oligo")
                    print(construct[0] + '\t' + str(i))
                    print(newoligo)
                        
                elif i==(len(construct)-1):
                    #last oligo
                    print("\naddBtsI: Bad number of restriction sites in last oligo")
                    print(construct[0] + '\t' + str(i))
                    print(newoligo)

                elif i==(len(construct)-2):
                    #second to last oligo
                    print("\naddBtsI: Bad number of restriction sites in second to last oligo oligo")
                    print(construct[0] + '\t' + str(i))
                    print(newoligo)

                else:
                    #middle oligos
                    print("\naddBtsI: Bad number of restriction sites in middle oligo")
                    print(construct[0] + '\t' + str(i))
                    print(newoligo)

            newconstruct.append(newoligo)
        newconstructs.append(newconstruct)
    return newconstructs

def plotOligo_length(constructs,filename):
    #plot histogram of lengths
    construct_lengths = []
    for construct in constructs:
        for i in range(1,len(construct)):
            construct_lengths.append(len(construct[i]))
    generate_histogram = True
    if generate_histogram:
        data = construct_lengths
        print("min Lev dist: " +str(min(data)) + " max Lev dist: " +str(max(data)))
        pylab.hist(data, bins=(max(data)-min(data)))
        pylab.title("%i oligos length distribution\nfrom %i to %i" \
                    % (len(data),min(data),max(data)))
        pylab.xlabel("Oligo length (nt)")
        pylab.ylabel("Count")
        pylab.savefig(filename+'.png', bbox_inches='tight')
        pylab.savefig(filename+'.pdf', bbox_inches='tight')
        pylab.show()

#generate random DNA and limit max homopolymer length
#http://stackoverflow.com/questions/21205836/generating-random-sequences-of-dna
def randomDNA(length, REbatch):
    #print(REbatch)
    gsflag = 0
    while gsflag == 0:
        temp_str_nocontext = ''.join(random.choice('CGTA') for _ in range(length))
        #this goes between BspQI and BspQI
        #it will stay attached to the bead
        temp_str = str(Seq("CGAAGAGC"))+temp_str_nocontext+str(Seq("GCAGTG"))
        #make sure no homopolymers greater than 5 in length or cut sites
        #BtsI (GCAGTG - RC: CACTGC), BspQI (GCTCTTC - RC: GAAGAGC)
        if ((temp_str.find('AAAAA') == -1) and (temp_str.find('TTTTT') == -1) and
            (temp_str.find('GGGGG') == -1) and (temp_str.find('CCCCC') == -1)):
            #check others
            fail_flag = 0
            fail_flag = checkRestrictionSites(SeqRecord(Seq(temp_str),id=""), REbatch, False)
            #if (temp_str_nocontext.find('GAAGAGC') != -1):
            #print("found GAAGAGC")
            #print(temp_str)
            #print(temp_str_nocontext)
            #print(fail_flag)
            if fail_flag != 1:
                gsflag = 1
    return temp_str_nocontext

def padOligoCP(constructs, totallength, REbatch, REbatch_random):
    #add bases until oligo is correct length
    newconstructs = []
    for construct in constructs:
        newconstruct = []
        newconstruct.append(construct[0])
        fail_flag = dict()
        for i in range(1,len(construct)):
            #pad oligo random bases
            oligo_len = len(construct[i])
            length_to_add = totallength - oligo_len
            
            full_padding_seq = randomDNA(length_to_add, REbatch_random)
            newoligo =  Seq(full_padding_seq) + construct[i] # add padding between BspQI and BtsI
            #rb = RestrictionBatch([BtsI, BspQI, BsaI])#BseRI
            #seqsearch = rb.search(newoligo)
            fail_flag[i] = checkRestrictionSites(SeqRecord(newoligo, id=""), REbatch)
            #rb = RestrictionBatch([BtsI, BspQI, BsaI])
            #seqsearch = rb.search(newoligo)
            if fail_flag[i] == 1:
                if i==1:
                    #first oligo
                    print("\nPadOligo: Bad number of restriction sites in first oligo")
                    print("padding: "+full_padding_seq)
                    print(construct[0] + '\t' + str(i))
                    print(newoligo)
                        
                elif i==(len(construct)-1):
                    #last oligo
                    print("\nPadOligo: Bad number of restriction sites in last oligo")
                    print("padding: "+full_padding_seq)
                    print(construct[0] + '\t' + str(i))
                    print(newoligo)

                elif i==(len(construct)-2):
                    #second to last oligo
                    print("\nPadOligo: Bad number of restriction sites in second to last oligo oligo")
                    print(construct[0] + '\t' + str(i))
                    print("padding: "+full_padding_seq)
                    print(newoligo)

                else:
                    #middle oligos
                    print("\nPadOligo: Bad number of restriction sites in middle oligo")
                    print("padding: "+full_padding_seq)
                    print(construct[0] + '\t' + str(i))
                    print(newoligo)

            newconstruct.append(newoligo)
        newconstructs.append(newconstruct)
    return newconstructs

def addBarcodes(constructs, barcodes, REbatch):
    newconstructs = []
    for index in range(len(constructs)):
        newconstruct = []
        fail_flag = dict()
        #print(barcodes[index].id)
        #print(barcodes[index+1].id)
        newconstruct.append(constructs[index][0]+';'+barcodes[index].id)
        for i in range(1,len(constructs[index])):#includes header in length
            newoligo = Seq("GCTCTTCG") + barcodes[index].seq + Seq("CGAAGAGC") + constructs[index][i]
            #rb = RestrictionBatch([BtsI, BspQI, BsaI])
            #seqsearch = rb.search(newoligo)
            fail_flag[i] = checkRestrictionSites(SeqRecord(newoligo, id=""), REbatch)
            if fail_flag[i] == 1:
                if i==1:
                    #first oligo
                    print("\naddBarcodes: Bad number of restriction sites in first oligo")
                    print(constructs[index][0] + '\t' + str(i))
                    print(newoligo)
                        
                elif i==(len(constructs[index])-1):
                    #last oligo
                    print("\naddBarcodes: Bad number of restriction sites in last oligo")
                    print(constructs[index][0] + '\t' + str(i))
                    print(newoligo)

                elif i==(len(constructs[index])-2):
                    #second to last oligo
                    print("\naddBarcodes: Bad number of restriction sites in second to last oligo oligo")
                    print(constructs[index][0] + '\t' + str(i))
                    print(newoligo)

                else:
                    #middle oligos
                    print("\naddBarcodes: Bad number of restriction sites in middle oligo")
                    print(constructs[index][0] + '\t' + str(i))
                    print(newoligo)

            newconstruct.append(newoligo)
        newconstructs.append(newconstruct)
    return newconstructs

def addAmpPrimers(constructs, fwdprim, revprim, REbatch):
    newconstructs = []
    for index in range(len(constructs)):
        newconstruct = []
        newconstruct.append(constructs[index][0])
        fail_flag = dict()
        for i in range(1,len(constructs[index])):#lgnth includes header, that's why no +1
            newoligo = fwdprim + constructs[index][i] + revprim.reverse_complement()
            fail_flag[i] = checkRestrictionSites(SeqRecord(newoligo, id=""), REbatch)
            if fail_flag[i] == 1:
                if i==1:
                    #first oligo
                    print("\naddAmpPrimers: Bad number of restriction sites in first oligo")
                    print(constructs[index][0] + '\t' + str(i))
                    print(newoligo)
                        
                elif i==(len(constructs[index])-1):
                    #last oligo
                    print("\naddAmpPrimers: Bad number of restriction sites in last oligo")
                    print(constructs[index][0] + '\t' + str(i))
                    print(newoligo)

                elif i==(len(constructs[index])-2):
                    #second to last oligo
                    print("\naddAmpPrimers: Bad number of restriction sites in second to last oligo oligo")
                    print(constructs[index][0] + '\t' + str(i))
                    print(newoligo)

                else:
                    #middle oligos
                    print("\naddAmpPrimers: Bad number of restriction sites in middle oligo")
                    print(constructs[index][0] + '\t' + str(i))
                    print(newoligo)

            newconstruct.append(newoligo)
        newconstructs.append(newconstruct)
    return newconstructs


def printOligoLibrary(constructs):
    for construct in constructs:
        print(construct[0])
        for oligo in construct[1:]:
            print('\t' + oligo + '\t' + str(len(oligo)))
            
def parseREstring_for_specific_function(REsettings, sitesheader, target_step):

    #output is a dictionary containing
    #name, site, number of times to check
    REout = dict()

    index_of_target = sitesheader.index(target_step)-1
    #print(target_step)
    #print(index_of_target)
    #print(sitesheader)
    #print(REsettings)

    for key in REsettings:
        #print(key)
        #if value is positive check it
        if REsettings[key][index_of_target] > -1:
            #name, site, 
            REout[key] = [key, REsettings[key][0], REsettings[key][1], int(REsettings[key][index_of_target])]
            #print(REout[key])

    return REout

def outputOligoLibrary(constructs,filename,append_or_write,codonver):
    handle = open(filename,append_or_write)
    print("Writing out files "+filename+" "+append_or_write+" "+str(codonver))
    if codonver == 1:
        codon_str = "Codon1"
    elif codonver == 2:
        codon_str = "Codon2"
    for construct in constructs:
        basename = construct[0] + ';' + codon_str + ';'
        for index in range(1,len(construct)):
            handle.write(basename+str(index)+'\n')
            handle.write(str(construct[index])+'\n')
    handle.close()

def main(argv):
   configfile = ''
   try:
      opts, args = getopt.getopt(argv,"hc:",["configfile="])
   except getopt.GetoptError:
      print('oligobuffergen.py -c <configfile>')
      sys.exit(2)
   for opt, arg in opts:
      if opt == '-h':
         print('oligobuffergen.py -c <configfile>')
         sys.exit()
      elif opt in ("-c", "--config"):
         inputfile = arg
   print('Config file is: ', inputfile)
   return inputfile

######################################################################
if __name__ == '__main__':

    config_path = main(sys.argv[1:])
    config = configparser.ConfigParser()
    config.read_file(open("lib_config/"+config_path))

    #was this split?
    needs_split = config.get('Lib_splitting', 'needs_split') == "True"

    output_location = config.get('Lib_splitting', 'output_location')

    if needs_split == True:
        
        output_folder = output_location

        number_of_libs = int(config.get('Lib_splitting', 'number_of_libs'))

        libarray= range(number_of_libs)

        filename = config.get('lib_specific', 'output_file_names').replace('.fasta','')

        filenames = []
        filenames_out = []

        for li in libarray:
            filename_split_lib = output_location+"/"+filename+"_split_Lib"+str(li+1)+".oligos"
            filenames.append(filename_split_lib)
            filenames_out.append(filename_split_lib.replace('.oligos','-finaloligos.fasta'))


    elif needs_split == False:
    
        output_folder = config.get('lib_specific', 'output_folder')

        #the source file with all sequences to be split
        filenames_temp = output_folder+"/"+config.get('DS_Lib_settings', 'oligo_input_file')
        filenames = [filenames_temp]
        #print(filenames)

        filenames_out_temp = filenames_temp.replace('.oligos','-finaloligos.fasta').replace('db_oligo',output_location)
        filenames_out = [filenames_out_temp]

    else:
        print("Something went very wrong!!!")

    filenames_codon2 = []
    for tempname in filenames:
        filenames_codon2.append(tempname.replace('.oligos','.codon2.oligos'))
        
    filenames_out_codon2 = []
    for tempname in filenames_out:
        filenames_out_codon2.append(tempname.replace('-finaloligos.fasta','-finaloligos.codon2.fasta'))

    files_for_BLAT = []
    for fnm in filenames:
        files_for_BLAT.append(fnm.replace('.oligos','-finaloligos_noAmp.fasta'))

    files_for_BLAT2 = []
    for fnm in filenames_codon2:
        files_for_BLAT2.append(fnm.replace('.oligos','-finaloligos_noAmp.fasta'))

    constructs_per_lib = [int(config.get('DS_Lib_settings', 'constructs_per_lib'))]

    barcodes = obtainFastaSequences('barcodes/'+config.get('DS_Lib_settings', 'barcode_set'))

    #check barcodes are 384, 1536, or 6144
    if not (len(barcodes) == 384 or len(barcodes) == 1536 or len(barcodes) == 6144 or len(barcodes) == 12288):
        print("Bead Barcode Length PROBLEM!!!!!!!")
    if len(barcodes) != constructs_per_lib[0]:
        print("Bead Barcode - Library Length MISMATCH!!!!!!!")

    #reverse barcodes
    barcodes_rev = Reverse(barcodes)

    # number of oligos to use to split genes in each lib:
    num_oligos_temp = int(config.get('lib_specific', 'num_oligos'))
    num_oligos = []
    for i in filenames:
        num_oligos.append(num_oligos_temp)

    #generate two different codon optimizations for each protein?
    gen_two_codon = config.get('second_lib', 'gen_two_codon') == "True"

    #should different barcodes be used for the same protein in the second library
    #this can mitigate bead effects
    switch_barcodes_for_second_lib = config.get('second_lib', 'switch_barcodes_for_second_lib') == "True"

    
    
    #length of payload + BtsaI sites + buffer such that all oligos are full length
    length_padded_payload = int(config.get('DS_Lib_settings', 'length_padded_payload'))


    primers_to_use_splitstring = config.get('DS_Lib_settings', 'primers_to_use').split(",")
    primers_to_use = []
    for i in primers_to_use_splitstring:
        primers_to_use.append(int(i))

    #ampF from 00_primer_screen.py output: skpp15 F&R for amplification primers
    #these have been further edited manually
    ampprimf_file = config.get('DS_Lib_settings', 'ampprimf_file')

    #ampR (not RC) from 00_primer_screen.py output:
    ampprimr_file = config.get('DS_Lib_settings', 'ampprimr_file')

    ampprimersf = obtainFastaSequences(ampprimf_file)
    #ampprimersf12 = ampprimersf[0:30]

    ampyFdict = dict()
    for ampy in ampprimersf:
        #skpp15-1-F filt15-1651
        ampyFdict[int(ampy.id.split("-")[1])] = ampy.seq

    #ampylist = []
    #for ampy in ampprimersf:
    #    ampylist.append(str(ampy.seq))
    #print(ampylist)

    ampprimersr = obtainFastaSequences(ampprimr_file)
    #ampprimersr12 = ampprimersr[0:30]

    ampyRdict = dict()
    for ampy in ampprimersr:
        #skpp15-1-F filt15-1651
        ampyRdict[int(ampy.id.split("-")[1])] = ampy.seq

    #make a list of the chosen primers
    ampprimersf12 = []
    ampprimersr12 = []
    for i in primers_to_use:
        print("primer: "+str(i)+" fwd: "+ampyFdict[i]+" rev: "+ampyRdict[i])
        ampprimersf12.append(SeqRecord(ampyFdict[i]))
        ampprimersr12.append(SeqRecord(ampyRdict[i]))

    #ampylist = []
    #for ampy in ampprimersr:
    #    ampylist.append(str(ampy.seq))
    #print(ampylist)


    ####################
    REsettings = dict()
    REsites_file = config.get('lib_specific', 'REsites_file')
    sitesfile = open("REsites/"+REsites_file)
    csvreader = csv.reader(sitesfile)
    sitesheader = next(csvreader)
    #print(header)
    #name   site    self_overlap_possible   3gene_no_buffer  gene_with_buffer    gene_with_asmprim   6after_split 7add_btsi    8final   random_dna_gene_padding random_dna_oligo_padding    notes
    #0name,1site,2self_overlap_possible,3gene_no_buffer,4gene_with_buffer,5gene_with_asmprim,6after_split,7final,8notes
    for row in csvreader:
        REsettings[row[0]] = [row[1],row[2],int(row[3]),int(row[4]),int(row[5]),int(row[6]),int(row[7]),int(row[8]),int(row[9]),int(row[10]),row[11]]
        #print(row[0]+': '+str(REsettings[row[0]]))
        #print(REsettings[row[0]][0])

    #print(rows)
    sitesfile.close()
    #testdict = parseREstring_for_specific_function(REsettings, sitesheader, "gene_no_buffer")
    #print(testdict)
    #testdict = parseREstring_for_specific_function(REsettings, sitesheader, "after_split")
    #print(testdict)
    ############

    append_flag = 0

    print("filenames: "+str(filenames))

    for index in range(len(filenames)):
        print("Processing Library: "+filenames[index])
        if gen_two_codon == True:
            index_mod = 2*index
        else:
            index_mod = index
        
        #read all of the split oligos
        buildoligos = getBuildOligos(filenames[index],num_oligos[index])
        #plotOligo_length(buildoligos,"plots/Oligo_Lib"+str(index)+"_length_in_hist")

        REbatch = parseREstring_for_specific_function(REsettings, sitesheader, "add_btsi")

        #add BtsI and check that everything is still good
        btsIoligos = addBtsI(buildoligos, REbatch)
        #plotOligo_length(btsIoligos,"plots/Oligo_Lib"+str(index)+"_length_after_btsI_hist")

        REbatch_random = parseREstring_for_specific_function(REsettings, sitesheader, "random_dna_oligo_padding")

        #add padding to make sure each oligo is same length
        paddedOligos = padOligoCP(btsIoligos,length_padded_payload, REbatch, REbatch_random)
        #plotOligo_length(paddedOligos,"plots/Oligo_"+str(index)+"_length_after_padding_hist")

        REbatch = parseREstring_for_specific_function(REsettings, sitesheader, "final")

        #add barcodes
        barcodedoligos = addBarcodes(paddedOligos,barcodes, REbatch)
        
        
        #save files for BLAT analysis
        if append_flag == 0:
            outputOligoLibrary(barcodedoligos,files_for_BLAT[index],'w',1)
            append_flag = 1
        else:
            outputOligoLibrary(barcodedoligos,files_for_BLAT[index],'a',1)

        #add amplification primers
        finaloligos = addAmpPrimers(barcodedoligos,ampprimersf12[index_mod].seq,ampprimersr12[index_mod].seq, REbatch)
        #printOligoLibrary(finaloligos)
        print(filenames_out[index])

        #save final oligos
        outputOligoLibrary(finaloligos,filenames_out[index],'w',1)

        if gen_two_codon == True:
            print("Processing Library: "+filenames_codon2[index])
            #read all of the split oligos
            buildoligos = getBuildOligos(filenames_codon2[index],num_oligos[index])
            #plotOligo_length(buildoligos,"plots/Oligo_Lib"+str(index)+"_length_in_hist")

            REbatch = parseREstring_for_specific_function(REsettings, sitesheader, "add_btsi")

            #add BtsI and check that everything is still good
            btsIoligos = addBtsI(buildoligos, REbatch)
            #plotOligo_length(btsIoligos,"plots/Oligo_Lib"+str(index)+"_length_after_btsI_hist")

            REbatch_random = parseREstring_for_specific_function(REsettings, sitesheader, "random_dna_oligo_padding")

            #add padding to make sure each oligo is same length
            paddedOligos = padOligoCP(btsIoligos,length_padded_payload, REbatch, REbatch_random)
            #plotOligo_length(paddedOligos,"plots/Oligo_"+str(index)+"_length_after_padding_hist")

            REbatch = parseREstring_for_specific_function(REsettings, sitesheader, "final")

            if switch_barcodes_for_second_lib == True:
                #add reversed barcodes
                barcodedoligos = addBarcodes(paddedOligos,barcodes_rev, REbatch)
            else:
                #add barcodes in normal order
                barcodedoligos = addBarcodes(paddedOligos,barcodes, REbatch)

            #save files for BLAT analysis
            if append_flag == 0:
                outputOligoLibrary(barcodedoligos,files_for_BLAT[index],'w',2)
                append_flag = 1
            else:
                outputOligoLibrary(barcodedoligos,files_for_BLAT[index],'a',2)

            #add amplification primers
            finaloligos = addAmpPrimers(barcodedoligos,ampprimersf12[index_mod+1].seq,ampprimersr12[index_mod+1].seq, REbatch)
            #printOligoLibrary(finaloligos)
            print(filenames_out_codon2[index])

            #save final oligos
            outputOligoLibrary(finaloligos,filenames_out_codon2[index],'w',2)



