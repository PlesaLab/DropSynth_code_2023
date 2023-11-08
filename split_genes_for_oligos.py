#!/usr/bin/python

#Converts protein sequences into split oligos
#allows multiple attempts for oligo splitting
#python 3.7
#updated 05.20.2022 Calin Plesa

#to run use
#python split_genes_for_oligos.py -c <configfile>
#python split_genes_for_oligos.py -c HP_Lib_1.conf

#to profile
#python -m cProfile -o profile_08.txt split_genes_for_oligos.py -c <configfile>
#echo 'stats' | python3 -m pstats profile_08.txt > profile_08_humanreadable.txt


from Bio import SeqIO
import random
random.seed(15643243242342759)
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
#from Bio import Restriction
#from Bio.Restriction import *
from typing import List, Tuple, Union, TextIO
import Primerselectiontools_py3
import csv
from Bio import pairwise2
import time
import json
import sys, getopt
import configparser
import re
import difflib

##################################
#FUNCTIONS:
def screenRestrictionSites(seqrec, REbatch):
    temp_total_str = ""
    #if (temp_str.count('GAATTC') != 0): #EcoRI
    #    print(str(temp_str.count('GAATTC'))+' EcoRI site(s) found: ' + str(temp_str.find('GAATTC')))

    if ((temp_str.find('GCAGTG') != -1) or (temp_str.find('CACTGC') != -1)): #BtsI
        temp_total_str = str(temp_str.count('GCAGTG')+temp_str.count('CACTGC'))
        print(temp_total_str+' BtsI site(s) found: ' + str(temp_str.find('GCAGTG')) + " " + str(temp_str.find('CACTGC')))

    if ((temp_str.find('GCTCTTC') != -1) or (temp_str.find('GAAGAGC') != -1)): #BspQI
        temp_total_str = str(temp_str.count('GCTCTTC')+temp_str.count('GAAGAGC'))
        print(temp_total_str+' BspQI site(s) found: ' + str(temp_str.find('GCTCTTC')) + " " + str(temp_str.find('GAAGAGC')))

    #if ((temp_str.find('GAGTC') != -1) or (temp_str.find('GACTC') != -1)): #MlyI
    #    temp_total_str = str(temp_str.count('GAGTC')+temp_str.count('GACTC'))
    #    print(temp_total_str+' MlyI site(s) found: ' + str(temp_str.find('GAGTC')) + " " + str(temp_str.find('GACTC')))

    #if ((temp_str.find('GAGGAG') != -1) or (temp_str.find('CTCCTC') != -1)): #BseRI
    #    temp_total_str = str(temp_str.count('GAGGAG')+temp_str.count('CTCCTC'))
    #    print(temp_total_str+' BseRI site(s) found: ' + str(temp_str.find('GAGGAG')) + " " + str(temp_str.find('CTCCTC')))
        
#generate random DNA and limit max homopolymer length
#http://stackoverflow.com/questions/21205836/generating-random-sequences-of-dna
def randomDNA(length, REbatch_random):
    gsflag = 0
    while gsflag == 0:
        temp_str = ''.join(random.choice('CGTA') for _ in range(length))
        #make sure no homopolymers greater than 5 in length or cut sites
        #NdeI (CATATG), KpnI (GGTACC), EcoRI (GAATTC), BtsI (GCAGTG - RC: CACTGC), BspQI (GCTCTTC - RC: GAAGAGC), MlyI (GAGTC - RC: GACTC)
        #BseRI (GAGGAG)
        if ((temp_str.find('AAAAA') == -1) and (temp_str.find('TTTTT') == -1) and
            (temp_str.find('GGGGG') == -1) and (temp_str.find('CCCCC') == -1)):
            #check others
            fail_flag = 0
            fail_flag = checkRestrictionSites(SeqRecord(Seq(temp_str),id=""), REbatch_random)
            if fail_flag != 1:
                gsflag = 1
    return temp_str

def readCodonUsage(filename):

    codon_threshold=0.149

    handle = open(filename)
    strings = []
    for line in handle:
        splitline = line.split() # breaks by whitespace
        for i in [0,5,10,15]:
            strings.append(splitline[i:i+5])
    
    codons = dict()
    for string in strings:
        if string[1] in codons: #string[1] is triplet
            dnacodon = string[0].replace('U','T')
            codons[string[1]].append([dnacodon, float(string[3])]) #string[3] is probability
        else:
            dnacodon = string[0].replace('U','T')
            codons[string[1]] = [[dnacodon, float(string[3])]]
    
    #make a fractional accounting
    for aa in codons:
        runningaacount = 0
        for i in codons[aa]:
            runningaacount += i[1] #add all probabilities
        for i in codons[aa]:
            fraction = i[1]/runningaacount #normalize
            i.append(fraction)
            #if fraction<0.15:
            #    print i
    
    #drop codons less than or equal to 15%
    for aa in codons:
        newlist = []
        for i in codons[aa]:
            if i[2]>=codon_threshold:
                newlist.append(i)
            else:
                print("dropped codon ",i[0]," with ",str(i[2]))
        codons[aa] = newlist
    
    #do running total accounting
    for aa in codons:
        runningaacount = 0
        runningtotal = 0
        for i in codons[aa]:
            runningaacount += i[1]
        for i in codons[aa]:
            fraction = i[1]/runningaacount
            i[2] = fraction
            runningtotal += fraction
            i.append(runningtotal)

    #take last codon from each and make 1.0
    for aa in codons:
        codons[aa][-1][-1] = 1.0
    
    #make reverse codon dictionary
    revcodons = dict()
    for aa in codons:
        for i in codons[aa]:
            revcodons[i[0]] = aa
            
    
    #for aa in codons:
    #    print aa + '\t\t' + str(codons[aa])
    return codons, revcodons

def obtainFastaSequences(filename):
    handle = open(filename)
    records = []
    for seqrecord in SeqIO.parse(handle, "fasta"):
        records.append(seqrecord)
        #print seqrecord.id
        #print len(seqrecord)
    #print len(records)
    return records

#do a basic optimization; not checking restriction sites or other sites to avoid
def calc_codons_no_homopol(codons):
    #make a new table for K, G, F, P with no homopolymer codons
    codons_no_hp = dict()
    temp_list = []
    for i in codons['K']:
        if i[0] != 'AAA':
            temp_list.append(i)
    codons_no_hp['K'] = temp_list
    temp_list = []
    for i in codons['G']:
        if i[0] != 'GGG':
            temp_list.append(i)
    codons_no_hp['G'] = temp_list
    temp_list = []
    for i in codons['F']:
        if i[0] != 'TTT':
            temp_list.append(i)
    codons_no_hp['F'] = temp_list
    temp_list = []
    for i in codons['P']:
        if i[0] != 'CCC':
            temp_list.append(i)      
    codons_no_hp['P'] = temp_list

    return codons_no_hp


#do a basic optimization; not checking restriction sites or other sites to avoid
def codonOptimize_two(codons, codons_no_hp, record, gentwo, furthest_codon, REbatch):
    
    sequence = ""

    if gentwo:
        sequence2 = ""

    #print(codons)
    #example of a codons record:
    #'D': [['GAT', 32.2, 0.6276803118908382, 0.6276803118908382], ['GAC', 19.1, 0.3723196881091618, 1.0]]

    #what was the previous codon chosen?
    #use this to prevent homopolymer runs
    last_codon_used =""
    last_letter=""

    key_array = []
    for key in REbatch:
        siteseq=REbatch[key][1]
        key_array.append(siteseq) #add RE site to list

    letter_index = 0
    total_seq_len = len(record.seq)
    #for each amino acid:
    for letter in record.seq:
        #print(str(letter) + '\t' + str(codons[letter][0][0]))
        if letter_index+1 < total_seq_len:
            next_letter = record.seq[letter_index+1]
        else:
            next_letter=""
        #pick random number
        rn = random.random()
        for i in codons[letter]:
            if rn < i[-1]:
                #force choice for homopolymer codons
                if (i[0] == "AAA") and (last_codon_used == "AAA"): #block A homopolymer
                    alternate_codon = random.choice(codons_no_hp['K'])
                    sequence = sequence + alternate_codon[0]
                    last_codon_used = alternate_codon[0]
                    last_letter=letter.upper()
                    if gentwo:
                        sequence2 = sequence2 + furthest_codon[alternate_codon[0]]
                    break

                elif (i[0] == "TTT") and (last_codon_used == "TTT"): #block T homopolymer
                    alternate_codon = random.choice(codons_no_hp['F'])
                    sequence = sequence + alternate_codon[0]
                    last_codon_used = alternate_codon[0]
                    last_letter=letter.upper()
                    if gentwo:
                        sequence2 = sequence2 + furthest_codon[alternate_codon[0]]
                    break

                elif (i[0] == "GGG") and (last_codon_used == "GGG"): #block G homopolymer
                    alternate_codon = random.choice(codons_no_hp['G'])
                    sequence = sequence + alternate_codon[0]
                    last_codon_used = alternate_codon[0]
                    last_letter=letter.upper()
                    if gentwo:
                        sequence2 = sequence2 + furthest_codon[alternate_codon[0]]
                    break

                elif (i[0] == "CCC") and (last_codon_used == "CCC"): #block C homopolymer
                    alternate_codon = random.choice(codons_no_hp['P'])
                    sequence = sequence + alternate_codon[0]
                    last_codon_used = alternate_codon[0]
                    last_letter=letter.upper()
                    if gentwo:
                        sequence2 = sequence2 + furthest_codon[alternate_codon[0]]
                    break
                elif ("NdeI" in REbatch) and (i[0] == "CAT") and (next_letter.upper() == "M"): #dont choose CATATG because of NdeI
                    sequence = sequence + "CAC"
                    last_codon_used = "CAC"
                    last_letter=letter.upper()
                    if gentwo:
                        sequence2 = sequence2 + "CAC"
                    break
                #BsaI GGTCTC GAGACC
                #frame 1
                elif (letter.upper() == "G") and (next_letter.upper() == "L"):#dont choose GGT because of BsaI
                    picked_codon = random.choice(["GGC","GGA","GGG"])
                    sequence = sequence + picked_codon
                    last_codon_used = picked_codon
                    last_letter=letter.upper()
                    if gentwo:
                        sequence2 = sequence2 + picked_codon
                    break
                elif (letter.upper() == "E") and (next_letter.upper() == "T"):#dont choose GAG because of BsaI
                    sequence = sequence + "GAA"
                    last_codon_used = "GAA"
                    last_letter=letter.upper()

                    if gentwo:
                        sequence2 = sequence2 + "GAA"
                    break
                #BsaI GGTCTC
                #frame 3
                elif (letter.upper() == "V") and (next_letter.upper() == "S"):#dont choose GTC because of BsaI
                    picked_codon = random.choice(["GTT","GTA","GTG"])
                    sequence = sequence + picked_codon
                    last_codon_used = picked_codon
                    last_letter=letter.upper()
                    if gentwo:
                        sequence2 = sequence2 + picked_codon
                    break

                #GCAGTG and CACTGC #BtsI
                elif (letter.upper() == "A") and (next_letter.upper() == "V"):#dont choose GCA if next aa is V because of BtsI
                    picked_codon = random.choice(["GCT","GCC","GCG"])
                    sequence = sequence + picked_codon
                    last_codon_used = picked_codon
                    last_letter=letter.upper()

                    if gentwo:
                        sequence2 = sequence2 + picked_codon
                    break
                elif (letter.upper() == "H") and (next_letter.upper() == "C"):#dont choose CAC if next aa is C because of BtsI
                    sequence = sequence + "CAT"
                    last_codon_used = "CAT"
                    last_letter=letter.upper()

                    if gentwo:
                        sequence2 = sequence2 + "CAT"
                    break
                elif (i[0] == "CAG") and (next_letter.upper() == "C" or next_letter.upper() == "W"):#dont choose GCA if next aa is V because of BtsI
                    picked_codon = "CAA"
                    sequence = sequence + picked_codon
                    last_codon_used = picked_codon
                    last_letter=letter.upper()

                    if gentwo:
                        sequence2 = sequence2 + picked_codon
                    break
                #BspQI - GCTCTTC
                #BspQIrc - GAAGAGC
                elif ("BspQI" in REbatch) and (i[0] == "CTT") and (last_codon_used == "GCT"): #block 
                    picked_codon = random.choice(["CTC","CTA","CTG"])
                    sequence = sequence + picked_codon
                    last_codon_used = picked_codon
                    last_letter=letter.upper()
                    if gentwo:
                        sequence2 = sequence2 + furthest_codon["CTT"]
                    break
                elif ("BspQI" in REbatch) and (i[0] == "TTC") and (last_codon_used == "CTC"): #block 
                    picked_codon = "TTT"
                    sequence = sequence + picked_codon
                    last_codon_used = picked_codon
                    last_letter=letter.upper()
                    if gentwo:
                        sequence2 = sequence2 + furthest_codon[picked_codon]
                    break
                elif ("BspQIrc" in REbatch) and (i[0] == "GAG") and (last_codon_used == "GAA"): #block 
                    picked_codon = "GAA"
                    sequence = sequence + picked_codon
                    last_codon_used = picked_codon
                    last_letter=letter.upper()
                    if gentwo:
                        sequence2 = sequence2 + furthest_codon[picked_codon]
                    break
                elif ("BspQIrc" in REbatch) and (i[0] == "AGC") and (last_codon_used == "AAG"): #block 
                    picked_codon = random.choice(["TCC","TCG","AGT"])
                    sequence = sequence + picked_codon
                    last_codon_used = picked_codon
                    last_letter=letter.upper()
                    if gentwo:
                        sequence2 = sequence2 + furthest_codon[picked_codon]
                    break

                else:
                    sequence = sequence + i[0]
                    last_codon_used = i[0]
                    last_letter=letter.upper()
                    if gentwo:
                        sequence2 = sequence2 + furthest_codon[i[0]]
                    break

        letter_index=letter_index+1        

    newseq = Seq(sequence)
    newrecord = SeqRecord(newseq, id=record.id)

    if gentwo:
        newseq2 = Seq(sequence2)
        newrecord2 = SeqRecord(newseq2, id=record.id)
        
        #determine percent identity between the two
        alignment = pairwise2.align.globalxx(newrecord.seq, newrecord2.seq, one_alignment_only=1, score_only=1)
        pct_ident = alignment/len(newrecord.seq)
        
        return newrecord, newrecord2, pct_ident
    else:
        return newrecord

def removeSite(seqrecord, site, codons, revcodons, REbatch, firstinstance=True):
    #find site
    seqtoremove = site
    if firstinstance:
        #this will only return FIRST occurence
        index = seqrecord.seq.find(site)
    else:
        #this will only return LAST occurence
        index = seqrecord.seq.rfind(site)

        #if both cloning sites the same this should be switched to:
        #index = seqrecord.seq.find(site,6)
        #where 6 is length to skip first instance

    if index == -1:
        #nothing found check reverse complement of the site to change
        if firstinstance:
            #this will only return FIRST occurence of the RC site
            index = seqrecord.seq.find(site.reverse_complement())
            seqtoremove = site.reverse_complement()
        else:
            #this will only return LAST occurence of the RC site
            index = seqrecord.seq.rfind(site.reverse_complement())
            seqtoremove = site.reverse_complement()        
    
    if index == -1: 
        #something went wrong could not find site
        return False

    #how long is the RE site
    recognitionlength = len(seqtoremove)
    if index % 3 == 0:
        #first base in RE site is first base in codon
        #calculate how many codons can be changed
        
        codonstochange = (((recognitionlength-1) // 3)+1) #py3
        
        #for each of the codons that can be changed
        for i in range(codonstochange):
            #starting position of the codon
            codonstart = i*3 + index
            #get codon seq
            codontochange = seqrecord.seq[codonstart:codonstart+3]
            #this will loop over all possible codons for this particular amino acid
            #revcodons will return the aa
            #codons will return all of the corresponding codons
            for possiblecodon in codons[revcodons[str(codontochange)]]:
                #make sure it's not the same codon
                if possiblecodon[0] != str(codontochange):
                    #replace the codon
                    tempseq = seqrecord.seq[:codonstart] + possiblecodon[0] + seqrecord.seq[codonstart+3:]

                    #grab the modified local sequence around the modification
                    localseq = tempseq[codonstart-3:codonstart+3*codonstochange+3]
                    tempseqrec = SeqRecord(Seq(str(localseq)),id="temp",description="none")

                    #make sure the change doesnt introduce another illegal site, but only check local region!
                    if checkRestrictionSites(tempseqrec, REbatch, withStartEndREs=False) == 0:
                        #print("site:"+str(site))
                        #print("old seq:"+str(seqrecord.seq))
                        #print("mod seq:"+str(localseq))
                        #print("mod accepted")
                        seqrecord.seq = tempseq

                        #make sure it's changed
                        if str(seqtoremove) != seqrecord.seq[index:index+recognitionlength]:
                            return True
        print("Could not change restriction site in "+seqrecord.id+" frame 1: "+seqtoremove+ " at index: "+str(index))
        return False
    elif index % 3 == 1:
        #first base in RE site is second base in codon
        #calculate how many codons can be changed
        codonstochange = (((recognitionlength) // 3)+1) #py3
        for i in range(codonstochange):
            codonstart = i*3 + index - 1
            codontochange = seqrecord.seq[codonstart:codonstart+3]
            for possiblecodon in codons[revcodons[str(codontochange)]]:
                if possiblecodon[0] != str(codontochange):
                    #replace the codon
                    tempseq = seqrecord.seq[:codonstart] + possiblecodon[0] + seqrecord.seq[codonstart+3:]

                    #grab the modified local sequence around the modification
                    localseq = tempseq[codonstart-3:codonstart+3*codonstochange+3]
                    tempseqrec = SeqRecord(Seq(str(localseq)),id="temp",description="none")

                    #make sure the change doesnt introduce another illegal site
                    if checkRestrictionSites(tempseqrec, REbatch, withStartEndREs=False) == 0:
                        seqrecord.seq = tempseq
                        #make sure it's changed
                        if str(seqtoremove) != seqrecord.seq[index:index+recognitionlength]:
                            return True
        print("Could not change restriction site in "+seqrecord.id+" frame 2: "+seqtoremove+ " at index: "+str(index))
        return False
    elif index % 3 == 2:
        # first base in RE site is third base in codon
        #calculate how many codons can be changed
        codonstochange = (((recognitionlength+1) // 3)+1) #py3
        for i in range(codonstochange):
            codonstart = i*3 + index - 2
            codontochange = seqrecord.seq[codonstart:codonstart+3]
            for possiblecodon in codons[revcodons[str(codontochange)]]:
                if possiblecodon[0] != str(codontochange):
                    #replace the codon
                    tempseq = seqrecord.seq[:codonstart] + possiblecodon[0] + seqrecord.seq[codonstart+3:]

                    #grab the modified local sequence around the modification
                    localseq = tempseq[codonstart-3:codonstart+3*codonstochange+3]
                    tempseqrec = SeqRecord(Seq(str(localseq)),id="temp",description="none")

                    #make sure the change doesnt introduce another illegal site
                    if checkRestrictionSites(tempseqrec, REbatch, withStartEndREs=False) == 0:
                        seqrecord.seq = tempseq
                        #make sure it's changed
                        if str(seqtoremove) != seqrecord.seq[index:index+recognitionlength]:
                            return True
        print("Could not change restriction site in "+seqrecord.id+" frame 3: "+seqtoremove+ " at index: "+str(index))
        #print(site)
        #print(seqrecord.seq)
        return False

def changeHomopolymerSites(seqrecord, codons, revcodons, max_homopolymer_repeat_length):
    tempseq = Seq(str(seqrecord.seq))
    newseqrecord = SeqRecord(tempseq,id=seqrecord.id,description=seqrecord.description)
    temp_str = str(seqrecord.seq)
    #are there illegal homopolymer runs
    if not (((temp_str.find('A' * max_homopolymer_repeat_length) == -1) and (temp_str.find('T' * max_homopolymer_repeat_length) == -1) and (temp_str.find('G' * max_homopolymer_repeat_length) == -1) and (temp_str.find('C' * max_homopolymer_repeat_length) == -1))):
        
        hpanalysis = dict()
        letters = ['A','T','G','C']
        #for each of the four types
        for i in range(0,3):
            key = letters[i] * max_homopolymer_repeat_length
            #try to find it
            results = temp_str.find(key)
            if results == -1: #no results
                hpanalysis[key] = []
            else: #hp found
                hpanalysis[key] = results
                seqkey = Seq(key)
                #print("removing homopolymer run")
                #try to remove it
                removeSite(newseqrecord, seqkey, codons, revcodons, REbatch)
    return newseqrecord

def occurrencesoverlap(string, sub):
    count = start = 0
    while True:
        start = string.find(sub, start) + 1
        if start > 0:
            count+=1
        else:
            return count

def changeRestrictionSites(seqrecordin, codons, revcodons, REbatch, codonversion, cloning_fwd, cloning_rev, withStartEndREs=True):
    tempseq = Seq(str(seqrecordin.seq))
    seqrecord = SeqRecord(tempseq,id=seqrecordin.id,description=seqrecordin.description)
    
    cloning_buffer_length_fwd = len(cloning_fwd)
    cloning_buffer_length_rev = len(cloning_rev)

    #stopped using BioPython Restriction for Type IIS because of issue with finding sites when cut site is past the end
    #BtsI, BspQI, , BseRI

    #not used MlyI: GACTC:GAGTC

    #key_array = ['GCAGTG', 'CACTGC', 'GCTCTTC', 'GAAGAGC', 'GAGGAG', 'CTCCTC']
    key_array = []
    key_count_allowed = dict()
    key_overlap = dict()
    name_of_seq = dict()
    for key in REbatch:
        siteseq=REbatch[key][1]
        key_array.append(siteseq) #add RE site to list
        key_overlap[siteseq] = REbatch[key][2]
        key_count_allowed[siteseq] = REbatch[key][3]
        name_of_seq[siteseq]=REbatch[key][0]

    for key in key_array:
        #print(key)
        #print(key_count_allowed[key])
        #print(key_overlap[key])
        #loop over all matches, if any
        #this site is not allowed at all:
        if key_count_allowed[key] == 0:
            for m in re.finditer(key, str(seqrecord.seq)):
                if withStartEndREs == True:
                    #print("trying to remove any instances of "+key)
                    seqkey = Seq(key)
                    #print("found ",key, " ",name_of_seq[key]," changeRestrictionSites with sequences ("+str(codonversion)+"):")
                    #print(str(seqrecordin.seq))
                    tempseq2 = Seq(str(seqrecordin.seq)[cloning_buffer_length_fwd:-cloning_buffer_length_rev])
                    #print(str(tempseq2))
                    seqrecord2 = SeqRecord(tempseq2,id=seqrecordin.id,description=seqrecordin.description)
                    removeSite(seqrecord2, seqkey, codons, revcodons, REbatch)
                    seqrecord = SeqRecord(Seq(str(seqrecordin.seq)[0:cloning_buffer_length_fwd]+str(seqrecord2.seq)+str(seqrecordin.seq)[-cloning_buffer_length_rev:]),id=seqrecordin.id,description=seqrecordin.description)
                    
                else:
                    seqkey = Seq(key)
                    #print("found ",key, " ",name_of_seq[key], " changeRestrictionSites ("+str(codonversion)+")")
                    removeSite(seqrecord, seqkey, codons, revcodons, REbatch)
        elif key_count_allowed[key] == 1:
            if key_overlap == "Yes":
                #count occurances with overlap
                count_sites = occurrencesoverlap(str(seqrecord.seq),key)
                if count_sites != 1:
                    print("Found "+str(count_sites)+" sites of "+key+" instead of 1.")
            elif key_overlap == "No":
                count_sites = str(seqrecord.seq).count(key)
                if count_sites != 1:
                    print("Found "+str(count_sites)+" sites of "+key+" instead of 1.")
    
    return seqrecord

def checkRestrictionSites(seqrecord, REbatch, withStartEndREs=True):
    if isinstance(seqrecord, (SeqRecord)):
        temp_str = str(seqrecord.seq)
    elif isinstance(seqrecord, (Seq)):
        temp_str = str(seqrecord)
    else:
        print("Variable type error checkRestrictionSites")
    fail_flag = 0

    key_array = []
    key_count_allowed = dict()
    key_overlap = dict()
    for key in REbatch:
        siteseq=REbatch[key][1] #site sequence
        key_array.append(siteseq) #add RE site to list
        key_overlap[siteseq] = REbatch[key][2] # Yes or No
        key_count_allowed[siteseq] = REbatch[key][3] #Number allowed

    for key in key_array:
        #loop over all matches, if any
        #this site is not allowed at all:
        if key_count_allowed[key] == 0:
            for m in re.finditer(key, temp_str):
                print("Found "+key+" sites in:")
                print(temp_str)
                fail_flag = 1
        #these sites are allowed once
        elif key_count_allowed[key] > 0:
            if key_overlap[key] == "Yes":
                #count occurances with overlap
                count_sites = occurrencesoverlap(temp_str,key)
                if count_sites != key_count_allowed[key]:
                    print("Found "+str(count_sites)+" sites of "+key+" instead of "+str(key_count_allowed[key]))
                    fail_flag = 1
            elif key_overlap[key] == "No":
                count_sites = temp_str.count(key)
                if count_sites != key_count_allowed[key]:
                    print("Found "+str(count_sites)+" sites of "+key+" instead of "+str(key_count_allowed[key]))
                    fail_flag = 1

    return fail_flag

def checkCodingRegions(aarecord, seqrecord):
    numberOfMatchedProteins = 0
    newproteinseq = seqrecord.seq.translate()
    oldproteinseq = aarecord.seq
    if str(newproteinseq) == str(oldproteinseq):
        numberOfMatchedProteins += 1
    else:
        print(seqrecord.id+" bad codon assignment, translation mismatch")
    #print str(numberOfMatchedProteins) + " proteins with correct aa sequence redesigned"
    return numberOfMatchedProteins

def addBufferREsites(seqrecordin,add_stop_codon,cloning_fwd,cloning_rev,clonging_fwd_includes_ATG, check_starting_ATG):
    forward_buffer = cloning_fwd # before start codon ATG = NdeI
    #print("fwd buff:"+str(forward_buffer))
    if add_stop_codon:
        reverse_buffer = "TAA"+cloning_rev # KpnI
    else:
        reverse_buffer = cloning_rev # KpnI
    #print(str(seqrecordin.seq[0:3]).upper())
    if (not ("ATG" == str(seqrecordin.seq[0:3]).upper()) and check_starting_ATG):
        print("Error! Sequence does not start with ATG. "+str(seqrecordin.id))
    newseq = Seq(forward_buffer + str(seqrecordin.seq) + reverse_buffer)
    record = SeqRecord(newseq,id=seqrecordin.id,description=seqrecordin.description)
    
    return record

def addAssemblyPrimers(seqrecordin, forward_buffer, reverse_buffer, padding_var, padding_length, REbatch_random):
    seqlength = len(seqrecordin.seq)
    #print("Length of seq before padding: "+str(seqlength)+", Padding variable:"+str(padding_var)+", Padding length:"+str(padding_length))
    pad_added_length = 0
    #padding is turned off OR padding is not needed because sequence is already long enough
    if (padding_var == False) or ((padding_var == True) and (seqlength >= padding_length)):
        tempseq = forward_buffer.seq + seqrecordin.seq + reverse_buffer.seq.reverse_complement()
        #tempseq = Seq(str(forward_buffer + seqrecordin.seq + reverse_buffer))
    #padding is on and needed
    else:
        #generate random sequence
        pad_random_seq = randomDNA(padding_length-len(seqrecordin.seq), REbatch_random)
        #the padding goes in front of the reverse primer like this
        tempseq = forward_buffer.seq + seqrecordin.seq + Seq(pad_random_seq) + reverse_buffer.seq.reverse_complement()
        #tempseq = Seq(str(forward_buffer + seqrecordin.seq + pad_random_seq + reverse_buffer))
        pad_added_length = len(pad_random_seq)
        #print("Length of seq after padding: "+str(len(tempseq)))
    seqrecord = SeqRecord(tempseq,id=seqrecordin.id,description=seqrecordin.description)
    return seqrecord, pad_added_length


def count_with_overlap(seqin, site):
  
    numsites = 0
    index = 0

    while index < len(seqin):
  
        # this will find the first occurence of site in seqin
        pos = seqin.find(site, index)
  
        if pos != -1: #if there is a match
            # move index past start of first site
            index = pos + 1
            numsites += 1

        else: #nothing found
            break

    return numsites

def parseREstring_for_specific_function(REsettings, sitesheader, target_step):

    #output is a dictionary containing
    #name, site, number of times to check
    REout = dict()

    index_of_target = sitesheader.index(target_step)-1

    for key in REsettings:
        #if value is positive check it
        if REsettings[key][index_of_target] > -1:
            #name, site, 
            REout[key] = [key, REsettings[key][0], REsettings[key][1], int(REsettings[key][index_of_target])]

    return REout


# #https://stackoverflow.com/questions/14128763/how-to-find-the-overlap-between-2-sequences-and-return-it
# def LongestCommonSubstring(s1, s2):
#     s = difflib.SequenceMatcher(None, s1, s2)
#     pos_a, pos_b, size = s.find_longest_match(0, len(s1), 0, len(s2)) 
#     return s1[pos_a:pos_a+size]

def checkOligos(oligosin, num_oligos, REbatch):

    #very important notes:
    #if you use biopython restriction, it will not find sites for type IIS cutters with cut outside sequence
    #if you use .count() it will not find overlapping sites
    # .find() will return index of first match or -1 otherwise
    #use count_with_overlap(seqin, site) to count sites that can overlap
    
    maxolgiolength = 0
    length_in_case_of_error = 100000
    #find length of longest oligo after split:
    for oligo in oligosin:
        if len(oligo) > maxolgiolength:
            maxolgiolength = len(oligo)
            
    checkOligos_msg = ""
    
    #check that adding BtsaI sites doesn't introduce new restriction sites
    #this loop will start by checking the LAST oligo and then decrement until the first oligo

    num_oligosin = len(oligosin)

    fail_flag = dict()
    for i in range(0,num_oligosin):
        fail_flag[i] = 0

    for i in reversed(range(0,num_oligosin)):
        newoligo = Seq("GCAGTG") + oligosin[i] + Seq("CACTGC")

        #fixed old bug here

        fail_flag[i] = checkRestrictionSites(newoligo, REbatch)
        if fail_flag[i] == 1:
            maxolgiolength = length_in_case_of_error
            if num_oligos == 1:
                if checkOligos_msg == "":
                    checkOligos_msg = "Bad number of restriction sites in only oligo"
                else:
                    checkOligos_msg = checkOligos_msg + ", only oligo"
            elif num_oligos == 2:
                if i==1:
                    if checkOligos_msg == "":
                        checkOligos_msg = "Bad number of restriction sites in last oligo (2 oligo split)"
                    else:
                        checkOligos_msg = checkOligos_msg + ", last oligo"
                elif i==0:
                    if checkOligos_msg == "":
                        checkOligos_msg = "Bad number of restriction sites in first oligo (2 oligo split)"
                    else:
                        checkOligos_msg = checkOligos_msg + ", first oligo"
            elif num_oligos >2:
                if i==0: #first
                    if checkOligos_msg == "":
                        checkOligos_msg = "Bad number of restriction sites in first oligo"
                    else:
                        checkOligos_msg = checkOligos_msg + ", first oligo"
                elif (i==(len(oligosin)-1)): #last
                    if checkOligos_msg == "":
                        checkOligos_msg = "Bad number of restriction sites in last oligo"
                    else:
                        checkOligos_msg = checkOligos_msg + ", last oligo"
                elif ((i==(num_oligosin-2)) and (num_oligosin > 2)):#second to last
                    if checkOligos_msg == "":
                        checkOligos_msg = "Bad number of restriction sites in second to last oligo"
                    else:
                        checkOligos_msg = checkOligos_msg + ", second to last oligo"

                else:#middle oligos
                    if checkOligos_msg == "":
                        checkOligos_msg = "Bad number of restriction sites in middle oligo"
                    else:
                        checkOligos_msg = checkOligos_msg + ", middle oligo"
        temp_str = str(oligosin[i])
        if temp_str.startswith("CAGTG") or temp_str.startswith("AAGAGC") or temp_str.startswith("CTCTTC") or temp_str.endswith("CACTG") or temp_str.endswith("GAAGAG") or temp_str.endswith("GCTCTT"):
            maxolgiolength = length_in_case_of_error
            if checkOligos_msg == "":
                checkOligos_msg = "Self overlapping BtsI or BspQI site over split"
            else:
                checkOligos_msg = checkOligos_msg + ", Self overlapping BtsI or BspQI site over split"

    return maxolgiolength, checkOligos_msg
    

class Settings():
    
    #codon usage table
    codons: dict
    revcodons: dict

    #make a table with non homopolymer codons
    codons_no_hp: dict

    # read table for each codon which other codon encodes same aa but is most distant
    furthest_codon: dict
    
    #these are the output file base names for each library, one per 384 lib:
    output_file_names: str
    output_folder: str

    #the source file with all sequences to be split
    input_file: str
    
    #are these proteins or DNA sequences
    input_type: str
    
    # number of oligos to use to split genes in each lib:
    num_oligos: int
    
    # number of oligos to use to split genes in each lib:
    cloning_fwd: str

    # 
    check_starting_ATG: bool

    # number of oligos to use to split genes in each lib:
    clonging_fwd_includes_ATG: bool

    # number of oligos to use to split genes in each lib:
    cloning_rev: str

    #total length of oligos to be used
    total_oligo_length: int

    #how much of the oligo is needed for processing (non-payload)
    seq_req_for_processing: int

    # max space on oligo between BtsaI cut sites (160 for 230mer):
    max_payload_length: int

    #maximum tries to split each gene (default 14):
    max_num_attempts: int
    
    #generate two different codon optimizations for each protein?
    gen_two_codon: bool
    
    #add a stop codon before KpnI?
    add_stop_codon: bool
    
    #use padding (add random sequence between final BtsI site and the reverse assembly primer)
    #this makes gel isolation easier since the distribution of lengths is reduced
    padding_var: bool
    
    #padding lengths for given number of oligos
    #sequences with lengths below this will get padding
    padding_options: List[str]
    padding_length_dict: dict

    padding_length: int
    
    #block any homopolymer repeats of this length or longer:
    max_homopolymer_repeat_length: int
    
    #this is a list of all seq IDs which don't work
    IDs_to_remove: List[str]
    
    #asmF primer file (only primer 504 used)
    #DO NOT USE skpp20-511 this is used in the anchor oligo!
    #Primers for alternate codon versions are offset by len(num_oligos): 15
    assemblyprimf_file: str
    assemblyprimf: List[SeqRecord]

    #enter the reverse primers (only primer 504 used)
    #Primers for alternate codon versions are offset by len(num_oligos): 15
    assemblyprimr_file: str
    assemblyprimr: List[SeqRecord]
    

    #codon2 forward assembly primers
    assemblyprimf_codon2: List[SeqRecord]

    #codon2 reverse assembly primers
    assemblyprimr_codon2: List[SeqRecord]

    seeds_for_libs: int

    #Primerselectiontools_settings
    lengthleeway: int
    positionleeway: int
    avgoverlapsize: int
    
    overlaptemps: List[int]
    deltaGThreshold: int
    selfDimersThreshold: int

    REsettings: dict
    sitesheader: List[str]


def parseSettings(config: dict) -> Settings:
    settings = Settings()

    #codon usage table
    settings.codons, settings.revcodons = readCodonUsage(config.get('lib_specific', 'codon_usage_file'))

    print(settings.codons)

    print(settings.revcodons)

    #make a table with non homopolymer codons
    settings.codons_no_hp = calc_codons_no_homopol(settings.codons)

    # read table for each codon which other codon encodes same aa but is most distant
    settings.furthest_codon = json.load(open(config.get('second_lib', 'furthest_codon_file')))
    
    #these are the output file base names for each library, one per 384 lib:
    settings.output_file_names = config.get('lib_specific', 'output_file_names')
    settings.output_folder = config.get('lib_specific', 'output_folder')

    #the source file with all sequences to be split
    settings.input_file = config.get('lib_specific', 'input_file')
    
    #are these proteins or DNA sequences
    settings.input_type = config.get('lib_specific', 'input_type')
    
    # number of oligos to use to split genes in each lib:
    settings.num_oligos = int(config.get('lib_specific', 'num_oligos'))
    
    # number of oligos to use to split genes in each lib:
    settings.cloning_fwd = config.get('lib_specific', 'cloning_fwd')

    # 
    settings.check_starting_ATG = config.get('lib_specific', 'check_starting_ATG') == "True"

    # number of oligos to use to split genes in each lib:
    settings.clonging_fwd_includes_ATG = config.get('lib_specific', 'clonging_fwd_includes_ATG') == "True"

    # number of oligos to use to split genes in each lib:
    settings.cloning_rev = config.get('lib_specific', 'cloning_rev')

    #total length of oligos to be used
    settings.total_oligo_length = int(config.get('general_settings', 'total_oligo_length'))

    #how much of the oligo is needed for processing (non-payload)
    settings.seq_req_for_processing = int(config.get('general_settings', 'seq_req_for_processing'))

    # max space on oligo between BtsaI cut sites (160 for 230mer):
    settings.max_payload_length = settings.total_oligo_length - settings.seq_req_for_processing

    #maximum tries to split each gene (default 14):
    settings.max_num_attempts = int(config.get('general_settings', 'max_num_attempts'))
    
    #generate two different codon optimizations for each protein?
    settings.gen_two_codon = config.get('second_lib', 'gen_two_codon') == "True"
    
    #add a stop codon before KpnI?
    settings.add_stop_codon = config.get('lib_specific', 'add_stop_codon') == "True"
    
    #use padding (add random sequence between final BtsI site and the reverse assembly primer)
    #this makes gel isolation easier since the distribution of lengths is reduced
    settings.padding_var = config.get('general_settings', 'padding_var') == "True"
    
    #padding lengths for given number of oligos
    #sequences with lengths below this will get padding
    settings.padding_options = config.get('general_settings', 'padding_length').split(",")
    settings.padding_length_dict = dict()
    for i in settings.padding_options:
        po_spvals = i.split(":")
        settings.padding_length_dict[po_spvals[0].replace(" ", "")] = int(po_spvals[1].replace(" ", ""))
    settings.padding_length = settings.padding_length_dict[str(settings.num_oligos)]
    
    #block any homopolymer repeats of this length or longer:
    settings.max_homopolymer_repeat_length = int(config.get('general_settings', 'max_homopolymer_repeat_length'))
    
    #this is a list of all seq IDs which don't work
    settings.IDs_to_remove = config.get('lib_specific', 'IDs_to_remove')
    
    #asmF primer file (only primer 504 used)
    #DO NOT USE skpp20-511 this is used in the anchor oligo!
    #Primers for alternate codon versions are offset by len(num_oligos) = 15
    assemblyprimf_file = config.get('general_settings', 'assemblyprimf_file')
    settings.assemblyprimf = obtainFastaSequences(assemblyprimf_file)[0]

    #enter the reverse primers (only primer 504 used)
    #Primers for alternate codon versions are offset by len(num_oligos) = 15
    assemblyprimr_file = config.get('general_settings', 'assemblyprimr_file')
    settings.assemblyprimr = obtainFastaSequences(assemblyprimr_file)[0]
    
    #use different assembly primers for the codon2 library?
    settings.diff_assembly_pim_on_two_codon = config.get('general_settings', 'diff_assembly_pim_on_two_codon') == "True"
    if settings.diff_assembly_pim_on_two_codon:
        #codon2 forward assembly primers
        assemblyprimf_file_codon2 = config.get('general_settings', 'assemblyprimf_file_codon2')
        settings.assemblyprimf_codon2 = obtainFastaSequences(assemblyprimf_file_codon2)[0]

        #codon2 reverse assembly primers
        assemblyprimr_file_codon2 = config.get('general_settings', 'assemblyprimr_file_codon2')
        settings.assemblyprimr_codon2 = obtainFastaSequences(assemblyprimr_file_codon2)[0]


    settings.seeds_for_libs = int(config.get('lib_specific', 'seeds_for_libs'))

    #Primerselectiontools_settings
    settings.lengthleeway = int(config.get('Primerselectiontools_settings', 'lengthleeway'))
    settings.positionleeway = int(config.get('Primerselectiontools_settings', 'positionleeway'))
    settings.avgoverlapsize = int(config.get('Primerselectiontools_settings', 'avgoverlapsize'))
    
    overlaptemps_splitstring = config.get('Primerselectiontools_settings', 'overlaptemps').split(",")
    settings.overlaptemps = [int(overlaptemps_splitstring[0]), int(overlaptemps_splitstring[1])]
    settings.deltaGThreshold = int(config.get('Primerselectiontools_settings', 'deltaGThreshold'))
    settings.selfDimersThreshold = int(config.get('Primerselectiontools_settings', 'selfDimersThreshold'))

    ####################
    settings.REsettings = dict()
    REsites_file = config.get('lib_specific', 'REsites_file')
    sitesfile = open("REsites/"+REsites_file)
    csvreader = csv.reader(sitesfile)
    settings.sitesheader = next(csvreader)
    #print(header)
    #0name,1site,2self_overlap_possible,3gene_no_buffer,4gene_with_buffer,5gene_with_asmprim,6after_split,7final,8notes
    for row in csvreader:
        settings.REsettings[row[0]] = [row[1],row[2],int(row[3]),int(row[4]),int(row[5]),int(row[6]),int(row[7]),int(row[8]),int(row[9]),int(row[10]),row[11]]
        #print(row[0]+': '+str(REsettings[row[0]]))
        #print(REsettings[row[0]][0])

    #print(rows)
    sitesfile.close()
    #testdict = parseREstring_for_specific_function(REsettings, sitesheader, "gene_no_buffer")
    #print(testdict)
    #testdict = parseREstring_for_specific_function(REsettings, sitesheader, "after_split")
    #print(testdict)
    ############


    return settings

def main(argv):
   configfile = ''
   try:
      opts, args = getopt.getopt(argv,"hc:",["configfile="])
   except getopt.GetoptError:
      print('split_genes_for_oligos.py -c <configfile>')
      sys.exit(2)
   for opt, arg in opts:
      if opt == '-h':
         print('split_genes_for_oligos.py -c <configfile>')
         sys.exit()
      elif opt in ("-c", "--config"):
         inputfile = arg
   print('Config file is: ', inputfile)
   return inputfile


######################################################################
if __name__ == '__main__':
    start_time = time.time()

    config_path = main(sys.argv[1:])
    config = configparser.ConfigParser()
    config.read_file(open("lib_config/"+config_path))

    settings = parseSettings(config)

    #seq fold = 1, unafold = 2
    hybridization_method = 1

    badIDsarray = []
    badIDfasta = []
    with open(settings.IDs_to_remove, 'r') as f:
            reader = csv.reader(f)
            lib_ident_list = list(reader)
    
    if len(lib_ident_list) > 0:
        badIDsarray = lib_ident_list[0]
    else:
        badIDsarray = []

    aarec_counter = 0
    bad_rec_counter = 0

    #obtain the protein records
    aarecords = obtainFastaSequences(settings.input_file)

    num_records = len(aarecords)

    print("Generating "+str(settings.num_oligos)+" oligo splits for ("+str(num_records)+" records):")
    random.seed(settings.seeds_for_libs)
    
    final_records_oligos = []
    final_records_gene_level = []
    final_records_gene_level_no_primer_or_RE = []
    final_records_gene_level_no_primer = []
    if settings.input_type == 'AA':
        final_records_protein_level = []
    if settings.input_type == 'DNA':
        final_records_dna_level = []

    updatedseqrecord = []
    seqrecord = []
    
    pct_ident_array = []
    padding_length_array = []
    write_extra_oligo_files = False
    
    if settings.gen_two_codon and settings.input_type == 'AA':
        final_records_oligos2 = []
        final_records_gene_level2 = []
        final_records_gene_level_no_primer_or_RE2 = []
        final_records_gene_level_no_primer2 = []
        final_records_protein_level2 = []
    
    final_records_oligos_extra_oligo = []
    final_records_gene_level_extra_oligo = []
    final_records_gene_level_no_primer_or_RE_extra_oligo = []
    final_records_gene_level_no_primer_extra_oligo = []
    if settings.input_type == 'AA':
        final_records_protein_level_extra_oligo = []
    if settings.input_type == 'DNA':
        final_records_dna_level_extra_oligo = []
    if settings.gen_two_codon and settings.input_type == 'AA':
        final_records_oligos_extra_oligo2 = []
        final_records_gene_level_extra_oligo2 = []
        final_records_gene_level_no_primer_or_RE_extra_oligo2 = []
        final_records_gene_level_no_primer_extra_oligo2 = []
        final_records_protein_level_extra_oligo2 = []
        
    total_rec_count = len(aarecords)
    
    ############ Protein Record Loop Start ####################
    #for each protein record:
    for ii in range(num_records):
        #get the next protein sequence
        caarec = aarecords[ii]
        #print(caarec.id)
        num_tries = 1
        while (num_tries < settings.max_num_attempts):
            #only translate for AA not DNA
            if settings.input_type == 'AA':
                #print("Starting "+caarec.id+" (construct # "+str(rec_count)+" of "+str(total_rec_count)+", Lib "+") try # "+str(num_tries))
                REbatch = parseREstring_for_specific_function(settings.REsettings, settings.sitesheader, "gene_no_buffer")
                #generate random weighted codon usage nt record:
                if settings.gen_two_codon:
                    seqrecord,seqrecord2,pct_ident = codonOptimize_two(settings.codons, settings.codons_no_hp, caarec, settings.gen_two_codon, settings.furthest_codon, REbatch)
                else:
                    seqrecord = codonOptimize_two(settings.codons, settings.codons_no_hp, caarec, settings.gen_two_codon, settings.furthest_codon, REbatch)
                
                #remove NdeI, EcoRI, KpnI, BtsI, BspQI:
                REbatch = parseREstring_for_specific_function(settings.REsettings, settings.sitesheader, "gene_no_buffer")
                seqrecord_noRE = changeRestrictionSites(seqrecord, settings.codons, settings.revcodons, REbatch, 1, settings.cloning_fwd, settings.cloning_rev, withStartEndREs=False)
                if settings.gen_two_codon:
                    seqrecord_noRE2 = changeRestrictionSites(seqrecord2, settings.codons, settings.revcodons, REbatch, 2, settings.cloning_fwd, settings.cloning_rev, withStartEndREs=False)
                
                #remove homopolymers:
                updatedseqrecord = changeHomopolymerSites(seqrecord_noRE, settings.codons, settings.revcodons, settings.max_homopolymer_repeat_length)
                if settings.gen_two_codon:
                    updatedseqrecord2 = changeHomopolymerSites(seqrecord_noRE2, settings.codons, settings.revcodons, settings.max_homopolymer_repeat_length)

                #make sure the translation is still the same:
                numMatch = checkCodingRegions(caarec, updatedseqrecord)
                if settings.gen_two_codon and numMatch == 1:
                    numMatch = checkCodingRegions(caarec, updatedseqrecord2)
            
            if settings.input_type == 'DNA':
                REbatch = parseREstring_for_specific_function(settings.REsettings, settings.sitesheader, "gene_no_buffer")
                updatedseqrecord = caarec
                #make sure there are no bad restriction sites
                numMatch = checkRestrictionSites(updatedseqrecord, REbatch, withStartEndREs=False)
                #flip for DNA
                if numMatch == 1:
                    numMatch = 0
                else:
                    numMatch = 1
            
            #make sure no long homopolymers, 
            temp_str = str(updatedseqrecord.seq)
            if settings.gen_two_codon and settings.input_type == 'AA':
                #concatenate both records together
                temp_str = str(updatedseqrecord.seq)+str(updatedseqrecord2.seq)
            
            gsflag = 0
            
            if ((temp_str.find('A' * settings.max_homopolymer_repeat_length) == -1) and (temp_str.find('T' * settings.max_homopolymer_repeat_length) == -1) and
                (temp_str.find('G' * settings.max_homopolymer_repeat_length) == -1) and (temp_str.find('C' * settings.max_homopolymer_repeat_length) == -1)):
                gsflag = 1

            #if everything is good:
            if numMatch == 1 and gsflag == 1:
                
                #add NdeI at start and Stop codon + KpnI at end:
                bufferedupdatedseqrecord = addBufferREsites(updatedseqrecord, settings.add_stop_codon, settings.cloning_fwd, settings.cloning_rev, settings.clonging_fwd_includes_ATG, settings.check_starting_ATG)
                if settings.gen_two_codon and settings.input_type == 'AA':
                    bufferedupdatedseqrecord2 = addBufferREsites(updatedseqrecord2, settings.add_stop_codon, settings.cloning_fwd, settings.cloning_rev, settings.clonging_fwd_includes_ATG, settings.check_starting_ATG)
                
                if settings.input_type == 'AA':
                    REbatch = parseREstring_for_specific_function(settings.REsettings, settings.sitesheader, "gene_with_buffer")
                    
                    #change codons if BtsaI, BspQI, or EcoRI are present. Change codons if >1 NdeI or KpnI site present.
                    finalupdatedseqrecord = changeRestrictionSites(bufferedupdatedseqrecord, settings.codons, settings.revcodons, REbatch, 1, settings.cloning_fwd, settings.cloning_rev, withStartEndREs=True)
                    if settings.gen_two_codon and settings.input_type == 'AA':
                        finalupdatedseqrecord2 = changeRestrictionSites(bufferedupdatedseqrecord2, settings.codons, settings.revcodons, REbatch, 2, settings.cloning_fwd, settings.cloning_rev, withStartEndREs=True)
                        
                else: #DNA case
                    finalupdatedseqrecord = bufferedupdatedseqrecord
                
                REbatch = parseREstring_for_specific_function(settings.REsettings, settings.sitesheader, "random_dna_gene_padding")
                pad_length_added = 0
                #add the assembly primers and also add padding with random sequence if necessary
                finalrecwithAssPrimers, pad_length_added = addAssemblyPrimers(finalupdatedseqrecord, settings.assemblyprimf, settings.assemblyprimr, settings.padding_var, settings.padding_length, REbatch)
                
                if settings.gen_two_codon and settings.input_type == 'AA':
                    pad_length_added2 = 0
                    #use different assembly primers for codon2 than codon1?
                    if settings.diff_assembly_pim_on_two_codon == True:
                        #different primers
                        finalrecwithAssPrimers2, pad_length_added2 = addAssemblyPrimers(finalupdatedseqrecord2, settings.assemblyprimf_codon2, settings.assemblyprimr_codon2, settings.padding_var, settings.padding_length, REbatch)
                    else:
                        #same primers
                        finalrecwithAssPrimers2, pad_length_added2 = addAssemblyPrimers(finalupdatedseqrecord2, settings.assemblyprimf, settings.assemblyprimr, settings.padding_var, settings.padding_length, REbatch)
                
                RE_fail_flag2 = 0

                REbatch = parseREstring_for_specific_function(settings.REsettings, settings.sitesheader, "gene_with_asmprim")
                # make sure 1 NdeI site, 1 KpnI site, and none of EcoRI, BtsaI, BspQI:
                RE_fail_flag2 = checkRestrictionSites(finalrecwithAssPrimers, REbatch, withStartEndREs=True)
                if settings.gen_two_codon and RE_fail_flag2 == 0 and settings.input_type == 'AA':
                    RE_fail_flag2 = checkRestrictionSites(finalrecwithAssPrimers2, REbatch, withStartEndREs=True)
                
                #if no extra restriction sites found:
                if RE_fail_flag2 == 0:
                    screenRestrictionSites(finalrecwithAssPrimers, REbatch)
                    if settings.gen_two_codon and settings.input_type == 'AA':
                        screenRestrictionSites(finalrecwithAssPrimers2, REbatch)
                    
                    if settings.num_oligos > 1:
                        #this is where genes are split into oligos:
                        #input parameters are:
                        #seq, oligosizemax, lengthleeway, positionleeway, avgoverlapsize, overlaptemps, deltaGThreshold, selfDimersThreshold, num_of_oligos
                        oligos = Primerselectiontools_py3.optimizedSplit(finalrecwithAssPrimers.seq, settings.max_payload_length, settings.lengthleeway, settings.positionleeway, settings.avgoverlapsize, settings.overlaptemps, settings.deltaGThreshold, settings.selfDimersThreshold, settings.num_oligos, hybridization_method)
                        if settings.gen_two_codon and settings.input_type == 'AA':
                            oligos2 = Primerselectiontools_py3.optimizedSplit(finalrecwithAssPrimers2.seq, settings.max_payload_length, settings.lengthleeway, settings.positionleeway, settings.avgoverlapsize, settings.overlaptemps, settings.deltaGThreshold, settings.selfDimersThreshold, settings.num_oligos, hybridization_method)
                    else: #dont split, this is used for small proteins which can fit on one oligo
                        oligos = [finalrecwithAssPrimers.seq]
                        if settings.gen_two_codon and settings.input_type == 'AA':
                            oligos2 = [finalrecwithAssPrimers2.seq]
                    
                    REbatch = parseREstring_for_specific_function(settings.REsettings, settings.sitesheader, "add_btsi") #changed

                    checkOligos_msg = ""

                    maxolgiolength, checkOligos_msg = checkOligos(oligos, settings.num_oligos, REbatch)
                    
                    #make sure max oligo size and number of oligos is good:
                    #note that num oligos +1 is allowed and saved into a special file
                    #genes close to the length cut off will often split into one more oligo
                    if (not ((len(oligos) == settings.num_oligos) or (len(oligos) == settings.num_oligos+1)) or (maxolgiolength > settings.max_payload_length)):
                        #you will get to this section if:
                        #you have incorrect number of oligos
                        #or extra oligos
                        #or the maxoligolength variable is off
                        error_msg = ""

                        if (not ((len(oligos) == settings.num_oligos)) or (len(oligos) == settings.num_oligos+1)):
                            error_msg = "Incorrect number of oligos: " + str(len(oligos))
                        if maxolgiolength > settings.max_payload_length:
                            error_msg = error_msg + " Max oligo length: " + str(maxolgiolength) + " of payload limit: "+ str(settings.max_payload_length)

                        print(caarec.id+" failed oligo splitting on try # "+str(num_tries)+". "+checkOligos_msg +error_msg)
                    else:
                        #you will enter here if everything is good with the codon1 oligos
                        codon2_good_flag = False
                        if settings.gen_two_codon and settings.input_type == 'AA':
                            #you are making two codon libraries of proteins
                            checkOligos_msg = ""
                            maxolgiolength,checkOligos_msg = checkOligos(oligos2, settings.num_oligos, REbatch)
                            #print(maxolgiolength)
                            #this checks if the codon2 oligos worked
                            if ((not ((len(oligos2) == settings.num_oligos) or #wrong number of oligos
                                (len(oligos2) == settings.num_oligos+1))) or #needs more oligos
                                (len(oligos2) != len(oligos)) or #number of oligos between the two codon versions doesnt match
                                (maxolgiolength > settings.max_payload_length)):
                                
                                error_msg2 = ""

                                if not ((len(oligos2) == settings.num_oligos) or (len(oligos2) == settings.num_oligos+1)):
                                    error_msg2 = "Incorrect number of oligos: " + str(len(oligos2))
                                if len(oligos2) != len(oligos):
                                    error_msg2 = error_msg2 + " Number of oligos dont match codon1: " + str(len(oligos2)) + " codon2: " + str(len(oligos))
                                if maxolgiolength > settings.max_payload_length:
                                    error_msg2 = error_msg2 + " Max oligo length: " + str(maxolgiolength) + " of payload limit: "+ str(settings.max_payload_length)

                                print(caarec.id+" failed oligo codon2 splitting on try # "+str(num_tries)+". "+checkOligos_msg +error_msg2)

                            else:
                                #everything worked
                                codon2_good_flag = True

                                ########################################
                                #split into correct number of oligos 
                                if len(oligos2) == settings.num_oligos:
                                    
                                    #save all correct records:
                                    final_records_oligos2.append([finalrecwithAssPrimers2.id, oligos2])
                                    final_records_gene_level2.append(finalrecwithAssPrimers2)
                                    final_records_protein_level2.append(caarec)
                                    final_records_gene_level_no_primer_or_RE2.append(updatedseqrecord2)
                                    final_records_gene_level_no_primer2.append(finalupdatedseqrecord2)
                                    
                                    pct_ident_array.append(pct_ident)
                                    padding_length_array.append(pad_length_added)

                                ########################################
                                elif len(oligos2) == settings.num_oligos+1:#split has extra oligo

                                    print(caarec.id+ " split needs extra oligo! Saved in extra_split files")

                                    write_extra_oligo_files = True
                                    
                                    #save all correct records:
                                    final_records_oligos_extra_oligo2.append([finalrecwithAssPrimers2.id, oligos2])
                                    final_records_gene_level_extra_oligo2.append(finalrecwithAssPrimers2)
                                    final_records_protein_level_extra_oligo2.append(caarec)
                                    final_records_gene_level_no_primer_or_RE_extra_oligo2.append(updatedseqrecord2)
                                    final_records_gene_level_no_primer_extra_oligo2.append(finalupdatedseqrecord2)
                                    
                                    #pct_ident_array.append(pct_ident)
                                    #padding_length_array.append(pad_length_added)

                                ########################################
                        
                        if (codon2_good_flag and settings.gen_two_codon and settings.input_type == 'AA') or (settings.input_type == 'DNA') or ((not settings.gen_two_codon) and settings.input_type == 'AA'):
                            #everything worked
                            aarec_counter += 1

                            ########################################
                            #split into correct number of oligos 
                            if len(oligos) == settings.num_oligos:

                                #save all correct records:
                                final_records_oligos.append([finalrecwithAssPrimers.id, oligos])
                                final_records_gene_level.append(finalrecwithAssPrimers)
                                final_records_gene_level_no_primer.append(finalupdatedseqrecord)
                                final_records_gene_level_no_primer_or_RE.append(updatedseqrecord)

                                #these only needed when starting from proteins
                                if settings.input_type == 'AA':
                                    final_records_protein_level.append(caarec)
                                    
                                if settings.input_type == 'DNA':
                                    final_records_dna_level.append(caarec)

                            ########################################
                            elif len(oligos) == settings.num_oligos+1: #split has an extra oligo
                                
                                print(caarec.id+ " split needs extra oligo! Saved in extra_split files")

                                write_extra_oligo_files = True

                                #save all correct records:
                                final_records_oligos_extra_oligo.append([finalrecwithAssPrimers.id, oligos])
                                final_records_gene_level_extra_oligo.append(finalrecwithAssPrimers)
                                final_records_gene_level_no_primer_extra_oligo.append(finalupdatedseqrecord)
                                final_records_gene_level_no_primer_or_RE_extra_oligo.append(updatedseqrecord)

                                #these only needed when starting from proteins
                                if settings.input_type == 'AA':
                                    final_records_protein_level_extra_oligo.append(caarec)
                                    
                                if settings.input_type == 'DNA':
                                    final_records_dna_level_extra_oligo.append(caarec)

                            ########################################
                            else:
                                print("Something has gone terribly wrong. ID: "+caarec.id+" try # "+str(num_tries))

                            break

                else:
                    print(caarec.id+ " failed checkRestrictionSites after adding RE, AsmPri, and buffer seq on try # "+str(num_tries))
                    #print(finalrecwithAssPrimers.seq)
            else:
                if (gsflag==0):#long homopolymer detected

                    print(caarec.id+ " failed, error: long stretch of homopolymer on try # "+str(num_tries))

                if (numMatch==0 and settings.input_type == 'AA'):

                    print(caarec.id+ " failed error: translation mismatch on try # "+str(num_tries))
                    print("caarec")
                    print(caarec.seq)
                    print("updatedseqrecord")
                    print(updatedseqrecord.translate().seq)
                else:
                    print(caarec.id+ " failed error: DNA numMatch==0 on try # "+str(num_tries))

            #this is the last code in the while loop for num_tries
            num_tries += 1
        #too many tries
        if (num_tries > settings.max_num_attempts-1):
            aarec_counter += 1
            bad_rec_counter += 1
            print("*** "+caarec.id+ " reached max number of attempts. Adding seq to badIDsarray.")
            #this seq cant be split, add it to list of files to remove
            if not (caarec.id in badIDsarray):
                badIDsarray.append(caarec.id)
                badIDfasta.append(caarec)
        
    ############ Protein Record Loop End ####################    
    
    print("Processed "+str(aarec_counter)+" records.")
    print(str(bad_rec_counter)+" records failed splitting.")
    print(str(len(final_records_oligos_extra_oligo))+" records need more oligos (extra_split files)")

    #write the split oligos here:
    fileout = open(settings.output_folder+'/'+settings.output_file_names.split('/')[-1].replace('fasta','oligos'),'w')
    for thisrec in final_records_oligos:
        fileout.write('>' + thisrec[0] + '\n')
        for currentoligo in thisrec[1]:
            #write split oligos:
            fileout.write(str(currentoligo) + '\n')
    fileout.close()
    
    #save file for each lib (codon 1) with full gene seq before splitting (with restriction sites (NdeI & KpnI) and assembly primers):
    gene_nt_fileout = open(settings.output_folder+'/'+settings.output_file_names.split('/')[-1].replace('fasta','full_wRE_wPrim.genes'),'w')
    SeqIO.write(final_records_gene_level, gene_nt_fileout, "fasta")
    gene_nt_fileout.close()
    
    #save file for each lib (codon 1) with full gene seq before splitting (with no restriction sites (NdeI & KpnI) and no assembly primers):
    #these seqs have start codon but no stop codon
    gene_nt_no_primer_or_RE_fileout = open(settings.output_folder+'/'+settings.output_file_names.split('/')[-1].replace('fasta','full_nRE_nPrim.genes'),'w')
    SeqIO.write(final_records_gene_level_no_primer_or_RE, gene_nt_no_primer_or_RE_fileout, "fasta")
    gene_nt_no_primer_or_RE_fileout.close()
    
    #save file for each lib (codon 1) with full gene seq before splitting (with restriction sites (NdeI & KpnI) and no assembly primers):
    #these seqs have start codon (part of NdeI) and stop codon (TAA) before KpnI
    gene_nt_no_primer_fileout = open(settings.output_folder+'/'+settings.output_file_names.split('/')[-1].replace('fasta','full_wRE_noPrim.genes'),'w')
    SeqIO.write(final_records_gene_level_no_primer, gene_nt_no_primer_fileout, "fasta")
    gene_nt_no_primer_fileout.close()
    
    if settings.input_type == 'AA':
        #save file for each lib (codon 1 & 2) with full protein seq (a.a.) (these have starting Met)
        protein_nt_fileout = open(settings.output_folder+'/'+settings.output_file_names.split('/')[-1].replace('fasta','proteins'),'w')
        SeqIO.write(final_records_protein_level, protein_nt_fileout, "fasta")
        protein_nt_fileout.close()
    
    if settings.input_type == 'DNA':
        #save file for each lib with input DNA seq
        dna_fileout = open(settings.output_folder+'/'+settings.output_file_names.split('/')[-1].replace('fasta','dna'),'w')
        SeqIO.write(final_records_dna_level, dna_fileout, "fasta")
        dna_fileout.close()

    #save for codon 2 also:
    if settings.gen_two_codon and settings.input_type == 'AA':

        fileout2 = open(settings.output_folder+'/'+settings.output_file_names.split('/')[-1].replace('fasta','codon2.oligos'),'w')
        for thisrec in final_records_oligos2:
            fileout2.write('>' + thisrec[0] + '\n')
            for currentoligo in thisrec[1]:
                #write split oligos:
                fileout2.write(str(currentoligo) + '\n')
        fileout2.close()

        #save file for each lib (codon 2) with full gene seq before splitting (with restriction sites (NdeI & KpnI) and assembly primers):
        gene_nt_fileout2 = open(settings.output_folder+'/'+settings.output_file_names.split('/')[-1].replace('fasta','codon2.full_wRE_wPrim.genes'),'w')
        SeqIO.write(final_records_gene_level2, gene_nt_fileout2, "fasta")
        gene_nt_fileout2.close()
        
        #save file for each lib (codon 2) with full gene seq before splitting (with no restriction sites (NdeI & KpnI) and no assembly primers):
        #these seqs have start codon but no stop codon
        gene_nt_no_primer_or_RE_fileout2 = open(settings.output_folder+'/'+settings.output_file_names.split('/')[-1].replace('fasta','codon2.full_nRE_nPrim.genes'),'w')
        SeqIO.write(final_records_gene_level_no_primer_or_RE2, gene_nt_no_primer_or_RE_fileout2, "fasta")
        gene_nt_no_primer_or_RE_fileout2.close()
        
        #save file for each lib (codon 2) with full gene seq before splitting (with restriction sites (NdeI & KpnI) and no assembly primers):
        #these seqs have start codon (part of NdeI) and stop codon (TAA) before KpnI
        gene_nt_no_primer_fileout2 = open(settings.output_folder+'/'+settings.output_file_names.split('/')[-1].replace('fasta','codon2.full_wRE_noPrim.genes'),'w')
        SeqIO.write(final_records_gene_level_no_primer2, gene_nt_no_primer_fileout2, "fasta")
        gene_nt_no_primer_fileout2.close()
        
        #save file for each lib (codon 1 & 2) with full protein seq (a.a.) (these have starting Met)
        #redundant
        protein_nt_fileout2 = open(settings.output_folder+'/'+settings.output_file_names.split('/')[-1].replace('fasta','codon2.proteins'),'w')
        SeqIO.write(final_records_protein_level2, protein_nt_fileout2, "fasta")
        protein_nt_fileout2.close()
        
        #save file for each lib with nt identity between two codon versions
        proteinnt_ident_fileout = settings.output_folder+'/'+settings.output_file_names.split('/')[-1].replace('fasta','lib_ident.csv')
        resultFile = open(proteinnt_ident_fileout,'w')
        wr = csv.writer(resultFile, dialect='excel')
        wr.writerow(pct_ident_array)

    #in the case when the split is one more than expected write the split oligos here:
    if write_extra_oligo_files:
        fileout_extra_oligo = open(settings.output_folder+'/'+settings.output_file_names.split('/')[-1].replace('fasta','extra_split.oligos'),'w')
        for thisrec in final_records_oligos_extra_oligo:
            fileout_extra_oligo.write('>' + thisrec[0] + '\n')
            for currentoligo in thisrec[1]:
                #write split oligos:
                fileout_extra_oligo.write(str(currentoligo) + '\n')
        fileout_extra_oligo.close()
        if settings.gen_two_codon and settings.input_type == 'AA':
            fileout_extra_oligo2 = open(settings.output_folder+'/'+settings.output_file_names.split('/')[-1].replace('fasta','extra_split.codon2.oligos'),'w')
            for thisrec in final_records_oligos_extra_oligo2:
                fileout_extra_oligo2.write('>' + thisrec[0] + '\n')
                for currentoligo in thisrec[1]:
                    #write split oligos:
                    fileout_extra_oligo2.write(str(currentoligo) + '\n')
            fileout_extra_oligo2.close()

        #save file for each lib (codon 1) with full gene seq before splitting (with restriction sites (NdeI & KpnI) and assembly primers):
        gene_nt_fileout_extra_oligo = open(settings.output_folder+'/'+settings.output_file_names.split('/')[-1].replace('fasta','full_wRE_wPrim_extra_split.genes'),'w')
        SeqIO.write(final_records_gene_level_extra_oligo, gene_nt_fileout_extra_oligo, "fasta")
        gene_nt_fileout_extra_oligo.close()
        
        #save file for each lib (codon 1) with full gene seq before splitting (with no restriction sites (NdeI & KpnI) and no assembly primers):
        #these seqs have start codon but no stop codon
        gene_nt_no_primer_or_RE_fileout_extra_oligo = open(settings.output_folder+'/'+settings.output_file_names.split('/')[-1].replace('fasta','full_nRE_nPrim_extra_split.genes'),'w')
        SeqIO.write(final_records_gene_level_no_primer_or_RE_extra_oligo, gene_nt_no_primer_or_RE_fileout_extra_oligo, "fasta")
        gene_nt_no_primer_or_RE_fileout_extra_oligo.close()
        
        #save file for each lib (codon 1) with full gene seq before splitting (with restriction sites (NdeI & KpnI) and no assembly primers):
        #these seqs have start codon (part of NdeI) and stop codon (TAA) before KpnI
        gene_nt_no_primer_fileout_extra_oligo = open(settings.output_folder+'/'+settings.output_file_names.split('/')[-1].replace('fasta','full_wRE_noPrim_extra_split.genes'),'w')
        SeqIO.write(final_records_gene_level_no_primer_extra_oligo, gene_nt_no_primer_fileout_extra_oligo, "fasta")
        gene_nt_no_primer_fileout_extra_oligo.close()
        
        if settings.input_type == 'AA':
            #save file for each lib (codon 1 & 2) with full protein seq (a.a.) (these have starting Met)
            protein_nt_fileout_extra_oligo = open(settings.output_folder+'/'+settings.output_file_names.split('/')[-1].replace('fasta','extra_split.proteins'),'w')
            SeqIO.write(final_records_protein_level_extra_oligo, protein_nt_fileout_extra_oligo, "fasta")
            protein_nt_fileout_extra_oligo.close()
        
        if settings.input_type == 'DNA':
            #save file for each lib with input DNA seq
            dna_fileout_extra_oligo = open(settings.output_folder+'/'+settings.output_file_names.split('/')[-1].replace('fasta','extra_split.dna'),'w')
            SeqIO.write(final_records_dna_level_extra_oligo, dna_fileout_extra_oligo, "fasta")
            dna_fileout_extra_oligo.close()

        #save for codon 2 also:
        if settings.gen_two_codon and settings.input_type == 'AA':
            #save file for each lib (codon 2) with full gene seq before splitting (with restriction sites (NdeI & KpnI) and assembly primers):
            gene_nt_fileout_extra_oligo2 = open(settings.output_folder+'/'+settings.output_file_names.split('/')[-1].replace('fasta','codon2.full_wRE_wPrim_extra_split.genes'),'w')
            SeqIO.write(final_records_gene_level_extra_oligo2, gene_nt_fileout_extra_oligo2, "fasta")
            gene_nt_fileout_extra_oligo2.close()
            
            #save file for each lib (codon 2) with full gene seq before splitting (with no restriction sites (NdeI & KpnI) and no assembly primers):
            #these seqs have start codon but no stop codon
            gene_nt_no_primer_or_RE_fileout_extra_oligo2 = open(settings.output_folder+'/'+settings.output_file_names.split('/')[-1].replace('fasta','codon2.full_nRE_nPrim_extra_split.genes'),'w')
            SeqIO.write(final_records_gene_level_no_primer_or_RE_extra_oligo2, gene_nt_no_primer_or_RE_fileout_extra_oligo2, "fasta")
            gene_nt_no_primer_or_RE_fileout_extra_oligo2.close()
            
            #save file for each lib (codon 2) with full gene seq before splitting (with restriction sites (NdeI & KpnI) and no assembly primers):
            #these seqs have start codon (part of NdeI) and stop codon (TAA) before KpnI
            gene_nt_no_primer_fileout_extra_oligo2 = open(settings.output_folder+'/'+settings.output_file_names.split('/')[-1].replace('fasta','codon2.full_wRE_noPrim_extra_split.genes'),'w')
            SeqIO.write(final_records_gene_level_no_primer_extra_oligo2, gene_nt_no_primer_fileout_extra_oligo2, "fasta")
            gene_nt_no_primer_fileout_extra_oligo2.close()
            
            #save file for each lib (codon 1 & 2) with full protein seq (a.a.) (these have starting Met)
            #redundant
            protein_nt_fileout_extra_oligo2 = open(settings.output_folder+'/'+settings.output_file_names.split('/')[-1].replace('fasta','codon2.extra_split.proteins'),'w')
            SeqIO.write(final_records_protein_level_extra_oligo2, protein_nt_fileout_extra_oligo2, "fasta")
            protein_nt_fileout_extra_oligo2.close()
            

    if badIDsarray:
        #save all of the IDs for sequences that failed completely:
        resultFile2 = open(settings.IDs_to_remove,'w')
        wr = csv.writer(resultFile2, dialect='excel')
        wr.writerow(badIDsarray)

    if badIDfasta:
        #make a FASTA file with them
        bad_fasta_fileout = open(settings.output_folder+'/'+settings.output_file_names.split('/')[-1].replace('fasta','fail_split.fasta'),'w')
        SeqIO.write(badIDfasta, bad_fasta_fileout, "fasta")
        bad_fasta_fileout.close()

    
    print("--- %s seconds ---" % (time.time() - start_time))
    