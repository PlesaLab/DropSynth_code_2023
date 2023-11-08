import os
os.chdir('/Users/Angus/Documents/Chip-4/DropSynth-master-final')

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.Alphabet import generic_dna
from Bio.Restriction import *
import pylab

def obtainFastaSequences(filename):
    handle = open(filename)
    records = []
    for seqrecord in SeqIO.parse(handle, "fasta"):
        records.append(seqrecord)
        #print(seqrecord.id)
        #print(len(seqrecord))
    #print(len(records))
    return records
        
barcodes = obtainFastaSequences('barcodes/filt_prim_12nt_Lev_3_Tm_38_44_GC_45_55_SD_2.fasta')

BColigolist1 = []
for i in range(0,384):
            BColigo = barcodes[i].reverse_complement() + Seq("AGGCATCTAATAGCACGTGT",generic_dna)
            BColigolist1.append(BColigo)
BColigolist2 = []
for i in range(384,768):
            BColigo = barcodes[i].reverse_complement() + Seq("AGGCATCTAATAGCACGTGT",generic_dna)
            BColigolist2.append(BColigo)  
BColigolist3 = []
for i in range(768,1152):
            BColigo = barcodes[i].reverse_complement() + Seq("AGGCATCTAATAGCACGTGT",generic_dna)
            BColigolist3.append(BColigo) 
BColigolist4 = []
for i in range(1152,1536):
            BColigo = barcodes[i].reverse_complement() + Seq("AGGCATCTAATAGCACGTGT",generic_dna)
            BColigolist4.append(BColigo) 
          
filename = "BColigolist1.xls"
fileout = open(filename,'w')
for i in range(0,384):
    fileout.write("bcoligo-cp12mer-"+str(i+1)+"\t")
    fileout.write(str(BColigolist1[i].seq)+"\n")
fileout.close()

filename = "BColigolist2.xls"
fileout = open(filename,'w')
for i in range(0,384):
    fileout.write("bcoligo-cp12mer-"+str(i+385)+"\t")
    fileout.write(str(BColigolist2[i].seq)+"\n")
fileout.close()

filename = "BColigolist3.xls"
fileout = open(filename,'w')
for i in range(0,384):
    fileout.write("bcoligo-cp12mer-"+str(i+769)+"\t")
    fileout.write(str(BColigolist3[i].seq)+"\n")
fileout.close()

filename = "BColigolist4.xls"
fileout = open(filename,'w')
for i in range(0,384):
    fileout.write("bcoligo-cp12mer-"+str(i+1153)+"\t")
    fileout.write(str(BColigolist4[i].seq)+"\n")
fileout.close()