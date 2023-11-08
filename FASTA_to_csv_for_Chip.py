from Bio import SeqIO
from Bio import Entrez
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Restriction import *
import csv

##################################
#FUNCTIONS:
def getOligos(filename):
    constructs = []
    currentconstruct = 'foo'
    for seqrec in SeqIO.parse(filename, "fasta"):
        #print(currentconstruct)
        if not seqrec.id.startswith(currentconstruct):
            if seqrec.id.count(';') == 0:
                currentconstruct = seqrec.id
            else:
                currentconstruct = seqrec.id[0:seqrec.id.rfind(';')]
            constructs.append([currentconstruct])
            constructs[-1].append(seqrec)
        else:
            constructs[-1].append(seqrec)
    return constructs

##################################
#INPUTS:
inputfiles = ['db_libs/HK4_He_out-finaloligos.fasta','db_libs/HK4_He_out-finaloligos.codon2.fasta',
              'db_libs/HK5_He_out-finaloligos.fasta','db_libs/HK5_He_out-finaloligos.codon2.fasta']


file_out_1 = "chip_out/300mer_chip8.csv"
file_out_2 = "chip_out/300mer_chip8_raw.csv"

##################################

csvfile = open(file_out_1, 'w')
fieldnames = ['Lib#','Construct#','Oligo#','Barcode','Accession','Source','Definition','Sequence']
writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
writer.writeheader()

for fileindex in range(len(inputfiles)):
    constructs = getOligos(inputfiles[fileindex].replace('.oligos','-finaloligos.fasta'))
    first_oligo_counter = 0
    middle_oligo_counter = 0
    last_oligo_counter = 0
    for index in range(len(constructs)):
        oligos = constructs[index][1:]
        oligo_counter=1
        #temp_str = constructs[index][0]
        temp_str = constructs[index][0].split(';')
        seqid = temp_str[0]
        #print(temp_str)
        #print(len(oligos))
        if len(oligos) == 1:
            barcode_id = "None"
        else:
            barcode_id = temp_str[1]
        for oligoseqrec in oligos:
            if oligo_counter == 1:
                writer.writerow({'Lib#':str(fileindex+1),'Construct#':str(index+1),'Oligo#':str(oligo_counter),'Barcode':barcode_id,'Accession': seqid, 'Sequence': oligoseqrec.seq})
            else:
                writer.writerow({'Lib#':'', 'Construct#':'','Oligo#':str(oligo_counter),'Barcode': '','Accession':'', 'Sequence': oligoseqrec.seq})
            oligo_counter += 1

csvfile.close()

csvfile = open(file_out_2, 'w')
fieldnames = ['Name','Sequence']
writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
#writer.writeheader()

for fileindex in range(len(inputfiles)):
    constructs = getOligos(inputfiles[fileindex].replace('.oligos','-finaloligos.fasta'))
    first_oligo_counter = 0
    middle_oligo_counter = 0
    last_oligo_counter = 0
    for index in range(len(constructs)):
        oligos = constructs[index][1:]
        oligo_counter=1
        for oligoseqrec in oligos:
            temp_str = str(fileindex+1)+'_'+constructs[index][0]+';'+str(oligo_counter)
            writer.writerow({'Name':temp_str, 'Sequence': oligoseqrec.seq})
            oligo_counter += 1

csvfile.close()