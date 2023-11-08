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
        #print seqrecord.id
        #print len(seqrecord)
    #print len(records)
    return records

def plotOligo_length(constructs,filename):
    #plot histogram of lengths
    construct_lengths = []

    for construct in constructs:
    	construct_lengths.append(len(construct.seq))

    #    for i in range(1,len(construct.seq)):
    #        construct_lengths.append(len(construct[i]))
    generate_histogram = True
    if generate_histogram:
        data = construct_lengths
        print("min length: " +str(min(data)) + " max length: " +str(max(data)))
        pylab.hist(data, bins=(max(data)-min(data)))
        pylab.title("%i gene length distribution\nfrom %i to %i" \
                    % (len(data),min(data),max(data)))
        pylab.xlabel("Gene length (nt)")
        pylab.ylabel("Count")
        pylab.savefig(filename+'.png', bbox_inches='tight')
        pylab.savefig(filename+'.pdf', bbox_inches='tight')
        pylab.show()

filenames = ['db_oligo/DHFR_Lib01_3oligo.full_wRE_wPrim.genes', 'db_oligo/DHFR_Lib01_3oligo.codon2.full_wRE_wPrim.genes',
             'db_oligo/DHFR_Lib02_4oligo.full_wRE_wPrim.genes', 'db_oligo/DHFR_Lib02_4oligo.codon2.full_wRE_wPrim.genes',
             'db_oligo/DHFR_Lib03_4oligo.full_wRE_wPrim.genes', 'db_oligo/DHFR_Lib03_4oligo.codon2.full_wRE_wPrim.genes']

for index in range(len(filenames)):
    print("Processing Lib"+str(index+1))
    buildoligos = obtainFastaSequences(filenames[index])
    plotOligo_length(buildoligos,"plots/Oligo_Lib"+str(index+1)+"_length_in_hist")