import Bio
from Bio import Seq
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

def getFastaSeqs(filename):
    fastqseqs = []
    handle = open(filename)
    for seqrec in SeqIO.parse(handle,"fasta"):
        fastqseqs.append(seqrec)
    handle.close()
    return fastqseqs

#parse through the results of the BLAT alignments and see which amp primers are best
#generate primer files based on this

outF_file = "ampprimers-skpp15/skpp15-forward_select_mod_Chip6.faa"
outR_file = "ampprimers-skpp15/skpp15-reverse_select_mod_Chip6.faa"
skpp_15_F_file = "ampprimers-skpp15/skpp15-forward.faa"
skpp_15_R_file = "ampprimers-skpp15/skpp15-reverse.faa"

print_all_even_below = False

#output file of BLAT alignment of first 200 skpp-15 primer pairs against all oligos in library (without amp primer regions):
hits_file = 'BLAT_psl/ampAll_200_hits.psl'
primerhits = dict()
numhitsbig = dict()

#primers that contain bad RE sites and other stuff
bad_primer_IDs = [0,5,8,14,47,52,61,66,73,81,82,87,90,107,109,131,134,137,141,144,173,177,180,190]

cutoff_for_hit_size = 12 #the most # of matching bases allowed (big hits)
selection_hitsizenum = 1 #how many hits above cutoff to allow
totalhitscutoff = 24 #the most total number of hits both F and R sum allowed

#open alignments file
with open(hits_file) as fp:
    #for each alignment
    for line in fp:
        #split into output values
        tmpline = line.split("\t")
        #primer reference id is here (eg: skpp15-1-R):
        primer = tmpline[9]
        #number of matches in alignment:
        hitsize = int(tmpline[0])
        #print(tmpline)
        #print(primer)
        temp = primer.split("-")
        #print(temp)
        #the primer name (eg 1F 1R 2F 2R etc...):
        primer = temp[1]+temp[2]

        #print(primer+" "+str(hitsize))

        #have we seen alignments for this primer before?
        if primer in primerhits:
            #count total alignments:
            primerhits[primer] += 1
            if hitsize > cutoff_for_hit_size:
                #count total alignments above cutoff:
                numhitsbig[primer] += 1
        else:
            primerhits[primer] = 1
            if hitsize > cutoff_for_hit_size:
                numhitsbig[primer] = 1
            else:
                numhitsbig[primer] = 0

runcounter = 0
numselected = 0
num_to_output = 30

skpp15F = getFastaSeqs(skpp_15_F_file)
skpp15R = getFastaSeqs(skpp_15_R_file)

skpp15Fmod = []
skpp15Rmod = []

print("index FWDhits REVhits totalBIGhits totalHits nameF nameR")
for index in range(1,200):
    primernameF = str(index)+'F'
    if primernameF in primerhits:
        primerhitsFcount=primerhits[primernameF]
        hitscountF = numhitsbig[primernameF]
    else:
        primerhitsFcount=0
        hitscountF=0
    primernameR = str(index)+'R'
    if primernameR in primerhits:
        primerhitsRcount=primerhits[primernameR]
        hitscountR = numhitsbig[primernameR]
    else:
        primerhitsRcount=0
        hitscountR=0
    totalBIGhits = hitscountF+hitscountR
    totalhits = primerhitsFcount+primerhitsRcount
    
    if print_all_even_below:
        print(str(index)+" "+str(primerhitsFcount)+" "+str(primerhitsRcount)+" "+str(totalBIGhits)+" "+str(totalhits)+" "+primernameF+" "+primernameR)

    #print only best ones
    if (totalBIGhits < selection_hitsizenum+1) and (totalhits < totalhitscutoff) and (index not in bad_primer_IDs):
        print(str(index)+" "+str(primerhitsFcount)+" "+str(primerhitsRcount)+" "+str(totalBIGhits)+" "+str(totalhits)+" "+primernameF+" "+primernameR)
        numselected += 1
        runcounter += totalhits
        skpp15Fmod.append(skpp15F[index-1])
        skpp15Rmod.append(skpp15R[index-1])

print("\n")
print("total_hits: "+str(runcounter))
print("numselected: "+str(numselected))

output_handleF = open(outF_file, "w")
SeqIO.write(skpp15Fmod, output_handleF, "fasta")
output_handleF.close()

output_handleR = open(outR_file, "w")
SeqIO.write(skpp15Rmod, output_handleR, "fasta")
output_handleR.close()

for index in range(len(skpp15Fmod)):
    print(skpp15Fmod[index].id +" "+ skpp15Fmod[index].seq +" "+ skpp15Rmod[index].id +" "+ skpp15Rmod[index].seq)

