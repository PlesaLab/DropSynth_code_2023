##
## DO NOT COPY EXTERAL CODE INTO HERE
##
## This must be completely independent to avoid propagation of errors.
##

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Restriction import *

##################################
#FUNCTIONS:
def obtainFastaSequences(filename):
    handle = open(filename)
    records = []
    for seqrecord in SeqIO.parse(handle, "fasta"):
        records.append(seqrecord)
        #print(seqrecord.id)
        #print(len(seqrecord))
    #print(len(records))
    return records

def getOligos(filename):
    print("opening: "+ filename)
    constructs = []
    currentconstruct = 'foo'
    for seqrec in SeqIO.parse(filename, "fasta"):
        if not seqrec.id.startswith(currentconstruct):
            currentconstruct = seqrec.id[0:seqrec.id.rfind(';')]
            constructs.append([currentconstruct])
            constructs[-1].append(seqrec)
        else:
            constructs[-1].append(seqrec)
    return constructs

# Reversing a list using reversed()
def Reverse(lst):
    return [ele for ele in reversed(lst)]

def checkConstruct(construct,lengthmax,filenum,ampprimf,ampprimr,barcode,assemblyprimf,assemblyprimr,padding_between_btsaI_ampR,allow_KpnI_in_middle,cloning_start,cloning_end):
    #print(construct)
    
    oligos = construct[1:]
    count_first = 0
    count_mid = 0
    count_last = 0
    BspQI_length = 8
    BtsaI_length = 6
    amppri_length = len(ampprimf)
    asmpri_length = len(assemblyprimf)
    bc_length = len(barcode)
    btsai_from_end_pos = amppri_length + BtsaI_length - 1
    bspqi_first = amppri_length + BspQI_length + 1
    bspqi_second = amppri_length + BspQI_length + bc_length -2
    asmpriFsearch_min_pos = amppri_length + 2*BspQI_length + bc_length +BtsaI_length-1
    asmpriRsearch_pos = lengthmax - amppri_length - BtsaI_length - asmpri_length
    
    rb = RestrictionBatch([BtsI, BspQI, BsaI, EcoRI])

    KpnI_min_from_end = lengthmax - amppri_length - BtsaI_length - asmpri_length
    KpnI_min_from_start = len(ampprimf)+ 2*BspQI_length + bc_length + BtsaI_length
    
    #check that lengths all pass
    for oligoseqrec in oligos:
        if len(oligoseqrec.seq)>lengthmax:
            print('Oligo length out of range:')
            print(oligoseqrec.id)
            print(str(oligoseqrec.seq) + '\t' + str(len(oligoseqrec.seq)))
    
    #check that amplification primers are correct
    for oligoseqrec in oligos:
        if ampprimf != str(oligoseqrec.seq[0:len(ampprimf)]) or ampprimr != str(oligoseqrec.seq[-1*len(ampprimr):].reverse_complement()):
            print('Amp primers not found:')
            print(oligoseqrec.id)
            print(str(oligoseqrec.seq))
    
    #check that barcode and BspQI sites are correct
    for oligoseqrec in oligos:
        #first site is 15(amp)+8(BspQI)+1=24
        #second site is 15(amp)+8(BspQI)+12(barcode)-2
        seqsearch = rb.search(oligoseqrec.seq)

        if len(seqsearch[BspQI])!=2:
            print("\naBad number of BspQI restriction sites")
            print(oligoseqrec.id)
            print(oligoseqrec.seq)
            print(BspQI.search(oligoseqrec.seq))
        
        if BspQI.search(oligoseqrec.seq) != [bspqi_first, bspqi_second]:
            print("BspQI SITES ARE WRONG")
            print(oligoseqrec.id)
            print(oligoseqrec.seq)
            print(BspQI.search(oligoseqrec.seq))
        #barcode site is 15amp + 8(BspQI)
        barcode_pos = amppri_length + BspQI_length
        if str(oligoseqrec.seq).find(barcode) != barcode_pos:
            print("BARCODE NOT FOUND")
            print("Barcode is "+oligoseqrec.seq[23:35]+" but expected "+barcode)
            print(oligoseqrec.id)
            print(oligoseqrec.seq)
            print(barcode)
    
    #check that BtsI sites are correct
    for oligoseqrec in oligos:
        seqsearch = rb.search(oligoseqrec.seq)

        btssearch = BtsI.search(oligoseqrec.seq)
        
        
        asmpriFsearch = str(oligoseqrec.seq).find(assemblyprimf)
        assemblyprimr_seq = Seq(assemblyprimr)
        asmpriRsearch = str(oligoseqrec.seq).find(str(assemblyprimr_seq.reverse_complement()))

        startREsearch = oligoseqrec.seq.find(cloning_start, asmpriFsearch)

        endREsearch = oligoseqrec.seq.find(cloning_end, asmpriFsearch)

        #pos end btsaI from end = amp length + BtsaI_length -1
        if len(btssearch) !=2 or len(oligoseqrec)-btssearch[1] != btsai_from_end_pos: # btssearch[0] != 43:#
            print("End BtsI site bad")
            print(oligoseqrec.id)
            print(oligoseqrec.seq)
            print(barcode)
            print(str(btssearch))
        
        first_oligo = 0
        if asmpriFsearch != -1:
            #this should be the first oligo in assembly
            first_oligo = 1
            count_first += 1
            #min asm search = amplength + 2*bspqilength + bclength+btsai-1
            if asmpriFsearch < asmpriFsearch_min_pos or asmpriFsearch > 100:
                print("Frist oligo Assembly FWD primer wrong")
                print(oligoseqrec.id)
                print(oligoseqrec.seq)
            

            if startREsearch == -1 or startREsearch != asmpriFsearch + asmpri_length:
                print("First oligo CloningStart site bad")
                print(oligoseqrec.id)
                print(oligoseqrec.seq)
                print(startREsearch)
                #print(asmpriFsearch)
                #print(asmpri_length)
            
            if endREsearch != -1:
                if not (endREsearch < asmpriRsearch_pos):
                    print("First oligo has a CloningEnd site")
                    print(oligoseqrec.id)
                    print(oligoseqrec.seq)
                    print(endREsearch)
            
            if btssearch[0] != asmpriFsearch + 3:
                print("First oligo BtsI site bad")
                print(oligoseqrec.id)
                print(oligoseqrec.seq)
        
        last_oligo = 0
        
        if asmpriRsearch != -1: #this should be the last oligo in assembly
            
            last_oligo = 1
            count_last += 1
            
            if first_oligo == 1:
                print("This oligo contains both FWD and REV assembly primers")
                print(oligoseqrec.id)
                print(oligoseqrec.seq)
            
            if asmpriRsearch != asmpriRsearch_pos:#200mer: 159 #230mer:189 #calc: total length - 15ampPriLength - 6BtsaI -20asmPriLength
                print("Last oligo assembly primer REV wrong")
                print(oligoseqrec.id)
                print(oligoseqrec.seq)
            
            if startREsearch != -1:
                if not (startREsearch < asmpriRsearch_pos):
                    print("Last oligo has an CloningStart site")
                    print(oligoseqrec.id)
                    print(oligoseqrec.seq)
            
            #disabled because not always true here
            #if len(endREsearch) != 1 or endREsearch[0] > KpnI_min_from_end or endREsearch[0] < KpnI_min_from_start:
            #    print("Last oligo KpnI site bad")
            #    print(oligoseqrec.id)
            #    print(oligoseqrec.seq)
            
        elif first_oligo == 0: # this is middle oligo
            
            count_mid += 1
            if startREsearch != -1:
                if not (startREsearch < asmpriRsearch_pos):
                    print("Middle oligo has an CloningStart site")
                    print(oligoseqrec.id)
                    print(oligoseqrec.seq)
            
            #if len(endREsearch) != 0 and not allow_KpnI_in_middle:
            #    print("Middle oligo has a KpnI site. Lib_num:"+str(filenum))
            #    print(oligoseqrec.id)
            #    print(oligoseqrec.seq)
                
            if asmpriFsearch != -1:
                print("Middle oligo has an Assembly primer FWD site")
                print(oligoseqrec.id)
                print(oligoseqrec.seq)
            	
            if asmpriRsearch != -1:
                print("Middle oligo has an Assembly primer REV site")
                print(oligoseqrec.id)
                print(oligoseqrec.seq)
                
    return count_first, count_mid, count_last


def checkTranslation(construct,dprotc,assemblyprimf,assemblyprimr,cloning_start,cloning_end):
    #print(construct)

    

    error_flag=0
    maxoverlap=50
    oligos = construct[1:]
    proteinname=construct[0].split(";")[0]
    #oligos_between_btsI=[]
    just_the_gene=[]
    cloning_start_skip=len(cloning_start)
    assemblyprimr_seq = Seq(assemblyprimr)
    kpni_in_middle_flag=0
    for oligoseqrec in oligos:
        seqstr=str(oligoseqrec.seq)
        #print("original_full "+seqstr)
        startpos=seqstr.find('GCAGTG')
        endpos=seqstr.find('CACTGC',startpos+6,len(seqstr))
        #print(str(startpos)+" "+str(endpos))
        temp_str=seqstr[startpos+6:endpos]
        #oligos_between_btsI.append(temp_str)
        #print(temp_str)
        asmpriFsearch = str(temp_str).find(assemblyprimf)
        asmpriRsearch = str(temp_str).find(str(assemblyprimr_seq.reverse_complement()))
        if asmpriFsearch != -1:
            #this is first oligo
            temp_str_skip_clnstrt=temp_str[len(assemblyprimf)+cloning_start_skip:]
            just_the_gene.append(temp_str_skip_clnstrt)
        elif asmpriRsearch == -1:
            #middle oligo

            just_the_gene.append(temp_str)#note not the same as above, nothing skipped
        elif asmpriRsearch != -1 and kpni_in_middle_flag == 0:
            #this is last oligo
            temp_str_skip_endPrim=temp_str[:len(temp_str)-len(assemblyprimr)]
            just_the_gene.append(temp_str_skip_endPrim)
    
    num_fragments=len(just_the_gene)

    if num_fragments>1:
        geneseq=just_the_gene[0]
        #print("firstfrag: "+geneseq)
        for i in range(1,len(just_the_gene)):
            #print(i)
            overlapseq=LongestCommonSubstring(just_the_gene[i-1],just_the_gene[i],maxoverlap)
            overlen=len(overlapseq)
            if overlen<5 or overlen>maxoverlap:
                print("Error in "+construct[0])
                print(str(overlen)+" "+overlapseq)
                print("In1: "+just_the_gene[i-1])
                print("In2: "+just_the_gene[i])
                print(just_the_gene)
                for oligoseqrec in oligos:
                    seqstr=str(oligoseqrec.seq)
                    print("original_full "+seqstr)
                    #print("original_full "+seqstr)
                    startpos=seqstr.find('GCAGTG')
                    endpos=seqstr.find('CACTGC',startpos+6,len(seqstr))
                    #print(str(startpos)+" "+str(endpos))
                    temp_str=seqstr[startpos+6:endpos]
                    print("start:"+str(startpos)+" end: "+str(endpos))
                error_flag=1
            else:
                geneseq=geneseq+just_the_gene[i][overlen:]
                #print("overlap: "+just_the_gene[i][:overlen])
                #print("no overlap: "+just_the_gene[i][overlen:])
                #print(str(i+1)+" frag: "+geneseq)

        geneseq_until_cloning_end = Seq(geneseq.partition(cloning_end)[0])
        if geneseq_until_cloning_end.translate(to_stop=True) != dprotc[proteinname]:
            print("BAD translation!!!!!! "+construct[0])
            #print(construct)
            #print("got dna: "+geneseq_until_cloning_end)
            #print("got dna no partition: "+Seq(geneseq))
            #print("got prot: "+Seq(geneseq.partition(cloning_end)[0]).translate())
            #print("got prot no partition: "+Seq(geneseq).translate())

            print("got prot: "+geneseq_until_cloning_end.translate())
            print("expected: "+dprotc[proteinname])
            print("last overlap: "+overlapseq)
            for oligoseqrec in oligos:
                seqstr=str(oligoseqrec.seq)
                #print("original_full "+seqstr)
            print(just_the_gene)
            error_flag=1

        #check for illegal sites in just gene
        rb = RestrictionBatch([EcoRI, BspQI])
        REgenesearch = rb.search(Seq(geneseq))
        if len(REgenesearch[EcoRI])!= 0 or len(REgenesearch[BspQI])!=0:# or len(REgenesearch[BsaI])!=0 or len(REgenesearch[NdeI])!=0 or len(REgenesearch[KpnI])!=0 or :
            print("\nBad number of restriction sites in assembled gene (not including primers)")
            print(construct[0])
            print(geneseq)
            print(REgenesearch)



    elif num_fragments==1:
        if just_the_gene[0].partition(cloning_end)[0].translate(to_stop=True) != dprotc[proteinname]:
            print("BAD translation!!!!!! "+construct[0])
            print("got: "+just_the_gene[0].partition(cloning_end)[0].translate(to_stop=True))
            print("expected: "+dprotc[proteinname])
            #for i in range(1,len(just_the_gene)):
            #    LongestCommonSubstringprint(S1, S2, maxoverlap)
            error_flag=1
    else:
        print("BAD!!!!!!")
        error_flag=1



    #print(just_the_gene)
    return error_flag

def LongestCommonSubstring(S1, S2, maxoverlap):
    overlap_tracker=0
    #print("inputs1: "+S1)
    #print("inputs2: "+S2)
    for i in range(5,maxoverlap):
        #print(i)
        #print(S1[len(S1)-i:])
        #print(S2[:i])
        if S1[len(S1)-i:]==S2[:i]:
            overlap_tracker=i
            #print("match!" +str(i))

            #print("match: "+str(i)+" in "+S1+" "+S2)
    
    if overlap_tracker==0:
        print("error in overlap trecker")
        S2=""
        # for i in range(5,maxoverlap):
        #     print(i)
        #     print(S1[len(S1)-i:])
        #     print(S2[:i])
        #     if S1[len(S1)-i:]==S2[:i]:
        #         overlap_tracker=i
        #         print("match!" +str(i))
    else:
        S2=S2[0:overlap_tracker]

    #print("overlap: "+S2)
    return S2

def LongestCommonSubstringprint(S1, S2, maxoverlap):
    overlap_tracker=0
    print("inputs1: "+S1)
    print("inputs2: "+S2)
    for i in range(5,maxoverlap):
        print(i)
        print(S1[len(S1)-i:])
        print(S2[:i])
        if S1[len(S1)-i:]==S2[:i]:
            overlap_tracker=i
            print("match!" +str(i))

            #print("match: "+str(i)+" in "+S1+" "+S2)
    
    if overlap_tracker==0:
        print("error in overlap trecker")
    else:
        S2=S2[0:overlap_tracker]

    #print("overlap: "+S2)
    return S2

    # S1=S1[len(S1)-maxoverlap:]
    # S2=S2[0:maxoverlap]
    # M = [[0]*(1+len(S2)) for i in range(1+len(S1))]
    # longest, x_longest = 0, 0
    # for x in range(1,1+len(S1)):
    #     for y in range(1,1+len(S2)):
    #         if S1[x-1] == S2[y-1]:
    #             M[x][y] = M[x-1][y-1] + 1
    #             if M[x][y]>longest:
    #                 longest = M[x][y]
    #                 x_longest  = x
    #         else:
    #             M[x][y] = 0
    # return S1[x_longest-longest: x_longest]


#####################################
########## OPTIONS ##################
#####################################

inputfiles = ['db_libs/HK4_He_out-finaloligos.fasta','db_libs/HK4_He_out-finaloligos.codon2.fasta',
              'db_libs/HK5_He_out-finaloligos.fasta','db_libs/HK5_He_out-finaloligos.codon2.fasta']


#asmF skpp20 from 00_primer_screen.py output #skpp504F
#Primers for alternate codon versions are offset by len(num_oligos) = 15
assemblyprimf = ['ATCGGGGATGGTAACTAACG','ATCGGGGATGGTAACTAACG','ATCGGGGATGGTAACTAACG','ATCGGGGATGGTAACTAACG',
                 'ATCGGGGATGGTAACTAACG','ATCGGGGATGGTAACTAACG','ATCGGGGATGGTAACTAACG','ATCGGGGATGGTAACTAACG']
#enter the reverse primers #skpp504R-rc
#Primers for alternate codon versions are offset by len(num_oligos) = 15
assemblyprimr = ['ATAGCTGATTGTCCGTTGGT', 'ATAGCTGATTGTCCGTTGGT','ATAGCTGATTGTCCGTTGGT','ATAGCTGATTGTCCGTTGGT',
                 'ATAGCTGATTGTCCGTTGGT', 'ATAGCTGATTGTCCGTTGGT','ATAGCTGATTGTCCGTTGGT','ATAGCTGATTGTCCGTTGGT']

#ampF from 00_primer_screen.py output: skpp15 5##F for amplification primers
ampprimersf =  ['CGATCGTGCCCACCT', 'GGGTTCGAGCGGGAG', 
                'GCGGCACCACAAACT', 'GACTGCGGCGTTGGT']

#ampR (not RC) from 00_primer_screen.py output:
ampprimersr =  ['GTGCGGGCTCCAACT', 'TAGCGCGCAGAGAGG', 
                'CGTGGCCTCTGTCCT', 'TACGCCCGGGACAGA']

barcodes = obtainFastaSequences('barcodes/filt_prim_12nt_Lev_3_Tm_38_44_GC_45_55_SD_2_trim.fasta')
barcodes_revorder=Reverse(barcodes)

cloning_start = 'GGTCTCACAT'
cloning_end = 'GACGTGAGACC'

oligos_per_construct = [4,4,5,5]
constructs_per_lib = [1536, 1536, 1525, 1525]
barcode_reversed = [False, True, False, True]

allow_KpnI_in_middle = True
padding_between_btsaI_ampR = True

oligo_length = 300

input_proteins = []
for fnm in inputfiles:
    input_proteins.append(fnm.replace('-finaloligos.fasta','.proteins').replace('-finaloligos.codon2.fasta','.proteins'))

#####################################
######### / OPTIONS #################
#####################################

total_translation_errors = 0

for fileindex in range(len(inputfiles)):
    constructs = getOligos(inputfiles[fileindex])
    proteins = obtainFastaSequences(input_proteins[fileindex].replace('db_libs','db_oligo'))
    dprotc=dict()
    for prot in proteins:
        dprotc[prot.id]=prot.seq

    first_oligo_counter = 0
    middle_oligo_counter = 0
    last_oligo_counter = 0
    translate_error_counter = 0
    for index in range(len(constructs)):
        if barcode_reversed[fileindex] == True:
            this_bc = str(barcodes_revorder[index].seq)

        else:
            this_bc = str(barcodes[index].seq)

        add_first, add_mid, add_last = checkConstruct(constructs[index],oligo_length,fileindex+1,ampprimersf[fileindex],ampprimersr[fileindex],this_bc,assemblyprimf[fileindex],assemblyprimr[fileindex],padding_between_btsaI_ampR,allow_KpnI_in_middle,cloning_start,cloning_end)
        translate_error = checkTranslation(constructs[index],dprotc,assemblyprimf[fileindex],assemblyprimr[fileindex],cloning_start,cloning_end)
        translate_error_counter+=translate_error
        first_oligo_counter += add_first
        middle_oligo_counter += add_mid
        last_oligo_counter += add_last
    print("Lib: " + inputfiles[fileindex])
    if translate_error_counter>0:
        print(str(translate_error_counter)+" translation errors")
        total_translation_errors+=translate_error_counter
    if first_oligo_counter != constructs_per_lib[fileindex] or middle_oligo_counter != ((oligos_per_construct[fileindex]-2)*constructs_per_lib[fileindex]) or last_oligo_counter != constructs_per_lib[fileindex]:
        print(str(first_oligo_counter) + " start oligos. Expect to have " + str(constructs_per_lib[fileindex]))
        print(str(middle_oligo_counter) + " middle oligos or " + str(middle_oligo_counter/constructs_per_lib[fileindex]) + " per construct. Expect to have "+ str((oligos_per_construct[fileindex]-2)*constructs_per_lib[fileindex])+ " total middle oligos.")
        print(str(last_oligo_counter) + " end oligos. Expect to have "+ str(constructs_per_lib[fileindex]))
    else:
        print("Everything good.")
    

if total_translation_errors>0:
    print(str(total_translation_errors)+" TOTAL translation errors")

