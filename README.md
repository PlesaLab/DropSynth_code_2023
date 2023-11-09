## Warning: this software is still in beta, some features have not been tested, use at own risk

## DropSynth oligo generation script

This set of scripts provides an example for how to generate DropSynth oligos for a set of libraries. We will use 4 libraries from this paper. For updated versions of this code please check [https://github.com/PlesaLab](https://github.com/PlesaLab)

#### Requirements:
* anaconda (see the ds063022.yml file)
 * python 3.5 environment
 * python-levenshtein
 * biopython
 * numpy
 * cython
 * (_optional_) matplotlib 
* [**seqfold**](https://github.com/Lattice-Automation/seqfold)
* [BLAT](https://genome.ucsc.edu/goldenpath/help/blatSpec.html) (faToTwoBit, pslPretty, gfClient, gfServer)
* (_optional_) [unafold](http://www.unafold.org/Dinamelt/software/obtaining-unafold.php) (only require hybrid-ss-min), this requires the appropriate [licence](http://www.unafold.org/Dinamelt/software/obtaining-unafold.php)

#### New features

* switched to a “recipe” based workflow with all parameters in a single file
* ability to use any cloning sites
* uses Lattice-Automation’s seqfold python library for minimum free energy structure calculations, much faster than unafold (hybrid-ss-min)
* implemented a programmable database for handling all restriction enzyme sites required
* implements the ability to split as many genes as necessary in first step with subsequent (384x, 1536x) library splitting
* virtual assembly and translation to verify oligo designs
* the option to do barcode reversal between libraries to offset barcoded bead effects
* improved codon optimization with lower split failures through the use of hardcoded rules
* added in the ability to require certain sequences (controls) in each library
* (Beta) Single oligo processing for very small genes
* (Beta) support for DNA (non-protein) constructs
* improved oligo junction length handling, with genes that fail due to length placed into a special file for input into higher oligo splits.


#### Sequence input

The input sequences can be found in the folders:

* fasta_input
* controls_input

You should sort your input sequences and place them into these folders. We typically split the inputs by lengths, corresponding to how many oligos are required to assemble them.

#### Configuration files

Each library design is controlled by a configuration file which can be found in the folders:

* lib_config

Each configuration file has the following variables:

1. codon\_usage\_file - codon frequencies to use for target organism. These are stored in codon_data for E. coli (W3110, Kazusa ID 316407), Human, and Yeast. Data is collected from [https://www.kazusa.or.jp/codon/](https://www.kazusa.or.jp/codon/). Eg. for [E. coli (W3110)](https://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?species=316407&aa=5&style=N)
2. output\_folder - where output files will go
3. input\_type - can be AA or DNA. DNA is still in Beta testing.
4. input\_file - source fasta file with input sequences
5. num\_oligos - (integer) the number of oligos to split these into. See the table at the bottom of [dropsynth.org](http://dropsynth.org/)
6. add\_stop\_codon - add stop (TAA) codon before end cloning site?
7. seeds\_for\_libs - random seed used for codon optimization
8. IDs\_to\_remove - location of file to store IDs of failed sequences
9. cloning\_fwd - the cloning sequence at the start of the gene
10. clonging\_fwd\_includes\_ATG - (True/False) flag if ATG is included
11. check\_starting\_ATG - (True/False) flag
12. cloning\_rev - the cloning sequence at the end of the gene
13. REsites\_file - the location of the file with restriction site definitions (see below)
14. gen\_two\_codon - (True/False) make two versions of each library?
15. furthest\_codon\_file - location of json lookup table for making codon 2
16. switch\_barcodes\_for\_second\_lib - (True/False) flag. If this is set to true the codon 2 will use microbead barcodes in the reverse order. This mitigates potential barcoded microbead effects, since the same gene no longer has the same barcode in each version.
17. total\_oligo\_length - (integer) oligo length to be ordered (typically 200 to 300)
18. seq\_req\_for\_processing - (integer) do not modify unless changing DropSynth process, always 70
19. max\_num\_attempts - (integer) how many times to retry making a gene
20. padding\_var - (True/False) use padding? (add random sequence between final cloning site and the reverse assembly primer). Highly recommended, as this makes gel isolation far easier since the distribution of lengths is tight.
21. padding\_length - (integer) what length to pad to?
22. max\_homopolymer\_repeat\_length - (integer) block any homopolymer repeats of this length or longer
23. assemblyprimf\_file - asmF primer file. Usually only 504 is used.
24. assemblyprimr\_file - asmR primer file. Usually only 504 is used.
25. diff\_assembly\_pim\_on\_two\_codon - (True/False) use different assembly primers for codon 2? This is a hack for one oligo constructs.
26. assemblyprimf\_file\_codon2 - asmF primer file for codon 2 if diff\_assembly\_pim\_on\_two\_codon set to True. Usually n/a.
27. assemblyprimr\_file\_codon2 - asmR primer file for codon 2 if diff\_assembly\_pim\_on\_two\_codon set to True. Usually n/a.
28. lengthleeway - (integer) 
29. positionleeway - (integer) 
30. avgoverlapsize - (integer) 
31. overlaptemps - 
32. deltaGThreshold - (integer) 
33. selfDimersThreshold - (integer) 
34. oligo\_input\_file - 
35. ampprimf\_file - 
36. ampprimr\_file - 
37. primers\_to\_use - integers separated by commas indicating which pair of skpp15 primers to use for subpool amplifcation
38. constructs\_per\_lib - how many total constructs in this library. Usually 1536 or 384.
39. barcode\_set - location of file with microbead barcodes. Make sure you're using the correct file for 384, or 1536.
40. length\_padded\_payload - integer. The length of payload + BtsaI sites + buffer such that all oligos are full length.
41. needs\_split - (True/False)
42. output\_location - which folder to place the split files
43. div\_lib\_size - (integer) 
44. number\_of\_libs - (integer) how many smaller libs to split the input into
45. seqs\_to\_require\_file - location of fasta file. Force these sequences into each of the split libraries. Useful for controls.

Each library design's controlled restriction site sequences are defined in a file (see REsites/sites.csv). This file is set in the REsites\_file variable of the lib config. Note that non-palindromic RE sites need two rows, one for each strand. This file has the following options:

1. name - the name to give this RE site
2. site - the RE sequence
3. self\_overlap\_possible - set this to Yes if this site can overlap itself
4. gene\_no\_buffer - How to treat this site when in just the raw gene seqeunce, without buffer. Use 0 to block. Use 1 to require. Use -1 to ignore.
5. gene\_with\_buffer - How to treat this site when in gene seqeunce with buffer and cloning sites. Use 0 to block. Use 1 to require. Use -1 to ignore.
6. gene\_with\_asmprim - How to treat this site when in gene seqeunce with buffer, cloning sites, and assembly primers added. Use 0 to block. Use 1 to require. Use -1 to ignore.
7. after\_split - How to treat this site after the gene is split into fragments. Use 0 to block. Use 1 to require. Use -1 to ignore.
8. add\_btsi - How to treat this site after the gene is split into fragments and BtsI is added to either side. Use 0 to block. Use 1 to require. Use -1 to ignore.
9. final - How to treat this site in the final oligo. Use 0 to block. Use 1 to require. Use -1 to ignore.
10. random\_dna\_gene\_padding - How to treat this site in the random DNA use to pad gene length. Use 0 to block. Use 1 to require. Use -1 to ignore.
11. random\_dna\_oligo\_padding - How to treat this site in the random DNA use to pad oligo length. Use 0 to block. Use 1 to require. Use -1 to ignore.
12. notes - for you own notes

#### Procedure

A general overview of the oligo generation procedure can be seen in the shell script _protocol.sh_. This file calls a number of scripts:

1. **_primer\_screen.py_**  
Used to screen assembly and amplification primers. Check which skpp-15 primers will work for amp and which skpp-20 primers for asm. For the skpp-15 check the first 200 primer pairs. The asm primers (skpp-20) are offset by 500. This also writes the Amp pairs into one single fasta for use with BLAT later on. Also outputs the individual Asm pairs.

2. **_split\_genes\_for\_oligos.py -c Library\_A.conf > logs/Library\_A.log_**  
This is the main script which generates the gene sequences and splits them. Before running check the following:  
 * Input library files.
 * Number of oligos to split for each file.
 * Generate multiple codon versions for each sequence?
 * Maximum length of payload.
 * Buffer size for padding. This padding is random sequence added between KpnI and the reverse assembly primer. This is used to improve gel based size selection after assembly.
 * Max attempts in splitting.
 * Stop codon present or added.
 * Restriction sites to be used.

 This will follow these processing steps:  

 * Open the big library protein fasta file (one 4-oligo file and another 5 oligo file)
 * Loop and retry oligo generation until max\_num\_attempts is reached, if it fails choose the next sequence. 
 * For each protein sequence:
     * Generate one or two codon sequences.
     * Remove illegal restriction sites (1).
     * Verify that the DNA sequences still translate correctly.
     * Add gene cloning restriction sites (NdeI & KpnI) and stop codon.
     * Add random sequence buffer, to pad sequence to set length. This should not have long homopolymer repeats.
     * Remove illegal restriction sites (2).
     * Add the assembly primers.
     * Check for illegal restriction sites (3).
     * Split the DNA sequence into pieces.
     * Make sure the max oligo length and number of pieces is good.
     * Make sure no restriction sites are present after BtsaI sites are added on.
     * Write files with:
         * Split oligo payloads.
         * Protein sequences
         * Final DNA sequences
         * DNA sequences without RE or primers
         * DNA sequences without primers
  * Do this until all of records in the library are processed. If any sequences fail these get added to _badIDarray_ and will not be used the next time the script is run.

  Most common splitting fail errors:

  * sequence is too big for maxpayload limit.
  * sequence is too small for design parameters (if buffering is off)
  * seq is partial and doesn't start with Met
  * restriction site is introduced which can't be removed by codon mutation

3. **_divide\_into\_libs.py -c Library\_A.conf > logs/Library\_A\_divide.log_**  
This will subdivide the library into multiple smaller libraries.


4. **_python oligobuffergen.py -c Library\_A.conf > logs/Library\_A\_buffergen.log_**  
Finalize DropSynth oligos by adding amp primers, microbead barcodes, and nt.BspQI nicking sites.

5. **_BLAT\_primers.sh_**  
Use BLAT to check for alignments between the DNA sequences and the assembly primers. If there are large matches change the primers and regenerate the corresponding oligo libraries. First check each assembly primer against each Lib. Second check each amplification primer against all Libs.

6. **_BLAT\_hits\_parse.py_**  
Find good amp primer pairs. This will generate skpp15-forward_select_mod.faa and skpp15-reverse_select_mod.faa in the ampprimers-skpp15 directory.

7. **_final\_oligo\_check.py_**  
Run a bunch of design rule checks on the final oligos.

8. **_FASTA\_to\_csv\_proteins.py_**  
make a csv file with all protein seqs (also done in make_master_CSV.py).

9. **_FASTA\_to\_csv\_for\_Chip.py_**  
Make a csv file with all oligo seqs.

10. **_homopolymer\_stats.py_**  
Check for long homopolymers and display stats.


