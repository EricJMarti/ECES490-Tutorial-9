# ECES 490/690 Tutorial 9: FragGeneScan & HMMer

### FragGeneScan
For a full ReadMe, click [here](http://omics.informatics.indiana.edu/FragGeneScan/README).

__Basic Syntax:__
```bash
run_FragGeneScan.pl -genome=< FASTA file > -out=< output path > -complete=< whole genome? > -train=< sequencing type >
```

The __-genome__ flag specifies the input sequence path.

The __-out__ flag tells FragGeneScan where to place the output files and what to name them.

The __-complete__ flag is a boolean that specifies if the input is a complete genomic sequence or a composition of short reads.

The __-train__ flag tells FragGeneScan to use a particular set of sequencing parameters for scanning the input genome. The options are:
* Complete
* Sanger Sequencing (0.5% & 1% error)
* 454 Pyrosequencing (1% & 3% error)
* Illumina Sequencing (0.5% & 1% error)

__Output:__
FragGeneScan generates three output files:
* SeqName.out
  * This file holds the locations of various genes in the input sequence. It details the start position, end position, strand, frame, and score.
* SeqName.faa
  * This file holds the amino acid sequences that correspond to the genes described in the .out file. We will use this file for HMMer.
* SeqName.ffn
  * This file holds the nucleotide sequences that correspond to the genes described in the .out file.

__HUGE Issue with FragGeneScan:__ No multithreaded support (this package is REALLY SLOW). Some of our datasets took 9-12 hours to analyze completely. A newer version has been developed with support for cluster computing called FragGeneScan-Plus, but their code would not compile on Proteus. FragGeneScan-Plus claims to analyze data 5x faster than FragGeneScan on a single core, and about 50x faster using 8 cores. This package can be found [here](https://github.com/hallamlab/FragGeneScanPlus).

### HMMer
For the offical (120 page) HMMer manual, click [here](http://hmmer.janelia.org) to visit the offical HMMer website and click the link to download the documentation PDF.

HMMer has many functions:
* __hmmbuild__ - Build profile HMM from input multiple alignment
* __hmmalign__ - Make multiple aligment to a common profile HMM
* __phmmer__ - Search single protein sequence against protein sequence database (BLASTP-like)
* __jackhmmer__ - Iteratively search protein sequence against protein sequence database (PSIBLAST-like)
* __hmmsearch__ - Search protein profile HMM against a protein sequence database
* __hmmscan__ - Seach protein sequence against a protein profile HMM database
* __hmmpgmd__ - Search daemon for submitting jobs to hmmer.org
* __nhmmer__ - Search DNA sequence, alignment, or profile HMM against a DNA sequence database (BLASTN-like)
* __nhmmscan__ - Search DNA sequence against a DNA profile HMM database
* __hmmpress__ - Format HMM database into binary format for hmmscan

This package encourages the user to create an HMM profile for the data being analyzed. A sample HMM profile looks like this (taken directly from the HMMer3.1 user manual):
```
HMMER3/f [3.1 | February 2013] NAME fn3
ACC PF00041.13
DESC Fibronectin type III domain LENG 86
ALPH amino
RF no
MM no
CONS yes
CS yes
MAP yes
DATE Fri
NSEQ 106
EFFN 11.415833
CKSUM 3564431818
GA 8.00 7.20
TC 8.00 7.20
NC 7.90 7.90
STATS LOCAL MSV
STATS LOCAL VITERBI -9.7737 0.71847 STATS LOCAL FORWARD -3.8341 0.71847
HMM
-9.4043 0.71847
A C D E F G H I (...) Y
(...)
//
85 2.48488 5.72055
2.68618 4.42225
0.00338 6.08833
86 3.03720 5.94099
2.68618 4.42225 0.00227 6.08723
Feb 15 06:04:13 2013
m->m m->i COMPO 2.70330 4.91262 2.68618 4.42225 0.00338 6.08833
1 3.16986 5.21447 2.68629 4.42236 0.09796 2.38361
2 2.70230 5.97353 2.68618 4.42225 0.00338 6.08833
m->d i->m i->i d->m 3.03272 2.64079 3.60307 2.84344 2.77519 2.73123 3.46354 2.40513 6.81068 0.61958 0.77255 0.00000 4.52134 3.29953 4.34285 4.18764 2.77530 2.73088 3.46365 2.40512 6.81068 0.10064 2.34607 0.48576 2.24744 2.62947 5.31433 2.60356 2.77519 2.73123 3.46354 2.40513 6.81068 0.61958 0.77255 0.48576
3.87501 1.97538 3.04853 3.48010
2.77519 2.73123 3.46354 2.40513 6.81068 0.61958 0.77255 0.48576 3.75455 2.96917 5.26587 2.91682
2.77519 2.73123 3.46354 2.40513 * 0.61958 0.77255 0.00000
d->d
3.74204 3.07942 (...) 3.21526 3.72494 3.29354 (...) 3.61503
*
4.30886 3.35801 (...) 3.93889 3.72505 3.29365 (...) 3.61514 0.95510
4.43584 4.79731 (...) 4.25623 3.72494 3.29354 (...) 3.61503 0.95510
4.51877 3.51898 (...) 3.43366
3.72494 3.29354 (...) 3.61503 0.95510
3.66571 4.11840 (...) 4.99111
3.72494 3.29354 (...) 3.61503 *
1p---
3s---
120 e - - B
121 s - - E
```

Creating an HMM profile is not always possible because the parser in __hmmbuild__ is incredibly selective of what data format it accepts. The recommended data format is Stockholm, which looks like this (again, taken from the HMMer3.1 user manual): 
```
# STOCKHOLM 1.0
       seq1  ACDEF...GHIKL
       seq2  ACDEF...GHIKL
       seq3  ...EFMNRGHIKL
       seq1  MNPQTVWY
       seq2  MNPQTVWY
       seq3  MNPQT...
       //
```

In this tutorial, we revert to using __hmmpress__ and __hmmscan__ to analyze the protein sequences produced by FragGeneScan (since HMM profiles cannot be created with our data). The __hmmscan__ program will accept a protein sequence and compare it against an HMM database. Before __hmmscan__ can be used, the HMM database must be converted into a binary format to optimize the speed of the scan. The __hmmpress__ program prepares an HMM profile for __hmmscan__, and its basic syntax is as follows:
```bash
hmmpress < HMM profile >
```
The basic syntax for __hmmscan__ is as follows:
```bash
hmmscan <options> < HMM profile > < FAA file >
```
For the complete set of options that go with every command of HMMer, consult the user manual or the man pages using the __-h__ flag after the HMMer command.

__Speed of HMMer:__ The user manual states that HMMer 3.1 should be as fast as BLAST for "most protein queries." However, the older FASTA format does not fall into this category. Even with 32 cores available, our first dataset is STILL being processed as we present this tutorial!

### Pfam_Scan
EMBL provides an alternate method of scanning protein sequences against their Pfam database. They provide a tool called Pfam_Scan, which is essentially a Perl script that runs HMMer against the Pfam-A and Pfam-B databases. It also offers a "pretty" output format for easier interpretation of the data, which includes Pfam abundance info. We were not able to get this script to install on Proteus due to permission issues around modifying the Perl configuration. However, if you were to run this script, the general syntax is:
```bash
${PFAMPATH}/pfam_scan.pl -fasta < Protein Sequence > -out < Pfam Database Directory >
```
For more on Pfam, click [here](http://pfam.xfam.org).

## Our Workflow:
General settings for Proteus:
```bash
#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -M ejm335@drexel.edu
#$ -P nsftuesPrj
#$ -pe shm 32
#$ -l h_rt=48:00:00
#$ -q all.q@@amdhosts

. /etc/profile.d/modules.sh
module load shared
module load proteus
module load sge/univa
module load gcc/4.8.1
```
First, we need to download our data. This command downloads our dataset (we must then move it all into one folder):
```bash
wget -r --no-parent ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX/SRX144/SRX144807/
```
Once all of the data is in one place, we need to unzip the files. 
```bash
gunzip *.gz
```
We need to convert all of the SRA data files into a FASTA format.
```bash
# Define the path for SRA Toolkit & HMMer
PATH=/mnt/HA/groups/nsftuesGrp/.local/bin:$PATH

# Convert from SRA to FASTA
fastq-dump --fasta ${DATAPATH}/*.sra
```
Now we have all of the files we need to run FragGeneScan.
```bash
# Define Data/Package Locations
FRAGPATH=/home/ejm335/fraggenescan-code
TMPPATH=${TMP}/FragGeneScanOutput
DATAPATH=/home/ejm335/TutorialData
OUTPUTPATH=/home/ejm335/TutorialData/output

# Define data array
datafiles=(SRR492065 SRR492066 SRR492182 SRR492183 SRR492184 SRR492185 SRR492186 SRR492187 SRR492188 SRR492189 SRR492190 SRR492191 SRR492192 SRR492193 SRR492194 SRR492195 SRR492196 SRR492197)

# Make a folder to hold output
mkdir ${TMP}/FragGeneScanOutput

# Run FragGeneScan (we ended up running this separately per file, each file takes ~9 hours to analyze)
for file in ${datafiles[@]}
do
${FRAGPATH}/run_FragGeneScan.pl -genome=${DATAPATH}/${file}.fasta -out=${TMPPATH}/${file} -complete=0 -train=454_30
done

# Move output from scatch to home
mv ${TMPPATH}/* ${OUTPUTPATH}/
```
We now have three output files from each FASTA file. Sample of each:
__SRR492065.out__
```
>SRR492065.1 HWI-EAS385_0095_FC:2:1:6702:1434 length=200
1	72	-	3	1.435572	I:10,	D:
>SRR492065.2 HWI-EAS385_0095_FC:2:1:6931:1435 length=200
1	86	-	3	1.358242	I:	D:
>SRR492065.3 HWI-EAS385_0095_FC:2:1:9984:1431 length=200
1	120	-	1	1.331453	I:	D:
>SRR492065.4 HWI-EAS385_0095_FC:2:1:11577:1434 length=200
>SRR492065.5 HWI-EAS385_0095_FC:2:1:15121:1434 length=200
1	84	-	1	1.379883	I:	D:
```
__SRR492065.ffn__
```
>SRR492065.1_1_72_-
GTAACGACAATCGTGACAAGTCGCGATTGGAGTGGTCTTCACAGTGAAGCCAAGCATAGGATGGCT
>SRR492065.2_1_86_-
ATNNNNNTTTCTTCCGACAAAGTTGATCAAGTCGCTGAGTTTGGAAATTCTAGTAAAATCACAGTCGGTGAGCCTGCTATT
>SRR492065.3_1_120_-
GTAAAACAANNNNNGAAAGGCGGCGATTGGCGTGCTNNNNNTTTGAATATTGTTTTAATTCTAGTGGCCATCTTGATTTACTATCCGTTCTTTGTAGCTTATGATAAAAATGAGCTT
>SRR492065.5_1_84_-
NNNNNTAAATCAGCCCATGACTCAGGTGCTGCTCTACTATGTAAACATGACGGTGTCGTAGAATTCGTCGATGCCAAAGAA
```
__SRR492065.faa__
```
>SRR492065.1_1_72_-
VTTIVTSRDWSGLHSEAKHRMA
>SRR492065.2_1_86_-
XXXSSDKVDQVAEFGNSSKITVGEPAI
>SRR492065.3_1_120_-
VKQXXKGGDWRAXXLNIVLILVAILIYYPFFVAYDKNEL
>SRR492065.5_1_84_-
XXKSAHDSGAALLCKHDGVVEFVDAKE
>SRR492065.6_1_200_+
NYPCKIYIYDGNVKLENLNIDEDFIXXKEEIWKAXXTSXXXXXVPLNHXXILLQQLLLQTHFYAL
```
Now that we have run FragGeneScan, we now need to download the Pfams database in order to anotate each dataset's protein sequences (.faa). Navigate to the desired directory and run this command: 
```bash
wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam//releases/Pfam27.0/Pfam-A.hmm.gz
```
HMMer offers many commands to analyze data, but the most straightfoward one to use is hmmscan. This is because it takes in a protein sequence and compares it to a protein database in HMM format. In the previous command, we downloaded the Pfam-A database in the HMM format. Before we can use hmmscan, we need to compress this database into a binary format using hmmpress:
```bash
# Define the path for SRA Toolkit & HMMer
PATH=/mnt/HA/groups/nsftuesGrp/.local/bin:$PATH

# Define other paths
DBPATH=/home/ejm335/TutorialData
DATAPATH=/home/ejm335/TutorialData/output/
OUTPATH=/home/ejm335/TutorialData/output/HMMer

# Compress Pfam Database
hmmpress ${DBPATH}/Pfam-A.hmm
```
This generates four binary files in the ${DATAPATH} path we defined earlier. We are now all set to run hmmscan (and wait a VERY LONG time):
```bash
## Run HMMer for each protein sequence
for file in ${datafiles[@]}
do
hmmscan -o ${OUTPATH}/${file} --cpu 32 ${DBPATH}/Pfam-A.hmm ${DATAPATH}/${file}.faa
done
```
This creates a long output file in the ${OUTPATH} directory that looks like this:
```
# hmmscan :: search sequence(s) against a profile database
# HMMER 3.1b1 (May 2013); http://hmmer.org/
# Copyright (C) 2013 Howard Hughes Medical Institute.
# Freely distributed under the GNU General Public License (GPLv3).
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# query sequence file:             /home/ejm335/TutorialData/output//SRR492190.faa
# target HMM database:             /home/ejm335/TutorialData/Pfam-A.hmm
# number of worker threads:        32
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

Query:       SRR492190.4_1_200_+  [L=65]
Scores for complete sequence (score includes all domains):
   --- full sequence ---   --- best 1 domain ---    -#dom-
    E-value  score  bias    E-value  score  bias    exp  N  Model    Description
    ------- ------ -----    ------- ------ -----   ---- --  -------- -----------
    1.3e-06   28.2   0.1    1.7e-06   27.8   0.1    1.2  1  IF-2      Translation-initiation factor 2


Domain annotation for each model (and alignments):
>> IF-2  Translation-initiation factor 2
   #    score  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc
 ---   ------ ----- --------- --------- ------- -------    ------- -------    ------- -------    ----
   1 !   27.8   0.1   1.2e-10   1.7e-06      15      48 ..       2      35 ..       1      42 [. 0.90

  Alignments for each domain:
  == domain 1  score: 27.8 bits;  conditional E-value: 1.2e-10
                 IF-2 15 vkelnvivkaDvqGsleAlkesLeklsteevkvk 48
                          ke+n+ivkaDvqG  eA+++sL+  + e v++ 
  SRR492190.4_1_200_+  2 FKEVNIIVKADVQGXDEAVSASLQXXDVEGVRLX 35
                         589************************9999986 PP
```
We can use the __--pfamtblout__ flag to produce a formatted output, and we can limit the inclusion threshold by specifying an inclusion e-value of 1e-5 with the __-incE__ flag. 
```bash
# HMMScan with formatted output
hmmscan --pfamtblout ${OUTPATH}/${file} --cpu 32 --incE 0.00001 ${DBPATH}/Pfam-A.hmm ${DATAPATH}/${file}.faa
```
This produces output that looks like this (first 25 matches in the SRR492190 dataset):
```
# name                  bits   E-value   n   exp  bias    description
# ------------------- ------ --------- --- ----- -----    ---------------------
IF-2                    28.2   1.3e-06   1   1.2   0.1    Translation-initiation factor 2
Asp23                   16.5    0.0074   1   1.1   0.0    Asp23 family
RIX1                    32.3   6.6e-08   1   1.1   0.2    rRNA processing/ribosome biogenesis
DUF2969                 23.2   5.1e-05   1   1.1   0.2    Protein of unknown function (DUF2969)
Malic_M                 13.6     0.036   1   1.1   0.0    Malic enzyme, NAD binding domain
Glyco_hydro_35          14.1     0.022   1   1.0   0.0    Glycosyl hydrolases family 35
ThiW                    35.8   6.3e-09   1   1.0   0.8    Thiamine-precursor transporter protein (ThiW)
YbaB_DNA_bd             15.5     0.011   1   1.0   0.1    YbaB/EbfC DNA-binding family
TPR_2                   24.8   1.2e-05   1   1.3   0.0    Tetratricopeptide repeat
TPR_1                   24.3   1.5e-05   1   1.4   0.0    Tetratricopeptide repeat
TPR_8                   24.3   1.6e-05   1   1.2   0.0    Tetratricopeptide repeat
TPR_12                  23.3   4.4e-05   1   1.1   0.0    Tetratricopeptide repeat
TPR_7                   22.9   4.6e-05   1   1.4   0.0    Tetratricopeptide repeat
TPR_11                  18.6     0.001   1   1.1   0.0    TPR repeat
TPR_14                  18.1    0.0028   1   1.3   0.0    Tetratricopeptide repeat
TPR_10                  16.9    0.0041   1   1.3   0.0    Tetratricopeptide repeat
TPR_16                  16.9    0.0076   1   1.2   0.1    Tetratricopeptide repeat
TPR_6                   15.9     0.013   1   1.2   0.0    Tetratricopeptide repeat
CbiQ                    16.8    0.0037   1   1.1   0.0    Cobalt transport protein
DUF1129                 13.8     0.026   1   1.1   1.0    Protein of unknown function (DUF1129)
SHMT                    19.6   0.00024   1   1.0   1.0    Serine hydroxymethyltransferase
SAICAR_synt             17.3    0.0018   1   1.1   0.1    SAICAR synthetase
H_kinase_N              19.8   0.00046   1   1.1   0.0    Signal transduction histidine kinase
Cna_B                   19.5   0.00059   1   1.2   0.5    Cna protein B-type domain
OTCace                  15.4     0.012   1   1.1   0.0    Aspartate/ornithine carbamoyltransferase, Asp/Orn binding domain
MoaE                    25.4   9.4e-06   1   1.1   0.0    MoaE protein
```
