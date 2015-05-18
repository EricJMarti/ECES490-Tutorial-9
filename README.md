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

__Misc Notes:__
__Biggest Issue with FragGeneScan:__ No multithreaded support (this package is REALLY SLOW). Some of our datasets took 9-12 hours to analyze completely. A newer version has been developed with support for cluster computing called FragGeneScan-Plus, but their code would not compile on Proteus. FragGeneScan-Plus claims to analyze data 5x faster than FragGeneScan on a single core, and about 50x faster using 8 cores. This package can be found [here](https://github.com/hallamlab/FragGeneScanPlus).

### HMMer
For the offical (120 page) HMMer manual, click [here](http://hmmer.janelia.org) to visit the offical HMMer website and click the link to download the documentation PDF.

HMMer has many functions:
* hmmbuild - Build profile HMM from input multiple alignment
* hmmalign - Make multiple aligment to a common profile HMM
* phmmer - Search single protein sequence against protein sequence database (BLASTP-like)
* jackhmmer - Iteratively search protein sequence against protein sequence database (PSIBLAST-like)
* hmmsearch - Search protein profile HMM against a protein sequence database
* hmmscan - Seach protein sequence against a protein profile HMM database
* hmmpgmd - Search daemon for submitting jobs to hmmer.org
* nhmmer - Search DNA sequence, alignment, or profile HMM against a DNA sequence database (BLASTN-like)
* nhmmscan - Search DNA sequence against a DNA profile HMM database
* hmmpress - Format HMM database into binary format for hmmscan

In this tutorial, we use __hmmpress__ and __hmmscan__ to analyze the protein sequences produced by FragGeneScan. The basic syntax for hmmpress is as follows:
```bash
hmmpress < HMM profile >
```
The basic syntax for hmmscan is as follows:
```bash
hmmscan <options> < HMM profile > < FAA file >
```

### Pfam_Scan
EMBL provides an alternate method of scanning protein sequences against their Pfam database. They provide a tool called Pfam_Scan, which is essentially a Perl script that runs HMMer against the Pfam-A and Pfam-B databases. It also offers a "pretty" output format for easier interpretation of the data. We were not able to get this script to install on Proteus due to permission issues around modifying the Perl configuration. If you were to run this script, the general syntax is:
```bash
${PFAMPATH}/pfam_scan.pl -fasta < Protein Sequence > -out < Pfam Database Directory >
```
For more on Pfam, click (here)[http://pfam.xfam.org].

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
First we need to download our data. This command downloads our dataset (we must then move it all into one folder):
```bash
wget -r --no-parent ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX/SRX144/SRX144807/
```
Once all of the data is in one place, we need to convert all of these data files into FASTA format.
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
