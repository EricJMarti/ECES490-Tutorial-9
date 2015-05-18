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
__Biggest Issue with FragGeneScan:__ No multithreaded support (this package is REALLY SLOW). Some of our datasets took 9-12 hours to analyze completely. A newer version has been developed with support for cluster computing called FragGeneScan-Plus, but their code would not compile on Proteus. FragGeneScan-Plus claims to analyze data 5x faster than FragGeneScan on a single core, and about 50x faster using 8 cores. This package can be found [here][https://github.com/hallamlab/FragGeneScanPlus].

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

__Our Workflow:__
First, we need to download the Pfams database in order to run HMMer.
