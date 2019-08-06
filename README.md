This is the version 1.4 of the HextractoR toolbox
-------------------------------------------------
Simple and integrated tool that automatically extracts and folds all hairpin
sequences from raw genome-wide data. It predicts the secondary structure of
several overlapped segments, with longer length than the mean length of
sequences of interest for the species under processing, ensuring that no one is
lost nor inappropriately cut.

The last version can be found at:
http://sourceforge.net/projects/sourcesinc/files/HextractoR/

Contact
- Cristian Yones <cyones@sinc.unl.edu.ar>
- sinc(i):  http://fich.unl.edu.ar/sinc/

Installation
------------
HextractoR was tested in R > 3.4, but it probably work in older versions as well.
For folding and alignment of sequences the toolbox uses some external well-known
packages. The following software must be installed:

- Vienna RNA: this package is used to fold sequences and alignments.
  Download: https://www.tbi.univie.ac.at/RNA/

- NCBI blast+: this package is used to identify the known miRNAs and the
  sequence to discard.
  Download: ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/

To install the R package simply execute from R:

> install.packages("HextractoR", dependencies=T)

How to preprocess a genome
--------------------------
To preprocess a genome, you need a file containing the raw genome in fasta format.
To run HExtractor, simply call the main function. This function creates 2 files
in the "out" folder and automatically names them. The input parameters are:
- input_file filename of the fasta file to proccess
- window_size Number of bases in the windows.
- window_step Window step. This number defines indirectly the overlap:
  window_overlap=window_size-window_step
- min_length Minimum sequence length. Shorter sequences are discarded.
- min_bp Minimum number of base-pairs that must form a sequence.
- margin_bp When the sequence is trimmed, at least min_bp+margin_bp base-pairs
  are left.
- min_valid_nucleotides Each input sequence must have this quantity of valid
  nucleotides (not 'N') to be processed.
- only_sloop Only extract single loop sequence.
- trim_sequences Use some heuristics to trim the hairpins.
- blast_evalue e-value used in blast to match the extracted sequences with the
  sequences from the filter files.
- identity_threshold Identity threshold used to match sequences with the
  sequences from the filter files.
- nthreads Allows using more than one thread in the execution.
- nworks Split each sequence in nworks to use less RAM memory.
- filter_files Fasta files with known sequences to separate the output stems.

HextractoR function returns a list with the path of the output files and the
result of the proccessing of each sequence (if it was succesful or failed)

Examples
--------
# Small example without filter files library(HextractoR)
# First we get the path of the example FASTA file
fpath <- system.file("Example_tiny.fasta", package="HextractoR")
# To run HextractoR, simply call the main function
HextractoR(input_file = fpath)

# Other example with filter files and bigger input file
fpath1 <- system.file("Example_human.fasta", package="HextractoR")
fpath2 <- system.file("Example_pre-miRNA.fasta", package="HextractoR")
HextractoR(input_file = fpath1, filter_files = {fpath2})
# This function creates 2 files in the working directory and automatically
# names them.

