#' HextractoR: Integrated Tool for Hairping Extraction of RNA Sequences
#'
#' To preprocess a genome, you need a file containing the raw genome in fasta
#' format. To run HExtractor, simply call the main function. This function
#' creates 2 files in the "out" folder and automatically names them.
#' @param input_file filename of the fasta file to proccess
#' @param window_size Number of bases in the windows.
#' @param window_step Window step. This number defines indirectly the overlap:
#' window_overlap=window_size-window_step
#' @param min_length Minimum sequence length. Shorter sequences are discarded.
#' @param min_bp Minimum number of base-pairs that must form a sequence.
#' @param margin_bp When the sequence is trimmed, at least min_bp+margin_bp
#' base-pairs are left.
#' @param min_valid_nucleotides Each input sequence must have this quantity of
#' valid nucleotides (not 'N') to be processed.
#' @param only_sloop Only extract single loop sequence.
#' @param trim_sequences Use some heuristics to trim the hairpins.
#' @param blast_evalue e-value used in blast to match the extracted sequences
#' with the sequences from the filter files.
#' @param identity_threshold Identity threshold used to match sequences with the
#' sequences from the filter files.
#' @param nthreads Allows using more than one thread in the execution.
#' @param nworks Split each sequence in nworks to use less RAM memory.
#' @param filter_files Fasta files with known sequences to separate the output
#' stems.
#' @return A list with the path of the output files and the result of the
#' proccessing of each sequence (if it was succesful or failed)
#' @examples
#' # Small example without filter files
#' library(HextractoR)
#' # First we get the path of the example FASTA file
#' fpath <- system.file("Example_tiny.fasta", package="HextractoR")
#' # To run HextractoR, simply call the main function
#' HextractoR(input_file = fpath)
#' \donttest{# Other example with filter files and bigger input file
#' fpath1 <- system.file("Example_human.fasta", package="HextractoR")
#' fpath2 <- system.file("Example_pre-miRNA.fasta", package="HextractoR")
#' HextractoR(input_file = fpath1, filter_files = {fpath2})
#' # This function creates 2 files in the working directory and automatically
#' # names them.}
#'
#' @import seqinr
#' @import parallel
#' @import doParallel
#' @import foreach
#' @import utils
#' @export
HextractoR <- function(input_file,
		       min_valid_nucleotides = 500,
		       window_size = 160,
		       window_step = 30,
		       only_sloop = T,
		       min_length = 60,
		       min_bp = 16,
		       trim_sequences = T,
		       margin_bp = 6,
		       blast_evalue = 1,
		       identity_threshold = 90,
		       nthreads = 4,
		       nworks = 4,
		       filter_files = {}) {
	if(!checkRequeriments())
		return()
	if(nthreads > nworks)
		nworks <- nthreads
	workers <- makeCluster(nthreads,type="SOCK")
	registerDoParallel(workers)
	numItems <- splitFasta(input_file)
	status <- rep(0, numItems)
	error.log <- rep("Ok", numItems)
	for (iseq in 1:numItems){
		partFile <- paste0(input_file, '.', iseq)
		rnaseq <- read.fasta(partFile, as.string=T, forceDNAtolower=T)
		nns <- nchar(rnaseq[[1]]) - length(gregexpr(rnaseq[[1]], pattern='n')[[1]])
		if (nns < min_valid_nucleotides) {
			error.log[iseq] <- paste0('Skipping sequence ', names(rnaseq),
						  ' due to lack of valid nucleotides.')
			next
		}
		files <- windowSequence(rnaseq,
					window_size = window_size,
					window_step = window_step,
					min_length = min_length,
					nworks = nworks)
		files <- foldSequences(files, nworks = nworks)
		files <- extractStems(files, min_length = min_length,
				      min_bp = min_bp, margin = margin_bp,
				      trim_sequences = trim_sequences, nworks = nworks)
		files <- uniqueSequences(files, nworks = nworks)
	}
	filter_files <- execBlast(files, filter_files = filter_files,
				  blast_evalue = blast_evalue,
				  nthreads = nthreads)
	files <- separateSequences(files, filter_files = filter_files,
				   identity_threshold = identity_threshold)
	stopCluster(workers)
	return(list(result_files = files, error.log = error.log))
}


checkRequeriments <- function() {
	rnafold <- system('RNAfold --version', ignore.stdout = T, ignore.stderr = T)
	blast <- system('formatdb --help', ignore.stdout = T, ignore.stderr = T)
	if(rnafold == 127) {
		cat("RNAfold not installed. Download the last version from:\n")
		cat("https://www.tbi.univie.ac.at/RNA/\n")
	}
	if(blast  == 127) {
		cat("BLAST not installed. Download the last version from:\n")
		cat("ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/\n")
	}
	return(blast != 127 && rnafold != 127)
}
