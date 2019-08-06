windowSequence <- function(rnaseq, window_size, window_step, min_length, nworks) {
	cat('- Windowing...\n')
	tstart = Sys.time()
	res <- {}
	seqName <- names(rnaseq)
	sq <- gsub('t', 'u', rnaseq)
	sc <- sq
	sc <- gsub('a', '1', sc)
	sc <- gsub('c', '2', sc)
	sc <- gsub('g', '3', sc)
	sc <- gsub('u', '4', sc)
	sc <- gsub('1', 'u', sc)
	sc <- gsub('2', 'g', sc)
	sc <- gsub('3', 'c', sc)
	sc <- gsub('4', 'a', sc)
	to_add <- data.frame(beg={}, end={})
	n_add <- 0
	for (i in seq(1, nchar(sq)-window_size, window_step)) {
		ssq <- substr(sq, start=i, stop=i+window_size)
		ibeg <- gregexpr(ssq, pattern = '[acgu]+')[[1]]
		iend <- ibeg + attr(ibeg, 'match.length')
		for (j in 1:length(ibeg)) {
			if (iend[j]-ibeg[j] > min_length) {
				n_add <- n_add+1
				to_add[n_add, 1] <- ibeg[j]-1+i
				to_add[n_add, 2] <- iend[j]-1+i
			}
		}
	}
	to_add <- unique(to_add)
	n_add <- 2 * nrow(to_add)
	for (i in 1:nworks) {
		res <- c(res, tempfile())
		fid <- file(tail(res, n=1), 'w')
		fseq <- 1+round((i-1) * nrow(to_add) / nworks)
		lseq <- round(i * nrow(to_add) / nworks)
		for (j in fseq:lseq){
			aux <- sprintf('>%s_%i-%i\n%s\n',
				       seqName, to_add[j,1], to_add[j,2],
				       substr(sq, to_add[j,1], to_add[j,2]))
			write(fid, x = aux)
		}
		for (j in fseq:lseq){
			aux <- sprintf('>%s_%i-%i_InvCom\n%s\n',
				       seqName, to_add[j,1], to_add[j,2],
				       substr(sc, to_add[j,1], to_add[j,2]))
			write(fid, x = aux)
		}
		close(fid)
	}
	cat('   -> Number of sequences: ', n_add, '\n')
	cat('   -> Elapsed time: ', Sys.time() - tstart, ' sec', '\n')
	return(res)
}

