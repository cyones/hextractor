foldSequences <- function(files, nworks) {
	cat('Folding...\n')
	tstart = Sys.time()
	i = NULL
	files <- foreach(i = 1:nworks, .inorder = F, .combine='c') %dopar% {
		ofile <- tempfile()
		system(paste0('RNAfold --noPS < ', files[i], ' > ', ofile))
		ofile
	}
	cat('   - Elapsed time: ', (Sys.time() - tstart) / 60, ' min\n')
	return(files)
}

