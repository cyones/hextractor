#' @importFrom utils unzip
execBlast  <- function(seq_file, filter_files, blast_evalue, nthreads) {
	if(length(filter_files) == 0)
		return(filter_files)
	cat('Executing blast...\n')
	tstart = Sys.time()
	cmd <- paste0('makeblastdb -max_file_sz 2GB -in ', seq_file,
		      ' -dbtype nucl > /dev/null')
	system(cmd)
	for(i in 1:length(filter_files)) {
		ffile <- filter_files[i]
		uncompressedFile = {}
		if(substr(ffile, nchar(ffile)-2, nchar(ffile)) == 'zip') {
			uncompressedFile = unzip(ffile)
			ffile = uncompressedFile$Name[1]
		}
		filter_files[i] <- tempfile()
		cmd <- paste0('blastn -query ', ffile,
			      ' -outfmt 10 -db ', seq_file,
			      ' -strand plus -evalue ', blast_evalue,
			      ' -num_threads ', nthreads,
			      ' -out ', filter_files[i])
		system(cmd)
		for(uf in uncompressedFile)
			file.remove(uf)
	}
	cat('   - Elapsed time: ', (Sys.time() - tstart) / 60, ' min\n')
	return(filter_files)
}
