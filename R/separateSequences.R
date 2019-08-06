separateSequences  <- function(in.file, filter_files, identity_threshold) {
	cat('Separating sequences according to the filter files...\n')
	tstart <- Sys.time()
	fid.in <- file(in.file,'rb')
	blockSize = 2^24
	hNdx <- {}
	pos <- 0
	cnt <- blockSize
	while(cnt == blockSize) {
		aux  <- readChar(fid.in, nchars=blockSize, useBytes=T)
		cnt<- nchar(aux)
		hNdx <- c(hNdx, pos + gregexpr(aux, pattern='>')[[1]])
		pos  <- seek(fid.in)
	}
	hNdx <- c(hNdx, pos + cnt)
	seek(fid.in, where=0)
	numItems <- length(hNdx)-1
	headers <- rep('', numItems);
	for(i in 1:numItems) {
		seek(fid.in, hNdx[i]);
		headers[i] = readLines(con = fid.in, n = 1)
	}
	label = rep(1, numItems)
	if(length(filter_files) > 0) {
		for(i in 1:length(filter_files)) {
			d <- tryCatch(read.table(filter_files[i], sep = ',', header = F),
				      error = function(e) data.frame(),
				      finally = function(e) data.frame())
			if(nrow(d) > 0) {
				d <- d[d[[3]] > identity_threshold,]
				d <- d[order(d[[1]]),]
				label[is.element(headers, d[[1]])] <- i+1
			}
		}
	}
	tstamp = Sys.time()
	out.files = {}
	for(i in unique(label)) {
		fid.out <- NA
		if(i == 1) {
			filename <- paste0(tstamp, '_filtered-unlabeled.fasta')
			fid.out <- file(filename, 'wb')
			cat('   -  unlabeled: ', sum(label == 1), '\n')
		} else {
			filename <- basename(filter_files[[i]])
			filename <- substr(filename, 1, nchar(filename)-4)
			filename <- paste0(tstamp, '_filtered-', filename, '.fasta')
			fid.out <- file(filename, 'wb')
			cat('   -  ', filename, ': ', sum(label==i))
		}
		out.files <- c(out.files, filename)
		for(j in which(label == i)) {
			chunkStart <- hNdx[j]-1
			seek(fid.in, chunkStart)
			chunk <- ""
			if(j < numItems) {
				chunkEnd <- hNdx[j+1]-1
				chunk <- readChar(fid.in, nchars = chunkEnd - chunkStart)
			} else
				chunk <- readChar(fid.in, nchars = 1e6)
			writeChar(chunk, con=fid.out, eos=NULL)
		}
		close(fid.out)
	}
	close(fid.in)
	cat('   - Elapsed time: ', Sys.time() - tstart, ' sec\n')
	return(out.files)
}
