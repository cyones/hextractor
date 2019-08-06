splitFasta <- function(filename) {
	fid <- file(filename, 'rb')
	blockSize <- 2^24
	hNdx <- {}
	pos <- 0
	cnt <- blockSize
	while (cnt == blockSize) {
		aux  <- readChar(fid, nchars=blockSize, useBytes=T)
		cnt <- nchar(aux)
		mtc <- gregexpr(aux, pattern='>')[[1]]
		hNdx <- c(hNdx, pos + mtc[mtc >= 0])
		pos  <- seek(fid)
	}
	hNdx <- c(hNdx, pos + cnt)
	seek(fid, where=0)
	numItems <- length(hNdx)-1
	for (fileNum in 1:numItems) {
		outName<- sprintf('%s.%d', filename, fileNum)
		fidOut <- file(outName, 'w')
		chunkStart <- hNdx[fileNum]-1
		chunkEnd <- hNdx[fileNum+1]-1
		chunkSize <- chunkEnd - chunkStart
		seek(fid, chunkStart)
		chunk <- readChar(fid, nchars=chunkSize, useBytes=T)
		writeChar(chunk, fidOut, useBytes=T, eos=NULL)
		close(fidOut)
	}
	close(fid)
	return(numItems)
}
