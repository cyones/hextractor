extractStems <- function(files, min_length, min_bp, margin, trim_sequences, nworks) {
	cat('Extracting stems...\n')
	tstart <- Sys.time()
	out_files = {}
	for (i in 1:nworks)
		out_files <- c(out_files, tempfile())
	foreach (i = 1:nworks, .inorder = F) %dopar% {
		in_seqs <- read.fasta(files[i], as.string=T)
		fid <- file(out_files[i], 'w')
		for (i in 1:length(in_seqs)) {
			rfo <- parseRNAfoldOut(in_seqs[[i]])
			if (BPNumber(rfo$str) < 5)
				next
			reg <- proccessWindow(rfo, min_length, min_bp,
					      margin, trim_sequences)
			for(r in reg)
				write(r, file = fid)
		}
		close(fid)
	}
	cat('   - Elapsed time: ', Sys.time() - tstart, ' sec\n')
	return(out_files)
}

BPNumber <- function(sq) {
	sum(attr(gregexpr(sq, pattern='\\(')[[1]], 'match.length'))
}

parseRNAfoldOut <- function(string) {
	rfo <- {}
	re <- gregexpr(text = string, pattern = '[acgut]+')[[1]]
	rfo$seq <- substr(string, re, re + attr(re, 'match.length') - 1)
	re <- gregexpr(text = string, pattern = '[\\(\\.\\)]+')[[1]]
	sel <- attr(re, 'match.length') > 1
	rfo$str <- substr(string, re[sel], re[sel] + attr(re, 'match.length')[sel] - 1)
	rfo$par <- rep(NA, nchar(rfo$str))
	pila <- rep(0, nchar(rfo$str))
	nadded <- 0
	for (j in 1:length(rfo$par)) {
		ss <- substr(rfo$str, j, j)
		if (ss =='(') {
			nadded <- nadded+1
			pila[nadded] <- j
		} else if (ss == ')') {
			rfo$par[pila[nadded]] <- j
			rfo$par[j] <- pila[nadded]
			nadded <- nadded-1
		}
	}
	return(rfo)
}

proccessWindow <- function(rfo, min_length, min_bp, margin, trim_sequences) {
	reg <- data.frame(name={}, sseq={}, ibeg={}, iend={})
	idx <- getIndexs(rfo)
	for (j in 1:length(idx$b)) {
		sstr <- substr(rfo$str, idx$b[j], idx$e[j])
		if (idx$e[j]-idx$b[j] > min_length && BPNumber(sstr) >= min_bp) {
			sseq <- substr(rfo$seq, idx$b[j], idx$e[j])
			if (trim_sequences) {
				spar <- rfo$par[idx$b[j]:idx$e[j]] - idx$b[j] + 1
				sseq <- tryCatch({trimSequence(sseq, sstr, spar, min_bp, margin)},
						 error=function(cond) {return(sseq)})
			}
			sname <- attr(rfo$seq,'name')
			reg[j,'name'] <- gsub(pattern = ' ', replacement = '_', x = sname)
			reg[j,'sseq'] <- sseq
			reg[j,'ibeg'] <- idx$b[j]
			reg[j,'iend'] <- idx$e[j]
		}
	}
	if(nrow(reg) > 0 && sum(!is.na(reg$name) > 0)) {
		reg = reg[!is.na(reg$name),]
		reg = sprintf('>%s_stem-%d-%d\n%s\n',
			      reg$name, reg$ibeg, reg$iend, reg$sseq)
	} else {
		reg=""
	}
	return(reg)
}

trimSequence <- function(sseq, sstr, spar, minbp, margin){
	if (tail(gregexpr('\\(', sstr)[[1]], n=1) > gregexpr('\\)', sstr)[[1]][1])
		return(sseq)
	auxpr = regexpr('\\(\\.*\\)', sstr)
	mid <-  auxpr[1] + round( attr(auxpr, 'match.length') / 2 )
	flen <- rep(0, nchar(sseq))
	dst <- 0
	for (i in mid:nchar(sseq)) {
		flen[i] <- dst
		dst <- dst + as.numeric(substr(sstr, i, i) != '.')
	}
	dst <- 0
	for (i in seq(mid, 1, -1)) {
		flen[i] <- dst
		dst <- dst + as.numeric(substr(sstr, i, i) != '.')
	}
	p1 <- (flen-minbp-margin) > 0
	drb <- gregexpr('\\.*', sstr)[[1]]
	dre <- drb + attr(drb, 'match.length') - 1
	flen <- rep(0, nchar(sseq))
	for (i in 1:length(dre)) {
		flen[drb[i]:dre[i]] <- dre[i] - drb[i] + 1
	}
	p2 <- flen
	flen <- -abs(1:length(p2) - mid)
	p3 <- (1 - flen / min(flen))^2
	p <- p1 * p2 * p3
	if (any(p > 0)) {
		idx <- which.max(p)
		sbeg <- idx(which.min(abs(idx-mid)))[1]
		send <- max(spar(sbeg+1), spar(sbeg-1))
		if (sbeg < send) {
			sseq <- substr(sseq, sbeg, send)
		} else {
			sseq <- substr(sseq, send, sbeg)
		}
	}
	return(sseq)
}

getIndexs <- function(rfo) {
	auxpr <- gregexpr('\\)\\.*\\(', rfo$str)[[1]]
	if(auxpr[1] == -1) {
		marks <- c(regexpr('\\(', rfo$str),
			   tail(gregexpr('\\)', rfo$str)[[1]], n=1))
	} else {
		marks <- c(regexpr('\\(', rfo$str),
			   auxpr,
			   tail(gregexpr('\\)', rfo$str)[[1]], n=1),
			   attr(auxpr,"match.length") + auxpr - 1)
	}
	cuts <- rep(F,nchar(rfo$str))
	cuts[marks] <- T
	cuts[rfo$par[marks]] <- T
	idx <- data.frame(b = rep(-1, sum(cuts) / 2), e = rep(-1, sum(cuts) / 2))
	nadded <- 0
	last_cut <- -1
	for (i in which(cuts)){
		if (substr(rfo$str, i, i) == '(') {
			last_cut <- i
		} else if (substr(rfo$str, i, i) == ')' && last_cut>0) {
			nadded <- nadded+1
			idx$b[nadded] <- last_cut-1
			while (idx$b[nadded]>0 &&
			       substr(rfo$str, idx$b[nadded], idx$b[nadded]) == '.')
				idx$b[nadded] <- idx$b[nadded]-1
			idx$b[nadded] <- idx$b[nadded]+1
			idx$e[nadded] <- i+1
			while (idx$e[nadded] <= nchar(rfo$str) &&
			       substr(rfo$str, idx$e[nadded], idx$e[nadded]) == '.')
				idx$e[nadded] <- idx$e[nadded]+1
			idx$e[nadded] <- idx$e[nadded]-1
			last_cut <- -1
		}
	}
	return(idx[1:nadded,])
}
