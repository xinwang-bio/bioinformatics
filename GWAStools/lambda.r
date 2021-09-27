		Args <- commandArgs()
		df=1
		data <- read.table(Args[6], header = T)
		data <- data[,4]
		data <- data[which(!is.na(data))]
		ntp <- round( 1 * length(data) )
		if ( max(data)<=1 ) {
			data <- qchisq(data, 1, lower.tail=FALSE)
		}
		data[which(abs(data)< -1)] <- NA
		data <- sort(data)
		ppoi <- ppoints(data)
		ppoi <- sort(qchisq(ppoi, df=df, lower.tail=FALSE))
		data <- data[1:ntp]
		ppoi <- ppoi[1:ntp]
		
		out <- list()
		s <- summary( lm(data~0+ppoi) )$coeff
		out$estimate <- s[1,1]
		out$se <- s[1,2]
		write.table(out$estimate, file = Args[7], row.names = FALSE, col.names = FALSE)

