eacf_pic<-function(x, ar.max=7, ma.max = 13, text_pvals=F){

        n<-length(x)

        eacf_pvalue<-function(mat1, n, ar.max =ar.max, ma.max = ma.max){ 
                tstat<-matrix(c(rep(1,(ar.max+1)*(ma.max+1))), nrow = ar.max+1, ncol = ma.max+1)
                se<-matrix(c(rep(1,(ar.max+1)*(ma.max+1))), nrow = ar.max+1, ncol = ma.max+1)  
                pvalue<-matrix(c(rep(1,(ar.max+1)*(ma.max+1))), nrow = ar.max+1, ncol = ma.max+1)
                symbol<-matrix(c(rep(1,(ar.max+1)*(ma.max+1))), nrow = ar.max+1, ncol = ma.max+1)
                for (i in 1:(ar.max+1))
                        {
                        for (j in 1:(ma.max+1)) { 
                                se[i,j]<- 1/sqrt((n-i-j+2))
                                tstat[i,j]<-abs( mat1[i,j]/se[i,j])  
                                pvalue[i,j]<-2*(1-pnorm(tstat[i,j]))
                                }
                        }
                symbol<-ifelse(pvalue>.05,"o","x")
                rownames(pvalue)<-0:ar.max
                pvalue<-data.frame(pvalue)
                names(pvalue)<-0:ma.max
                pvalue
        }
        
        eacf<-function (z, ar.max = 7, ma.max = 13) 
        {
                lag1 <- function(z, lag = 1) {
                        c(rep(NA, lag), z[1:(length(z) - lag)])
                }
                reupm <- function(m1, nrow, ncol) {
                        k <- ncol - 1
                        m2 <- NULL
                        for (i in 1:k) {
                                i1 <- i + 1
                                work <- lag1(m1[, i])
                                work[1] <- -1
                                temp <- m1[, i1] - work * m1[i1, i1]/m1[i, i]
                                temp[i1] <- 0
                                m2 <- cbind(m2, temp)
                        }
                        m2
                }
                ceascf <- function(m, cov1, nar, ncol, count, ncov, z, zm) {
                        result <- 0 * seq(1, nar + 1)
                        result[1] <- cov1[ncov + count]
                        for (i in 1:nar) {
                                temp <- cbind(z[-(1:i)], zm[-(1:i), 1:i]) %*% c(1, 
                                                                                -m[1:i, i])
                                result[i + 1] <- acf(temp, plot = FALSE, lag.max = count, 
                                                     drop.lag.0 = FALSE)$acf[count + 1]
                        }
                        result
                }
                ar.max <- ar.max + 1
                ma.max <- ma.max + 1
                nar <- ar.max - 1
                nma <- ma.max
                ncov <- nar + nma + 2
                nrow <- nar + nma + 1
                ncol <- nrow - 1
                z <- z - mean(z)
                zm <- NULL
                for (i in 1:nar) zm <- cbind(zm, lag1(z, lag = i))
                cov1 <- acf(z, lag.max = ncov, plot = FALSE, drop.lag.0 = FALSE)$acf
                cov1 <- c(rev(cov1[-1]), cov1)
                ncov <- ncov + 1
                m1 <- matrix(0, ncol = ncol, nrow = nrow)
                for (i in 1:ncol) m1[1:i, i] <- ar.ols(z, order.max = i, 
                                                       aic = FALSE, demean = FALSE, intercept = FALSE)$ar
                eacfm <- NULL
                for (i in 1:nma) {
                        m2 <- reupm(m1 = m1, nrow = nrow, ncol = ncol)
                        ncol <- ncol - 1
                        eacfm <- cbind(eacfm, ceascf(m2, cov1, nar, ncol, i, 
                                                     ncov, z, zm))
                        m1 <- m2
                }
                work <- 1:(nar + 1)
                work <- length(z) - work + 1
                symbol <- NULL
                for (i in 1:nma) {
                        work <- work - 1
                        symbol <- cbind(symbol, ifelse(abs(eacfm[, i]) > 2/work^0.5, 
                                                       "x", "o"))
                }
                rownames(symbol) <- 0:(ar.max - 1)
                colnames(symbol) <- 0:(ma.max - 1)
                invisible(list(eacf = eacfm, ar.max = ar.max, ma.ma = ma.max, 
                               symbol = symbol))
        }

        out=invisible(eacf(x, ar.max=ar.max, ma.max=ma.max))
        epo<-eacf_pvalue(out$eacf, n, ar.max=ar.max, ma.max=ma.max)
        m <- melt(epo)
        m$pvalue<-round(m$value, 2)
        m<-cbind(m, rep(0:ar.max, ma.max+1))
        names(m)<-c("rows","value","pvalue", "columns")
        Spectrum<-c()
        for(i in 1:nrow(m)){
                if(m$value[i]<=.001){
                        Spectrum[i]<-"P </= .001"
                } else if(m$value[i]<=.01){
                        Spectrum[i]<-".001 < P </= .01"
                } else if(m$value[i]<=.05){
                        Spectrum[i]<-".01 < P </= .05"
                } else if(m$value[i]<=.1){
                        Spectrum[i]<-".05 < P </= .10"
                } else if(m$value[i]<=.2){
                        Spectrum[i]<-".10 < P </= .2"
                } else if(m$value[i]>.2){
                        Spectrum[i]<-".P > .02"
                }
        }
        Spectrum <- factor(Spectrum, levels = c("P </= .001", ".001 < P </= .01",
                                                ".01 < P </= .05",  ".05 < P </= .10",  ".10 < P </= .20",  ".P > .20"))
        m<-cbind(m, Spectrum)
        rm(Spectrum)
        
        p<-ggplot(m, aes(x = rows, y = reorder(columns, desc(columns)))) +
                geom_tile(aes(fill = Spectrum), colour="black", show.legend=T) +
                scale_fill_manual(values= c("#EF3B2C", "#FB6A4A", "#FC9272", "#BDBDBD", "#737373", "#525252"))+
                scale_x_discrete(position = "top")+
                labs(title="Extended Autocorrelation Function (EACF)", caption = "Smaller P-values are shown in red. P-values above .05 are shown in grey.")+
                xlab("Moving Average (MA)")+
                ylab("Autoregressive (AR)")+
                theme_bw()
        
        if(text_pvals){
                p + geom_text(aes(label=pvalue), color="black", size=rel(4.5))
        } else {
                p
        }
}



