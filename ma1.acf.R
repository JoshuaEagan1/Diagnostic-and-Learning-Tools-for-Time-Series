ma1.acf<-function(y, true_ma, lag.max=NULL, pvalue=.05){
        acf.out<-as.vector(acf(y, lag.max=lag.max, plot=F)$acf) #getting the sample acf values
        lag.max<-length(acf.out)
        p1<-acf.out[1]
        var<-(1/length(y))*(1-(3*p1^2)+(4*p1^4)) #calculating the variance for roh(1)
        var_oth<-(1/length(y))*(1+(2*p1^2)) #calculating the variance for roh > 1
        var<-append(var, rep(var_oth, lag.max-1))
        z <- qnorm(1-(pvalue/2))
        k <- 1:lag.max
        true_vals<-ARMAacf(ma=true_ma, lag.max = lag.max)[1:lag.max+1] #finding the theoretical autocorrelation function
        confband <- matrix(NA,2,lag.max)
        for (i in 1:2){
                confband[i,] <- (acf.out + c(-1,1)[i]*z*sqrt(var))
        }
        matplot(k,t(confband),type="l",lty=c(2,2,1),col="blue",ylab=expression(rho[k]), ylim=c(min(confband)-.1, max(confband)+.1))
        abline(0,0)
        for (i in 1:lag.max){
                points(i, true_vals[i])}
        title(paste0((1-pvalue)*100, "% ACF Confidence Intervals Against the True MA(1)"))
}