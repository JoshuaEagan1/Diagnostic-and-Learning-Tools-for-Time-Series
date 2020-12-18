ar1.est<-function(x, true_ar, lag.max=NULL, pvalue=.05){
        z <- qnorm(1-(pvalue/2))
        acf_vals<-as.vector(TSA::acf(x, lag.max=lag.max, plot=F)$acf)
        lag.max<-length(acf_vals)
        k<-1:lag.max
        var<-((1 - (acf_vals^(2*k)))*(1 + (acf_vals^2))*((1 - (acf_vals^2))^-1))-(2*k*(acf_vals^(2*k)))
        confband <- matrix(NA,2,lag.max)
        for (i in 1:2){
                confband[i,] <- (acf_vals + c(-1,1)[i]*z*sqrt(var/length(x)))}
        matplot(k,t(confband),type="l",lty=c(2,2,1),col="blue",ylab=expression(rho[k]), ylim=c(min(confband)-.1, max(confband)+.1), xlab="lags")
        abline(0,0)
        true_vals<-as.vector(ARMAacf(ar=true_ar, lag.max=lag.max)[2:lag.max])
        for (i in 1:length(acf_vals)){
                points(i, true_vals[i])}
        title(paste0((1-pvalue)*100, "% ACF Confidence Intervals Against the True AR(1)"))
}