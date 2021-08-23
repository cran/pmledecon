### R code from vignette source 'pmledecon.Rnw'

###################################################
### code chunk number 1: pmledecon.Rnw:45-60
###################################################
library(pmledecon)
sz=esz=30
set.seed(45217)
truth=rnorm(sz,0,1)
error=rnorm(esz,0,2)
ob=truth+error
error1=rnorm(esz,0,2)
est=pmledecon(ob,error1)
plot(density(ob,n=1000),col="red",lwd=2,lty=3,type="l",ylim=c(0,0.4),
     xlab="",main="unknown error")
lines(seq(-10,10,length.out=1000),dnorm(seq(-10,10,length.out=1000),0,1),
      lwd=2,lty=4,col="green")
lines(est$sup,est$f,lwd=2) 
legend("topright", lty=c(1,4,3),col=c("black","green","red"),lwd=2, 
       legend=c("Pmle","true density","kernel density of data"))


