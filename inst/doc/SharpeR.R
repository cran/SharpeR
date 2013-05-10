
## @knitr 'preamble', include=FALSE, warning=FALSE, message=FALSE
library(knitr)

# set the knitr options ... for everyone!
# if you unset this, then vignette build bonks. oh, joy.
#opts_knit$set(progress=TRUE)
opts_knit$set(eval.after='fig.cap')
# for a package vignette, you do want to echo.
# opts_chunk$set(echo=FALSE,warning=FALSE,message=FALSE)
opts_chunk$set(warning=FALSE,message=FALSE)
#opts_chunk$set(results="asis")
opts_chunk$set(cache=TRUE,cache.path="cache/")
opts_chunk$set(fig.path="figure/",dev=c("pdf","cairo_ps"))
# for figures? this is sweave-specific?
#opts_knit$set(eps=TRUE)

# this would be for figures:
#opts_chunk$set(out.width='.8\\textwidth')
# for text wrapping:
options(width=64,digits=2)
opts_chunk$set(size="small")
opts_chunk$set(tidy=TRUE,tidy.opts=list(width.cutoff=50,keep.blank.line=TRUE))

compile.time <- Sys.time()

# compiler flags!

# not used yet
LONG.FORM <- TRUE

library(quantmod)
options("getSymbols.warning4.0"=FALSE)


## @knitr 'babysteps'
library(SharpeR)
# suppose you computed the Sharpe of your strategy to
# be 1.3 / sqrt(yr), based on 1200 daily observations.
# store them as follows:
my.sr <- sr(sr=1.3,df=1200-1,ope=252,epoch="yr")
print(my.sr)
# multiple strategies can be tracked as well.
# one can attach names to them.
srstats <- c(0.5,1.2,0.6)
dim(srstats) <- c(3,1)
rownames(srstats) <- c("strat. A","strat. B","benchmark")
my.sr <- sr(srstats,df=1200-1,ope=252,epoch="yr")
print(my.sr)


## @knitr 'showoff'
set.seed(as.integer(charToRaw("set the seed")))
# Sharpe's 'model': just given a bunch of returns.
returns <- rnorm(253*8,mean=3e-4,sd=1e-2)
asr <- as.sr(returns,ope=253,epoch="yr")
print(asr)
# a data.frame with a single strategy
asr <- as.sr(data.frame(my.strategy=returns),ope=253,epoch="yr")
print(asr)


## @knitr 'more_data_frame'
# a data.frame with multiple strategies
asr <- as.sr(data.frame(strat1=rnorm(253*8),strat2=rnorm(253*8,mean=4e-4,sd=1e-2)),
	ope=253,epoch="yr")
print(asr)


## @knitr 'some_stocks'
library(quantmod)
get.ret <- function(sym,warnings=FALSE,...) {
	# getSymbols.yahoo will barf sometimes; do a trycatch
	ntry <- 3
  trynum <- 0
	while (!exists("OHCLV") && (trynum < ntry)) {
		trynum <- trynum + 1
		try(OHLCV <- getSymbols(sym,auto.assign=FALSE,warnings=warnings,...),silent=TRUE)
  }
	adj.names <- paste(c(sym,"Adjusted"),collapse=".",sep="")
	lrets <- diff(log(OHLCV[,adj.names]))
	#chop first
	lrets[-1,]
}
get.rets <- function(syms,...) { some.rets <- do.call("cbind",lapply(syms,get.ret,...)) }
# quantmod::periodReturn does not deal properly with multiple
# columns, and the straightforward apply(mtms,2,periodReturn) barfs
my.periodReturn <- function(mtms,...) {
	per.rets <- do.call(cbind,lapply(mtms,
		function(x) {
			retv <- periodReturn(x,...)
			colnames(retv) <- colnames(x)
			return(retv) 
		}))
}
some.rets <- get.rets(c("AAPL","IBM","A","C"),from="2001-01-01")
print(as.sr(some.rets))


## @knitr 'reannualize'
some.rets <- get.rets(c("XOM"),from="2001-01-01")
yearly <- as.sr(some.rets)
monthly <- reannualize(yearly,new.ope=21,new.epoch="mo.")
print(yearly)
# significance should be the same, but units changed.
print(monthly)


## @knitr 'BRKB'
# get the returns (see above for the function)
some.rets <- get.rets(c("BRK.B","SPY"),from="1996-05-09")
# make them monthly:
mo.rets <- my.periodReturn(exp(cumsum(some.rets)),period='monthly',type='arithmetic')
# look at both of them together:
both.sr <- as.sr(mo.rets)
print(both.sr)
# confindence intervals on the Sharpe:
print(confint(both.sr))
# perform a CAPM attribution, using SPY as 'the market'
linmod <- lm(BRK.B.Adjusted ~ SPY.Adjusted,data=mo.rets)
# convert attribution model to Sharpe
CAPM.sr <- as.sr(linmod,ope=both.sr$ope,epoch="yr")
# statistical significance does not change (though note the sr summary
# prints a 1-sided p-value)
print(summary(linmod))
print(CAPM.sr)
# the confidence intervals tell the same story, but in different units:
print(confint(linmod,'(Intercept)'))
print(confint(CAPM.sr))


## @knitr 'SPDRcheck'
# get the sector 'spiders'
some.rets <- get.rets(c("XLY","XLE","XLP","XLF","XLV","XLI","XLB","XLK","XLU"),from="1999-01-01")
# make them monthly:
mo.rets <- my.periodReturn(exp(cumsum(some.rets)),period='monthly',type='arithmetic')
# one-sample test on utilities:
XLU.monthly <- mo.rets[,"XLU.Adjusted"]
print(sr_test(XLU.monthly),alternative="two.sided")

# test for equality of Sharpe among the different spiders
print(sr_equality_test(some.rets))
# perform a paired two-sample test via sr_test:
XLF.monthly <- mo.rets[,"XLF.Adjusted"]
print(sr_test(x=XLU.monthly,y=XLF.monthly,ope=12,paired=TRUE))


## @knitr 'power_thing',fig.cap="The percent error of the power mnemonic $e\\approx\\ssiz \\psnrsq$ is plotted versus \\psnr."
ope <- 253
 zetas <- seq(0.1,2.5,length.out=101)
ssizes <- sapply(zetas,function(zed) { 
 x <- power.sr_test(n=NULL,zeta=zed,sig.level=0.05,power=0.5,ope=ope)
 x$n / ope 
})
plot(zetas,100 * ((exp(1) / zetas^2) - ssizes)/ssizes, ylab="error in mnemonic rule (as %)")


## @knitr 'sobering',include=FALSE
foo.power <- power.sr_test(n=253,zeta=NULL,sig.level=0.05,power=0.5,ope=253)


## @knitr 'sropt_basics'
# from a matrix object:
ope <- 253
n.stok <- 7
n.yr <- 8
# somewhat unrealistic: the returns are independent.
some.rets <- matrix(rnorm(n.yr * ope * n.stok),ncol=n.stok)
asro <- as.sropt(some.rets,ope=ope)
print(asro)
# from an xts object
some.rets <- get.rets(c("IBM","AAPL","A","C","SPY","XOM"),from="2001-01-01")
asro <- as.sropt(some.rets)
print(asro)


## @knitr 'sropt_estim'
# confidence intervals:
print(confint(asro,level.lo=0.05,level.hi=1))
# estimation
print(inference(asro,type="KRS"))
print(inference(asro,type="MLE"))


## @knitr 'MLE_rule'
ope <- 253
zeta.s <- 0.8
n.check <- 1000
df1 <- 10
df2 <- 6 * ope
rvs <- rsropt(n.check,df1,df2,zeta.s,ope,drag=0)
roll.own <- sropt(z.s=rvs,df1,df2,drag=0,ope=ope,epoch="yr")
MLEs <- inference(roll.own,type="MLE")
zerMLE <- MLEs <= 0
crit.value <- 0.5 * (max(rvs[zerMLE])^2 + min(rvs[!zerMLE])^2)
aspect.ratio <- df1 / (df2 / ope)
cat(sprintf("empirical cutoff for zero MLE is %2.2f yr^{-1}\n", crit.value))
cat(sprintf("the aspect ratio is %2.2f yr^{-1}\n",aspect.ratio))


## @knitr 'haircutting',fig.cap=paste("Q-Q plot of",n.sim,"simulated haircut values versus the approximation given by \\eqnref{hcut_apx} is shown.")
require(MASS)

# simple markowitz.
simple.marko <- function(rets) {
	mu.hat <- as.vector(apply(rets,MARGIN=2,mean,na.rm=TRUE))
	Sig.hat <- cov(rets)
	w.opt <- solve(Sig.hat,mu.hat)
	retval <- list('mu'=mu.hat,'sig'=Sig.hat,'w'=w.opt)
	return(retval)
}
# make multivariate pop. & sample w/ given zeta.star
gen.pop <- function(n,p,zeta.s=0) {
	true.mu <- matrix(rnorm(p),ncol=p)
	#generate an SPD population covariance. a hack.
	xser <- matrix(rnorm(p*(p + 100)),ncol=p)
	true.Sig <- t(xser) %*% xser
	pre.sr <- sqrt(true.mu %*% solve(true.Sig,t(true.mu)))
	#scale down the sample mean to match the zeta.s
	true.mu <- (zeta.s/pre.sr[1]) * true.mu 
  X <- mvrnorm(n=n,mu=true.mu,Sigma=true.Sig)
	retval = list('X'=X,'mu'=true.mu,'sig'=true.Sig,'SNR'=zeta.s)
	return(retval)
}
# a single simulation
sample.haircut <- function(n,p,...) {
	popX <- gen.pop(n,p,...)
	smeas <- simple.marko(popX$X)
	# I have got to figure out how to deal with vectors...
	ssnr <- (t(smeas$w) %*% t(popX$mu)) / sqrt(t(smeas$w) %*% popX$sig %*% smeas$w)
	hcut <- 1 - (ssnr / popX$SNR)
	# to compute the plugin estimator, estimate zeta.star
	asro <- sropt(z.s=sqrt(t(smeas$w) %*% smeas$mu),df1=p,df2=n)
	zeta.hat.s <- inference(asro,type="KRS")  # or 'MLE', 'unbiased'
	return(c(hcut,zeta.hat.s))
}

# set everything up
set.seed(as.integer(charToRaw("496509a9-dd90-4347-aee2-1de6d3635724")))
ope <- 253
n.sim <- if (LONG.FORM) 2048 else 512
n.stok <- 8
n.yr <- 4
n.obs <- ceiling(ope * n.yr)
zeta.s <- 1.20 / sqrt(ope)   # optimal SNR, in daily units

# run some experiments
system.time(experiments <- replicate(n.sim,sample.haircut(n.obs,n.stok,zeta.s)))
hcuts <- experiments[1,]
print(summary(hcuts))
# haircut approximation in the equation above
qhcut <- function(p, df1, df2, zeta.s, lower.tail=TRUE) {
	1 - sin(atan((1/sqrt(df1-1)) * qt(p,df=df1-1,ncp=sqrt(df2)*zeta.s,lower.tail=!lower.tail)))
}
# if you wanted to look at how bad the plug-in estimator is, then
# uncomment the following (you are warned):
# zeta.hat.s <- experiments[2,];                                   
# qqplot(qhcut(ppoints(length(hcuts)),n.stok,n.obs,zeta.hat.s),hcuts,
# 			 xlab = "Theoretical Approximate Quantiles", ylab = "Sample Quantiles");
# qqline(hcuts,datax=FALSE,distribution = function(p) { qhcut(p,n.stok,n.obs,zeta.hat.s) },
# 			 col=2)

# qqplot;
qqplot(qhcut(ppoints(length(hcuts)),n.stok,n.obs,zeta.s),hcuts,
			 xlab = "Theoretical Approximate Quantiles", ylab = "Sample Quantiles")
qqline(hcuts,datax=FALSE,distribution = function(p) { qhcut(p,n.stok,n.obs,zeta.s) },
			 col=2)


## @knitr 'hcut_med',include=FALSE
medv.true <- median(hcuts)
med.snr.true <- zeta.s * (1 - medv.true)


## @knitr 'estimate_overfit',fig.cap=paste("Q-Q plot of",n.sim,"achieved optimal \\txtSR values from brute force search over both windows of a Moving Average Crossover under the null of driftless log returns with zero autocorrelation versus the approximation by a 2-parameter optimal \\txtSR distribution is shown.")
require(TTR)
# brute force search two window MAC
brute.force <- function(lrets,rrets=exp(lrets)-1,win1,win2=win1) {
	mtms <- c(1,exp(cumsum(lrets)))  # prepend a 1.
  # do all the SMAs;
  SMA1 <- sapply(win1,function(n) { SMA(mtms,n=n) }) 
  symmetric <- missing(win2)
  if (!symmetric)
  	SMA2 <- sapply(win2,function(n) { SMA(mtms,n=n) }) 

  mwin <- max(c(win1,win2))
  zeds <- matrix(NaN,nrow=length(win1),ncol=length(win2))
	upb <- if (symmetric) length(win1) - 1 else length(win1)
	# 2FIX: vectorize this!
	for (iidx in 1:upb) {
		SM1 <- SMA1[,iidx]
		lob <- if (symmetric) iidx + 1 else 1
		for (jidx in lob:length(win2)) {
			SM2 <- if (symmetric) SMA1[,jidx] else SMA2[,jidx]
			trades <- sign(SM1 - SM2)
			dum.bt <- trades[mwin:(length(trades)-1)] * rrets[mwin:length(rrets)]  # braindead backtest.
			mysr <- as.sr(dum.bt)
			zeds[iidx,jidx] <- mysr$sr
			if (symmetric) zeds[jidx,iidx] <- - zeds[iidx,jidx]  # abuse symmetry of arithmetic returns
		}
	}
	retv <- max(zeds,na.rm=TRUE) 
	return(retv)
}
# simulate one.
sim.one <- function(nbt,win1,...) {
	lrets <- rnorm(nbt+max(win1),sd=0.01)
	retv <- brute.force(lrets,win1=win1,...)
	return(retv)
}
# set everything up
set.seed(as.integer(charToRaw("e23769f4-94f8-4c36-bca1-28c48c49b4fb")))
ope <- 253
n.yr <- 4
n.obs <- ceiling(ope * n.yr)
n.sim <- if (LONG.FORM) 2048 else 512
win1 <- c(2,4,8,16,32,64,128,256)

# run them
system.time(max.zeds <- replicate(n.sim,sim.one(n.obs,win1)))
# qqplot;
qqplot(qsropt(ppoints(length(max.zeds)),df1=2,df2=n.obs),max.zeds,
			 xlab = "Theoretical Approximate Quantiles", ylab = "Sample Quantiles")
qqline(max.zeds,datax=FALSE,distribution = function(p) { qsropt(p,df1=2,df2=n.obs) },
			 col=2)


## @knitr 'now_on_spy'
# is MAC on SPY significant?
SPY.lret <- get.ret('SPY',from="1995-01-01")
mysr <- as.sr(SPY.lret)  # just to get the ope
print(mysr)
# try a whole lot of windows:
win1 <- seq(4,204,by=10)
zeds <- brute.force(SPY.lret,win1=win1)
asro <- sropt(z.s=zeds,df1=2,df2=length(SPY.lret) - max(win1),ope=mysr$ope)
print(asro)
print(inference(asro,type="KRS"))


