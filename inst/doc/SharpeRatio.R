
## ----'preamble', include=FALSE, warning=FALSE, message=FALSE-------------
library(knitr)

# set the knitr options ... for everyone!
# if you unset this, then vignette build bonks. oh, joy.
#opts_knit$set(progress=TRUE)
opts_knit$set(eval.after='fig.cap')
# for a package vignette, you do want to echo.
# opts_chunk$set(echo=FALSE,warning=FALSE,message=FALSE)
opts_chunk$set(warning=FALSE,message=FALSE)
#opts_chunk$set(results="asis")
opts_chunk$set(cache=TRUE,cache.path="cache/SharpeRatio")

#opts_chunk$set(fig.path="figure/",dev=c("pdf","cairo_ps"))
opts_chunk$set(fig.path="figure/SharpeRatio",dev=c("pdf"))
opts_chunk$set(fig.width=5,fig.height=4,dpi=64)

# doing this means that png files are made of figures;
# the savings is small, and it looks like shit:
#opts_chunk$set(fig.path="figure/",dev=c("png","pdf","cairo_ps"))
#opts_chunk$set(fig.width=4,fig.height=4)
# for figures? this is sweave-specific?
#opts_knit$set(eps=TRUE)

# this would be for figures:
#opts_chunk$set(out.width='.8\\textwidth')
# for text wrapping:
options(width=64,digits=2)
opts_chunk$set(size="small")
opts_chunk$set(tidy=TRUE,tidy.opts=list(width.cutoff=50,keep.blank.line=TRUE))

compile.time <- Sys.time()

# from the environment

# only recompute if FORCE_RECOMPUTE=True w/out case match.
FORCE_RECOMPUTE <- 
	(toupper(Sys.getenv('FORCE_RECOMPUTE',unset='False')) == "TRUE")

# compiler flags!

# not used yet
LONG.FORM <- FALSE

mc.resolution <- ifelse(LONG.FORM,1000,200)
mc.resolution <- max(mc.resolution,100)

library(quantmod)
options("getSymbols.warning4.0"=FALSE)

library(SharpeR)

gen_norm <- rnorm
lseq <- function(from,to,length.out) { 
	exp(seq(log(from),log(to),length.out = length.out))
}


## ----'generate_tbias',include=FALSE--------------------------------------
# the bias of the t
f_tbias <- function(n) { 
	sqrt((n-1) / 2) * exp(lgamma((n-2)/2) - lgamma((n-1)/2))
}
# approximate tbias
f_apx_tbias1 <- function(n) { 
	1 + (0.75 / n)
}
f_apx_tbias2 <- function(n) { 
	#1 + (0.75 / n) + (2 / n**2)
	1 + (0.75 / (n - 1)) + (32 / (25 * ((n-1) ** 2)))
}
n <- 12
tbias <- f_tbias(n)


## ----'example_usage_Fpower',include=FALSE--------------------------------
ex_p <- 30
ex_n <- 1000
ex_snr <- 1.5
dpy <- 253
ex_snr_d <- ex_snr / sqrt(dpy)
c <- ex_p / ex_n
cplus <- (c + ex_snr_d ** 2)
meanv <- sqrt(cplus / (1 - c))
stdv <- sqrt((1 / ex_n) * (ex_snr_d ** 4 + 2 * ex_snr_d ** 2 + c) / (2 * cplus * (1 - c) ** 2))


## ----'pos_cons_example', include=FALSE, warning=FALSE, message=FALSE-----
base.opt.1 <- function(rho,c1,c2) {
# compute v1, v2 that solve
#
#   min      v1^2 + v2^2
#   st.         v1 >= c1
#     -rho v1 + v2 >= c2 
# 
# return as a list?
# the prospective solutions are
# (0,0), (c1,0), (-2rho c2/(1+rho^2),c2(1-rho^2)/(1+rho^2)),
# and (c1,c2+rho*c1)
	rho2 <- rho^2
	oneprho2 <- 1 + rho2;
	solns <- matrix(c(0,0, c1,0, -2*rho*c2/(oneprho2),c2*(1-rho^2)/(oneprho2),
c1,c2+rho*c1),nrow=2)
	feasible <- solns[1,] >= c1
	solns <- solns[,feasible]
	feasible <- -rho * solns[1,] + solns[2,] >= c2
	solns <- solns[,feasible]
	solns <- matrix(solns,nrow=2)

	vals <- colSums(solns ^ 2)
	retval <- solns[,which.min(vals)]
	return(retval)
}
base.opt.2 <- function(rho,psi1,psi2) {
# compute lam1, lam2 that solve
# 
#   min     [lam1,lam2] * R^-1 [lam1,lam2]'
#   st.     R^-1 ([psi1,psi2]' + [lam1,lam2]') >= 0
#
#   where   R = [1,rho]
#               [rho,1]
#    so  R^-1 = [1, -rho]
#               [-rho, 1] / (1-rho^2)
#
# and then return R^-1 ([psi1,psi2]' + [lam1,lam2]')
	c1 <- rho * psi2 - psi1;
	c2 <- rho * psi1 - psi2;
	rv <- base.opt.1(rho,c1,c2)
	lam2 <- rv[2] / (1-rho^2)
	lam1 <- rv[1] + rho * lam2
	nmu1 <- psi1 + lam1
	nmu2 <- psi2 + lam2
	# deal with roundoff
	w1 <- max(0,nmu1 - rho * nmu2)
	w2 <- max(0,-rho * nmu1 + nmu2)
	return(c(w1,w2))
}
base.opt.3 <- function(mu1,mu2,sig1,sig12,sig2) {
# compute w1, w2 to solve
#
#   max   [w1,w2] * [mu1,mu2]'  / sqrt([w1,w2] * Sig * [w1,w2]')
#   st.   w1 >= 0
#         w2 >= 0
#
#  where  Sig = [sig1   sig12]
#               [sig12   sig2]
	rho <- sig12 / sqrt(sig1 * sig2)
	psi1 <- mu1 / sqrt(sig1)
	psi2 <- mu2 / sqrt(sig2)
	rv <- base.opt.2(rho,psi1,psi2)
	return(rv)
}
# test it
base.opt.3(1,0.5,1,0.4,1);
base.opt.3(1,0.5,1,0.5,1);
base.opt.3(1,0.5,1,0.6,1);
base.opt.3(1,0.5,1,0.8,1);

base.opt.3(1,0.5,1,-0.8,1);
opt.pos.T2 <- function(mu1,mu2,sig1,sig12,sig2) {
# compute the maximum value of 
#
#   max   [w1,w2] * [mu1,mu2]'  / sqrt([w1,w2] * Sig * [w1,w2]')
#   st.   w1 >= 0
#         w2 >= 0
	rv <- base.opt.3(mu1,mu2,sig1,sig12,sig2)
	w <- as.vector(rv)
	mu <- as.vector(c(mu1,mu2))
	Sig <- matrix(c(sig1,sig12,sig12,sig2),nrow=2)
	T <- (mu %*% w) / sqrt(t(w) %*% Sig %*% w)
	return(T^2)
}

	




## ----'leverage_foo',include=FALSE----------------------------------------
Elevi <- 1.5 
Elevis <- Elevi ** 2
Vlevi <- 1 / 12
ratlevi <- Vlevi / Elevis
vixlevi <- (0.4)^2


## ----'haircutting',fig.cap=paste("Q-Q plot of",n.sim,"simulated haircut values versus the approximation given by \\eqnref{hcut_apx} is shown.")----
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
	# for plugin estimator, estimate zeta.star
	asro <- sropt(z.s=sqrt(t(smeas$w) %*% smeas$mu),df1=p,df2=n)
	zeta.hat.s <- inference(asro,type="KRS")  # or 'MLE', 'unbiased'
	return(c(hcut,zeta.hat.s))
}

# set everything up
set.seed(as.integer(charToRaw("496509a9-dd90-4347-aee2-1de6d3635724")))
ope <- 253
LONG.FORM <- FALSE
n.sim <- if (LONG.FORM) 2048 else 512
n.stok <- if (LONG.FORM) 8 else 6
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


## ----'hcut_med',include=FALSE--------------------------------------------
medv.true <- median(hcuts)
med.snr.true <- zeta.s * (1 - medv.true)


