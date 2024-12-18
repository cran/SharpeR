% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/estimation.r
\name{asnr_confint}
\alias{asnr_confint}
\alias{asnr_confint.sropt}
\alias{asnr_confint.del_sropt}
\title{Confidence intervals on achieved SnR}
\usage{
asnr_confint(z.s, level = 0.95, level.lo = (1 - level), level.hi = 1)

\method{asnr_confint}{sropt}(z.s, level = 0.95, level.lo = (1 - level), level.hi = 1)

\method{asnr_confint}{del_sropt}(z.s, level = 0.95, level.lo = (1 - level), level.hi = 1)
}
\arguments{
\item{z.s}{an object of type \code{sropt}, or \code{del_sropt}}

\item{level}{the confidence level required.}

\item{level.lo}{the lower confidence level required.}

\item{level.hi}{the upper confidence level required.}
}
\value{
an estimate of the non-centrality parameter, which is
the maximal population Sharpe ratio.
}
\description{
Computes approximate bounds on the achieved signal-noise ratio of the
Markowitz portfolio built on sample data.
}
\details{
Provides an approximate bound on the achieved Signal-noise ratio of the
sample Markowitz portfolio. That is if \eqn{\mu}{mu} and \eqn{\Sigma}{Sigma}
are the unknown mean and covariance of returns, and \eqn{w}{w} is the 
sample Markowitz portfolio, then the probability that 
\deqn{w^{\top}\mu / \sqrt{w^{\top}\Sigma w} \ge b}{w'mu/sqrt(w'Sigma w) >= b}
is the given probability level.
See section 8.3.1 of \sQuote{The Sharpe Ratio: Statistics and Applications}.
Plugs in the \eqn{\delta_2}{delta_2} estimator.
}
\examples{
# generate some sropts
nfac <- 3
nyr <- 5
ope <- 253
# simulations with no covariance structure.
# under the null:
set.seed(as.integer(charToRaw("determinstic")))
Returns <- matrix(rnorm(ope*nyr*nfac,mean=0,sd=0.0125),ncol=nfac)
asro <- as.sropt(Returns,drag=0,ope=ope)
asnr_confint(asro)

# for del_sropt:
nfac <- 5
nyr <- 10
ope <- 253
set.seed(as.integer(charToRaw("fix seed")))
Returns <- matrix(rnorm(ope*nyr*nfac,mean=0.0005,sd=0.0125),ncol=nfac)
# hedge out the first one:
G <- matrix(diag(nfac)[1,],nrow=1)
asro <- as.del_sropt(Returns,G,drag=0,ope=ope)
asnr_confint(asro)

}
\references{
Pav, S. E. "The Sharpe Ratio: Statistics and Applications." CRC Press, 2021.

Pav, S. E. "Inference on achieved signal noise ratio."
2020 \url{https://arxiv.org/abs/2005.06171}
}
\seealso{
Other sropt Hotelling: 
\code{\link{inference}()},
\code{\link{sric}()}
}
\author{
Steven E. Pav \email{shabbychef@gmail.com}
}
\concept{sropt Hotelling}
