% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tests.r
\name{sr_equality_test}
\alias{sr_equality_test}
\title{Paired test for equality of Sharpe ratio}
\usage{
sr_equality_test(X,type=c("chisq","F","t"),
                 alternative=c("two.sided","less","greater"),
                 contrasts=NULL,
                 vcov.func=vcov)
}
\arguments{
\item{X}{an \eqn{n \times p}{n x p} matrix of paired observations.}

\item{type}{which approximation to use. \code{"chisq"} is preferred when
the returns are non-normal, but the approximation is asymptotic.
the \code{"t"} test is only supported when \eqn{k = 1}{k = 1}.}

\item{alternative}{a character string specifying the alternative hypothesis,
must be one of \code{"two.sided"} (default), \code{"greater"} or
\code{"less"}. You can specify just the initial letter.
This is only relevant for the \code{"t"} test.
\code{"greater"} corresponds to \eqn{H_a: E s > 0}{Ha: E s > 0}.}

\item{contrasts}{an \eqn{k \times p}{k x p} matrix of the contrasts}

\item{vcov.func}{a function which takes a model of class lm (one of
the form x ~ 1), and produces a variance-covariance matrix.
The default is \code{\link{vcov}}, which produces a 'vanilla'
estimate of covariance. Other sensible options are
\code{vcovHAC} from the \code{sandwich} package.}
}
\value{
Object of class \code{htest}, a list of the test statistic,
the size of \code{X}, and the \code{method} noted.
}
\description{
Performs a hypothesis test of equality of Sharpe ratios of p assets
given paired observations.
}
\details{
Given \eqn{n} \emph{i.i.d.} observations of the excess returns of
\eqn{p} strategies, we test
\deqn{H_0: \frac{\mu_i}{\sigma_i} = \frac{\mu_j}{\sigma_j}, 1 \le i < j \le p}{H0: sr1 = sr2 = ...}
using the method of Wright, et. al. 

More generally, a matrix of constrasts, \eqn{E}{E} can be given, and we can
test
\deqn{H_0: E s = 0,}{H0: E s = 0,}
where \eqn{s}{s} is the vector of Sharpe ratios of the \eqn{p} strategies.

When \eqn{E}{E} consists of a single row (a single contrast), as is the
case when the default contrasts are used and only two strategies are
compared, then an approximate t-test can be performed against the
alternative hypothesis \eqn{H_a: E s > 0}{Ha: E s > 0}

Both chi-squared and F- approximations are supported; the former is
described by Wright. \emph{et. al.}, the latter by Leung and Wong.

See \sQuote{The Sharpe Ratio: Statistics and Applications},
section 3.3.1.
}
\examples{
# under the null 
set.seed(1234)
rv <- sr_equality_test(matrix(rnorm(500*5),ncol=5))

# under the alternative (but with identity covariance)
ope <- 253
nyr <- 10
nco <- 5
set.seed(909)
rets <- 0.01 * sapply(seq(0,1.7/sqrt(ope),length.out=nco),
  function(mu) { rnorm(ope*nyr,mean=mu,sd=1) })
rv <- sr_equality_test(rets)

# using real data
if (require(xts)) {
 data(stock_returns)
 pvs <- sr_equality_test(stock_returns)
}

# test for uniformity
pvs <- replicate(1024,{ x <- sr_equality_test(matrix(rnorm(400*5),400,5),type="chisq")
                       x$p.value })
plot(ecdf(pvs))
abline(0,1,col='red') 

\donttest{
if (require(sandwich)) {
  set.seed(as.integer(charToRaw("0b2fd4e9-3bdf-4e3e-9c75-25c6d18c331f")))
  n.manifest <- 10
  n.latent <- 4
  n.day <- 1024
  snr <- 0.95
  la_A <- matrix(rnorm(n.day*n.latent),ncol=n.latent)
  la_B <- matrix(runif(n.latent*n.manifest),ncol=n.manifest)
  latent.rets <- la_A \%*\% la_B
  noise.rets <- matrix(rnorm(n.day*n.manifest),ncol=n.manifest)
  some.rets <- snr * latent.rets + sqrt(1-snr^2) * noise.rets
  # naive vcov
  pvs0 <- sr_equality_test(some.rets)
  # HAC vcov
  pvs1 <- sr_equality_test(some.rets,vcov.func=vcovHAC)
  # more elaborately:
  pvs <- sr_equality_test(some.rets,vcov.func=function(amod) {
	vcovHAC(amod,prewhite=TRUE) })
}
}

}
\references{
Sharpe, William F. "Mutual fund performance." Journal of business (1966): 119-138.
\url{https://ideas.repec.org/a/ucp/jnlbus/v39y1965p119.html}

Wright, J. A., Yam, S. C. P., and Yung, S. P. "A note on the test for the
equality of multiple Sharpe ratios and its application on the evaluation
of iShares." J. Risk. to appear. 
\url{https://www.risk.net/journal-risk/2340067/test-equality-multiple-sharpe-ratios}

Leung, P.-L., and Wong, W.-K. "On testing the equality of multiple Sharpe ratios, with 
application on the evaluation of iShares." J. Risk 10, no. 3 (2008): 15--30.
\url{https://papers.ssrn.com/sol3/papers.cfm?abstract_id=907270}

Memmel, C. "Performance hypothesis testing with the Sharpe ratio." Finance
Letters 1 (2003): 21--23.

Ledoit, O., and Wolf, M. "Robust performance hypothesis testing with the 
Sharpe ratio." Journal of Empirical Finance 15, no. 5 (2008): 850-859.
\doi{10.1016/j.jempfin.2008.03.002}

Lo, Andrew W. "The statistics of Sharpe ratios." Financial Analysts Journal 58, no. 4 
(2002): 36-52. \url{https://www.ssrn.com/paper=377260}

Pav, S. E. "The Sharpe Ratio: Statistics and Applications." CRC Press, 2021.
}
\seealso{
\code{\link{sr_test}}

Other sr: 
\code{\link{as.sr}()},
\code{\link{confint.sr}()},
\code{\link{dsr}()},
\code{\link{is.sr}()},
\code{\link{plambdap}()},
\code{\link{power.sr_test}()},
\code{\link{predint}()},
\code{\link{print.sr}()},
\code{\link{reannualize}()},
\code{\link{se}()},
\code{\link{sr_test}()},
\code{\link{sr_unpaired_test}()},
\code{\link{sr_vcov}()},
\code{\link{sr}},
\code{\link{summary.sr}}
}
\author{
Steven E. Pav \email{shabbychef@gmail.com}
}
\concept{sr}
\keyword{htest}
