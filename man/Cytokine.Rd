% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/StatComp21019R.R
\name{Cytokine}
\alias{Cytokine}
\title{cytokine Response}
\usage{
Cytokine(
  t,
  d,
  t0 = c(1, 5, 9, 13, 17, 21, 25),
  EMAX = 359000,
  EC50 = 10000,
  H = 0.92 * 0.9,
  IMAX = 0.995,
  IC = 18200,
  KDEG = 3.6,
  K = 2.83
)
}
\arguments{
\item{t}{the time}

\item{d}{the dose of drug}

\item{t0}{Administration time}

\item{EMAX}{Maximum cytokine release rate}

\item{EC50}{Drug exposure for half-maximum release}

\item{H}{Hill coefficient for cytokine release}

\item{IMAX}{Maximum inhibition of cytokine release}

\item{IC}{Cytokine exposure for half-maximum inhibition}

\item{KDEG}{Degradation rate for cytokine}

\item{K}{Priming factor for cytokine release}
}
\value{
A list of the cytokine response , its derivatives and its AUC.
}
\description{
A cytokine response function for the PD model
}
\examples{
\dontrun{
t<-seq(from=0,to=35,by=0.02)
plot(t,Cytokine(t=t,d=rep(25,7))$EE,type="l",lty=2,ylab="E(t)")
abline(h=max(Cytokine(t=t,d=rep(25,7))$EE),lty=2)
lines(t,Cytokine(t=t,d=c(1,5,10,rep(25,4)))$EE,type="l")
abline(h=max(Cytokine(t=t,d=c(1,5,10,rep(25,4)))$EE),lty=1)
}
}
