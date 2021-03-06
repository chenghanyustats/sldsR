% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/update_x.R
\name{update_x}
\alias{update_x}
\title{Update continous latent states using FFBS algorithm}
\usage{
update_x(m0, C0, y, v, w, phi, st, k2)
}
\arguments{
\item{m0}{mean of the normal prior}

\item{C0}{(co)variance of the normal prior}

\item{y}{observation vector}

\item{v}{mcmc sample of observational variance}

\item{w}{mcmc sample of AR(1) coefficient}

\item{phi}{mcmc sample of AR(1) coefficient}

\item{st}{mcmc sample of \code{s_t}, the discrete hidden states}

\item{k2}{observational variance when \code{s} is 2}
}
\value{
A vector of continous latent state
}
\description{
Sample new continous latent states using forward filtering
    backward sampling (FFBS) algorithm with the normal prior \code{N(m0, C0)}
}
