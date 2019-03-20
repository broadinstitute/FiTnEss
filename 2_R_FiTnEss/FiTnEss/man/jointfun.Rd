\name{jointfun}
\alias{jointfun}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Construct distribution using different parameter combinations
%%  ~~function to do ... ~~
}
\description{
For each pair of parameters, this function constructs theoretical distributions which could be used in \code{my_cvm} function to identify optimized parameters.
}
\usage{
jointfun(data_input,lambda,lp,sigma,nta)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data_input}{vector of data input to help fit parameters}
  \item{lambda}{proportion of non-essential genes/whole genes}
  \item{lp}{parameter for negative binomial distribution}
  \item{sigma}{parameter for negative binomial distribution}
  \item{nta}{number of TA sites in gene}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
jointfun(data_vector,0.7,5,0,5,10)
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~fit_distribution }
