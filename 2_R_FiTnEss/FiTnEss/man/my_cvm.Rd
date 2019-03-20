\name{my_cvm}
\alias{my_cvm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Modified version of CVM function from package goftest
%%  ~~function to do ... ~~
}
\description{
Modified version of \code{cvm.test} function from package goftest, introduced \code{Heaviside} function to truncate lower 0.25 percent of data.

Parameters that construct a theoretical distribution with minimum CVM distance between actual data are optimized parameters.
}
\usage{
my_cvm(data_vector,"function",parameters,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data_input}{vector of actual data used to calculate cvm distance between theoretical distribution and actual one}
  \item{function}{function to construct theoretical distribution}
  \item{parameters}{parameters used to construct theoretical distribution}
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
my.cvm(data_vector,"jointfun",0.7,5,0.5,10)
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~cvm_distance }
