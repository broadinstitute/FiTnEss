\name{callessfun}
\alias{callessfun}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Identify essential genes using calculated parameters
%%  ~~function to do ... ~~
}
\description{
Call essential genes based on calculated parameters from function \code{calcparafun}
}
\usage{
callessfun(file_location,usable_tally_list,parameter_list)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{file_location}{full path of all raw tally files, replicates are separated by \code{","} without space}
  \item{usable_tally_list}{output of \code{self_tallyprepfun}, usable tally files that prepared using raw tally and supporting files}
  \item{parameter_list}{output parameter value from \code{calcparafun} function. Optimized parameters are used to construct null distribution.}
}
\details{
This function constructs potential non-essential distribution using calculated parameters, and calculate p-value for each gene to identify essential genes.
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
callessfun("/home/TnSeq/data/test_data/PA14_M9_rep1_tally.txt",
           usable_tally,parameter_list)
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~calculate_essential }
