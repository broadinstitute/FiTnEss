\name{denote_coreTA}
\alias{denote_coreTA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Identify core TA sites based on TA location
%%  ~~function to do ... ~~
}
\description{
Identify core TA sites, and mark TA sites that locate within 50bp of gene ends. 
}
\usage{
denote_coreTA(tally,bp=50)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tally}{raw tally file, path and name of tally file for run, e.g. `/home/your_folder/your_tally.txt`}
  \item{bp}{define core TA sites by defining distance from gene ends, default is 50bp}
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
denote_coreTA("/home/TnSeq/data/test_data/PA14_M9_rep1_tally.txt",
              50bp)
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~core_TA }
