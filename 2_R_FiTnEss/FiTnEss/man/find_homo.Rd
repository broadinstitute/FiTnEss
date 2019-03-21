\name{find_homo}
\alias{find_homo}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Find homologous TA sites from raw tally file
%%  ~~function to do ... ~~
}
\description{
This function helps identify homologous TA sites from raw tally file.
}
\usage{
find_homo(tally,homofile)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tally}{raw tally file, path and name of tally file for run, e.g. `/home/your_folder/your_tally.txt`}
  \item{homofile}{path and name of homologous marker file generated from genomic pre-processing step}
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
find_homo("/home/TnSeq/data/test_data/PA14_M9_rep1_tally.txt",
          "/home/TnSeq/data/test_data/homologous_TA_sites.txt")
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~homologous }
