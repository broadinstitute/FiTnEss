\name{calc_TApos}
\alias{calc_TApos}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Map TA sites onto corresponding genes
}
\description{
This function is part of the \code{tallyprepfun} function, which pre-process raw tally file for analysis.

In this step, the function uses raw tally .txt file (imported from function \code{import_raw_tally})
and genelist file (internal data called from function \code{tallyprepfun}) to map TA sites into corresponding genes.
}
\usage{
calc_TApos(tally,genelist)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tally}{raw tally.txt file that imported from \code{import_raw_tally} function}
  \item{genelist}{internal data in the package, strain-specific, containing gene name, TA site positions of each gene}
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
calc_TApos(unique_map_tally_file,genelist)
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ TA position }
\keyword{ TA site mapping }% __ONLY ONE__ keyword per line
