\name{self_tallyprepfun}
\alias{self_tallyprepfun}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Preparing usable tally files from raw data
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
tallyprepfun(strain, file_location)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{strain}{strain name}
  \item{file_location}{full path of all raw tally files, replicates are separated by \code{","} without space}
}
\details{
This is a function for raw tally data processing. Input data are raw tally files and strain type, output data are usable tally files for downstream analysis.

Current available \code{strain} options include \code{c("PA14","X13273","PS75","CF77","BWH015","BWH013","BWH005","BL23","19660")}.

Raw data are usually \code{.txt} files, which contain two columns: \code{}. Example of raw data see \code{data("PA14_raw_tally")}
}
\value{
A list containing all usable tally files of each replicate will be returned.

For example, \code{list(usable_replicate1,usable_replicate2)}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{

}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{

\code{\link{import_raw_tally}}

\code{\link{calc_TApos}}

\code{\link{find_homo}}

\code{\link{denote_coreTA}}

}
\examples{
data("PA14_raw_tally")
usable_list<-tallyprepfun("PA14","PA14_raw_tally_rep1,PA14_raw_tally_rep2")

own_data<-tallyprepfun("CF77","/folder1/CF77_rep1.txt,/folder1/CF77_rep2.txt,/folder1/CF77_rep3.txt")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~tally }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
