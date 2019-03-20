\name{calcparafun}
\alias{calcparafun}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Calculating global parameters
}
\description{
This function is part of \code{fitnessRun} function, where global parameters describing non-essential genes are being calculated.

In this step, usable tally files prepared in \code{tallyprepfun} were imported, and global parameters were estimated with minimized Cramer-von Mises value.
}
\usage{
calcparafun(strain,usable_tally_list,save_location,rep_time)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{strain}{strain name}
  \item{usable_tally_list}{output of \code{self_tallyprepfun}, usable tally files that prepared using raw tally and supporting files}
  \item{save_location}{path and name of where to save final results file, e.g. \code{"/home/results_folder/results.xlsx"}}
  \item{rep_time}{same argument of \code{repeat_time} from main function. Define how many times to run the pipeline in order to obtain best results: by default, we run 3 times.}
  
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
calcparafun("PA14",usable_tally,
            "/home/TnSeq/results/results.xlsx",rep_time=5)
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
