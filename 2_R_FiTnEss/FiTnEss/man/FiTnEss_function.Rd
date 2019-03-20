\name{FiTnEss_function}
\alias{FiTnEss_function}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Run FiTnEss
%%  ~~function to do ... ~~
}
\description{
Single function to run whole FiTnEss: from tally files to final results
}
\usage{
FiTnEss_Run(strain,file_location,permissive_file,
            homologous_file,gene_file,save_location,
            repeat_time=3)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{strain}{name of the strain}
  \item{file_location}{path and name of tally file for run, e.g. `/home/your_folder/your_tally.txt`}
  \item{permissive_file}{path and name of non-permissive TA site file that generated from genomic pre-processing step}
  \item{homologous_file}{path and name of homologous TA site file that generated from pre-processing step}
  \item{gene_file}{path and name of GFF3 gene file that downloaded from [Pseudomonas Genome Database](http://www.pseudomonas.com/strain/show?id=109)}
  \item{save_location}{path and name of where to save final results file}
  \item{repeat_time}{how many times to run the pipeline in order to obtain best results: by default, we run 3 times}
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

1. install FiTnEss package from github using devtools

devtools::install_github("ruy204/FiTnEss")

#2. load packages

require(FiTnEss)
library(dplyr)
library(tidyr)
library(openxlsx)

#3. run FiTnEss

FiTnEss_Run("UCBPP",
            "/home/TnSeq/data/test_data/PA14_M9_rep1_tally.txt",
            "/home/TnSeq/data/test_data/nonpermissive_TA_sites.txt",
            "/home/TnSeq/data/test_data/homologous_TA_sites.txt",
            "/home/TnSeq/data/test_data/PA14_gene_file.txt",
            "/home/TnSeq/test_result/test_results.xlsx",
            repeat_time = 3)
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~FiTnEss_Run }
