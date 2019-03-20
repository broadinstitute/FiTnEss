library(dplyr)
require(graphics)
require(pastecs)
require(ggplot2)

#4. denotate core TAs

denote_coreTA <- function(tally, bp = 50) {
  tally$coreTA <- ((tally$TA.gene.pos > bp) &
                     (tally$TA.gene.pos < (tally$gene.size.CIA-bp + 1)))
  return(tally)
}
