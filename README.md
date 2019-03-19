# _FiTnEss_
Finding Tn-Seq Essential genes (_FiTnEss_)

_FiTnEss_ is a package using Transposon insertion sequencing data to identify essential genes in the genome. 

Preprint manuscript on bioRxiv: [Defining the core essential genome of Pseudomonas aeruginosa](https://www.biorxiv.org/content/early/2019/01/12/396689)

Published manuscript:

#### Quick start

```
fitnessRun(strain,file_location,save_location,repeat_time)

#Example
fitnessRun(strain = "PA14",
           file_location = "/home/documents/raw_tally/PA14_rep1.txt,
                            /home/documents/raw_tally/PA14_rep2.txt,
                            /home/documents/raw_tally/PA14_rep3.txt",
           save_location = "/home/documents/results/PA14_results.xlsx",
           repeat_time = 3)
```















