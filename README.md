# _FiTnEss_
Finding Tn-Seq Essential genes (_FiTnEss_)

_FiTnEss_ is a package to analyze transposon insertion sequencing (Tn-Seq) data to identify essential genes in the genome. 

Manuscript (preprint on bioRxiv): [Defining the core essential genome of Pseudomonas aeruginosa](https://www.biorxiv.org/content/early/2019/01/12/396689)


### Quick start

After installing the FiTnEss package as described below, run the main FiTnEss function by ```FiTnEss_Run```

#### Function arguments

Arguments in this function include: 
- **_strain_**
- **_file_location_**: path and name of tally file for run: 
e.g. `"/home/your_folder/your_tally.txt"`
- **_permissive_file_**: path and name of non-permissive TA site file that generated from genomic pre-processing step: 
e.g. `"/home/your_folder/non_permissive_TA_sites.txt"`
- **_homologous_file_**: path and name of homologous TA site file that generated from pre-processing step: 
e.g. `"/home/your_folder/homologous_TA_sites.txt"`
- **_gene_file_**: path and name of GFF3 gene file that downloaded from [Pseudomonas Genome Database](http://www.pseudomonas.com/strain/show?id=109): 
e.g. `"/home/your_folder/PA14_gene_file.txt"`
- **_save_location_**: path and name of where to save final results file: 
e.g. `"/home/results_folder/results.xlsx"`
- **_repeat_time_**: how many times to run the pipeline in order to obtain best results: by default, we run 3 times.

#### Step 1. install dependent packages

```
install.packages(c("devtools","dplyr","fBasics","goftest","openxlsx","scales","stats","tidyr"))
```

#### Step 2. install FiTnEss package from github

```
devtools::install_github("ruy204/FiTnEss")
```

#### Step 3. load FiTnEss and dependent packages

```
Packages <- c("devtools","dplyr","fBasics","goftest","openxlsx","scales","stats","tidyr")
lapply(Packages, library, character.only = TRUE)

require(FiTnEss)
```
#### Step 4. run FiTnEss

```
FiTnEss_Run("UCBPP",
            "/home/TnSeq/data/test_data/PA14_M9_rep1_tally.txt",
            "/home/TnSeq/data/test_data/nonpermissive_TA_sites.txt",
            "/home/TnSeq/data/test_data/homologous_TA_sites.txt",
            "/home/TnSeq/data/test_data/PA14_gene_file.txt",
            "/home/TnSeq/test_result/test_results.xlsx",
            repeat_time = 3)
```

