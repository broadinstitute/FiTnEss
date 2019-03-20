# _FiTnEss_
Finding Tn-Seq Essential genes (_FiTnEss_)

_FiTnEss_ is a package to analyze transposon insertion sequencing (Tn-Seq) data to identify essential genes in the genome. 

Manuscript (preprint on bioRxiv): [Defining the core essential genome of Pseudomonas aeruginosa](https://www.biorxiv.org/content/early/2019/01/12/396689)


### Quick start

FiTnEss requires input files that can be generated using your genome and sequencing data and the provided UNIX-based scripts. The instructions are described in the 1_UNIX_genome-preprocessing_mapping folder.

After installing the FiTnEss package as described below, run the main FiTnEss function by ```FiTnEss_Run```

#### Function arguments

Arguments in this function include: 
- **_strain_**
- **_file_location_**: path and name of tally file for run: 
e.g. `"/your/folder/location/your_data_tally.txt"`
- **_permissive_file_**: path and name of non-permissive TA site file that generated from genomic pre-processing step: 
e.g. `"/your/folder/location/non_permissive_TA_sites.txt"`
- **_homologous_file_**: path and name of homologous TA site file that generated from pre-processing step: 
e.g. `"/your/folder/location/homologous_TA_sites.txt"`
- **_gene_file_**: path and name of GFF3 gene annotation file.
e.g. `"/your/folder/location/your_gff3_file.txt"`
- **_save_location_**: path and name of where to save final results file: 
e.g. `"/your/folder/location/results.xlsx"`
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
FiTnEss_Run("PA14",
            "/your/folder/location/Test_set_P_aeruginosa/sample_data/PA14_M9_rep1_tally.txt",
            "/your/folder/location/Test_set_P_aeruginosa/TAsite_info/nonpermissive_TA_sites.txt",
            "/your/folder/location/Test_set_P_aeruginosa/TAsite_info/homologous_TA_sites.txt",
            "/your/folder/location/Test_set_P_aeruginosa/genome_info/PA14.gff.txt",
            "/your/folder/location/Test_set_P_aeruginosa/sample_data/test_results.xlsx",
            repeat_time = 3)
```

