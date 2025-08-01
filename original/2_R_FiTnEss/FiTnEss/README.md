# _FiTnEss_
Finding Tn-Seq Essential genes (_FiTnEss_)

_FiTnEss_ is a package using Transposon insertion sequencing data to identify essential genes in the genome. 

Original paper on bioRxiv: [Defining the core essential genome of Pseudomonas aeruginosa](https://www.biorxiv.org/content/early/2019/01/12/396689)


### Quick start

After installing FiTnEss package, run main FiTnEss function by ```FiTnEss_Run```

#### Function arguments

Arguments in this function include: 
- **_strain_**
- **_file_location_**: path and name of tally file for run: 
e.g. `"/home/your_folder/your_tally.txt"`
- **_permissive_file_**: path and name of non-permissive TA site file that generated from genomic pre-processing step: 
e.g. `"/home/your_folder/non_permissive_TA_sites.txt"`
- **_homologous_file_**: path and name of homologous TA site file that generated from pre-processing step: 
e.g. `"/home/your_folder/homologous_TA_sites.txt"`
- **_gene_file_**: path and name of GFF3 gene annotation file. For example, GFF3 file could be downloaded from [Pseudomonas Genome Database](http://www.pseudomonas.com/strain/show?id=109): 
e.g. `"/your/folder/location/your_gff3_file.txt"`
- **_save_location_**: path and name of where to save final results file: 
e.g. `"/home/results_folder/results.xlsx"`
- **_repeat_time_**: how many times to run the pipeline in order to obtain best results: by default, we run 3 times.

#### Step 1. install devtools package

```
install.packages("devtools")
```

#### Step 2. install FiTnEss package from github

```
devtools::install_github("ruy204/FiTnEss")
```

#### Step 3. load FiTnEss and dependent packages

```
Packages <- c("dplyr","fBasics","goftest","openxlsx","scales","stats","tidyr")
lapply(Packages, library, character.only = TRUE)

require(FiTnEss)
```
#### Step 4. run FiTnEss

```
FiTnEss_Run("PA14",
            "/your/folder/location/Test_set_P_aeruginosa/sample_data/PA14_M9_rep1_tally.txt",
            "/your/folder/location/Test_set_P_aeruginosa/TAsite_info/nonpermissive_TA_sites.txt",
            "/your/folder/location/Test_set_P_aeruginosa/TAsite_info/homologous_TA_sites.txt",
            "/your/folder/location/Test_set_P_aeruginosa/genome_info/PA14_gff.txt",
            "/your/folder/location/Test_set_P_aeruginosa/sample_data/test_results.xlsx",
            repeat_time = 3)
```

#### Step 5. retrieve results

|Locus.CIA |gtot|Nta|pvalue  |padj|Ess_fwer|pfdr    |Ess_fdr|
|----------|----|---|--------|----|--------|--------|-------|
|PA14_00410|   5|  1|0.015989|   1| NE_fwer|0.093033| NE_fdr|

Each tab in the .xlsx file saves results from each replicate.
Within each results table, there are 8 columns: 
- Locus.CIA: gene index
- gtot: total reads for the gene
- Nta: number of TA sites in this gene
- pvalue: unadjusted p-value of being essential
- padj: FWER-adjusted p-value
- Ess_fwer: confident essential category
- pfdr: FDR-adjusted p-value
- Ess_fdr: candidate essential category




















