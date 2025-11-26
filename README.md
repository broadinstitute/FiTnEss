# _FiTnEss_
Finding Tn-Seq Essential genes (_FiTnEss_)

_FiTnEss_ is a package to analyze transposon insertion sequencing (Tn-Seq) data to identify essential genes in the genome. 

Manuscript (preprint on bioRxiv): [Defining the core essential genome of Pseudomonas aeruginosa](https://www.biorxiv.org/content/early/2019/01/12/396689)

Modified on<br>
11/26/2025:<br>
* allow user to input names of contigs/plasmids that should be removed from analysis (after finding regions of homology and before fitting the distribution of reads). This may be useful when some plasmids are in high copy number, which would result in ambiguity concerning essentiality for those genes and potentially skew the count distribution.<br>
8/7/2025:<br>
* allow for genomes with multiple contigs
* allow user to input gff tag to label genes
* updated to python3

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
- **_gff_name_tag_**: tag in GFF file that uniquely identifies each gene, usually "locus_tag": 
e.g. `"locus_tag"`
- **_remove_multicopy_plasmid_names_**: a vector of strings naming contigs (matching the gff file) that should be removed from analysis, leave as NA (default) if there are none to be removed
- **_repeat_time_**: how many times to run the pipeline in order to obtain best results: by default, we run 3 times.

#### Step 1. install and load dependent packages

```
packages <- c("devtools", "dplyr", "fBasics", "goftest", "openxlsx", "scales", "stats", "tidyr")

for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}
```

#### Step 2. install and load FiTnEss package from github

```
if (!require("FiTnEss", character.only = TRUE)) {
  devtools::install_github("broadinstitute/FiTnEss", subdir = "2_R_FiTnEss/FiTnEss")
}

library(FiTnEss)

```
#### Step 3. run FiTnEss

```
FiTnEss_Run(strain = "PA14",
            file_location = "/your/folder/location/Test_set_P_aeruginosa/sample_data/PA14_M9_rep1_tally.txt",
            permissive_file = "/your/folder/location/Test_set_P_aeruginosa/TAsite_info/nonpermissive_TA_sites.txt",
            homologous_file = "/your/folder/location/Test_set_P_aeruginosa/TAsite_info/homologous_TA_sites.txt",
            gene_file = "/your/folder/location/Test_set_P_aeruginosa/genome_info/PA14.gff",
            save_location = "/your/folder/location/Test_set_P_aeruginosa/sample_data/test_results.xlsx",
            gff_name_tag = "locus_tag", #tag in gff file that uniquely labels genes
            remove_multicopy_plasmid_names = NA,
            repeat_time = 3)
```


