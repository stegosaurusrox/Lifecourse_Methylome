# Lifecourse_Methylome

The aim of this project is to characterize the Human methylome across the lifecourse.
We aggregated data from several UK-based cohorts with DNA methylation data (Illumina 450k or EPIC array):

+ ARIES (mums, dads, and children)
+ Lothian Birth Cohort (1921 and 1936)
+ Understanding Society
+ 1958 cohort
+ 1970 cohort
+ TwinsUK
+ SABRE (two time points at age 50y and age 70y)

Analysis flow of the project was as follows:

a) harmonize methylation data across all cohorts
  + (re-)process raw IDAT files using the meffil R package (done for all cohorts, apart from Understanding Society, for which idats were not available)
  + merge qc.object files from all cohorts to harmonize data
  + functional normalization of methylation data using cross-harmonized data
 
 b) generate summary statistics per CpG and for demographics (sex, age, cell type)
 
 For cohorts spanning a large age range (>10years), create 10-year age bins and run a) and b) for each bin.
  

## 1) Harmonization

To achieve harmonization without sharing individual level data, we used the meffil R package.
We don't have the idats for Understanding Society, so this cohort is not part of this harmonization protocol.
To harmonize DNA methylation data for the remaining cohorts, we need to go throught three steps using the meffil R package:

- step 1: generate `qc.object` files for each cohort (done at each site and shared centrally)
- step 2: merge `qc.object` files, quantile normalize together, send `norm.object` files back (done centrally) 
- step 3: use `qc.object` files to normalize samples (done at each site)

Importantly, step 1 need to be run twice for common CpG sites (shared across 450k and EPIC) and EPIC-only CpGs. In ARIES, we need to run this three times (450k-only, EPIC-only and common probes), because we have both EPIC and 450k arrays.
Step 2 also needs to be run three times, for 450k, EPIC and common probes.

For a detailed protocol, see here.

## 2) extract CpG-specific summary statistics
