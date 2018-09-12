# QTIGCmeta
### Data and code relevant to the QT-IGC meta-analysis of QT interval

## Data for running meta-analysis

(1) meta_qtigc.rev.17Jan2012.params    
Params file to be passed to mantel.pl. Lists the cohort name, lambda from GWAS of that individual cohort, sample size, phenotype adjustment factor (if necessary) and name of the file containing summary-level data for that cohort for the meta-analysis.

(2) Statistics/Distributions.pm   
Perl module for calculating probabilities from distributions (used by mantel.pl)

## Scripts related to running the QT-IGC meta-analysis
(1) mantel.pl and meta_qt_03092011.pl run inverse variance-weighted fixed effects meta-analysis

(2) qqplot_by_INFO.R and qqplot_by_maf.R plot a QQ plot, stratifying SNPs by imputation info score and minor allele frequency, respectively.

(3) parse_loci.pl and parse_clumps.pl for processing post-meta data.
