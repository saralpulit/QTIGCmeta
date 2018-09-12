#!/usr/bin/perl
####################################################################################################
#
# Author:           Paul de Bakker
#                   Division of Genetics, Brigham and Women's Hospital
#                   Program in Medical and Population Genetics, Broad Institute of MIT and Harvard
# 
# Email:            pdebakker@rics.bwh.harvard.edu  or  debakker@broadinstitute.org
#
# Web:              http://debakker.med.harvard.edu
#        
# Acknowledgements: Thanks to Nikolaos Patsopoulos for help with the heterogeneity and 
#                   random-effects code.   Thanks to Chris Cotsapas and Ben Voight for sharing
#                   code for the exponential decay model.
#
# Description:      This Perl script performs a meta-analysis across an arbitrary number of 
#                   genome-wide studies, where for each study association statistics are computed 
#                   for a predefined set of SNPs.  The meta-analysis is based on BETA and SE 
#                   estimates from a linear or logistic regression analysis for each study.  The 
#                   meta-analysis statistic is based on inverse-variance weighting as well as a 
#                   sample-size weighting, where we explicitly correct for imputation quality. We
#                   compute both fixed and random-effects models.  Fisher's rule for combining 
#                   p-values and a novel exponential-decay test as well as conventional tests for 
#                   heterogeneity (Cochran's Q and I-squared) are included.
#
####################################################################################################
#
# Required input:
#
# 1) A plain-text file with association analysis results for all studies combined into a single file
#    (one line per SNP).
#
#    Expected file format (where the number indicates the column):
#
#    SNP CHR POS BETA SE PVALUE CA A2 CAF RATIO
#      1   2   3    4  5      6  7  8   9    10
#
#    where CHR is the chromosome number, POS is the chromosomal position of the SNP,
#    BETA is the computed estimate parameter of the effect of the given SNP,
#    SE is the standard error around that BETA estimate,
#    CA represents the coded allele that the BETA and SE are referring to, A2 is the other allele,
#    CAF refers to the allele frequency of the coded allele (CA), and
#    RATIO is the ratio of the observed variance of the dosage to the expected (binomial) variance.
#
#    Columns 2-10 are repeated (on the same line) for every additional GWAS that is part of the
#    meta-analysis.
#
#    The RATIO is used to correct the weight of the contribution of each individual study depending
#    on the imputation quality of the SNP.  This is only necessary when imputation was used (set to
#    1 if all SNPs are genotyped experimentally).   See de Bakker et al., Human Molecular Genetics,
#    2008 for more background information on this topic.
#
# 2) A plain-text file that contains study-specific parameters
#
#    Example file format:
#
#    FHS    1.023     7650      1
#    CHS    1.040      854     17.2755
#    ERGO   1.034     4606     18.1075
#
#    Column 1 : contains an alphanumeric name to identify the studies listed in the file specified
#               above
#    Column 2 : lists the genomic inflation factor (lambda) -- used for adjusting the SE on the fly
#               (set to 1 if the association results are already adjusted)
#    Column 3 : lists the sample size of each study -- used for meta-analysis based on sample size-
#               weighted z-scores
#    Column 4 : lists the correction factor to standardize the BETA and SE estimates across all
#               studies to ensure the scale and units of the BETA and SE are identical (the BETA
#               and SE are divided by the given factor)
#
#    Note that the order of the studies in this file is important -- it reflects the order of the
#    association results as they appear in the data file specified above.
#
# 3) A PLINK BIM file -- for SNP positions
#
#    Expected file format (note that a PLINK BIM file has NO header line!):
#
#    CHR  SNP     MORGAN     POS  A    B
#      1    2          3       4  5    6
#
# 4) A PLINK generated file (--freq) for HapMap which is used as the reference to resolve
#    ambiguities in allele coding.
#
#    Expected file format:
#
#    CHR          SNP   A1   A2          MAF  NCHROBS
#      1            2    3    4            5        6
#
# 5) Gene annotations
#
#    Expected file format:
#
#    CHR START STOP GENE_SYMBOL
#
# 6) Maximal distance to a gene (in kilobase units) -- this is used in the final output as for
#    every SNP genes within the specified distance are listed.  A reasonable choice is 100-200
#    kilobases.
#
#
####################################################################################################

### VERSION
my $date_version = "08 April 2010";

use lib '.';
use strict;
use FileHandle;
use Getopt::Long;
use Statistics::Distributions;

my $paramsFile;
my $snpFile;
my $extractFile;
my $dbsnpFile;
my $hmfreqFile;
my $genesFile;
my $extension = "";
my $gene_dist = 200;
#my $no_header = '';
my $verbose_p = '';
my $outFile = "meta.out";
my $cutoff1 = 1e-2;
my $cutoff2 = 1e-5;
my $cutoff3 = 1e-8;

GetOptions(
           "params=s"       => \$paramsFile,
	   "snps=s"         => \$snpFile,
	   "extract=s"      => \$extractFile,
	   "dbsnp=s"        => \$dbsnpFile,
	   "freq=s"         => \$hmfreqFile,
	   "genes=s"        => \$genesFile,
	   "dist=i"         => \$gene_dist,
	   "out=s"          => \$outFile,
	   "ext=s"          => \$extension,
	   "verbose"        => \$verbose_p,
           );

if ( ! $paramsFile || ! $snpFile || ! $dbsnpFile || ! $hmfreqFile || ! $genesFile ) {
  die "mantel.pl by Paul de Bakker, last update: $date_version\n\nUsage: mantel.pl --params params_file --snps snps_file --dbsnp dbsnp_file [--freq freq_file] [--genes genes_file] [--verbose] [--dist distance] [--out outfile] [--ext extension]\n";
}

print STDOUT "mantel.pl by Paul de Bakker, last update: $date_version\n";
print STDOUT "running with the following parameters:\n";
print STDOUT "  --params              = $paramsFile\n";
print STDOUT "  --snps                = $snpFile\n";
print STDOUT "  --dbsnp               = $dbsnpFile\n";
print STDOUT "  --freq                = $hmfreqFile\n";
print STDOUT "  --genes               = $genesFile\n";
print STDOUT "  --dist                = $gene_dist\n";
print STDOUT "  --out                 = $outFile\n";

if ( $extension ne "" ) {
  print STDOUT "  --ext                 = $extension\n";
}

if ( $extractFile ne "" ) {
  print STDOUT "  --extract             = $extractFile\n";
}

#if ( $no_header ) {
#  print STDOUT "  --no-header\n";
#}

if ( $verbose_p ) {
  print STDOUT "  --verbose\n";
}

print STDOUT "\n";



################################################################################
################################################################################
###
### read in SNP list
###
################################################################################
################################################################################

my @snp_name = ();
my %snplist = ();
my $n_total_snps = 0;
#my $header = 0;

my $snpFilename = $snpFile;
if ( $extension ne "" ) { $snpFilename .= ".$extension"; }
open (SNP, $snpFilename) or die "cannot open $snpFilename\n";
print STDOUT "reading SNP list from: $snpFilename\n";
while(my $c = <SNP>){
  chomp $c;
  $c =~ s/^\s+//;
  my @fields = split /\s+/, $c;
 
#  if ( $header == 0 && ! $no_header ) { $header = 1; next; }

  $snp_name[$n_total_snps] = $fields[0];
  $snplist{$fields[0]} = 1;

  $n_total_snps++;
}
close(SNP);

print STDOUT "number of SNPs expected in meta-analysis: $n_total_snps\n";



################################################################################
################################################################################
###
### read in SNP extract list ( optional )
###
################################################################################
################################################################################

my %extract = ();
my $n_extract_snps = 0;

if ( $extractFile ) {
  open (SNP, $extractFile) or die "cannot open $extractFile\n";
  print STDOUT "reading SNP extract list from: $extractFile\n";
  while(my $c = <SNP>){
    chomp $c;
    $c =~ s/^\s+//;
    my @fields = split /\s+/, $c;

    $extract{$fields[0]} = 1;
    $n_extract_snps++;
  }
  close(SNP);
  print STDOUT "extracting SNPs: $n_extract_snps\n";
}


################################################################################
################################################################################
###
### read in study parameters
###
################################################################################
################################################################################

my @study_name = ();
my @filename = ();
my @lambda = ();
my @sample_size = ();
my @correction_factor = ();
my @n_strand_flips = ();
my @n_sign_flips = ();
my @n_informative_snps = ();

my @fh = ();
my $n_studies = 0;
my $total_sample_size = 0;

my @col_snp = ();
my @col_beta = ();
my @col_se = ();
my @col_pval = ();
my @col_ca = ();
my @col_nca = ();
my @col_caf = ();
my @col_ratio = ();

my $beta_p = 1; ## assume we will do calcs using beta and SE etc.
my $pval_p = 1; ## assume we will do fishers etc.
my $freq_p = 1; ## assume we will use frequency etc.
my $alleles_p = 1; ## assume we will have coded and non-coded alleles

open (PARAMS, $paramsFile) or die "cannot open $paramsFile\n";

print STDOUT "reading parameter file: $paramsFile\n";
print STDOUT "\n";
print STDOUT "          study name     lambda     sample size     correction factor     number of SNPs     file name\n";
print STDOUT "          ----------     ------     -----------     -----------------     --------------     ---------\n";
 
while(my $c = <PARAMS>){
  chomp $c;
  $c =~ s/^\s+//;
  my @fields = split /\s+/, $c;

  if ( $#fields != 4 ) { die "number of columns in the $paramsFile must be 5!\n"; }

  $study_name[$n_studies] = $fields[0];
  $lambda[$n_studies] = $fields[1];
  if ( $fields[1] < 1 ) { 
    die "lambda $fields[0] cannot be less than 1\n";
  }
  $sample_size[$n_studies] = $fields[2];
  $total_sample_size += $fields[2];
  $correction_factor[$n_studies] = $fields[3];
  $filename[$n_studies] = $fields[4];
  if ( $extension ne "" ) {
    $filename[$n_studies] .= ".$extension";
  }

  # set counters to zero 
  $n_strand_flips[$n_studies] = 0;  
  $n_sign_flips[$n_studies] = 0;  
  $n_informative_snps[$n_studies] = 0;  
  
  $col_snp[$n_studies] = -1;
  $col_beta[$n_studies] = -1;
  $col_se[$n_studies] = -1;
  $col_pval[$n_studies] = -1;
  $col_ca[$n_studies] = -1;
  $col_nca[$n_studies] = -1;
  $col_caf[$n_studies] = -1;
  $col_ratio[$n_studies] = -1;

  # make filehandle for simultaneous access
  $fh[$n_studies] = new FileHandle;
  $fh[$n_studies]->open($filename[$n_studies]) || die "$filename[$n_studies] did not open\n";
  my $FILE = $fh[$n_studies];

  my $counter = 0;
  my $header = 0;
  
  while(<$FILE>) {
    my $line = $_;
    chomp($line); 
    $line =~ s/^\s+//;
    my @cols = split /\s+/, $line;
   
    #if ( $header == 0 && ! $no_header ) { $header = 1; next; }
    #if ( $#cols != 9 ) { die "number of columns in $filename[$n_studies] must be 10\n"; }
    #if ( $cols[0] ne $snp_name[$counter] ) { die "$filename[$n_studies] not in order with $snpFile\n"; }

    if ( $header == 0 ) { 
      for (my $i = 0; $i < $#cols; $i++) { 
        if ( uc($cols[$i]) eq "SNP" ) { $col_snp[$n_studies] = $i; }
        if ( ( $cols[$i] =~ m/beta/i ) || ( $cols[$i] =~ m/effect/i ) ) { $col_beta[$n_studies] = $i; }
        if ( ( $cols[$i] =~ m/stderr/i ) || ( $cols[$i] =~ m/se/i ) ) { $col_se[$n_studies] = $i; }
        if ( uc($cols[$i]) eq "P_FIXED" || ( uc($cols[$i]) =~ m/pval/i ) || uc($cols[$i]) eq "P.VALUE" ) { $col_pval[$n_studies] = $i; }
        if ( uc($cols[$i]) eq "CODED_ALLELE" || uc($cols[$i]) eq "CA" || $cols[$i] eq "A1" ) { $col_ca[$n_studies] = $i; }
        if ( uc($cols[$i]) eq "NONCODED_ALLELE" || uc($cols[$i]) eq "NON_CODED_ALLELE" || uc($cols[$i]) eq "NCA" || $cols[$i] eq "A2" ) { $col_nca[$n_studies] = $i; }
        if ( ( $cols[$i] =~ m/freq/i ) || uc($cols[$i]) eq "CAF" ) { $col_caf[$n_studies] = $i; }
        if ( ( $cols[$i] =~ m/qual/i ) || ( $cols[$i] =~ m/rsq/i ) || ( $cols[$i] =~ m/oevar/i ) || ( $cols[$i] =~ m/ratio/i ) || ( $cols[$i] =~ m/var-quot/i ) ) { $col_ratio[$n_studies] = $i; }
      }
      $header = 1;
    }

    $counter++;
  }
  close($FILE);

  ### reopen it and skip first line
  $fh[$n_studies]->open($filename[$n_studies]) || die "$filename[$n_studies] did not open\n";
#  if ( ! $no_header ) { 
    my $FILE = $fh[$n_studies];
    my $ignore = <$FILE>; 
#  }

  printf STDOUT "%20s %10.5f %15d %21.5f %18d     %s [", $study_name[$n_studies], $lambda[$n_studies], $sample_size[$n_studies], $correction_factor[$n_studies], $counter, $filename[$n_studies];
  if ( $col_beta[$n_studies] != -1 ) { printf STDOUT " BETA"; } else { $beta_p = 0; }
  if ( $col_se[$n_studies]   != -1 ) { printf STDOUT " SE"; } else { $beta_p = 0; }
  if ( $col_pval[$n_studies] != -1 ) { printf STDOUT " PVALUE"; } else { $pval_p = 0; }
  if ( $col_caf[$n_studies]  != -1 ) { printf STDOUT " FREQ"; } else { $freq_p = 0; }
  if ( $col_ca[$n_studies]   != -1 && $col_nca[$n_studies] != -1 ) { printf STDOUT " ALLELES"; } else { $alleles_p = 0; }
  printf STDOUT " ]\n";

  #printf STDOUT " beta=%d se=%d pval=%d caf=%d ca=%d nca=%d\n", $col_beta[$n_studies], $col_se[$n_studies], $col_pval[$n_studies], $col_caf[$n_studies], $col_ca[$n_studies], $col_nca[$n_studies];

  $n_studies++;
}
close(PARAMS);

print STDOUT "                                    ===========\n";
printf STDOUT "                                    %11d\n", $total_sample_size;
print STDOUT "\n";
print STDOUT "total number of studies: $n_studies\n";
  
if ( $alleles_p == 0 && $beta_p == 1 ) { 
  printf STDOUT "Cannot perform beta/se based meta-analysis if coded/noncoded alleles are not given\n";
  $beta_p = 0;
}


################################################################################
################################################################################
###
### read in dbSNP file to get inventory of markers, positions and annotation
###
################################################################################
################################################################################

my $dbsnp_p = 0;
my $n_dbsnp_annotations = 0;
my %skip_list = ();
my %dbsnp_chr = ();
my %dbsnp_pos = ();
my %dbsnp_alleles = ();
my %dbsnp_function = ();
my %dbsnp_caveat = ();

if ( $dbsnpFile ) { 
  open (DBSNP, $dbsnpFile) or die "cannot open $dbsnpFile\n";
  print STDOUT "reading SNP annotations file: $dbsnpFile\n";

  #chrom	chromStart	chromEnd	name	strand	observed	class	func
  #chr1	6946796	6946821	rs57898978	+	A/G	single	intron
  ##chr1	8912885	8912910	rs57188530	+	C/G	single	unknown
  # chr1	34340790	34340885	rs6143185	+	(LARGEDELETION)/-	named	intron
  #chr1	102891517	102891542	rs56752146	+	A/G	single	unknown

  while(my $c = <DBSNP>){
    chomp $c;
    $c =~ s/^\s+//;
    my @fields = split /\s+/, $c;
    my $snp = $fields[3];

    ### skipping SNPs that map to alternate chromosome (e.g. chr6_qbl) 
    if ( $fields[0] =~ m/_/ ) { next; }
  
    if ( defined( $snplist{$snp} ) && ( ( ! $extractFile ) || defined( $extract{$snp} ) ) ) {

      if ( defined( $dbsnp_alleles{$snp} ) ) {
#        print STDERR "$snp appears more than once -- skipping it\n";
        $dbsnp_caveat{$snp} = "MULTIPLE_MAPPINGS";
#        $skip_list{$snp} = 1;
#        next;
      }
    
      my @alleles = split /\//, $fields[5];
   
#      if ( $#alleles > 1 ) { 
#        print STDERR "$snp has more than 2 alleles [" . $fields[5] . "] -- skipping it\n";
#        $skip_list{$snp} = 1;
#        next;
#       } 
      if ( $fields[5] =~ m/lengthTooLong/ ) {
        print STDERR "$snp has alleles with lengthTooLong -- skipping it\n";
        $skip_list{$snp} = 1;
        next;
      }
      if ( $#alleles == 0 ) { 
        print STDERR "$snp has only 1 allele [" . $fields[5] . "] -- skipping it\n";
        $skip_list{$snp} = 1;
        next;
      } 

      $fields[0] =~ s/chr//;
      $dbsnp_chr{$snp} = $fields[0];
      $dbsnp_pos{$snp} = $fields[1] + 1;
      $dbsnp_function{$snp} = $fields[7];
      $dbsnp_alleles{$snp} = [ @alleles ];
#      print "from dbsnp read $snp with $dbsnp_alleles{$snp}[0] $dbsnp_alleles{$snp}[1] alleles strand = $strand\n";
    
      my $strand = $fields[4]; 
      if ( $strand eq "-" ) { 
        @{$dbsnp_alleles{$snp}} = ();
        foreach my $allele ( @alleles ) { 
          push @{$dbsnp_alleles{$snp}}, strand_flip( $allele );
        }
      }
#      print "from dbsnp read $snp with $dbsnp_alleles{$snp}[0] $dbsnp_alleles{$snp}[1] alleles strand = $strand function = $dbsnp_function{$snp}\n";
 
      $n_dbsnp_annotations++;
    }
  }
  close (DBSNP);
  $dbsnp_p = 1;
  print STDOUT "number of annotated SNPs: $n_dbsnp_annotations\n";
} 


################################################################################
################################################################################
###
### check all SNPs on the list if they are in dbSNP 
###
################################################################################
################################################################################

for (my $nsnp; $nsnp < $n_total_snps; $nsnp++) {
  my $snp = $snp_name[$nsnp];
  if ( ! defined( $skip_list{$snp} ) && ( ( ! $extractFile ) || defined( $extract{$snp} ) ) && ! defined( $dbsnp_chr{$snp} ) ) {
    print STDERR "$snp in $snpFile is not annotated in dbSNP -- skipping it\n";
    $skip_list{$snp} = 1;
  }
}



################################################################################
################################################################################
###
### read in the minor and major alleles from the HapMap frequency file 
### (nothing is done currently with those - but may need them later for A/T and C/G SNPs)
###
################################################################################
################################################################################

my %hapmap_a1 = ();
my %hapmap_a2 = ();
my %hapmap_maf = ();
my $n_hapmap_snps = 0;
my $hapmap_p = 0;

if ( $hmfreqFile ) {
  open (HM, $hmfreqFile) or die "cannot open $hmfreqFile\n";
  print STDOUT "reading HapMap frequency file: $hmfreqFile\n";

  while(my $c = <HM>){
    chomp $c;
    $c =~ s/^\s+//;
    my @fields = split /\s+/, $c;
    my $snp = $fields[1]; 
    if ( $snp eq "SNP" ) { next; }

    if ( ( ! defined( $skip_list{$snp} ) ) && defined( $snplist{$snp} ) && ( ( ! $extractFile ) || defined( $extract{$snp} ) ) ) {

      my $a1 = $fields[2];
      my $a2 = $fields[3]; 
      my $alleles_present = 0;

      foreach my $allele ( @{$dbsnp_alleles{$snp}} ) {
        if ( $allele eq $a1 || $allele eq $a2 ) { $alleles_present++; }
      } 

      if ( $alleles_present == 2 || ( $alleles_present == 1 && ( $a1 eq "0" || $a2 eq "0" ) ) ) {
        $hapmap_a1{$snp} = $a1;
        $hapmap_a2{$snp} = $a2;
        $hapmap_maf{$snp} = $fields[4];
        $n_hapmap_snps++;
      }
      else {
        $alleles_present = 0;

        foreach my $allele ( @{$dbsnp_alleles{$snp}} ) {
          if ( $allele eq strand_flip( $a1 ) || $allele eq strand_flip( $a2 ) ) { $alleles_present++; }
        } 

        if ( $alleles_present == 2 || ( $alleles_present == 1 && ( $a1 eq "0" || $a2 eq "0" ) ) ) {
          print STDERR "for $snp, HapMap alleles $a1 and $a2 are not consistent with annotated alleles [" . join("/", @{$dbsnp_alleles{$snp}}) . "] -- but flipping them appears to fix it\n";
          $hapmap_a1{$snp} = strand_flip( $a1 );
          $hapmap_a2{$snp} = strand_flip( $a2 );
          $hapmap_maf{$snp} = $fields[4];
          $n_hapmap_snps++;
        }
        else {
          print STDERR "for $snp, HapMap alleles $a1 and $a2 are not consistent with annotated alleles [" . join("/", @{$dbsnp_alleles{$snp}}) . "] -- even after flipping them -- skipping it\n";
          $skip_list{$snp} = 1;
        }
      }
    }
  }
  close (HM);
  $hapmap_p = 1;
  print STDOUT "number of HapMap SNPs: $n_hapmap_snps\n";
}



################################################################################
################################################################################
###
### read in the genes
### 
################################################################################
################################################################################

my $n_genes = 0;
my @gene_name = ();
my @gene_chr = ();
my @gene_start = ();
my @gene_stop = ();
my $genes_p = 0;

if ( $genesFile ) {
  open (GENE, $genesFile) or die "cannot open $genesFile\n";
  print STDOUT "reading gene annotations file: $genesFile\n";

  while(my $c = <GENE>){
    chomp $c;
    my @fields = split /\s+/, $c;

    $gene_chr[$n_genes] = $fields[0];
    $gene_start[$n_genes] = $fields[1];
    $gene_stop[$n_genes] = $fields[2];
    $gene_name[$n_genes] = $fields[3];

    $n_genes++;
  }
  close (GENE);
  $genes_p = 1;
  print STDOUT "number of annotated genes: $n_genes\n";
  print STDOUT "maximal distance to genes: $gene_dist KB\n";
}




################################################################################
################################################################################
###
### prepare output file -- write out the header line
###
################################################################################
################################################################################

open (OUT, ">$outFile") or die "cannot open $outFile\n";
print OUT "SNP";

if ( $dbsnp_p ) {
  print OUT " CHR POS ALLELES FUNCTION CAVEAT";
}

if ( $hapmap_p ) {
  print OUT " HM_MINOR HM_MAJOR HM_MAF";
}

if ( $verbose_p ) {
  for (my $i=0; $i < $n_studies; $i++) {
    my $name = $study_name[$i];
    print OUT " CA_$name NCA_$name STRAND_FLIPPED_$name SIGN_FLIPPED_$name CAF_$name BETA_$name SE_$name P_$name RATIO_$name NEFF_$name";
  }
}

if ( $alleles_p ) {
  print OUT " CODED_ALLELE NONCODED_ALLELE";
  if ( $freq_p ) {
    print OUT " CODED_ALLELE_FREQ";
  }
}

if ( $pval_p ) {
  print OUT " CHISQ_FISHER DF_FISHER P_FISHER CHISQ_EXPDECAY P_EXPDECAY";
} 

if ( $beta_p ) {
  print OUT " Z_SQRTN P_SQRTN";
  print OUT " BETA_FIXED SE_FIXED Z_FIXED P_FIXED BETA_LOWER_FIXED BETA_UPPER_FIXED";
  print OUT " BETA_RANDOM SE_RANDOM Z_RANDOM P_RANDOM BETA_LOWER_RANDOM BETA_UPPER_RANDOM";
  print OUT " COCHRAN_Q DF_COCHRAN P_COCHRAN I_SQUARED TAU_SQUARED DIRECTIONS";
}

print OUT " K N_EFF";

if ( $genes_p ) {
  print OUT " GENES_" . $gene_dist . "KB NEAREST_GENE";
}

print OUT "\n";


################################################################################
################################################################################
###
### now loop over all SNPs from the list - and do some work
###
################################################################################
################################################################################

print STDOUT "now running meta-analysis...\n";

my $n_snps_in_meta = 0;
my $n_snps_not_on_hapmap = 0;
my $n_uninformative_snps = 0;
my $n_snps_skipped = 0;

for (my $nsnp; $nsnp < $n_total_snps; $nsnp++) {

  my $snp = $snp_name[$nsnp];
  
  ### skip this SNP if it is not on the extract hash 
  if ( $skip_list{$snp} || ( $extractFile && ! defined( $extract{$snp} ) ) ) { $n_snps_skipped++; next; }
  
  my $refchr = defined($dbsnp_chr{$snp}) ? $dbsnp_chr{$snp} : "NA";
  my $refpos = defined($dbsnp_pos{$snp}) ? $dbsnp_pos{$snp} : "NA";
  my $ref1 = "0";
  my $ref2 = "0";
  my $coded_allele = "0";
  my $noncoded_allele = "0";
  
  my @beta = ();
  my @se = ();
  my @pval = ();
  my @ca = ();
  my @nca = ();
  my @caf = ();
  my @ratio = ();
  
  my $skip_freq_p = 0;
  my @study_ok = ();
  my @flip_alleles = ();
  my @effective_size = ();
  my $n_ok_studies = 0;

  my $total_weight = 0;
  my $total_weight_squared = 0;
  my $total_weighted_beta = 0;
  my $total_weighted_beta_squared = 0;
  my $n_eff = 0;
  my $z_sqrtn = 0;
  my $af_weighted = 0;
  my $chisq_fisher = 0;
  my $df_fisher = 0;
  my $log_pvals_sum = 1e-8;

#  my $bin1 = 0; # p > 0.01
#  my $bin2 = 0; # 1e-4 < p < 0.01
#  my $bin3 = 0; # 1e-7 < p
#  my $bin4 = 0; # 1e-7 < p

 
  ###
  ### first read in the data
  ###
  for ( my $study = 0; $study < $n_studies; $study++ ) {
    my $FILE = $fh[$study]; 
    for ( my $a = 0; $a < $n_snps_skipped; $a++ ) { my $void = <$FILE>; }
    my $c = <$FILE>;
    chomp $c;
    my @fields = split /\s+/, $c;

    if ( $fields[0] ne $snp ) { die "$filename[$study] not in order with $snpFile -- was expecting $snp but got $fields[0]\n"; }

    if ( $beta_p ) {
      $beta[$study] = $fields[ $col_beta[$study] ];
      $se[$study]   = $fields[ $col_se[$study] ];
    }

    if ( $pval_p ) {
      $pval[$study] = $fields[ $col_pval[$study] ];
    }

    if ( $alleles_p ) {
      $ca[$study]  = allele_1234_to_ACGT( $fields[ $col_ca[$study] ] );
      $nca[$study] = allele_1234_to_ACGT( $fields[ $col_nca[$study] ] );
    }
    
    if ( $freq_p ) {
      $caf[$study] = $fields[ $col_caf[$study] ];
      if ( $caf[$study] eq "NA" || $caf[$study] == 0 || $caf[$study] == 1 ) { $skip_freq_p = 1; }
    }

    if ( $col_ratio[$study] != -1 ) {
      $ratio[$study] = $fields[ $col_ratio[$study] ];
    } else {
      $ratio[$study] = 1;
    }

    #print " beta = $beta[$study]  se = $se[$study]  pval = $pval[$study]  ca = $ca[$study]  nca = $nca[$study]  ratio = $ratio[$study]\n";
  }

  ### reset skip 
  $n_snps_skipped = 0;
   
  ### check if this SNP has a HapMap frequency 
  if ( ! defined( $hapmap_maf{$snp} ) ) {
    $n_snps_not_on_hapmap++;
  }
  else {
    ### if SNP is not monomorphic -- this is only so that A1 is listed as the HapMap minor allele
    if ( $hapmap_a1{$snp} ne "0" && $hapmap_a2{$snp} ne "0" ) {
      $ref1 = $hapmap_a1{$snp};
      $ref2 = $hapmap_a2{$snp};
    }
  }
  
  ###
  ### walk through studies to see who's present
  ###
  for ( my $study=0; $study < $n_studies; $study++ ) {  
      
    ### first assume the study is okay
    $study_ok[$study] = 1;
    
    if ( ( $beta_p && ( $beta[$study] eq "NA" || $se[$study] eq "NA" || $se[$study] == 0 ) ) || ( $pval_p && ( $pval[$study] eq "NA" || $pval[$study] == 0 ) ) || ( $alleles_p && ( $ca[$study] eq "NA" || $nca[$study] eq "NA" ) ) ) {

      $study_ok[$study] = 0;

    } else {
  
      if ( $alleles_p ) {
        ### are the reference alleles set yet?
        if ( $ref1 eq "0" && $ref2 eq "0" ) { 

          ### first check if they are consistent with dbSNP
          my $alleles_present = 0;
          foreach my $allele ( @{$dbsnp_alleles{$snp}} ) {
            if ( $ca[$study] eq $allele || $nca[$study] eq $allele ) { $alleles_present++; }
          }

          if ( $alleles_present == 2 ) {
            $ref1 = $ca[$study];
            $ref2 = $nca[$study];
            $study_ok[$study] = 1;
            $flip_alleles[$study] = 0;
          }
          else {
            $alleles_present = 0;
            foreach my $allele ( @{$dbsnp_alleles{$snp}} ) {
              if ( strand_flip( $ca[$study] ) eq $allele || strand_flip( $nca[$study] ) eq $allele ) { $alleles_present++; }
            }
            if ( $alleles_present == 2 ) {
              $ref1 = strand_flip( $ca[$study] );
              $ref2 = strand_flip( $nca[$study] );
              $study_ok[$study] = 1;
              $flip_alleles[$study] = 1;
            }
          }
        }
        ### allele a1 and allele a2 match the two reference alleles 
        elsif ( ( $ca[$study] eq $ref1 && $nca[$study] eq $ref2 ) || ( $ca[$study] eq $ref2 && $nca[$study] eq $ref1 ) ) { 
          $flip_alleles[$study] = 0;

          ###
          ### here we could do a frequency-based test for A/T or C/G SNPs (to be done)
          ###
        }
        ### allele a1 and allele a2 do not match the two reference alleles 
        elsif ( $ca[$study] ne $ref1 && $ca[$study] ne $ref2 && $nca[$study] ne $ref1 && $nca[$study] ne $ref2 ) { 
          $flip_alleles[$study] = 1;
        } 
        ### coded and noncoded alleles do not match the two reference alleles -- even after flipping
        else {
          print STDERR "in $study_name[$study], $snp has alleles $ca[$study] $nca[$study] inconsistent with reference alleles $ref1 $ref2 -- skipping this SNP for this study\n"; 
          $study_ok[$study] = 0;
        }
      }
    } 

    if ( $study_ok[$study] == 1 ) {
      $effective_size[$study] = $sample_size[$study] * ( $ratio[$study] > 1 ? 1 : $ratio[$study] ) ;
      $n_eff += $effective_size[$study];
      $n_ok_studies++;
    }

    #print " $study_name[$study]  ok=$study_ok[$study]  flip=$flip_alleles[$study] n_eff=$n_eff  n_ok_studies=$n_ok_studies\n";
  }

  if ( $n_ok_studies == 0 || $n_eff == 0 ) { 
    print STDERR "for $snp , number of studies is $n_ok_studies and effective sample size is $n_eff -- so skipping this SNP\n";
    $n_uninformative_snps++;
  }
  
  if ( $n_ok_studies > 0 && $n_eff > 0 ) {  

    ###
    ### we can now print out SNP information (based on HapMap, if available)
    ###
    print OUT "$snp";
   
    if ( $dbsnp_p ) {
      print OUT " $refchr $refpos " . join("/", @{$dbsnp_alleles{$snp}}) . " $dbsnp_function{$snp}";
      if ( defined( $dbsnp_caveat{$snp} ) ) { 
        print OUT " $dbsnp_caveat{$snp}";
      } else { 
        print OUT " N";
      }
    }

    if ( $hapmap_p ) {
      if ( defined( $hapmap_maf{$snp} ) ) { print OUT " $hapmap_a1{$snp} $hapmap_a2{$snp} $hapmap_maf{$snp}"; } else { print OUT " NA NA NA"; }
    }
    
    my @signed_beta = ();
    my @weight = (); 

    ###
    ### iterate over the association results from all studies
    ###
    for ( my $study = 0; $study < $n_studies; $study++ ) {  

      ### if everything is really okay, proceed
      if ( $study_ok[$study] == 1 ) {

        $n_informative_snps[$study] += 1;

        my $sign = 1;
        my $strand_flipped = "N";
        my $sign_flipped = "N";

        if ( $flip_alleles[$study] == 1 ) {
          $ca[$study] = strand_flip( $ca[$study] );
          $nca[$study] = strand_flip( $nca[$study] );
          $strand_flipped = "Y";
          $n_strand_flips[$study] += 1;
        }

        ### the first study sets the coded/noncoded alleles and is checked against for the other studies
        if ( $coded_allele eq "0" && $noncoded_allele eq "0" ) {
          $coded_allele = $ca[$study];
          $noncoded_allele = $nca[$study];
        } 
        else {
          if ( $ca[$study] eq $coded_allele && $nca[$study] eq $noncoded_allele ) {
            #$sign = 1;
          }
          elsif ( $ca[$study] eq $noncoded_allele && $nca[$study] eq $coded_allele ) {
            # change the sign of the beta if the coded/noncoded alleles are reversed compared to the first study
            #print STDERR "flipped sign for $a1 $a2\n";
            $sign = -1; 
            $sign_flipped = "Y";
            $n_sign_flips[$study] += 1;
          } 
          else {
            die "internal error: coded and non-coded alleles do not match reference alleles\n";
          }
        }

        my $z = 0;
 
        if ( $beta_p ) {
        
          ### put BETA and SE on same scale across studies and correct SE for inflation 
          $beta[$study] = $beta[$study] * $correction_factor[$study];
          $se[$study] = $se[$study] * sqrt($lambda[$study]) * $correction_factor[$study];

          ### inverse variance weighting
          $signed_beta[$study] = $sign * $beta[$study];
          $weight[$study] = 1 / ( $se[$study] * $se[$study] );
          my $weighted_beta = $signed_beta[$study] * $weight[$study];
          $total_weighted_beta += $weighted_beta;
          $total_weight += $weight[$study];
          $total_weight_squared += $weight[$study] * $weight[$study];
          
          ### sample-size weighted z-score
          my $z_weight = sqrt( $effective_size[$study] / $n_eff ); 
          $z = $signed_beta[$study] / $se[$study];
          $z_sqrtn += $z * $z_weight;
        }

        ### sample-size weighted allele frequency
        if ( $freq_p && ! $skip_freq_p ) {
          my $af_weight = $effective_size[$study] / $n_eff; 
          if ( $sign == -1 ) { $caf[$study] = 1 - $caf[$study]; }
          $af_weighted += $caf[$study] * $af_weight; 
        }

        ### Fisher's
        if ( $pval_p ) {
          ## first apply lambda - use z-score if available
          my $chisq = 1;
          if ( $beta_p ) {
            $chisq = $z * $z;
          }
          else {
            $chisq = Statistics::Distributions::chisqrdist( 1, $pval[$study] ) / $lambda[$study];
          }
          $pval[$study] = Statistics::Distributions::chisqrprob( 1, $chisq );
          if ( $pval[$study] == 0 ) { $pval[$study] = 1e-100; }

          $chisq_fisher += -2 * log( $pval[$study] );
          $df_fisher += 2;

          ### log test
          $log_pvals_sum += -log( $pval[$study] );

          ### binned pvalues
#          if ( $pval[$study] >= $cutoff1 ) { $bin1++; } 
#          elsif ( $pval[$study] < $cutoff1 && $pval[$study] >= $cutoff2 ) { $bin2++; }
#          elsif ( $pval[$study] < $cutoff2 && $pval[$study] >= $cutoff3 ) { $bin3++; }
#          elsif ( $pval[$study] < $cutoff3 ) { $bin4++; }
        }

        ### dump out per-study statistics (only if verbose_p)
        if ( $verbose_p ) {
          print OUT sprintf(" %s %s %s %s %.4f %.4f %.4f %.2e %.4f %.1f", $ca[$study], $nca[$study], $strand_flipped, $sign_flipped, $caf[$study], $signed_beta[$study], $se[$study], $pval[$study], $ratio[$study], $effective_size[$study]);
        } 
      }
      else {
        if ( $verbose_p ) {
          print OUT " NA NA NA NA NA NA NA NA NA NA";
        }
      }
    }

    ###
    ### print out the overall summary and annotated genes
    ###
    
    if ( $alleles_p ) {
      printf OUT " %s %s", $coded_allele, $noncoded_allele;
      if ( $freq_p ) { 
        if ( $skip_freq_p ) { print OUT " NA"; } else { printf OUT " %.3f", $af_weighted; }
      }
    }

    if ( $pval_p ) {
      ### Fisher's
      my $p_fisher = Statistics::Distributions::chisqrprob( $df_fisher, $chisq_fisher ); 
      if ( $p_fisher == 0 ) { 
        printf OUT " %.3f %d NA", $chisq_fisher, $df_fisher;
      }
      else {
        printf OUT " %.3f %d %.4e", $chisq_fisher, $df_fisher, $p_fisher;
      }
     
      ### binned pvalues
#      my $exp_bin4 = $n_ok_studies * $cutoff3;
#      my $exp_bin3 = ( $n_ok_studies * $cutoff2 ) - $exp_bin4;
#      my $exp_bin2 = ( $n_ok_studies * $cutoff1 ) - $exp_bin4 - $exp_bin3;
#      my $exp_bin1 = $n_ok_studies - $exp_bin4 - $exp_bin3 - $exp_bin2;

      ### exp-decay test
      ## return the log likelihood ratio of -log(p) being exponentially
      ## distributed (cf biased exponential). Effectively, test if
      ## exponential decay parameter == 1

      ## this should be chi-sq,df=1
      ## abs(-2 * log(    prod(int.exp.fn(pvals,1/mean(-log(pvals)))) / prod(int.exp.fn(pvals,1)) ))
      my $log_pvals_mean = $log_pvals_sum / $n_ok_studies;
      my $lambda = 1.0 / $log_pvals_mean;
      my $exp_alt = 1;
      my $exp_null = 1;

      for (my $study = 0; $study < $n_studies; $study++) {
        if ( $study_ok[$study] == 1 ) {
          $exp_alt *= int_exp_fn($pval[$study], $lambda);
          $exp_null *= int_exp_fn($pval[$study], 1.0);
        }
      }
      my $chisq_expdecay = abs(-2 * log( $exp_alt / $exp_null ));
      my $p_expdecay = Statistics::Distributions::chisqrprob(1, $chisq_expdecay);
      if ( $p_expdecay == 0 ) {
        printf OUT " %.3f NA", $chisq_expdecay;
      }
      else {
        printf OUT " %.3f %.4e", $chisq_expdecay, $p_expdecay;
      }
    }
    
    if ( $beta_p ) {
      my $p_sqrtn = Statistics::Distributions::chisqrprob( 1, $z_sqrtn * $z_sqrtn );
      if ( $p_sqrtn == 0 ) {
        printf OUT " %.4f NA", $z_sqrtn;
      } 
      else {
        printf OUT " %.4f %.4e", $z_sqrtn, $p_sqrtn;
      }
 
      ### fixed effects
      my $weighted_mean_beta = $total_weighted_beta / $total_weight;
      my $se_mean_beta = sqrt( 1 / $total_weight );
      my $z = $weighted_mean_beta / $se_mean_beta;
      my $beta_lower = $weighted_mean_beta - 1.96 * ( $se_mean_beta );
      my $beta_upper = $weighted_mean_beta + 1.96 * ( $se_mean_beta );
      my $p = Statistics::Distributions::chisqrprob( 1, $z * $z );
  
      if ( $p == 0 ) { 
        printf OUT " %.3f %.3f %.3f NA %.3f %.3f ", $weighted_mean_beta, $se_mean_beta, $z, $p, $beta_lower, $beta_upper;
      }
      else {
        printf OUT " %.3f %.3f %.3f %.4e %.3f %.3f ", $weighted_mean_beta, $se_mean_beta, $z, $p, $beta_lower, $beta_upper;
      }

      my $df = $n_ok_studies - 1;

      my $cochran_q = 0;
      my $p_cochran = 1;
      my $i_squared = 0;
      my $tau_squared = 0;

      ### by default, set the random-effects results equal to the fixed-effects results
      my $weighted_mean_beta_random = $weighted_mean_beta;
      my $se_mean_beta_random = $se_mean_beta;
      my $z_random = $z;
      my $beta_lower_random = $beta_lower;
      my $beta_upper_random = $beta_upper;
      my $p_random = $p;

      ### 
      ### only test heterogeneity when >2 studies
      ###
      if ( $n_ok_studies > 2 ) {

        ### Cochran's Q
        for (my $study = 0; $study < $n_studies; $study++) {
          if ( $study_ok[$study] == 1 ) {
            my $tmp = $signed_beta[$study] - $weighted_mean_beta;
            $cochran_q += $weight[$study] * $tmp * $tmp;
          }
        }
        $p_cochran = Statistics::Distributions::chisqrprob($df, $cochran_q);

        ### I-squared
        $i_squared = 100.0 * ( $cochran_q - $df ) / $cochran_q;
        if ( $i_squared < 0 ) { $i_squared = 0; }

        ### random effects
        my $total_weighted_beta_random = 0;
        my $total_weight_random = 0;

        my $mean_weight = $total_weight / $n_ok_studies;
        my $variance_weights = ( 1 / $df ) * ( $total_weight_squared - ( $n_ok_studies * ( $mean_weight * $mean_weight ) ) ); # this is the variance of the FE weights
        my $U = $df * ( $mean_weight - ( $variance_weights / ( $n_ok_studies * $mean_weight ) ) ); # U is used in the calculation of tau2

        ### tau-squared 
        $tau_squared = ( $cochran_q - $df ) / $U;
        if ( $cochran_q <= $df ) { $tau_squared = 0; }

        for (my $study = 0; $study < $n_studies; $study++) {
          if ( $study_ok[$study] == 1 ) {
            my $weight_random = 1 / ( ( 1 / $weight[$study] ) + $tau_squared ); 
            $total_weighted_beta_random += $signed_beta[$study] * $weight_random;
            $total_weight_random += $weight_random;
          }
        }
 
        $weighted_mean_beta_random = $total_weighted_beta_random / $total_weight_random;
        $se_mean_beta_random = sqrt( 1 / $total_weight_random );
        $z_random = $weighted_mean_beta_random / $se_mean_beta_random;
        $beta_lower_random = $weighted_mean_beta_random - 1.96 * ( $se_mean_beta_random );
        $beta_upper_random = $weighted_mean_beta_random + 1.96 * ( $se_mean_beta_random );
        $p_random = Statistics::Distributions::chisqrprob(1, $z_random * $z_random);  
      }

      if ( $p_random == 0 ) {
        printf OUT "%.3f %.3f %.3f NA %.3f %.3f ", $weighted_mean_beta_random, $se_mean_beta_random, $z_random, $beta_lower_random, $beta_upper_random;
      }
      else {
        printf OUT "%.3f %.3f %.3f %.4e %.3f %.3f ", $weighted_mean_beta_random, $se_mean_beta_random, $z_random, $p_random, $beta_lower_random, $beta_upper_random;
      }
   
      if ( $p_cochran == 0 ) {
        printf OUT "%.3f %d NA %.1f %.3f ", $cochran_q, $df, $i_squared, $tau_squared;
      }
      else {
        printf OUT "%.3f %d %.4e %.1f %.3f ", $cochran_q, $df, $p_cochran, $i_squared, $tau_squared;
      }
   
      ### print out the directions for each study
      for (my $study = 0; $study < $n_studies; $study++) {
        if ( $study_ok[$study] == 1 ) { 
          if ( $signed_beta[$study] > 0 ) { printf OUT "+"; } else { printf OUT "-"; }  
        } 
        else { 
          printf OUT ".";
        }
      }
    }

    printf OUT " %d %.1f", $n_ok_studies, $n_eff;
  
    if ( $genes_p == 1 ) {
      ### print out any nearby genes
      my $yes_genes = 0;
      my %listed_genes = ();
      my $nearest_gene = "NA";
      my $nearest_distance = 10000000000;
      my $gene_length_temp = 10000000000;
      my $left_most = $refpos - ($gene_dist * 1000);
      my $right_most = $refpos + ($gene_dist * 1000);

      for (my $i = 0; $i < $n_genes; $i++) {
        if ( $gene_chr[$i] eq $refchr && ! defined( $listed_genes{ $gene_name[$i] } ) ) {
          if ( ( $gene_start[$i] > $left_most && $gene_start[$i] < $right_most ) || ( $gene_stop[$i] > $left_most && $gene_stop[$i] < $right_most ) ) {
            if ( $yes_genes == 1 ) { print OUT ","; } else { print OUT " "; }
            print OUT "$gene_name[$i]"; 
            $yes_genes = 1;
            $listed_genes{$gene_name[$i]} = 1;
        
            my $gene_length = $gene_stop[$i] - $gene_start[$i];
            my $dist_left = $refpos - $gene_stop[$i]; 
            my $dist_right = $gene_start[$i] - $refpos;

            if ( $refpos > $gene_start[$i] && $refpos < $gene_stop[$i] && $gene_length < $gene_length_temp ) { 
              $nearest_gene = $gene_name[$i];
              $gene_length_temp = $gene_length;
              $nearest_distance = 0;
            }
            elsif ( $dist_left > 0 && $dist_left < $nearest_distance ) {
              $nearest_gene = $gene_name[$i];
              $nearest_distance = $dist_left;
            }
            elsif ( $dist_right > 0 && $dist_right < $nearest_distance ) {
              $nearest_gene = $gene_name[$i];
              $nearest_distance = $dist_right;
            }
          }
        } 
      }

      if ( $yes_genes == 0 ) { print OUT " NA"; }
      print OUT " $nearest_gene";
    }
    
    print OUT "\n";

    $n_snps_in_meta++;
  }
}
close(OUT);

 
################################################################################
################################################################################
###
### print summary and close everything 
###
################################################################################
################################################################################

print STDOUT "number of SNPs in meta-analysis: $n_snps_in_meta\n";
print STDOUT "number of SNPs not on HapMap: $n_snps_not_on_hapmap\n";
print STDOUT "number of SNPs uninformative and skipped: $n_uninformative_snps\n";
print STDOUT "\n";
print STDOUT "          study name     allele flips     sign flips     informative SNPs\n";
print STDOUT "          ----------     ------------     ----------     ----------------\n";

for (my $study = 0; $study < $n_studies; $study++) {
  close $fh[$study]; 
  printf STDOUT "%20s %16d %14d %20d\n", $study_name[$study], $n_strand_flips[$study], $n_sign_flips[$study], $n_informative_snps[$study];
}

print STDOUT "\n";
print STDOUT "mantel successfully finished!!!\n";






sub strand_flip($)
{
  my $allele = shift;
  my $flipped_allele = "";
  if ( $allele eq "(LARGEDELETION)" || $allele eq "lengthTooLong" ) { return $allele; }

  for (my $i=0; $i < length($allele); $i++) {
    my $current_base = substr $allele, $i, 1;
    if ( $current_base eq "A" ) { $flipped_allele .= "T"; }
    elsif ( $current_base eq "C" ) { $flipped_allele .= "G"; }
    elsif ( $current_base eq "G" ) { $flipped_allele .= "C"; }
    elsif ( $current_base eq "T" ) { $flipped_allele .= "A"; }
    else { $flipped_allele .= $current_base; }    
  }

  return $flipped_allele;
}


sub allele_1234_to_ACGT($)
{
  my $allele = shift;
  if ( $allele eq "1" ) { return "A"; }
  if ( $allele eq "2" ) { return "C"; }
  if ( $allele eq "3" ) { return "G"; }
  if ( $allele eq "4" ) { return "T"; }
  return $allele;
}


sub int_exp_fn($)
{
  my $x = $_[0];
  my $lambda = $_[1];
  my $epsilon = 0.001;

  return(exp(-$lambda * (-log($x)-$epsilon)) - exp(-$lambda * (-log($x)+$epsilon)))
}

