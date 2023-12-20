
#=============================================================================#
# This script create .ecf files for the software EasyQC. (BINARY TRAITS ONLY) 
#=============================================================================#

#=====SETUP====================================================================

# install.packages(pacman)
library(pacman)
p_load(glue,
       rio)

#=EDIT THESE STRINGS===========================================================
ANCEST      <- "EUR" #ALL CAPS
PHENO       <- "ADHD" #ALL CAPS
VERSION     <- "" #Specific version (if multiple versions exist)
MAF         <- "005" #MAF filter level ("001" for AFR; "005" for EUR)

#=DO NOT EDIT THESE============================================================

#===INPUT - Sumstats
#conditional statement for file nickname
if (VERSION != "") {
  FILE_NAME <- tolower(glue("{ANCEST}_{PHENO}_{VERSION}"))
} else {
  FILE_NAME <- tolower(glue("{ANCEST}_{PHENO}"))
}

FILE_IN     <- glue("../temp/{ANCEST}/{PHENO}/{FILE_NAME}_PRE_EASYQC.txt") #Pathname of file

COLS_IN    <- paste(c("rsID","coded_all","noncoded_all","Chr","position",
                      "N","Beta","SE","Pval","AF_coded_all", 
                      "imputed","oevar_imp","HWE_pval", "OR", "Neff"), collapse = ";") #Column names
CLASSES_IN <- paste(c(rep("character", 5), rep("numeric", 10)), collapse=";") #Column types (CHARACTER;INTEGER;NUMERIC)
COLS_NEW   <- paste(c("SNP","EFFECT_ALLELE","OTHER_ALLELE","CHR","BP",
                      "N","BETA","SE","PVAL","EAF", 
                      "IMPUTED","IMPUTATION","HWE", "OR", "Neff"), collapse = ";") #Column names used in QC (DO NOT CHANGE)

#===INPUT - Reference panel 
REF_FILE_RSID   <- glue("../temp/ref/1KG_{ANCEST}_AF_maf{MAF}_rsid.txt"  ) #Pathname of reference panel
REF_FILE_ALLELE <- glue("../temp/ref/1KG_{ANCEST}_AF_maf{MAF}_allele.txt") #Pathname of reference panel
REF_COLS    <- paste(c("ChrPosID","REF","ALT",glue("AF_ALT_{ANCEST}_1000G")), collapse = ";") #Column names
REF_CLASSES <- paste(c(rep("character", 3), "numeric"), collapse=";") #Column types (CHARACTER;INTEGER;NUMERIC)
REF_MARKER  <- "ChrPosID"
REF_RSID    <- "RSID" #Column name of RSID column
REF_CHR     <- "Chr" #Column name of chromosome column
REF_POS     <- "Pos" #Column name of position column
REF_AF      <- glue("AF_ALT_{ANCEST}_1000G") #Column name of allele freq

#===OUTPUT
PATH_OUT    <- glue("../output/{ANCEST}/{PHENO}") # Destination
COLS_OUT    <- paste(c("SNP","cptid","Chr","position","EFFECT_ALLELE",
                       "OTHER_ALLELE","EAF",
                       "IMPUTED","IMPUTATION","BETA",
                       "SE","Z","PVAL","N","HWE",
                       "OR","Neff"), collapse = ";") 


#=Generate ECF File============================================================

ecf_file_binary <- glue(
  "#################################################################################################################
##### EasyQC-script to perform study-level and meta-level QC on imputed 1000G data
##### EasyQC version: 9.0
##### Programmer: Thomas Winkler, 2014-09-22
##### Contact: thomas.winkler@klinik.uni-regensburg.de
#################################################################################################################

### Please DEFINE here format and input columns of the following EASYIN files
DEFINE		--pathOut {PATH_OUT}
		--strMissing NA
		--strSeparator TAB
		--acolIn {COLS_IN}
		--acolInClasses {CLASSES_IN}
		--acolNewName {COLS_NEW}

## Please DO NOT CHANGE --acolNewName values because these reflect the column names used throughout the script
## If the study used different column names, please amend the respective value at --acolIn, the column will then
## be automatically renamed to the respective --acolNewName value

### Please define here the input files of the study:
EASYIN	--fileIn {FILE_IN}
	--fileInShortName {FILE_NAME}
	--astrSetNumCol MAF_THRESHOLD=0.{MAF};INFO_THRESHOLD=0.9;SDY=1

##############################################################
## 0. Start EasyQC scripting interface                      ##
##############################################################

START EASYQC

##############################################################
## 1. Sanity checks                                         ##
##############################################################

## Filtering  SNPs with missing values
CLEAN --rcdClean !(EFFECT_ALLELE%in%c('A','C','G','T')) --strCleanName numDrop_invalid_EA --blnWriteCleaned 0
CLEAN --rcdClean !(OTHER_ALLELE%in%c('A','C','G','T')) --strCleanName numDrop_invalid_OA --blnWriteCleaned 0
CLEAN --rcdClean is.na(SNP) --strCleanName numDrop_Missing_SNP --blnWriteCleaned 0
CLEAN --rcdClean is.na(EFFECT_ALLELE)&is.na(OTHER_ALLELE) --strCleanName numDrop_Missing_Both_Alleles --blnWriteCleaned 0
CLEAN --rcdClean is.na(PVAL) --strCleanName numDrop_Missing_P --blnWriteCleaned 0
CLEAN --rcdClean is.na(BETA) --strCleanName numDrop_Missing_BETA --blnWriteCleaned 0
CLEAN --rcdClean is.na(SE) --strCleanName numDrop_Missing_SE --blnWriteCleaned 0
CLEAN --rcdClean is.na(EAF) --strCleanName numDrop_Missing_EAF --blnWriteCleaned 0
CLEAN --rcdClean is.na(N) --strCleanName numDrop_Missing_N --blnWriteCleaned 0
CLEAN --rcdClean is.na(IMPUTED) --strCleanName numDrop_Missing_IMPUTED --blnWriteCleaned 0
CLEAN	--rcdClean IMPUTED==1 & is.na(IMPUTATION) --strCleanName numDrop_Imputed_Missing_INFO --blnWriteCleaned 0


## Filtering  SNPs with nonsense values
CLEAN --rcdClean PVAL<0|PVAL>1 --strCleanName numDrop_invalid_PVAL --blnWriteCleaned 0
CLEAN --rcdClean SE<=0|SE==Inf --strCleanName numDrop_invalid_SE --blnWriteCleaned 0
CLEAN --rcdClean abs(BETA)==Inf --strCleanName numDrop_invalid_BETA --blnWriteCleaned 0
CLEAN --rcdClean EAF<0|EAF>1 --strCleanName numDrop_invalid_EAF --blnWriteCleaned 0
CLEAN --rcdClean N<=0 --strCleanName numDrop_invalid_N --blnWriteCleaned 0
CLEAN --rcdClean IMPUTATION<0 & is.na(IMPUTATION)==0 --strCleanName numDrop_invalid_IMPUTATION --blnWriteCleaned 0
CLEAN --rcdClean IMPUTED!=0 & IMPUTED!=1 & is.na(IMPUTED)==0 --strCleanName numDrop_invalid_IMPUTED --blnWriteCleaned 0

## Calculate descriptives for imputation quality prior to QC
EVALSTAT    --colStat IMPUTATION
            --strTag Stat.PreQC

##############################################################
## 2. Main QC: apply general and cohort-specific thresholds ##
##############################################################

## Exclude monomorphic SNPs:
CLEAN --rcdClean (EAF==0)|(EAF==1) --strCleanName numDrop_Monomorphic --blnWriteCleaned 0

## MAF
CLEAN --rcdClean (EAF<MAF_THRESHOLD)|(EAF>(1-MAF_THRESHOLD)) --strCleanName numDrop_MAF --blnWriteCleaned 0
CALCULATE --rcdCalc mean(MAF_THRESHOLD) --strCalcName MAF_threshold

## Print the number of SNPs that has been imputed, genotyped, and used for imputation
GETNUM --rcdGetNum IMPUTED==1 --strGetNumName num_Imputed
GETNUM --rcdGetNum IMPUTED==0 --strGetNumName num_Genotyped
GETNUM --rcdGetNum is.na(IMPUTED) --strGetNumName num_ImputedNA

## IMPUTATION QUALITY: DO NOT FILTER ON IMPUTATION QUALITY IN CASE THE SNP HAS NOT BEEN IMPUTED (i.e., IMPUTED!=0 (does not work properly) USE: IMPUTED==1 or NA)
CLEAN --rcdClean is.na(IMPUTATION) & (IMPUTED==1 | is.na(IMPUTED)) --strCleanName numDrop_Imputed_ImpQual_Missing --blnWriteCleaned 0
CLEAN --rcdClean IMPUTATION<INFO_THRESHOLD & IMPUTED==1 --strCleanName numDrop_Imputed_ImpQual --blnWriteCleaned 0
CALCULATE --rcdCalc mean(INFO_THRESHOLD) --strCalcName ImpQual_threshold

## HWE: FILTER ONLY IF HWE INFO IS AVAILABLE AND THE MARKER IS KNOWN TO HAVE BEEN GENOTYPED
CLEAN --rcdClean HWE<=0.001 & is.na(HWE)==0 & N<1000 & IMPUTED==0 --strCleanName numDrop_Genotyped_HWElet1EM3 --blnWriteCleaned 0
CLEAN --rcdClean HWE<=0.0001 & is.na(HWE)==0 & N>=1000 & N<2000 & IMPUTED==0 --strCleanName numDrop_Genotyped_HWElet1EM4 --blnWriteCleaned 0
CLEAN --rcdClean HWE<=0.00001 & is.na(HWE)==0 & N>=2000 & N<10000 & IMPUTED==0 --strCleanName numDrop_Genotyped_HWElet1EM5 --blnWriteCleaned 0


##Add cptid column (colMap[Marker|Chr|Pos] are default values)
CREATECPTID
		--colInMarker SNP
		--colInA1 EFFECT_ALLELE
		--colInA2 OTHER_ALLELE
		--colInChr CHR
		--colInPos BP
		--fileMap {REF_FILE_RSID}
		--colMapMarker {REF_RSID}
		--colMapChr {REF_CHR}
		--colMapPos {REF_POS}
		--blnUseInMarker 0

##############################################################
## 3. Harmonization of allele coding (I/D)                  ##
##############################################################

## Aim: compile uniform allele codes A/C/G/T or I/D from different versions
HARMONIZEALLELES 	--colInA1 EFFECT_ALLELE
			--colInA2 OTHER_ALLELE


##############################################################
## 4.Filter duplicate SNPs                                  ##
##############################################################

## Aim: For duplicate SNPs keep only the duplicate with highest N

CLEANDUPLICATES	--colInMarker cptid
			--strMode removeall
			--colN N

## Dropped duplicates written to '*duplicates.txt'


##############################################################
## 5. AF Checks                                             ##
##############################################################

## Merge with file containing AFs for 1kG
MERGE 		--colInMarker cptid
		--fileRef {REF_FILE_ALLELE}
		--acolIn {REF_COLS}
		--acolInClasses {REF_CLASSES}
		--strRefSuffix .ref
		--colRefMarker {REF_MARKER}
		--blnWriteNotInRef 1
		--blnInAll 0
		--blnRefAll 0

## Create column with expected allele frequency, given imputation ref. sample
ADDCOL --rcdAddCol {REF_AF}.ref --colOut EAF.ref

## Adjust alleles to all be on the forward strand and to match the reference sample
## (this automatically adjust EAF s.t. the frequency corresponds to A1 in the reference sample)
## All mismatches will be removed (e.g. A/T in input, A/C in reference)
ADJUSTALLELES
			--colInA1 EFFECT_ALLELE
			--colInA2 OTHER_ALLELE
			--colInFreq EAF
			--colInBeta BETA
			--colRefA1 ALT.ref
			--colRefA2 REF.ref
			--blnRemoveMismatch 1
			--blnRemoveInvalid 1
			--blnWriteInvalid 1

## Compute difference expected and observed allele frequency
ADDCOL --rcdAddCol abs(EAF-EAF.ref) --colOut DAF

## Compare allele frequencies to those in the reference sample
## blnPlotAll 0 causes that only outlying SNPs with |Freq-Freq.ref|>0.2 will be plotted (way less computational time)
AFCHECK 	--colInFreq EAF
		--colRefFreq EAF.ref
		--numLimOutlier 0.2
		--blnPlotAll 0


##############################################################
## 6. Compute Pval(Z) and filter on |reported P - Pval(Z)|   #
##############################################################

## Create column with Z-scores
ADDCOL --rcdAddCol BETA/SE --colOut Z

## Compute P-values based on Z-scores and the difference of the log p-values
ADDCOL --rcdAddCol 2*pnorm(-abs(Z)) --colOut PVALZ
ADDCOL --rcdAddCol abs(log10(PVAL)-log10(PVALZ)) --colOut DLP

##############################################################
## 7. Approximate expected SE and plot versus reported SE   ##
##############################################################

## Compute expected SEs if phenotype standardized and not standardized
ADDCOL --rcdAddCol 1/sqrt(2*N*EAF*(1-EAF)) --colOut SE_EXPECT_PHEN_STD
ADDCOL --rcdAddCol SDY/sqrt(2*N*EAF*(1-EAF)) --colOut SE_EXPECT
ADDCOL --rcdAddCol sample(length(SE)) --colOut PERMUTATION
ADDCOL --rcdAddCol SE_EXPECT_PHEN_STD[PERMUTATION] --colOut SE_EXPECT_PHEN_STD_PERMUTED
ADDCOL --rcdAddCol SE_EXPECT[PERMUTATION] --colOut SE_EXPECT_PERMUTED
ADDCOL --rcdAddCol SE[PERMUTATION] --colOut SE_PERMUTED

# Concatenate the expected SE and the expected SE assuming the phenotype has been standardized
# Ensure that the SEs for which standardization has been assumed are printed black
# The counterpart (Y-axis) are the reported SEs
SPLOT  	--rcdSPlotX c(SE_EXPECT_PERMUTED[1:min(50000,length(SE))],SE_EXPECT_PHEN_STD_PERMUTED[1:min(50000,length(SE))])
              --rcdSPlotY c(SE_PERMUTED[1:min(50000,length(SE))],SE_PERMUTED[1:min(50000,length(SE))])
		--strDefaultColour red
		--arcdColourCrit (c(rep(0,min(50000,length(SE))),rep(1,min(50000,length(SE))))==1)
		--astrColour black
		--arcdAdd2Plot abline(0,1,col='orange')
              --strXlab EXPECTED SE
              --strYlab REPORTED SE
		--strAxes zeroequal
		--strPlotName SE-PLOT
		--strMode subplot

## Compute the expected SD(Y) and the input/reported SD(Y)
CALCULATE --rcdCalc mean(SE*sqrt(2*N*EAF*(1-EAF))) --strCalcName expSDY
CALCULATE --rcdCalc mean(SDY) --strCalcName repSDY


##############################################################
## 8. Rearrange columns and Write CLEANED output            ##
##############################################################

##Split created Chr and BP columns after removal of duplicates etc.
STRSPLITCOL --colSplit cptid
	--strSplit :
	--numSplitIdx 1
	--colOut Chr

STRSPLITCOL --colSplit cptid
		--strSplit :
		--numSplitIdx 2
		--colOut position

## Retain only columns relevant for exporting cleaned file
GETCOLS --acolOut {COLS_OUT}

## Write cleaned file
WRITE	--strPrefix CLEANED.
	--strMissing NA
	--strMode gz


##############################################################
## 9. Plot Pval(Z) versus reported P                        ##
##############################################################

PZPLOT	--colBeta BETA
	--colSe SE
	--colPval PVAL


##############################################################
## 10. Generate QQ plot                                     ##
##############################################################

QQPLOT	--acolQQPlot PVAL
	--numPvalOffset 0.05
	--strMode subplot


##############################################################
## 11. Generate GC PLOT                                     ##
##############################################################

CALCULATE --rcdCalc max(N,na.rm=T) --strCalcName N_max
GC	--colPval PVAL --blnSuppressCorrection 1

RPLOT	--rcdRPlotX N_max
	--rcdRPlotY Lambda.PVAL.GC
	--arcdAdd2Plot abline(h=1,col='orange');abline(h=1.1,col='red')
	--strAxes lim(0,NULL,0,NULL)
	--strPlotName GC-PLOT


##############################################################
## 12. Generate SE-N Plot                                   ##
##############################################################

CALCULATE --rcdCalc median(SE,na.rm=T) --strCalcName SEmedian
CALCULATE --rcdCalc median(1/sqrt(2*EAF*(1-EAF)), na.rm=T) --strCalcName c_trait_transf

RPLOT 	--rcdRPlotX sqrt(N_max)
	--rcdRPlotY c_trait_transf/SEmedian
	--arcdAdd2Plot abline(0,1,col='orange')
	--strAxes zeroequal
	--strPlotName SEN-PLOT


##############################################################
## 13. Calculate Descriptive statistics                     ##
##############################################################

## Get number of SNPs with large Beta and SE.
GETNUM --rcdGetNum SE>=5 --strGetNumName num_SEget5
GETNUM --rcdGetNum abs(BETA)>=5 --strGetNumName num_absBETAget5

## Get descriptives for BETA, SE, PVAL, and IMPUTATION
EVALSTAT	--colStat BETA
		--strTag Stat
EVALSTAT	--colStat SE
		--strTag Stat
EVALSTAT	--colStat PVAL
		--strTag Stat
EVALSTAT      --colStat IMPUTATION
              --strTag Stat.PostQC


STOP EASYQC")

#=Write text to .ecf file======================================================
system(glue("mkdir -p {PATH_OUT}"))
sink(glue("../temp/ecf_files/easyqc_ext2.0_{FILE_NAME}.ecf"))
cat(ecf_file_binary)
sink()

