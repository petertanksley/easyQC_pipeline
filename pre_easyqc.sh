#!/bin/bash
#SBATCH -J tt_pgs         # Job name
#SBATCH -o tt_pgs.o%j     # Name of stdout output file
#SBATCH -N 1		  # Total # of nodes (must be 1 for serial)
#SBATCH -n 24		  # Total # of mpi tasks (should be 1 for serial)
#SBATCH -p normal   	  # Queue (partition) name
#SBATCH -t 00:30:00	  # Run time (hh:mm:ss)
#SBATCH -A OTH21060	  # Project/Allocation name (req'd if you have more than 1)
#SBATCH --mail-type=all   # Send email at begin and end of job
#SBATCH --mail-user=peter.tanksley@austin.utexas.edu


#========================================================================================#
# SET WORKING ENVIRONMENT
#========================================================================================#

export PATH=$PATH:/work/07624/tankslpr/ls6/TOOLS

#=====EasyQC files===============================================#

#EasyQC expected headers
#rsID coded_all noncoded_all Chr position N Beta SE Pval AF_coded_all imputed oevar_imp HWE_pval OR Neff

#================================================================#
# AFR SUMSTATS
#================================================================#

#=========================================ALCP (Zhou et al., 2023)

#head -n2 ../input/AFR/ALCP/AUD_AFR_Aug2023.txt
#SNP_ID	Chrosome	Position	Allele1	Allele2	EA	EAF	SampleSize	Effect	SE	PValue	Direction
#rs3131972	1	752721	A	G	A	0.6529	105330.5	-0.00543256610890553	0.00457671955257416	0.2351	---

#assign variables
AFR_ALCP_RAW="../input/AFR/ALCP/AUD_AFR_Aug2023.txt"
AFR_ALCP_PRE="../temp/AFR/ALCP"
mkdir -p "$AFR_ALCP_PRE"

## REMOVE NON-SNPs, SNPS NOT INCLUDED IN GWAS, AND X-CHR ##
awk -F"\t" 'NR==1 {print} NR>1 && ($4 != "D" && $4 != "I") &&
	substr($1, 1, 2) == "rs" && ($9 != "" && $9 != "NA") &&
	($2 != "X") {print}' OFS="\t" \
	"$AFR_ALCP_RAW" > "$AFR_ALCP_PRE/afr_alcp_rsid_snps.txt"

#=NOTES
#combine chr and pos to form ChrPosID
#set oever_imp (INFO) to "NA"
#exponentiate beta to get OR
#divide beta by SE to get Z
#set N to "122571" (max sample size)
#set HWE_pval to "NA" (no data)
#set imputed to "0" (default to assuming no imputation)

## ADD FINAL HEADERS AND MAKE ADJUSTMENTS TO FIT STRUCTURE ##
awk -F"\t" '{if(NR == 1) {print "ChrPosID", "Chr", "rsID", "position", "coded_all", "noncoded_all",
	"AF_coded_all", "oevar_imp", "OR", "Beta", "SE", "Z", "Pval", "N", "Neff", "HWE_pval", "imputed"} \
	else {print $2":"$3, $2, $1, $3, $4, $5, $7, "NA", exp($9), $9, $10, $9/$10, $11, "122571", $8, "NA", 0}}' OFS="\t" \
	"$AFR_ALCP_PRE/afr_alcp_rsid_snps.txt" > "$AFR_ALCP_PRE/afr_alcp_PRE_EASYQC.txt"


#=======================================CUD (Johnson et al., 2020)

#head -n2 ../input/AFR/CUD/CUD_AFR_full_public_11.14.2020
#CHR SNP BP A1 A2 Z P N_CAS N_CON N
#1 rs111906919 16160196 T C 0.532 0.5944 2505 3751 6256

#1 CHR
#2 SNP
#3 BP
#4 A1
#5 A2
#6 Z
#7 P
#8 N_CAS
#9 N_CON
#10 N

#AFR_CUD_RAW="../input/AFR/CUD/CUD_AFR_full_public_11.14.2020"
#AFR_CUD_PRE="..temp/AFR/CUD"
#mkdir -p $AFR_CUD_PRE

## REMOVE NON-SNPs, SNPS NOT INCLUDED IN GWAS, AND X-CHR ##
#awk -F"\t" 'NR==1 {print} NR>1 && ($4 != "D" && $4 != "I") &&
#        substr($2, 1, 2) == "rs" && ($7 != "" && $7 != "NA") &&
#        ($1 != "X") {print}' OFS="\t" \
#        $AFR_CUD_RAW > $AFR_CUD_PRE/afr_cud_rsid.snps.txt

#NOTES
#combine chr and pos to form ChrPosID
#set oever_imp (INFO) to "NA"
#Calculate beta from Z: Z/sqrt((2MAF)*(1-MAF)*(N+Z^2))
#Calculate se from beta and Z: se=beta/Z
#exponentiate beta to get OR
#calculate Neff: SAMP_PREV=N_CAS/(N_CAS+N_CON); Neff=4*SAMP_PREV*(1-SAMP_PREV)*(N_CAS+N_CON)
#set HWE_pval to "NA" (no data)
#set imputed to "0" (default to assuming no imputation)

#awk -F"\t" '{if(NR == 1) {print "ChrPosID", "Chr", "rsID", "position", "coded_all", "noncoded_all",
#	"oevar_imp", "OR", "Beta", "SE", "Z", "Pval", "N", "Neff", "HWE_pval", "imputed"}
#	else if {print $1":"$3, $1, $2, $3, $4, $5, "NA", exp($9), $9, $10, $9/$10, $11, "122571", $8, "NA", 1}}' OFS="\t" \
#	$AFR_CUD_RAW > $AFR_CUD_PRE/afr_cud_PRE_EASYQC.txt


#====================================SMOK (Saunders et al., 2022)
#1-23andMe (females)
#2-23andMe (males)
#3-public (combined)

#23andMe - FEMALES

#head -n2 ../input/AFR/SMOK/GSCAN_SMOK_INI_FEMALE_AFR_23andme.txt
#rsid	chr	position	strand	A.OA	B.EA	B.EAF	effect	stderr	pvalue	N	imputed	gt.rate	hw.p.value
#rs28544273	1	751343	+	A	T	0.59031	0.0274084	0.0213779	0.199824	54757	1	NA	NA

AFR_SMOK_23F_RAW="../input/AFR/SMOK/GSCAN_SMOK_INI_FEMALE_AFR_23andme.txt"
AFR_SMOK_23F_PRE="../temp/AFR/SMOK"
mkdir -p "$AFR_SMOK_23F_PRE"

#Fix imputation field (currently, all rows show imputation)
awk -F "\t" 'BEGIN {OFS=FS} NR==1 || ($13 != "NA") { $12="0"} 1' \
	"$AFR_SMOK_23F_RAW" > "$AFR_SMOK_23F_PRE/afr_smok_23f_impfix.txt"

# REMOVE NON-SNPs, SNPS NOT INCLUDED IN GWAS, AND X-CHR ##
awk -F"\t" 'NR==1 {print} NR>1 && ($6 != "D" && $6 != "I") &&
        substr($1, 1, 2) == "rs" && ($10 != "" && $10 != "NA") &&
        ($2 != "X") {print}' OFS="\t" \
        "$AFR_SMOK_23F_PRE/afr_smok_23f_impfix.txt" > "$AFR_SMOK_23F_PRE/afr_smok_23f_impfix_rsid_snps.txt"

#NOTES
#combine chr and pos to form ChrPosID
#set oever_imp (INFO) to "NA"
#Calculate se from beta and Z: se=beta/Z
#exponentiate beta to get OR
#
#Calculate Neff:
#=N_CAS: 54757*.23=12594
#=N_CON: 54757-12594=42163
#=SAMP_PREV: 12594/54757=.23
#=Neff=4*.23*(1-.23)*(54757)=38789.86

awk -F"\t" '{if(NR == 1) {print "ChrPosID", "Chr", "rsID", "position", "coded_all", "noncoded_all",
	"oevar_imp", "OR", "Beta", "SE", "Z", "Pval", "N", "Neff", "HWE_pval", "imputed"} \
	else {print $2":"$3, $2, $1, $3, $5, $6, "NA", exp($8), $8, $9, $8/$9, $10, $11, "38789.86", $14, $12}}' OFS="\t" \
	"$AFR_SMOK_23F_PRE/afr_smok_23f_impfix_rsid_snps.txt" > "$AFR_SMOK_23F_PRE/afr_smok_23f_PRE_EASYQC.txt"

#23andMe - MALES

#head -n2 ../input/AFR/SMOK/GSCAN_SMOK_INI_MALE_AFR_23andme.txt
#rsid	chr	position	strand	A.OA	B.EA	B.EAF	effect	stderr	pvalue	N	imputed	gt.rate	hw.p.value
#rs28544273	1	751343	+	A	T	0.59031	0.019607	0.0228292	0.390432882	40554	1	NA	NA

AFR_SMOK_23M_RAW="../input/AFR/SMOK/GSCAN_SMOK_INI_MALE_AFR_23andme.txt"
AFR_SMOK_23M_PRE="../temp/AFR/SMOK"
mkdir -p "$AFR_SMOK_23M_PRE"

#Fix imputation field (currently, all rows show imputation)
awk -F "\t" 'BEGIN {OFS=FS} NR==1 || ($13 != "NA") { $12="0"} 1' \
	"$AFR_SMOK_23M_RAW" > "$AFR_SMOK_23M_PRE/afr_smok_23m_impfix.txt"

# REMOVE NON-SNPs, SNPS NOT INCLUDED IN GWAS, AND X-CHR ##
awk -F"\t" 'NR==1 {print} NR>1 && ($6 != "D" && $6 != "I") &&
        substr($1, 1, 2) == "rs" && ($10 != "" && $10 != "NA") &&
        ($2 != "X") {print}' OFS="\t" \
	"$AFR_SMOK_23M_PRE/afr_smok_23m_impfix.txt" > "$AFR_SMOK_23M_PRE/afr_smok_23m_impfix_rsid_snps.txt"

#NOTES
#combine chr and pos to form ChrPosID
#set oever_imp (INFO) to "NA"
#Calculate se from beta and Z: se=beta/Z
#exponentiate beta to get OR
#
#Calculate Neff:
#=N_CAS: 40554*.29=11720
#=N_CON: 40554-11720=28834
#=SAMP_PREV: 11720/40554=.29
#=Neff=4*.29*(1-.29)*(40554)=33400.27

awk -F"\t" '{if(NR == 1) {print "ChrPosID", "Chr", "rsID", "position", "coded_all", "noncoded_all",
	"oevar_imp", "OR", "Beta", "SE", "Z", "Pval", "N", "Neff", "HWE_pval", "imputed"} \
        else {print $2":"$3, $2, $1, $3, $5, $6, "NA", exp($8), $8, $9, $8/$9, $10, $11, "33400.27", $14, $12}}' OFS="\t" \
        "$AFR_SMOK_23M_PRE/afr_smok_23m_impfix_rsid_snps.txt" > "$AFR_SMOK_23M_PRE/afr_smok_23m_PRE_EASYQC.txt"


#Public - COMBINED (GRCH38)

#head -n2 ../input/AFR/SMOK/GSCAN_SmkInit_2022_GWAS_SUMMARY_STATS_AFR.txt
#CHR	POS	RSID	EFFECT_ALLELE	OTHER_ALLELE	AF_1000G	BETA	SE	P	N
#chr10	100000222	rs138880521	A	G	0.0189107	-0.0355	0.03	0.237594	24278

AFR_SMOK_PUB_RAW="../input/AFR/SMOK/GSCAN_SmkInit_2022_GWAS_SUMMARY_STATS_AFR.txt"
AFR_SMOK_PUB_PRE="../temp/AFR/SMOK"
mkdir -p "$AFR_SMOK_PUB_PRE"

# REMOVE NON-SNPs, SNPS NOT INCLUDED IN GWAS, AND X-CHR ##
awk -F"\t" 'NR==1 {print} NR>1 && ($4 != "D" && $4 != "I") &&
        substr($3, 1, 2) == "rs" && ($9 != "" && $9 != "NA") &&
        ($2 != "X") {print}' OFS="\t" \
        "$AFR_SMOK_PUB_RAW" > "$AFR_SMOK_PUB_PRE/afr_smok_pub_rsid_snps.txt"

#NOTES
#combine chr and pos to form ChrPosID
#set oever_imp (INFO) to "NA"
#Calculate se from beta and Z: se=beta/Z
#exponentiate beta to get OR
#
#Calculate Neff (prev reported in paper):
#=N_CAS: 24278*.41=9915.72
#=N_CON: 24278-9915.72=14362.28
#=SAMP_PREV: 9915.72/24278=.41
#=Neff=4*.41*(1-.41)*(24278)=23491.39

awk -F"\t" '{if(NR == 1) {print "ChrPosID", "Chr", "rsID", "position", "coded_all", "noncoded_all",
	"oevar_imp", "OR", "Beta", "SE", "Z", "Pval", "N", "Neff", "HWE_pval", "imputed"} \
        else {print $1":"$2, $1, $3, $2, $4, $5, "NA", exp($7), $7, $8, $7/$8, $9, $10, "23491.39", "NA", "0"}}' OFS="\t" \
        "$AFR_SMOK_PUB_PRE/afr_smok_pub_rsid_snps.txt" > "$AFR_SMOK_PUB_PRE/afr_smok_pub_PRE_EASYQC.txt"




#================================================================#
# EUR SUMSTATS
#================================================================#

#====================================ADHD (Demontis et al., 2022)

#head -n2 ../input/EUR/ADHD/ADHD2022_iPSYCH_deCODE_PGC.meta
#CHR SNP BP A1 A2 FRQ_A_38691 FRQ_U_186843 INFO OR SE P Direction Nca Nco
#8 rs62513865 101592213 C T 0.925 0.937 0.981 0.99631 0.0175 0.8325 +---+++0-++-+ 38691 186843

EUR_ADHD_RAW="../input/EUR/ADHD/ADHD2022_iPSYCH_deCODE_PGC.meta"
EUR_ADHD_PRE="../temp/EUR/ADHD"
mkdir -p "$EUR_ADHD_PRE"

# REMOVE NON-SNPs, SNPS NOT INCLUDED IN GWAS, AND X-CHR ##

awk -F'[[:space:]]' 'BEGIN {OFS="\t"} NR==1 {print} NR>1 && ($4 != "D" && $4 != "I") && \
    substr($2, 1, 2) == "rs" && ($11 != "" && $11 != "NA") && \
    ($1 != "X") {print}' "$EUR_ADHD_RAW" > "$EUR_ADHD_PRE/eur_adhd_rsid_snps.txt"


#NOTES
#combine chr and pos to form ChrPosID
#Calculate Z: Z=log(OR)/se
#get beta from log(OR)
#
#Calculate Neff (prev reported in paper):
#=N_CAS: 38691
#=N_CON: 186843
#=SAMP_PREV: 38691/(38691+186843)=.17
#=Neff=4*.17*(1-.17)*(225534)=127291.39

awk -F"\t" 'BEGIN {print "ChrPosID", "Chr", "rsID", "position", "coded_all", "noncoded_all",
	"oevar_imp", "OR", "Beta", "SE", "Z", "Pval", "N", "Neff", "HWE_pval", "imputed"} \
	NR > 1 {Beta=log($9); \
	Z=($10 != 0) ? log($9)/$10 : "NA"; \
	print $1 ":" $3, $1, $2, $3, $4, $5, $8, $9, Beta, $10, Z, $11, "225534", "127291.39", "NA", "1"}' \
	"$EUR_ADHD_PRE/eur_adhd_rsid_snps.txt" > "$EUR_ADHD_PRE/eur_adhd_PRE_EASYQC.txt"


#========================================ALCP (Zhou et al., 2023)

#head -n2 ../input/EUR/ALCP/PAU_EUR_Aug2023.txt
#SNP_ID	Chrosome	Position	Allele1	Allele2	EA	EAF	SampleSize	Effect	SE	PValue	Direction
#rs144155419	1	717587	A	G	A	0.0100	183925.9	-0.0171839676620293	0.0165708463471835	0.2996	?-+-???????

#assign variables
EUR_ALCP_RAW="../input/EUR/ALCP/PAU_EUR_Aug2023.txt"
EUR_ALCP_PRE="../temp/EUR/ALCP"
mkdir -p "$EUR_ALCP_PRE"


## REMOVE NON-SNPs, SNPS NOT INCLUDED IN GWAS, AND X-CHR ##
awk -F"\t" 'NR==1 {print} NR>1 && ($4 != "D" && $4 != "I") && \
	substr($1, 1, 2) == "rs" && ($11 != "" && $11 != "NA") && \
	($2 != "X") {print}' OFS="\t" \
	"$EUR_ALCP_RAW" > "$EUR_ALCP_PRE/eur_alcp_rsid_snps.txt"

#=NOTES
#combine chr and pos to form ChrPosID
#set oever_imp (INFO) to "NA"
#exponentiate beta to get OR
#divide beta by SE to get Z
#set N to "122571" (max sample size)
#set HWE_pval to "NA" (no data)
#set imputed to "0" (default to assuming no imputation)

## ADD FINAL HEADERS AND MAKE ADJUSTMENTS TO FIT STRUCTURE ##
awk -F"\t" '{if(NR == 1) {print "ChrPosID", "Chr", "rsID", "position", "coded_all", "noncoded_all", \
	"oevar_imp", "OR", "Beta", "SE", "Z", "Pval", "N", "Neff", "HWE_pval", "imputed"} \
	else {print $2":"$3, $2, $1, $3, $4, $5, $7, "NA", exp($9), $9, $10, $9/$10, $11, "122571", $8, "NA", 0}}' OFS="\t" \
	"$EUR_ALCP_PRE/eur_alcp_rsid_snps.txt" > "$EUR_ALCP_PRE/eur_alcp_PRE_EASYQC.txt"


#======================================CUD (Johnson et al., 2020)

#head -n2 ../input/EUR/CUD/CUD_EUR_full_noAddHealth.chrbp.short
#CHR SNP BP A1 A2 Z P N N_CAS N_CON
#1 rs12562373 166367755 A G -2.001 0.04537 362864 15170 347694

#====================================SMOK (Saunders et al., 2022)

#1-23andMe (females)

#head -n2 ../input/EUR/SMOK/GSCAN_SMOK_INI_FEMALE_EUR_23andme.txt
#rsid	chr	position	strand	A.OA	B.EA	B.EAF	effect	stderr	pvalue	N	imputed	gt.rate	hw.p.value
#rs3131972	1	752721	+	A	G	0.82494	0.00127845	0.00399547	0.748967	997356	1	0.9908	0

EUR_SMOK_23F_RAW="../input/EUR/SMOK/GSCAN_SMOK_INI_FEMALE_EUR_23andme.txt"
EUR_SMOK_23F_PRE="../temp/EUR/SMOK"
mkdir -p "$EUR_SMOK_23F_PRE"

#Fix imputation field (currently, all rows show imputation)
awk -F "\t" 'BEGIN {OFS=FS} NR==1 || ($13 != "NA") { $12="0"} 1' \
        "$EUR_SMOK_23F_RAW" > "$EUR_SMOK_23F_PRE/eur_smok_23f_impfix.txt"

# REMOVE NON-SNPs, SNPS NOT INCLUDED IN GWAS, AND X-CHR ##
awk -F"\t" 'NR==1 {print} NR>1 && ($6 != "D" && $6 != "I") &&
        substr($1, 1, 2) == "rs" && ($10 != "" && $10 != "NA") &&
        ($2 != "X") {print}' OFS="\t" \
        "$EUR_SMOK_23F_PRE/eur_smok_23f_impfix.txt" > "$EUR_SMOK_23F_PRE/eur_smok_23f_impfix_rsid_snps.txt"

#NOTES
#combine chr and pos to form ChrPosID
#set oever_imp (INFO) to "NA"
#Calculate se from beta and Z: se=beta/Z
#exponentiate beta to get OR
#
#Calculate Neff (prev provided by article):
#=N_CAS: 997356*.41=408916
#=N_CON: 997356-408916=588440
#=SAMP_PREV: 408916/997356=.41
#=Neff=4*.41*(1-.41)*(997356)=965041.67

awk -F"\t" '{if(NR == 1) {print "ChrPosID", "Chr", "rsID", "position", "coded_all", "noncoded_all",
	"oevar_imp", "OR", "Beta", "SE", "Z", "Pval", "N", "Neff", "HWE_pval", "imputed"} \
        else {print $2":"$3, $2, $1, $3, $5, $6, "NA", exp($8), $8, $9, $8/$9, $10, $11, "965041.67", $14, $12}}' OFS="\t" \
        "$EUR_SMOK_23F_PRE/eur_smok_23f_impfix_rsid_snps.txt" > "$EUR_SMOK_23F_PRE/eur_smok_23f_PRE_EASYQC.txt"

#2-23andMe (males)


EUR_SMOK_23M_RAW="../input/EUR/SMOK/GSCAN_SMOK_INI_MALE_EUR_23andme.txt"
EUR_SMOK_23M_PRE="../temp/EUR/SMOK"
mkdir -p "$EUR_SMOK_23M_PRE"

#Fix imputation field (currently, all rows show imputation)
awk -F "\t" 'BEGIN {OFS=FS} NR==1 || ($13 != "NA") { $12="0"} 1' \
        "$EUR_SMOK_23M_RAW" > "$EUR_SMOK_23M_PRE/eur_smok_23m_impfix.txt"

# REMOVE NON-SNPs, SNPS NOT INCLUDED IN GWAS, AND X-CHR ##
awk -F"\t" 'NR==1 {print} NR>1 && ($6 != "D" && $6 != "I") &&
        substr($1, 1, 2) == "rs" && ($10 != "" && $10 != "NA") &&
        ($2 != "X") {print}' OFS="\t" \
        "$EUR_SMOK_23M_PRE/eur_smok_23m_impfix.txt" > "$EUR_SMOK_23M_PRE/eur_smok_23m_impfix_rsid_snps.txt"

#NOTES
#combine chr and pos to form ChrPosID
#set oever_imp (INFO) to "NA"
#Calculate se from beta and Z: se=beta/Z
#exponentiate beta to get OR
#
#Calculate Neff (prev provided by article):
#=N_CAS: 866242*.43=372484
#=N_CON: 866242-372484=493758
#=SAMP_PREV: 372484/866242=.43
#=Neff=4*.43*(1-.43)*(866242)=849263.66

awk -F"\t" '{if(NR == 1) {print "ChrPosID", "Chr", "rsID", "position", "coded_all", "noncoded_all",
	"oevar_imp", "OR", "Beta", "SE", "Z", "Pval", "N", "Neff", "HWE_pval", "imputed"} \
        else {print $2":"$3, $2, $1, $3, $5, $6, "NA", exp($8), $8, $9, $8/$9, $10, $11, "849263.66", $14, $12}}' OFS="\t" \
        "$EUR_SMOK_23M_PRE/eur_smok_23m_impfix_rsid_snps.txt" > "$EUR_SMOK_23M_PRE/eur_smok_23m_PRE_EASYQC.txt"

#3-public

#head -n2 ../input/EUR/SMOK/GSCAN_SmkInit_2022_GWAS_SUMMARY_STATS_EUR_without_UKB.txt
#CHR	POS	RSID	EFFECT_ALLELE	OTHER_ALLELE	AF_1000G	BETA	SE	P	N
#chr10	100000235	rs11596870	T	C	0.314115	-0.000605	0.003	0.811963	357235

EUR_SMOK_PUB_RAW="../input/EUR/SMOK/GSCAN_SmkInit_2022_GWAS_SUMMARY_STATS_EUR_without_UKB.txt"
EUR_SMOK_PUB_PRE="../temp/EUR/SMOK"
mkdir -p "$EUR_SMOK_PUB_PRE"

#Fix imputation field (currently, all rows show imputation)
awk 'BEGIN {FS=OFS="\t"} NR > 1 {$1 = substr($1, 4)} 1' \
        "$EUR_SMOK_PUB_RAW" > "$EUR_SMOK_PUB_PRE/eur_smok_pub_chrfix.txt"

# REMOVE NON-SNPs, SNPS NOT INCLUDED IN GWAS, AND X-CHR ##
awk -F"\t" 'NR==1 {print} NR>1 && ($4 != "D" && $4 != "I") &&
        substr($3, 1, 2) == "rs" && ($9 != "" && $9 != "NA") &&
        ($2 != "X") {print}' OFS="\t" \
        "$EUR_SMOK_PUB_PRE/eur_smok_pub_chrfix.txt" > "$EUR_SMOK_PUB_PRE/eur_smok_pub_chrfix_rsid_snps.txt"

#NOTES
#combine chr and pos to form ChrPosID
#set oever_imp (INFO) to "NA"
#Calculate se from beta and Z: se=beta/Z
#exponentiate beta to get OR
#
#Calculate Neff (prev reported in paper):
#=N_CAS: 24278*.41=9915.72
#=N_CON: 24278-9915.72=14362.28
#=SAMP_PREV: 9915.72/24278=.41
#=Neff=4*.41*(1-.41)*(24278)=23491.39

awk -F"\t" '{if(NR == 1) {print "ChrPosID", "Chr", "rsID", "position", "coded_all", "noncoded_all",
	"oevar_imp", "OR", "Beta", "SE", "Z", "Pval", "N", "Neff", "HWE_pval", "imputed"} \
        else {print $1":"$2, $1, $3, $2, $4, $5, "NA", exp($7), $7, $8, $7/$8, $9, $10, "23491.39", "NA", "0"}}' OFS="\t" \
        "$EUR_SMOK_PUB_PRE/eur_smok_pub_rsid_snps.txt" > "$EUR_SMOK_PUB_PRE/eur_smok_pub_PRE_EASYQC.txt"

