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

export PATH="$PATH:/work/07624/tankslpr/ls6/TOOLS"
CHAIN_38to19="/work/07624/tankslpr/ls6/TOOLS/liftOver_files/hg38ToHg19.over.chain"

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

if [ ! -f "$AFR_ALCP_PRE/afr_alcp_PRE_EASYQC.txt" ]; then

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
fi

#=======================================CUD (Johnson et al., 2020)


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

if [ ! -f "$AFR_SMOK_23F_PRE/afr_smok_23f_PRE_EASYQC.txt" ]; then

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

#CHANGED LAST TO COLUMNS TO ("NA", 1)
#THIS IGNORES THE "IMPFIX" FILE MODIFICATIONS

awk -F"\t" '{if(NR == 1) {print "ChrPosID", "Chr", "rsID", "position", "coded_all", "noncoded_all",
	"AF_coded_all", "oevar_imp", "OR", "Beta", "SE", "Z", "Pval", "N", "Neff", "HWE_pval", "imputed"} \
	else {print $2":"$3, $2, $1, $3, $5, $6, $7, "NA", exp($8), $8, $9, $8/$9, $10, $11, "38789.86", "NA", 0}}' OFS="\t" \
	"$AFR_SMOK_23F_PRE/afr_smok_23f_impfix_rsid_snps.txt" > "$AFR_SMOK_23F_PRE/afr_smok_23f_PRE_EASYQC.txt"

fi

#23andMe - MALES

#head -n2 ../input/AFR/SMOK/GSCAN_SMOK_INI_MALE_AFR_23andme.txt
#rsid	chr	position	strand	A.OA	B.EA	B.EAF	effect	stderr	pvalue	N	imputed	gt.rate	hw.p.value
#rs28544273	1	751343	+	A	T	0.59031	0.019607	0.0228292	0.390432882	40554	1	NA	NA

AFR_SMOK_23M_RAW="../input/AFR/SMOK/GSCAN_SMOK_INI_MALE_AFR_23andme.txt"
AFR_SMOK_23M_PRE="../temp/AFR/SMOK"
mkdir -p "$AFR_SMOK_23M_PRE"

if [ ! -f "$AFR_SMOK_23M_PRE/afr_smok_23m_PRE_EASYQC.txt" ]; then

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
	"AF_coded_all", "oevar_imp", "OR", "Beta", "SE", "Z", "Pval", "N", "Neff", "HWE_pval", "imputed"} \
        else {print $2":"$3, $2, $1, $3, $5, $6, $7, "NA", exp($8), $8, $9, $8/$9, $10, $11, "33400.27", "NA", 0}}' OFS="\t" \
        "$AFR_SMOK_23M_PRE/afr_smok_23m_impfix_rsid_snps.txt" > "$AFR_SMOK_23M_PRE/afr_smok_23m_PRE_EASYQC.txt"
fi

#Public - COMBINED (GRCH38)

#head -n2 ../input/AFR/SMOK/GSCAN_SmkInit_2022_GWAS_SUMMARY_STATS_AFR.txt
#CHR	POS	RSID	EFFECT_ALLELE	OTHER_ALLELE	AF_1000G	BETA	SE	P	N
#chr10	100000222	rs138880521	A	G	0.0189107	-0.0355	0.03	0.237594	24278

AFR_SMOK_PUB_RAW="../input/AFR/SMOK/GSCAN_SmkInit_2022_GWAS_SUMMARY_STATS_AFR.txt"
AFR_SMOK_PUB_PRE="../temp/AFR/SMOK"
mkdir -p "$AFR_SMOK_PUB_PRE"

if [ ! -f "$AFR_SMOK_PUB_PRE/afr_smok_pub_PRE_EASYQC.txt" ]; then

#prep file for liftOver (38->19)
awk 'BEGIN {OFS="\t"} NR>1 {print $1, $2-1, $2, $0}' "$AFR_SMOK_PUB_RAW" > "$AFR_SMOK_PUB_PRE/afr_smok_pub_liftOver.bed"

#liftOver to HG19
liftOver -bedPlus=3 \
	-tab \
	"$AFR_SMOK_PUB_PRE/afr_smok_pub_liftOver.bed" \
	"$CHAIN_38to19" \
	"$AFR_SMOK_PUB_PRE/afr_smok_pub_hg19.bed" \
	"$AFR_SMOK_PUB_PRE/afr_smok_pub_unmapped38.bed"

#head -n2 ../temp/AFR/SMOK/afr_smok_pub_hg19.bed
#chr10	101759978	101759979	chr10	100000222	rs138880521	A	G	0.0189107	-0.0355	0.03	0.237594	24278
#chr10	101759991	101759992	chr10	100000235	rs11596870	T	C	0.30938	0.00543	0.01	0.582511	24278

#remove extra columns (2&4) and remove "chr" prefix from column
awk 'BEGIN {OFS="\t"; print "CHR", "POS", "RSID", "EFFECT_ALLELE", "OTHER_ALLELE", "AF_1000G", "BETA", "SE", "P", "N"} \
	{sub(/^chr/, "", $1); print $1, $3, $6, $7, $8, $9, $10, $11, $12, $13}' \
	"$AFR_SMOK_PUB_PRE/afr_smok_pub_hg19.bed" > "$AFR_SMOK_PUB_PRE/afr_smok_pub_hg19.txt"

# REMOVE NON-SNPs, SNPS NOT INCLUDED IN GWAS, AND X-CHR ##
awk -F"\t" 'NR==1 {print} NR>1 && ($4 != "D" && $4 != "I") &&
        substr($3, 1, 2) == "rs" && ($9 != "" && $9 != "NA") &&
        ($2 != "X") {print}' OFS="\t" \
        "$AFR_SMOK_PUB_PRE/afr_smok_pub_hg19.txt" > "$AFR_SMOK_PUB_PRE/afr_smok_pub_hg19_rsid_snps.txt"

#ATTACH AF AND ALT ALLELE FROM REF PANEL
awk -v OFS='\t' 'NR==FNR{a[$3]=$6 "\t" $7;next} FNR==1{print $0, "ALT", "AF_ALT_AFR_1000G"; next} $3 in a{print $0, a[$3];}' \
    ../temp/ref/1KG_AFR_AF_maf001.txt "$AFR_SMOK_PUB_PRE/afr_smok_pub_hg19_rsid_snps.txt" \
    > "$AFR_SMOK_PUB_PRE/afr_smok_pub_hg19_rsid_snps_af.txt"

# ALIGN AF_coded_all FROM THE ATTACHED COLUMNS ##
awk -F"\t" 'BEGIN {OFS="\t"} {
    if (NR==1) {
        print $0, "AF_coded_all"
    } else if (NR > 1 && $4 == $11) {
        print $0, $12
    } else if (NR > 1 && $5 == $11) {
        print $0, 1 - $12
    }
}' "$AFR_SMOK_PUB_PRE/afr_smok_pub_hg19_rsid_snps_af.txt" > "$AFR_SMOK_PUB_PRE/afr_smok_pub_hg19_rsid_snps_af_aligned.txt"

#head -n2 ../temp/AFR/SMOK/afr_smok_pub_hg19_rsid_snps_af_aligned.txt
#CHR	POS	RSID	EFFECT_ALLELE	OTHER_ALLELE	AF_1000G	BETA	SE	P	N	ALT	AF_ALT_AFR_1000G	AF_coded_all
#10	101759979	rs138880521	A	G	0.0189107	-0.0355	0.03	0.237594	24278	A 0.0192012

#NOTES
#combine chr and pos to form ChrPosID
#set oever_imp (INFO) to "NA"
#Calculate se from beta and Z: se=beta/Z
#exponentiate beta to get OR
#
#Calculate Neff as sum of effective sample size (cohort-specific prev reported in supp)
#Neff==22109.62

awk -F"\t" '{if(NR == 1) {print "ChrPosID", "Chr", "rsID", "position", "coded_all", "noncoded_all",
	"AF_coded_all", "oevar_imp", "OR", "Beta", "SE", "Z", "Pval", "N", "Neff", "HWE_pval", "imputed"} \
        else {print $1":"$2, $1, $3, $2, $4, $5, $13, "NA", exp($7), $7, $8, $7/$8, $9, $10, "22109.62", "NA", 0}}' OFS="\t" \
        "$AFR_SMOK_PUB_PRE/afr_smok_pub_hg19_rsid_snps_af_aligned.txt" > "$AFR_SMOK_PUB_PRE/afr_smok_pub_PRE_EASYQC.txt"

fi



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

if [ ! -f "$EUR_ADHD_PRE/eur_adhd_PRE_EASYQC.txt" ]; then

#change delimiter from space to tab
cat "$EUR_ADHD_RAW" | tr ' ' '\t' > "$EUR_ADHD_PRE/eur_adhd_tab.txt"

# REMOVE NON-SNPs, SNPS NOT INCLUDED IN GWAS, AND X-CHR ##
awk -F'\t' 'BEGIN {OFS="\t"} NR==1 {print} NR>1 && ($4 != "D" && $4 != "I") && \
    substr($2, 1, 2) == "rs" && ($11 != "" && $11 != "NA") && \
    ($1 != "X") {print}' "$EUR_ADHD_PRE/eur_adhd_tab.txt" > "$EUR_ADHD_PRE/eur_adhd_tab_rsid_snps.txt"

#head -n2 ../temp/ref/1KG_EUR_AF_maf005.txt
#ChrPosID	Chr	RSID	Pos	REF	ALT	AF_ALT_EUR_1000G
#1:11008	1	rs575272151	11008	C	G	0.0884692

#ATTACH AF AND ALT ALLELE FROM REF PANEL
awk -v OFS='\t' 'NR==FNR{a[$3]=$6 "\t" $7;next} FNR==1{print $0, "ALT", "AF_ALT_AFR_1000G"; next} $2 in a{print $0, a[$2];}' \
    ../temp/ref/1KG_EUR_AF_maf005.txt "$EUR_ADHD_PRE/eur_adhd_tab_rsid_snps.txt" \
    > "$EUR_ADHD_PRE/eur_adhd_tab_rsid_snps_af.txt"

# ALIGN AF_coded_all FROM THE ATTACHED COLUMNS ##
awk -F"\t" 'BEGIN {OFS="\t"} {
    if (NR==1) {
        print $0, "AF_coded_all"
    } else if (NR > 1 && $4 == $15) {
        print $0, $16
    } else if (NR > 1 && $5 == $15) {
        print $0, 1 - $16
    }
}' "$EUR_ADHD_PRE/eur_adhd_tab_rsid_snps_af.txt" > "$EUR_ADHD_PRE/eur_adhd_tab_rsid_snps_af_aligned.txt"

#head -n2 ../temp/EUR/ADHD/eur_adhd_tab_rsid_snps.txt
#CHR	SNP	BP	A1	A2	FRQ_A_38691	FRQ_U_186843	INFO	OR	SE	P	Direction	Nca	Nco
#8	rs62513865	101592213	C	T	0.925	0.937	0.981	0.99631	0.0175	0.8325	+---+++0-++-+	38691	186843

#NOTES
#combine chr and pos to form ChrPosID
#Calculate Z: Z=log(OR)/se
#get beta from log(OR)
#
#Calculate Neff as sum of effective sample size (cohort-specific prev reported in supp):
#Neff==103135.53


awk -F"\t" 'BEGIN {OFS="\t"} NR==1 {print "ChrPosID", "Chr", "rsID", "position", "coded_all", "noncoded_all",
	"AF_coded_all", "oevar_imp", "OR", "Beta", "SE","Z","Pval", "N", "Neff", "HWE_pval", "imputed"}
	NR>1 {print $1 ":" $3, $1, $2, $3, $4, $5, $17, $8, $9, log($9), $10,log($9)/$10,$11, "225534", "103135.53", "NA", 1}' \
	"$EUR_ADHD_PRE/eur_adhd_tab_rsid_snps_af_aligned.txt" > "$EUR_ADHD_PRE/eur_adhd_PRE_EASYQC.txt"
fi

#========================================ALCP (Zhou et al., 2023)

#head -n2 ../input/EUR/ALCP/PAU_EUR_Aug2023.txt
#SNP_ID	Chrosome	Position	Allele1	Allele2	EA	EAF	SampleSize	Effect	SE	PValue	Direction
#rs144155419	1	717587	A	G	A	0.0100	183925.9	-0.0171839676620293	0.0165708463471835	0.2996	?-+-???????

#assign variables
EUR_ALCP_RAW="../input/EUR/ALCP/PAU_EUR_Aug2023.txt"
EUR_ALCP_PRE="../temp/EUR/ALCP"
mkdir -p "$EUR_ALCP_PRE"

if [ ! -f "$EUR_ALCP_PRE/eur_alcp_PRE_EASYQC.txt" ]; then

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
#set N to "359398.05" (max sample size)
#set HWE_pval to "NA" (no data)
#set imputed to "0" (default to assuming no imputation)

## ADD FINAL HEADERS AND MAKE ADJUSTMENTS TO FIT STRUCTURE ##
awk -F"\t" '{if(NR == 1) {print "ChrPosID", "Chr", "rsID", "position", "coded_all", "noncoded_all", \
	"AF_coded_all", "oevar_imp", "OR", "Beta", "SE", "Z", "Pval", "N", "Neff", "HWE_pval", "imputed"} \
	else {print $2":"$3, $2, $1, $3, $4, $5, $7, "NA", exp($9), $9, $10, $9/$10, $11, "903147", "359398.05", "NA", 0}}' OFS="\t" \
	"$EUR_ALCP_PRE/eur_alcp_rsid_snps.txt" > "$EUR_ALCP_PRE/eur_alcp_PRE_EASYQC.txt"
fi


#========================================ALCP (UKB - no sibs)

#head -n2 ../input/EUR/ALCP/CLEANED.UKB_L_AUDIT_P_log10
#SNP    cptid   Chr     position        EFFECT_ALLELE   OTHER_ALLELE    EAF     IMPUTED IMPUTATION      BETA    SE      Z       PVAL    N       HWE
#rs113190524    10:62093307     10      62093307        C       T       0.719602        1       0.994363        0.000261662     0.00111148      0.235417641343074       0.81    130999  NA

#assign variables
EUR_ALCP_UKB_RAW="../input/EUR/ALCP/CLEANED.UKB_L_AUDIT_P_log10"
EUR_ALCP_UKB_PRE="../temp/EUR/ALCP"
mkdir -p "$EUR_ALCP_UKB_PRE"

if [ ! -f "$EUR_ALCP_UKB_PRE/eur_alcp_ukb_PRE_EASYQC.txt" ]; then

## REMOVE NON-SNPs, SNPS NOT INCLUDED IN GWAS, AND X-CHR ##
awk -F"\t" 'NR==1 {print} NR>1 && ($5 != "D" && $5 != "I") && \
        substr($1, 1, 2) == "rs" && ($13 != "" && $13 != "NA") && \
        ($2 != "X") {print}' OFS="\t" \
        "$EUR_ALCP_RAW" > "$EUR_ALCP_PRE/eur_alcp_rsid_snps.txt"

#=NOTES

## ADD FINAL HEADERS AND MAKE ADJUSTMENTS TO FIT STRUCTURE ##
awk -F"\t" '{if(NR == 1) {print "ChrPosID", "Chr", "rsID", "position", "coded_all", "noncoded_all",
        "AF_coded_all", "oevar_imp", "OR", "Beta", "SE", "Z", "Pval", "N", "Neff", "HWE_pval", "imputed"} \
        else {print $2, $3, $1, $4, $5, $6, $7, $9, exp($10), $10, $11, $12, $13, $14, "399315.51", $15, $8}}' OFS="\t" \
        "$EUR_ALCP_UKB_PRE/eur_alcp_ukb_rsid_snps.txt" > "$EUR_ALCP_UKB_PRE/eur_alcp_ukb_PRE_EASYQC.txt"

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

if [ ! -f "$EUR_SMOK_23F_PRE/eur_smok_23f_PRE_EASYQC.txt" ]; then

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
	"AF_coded_all", "oevar_imp", "OR", "Beta", "SE", "Z", "Pval", "N", "Neff", "HWE_pval", "imputed"} \
        else {print $2":"$3, $2, $1, $3, $5, $6, $7, "NA", exp($8), $8, $9, $8/$9, $10, $11, "965041.67", "NA", 0}}' OFS="\t" \
        "$EUR_SMOK_23F_PRE/eur_smok_23f_impfix_rsid_snps.txt" > "$EUR_SMOK_23F_PRE/eur_smok_23f_PRE_EASYQC.txt"

fi

#2-23andMe (males)


EUR_SMOK_23M_RAW="../input/EUR/SMOK/GSCAN_SMOK_INI_MALE_EUR_23andme.txt"
EUR_SMOK_23M_PRE="../temp/EUR/SMOK"
mkdir -p "$EUR_SMOK_23M_PRE"

if [ ! -f "$EUR_SMOK_23M_PRE/eur_smok_23m_PRE_EASYQC.txt" ]; then

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
	"AF_coded_all", "oevar_imp", "OR", "Beta", "SE", "Z", "Pval", "N", "Neff", "HWE_pval", "imputed"} \
        else {print $2":"$3, $2, $1, $3, $5, $6, $7, "NA", exp($8), $8, $9, $8/$9, $10, $11, "849263.66", "NA", 0}}' OFS="\t" \
        "$EUR_SMOK_23M_PRE/eur_smok_23m_impfix_rsid_snps.txt" > "$EUR_SMOK_23M_PRE/eur_smok_23m_PRE_EASYQC.txt"

fi

#3-public

#head -n2 ../input/EUR/SMOK/GSCAN_SmkInit_2022_GWAS_SUMMARY_STATS_EUR_without_UKB.txt
#CHR	POS	RSID	EFFECT_ALLELE	OTHER_ALLELE	AF_1000G	BETA	SE	P	N
#chr10	100000235	rs11596870	T	C	0.314115	-0.000605	0.003	0.811963	357235

EUR_SMOK_PUB_RAW="../input/EUR/SMOK/GSCAN_SmkInit_2022_GWAS_SUMMARY_STATS_EUR_without_UKB.txt"
EUR_SMOK_PUB_PRE="../temp/EUR/SMOK"
mkdir -p "$EUR_SMOK_PUB_PRE"

if [ ! -f "$EUR_SMOK_PUB_PRE/eur_smok_pub_PRE_EASYQC.txt" ]; then

#prep file for liftOver (38->19)
awk 'BEGIN {OFS="\t"} NR>1 {print $1, $2-1, $2, $0}' "$EUR_SMOK_PUB_RAW" > "$EUR_SMOK_PUB_PRE/eur_smok_pub_liftOver.bed"

#liftOver to HG19
liftOver -bedPlus=3 \
        -tab \
	"$EUR_SMOK_PUB_PRE/eur_smok_pub_liftOver.bed" \
        "$CHAIN_38to19" \
        "$EUR_SMOK_PUB_PRE/eur_smok_pub_hg19.bed" \
        "$EUR_SMOK_PUB_PRE/eur_smok_pub_unmapped38.bed"

#head -n2 ../temp/EUR/SMOK/eur_smok_pub_hg19.bed
#chr10	101759991	101759992	chr10	100000235	rs11596870	T	C	0.314115	-0.000605	0.003	0.811963	357235
#chr10	101760699	101760700	chr10	100000943	rs11190359	A	G	0.0994036	-0.0034	0.004	0.385971	357235

awk 'BEGIN {OFS="\t"; print "CHR", "POS", "RSID", "EFFECT_ALLELE", "OTHER_ALLELE", "AF_1000G", "BETA", "SE", "P", "N"} \
	{sub(/^chr/, "", $1); print $1, $3, $6, $7, $8, $9, $10, $11, $12, $13}' \
	"$EUR_SMOK_PUB_PRE/eur_smok_pub_hg19.bed" > "$EUR_SMOK_PUB_PRE/eur_smok_pub_hg19.txt"

# REMOVE NON-SNPs, SNPS NOT INCLUDED IN GWAS, AND X-CHR ##
awk -F"\t" 'NR==1 {print} NR>1 && ($4 != "D" && $4 != "I") &&
	substr($3, 1, 2) == "rs" && ($9 != "" && $9 != "NA") &&
	($2 != "X") {print}' OFS="\t" \
	"$EUR_SMOK_PUB_PRE/eur_smok_pub_hg19.txt" > "$EUR_SMOK_PUB_PRE/eur_smok_pub_hg19_rsid_snps.txt"

# ATTACH AF AND ALT ALLELE FROM REF PANEL
awk -v OFS='\t' 'NR==FNR{a[$3]=$6 "\t" $7;next} FNR==1{print $0, "ALT", "AF_ALT_EUR_1000G"; next} $3 in a{print $0, a[$3];}' \
    ../temp/ref/1KG_EUR_AF_maf005.txt "$EUR_SMOK_PUB_PRE/eur_smok_pub_hg19_rsid_snps.txt" \
    > "$EUR_SMOK_PUB_PRE/eur_smok_pub_hg19_rsid_snps_af.txt"

# ALIGN AF_coded_all FROM THE ATTACHED COLUMNS ##
awk -F"\t" 'BEGIN {OFS="\t"} {
    if (NR==1) {
        print $0, "AF_coded_all"
    } else if (NR > 1 && $4 == $11) {
        print $0, $12
    } else if (NR > 1 && $5 == $11) {
        print $0, 1 - $12
    }
}' "$EUR_SMOK_PUB_PRE/eur_smok_pub_hg19_rsid_snps_af.txt" > "$EUR_SMOK_PUB_PRE/eur_smok_pub_hg19_rsid_snps_af_aligned.txt"

#NOTES
#combine chr and pos to form ChrPosID
#set oever_imp (INFO) to "NA"
#Calculate se from beta and Z: se=beta/Z
#exponentiate beta to get OR
#
#Calculate Neff sum of effective sample sizes (reported in supp)
#=Neff=4*sample_prev*(1-sample_prev)*(sample_total)=783473.64

awk -F"\t" '{if(NR == 1) {print "ChrPosID", "Chr", "rsID", "position", "coded_all", "noncoded_all",
        "AF_coded_all", "oevar_imp", "OR", "Beta", "SE", "Z", "Pval", "N", "Neff", "HWE_pval", "imputed"} \
        else {print $1":"$2, $1, $3, $2, $4, $5, $13, "NA", exp($7), $7, $8, $7/$8, $9, $10, "339759.6", "NA", 0}}' OFS="\t" \
        "$EUR_SMOK_PUB_PRE/eur_smok_pub_hg19_rsid_snps_af_aligned.txt" > "$EUR_SMOK_PUB_PRE/eur_smok_pub_PRE_EASYQC.txt"

fi



#4-UKB (no sibs)


EUR_SMOK_UKB_RAW="../input/EUR/SMOK/CLEANED.UKB_J_SMOKE_EVER"
EUR_SMOK_UKB_PRE="../temp/EUR/SMOK"
mkdir -p "$EUR_SMOK_UKB_PRE"

if [ ! -f "$EUR_SMOK_UKB_PRE/eur_smok_ukb_PRE_EASYQC.txt" ]; then

#head -n2 CLEANED.UKB_J_SMOKE_EVER
#SNP	cptid	Chr	position	EFFECT_ALLELE	OTHER_ALLELE	EAF	IMPUTED	IMPUTATION	BETA	SE	Z	PVAL	N	HWE
#rs113190524	10:62093307	10	62093307	C	T	0.719095	1	0.994363	0.000752903	0.00120016	0.627335521930409	0.53	403349	NA

# REMOVE NON-SNPs, SNPS NOT INCLUDED IN GWAS, AND X-CHR ##
awk -F"\t" 'NR==1 {print} NR>1 && ($5 != "D" && $5 != "I") &&
        substr($1, 1, 2) == "rs" && ($13 != "" && $13 != "NA") &&
        ($3 != "X") {print}' OFS="\t" \
        "$EUR_SMOK_UKB_RAW" > "$EUR_SMOK_UKB_PRE/eur_smok_ukb_rsid_snps.txt"

#NOTES
#Calculate Neff (prev from GSCAN2 article based on full UKB):
#=N_CAS: 403349*.45=181507
#=N_CON: 403349-181507=221842
#=SAMP_PREV: 181507/403349=.45
#=Neff=4*.45*(1-.45)*(403349)=399315.51

awk -F"\t" '{if(NR == 1) {print "ChrPosID", "Chr", "rsID", "position", "coded_all", "noncoded_all",
        "AF_coded_all", "oevar_imp", "OR", "Beta", "SE", "Z", "Pval", "N", "Neff", "HWE_pval", "imputed"} \
        else {print $2, $3, $1, $4, $5, $6, $7, $9, exp($10), $10, $11, $12, $13, $14, "399315.51", $15, $8}}' OFS="\t" \
        "$EUR_SMOK_UKB_PRE/eur_smok_ukb_rsid_snps.txt" > "$EUR_SMOK_UKB_PRE/eur_smok_ukb_PRE_EASYQC.txt"

fi

