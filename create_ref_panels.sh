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

#========================================================================================#
# OPERATIONS
#========================================================================================#

#set path variables to ref panel files (downloaded from Plink2 website)
PGEN_RAW="../input/ref/1kg_pgen/all_phase3"
PSAM_RAW="../input/ref/1kg_pgen/all_phase3.psam"
REL_IDS="../input/ref/1kg_pgen/deg2_phase3.king.cutoff.out.id"
TEMP_REF="../temp/ref"

mkdir -p $TEMP_REF

#set MAF thesholds
AFR_MAF="0.005"
AFR_MAF_NAME="maf005"

EUR_MAF="$AFR_MAF"
EUR_MAF_NAME="$AFR_MAF_NAME"

#=AFR REF PANEL==========================================================================#

if [ ! -f "$TEMP_REF/1KG_AFR_AF_${AFR_MAF_NAME}_rsid.txt" ]; then

#find all subject IDs that were not reported as 2nd degree relatives
#select only individuals belonging to the AFR superpopulation
grep -vFf "$REL_IDS" "$PSAM_RAW" | awk '$5 == "AFR" {print $1}' > "$TEMP_REF/all_phase3_king2d_afr.txt"

#keep individuals identified above
#keep only biallelic SNPs (no indels)
#restrict to SNPs with greater than 0.005 MAF
plink2 	--pfile "$PGEN_RAW" \
	--keep "$TEMP_REF/all_phase3_king2d_afr.txt" \
	--snps-only just-acgt \
	--maf "$AFR_MAF" \
	--make-pgen \
	--out "$TEMP_REF/all_phase3_king2d_afr_snps_${AFR_MAF_NAME}"

#get frequency file
plink2	--pfile "$TEMP_REF/all_phase3_king2d_afr_snps_${AFR_MAF_NAME}" \
	--freq \
	--out "$TEMP_REF/all_phase3_king2d_afr_snps_${AFR_MAF_NAME}"


#extract SNP info
awk 'BEGIN {OFS="\t"} !/^##/{print $2, $3}' "$TEMP_REF/all_phase3_king2d_afr_snps_${AFR_MAF_NAME}.pvar" > "$TEMP_REF/positions_afr.txt"

#check that rows were kept consistent
#cut -f2 ../temp/ref/all_phase3_king2d_afr_snps_maf01.afreq | head
#ID
#rs558604819
#rs575272151
#rs544419019
#rs561109771
#rs62635286
#rs200579949
#rs531730856
#rs571093408
#rs546169444

#cut -f2 ../temp/ref/positions.txt | head
#ID
#rs558604819
#rs575272151
#rs544419019
#rs561109771
#rs62635286
#rs200579949
#rs531730856
#rs571093408
#rs546169444

#Paste files together to add position data to freq file
paste "$TEMP_REF/all_phase3_king2d_afr_snps_${AFR_MAF_NAME}.afreq" "$TEMP_REF/positions_afr.txt" > "$TEMP_REF/all_phase3_king2d_afr_snps_${AFR_MAF_NAME}_freq.txt"

#double check that rows were correctly aligned

#awk '{ if ($2 != $8) { print } }' ../temp/ref/all_phase3_king2d_afr_snps_maf01_freq.txt
#[no rows returned]

#final formatting
awk 'BEGIN {OFS="\t"} \
	NR==1 {print "ChrPosID", "Chr", "RSID", "Pos", "REF", "ALT", "AF_ALT_AFR_1000G"} \
	NR>1 {print $1":"$7, $1, $2, $7, $3, $4, $5}' \
	"$TEMP_REF/all_phase3_king2d_afr_snps_${AFR_MAF_NAME}_freq.txt" > "$TEMP_REF/1KG_AFR_AF_${AFR_MAF_NAME}.txt"


#wc -l ../temp/ref/1KG_AFR_AF_MAF01.txt
#15482394 ../temp/ref/1KG_AFR_AF_MAF01.txt

#head -n2 ../temp/ref/1KG_AFR_AF_MAF01.txt
#ChrPosID       Chr     RSID    Pos     REF     ALT     AF_ALT_AFR_1000G
#1:10642        1       rs558604819     10642   G       A       0.0130568

#split reference file into RSID and ALLELE files
  # RSID file (RSID, CHROM, POS)
  awk -F'\t' -v OFS='\t' '{print $3, $2, $4}' "$TEMP_REF/1KG_AFR_AF_${AFR_MAF_NAME}.txt" > "$TEMP_REF/1KG_AFR_AF_${AFR_MAF_NAME}_rsid.txt"

  # ALLELE file (ChrPosID, REF, ALT, AF_ALT_AFR_1000G)
  awk -F'\t' -v OFS='\t' '{print $1, $5, $6, $7}' "$TEMP_REF/1KG_AFR_AF_${AFR_MAF_NAME}.txt" > "$TEMP_REF/1KG_AFR_AF_${AFR_MAF_NAME}_allele.txt"

fi






#=EUR REF PANEL==========================================================================#

if [ ! -f "$TEMP_REF/1KG_EUR_AF_${EUR_MAF_NAME}_rsid.txt" ]; then

#find all subject IDs that were not reported as 2nd degree relatives
#select only individuals belonging to the EUR superpopulation
grep -vFf "$REL_IDS" "$PSAM_RAW" | awk '$5 == "EUR" {print $1}' > "$TEMP_REF/all_phase3_king2d_eur.txt"

#keep individuals identified above
#keep only biallelic SNPs (no indels)
#restrict to SNPs with greater than 0.005 MAF
plink2  --pfile "$PGEN_RAW" \
        --keep "$TEMP_REF/all_phase3_king2d_eur.txt" \
        --snps-only just-acgt \
        --maf "$EUR_MAF" \
        --make-pgen \
        --out "$TEMP_REF/all_phase3_king2d_eur_snps_${EUR_MAF_NAME}"

#get frequency file
plink2  --pfile "$TEMP_REF/all_phase3_king2d_eur_snps_${EUR_MAF_NAME}" \
        --freq \
        --out "$TEMP_REF/all_phase3_king2d_eur_snps_${EUR_MAF_NAME}"


#extract SNP info
awk 'BEGIN {OFS="\t"} !/^##/{print $2, $3}' "$TEMP_REF/all_phase3_king2d_eur_snps_${EUR_MAF_NAME}.pvar" > "$TEMP_REF/positions_eur.txt"


#Paste files together to add position data to freq file
paste "$TEMP_REF/all_phase3_king2d_eur_snps_${EUR_MAF_NAME}.afreq" "$TEMP_REF/positions_eur.txt" > "$TEMP_REF/all_phase3_king2d_eur_snps_${EUR_MAF_NAME}_freq.txt"

#double check that rows were correctly aligned

#awk '{ if ($2 != $8) { print } }' "$TEMP_REF/all_phase3_king2d_eur_snps_${EUR_MAF_NAME}_freq.txt"
#[no rows returned]

#final formatting
awk 'BEGIN {OFS="\t"} \
        NR==1 {print "ChrPosID", "Chr", "RSID", "Pos", "REF", "ALT", "AF_ALT_EUR_1000G"} \
        NR>1 {print $1":"$7, $1, $2, $7, $3, $4, $5}' \
        "$TEMP_REF/all_phase3_king2d_eur_snps_${EUR_MAF_NAME}_freq.txt" > "$TEMP_REF/1KG_EUR_AF_${EUR_MAF_NAME}.txt"


#wc -l

#head -n2

#split reference file into RSID and ALLELE files
  # RSID file (RSID, CHROM, POS)
  awk -F'\t' -v OFS='\t' '{print $3, $2, $4}' "$TEMP_REF/1KG_EUR_AF_${EUR_MAF_NAME}.txt" > "$TEMP_REF/1KG_EUR_AF_${EUR_MAF_NAME}_rsid.txt"

  # ALLELE file (ChrPosID, REF, ALT, AF_ALT_AFR_1000G)
  awk -F'\t' -v OFS='\t' '{print $1, $5, $6, $7}' "$TEMP_REF/1KG_EUR_AF_${EUR_MAF_NAME}.txt" > "$TEMP_REF/1KG_EUR_AF_${EUR_MAF_NAME}_allele.txt"

fi

