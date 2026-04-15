#!/bin/bash

set -euo pipefail

#  Filter parameters to vcftools

MAF=0.025
MISS=0.875
QUAL=50
MIN_DEPTH=10
MAX_DEPTH=50

vcftools --gzvcf NC_058189.1.vcf.gz --remove-indels --maf $MAF --max-missing $MISS --minQ $QUAL --min-meanDP $MIN_DEPTH --max-meanDP $MAX_DEPTH --minDP $MIN_DEPTH --maxDP $MAX_DEPTH --recode --stdout | bgzip -c > FILTERED_CHR_VCFS/NC_058189.1.filtered.vcf.gz
tabix -p vcf FILTERED_CHR_VCFS/NC_058189.1.filtered.vcf.gz

vcftools --gzvcf NC_058190.1.vcf.gz --remove-indels --maf $MAF --max-missing $MISS --minQ $QUAL --min-meanDP $MIN_DEPTH --max-meanDP $MAX_DEPTH --minDP $MIN_DEPTH --maxDP $MAX_DEPTH --recode --stdout | bgzip -c > FILTERED_CHR_VCFS/NC_058190.1.filtered.vcf.gz
tabix -p vcf FILTERED_CHR_VCFS/NC_058190.1.filtered.vcf.gz

vcftools --gzvcf NC_058191.1.vcf.gz --remove-indels --maf $MAF --max-missing $MISS --minQ $QUAL --min-meanDP $MIN_DEPTH --max-meanDP $MAX_DEPTH --minDP $MIN_DEPTH --maxDP $MAX_DEPTH --recode --stdout | bgzip -c > FILTERED_CHR_VCFS/NC_058191.1.filtered.vcf.gz
tabix -p vcf FILTERED_CHR_VCFS/NC_058191.1.filtered.vcf.gz

vcftools --gzvcf NC_058192.1.vcf.gz --remove-indels --maf $MAF --max-missing $MISS --minQ $QUAL --min-meanDP $MIN_DEPTH --max-meanDP $MAX_DEPTH --minDP $MIN_DEPTH --maxDP $MAX_DEPTH --recode --stdout | bgzip -c > FILTERED_CHR_VCFS/NC_058192.1.filtered.vcf.gz
tabix -p vcf FILTERED_CHR_VCFS/NC_058192.1.filtered.vcf.gz

vcftools --gzvcf NC_058193.1.vcf.gz --remove-indels --maf $MAF --max-missing $MISS --minQ $QUAL --min-meanDP $MIN_DEPTH --max-meanDP $MAX_DEPTH --minDP $MIN_DEPTH --maxDP $MAX_DEPTH --recode --stdout | bgzip -c > FILTERED_CHR_VCFS/NC_058193.1.filtered.vcf.gz
tabix -p vcf FILTERED_CHR_VCFS/NC_058193.1.filtered.vcf.gz

vcftools --gzvcf NC_058194.1.vcf.gz --remove-indels --maf $MAF --max-missing $MISS --minQ $QUAL --min-meanDP $MIN_DEPTH --max-meanDP $MAX_DEPTH --minDP $MIN_DEPTH --maxDP $MAX_DEPTH --recode --stdout | bgzip -c > FILTERED_CHR_VCFS/NC_058194.1.filtered.vcf.gz
tabix -p vcf FILTERED_CHR_VCFS/NC_058194.1.filtered.vcf.gz

vcftools --gzvcf NC_058195.1.vcf.gz --remove-indels --maf $MAF --max-missing $MISS --minQ $QUAL --min-meanDP $MIN_DEPTH --max-meanDP $MAX_DEPTH --minDP $MIN_DEPTH --maxDP $MAX_DEPTH --recode --stdout | bgzip -c > FILTERED_CHR_VCFS/NC_058195.1.filtered.vcf.gz
tabix -p vcf FILTERED_CHR_VCFS/NC_058195.1.filtered.vcf.gz

vcftools --gzvcf NC_058196.1.vcf.gz --remove-indels --maf $MAF --max-missing $MISS --minQ $QUAL --min-meanDP $MIN_DEPTH --max-meanDP $MAX_DEPTH --minDP $MIN_DEPTH --maxDP $MAX_DEPTH --recode --stdout | bgzip -c > FILTERED_CHR_VCFS/NC_058196.1.filtered.vcf.gz

tabix -p vcf FILTERED_CHR_VCFS/NC_058196.1.filtered.vcf.gz

vcftools --gzvcf NC_058197.1.vcf.gz --remove-indels --maf $MAF --max-missing $MISS --minQ $QUAL --min-meanDP $MIN_DEPTH --max-meanDP $MAX_DEPTH --minDP $MIN_DEPTH --maxDP $MAX_DEPTH --recode --stdout | bgzip -c > FILTERED_CHR_VCFS/NC_058197.1.filtered.vcf.gz

tabix -p vcf FILTERED_CHR_VCFS/NC_058197.1.filtered.vcf.gz

vcftools --gzvcf NC_058198.1.vcf.gz --remove-indels --maf $MAF --max-missing $MISS --minQ $QUAL --min-meanDP $MIN_DEPTH --max-meanDP $MAX_DEPTH --minDP $MIN_DEPTH --maxDP $MAX_DEPTH --recode --stdout | bgzip -c > FILTERED_CHR_VCFS/NC_058198.1.filtered.vcf.gz

tabix -p vcf FILTERED_CHR_VCFS/NC_058198.1.filtered.vcf.gz