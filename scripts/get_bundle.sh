#!/bin/env bash

#############################
# get Broad hg38 BUNDLE data
#############################

LIST=$(cat <<'END_HEREDOC'
1000G_omni2.5.hg38.vcf.gz
1000G_omni2.5.hg38.vcf.gz.tbi
1000G_phase1.snps.high_confidence.hg38.vcf.gz
1000G_phase1.snps.high_confidence.hg38.vcf.gz.tbi
Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz
Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz.tbi
Homo_sapiens_assembly38.dict
Homo_sapiens_assembly38.fasta.64.alt
Homo_sapiens_assembly38.fasta.fai
Homo_sapiens_assembly38.fasta.gz
Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi
dbsnp_138.hg38.vcf.gz
dbsnp_138.hg38.vcf.gz.tbi
dbsnp_144.hg38.vcf.gz
dbsnp_144.hg38.vcf.gz.tbi
dbsnp_146.hg38.vcf.gz
dbsnp_146.hg38.vcf.gz.tbi
hapmap_3.3.hg38.vcf.gz
hapmap_3.3.hg38.vcf.gz.tbi
hapmap_3.3_grch38_pop_stratified_af.vcf.gz
hapmap_3.3_grch38_pop_stratified_af.vcf.gz.tbi
wgs_calling_regions.hg38.interval_list
END_HEREDOC
)

# create folder and get data
mkdir -p reference

# Broad bundle repo for hg38
bundleurl=ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38

for file in $LIST; do
echo "# downloading $url"
wget -P reference -np --ftp-user=gsapubftp-anonymous ${bundleurl}/$file"
done
