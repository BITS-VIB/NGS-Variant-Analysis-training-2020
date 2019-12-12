#!/bin/env bash
# author:Stephane Plaisance (VIB-NC), 2019-12-12

## Where are we?
workdir=/data/NC_projects/BITS_variants_2020

## remap reads from 300x chr22 extarct
cd ${workdir}

# create tmp folder
mkdir -p tmpfiles

# create date-tagged logs 
mkdir -p logs
log=logs/runlog.txt
cat /dev/null > ${log}

# record all script actions
set -x
#exec > ${log} 2>&1
exec &> >(tee -i ${log})

# multithreading, adapt to your cpu count
bwathr=84
samtoolsthr=24

# how much RAM can java use? - set these to less than your available RAM!
javaopts="-Xms24g -Xmx24g"

# my samtools is here
samtools=$BIOTOOLS/samtools/bin/samtools

#############################################
# BWA MEM mapping
#############################################

## create a bwa index when absent
mkdir -p bwaidx

reference_fa=reference/Homo_sapiens.UCSC.hg38.fa
bwaidx="bwaidx/Homo_sapiens.UCSC.hg38.fa"

if [ ! -f bwaidx/index_created ]; then
echo "# creating BWA index"
bwa index ${reference_fa} -p ${bwaidx} \
  && touch bwaidx/index_created
else
echo "# BWA index already exists"
fi

## map reads to reference
reads_1="reads/HG001.GRCh38_chr22_1.fq.gz"
reads_2="reads/HG001.GRCh38_chr22_2.fq.gz"
thr=84
outpfx="HG001_chr22"
samplename="NA12878"

outfolder=bwa_mappings
mkdir -p ${workdir}/${outfolder}

# map using BWA mem
cmd="bwa mem -t ${bwathr} \
	-M \
	-R '@RG\tID:HG001\tLB:NA12878_giab\tPU:unknown-0.0\tPL:Illumina\tSM:NA12878' \
	${bwaidx} \
	${reads_1} \
	${reads_2} | ${samtools} view -b - -o ${outfolder}/${outpfx}_rawmappings.bam"

echo "# ${cmd}"
eval ${cmd}

# get flagstats
${samtools} flagstat -@ ${samtoolsthr} \
	${outfolder}/${outpfx}_rawmappings.bam \
	> ${outfolder}/${outpfx}_rawmappings_flagstats.txt


#############################################
# PICARD Cleanup & MarkDuplicates
#############################################

# more records in RAM speeds up when enough RAM is present
recinram=10000000

# sort by queryname
java ${javaopts} -jar $PICARD/picard.jar \
	SortSam \
	I=${outfolder}/${outpfx}_rawmappings.bam \
	O=${outfolder}/${outpfx}_rawmappings_qrysrt.bam \
	SO=queryname \
	MAX_RECORDS_IN_RAM=${recinram} \
	TMP_DIR=tmpfiles/

# mark duplicates
pixdist=100
java ${javaopts} -jar $PICARD/picard.jar \
	MarkDuplicates \
	I=${outfolder}/${outpfx}_rawmappings_qrysrt.bam \
	O=${outfolder}/${outpfx}_mrkdup.bam \
	M=${outfolder}/${outpfx}_MarkDuplicates.txt \
	ASO=queryname \
	REMOVE_DUPLICATES=false \
	OPTICAL_DUPLICATE_PIXEL_DISTANCE=${pixdist} \
	MAX_RECORDS_IN_RAM=${recinram} \
	TMP_DIR=tmpfiles/

# sort by coordinate and index
java ${javaopts} -jar $PICARD/picard.jar \
	SortSam \
	I=${outfolder}/${outpfx}_mrkdup.bam \
	O=${outfolder}/${outpfx}_mrkdup_srt.bam \
	SO=coordinate \
	CREATE_INDEX=true \
	MAX_RECORDS_IN_RAM=${recinram} \
	TMP_DIR=tmpfiles/

# fix tags
java ${javaopts} -jar $PICARD/picard.jar \
	SetNmMdAndUqTags \
	I=${outfolder}/${outpfx}_mrkdup_srt.bam \
	O=${outfolder}/${outpfx}_mrkdup_srt-tags.bam \
	R=${reference_fa} \
	CREATE_INDEX=true \
	MAX_RECORDS_IN_RAM=${recinram} \
	TMP_DIR=tmpfiles/

# validate final BAM content
java ${javaopts} -jar $PICARD/picard.jar \
	ValidateSamFile \
	I=${outfolder}/${outpfx}_mrkdup_srt-tags.bam \
	O=${outfolder}/${outpfx}_mrkdup_srt-tags.bam_ValidateSamFile.txt \
	R=${reference_fa} \
	M=SUMMARY \
	MO=100 \
	IGNORE_WARNINGS=FALSE \
	MAX_RECORDS_IN_RAM=${recinram} \
	TMP_DIR=tmpfiles/

# extract chr22 reads for calling only chr22 variants
chr="chr22"
${samtools} view -@ ${samtoolsthr} -h -b \
	${outfolder}/${outpfx}_mrkdup_srt-tags.bam ${chr} \
	-o ${outfolder}/${outpfx}_mrkdup_srt_22only.bam && \
	${samtools} index ${outfolder}/${outpfx}_mrkdup_srt_22only.bam


#############################################
# GATK4 BAM PROCESSING
#############################################

outfolder=gatk_preprocessing
mkdir -p ${outfolder}

bamfile=bwa_mappings/${outpfx}_mrkdup_srt_22only.bam
recalbamfile=${outpfx}_mrkdup_srt_recal.bam

# base quality score recalibration
knownsites="reference/hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf.gz"
knownindels="reference/hg38_v0_Homo_sapiens_assembly38.known_indels.vcf.gz"

# compute table before
java ${javaopts} -jar $GATK/gatk.jar \
	BaseRecalibrator \
	-I ${bamfile} \
	-R ${reference_fa} \
	--known-sites ${knownsites} \
	--known-sites ${knownindels} \
	-O ${outfolder}/recal_data.table \
	--tmp-dir tmpfiles/

# apply recalibration table
java ${javaopts} -jar $GATK/gatk.jar \
	ApplyBQSR \
	-I ${bamfile} \
	-R ${reference_fa} \
	-bqsr ${outfolder}/recal_data.table \
	-O ${outfolder}/${recalbamfile} \
	--tmp-dir tmpfiles/

# compute table after
java ${javaopts} -jar $GATK/gatk.jar \
	BaseRecalibrator \
	-I ${outfolder}/${recalbamfile} \
	-R ${reference_fa} \
	--known-sites ${knownsites} \
	--known-sites ${knownindels} \
	-O ${outfolder}/recal_data_after.table \
	--tmp-dir tmpfiles/

# create plots from both tables
java ${javaopts} -jar $GATK/gatk.jar \
	AnalyzeCovariates \
	-before ${outfolder}/recal_data.table \
	-after ${outfolder}/recal_data_after.table \
	-plots ${outfolder}/BQSR_report.pdf \
	-csv ${outfolder}/BQSR-report.csv \
	--tmp-dir tmpfiles/

# Picard CollectMultipleMetrics on final BAM
java ${javaopts} -jar $PICARD/picard.jar \
	CollectMultipleMetrics \
	I=${outfolder}/${recalbamfile} \
	O=gatk_multiple_metrics/multiple_metrics \
	R=${reference_fa} \
	TMP_DIR=tmpfiles/


#############################################
# GATK4 VARIANT CALLING (chr22 only)
#############################################

outfolder=gatk_variantcalling
mkdir -p ${outfolder}

# call short variants to gvcf format and save supporting reads
java ${javaopts} -jar $GATK/gatk.jar \
	HaplotypeCaller  \
	--input gatk_preprocessing/${recalbamfile} \
	--output ${outfolder}/${samplename}.g.vcf.gz \
	--reference ${reference_fa} \
	--emit-ref-confidence GVCF \
	--sample-ploidy 2 \
	--intervals chr22 \
	--bam-output ${outfolder}/${samplename}_HC_aligned_reads.bam \
	--tmp-dir tmpfiles/


#############################################
# GATK4 VARIANT RECALIBRATION
#############################################

outfolder=gatk_variantrecalibration
mkdir -p ${outfolder}

# recalibration sources

## variants-sets
# True sites training resource: HapMap
truetraining15=reference/hg38_v0_hapmap_3.3.hg38.vcf.gz

# True sites training resource: Omni
truetraining12=reference/hg38_v0_1000G_omni2.5.hg38.vcf.gz

# Non-true sites training resource: 1000G
nontruetraining10=reference/hg38_v0_Mills_and_1000G_gold_standard.indels.hg38.vcf.gz

# Known sites resource, not used in training: dbSNP
knowntraining2=reference/hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf.gz

# indels True sites training resource: Mills
truetrainingindel12=reference/hg38_v0_Mills_and_1000G_gold_standard.indels.hg38.vcf.gz

# indels Axiom
axiom10=reference/Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz

# convert to VCF
java ${javaopts} -jar $GATK/gatk.jar \
	GenotypeGVCFs \
	--reference ${reference_fa} \
	--variant gatk_variantcalling/${samplename}.g.vcf.gz \
	--output ${outfolder}/${samplename}.vcf.gz \
	--tmp-dir tmpfiles/

# mark ExcessHet
threshold=54.69
java ${javaopts} -jar $GATK/gatk.jar \
	VariantFiltration \
	--variant ${outfolder}/${samplename}.vcf.gz \
	--filter-expression "ExcessHet > ${threshold}" \
	--filter-name "ExcessHet" \
	-O ${outfolder}/${samplename}_excesshet_filtered.vcf.gz \
	--tmp-dir tmpfiles/

# extract a 6-column version of the data for recalibration
java ${javaopts} -jar $GATK/gatk.jar \
	MakeSitesOnlyVcf \
	-I ${outfolder}/${samplename}_excesshet_filtered.vcf.gz \
	-O ${outfolder}/${samplename}_excesshet_sitesonly.vcf.gz

# Build the SNP recalibration model
java ${javaopts} -jar $GATK/gatk.jar \
	VariantRecalibrator \
	-R ${reference_fa} \
	-V ${outfolder}/${samplename}_excesshet_sitesonly.vcf.gz \
	-O ${outfolder}/${samplename}_recalibrate_SNP.recal.vcf.gz \
	--resource:hapmap,known=false,training=true,truth=true,prior=15.0 ${truetraining15} \
	--resource:omni,known=false,training=true,truth=true,prior=12.0 ${truetraining12} \
	--resource:1000G,known=false,training=true,truth=false,prior=10.0 ${nontruetraining10} \
	--resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ${knowntraining2} \
	--use-annotation DP \
	--use-annotation QD \
	--use-annotation FS \
	--use-annotation SOR \
	--use-annotation MQ \
	--use-annotation MQRankSum \
	--use-annotation ReadPosRankSum \
	--mode SNP \
	--tranches-file ${outfolder}/${samplename}_recalibrate_snp.tranches \
	--tranche 100.0 \
	--tranche 99.9 \
	--tranche 99.0 \
	--tranche 90.0 \
	--rscript-file ${outfolder}/${samplename}_recalibrate_snp_plots.R \
	--tmp-dir tmpfiles/

# Build the Indel recalibration model
java ${javaopts} -jar $GATK/gatk.jar \
	VariantRecalibrator \
	-R ${reference_fa} \
	-V ${outfolder}/${samplename}_excesshet_sitesonly.vcf.gz \
	-O ${outfolder}/${samplename}_recalibrate_indels.recal.vcf.gz \
	--resource:mills,known=false,training=true,truth=true,prior=12 ${truetrainingindel12} \
	--resource:axiomPoly,known=false,training=true,truth=false,prior=10 ${axiom10} \
	--resource:dbsnp,known=true,training=false,truth=false,prior=2 ${knowntraining2} \
	--trust-all-polymorphic \
	--use-annotation QD \
	--use-annotation DP \
	--use-annotation FS \
	--use-annotation SOR \
	--use-annotation MQRankSum \
	--use-annotation ReadPosRankSum \
	--mode INDEL \
	--max-gaussians 4 \
	--tranches-file ${outfolder}/${samplename}_recalibrate_indels.tranches \
	--tranche 100.0 \
	--tranche 99.95 \
	--tranche 99.9 \
	--tranche 99.5 \
	--tranche 99.0 \
	--tranche 97.0 \
	--tranche 96.0 \
	--tranche 95.0 \
	--tranche 94.0 \
	--tranche 93.5 \
	--tranche 93.0 \
	--tranche 92.0 \
	--tranche 91.0 \
	--tranche 90.0 \
	--rscript-file ${outfolder}/${samplename}_recalibrate_indels_plots.R \
	--tmp-dir tmpfiles/

# Apply the desired level of recalibration to the SNPs in the call set
java ${javaopts} -jar $GATK/gatk.jar \
	ApplyVQSR \
	-R ${reference_fa} \
	-V ${outfolder}/${samplename}_excesshet_filtered.vcf.gz \
	--mode SNP \
	--recal-file ${outfolder}/${samplename}_recalibrate_SNP.recal.vcf.gz \
	--tranches-file ${outfolder}/${samplename}_recalibrate_snp.tranches \
	--truth-sensitivity-filter-level 99.0 \
	-O ${outfolder}/${samplename}_recalibrated_snps_raw_indels.vcf.gz \
	--tmp-dir tmpfiles/

# Apply the desired level of recalibration to the Indels in the call set
java ${javaopts} -jar $GATK/gatk.jar \
	ApplyVQSR \
	-R ${reference_fa} \
	-V ${outfolder}/${samplename}_recalibrated_snps_raw_indels.vcf.gz \
	--recal-file ${outfolder}/${samplename}_recalibrate_indels.recal.vcf.gz \
	--tranches-file ${outfolder}/${samplename}_recalibrate_indels.tranches \
	--truth-sensitivity-filter-level 99.7 \
	--create-output-variant-index true \
	-mode INDEL \
	-O ${outfolder}/${samplename}_VQSR.vcf.gz \
	--tmp-dir tmpfiles/


#############################################
# SnpEff on VQSR
#############################################

outfolder="snpeff"
mkdir -p ${outfolder}

build="GRCh38.86"

java ${javaopts} -jar $SNPEFF/snpEff.jar \
	-htmlStats ${outfolder}/gatk_snpEff_summary.html \
	${build} \
	gatk_variantrecalibration/${samplename}_VQSR.vcf.gz | \
	bgzip -c > ${outfolder}/${samplename}_VQSR_snpeff.vcf.gz && \
	tabix -p vcf ${outfolder}/${samplename}_VQSR_snpeff.vcf.gz


#############################################
# GATK4 VARIANT HARD FILTERING
#############################################

# instructions and filters from:
# https://gatkforums.broadinstitute.org/gatk/discussion/23216/how-to-filter-variants-either-with-vqsr-or-by-hard-filtering
# Genepattern does the following: QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || SOR > 3.0 || QUAL < 30

outfolder=gatk_varianthardfiltering
mkdir -p ${outfolder}

# copy the GenotypeGVCFs vcf output with marked excesshet output from above
# create local copy of previous file
cp gatk_variantrecalibration/${samplename}_excesshet_filtered.vcf.gz* ${outfolder}/


###########################################
# 2) hard-Filter SNPs on multiple metrics
###########################################

# produces a VCF with records with SNP-type variants only.
java ${javaopts} -jar $GATK/gatk.jar \
    SelectVariants \
    -V ${outfolder}/${samplename}_excesshet_filtered.vcf.gz \
    --select-type-to-include SNP \
    -O ${outfolder}/${samplename}_snp.vcf.gz

# This produces a VCF with the same variant records now annotated with filter status. Specifically, if a record passes all the filters, it receives a PASS label in the FILTER column.
# A record that fails a filter #receives the filter name in the FILTER column, e.g. SOR3.
# If a record fails multiple filters, then each failing filter name appears in the FILTER column separated by semi-colons ; e.g. "MQRankSum-12.5;ReadPosRankSum-8".

java ${javaopts} -jar $GATK/gatk.jar \
    VariantFiltration \
    -V ${outfolder}/${samplename}_snp.vcf.gz \
    --filter "QD < 2.0" --filter-name "QD2" \
    --filter "QUAL < 30.0" --filter-name "QUAL30" \
    --filter "SOR > 3.0" --filter-name "SOR3" \
    --filter "FS > 60.0" --filter-name "FS60" \
    --filter "MQ < 40.0" --filter-name "MQ40" \
    --filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
    --filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
    -O ${outfolder}/${samplename}_snp_filtered.vcf.gz

###########################################
# 3) hard-Filter INDELs on multiple metrics
###########################################

# produces a VCF with records with INDEL-type variants only.
java ${javaopts} -jar $GATK/gatk.jar \
    SelectVariants \
    -V ${outfolder}/${samplename}_excesshet_filtered.vcf.gz \
    --select-type-to-include INDEL \
    -O ${outfolder}/${samplename}_indels.vcf.gz

java ${javaopts} -jar $GATK/gatk.jar \
    VariantFiltration \
    -V ${outfolder}/${samplename}_indels.vcf.gz \
    --filter "QD < 2.0" --filter-name "QD2" \
    --filter "QUAL < 30.0" --filter-name "QUAL30" \
    --filter "FS > 200.0" --filter-name "FS200" \
    --filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" \
    -O ${outfolder}/${samplename}_indels_filtered.vcf.gz

###########################################
## 4) merge SNP and Indel filtered calls    
###########################################

# combine filtered into one file using Picard
java ${javaopts} -jar $GATK/gatk.jar \
    MergeVcfs \
    -I ${outfolder}/${samplename}_snp_filtered.vcf.gz \
    -I ${outfolder}/${samplename}_indels_filtered.vcf.gz \
    -R ${reference_fa} \
    -O ${outfolder}/${samplename}_snp_indel_filtered.vcf.gz \
	CREATE_INDEX=true \
	TMP_DIR=tmpfiles/

#############################################
# SnpEff on HardFiltering
#############################################

outfolder="snpeff"
mkdir -p ${outfolder}

build="GRCh38.86"

java ${javaopts} -jar $SNPEFF/snpEff.jar \
	-htmlStats ${outfolder}/gatk_hardfiltering_snpEff_summary.html \
	${build} \
	gatk_varianthardfiltering/${samplename}_snp_indel_filtered.vcf.gz | \
	bgzip -c > ${outfolder}/${samplename}_hardfiltering_snpeff.vcf.gz && \
	tabix -p vcf ${outfolder}/${samplename}_hardfiltering_snpeff.vcf.gz


#############################################
# Compare calls to Golden
#############################################

gold=reference/NA12878_HG001-chr22_gold.vcf.gz

raw=gatk_variantrecalibration/NA12878.vcf.gz
vqsr=snpeff/NA12878_VQSR_snpeff.vcf.gz
hard=snpeff/NA12878_hardfiltering_snpeff.vcf.gz

# compare 3
vcf-compare -r chr22 \
    ${vqsr} \
    ${hard} \
    ${gold} | grep "^VN" > snpeff/3venn_counts.txt

#A.C VN	29	reference/NA12878_HG001-chr22_gold.vcf.gz (0.1%)	snpeff/NA12878_VQSR_snpeff.vcf.gz (0.0%)
#..C VN	51	reference/NA12878_HG001-chr22_gold.vcf.gz (0.1%)
#A.. VN	128	snpeff/NA12878_VQSR_snpeff.vcf.gz (0.1%)
#ABC VN	42110	reference/NA12878_HG001-chr22_gold.vcf.gz (99.8%)	snpeff/NA12878_VQSR_snpeff.vcf.gz (46.1%)	snpeff/NA12878_hardfiltering_snpeff.vcf.gz (46.2%)
#AB. VN	49095	snpeff/NA12878_VQSR_snpeff.vcf.gz (53.7%)	snpeff/NA12878_hardfiltering_snpeff.vcf.gz (53.8%)

# create venn plot
# use R to plot
# https://wiki.bits.vib.be/index.php/NGS_Exercise.6#Create_a_VENN_diagram_for_the_4_mappings

3DVenn.R -a 128 -A chr22_VQSR \
    -b 0 -B chr22_HardF \
    -c 51 -C chr22_GoldS \
    -d 49095 -e 29 -f 0 -i 42110 -t "NA12878 variant recall" -x 2 -o snpeff/recall3 -u 1
    
3DVenn.R -a 128 -A chr22_VQSR \
    -b 0 -B chr22_HardF \
    -c 51 -C chr22_GoldS \
    -d 49095 -e 29 -f 0 -i 42110 -t "NA12878 variant recall" -x 2 -o snpeff/recall3pc -u 1 -P 1


# compare 4 including raw calls
vcf-compare -r chr22 \
	${raw} \
    ${vqsr} \
    ${hard} \
    ${gold} | grep "^VN" > snpeff/4venn_counts.txt

#AB.D VN	29	gatk_variantrecalibration/NA12878.vcf.gz (0.0%)	reference/NA12878_HG001-chr22_gold.vcf.gz (0.1%)	snpeff/NA12878_VQSR_snpeff.vcf.gz (0.0%)
#...D VN	51	reference/NA12878_HG001-chr22_gold.vcf.gz (0.1%)
#AB.. VN	128	gatk_variantrecalibration/NA12878.vcf.gz (0.1%)	snpeff/NA12878_VQSR_snpeff.vcf.gz (0.1%)
#ABCD VN	42110	gatk_variantrecalibration/NA12878.vcf.gz (46.1%)	reference/NA12878_HG001-chr22_gold.vcf.gz (99.8%)	snpeff/NA12878_VQSR_snpeff.vcf.gz (46.1%)	snpeff/NA12878_hardfiltering_snpeff.vcf.gz (46.2%)
#ABC. VN	49095	gatk_variantrecalibration/NA12878.vcf.gz (53.7%)	snpeff/NA12878_VQSR_snpeff.vcf.gz (53.7%)	snpeff/NA12878_hardfiltering_snpeff.vcf.gz (53.8%)

4DVenn.R -a 128 -A chr22_raw \
    -b 0 -B chr22_VQSR \
    -c 0 -C chr22_HardF \
	-d 51 -D chr22_GoldS \
	-e 128 \
	-f 0 \
	-G 0 \
	-j 0 \
	-k 0 \
	-l 0 \
	-m 49095 \
	-n 29 \
	-p 0 \
	-q 0 \
	-i 42110 \
	-t "NA12878 variant recall" -x 2 -o snpeff/recall4 -u 1

4DVenn.R -a 128 -A chr22_raw \
    -b 0 -B chr22_VQSR \
    -c 0 -C chr22_HardF \
	-d 51 -D chr22_GoldS \
	-e 128 \
	-f 0 \
	-G 0 \
	-j 0 \
	-k 0 \
	-l 0 \
	-m 49095 \
	-n 29 \
	-p 0 \
	-q 0 \
	-i 42110 \
	-t "NA12878 variant recall" -x 2 -o snpeff/recall4pc -u 1 -P 1

