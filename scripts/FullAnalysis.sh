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
#set -x
#exec > ${log} 2>&1
exec &> >(tee -i ${log})

# multithreading, adapt to your cpu count
bwathr=84
samtoolsthr=24

# how much RAM can java use? - set these to less than your available RAM!
javaopts="-Xms24g -Xmx24g"

# my samtools is here
samtools=$BIOTOOLS/samtools/bin/samtools

function mytest {
	"$@"
	local status=$?
	if [ ${status} -ne 0 ]; then
		echo "error with $1" >&2
		exit
	fi
	return ${status}
}

#############################################
# Get GATK Bundle files
#############################################

# Please run the separate script get_bundle.sh before attempting to run this one

#############################################
# Get chr22 reads from GIAB ftp
# source: NIST_NA12878_HG001_HiSeq_300x
#############################################

# Note: the chr22 read subsets should be first prepared as described next (this takes time!!)
# mkdir -p reads
# # download bam chr22 mappings
# samtools view -b -F 4 \
#	ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/NA12878/NIST_NA12878_HG001_HiSeq_300x/NHGRI_Illumina300X_novoalign_bams/HG001.GRCh38_full_plus_hs38d1_analysis_set_minus_alts.300x.bam \
#	chr22 > NHGRI_Illumina300X_novoalign_chr22.bam
#
# # sort by names and extract pairs to fastq
# samtools sort -n \
#	NHGRI_Illumina300X_novoalign_chr22.bam \
#	> NHGRI_Illumina300X_novoalign_chr22_nmsrt.bam
#
# # convert to fastq (discard unmapped and singletons)
# samtools fastq \
#		-F 0x900 \
#		-@ 24 \
#		-0 /dev/null \
#		-1 HG001.GRCh38_chr22_1.fq \
#		-2 HG001.GRCh38_chr22_2.fq \
#		-s /dev/null \
#		-n \
#		NHGRI_Illumina300X_novoalign_chr22_nmsrt.bam
# [M::bam2fq_mainloop] discarded 1508933 singletons
# [M::bam2fq_mainloop] processed 66909581 reads
#
# # cleanup
# rm ./*.bam*

#############################################
# BWA index & MEM mapping
#############################################

## create a bwa index when absent
bwaidxfolder=bwa_index
mkdir -p ${bwaidxfolder}

reference_fa=reference/Homo_sapiens_assembly38.fasta
bwaidx="${bwaidxfolder}/$(basename ${reference_fa})"

if [ ! -f ${bwaidxfolder}/index_created ]; then
	echo "# creating BWA index"
	bwa index ${reference_fa} -p ${bwaidx} \
  	&& touch ${bwaidxfolder}/index_created
else
	echo "# BWA index already exists"
fi


## map reads to reference
reads_1="reads/HG001.GRCh38_chr22_1.fq.gz"
reads_2="reads/HG001.GRCh38_chr22_2.fq.gz"

# mapping settings
thr=84
outpfx="HG001_chr22"
samplename="NA12878"

outfolder=bwa_mappings
mkdir -p ${workdir}/${outfolder}

# map using BWA mem (only once)
if [ ! -f bwa_mappings/mapping_done ]; then
	echo "# mapping reads with BWA mem"
	cmd="bwa mem -t ${bwathr} \
		-M \
		-R '@RG\tID:HG001\tLB:NA12878_giab\tPU:unknown-0.0\tPL:Illumina\tSM:NA12878' \
		${bwaidx} \
		${reads_1} \
		${reads_2} | ${samtools} view -b - -o ${outfolder}/${outpfx}_rawmappings.bam"
	echo "# ${cmd}"
	eval ${cmd}
	# get flagstats
	# then flag record that mapping was done for next run
	${samtools} flagstat -@ ${samtoolsthr} \
		${outfolder}/${outpfx}_rawmappings.bam \
		> ${outfolder}/${outpfx}_rawmappings_flagstats.txt && \
		touch bwa_mappings/mapping_done
else
	echo "# BWA mapping already done"
fi


#############################################
# PICARD Cleanup & MarkDuplicates
#############################################

# more records in RAM speeds up when enough RAM is present
recinram=10000000

if [ ! -f bwa_mappings/preprocessing_done ]; then
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
	${samtools} index ${outfolder}/${outpfx}_mrkdup_srt_22only.bam && \
		touch bwa_mappings/preprocessing_done
else
	echo "# GATK preprocessing already done"
fi

#############################################
# GATK4 BAM RECALIBRATION
#############################################

outfolder=gatk_preprocessing
mkdir -p ${outfolder}

bamfile=bwa_mappings/${outpfx}_mrkdup_srt_22only.bam
recalbamfile=${outpfx}_mrkdup_srt_recal.bam

# base quality score recalibration
knownsites="reference/dbsnp_138.hg38.vcf.gz"
knownindels="reference/Homo_sapiens_assembly38.known_indels.vcf.gz"

if [ ! -f gatk_preprocessing/recalibration_done ]; then
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
	--interval-padding 100 \
	--add-output-sam-program-record \
	--use-original-qualities \
	--static-quantized-quals 10 \
	--static-quantized-quals 20 \
	--static-quantized-quals 30 \
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
	O=${outfolder}/${outpfx}_multiple_metrics \
	R=${reference_fa} \
	MAX_RECORDS_IN_RAM=${recinram} \
	TMP_DIR=tmpfiles/ && \
	touch gatk_preprocessing/recalibration_done
else
	echo "# GATK BAM recalibration already done"
fi


#############################################
# GATK4 VARIANT CALLING (chr22 only)
#############################################

outfolder=gatk_variantcalling
mkdir -p ${outfolder}

# dbSNP for ID-field annotation
dbsnp146=reference/dbsnp_146.hg38.vcf.gz

if [ ! -f gatk_variantcalling/calling_done ]; then
# call short variants to gvcf format and save supporting reads

# parallel threads
hmmt=4

java ${javaopts} -jar $GATK/gatk.jar \
	HaplotypeCaller  \
	--input gatk_preprocessing/${recalbamfile} \
	--output ${outfolder}/${samplename}.g.vcf.gz \
	--reference ${reference_fa} \
	--emit-ref-confidence GVCF \
	--sample-ploidy 2 \
	--native-pair-hmm-threads ${hmmt} \
	--intervals chr22 \
	--bam-output ${outfolder}/${samplename}_HC_aligned_reads.bam \
	--tmp-dir tmpfiles/ 
	
	# convert to VCF
# https://software.broadinstitute.org/gatk/documentation/article?id=11813
java ${javaopts} -jar $GATK/gatk.jar \
	GenotypeGVCFs \
	--reference ${reference_fa} \
	--variant ${outfolder}/${samplename}.g.vcf.gz \
	--output ${outfolder}/${samplename}.vcf.gz \
	--intervals chr22 \
	--dbsnp ${dbsnp146} \
	--use-new-qual-calculator \
	--tmp-dir tmpfiles/
	
# mark ExcessHet
threshold=54.69
java ${javaopts} -jar $GATK/gatk.jar \
	VariantFiltration \
	--variant ${outfolder}/${samplename}.vcf.gz \
	--filter-expression "ExcessHet > ${threshold}" \
	--filter-name "ExcessHet" \
	-O ${outfolder}/${samplename}_excesshet_filtered.vcf.gz \
	--tmp-dir tmpfiles/ && \
	touch gatk_variantcalling/calling_done
else
	echo "# GATK calling already done"
fi


#############################################
# GATK4 VARIANT RECALIBRATION
#############################################

outfolder=gatk_variantrecalibration
mkdir -p ${outfolder}

# recalibration sources
## variants-sets

# True sites training resource: HapMap
truetraining15=reference/hapmap_3.3.hg38.vcf.gz

# True sites training resource: Omni
truetraining12=reference/1000G_omni2.5.hg38.vcf.gz

# Non-true sites training resource: 1000G
nontruetraining10=reference/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz

# Known sites resource, not used in training: dbSNP
knowntraining2=reference/dbsnp_138.hg38.vcf.gz

# indels True sites training resource: Mills
truetrainingindel12=reference/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz

# indels Axiom
axiom10=reference/Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz

# intervals to restrict analysis to chr22
intervals=reference/chr22.bed

# maxgaussians default to 8 is too high for chr22-only calls
maxSNPgaussians=6
maxINDELgaussians=4

if [ ! -f gatk_variantrecalibration/variantrecalibation_done ]; then
# copy last files from previous step
cp gatk_variantcalling/${samplename}_excesshet_filtered.vcf.gz* ${outfolder}/

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
	--intervals chr22 \
	--resource:hapmap,known=false,training=true,truth=true,prior=15.0 ${truetraining15} \
	--resource:omni,known=false,training=true,truth=true,prior=12.0 ${truetraining12} \
	--resource:1000G,known=false,training=true,truth=false,prior=10.0 ${nontruetraining10} \
	--resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ${knowntraining2} \
	--trust-all-polymorphic \
	--use-annotation DP \
	--use-annotation QD \
	--use-annotation FS \
	--use-annotation SOR \
	--use-annotation MQ \
	--use-annotation MQRankSum \
	--use-annotation ReadPosRankSum \
	--mode SNP \
	--max-gaussians ${maxSNPgaussians} \
	--tranches-file ${outfolder}/${samplename}_recalibrate_snp.tranches \
	--tranche 100.0 \
	--tranche 99.95 \
	--tranche 99.9 \
	--tranche 99.8 \
	--tranche 99.6 \
	--tranche 99.5 \
	--tranche 99.4 \
	--tranche 99.3 \
	--tranche 99.0 \
	--tranche 98.0 \
	--tranche 97.0 \
	--tranche 90.0 \
	--rscript-file ${outfolder}/${samplename}_recalibrate_snp_plots.R \
	--tmp-dir tmpfiles/


# Build the Indel recalibration model
java ${javaopts} -jar $GATK/gatk.jar \
	VariantRecalibrator \
	-R ${reference_fa} \
	-V ${outfolder}/${samplename}_excesshet_sitesonly.vcf.gz \
	-O ${outfolder}/${samplename}_recalibrate_indels.recal.vcf.gz \
	--intervals chr22 \
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
	--max-gaussians ${maxINDELgaussians} \
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
	-O ${outfolder}/${samplename}_recalibrated_snps_raw_indels.vcf.gz \
	--intervals chr22 \
	--mode SNP \
	--recal-file ${outfolder}/${samplename}_recalibrate_SNP.recal.vcf.gz \
	--tranches-file ${outfolder}/${samplename}_recalibrate_snp.tranches \
	--truth-sensitivity-filter-level 99.0 \
	--tmp-dir tmpfiles/

# Apply the desired level of recalibration to the Indels in the call set
java ${javaopts} -jar $GATK/gatk.jar \
	ApplyVQSR \
	-R ${reference_fa} \
	-V ${outfolder}/${samplename}_recalibrated_snps_raw_indels.vcf.gz \
	-O ${outfolder}/${samplename}_VQSR.vcf.gz \
	--intervals chr22 \
	-mode INDEL \
	--recal-file ${outfolder}/${samplename}_recalibrate_indels.recal.vcf.gz \
	--tranches-file ${outfolder}/${samplename}_recalibrate_indels.tranches \
	--truth-sensitivity-filter-level 99.7 \
	--create-output-variant-index true \
	--tmp-dir tmpfiles/ && \
	touch gatk_variantrecalibration/variantrecalibation_done
else
	echo "# GATK variant recalibration already done"
fi


#############################################
# SnpEff on VQSR
#############################################

outfolder="snpeff"
mkdir -p ${outfolder}

build="GRCh38.86"

if [ ! -f snpeff/VQSR_annotation_done ]; then
java ${javaopts} -jar $SNPEFF/snpEff.jar \
	-htmlStats ${outfolder}/gatk_recal_snpEff_summary.html \
	${build} \
	gatk_variantrecalibration/${samplename}_VQSR.vcf.gz | \
	bgzip -c > ${outfolder}/${samplename}_VQSR_snpeff.vcf.gz && \
	tabix -p vcf ${outfolder}/${samplename}_VQSR_snpeff.vcf.gz && \
	touch snpeff/VQSR_annotation_done
else
	echo "# SNPEff annotation for VQSR variant already done"
fi


#############################################
# GATK4 VARIANT HARD FILTERING
#############################################

# instructions and filters from:
# https://gatkforums.broadinstitute.org/gatk/discussion/23216/how-to-filter-variants-either-with-vqsr-or-by-hard-filtering
# Genepattern: QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || SOR > 3.0 || QUAL < 30

outfolder=gatk_varianthardfiltering
mkdir -p ${outfolder}

if [ ! -f gatk_varianthardfiltering/hardFiltering_recalibation_done ]; then

# copy the GenotypeGVCFs vcf output with marked excesshet output from above
# create local copy of previous file
cp gatk_variantcalling/${samplename}_excesshet_filtered.vcf.gz* ${outfolder}/


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
	-O ${outfolder}/${samplename}_snp_filtered.vcf.gz \
	--filter "QD < 2.0" --filter-name "QD2" \
	--filter "QUAL < 30.0" --filter-name "QUAL30" \
	--filter "SOR > 3.0" --filter-name "SOR3" \
	--filter "FS > 60.0" --filter-name "FS60" \
	--filter "MQ < 40.0" --filter-name "MQ40" \
	--filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
	--filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8"


######################################################
# 3) hard-Filter INDELs and MIXED on multiple metrics
######################################################

# produces a VCF with records with INDEL-type variants only.
java ${javaopts} -jar $GATK/gatk.jar \
	SelectVariants \
	-V ${outfolder}/${samplename}_excesshet_filtered.vcf.gz \
	-O ${outfolder}/${samplename}_mixed_indels.vcf.gz \
	--select-type-to-include INDEL \
	--select-type-to-include MIXED

java ${javaopts} -jar $GATK/gatk.jar \
	VariantFiltration \
	-V ${outfolder}/${samplename}_mixed_indels.vcf.gz \
	-O ${outfolder}/${samplename}_mixed_indels_filtered.vcf.gz \
	--filter "QD < 2.0" --filter-name "QD2" \
	--filter "QUAL < 30.0" --filter-name "QUAL30" \
	--filter "FS > 200.0" --filter-name "FS200" \
	--filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20"


###########################################
## 4) merge SNP and Indel filtered calls	
###########################################

# combine filtered into one file using Picard
java ${javaopts} -jar $GATK/gatk.jar \
	MergeVcfs \
	-I ${outfolder}/${samplename}_snp_filtered.vcf.gz \
	-I ${outfolder}/${samplename}_mixed_indels_filtered.vcf.gz \
	-R ${reference_fa} \
	-O ${outfolder}/${samplename}_snp_indel_filtered.vcf.gz && \
	touch ${outfolder}/hardFiltering_recalibation_done
else
	echo "# GATK hardFiltering recalibration already done"
fi


#############################################
# SnpEff on HardFiltering
#############################################

outfolder="snpeff"
mkdir -p ${outfolder}

build="GRCh38.86"

if [ ! -f snpeff/hardFiltering_annotation_done ]; then
java ${javaopts} -jar $SNPEFF/snpEff.jar \
	-htmlStats ${outfolder}/gatk_hardfiltering_snpEff_summary.html \
	${build} \
	gatk_varianthardfiltering/${samplename}_snp_indel_filtered.vcf.gz | \
	bgzip -c > ${outfolder}/${samplename}_hardfiltering_snpeff.vcf.gz && \
	tabix -p vcf ${outfolder}/${samplename}_hardfiltering_snpeff.vcf.gz && \
	touch snpeff/hardFiltering_annotation_done
else
	echo "# SNPEff annotation for hard-filtered variant already done"
fi


#############################################
# Compare calls to Golden
#############################################

gold=reference/NA12878_HG001-chr22_gold.vcf.gz
vqsr=snpeff/NA12878_VQSR_snpeff.vcf.gz
hard=snpeff/NA12878_hardfiltering_snpeff.vcf.gz

if [ ! -f snpeff/comparisons_done ]; then

# compare 3 VCF on PASS and . rows only
vcf-compare -a \
	-r chr22 \
	${vqsr} \
	${hard} \
	${gold} | grep ^VN | cut -f 2- > snpeff/3venn_counts.txt

# a.c 26	reference/NA12878_HG001-chr22_gold.vcf.gz (0.1%)	snpeff/NA12878_VQSR_snpeff.vcf.gz (0.0%)
# a.. 484	snpeff/NA12878_VQSR_snpeff.vcf.gz (0.9%)
# ..c 1569	reference/NA12878_HG001-chr22_gold.vcf.gz (3.7%)
# .bc 3070	reference/NA12878_HG001-chr22_gold.vcf.gz (7.3%)	snpeff/NA12878_hardfiltering_snpeff.vcf.gz (4.6%)
# .b. 11541	snpeff/NA12878_hardfiltering_snpeff.vcf.gz (17.1%)
# ab. 15328	snpeff/NA12878_VQSR_snpeff.vcf.gz (28.7%)	snpeff/NA12878_hardfiltering_snpeff.vcf.gz (22.7%)
# abc 37525	reference/NA12878_HG001-chr22_gold.vcf.gz (88.9%)	snpeff/NA12878_VQSR_snpeff.vcf.gz (70.3%)	snpeff/NA12878_hardfiltering_snpeff.vcf.gz (55.6%)

# create venn plot
# use R to plot
# https://wiki.bits.vib.be/index.php/NGS_Exercise.6#Create_a_VENN_diagram_for_the_4_mappings

3DVenn.R -a 484 -A chr22_VQSR \
	-b 11541 -B chr22_HardF \
	-c 1569 -C chr22_GoldS \
	-d 15328 \
	-e 26 \
	-f 3070 \
	-i 37525 \
	-t "NA12878 variant recall" \
	-x 2 \
	-o snpeff/recall3 \
	-u 1

3DVenn.R -a 484 -A chr22_VQSR \
	-b 11541 -B chr22_HardF \
	-c 1569 -C chr22_GoldS \
	-d 15328 \
	-e 26 \
	-f 3070 \
	-i 37525 \
	-t "NA12878 variant recall" \
	-x 2 \
	-o snpeff/recall3pc \
	-u 1 \
	-P 1 && \
	touch snpeff/comparisons_done
else
	echo "# VCF-tools comparisons already done"
fi

# cleanup leftovers
rm tmpfiles/*
