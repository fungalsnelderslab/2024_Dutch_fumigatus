#this script will do the hard filtering of the vcf files, as well as produce various QC information

#the main goal is to filter low quality variants. A recent publication used the following filters, which can be a starting point for our purposes
#Barber et al.,2021: (SNP: QD < 2.0 || MQ < 40.0 || FS > 60.0 || MQRankSum < −12.5 || ReadPosRankSum < −8.0; indel: QD < 2.0 || FS > 200.0 || ReadPosRankSum < −20.0)
#however, these filter settings were used on individual VCFs, and not the combined VCF. Barber et al. then merged the filtered VCFs, which is not recommended by GATK
#so we will need to adapt.

#first we need to split the VCF

/mnt/LTR_userdata/auxie001/programs/gatk-4.2.6.1/gatk SplitVcfs -I dutch_vcf/dutch_raw.vcf.gz \
	--SNP_OUTPUT dutch_vcf/dutch_SNPs.vcf.gz --INDEL_OUTPUT dutch_vcf/dutch_INDELs.vcf.gz --STRICT false

echo "now counting hets"
#this command will count the heterozygotes per line:
#
bcftools query -f "%CHROM\t%POS\t[%GT\t]\n" dutch_vcf/dutch_raw.vcf.gz | sed "s/0|1/0\/1/g" | awk -v OFS="\t" '{print $1,$2,$2,gsub("0/1","")}' | awk '{if ($4 > 20) print $0}' > dutch_vcf/het_variants.bed
bedtools slop -i dutch_vcf/het_variants.bed -b 100 -g /mnt/LTR_userdata/auxie001/groeniii/genomes/Af293_combined.genome | bedtools merge > dutch_vcf/het_regions.bed
rm dutch_vcf/het_variants.bed
echo "finished counting hets"

echo "the heterozygous regions cover a total of :"
cat dutch_vcf/het_regions.bed  | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}'

echo "starting filtering SNPs"
bcftools filter -O u -s "QD"			-m + -e "QD < 2.0" dutch_vcf/dutch_SNPs.vcf.gz | \
	bcftools filter -O u -s "MQ"		-m + -e "MQ < 40.0" | \
	bcftools filter -O u -s "FS"		-m + -e "FS > 60.0" | \
	bcftools filter -O u -s "MQRank"	-m + -e "MQRankSum < -12.5" | \
	bcftools filter -O u -s "ReadPos"	-m + -e "ReadPosRankSum < -8.0" | \
	bcftools filter -O u -s "AN"		-m + -e "AN < 300" | \
	bcftools view -O z -f "PASS" -T ^dutch_vcf/het_regions.bed > dutch_vcf/dutch_SNPs_filtered.vcf.gz
echo "finished filtering SNPs"

bcftools index dutch_vcf/dutch_SNPs_filtered.vcf.gz

echo "starting filtering INDELs"
bcftools filter -O u -s "QD"                    -m + -e "QD < 2.0" dutch_vcf/dutch_INDELs.vcf.gz | \
        bcftools filter -O u -s "FS"         -m + -e "FS > 200.0" | \
        bcftools filter -O u -s "ReadPos"    -m + -e "ReadPosRankSum < -20.0" | \
        bcftools filter -O u -s "AN"         -m + -e "AN < 300" | \
        bcftools view -O z -f "PASS" -T ^dutch_vcf/het_regions.bed > dutch_vcf/dutch_INDELs_filtered.vcf.gz
echo "finished filtering INDELs"

bcftools index dutch_vcf/dutch_INDELs_filtered.vcf.gz

bcftools concat -O z -a dutch_vcf/dutch_INDELs_filtered.vcf.gz dutch_vcf/dutch_SNPs_filtered.vcf.gz > dutch_vcf/dutch_combined_allsamples.vcf.gz
bcftools stats -s - dutch_vcf/dutch_combined_allsamples.vcf.gz > dutch_vcf/dutch_stats_allsamples.txt

#this will print the lines of the VCF with more than 10% of variant sites are heterozygous
grep "PSC" dutch_vcf/dutch_stats_allsamples.txt |  awk '{if ($6/($5+1) > 0.1) print $0,$6/($5+1)}'

#so exclude these sites, as well as V273-05 which has a significant amount of missing data.
bcftools view -O u -s ^2108,27C7,50C32,513,72A19,72A8,76A17,V209-50,V215-11,V219-48,V224-66,V226-13,V228-47,V237-62,V265-56,V293-59,V273-05,V254-50,out.Afis_SRR11363404,out.Afis_SRR11363405,out.Aoer_SRR11363406 dutch_vcf/dutch_combined_allsamples.vcf.gz | \
	bcftools view -O z --min-ac 1 > dutch_vcf/dutch_combined_filtered.vcf.gz
bcftools index dutch_vcf/dutch_combined_filtered.vcf.gz
bcftools stats -s - dutch_vcf/dutch_combined_filtered.vcf.gz > dutch_vcf/dutch_stats.txt

bcftools annotate -x "INFO" dutch_vcf/dutch_combined_filtered.vcf.gz | bcftools annotate -x "FORMAT" > dutch_vcf/dutch_combined_filtered.vcf
