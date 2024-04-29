#================================= =========Stats=================== first we want to calculate some basic stats about the BAM files. Some samples are no good, we want to 
#exclude very low coverage we make a header using the stats information from the outgroup isolate, although any isolate would do just fine for this purpose rm bam_stats.tsv 
#echo -n -e "sample\t" >> bam_stats.tsv samtools stats bams/dutch.V254-50.sorted.bam | grep "^SN" | cut -f 2 | tr "\n" "\t" >> bam_stats.tsv echo >> bam_stats.tsv
  
#now we make a loop and walk through all the bams, keeping the sample name and the various stats
#for i in bams/*sorted.bam
#do
#echo -n -e ${i/\.sorted\.bam/}"\t" >> bam_stats.tsv
#samtools stats -@ 4 $i | grep "^SN" | cut -f 3 | tr "\n" "\t" >> bam_stats.tsv
#echo >> bam_stats.tsv
#done

## the same name is stored in column 1 the number of properly paired reads in 11 duplicated reads in 13 and length in 24 of bam_stats.tsv
## then we can compute the set of "good" samples by having 11-13*24/28000000, the proper pairs minus duplicates * length / genome length to get effective coverage

#=================================
#========Merging VCFs=============
#nothing fancy here, we just make a list of the files, and merge all of them
#but we remove the small ones, with low coverage

find bams/ -type f -name "*.gvcf.gz" | grep -v "barber.SRR10714231" | \
	grep -v "bams/dutch.V273-05" | grep -v "bams/etienne.SRR11785206" | \
	grep -v "bams/lofgren_pbio.SRR12763094" | grep -v "bams/lofgren_pbio.SRR12763095" | \
	grep -v "bams/lofgren_pbio.SRR12763102" | grep -v "bams/lofgren_pbio.SRR14584215" | \
	grep -v "bams/lofgren_pbio.SRR14584260" | grep -v "bams/lofgren_pbio.SRR14584267" | \
	grep -v "bams/lofgren_pbio.SRR14584274" | grep -v "bams/shelton.ERR12374166" | \
	grep -v "bams/shelton.ERR12374179" | grep -v "bams/lofgren_pbio.SRR14584249" > temp_names.list

mkdir global_vcf

/mnt/LTR_userdata/auxie001/programs/gatk-4.2.6.1/gatk CombineGVCFs \
       -O  global_vcf/global_raw.gvcf.gz \
       --variant temp_names.list \
       -R /mnt/LTR_userdata/auxie001/groeniii/genomes/Af293_combined.fna

/mnt/LTR_userdata/auxie001/programs/gatk-4.2.6.1/gatk GenotypeGVCFs \
       -O global_vcf/global_raw.vcf.gz \
       -R /mnt/LTR_userdata/auxie001/groeniii/genomes/Af293_combined.fna \
       -V global_vcf/global_raw.gvcf.gz

rm temp_names.list



#=================================
#========Filtering VCFs=============
##this script will do the hard filtering of the vcf files, as well as produce various QC information
#the main goal is to filter low quality variants. A recent publication used the following filters, which can be a starting point for our purposes
#Barber et al.,2021: (SNP: QD < 2.0 || MQ < 40.0 || FS > 60.0 || MQRankSum < −12.5 || ReadPosRankSum < −8.0; indel: QD < 2.0 || FS > 200.0 || ReadPosRankSum < −20.0)
#however, these filter settings were used on individual VCFs, and not the combined VCF. Barber et al. then merged the filtered VCFs, which is not recommended by GATK
#so we will need to adapt.

#first we need to split the VCF

/mnt/LTR_userdata/auxie001/programs/gatk-4.2.6.1/gatk SplitVcfs -I global_vcf/global_raw.vcf.gz \
	--SNP_OUTPUT global_vcf/global_SNPs.vcf.gz --INDEL_OUTPUT global_vcf/global_INDELs.vcf.gz --STRICT false

echo "now counting hets"
#this command will count the heterozygotes per line:
#it prints any site with more than 63 heterozygotes (more than 5% of sample) into a bed file, and then slops into windows where it also makes regions including the up- and downstream 50 basepairs
bcftools query -f "%CHROM\t%POS\t[%GT\t]\n" global_vcf/global_raw.vcf.gz | sed "s/0|1/0\/1/g" | awk -v OFS="\t" '{print $1,$2,gsub("0/1","")}' | awk -v OFS="\t" '{if ($3 > 63) print $1,$2,$2+1}' > global_vcf/het_variants.bed
bedtools slop -i global_vcf/het_variants.bed -b 50 -g /mnt/LTR_userdata/auxie001/groeniii/genomes/Af293_combined.genome | bedtools merge > global_vcf/het_regions.bed
rm global_vcf/het_variants.bed
echo "finished counting hets"

echo "the heterozygous regions cover a total of :"
cat global_vcf/het_regions.bed  | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}'

#when filtering, we also exculde sites where less than 90% of samples have called genotypes
echo "starting filtering SNPs"
bcftools filter -O u -s "QD"			-m + -e "QD < 2.0" global_vcf/global_SNPs.vcf.gz | \
	bcftools filter -O u -s "MQ"		-m + -e "MQ < 40.0" | \
	bcftools filter -O u -s "FS"		-m + -e "FS > 60.0" | \
	bcftools filter -O u -s "MQRank"	-m + -e "MQRankSum < -12.5" | \
	bcftools filter -O u -s "ReadPos"	-m + -e "ReadPosRankSum < -8.0" | \
	bcftools filter -O u -s "AN"		-m + -e "AN < 1800" | \
	bcftools view -O z -f "PASS" -T ^global_vcf/het_regions.bed > global_vcf/global_SNPs_filtered.vcf.gz
echo "finished filtering SNPs"

bcftools index global_vcf/global_SNPs_filtered.vcf.gz

echo "starting filtering INDELs"
bcftools filter -O u -s "QD"                    -m + -e "QD < 2.0" global_vcf/global_INDELs.vcf.gz | \
        bcftools filter -O u -s "FS"         -m + -e "FS > 200.0" | \
        bcftools filter -O u -s "ReadPos"    -m + -e "ReadPosRankSum < -20.0" | \
        bcftools filter -O u -s "AN"         -m + -e "AN < 1800" | \
        bcftools view -O z -f "PASS" -T ^global_vcf/het_regions.bed > global_vcf/global_INDELs_filtered.vcf.gz
echo "finished filtering INDELs"

bcftools index global_vcf/global_INDELs_filtered.vcf.gz

bcftools concat -O z -a global_vcf/global_INDELs_filtered.vcf.gz global_vcf/global_SNPs_filtered.vcf.gz > global_vcf/global_combined_allsamples.vcf.gz
bcftools stats -s - global_vcf/global_combined_allsamples.vcf.gz > global_vcf/global_stats_allsamples.txt

#========================================================
#==========Removing highly divergent samples=============
#========================================================

#this will print the samples of the global dataset with more than 10% of variant sites being heterozygous
grep "PSC" global_vcf/global_stats_allsamples.txt | grep -v "#" | awk '{if ($6/($5+$6) > 0.1) print $0,$6/($5+$6)}' | cut -f 3 > het_exclude_list.txt
#and now exclude them
bcftools view -S ^het_exclude_list.txt  global_vcf/global_combined_allsamples.vcf.gz | \
	bcftools view -O z --min-ac 1 > global_vcf/global_combined_partial_filtered_nohets.vcf.gz

#we also want to remove the highly divergent outgroups, which can be found in the stats file:
grep "PSC" global_vcf/global_stats_allsamples.txt | grep -v "#" | awk '{if ($5>60000) print $3}' > divergent_exclude_list.txt

#so exclude these sites, as well as V273-05 which has a significant amount of missing data.
#need to use --force-samples since some samples have been excluded already by the heterozygous filter
bcftools view --force-samples -S ^divergent_exclude_list.txt  global_vcf/global_combined_partial_filtered_nohets.vcf.gz | \
	bcftools view -O z --min-ac 1 > global_vcf/global_combined_filtered.vcf.gz
bcftools index global_vcf/global_combined_filtered.vcf.gz
bcftools stats -s - global_vcf/global_combined_filtered.vcf.gz > global_vcf/global_stats.txt

rm het_exclude_list.txt divergent_exclude_list.txt

