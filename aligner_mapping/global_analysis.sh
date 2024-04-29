mkdir global_analysis
cd global_analysis

#==================================================
#===============                ===================
#===============   Phylogeny    ===================
#===============                ===================
#==================================================

#first we want to make a phylogeny
#we need to exclude two samples of A. spinosus from the Celia-Sanchez et al samples.
/mnt/LTR_userdata/auxie001/programs/vcftools/bin/vcftools --gzvcf ../global_vcf/global_combined_partial_filtered_nohets.vcf.gz --thin 1000 --stdout --recode --recode-INFO-all | \
bcftools annotate -x ^FORMAT/GT -O z -o ../global_vcf/global_combined_partial_filtered_nohets_smaller.vcf.gz
python3 /mnt/scratch/auxie001/groeniii/vcf2phylip.py --input ../global_vcf/global_combined_partial_filtered_nohets_smaller.vcf.gz --nexus
iqtree -s global_combined_partial_filtered_nohets_smaller.min4.phy -fast -m JC

#making FST plots
#vcftools --gzvcf ../global_vcf/global_combined_filtered.vcf.gz --weir-fst-pop susceptible.txt --weir-fst-pop resistant.txt   --fst-window-size 10000 --fst-window-step 5000 --out fst_res_susc
#vcftools --gzvcf ../global_vcf/global_combined_filtered.vcf.gz --weir-fst-pop environmental_samples.txt --weir-fst-pop clinical_samples.txt --fst-window-size 10000 --fst-window-step 50000 --out fst_clin_env

#making splitstree data
#wget https://raw.githubusercontent.com/edgardomortiz/vcf2phylip/master/vcf2phylip.py
python3 /mnt/scratch/auxie001/groeniii/vcf2phylip.py --input ../global_vcf/global_combined_filtered.vcf.gz --nexus

#we want to make two subsetted Splitstree plots:
#first we want 50% resistance
grep "Susceptible" ../global_dataset_metadata.csv | shuf -n 200 | cut -f 2 > fifty_percent_list.txt
grep "Resistant" ../global_dataset_metadata.csv | shuf -n 200 | cut -f 2 >> fifty_percent_list.txt
#now extract these sampes, and subset to MAF > 0.05
bcftools view --samples-file fifty_percent_list.txt --force-samples -O z -o global_fifty_percent.vcf.gz -q 0.05:minor ../global_vcf/global_combined_filtered.vcf.gz
python3 /mnt/scratch/auxie001/groeniii/vcf2phylip.py --input global_fifty_percent.vcf.gz --nexus

#first we want 15% resistance
grep "Susceptible" ../global_dataset_metadata.csv | shuf -n 355 | cut -f 2 > fifteen_percent_list.txt
grep "Resistant" ../global_dataset_metadata.csv | shuf -n 45 | cut -f 2 >> fifteen_percent_list.txt
#now extract these sampes, and subset to MAF > 0.05
bcftools view --samples-file fifteen_percent_list.txt --force-samples -O z -o global_fifteen_percent.vcf.gz -q 0.05:minor ../global_vcf/global_combined_filtered.vcf.gz
python3 /mnt/scratch/auxie001/groeniii/vcf2phylip.py --input global_fifteen_percent.vcf.gz --nexus

#now make 5% which reflects airsampling of Shelton and Kortenbosch
grep "Susceptible" ../global_dataset_metadata.csv | shuf -n 385 | cut -f 2 > five_percent_list.txt
grep "Resistant" ../global_dataset_metadata.csv | shuf -n 15 | cut -f 2 >> five_percent_list.txt
#now extract these sampes, and subset to MAF > 0.05
bcftools view --samples-file five_percent_list.txt --force-samples -O z -o global_five_percent.vcf.gz -q 0.05:minor ../global_vcf/global_combined_filtered.vcf.gz
python3 /mnt/scratch/auxie001/groeniii/vcf2phylip.py --input global_five_percent.vcf.gz --nexus

rm *percent_list.txt

echo "starting plink steps"
#making PCA plot using plink
/mnt/LTR_userdata/auxie001/programs/plink --vcf ../global_vcf/global_combined_filtered.vcf.gz --make-bed --out global_plink --allow-extra-chr -recode
#now we remake the datafile but with no LD pruning to look for clones
#the resulting file plink.genome has the relatedness values inside of it
/mnt/LTR_userdata/auxie001/programs/plink --file global_plink --allow-extra-chr --genome
#thin the data using ld pruning for admixture
/mnt/LTR_userdata/auxie001/programs/plink --allow-extra-chr --bfile global_plink --maf 0.01 --indep-pairwise 50 5 0.2 --out global_data
/mnt/LTR_userdata/auxie001/programs/plink --allow-extra-chr --bfile global_plink --extract global_data.prune.in --make-bed --out global_data_thin
/mnt/LTR_userdata/auxie001/programs/plink --allow-extra-chr --bfile global_data_thin --pca --out global_data_thin
#echo "finished plink"

#now we want to produce the input file for hmmIBD, we need to replace missing data with -1, reference with 0, and alternate as 1,2,3 etc.
bcftools query -f "%CHROM\t%POS[\t%GT]\n" ../global_vcf/global_combined_filtered.vcf.gz | tr "|" "/" | \
	 sed "s/0\/0/0/g" | sed "s/1\/1/1/g" | sed "s/2\/2/2/g" | sed "s/3\/3/3/g" | sed "s/4\/4/4/g" | sed "s/\.\/\./-1/g" | \
         sed "s/NC_007194.1/1/g" | sed "s/NC_007195.1/2/g" | sed "s/NC_007196.1/3/g" | sed "s/NC_007197.1/4/g" | \
	sed "s/NC_007198.1/5/g" | sed "s/NC_007199.1/6/g" | sed "s/NC_007200.1/7/g" | sed "s/NC_007201.1/8/g" | sed "s/JQ346808.1/9/g" > global_vcf_for_hmmIBD.tsv

#echo "starting hmmIBD"
python hmmIBD/vcf2hmm.py ../global_vcf/global_combined_filtered.vcf.gz global_vcf_for_hmmIBD
python hmmIBD/thin_sites.py global_vcf_for_hmmIBD_freq.txt global_vcf_for_hmmIBD_thin_list.txt
python hmmIBD/vcf2hmm.py -l global_vcf_for_hmmIBD_thin_list.txt ../global_vcf/global_combined_filtered.vcf.gz global_vcf_for_hmmIBD_thinned
hmmIBD/hmmIBD_updated -i global_vcf_for_hmmIBD_thinned_seq.txt -o global_vcf_for_hmmIBD_thinned -f global_vcf_for_hmmIBD_thinned_freq.txt

echo "finished hmmIBD"

mkdir admixture


rm antifungal_mutations.tsv
touch antifungal_mutations.tsv
#extract known fungicide mutations based on sample names
bcftools query -l ../global_vcf/global_combined_partial_filtered_nohets.vcf.gz | tr "\n" "\t" > antifungal_mutations.tsv
echo "" >> antifungal_mutations.tsv
#cytBG213A
bcftools view ../global_vcf/global_combined_partial_filtered_nohets.vcf.gz JQ346808.1:428 | bcftools query -f "[%TGT\t]\n" >> antifungal_mutations.tsv
#sdhB H270Y is C->T
bcftools view ../global_vcf/global_combined_partial_filtered_nohets.vcf.gz NC_007198.1:2654913 | bcftools query -f "[%TGT\t]\n" >> antifungal_mutations.tsv
#benA F219Y 
bcftools view ../global_vcf/global_combined_partial_filtered_nohets.vcf.gz NC_007194.1:2849059 | bcftools query -f "[%TGT\t]\n" >> antifungal_mutations.tsv
#cyp51A TR34
bcftools view ../global_vcf/global_combined_partial_filtered_nohets.vcf.gz NC_007197.1:1782107 | bcftools query -f "[%TGT\t]\n" >> antifungal_mutations.tsv
#cyp51A TR46
bcftools view ../global_vcf/global_combined_partial_filtered_nohets.vcf.gz NC_007197.1:1782102 | bcftools query -f "[%TGT\t]\n" >> antifungal_mutations.tsv
#msh6 G240A
bcftools view ../global_vcf/global_combined_partial_filtered_nohets.vcf.gz NC_007197.1:2148956 | bcftools query -f "[%TGT\t]\n" >> antifungal_mutations.tsv
Rscript ../R_transposing.R


#now we want to annotate the variants:
/mnt/LTR_userdata/auxie001/programs/snpEff/exec/snpeff eff -ud 250 -csvStats global_combined_filtered.stats.csv Aspergillus_fumigatus ../global_vcf/global_combined_filtered.vcf.gz > global_combined_filtered.snpEFF.vcf.gz


#admixture doesn't like non-standard chromosome names, so rename
sed -i "s/NC_007194.1/1/g" global_data_thin.bim
sed -i "s/NC_007195.1/2/g" global_data_thin.bim
sed -i "s/NC_007196.1/3/g" global_data_thin.bim
sed -i "s/NC_007197.1/4/g" global_data_thin.bim
sed -i "s/NC_007198.1/5/g" global_data_thin.bim
sed -i "s/NC_007199.1/6/g" global_data_thin.bim
sed -i "s/NC_007200.1/7/g" global_data_thin.bim
sed -i "s/NC_007201.1/8/g" global_data_thin.bim
sed -i "s/JQ346808.1/9/g" global_data_thin.bim

cd admixture
for k in 2 3 4 5 6 7 8 9 10
do /mnt/LTR_userdata/auxie001/programs/admixture_linux-1.3.0/admixture -j24 -s 12345  --cv ../global_data_thin.bed $k | tee log{$k}.run1.out
done

for k in 2 3 4 5 6 7 8 9 10
do /mnt/LTR_userdata/auxie001/programs/admixture_linux-1.3.0/admixture -j8 -s 42  --cv global_data_prune.bed $k | tee log{$k}.run2.out
done

for k in 2 3 4 5 6 7 8 9 10
do /mnt/LTR_userdata/auxie001/programs/admixture_linux-1.3.0/admixture -j8 -s 82507  --cv global_data_prune.bed $k | tee log{$k}.run3.out
done

for k in 2 3 4 5 6 7 8 9 10
do /mnt/LTR_userdata/auxie001/programs/admixture_linux-1.3.0/admixture -j8 -s 451309  --cv global_data_prune.bed $k | tee log{$k}.run4.out
done

for k in 2 3 4 5 6 7 8 9 10
do /mnt/LTR_userdata/auxie001/programs/admixture_linux-1.3.0/admixture -j8 -s 441479  --cv global_data_prune.bed $k | tee log{$k}.run5.out
done

for k in 2 3 4 5 6 7 8 9 10
do /mnt/LTR_userdata/auxie001/programs/admixture_linux-1.3.0/admixture -j8 -s 65353  --cv global_data_prune.bed $k | tee log{$k}.run6.out
done
