/mnt/LTR_userdata/auxie001/programs/vcftools/bin/vcftools --gzvcf dutch_combined_filtered.vcf.gz --weir-fst-pop susceptible.txt --weir-fst-pop resistant.txt   --fst-window-size 10000 --fst-window-step 5000 --out fst_res_susc
/mnt/LTR_userdata/auxie001/programs/vcftools/bin/vcftools --gzvcf dutch_combined_filtered.vcf.gz --weir-fst-pop susceptible.txt --weir-fst-pop resistant.txt --out fst_res_susc.nowindowed

/mnt/LTR_userdata/auxie001/programs/vcftools/bin/vcftools --gzvcf dutch_combined_filtered.vcf.gz --weir-fst-pop environmental_samples.txt --weir-fst-pop clinical_samples.txt --fst-window-size 10000 --fst-window-step 5000 --out fst_clin_env


#making PCA plot

/mnt/LTR_userdata/auxie001/programs/plink --vcf dutch_combined_filtered.vcf.gz --make-bed --out dutch_plink --allow-extra-chr

/mnt/LTR_userdata/auxie001/programs/plink --allow-extra-chr --bfile dutch_plink --maf 0.01 --indep-pairwise 50 5 0.2 --out dutch_data_clean
/mnt/LTR_userdata/auxie001/programs/plink --allow-extra-chr --bfile dutch_plink --extract dutch_data_clean.prune.in --make-bed --out dutch_data_clean_prune

/mnt/LTR_userdata/auxie001/programs/plink --allow-extra-chr --bfile dutch_data_clean_prune --pca --out dutch_data_clean_prune

#wget https://raw.githubusercontent.com/edgardomortiz/vcf2phylip/master/vcf2phylip.py

python3 /mnt/scratch/auxie001/groeniii/vcf2phylip.py --input dutch_combined_filtered.vcf.gz --nexus

