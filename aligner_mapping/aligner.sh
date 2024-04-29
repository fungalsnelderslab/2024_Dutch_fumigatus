BAMS=/mnt/scratch/auxie001/groeniii/bams

READS=/mnt/scratch/auxie001/groeniii/barber
##now process the samples from Barber et al.
#for i in $READS/*_1.fastq.gz
#do
#BARBER_SAMPLE=$(basename ${i/_1.fastq.gz/})
#echo $BARBER_SAMPLE
#if test -f $BAMS/barber.$BARBER_SAMPLE.gvcf.gz; then
#        echo "already finished vcf file"
#fi
#if ! test -f $BAMS/barber.$BARBER_SAMPLE.gvcf.gz; then
#        echo "doesn't exist, now generating"
#	/mnt/LTR_userdata/auxie001/programs/bwa-mem2-2.0pre2_x64-linux/bwa-mem2 mem -t 4 -R '@RG\tID:'$BARBER_SAMPLE'\tSM:'$BARBER_SAMPLE'\tLB:barber' -v 1 /mnt/LTR_userdata/auxie001/groeniii/genomes/Af293_combined.fna $READS/$BARBER_SAMPLE\_1.fastq.gz $READS/$BARBER_SAMPLE\_2.fastq.gz | \
#	samtools view -@ 4 -b | samtools fixmate -@ 4 -m - - | samtools sort -m 3G -@ 4 - | samtools markdup - $BAMS/barber.$BARBER_SAMPLE.sorted.bam
#	samtools index -@ 4 $BAMS/barber.$BARBER_SAMPLE.sorted.bam
#        /mnt/LTR_userdata/auxie001/programs/gatk-4.2.6.1/gatk HaplotypeCaller --seconds-between-progress-updates 120 -ERC GVCF -R /mnt/LTR_userdata/auxie001/groeniii/genomes/Af293_combined.fna \
#         -I $BAMS/barber.$BARBER_SAMPLE.sorted.bam -O $BAMS/barber.$BARBER_SAMPLE.gvcf.gz
#fi
#sleep 1s
#done

READS=/mnt/scratch/auxie001/groeniii/rhodes
#now process the samples from Rhodes et al.
#for i in $READS/ERR9791746_1.fastq.gz
#do
#RHODES_SAMPLE=$(basename ${i/_1.fastq.gz/})
#echo $RHODES_SAMPLE
#/mnt/LTR_userdata/auxie001/programs/bwa-mem2-2.0pre2_x64-linux/bwa-mem2 mem -t 4 -R '@RG\tID:'$RHODES_SAMPLE'\tSM:'$RHODES_SAMPLE'\tLB:rhodes' -v 1 /mnt/LTR_userdata/auxie001/groeniii/genomes/Af293_combined.fna $READS/$RHODES_SAMPLE\_1.fastq.gz $READS/$RHODES_SAMPLE\_2.fastq.gz | \
# samtools view -@ 4 -b | samtools fixmate -@ 4 -m - - | samtools sort -m 3G -@ 4 - | samtools markdup - $BAMS/rhodes.$RHODES_SAMPLE.sorted.bam
#samtools index -@ 4 $BAMS/rhodes.$RHODES_SAMPLE.sorted.bam
#done


READS=/mnt/scratch/auxie001/groeniii/kang
##now process the samples from Kang et al.
#for i in $READS/*_1.fastq.gz
#do
#KANG_SAMPLE=$(basename ${i/_1.fastq.gz/})
#echo $KANG_SAMPLE
#bwa-mem2 mem -v 1 -t 12 -R '@RG\tID:'$KANG_SAMPLE'\tSM:'$KANG_SAMPLE'\tLB:kang_et_al' -v 1 /mnt/LTR_userdata/auxie001/groeniii/genomes/Af293_combined.fna $READS/$KANG_SAMPLE\_1.fastq.gz $READS/$KANG_SAMPLE\_2.fastq.gz | \
# samtools view -@ 6 -b | samtools fixmate -@ 6 -m - - | samtools sort -m 3G -@ 12 - | samtools markdup -@ 6 - $BAMS/kang.$KANG_SAMPLE.sorted.bam
#samtools index $BAMS/kang.$KANG_SAMPLE.sorted.bam
#done

READS=/mnt/scratch/auxie001/groeniii/other
#get A. oerlinhausenensis
#wget -O $READS/Aoer_SRR11363406_1.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR113/006/SRR11363406/SRR11363406_1.fastq.gz
#wget -O $READS/Aoer_SRR11363406_2.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR113/006/SRR11363406/SRR11363406_2.fastq.gz

#get some A. fischeri
#wget -O $READS/Afis_SRR11363404_1.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR113/004/SRR11363404/SRR11363404_1.fastq.gz
#wget -O $READS/Afis_SRR11363404_2.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR113/004/SRR11363404/SRR11363404_2.fastq.gz

#wget -O $READS/Afis_SRR11363405_1.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR113/005/SRR11363405/SRR11363405_1.fastq.gz
#wget -O $READS/Afis_SRR11363405_2.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR113/005/SRR11363405/SRR11363405_2.fastq.gz

#for i in $READS/Aoer_SRR11363406_1.fastq.gz $READS/Afis_SRR11363404_1.fastq.gz $READS/Afis_SRR11363405_1.fastq.gz
#do
#OTHER_SAMPLE=$(basename ${i/_1.fastq.gz/})
#echo $OTHER_SAMPLE
#bwa-mem2 mem -v 1 -t 12 -R '@RG\tID:out.'$OTHER_SAMPLE'\tSM:out.'$OTHER_SAMPLE'\tLB:outgroup' /mnt/LTR_userdata/auxie001/groeniii/genomes/Af293_combined.fna $READS/$OTHER_SAMPLE\_1.fastq.gz $READS/$OTHER_SAMPLE\_2.fastq.gz | \
# samtools view -@ 6 -b | samtools fixmate -@ 6 -m - - | samtools sort -m 3G -@ 12 - | samtools markdup -@ 6 - $BAMS/out.$OTHER_SAMPLE.sorted.bam
#samtools index $BAMS/out.$OTHER_SAMPLE.sorted.bam
#done


##now make VCFs in GVCF format first
#mkdir /mnt/scratch/auxie001/groeniii/vcfs
#/mnt/LTR_userdata/auxie001/programs/gatk-4.2.6.1/gatk CreateSequenceDictionary -R /mnt/LTR_userdata/auxie001/groeniii/genomes/Af293_combined.fna

#for i in $BAMS/out*.sorted.bam $BAMS/dutch*.sorted.bam $BAMS/rhodes*.sorted.bam $BAMS/kang*.sorted.bam $BAMS/barber*.sorted.bam $BAMS/kang*.sorted.bam
#for i in $BAMS/ZHAO*.sorted.bam $BAMS/lofgren_bird*.sorted.bam $BAMS/lofgren_pbio*.sorted.bam
#for i in $BAMS/barber.SRR13579407.sorted.bam
#do /mnt/LTR_userdata/auxie001/programs/gatk-4.2.6.1/gatk HaplotypeCaller --seconds-between-progress-updates 120 -ERC GVCF -R /mnt/LTR_userdata/auxie001/groeniii/genomes/Af293_combined.fna \
# -I $i -O ${i/sorted.bam/gvcf.gz}
#done

#change to using CombineGVCFs
find bams/ -type f -name "*.gvcf.gz" > temp_names.list
#grep -e "dutch" -e "out" temp_names.list > dutch_subset.list

#mkdir global_vcf
#mkdir dutch_vcf

/mnt/LTR_userdata/auxie001/programs/gatk-4.2.6.1/gatk CombineGVCFs \
	-O  global_vcf/global_raw.gvcf.gz \
	--variant temp_names.list \
	-R /mnt/LTR_userdata/auxie001/groeniii/genomes/Af293_combined.fna

#/mnt/LTR_userdata/auxie001/programs/gatk-4.2.6.1/gatk CombineGVCFs \
#	-O dutch_vcf/dutch_raw.gvcf.gz \
#	--variant dutch_subset.list \
#	-R /mnt/LTR_userdata/auxie001/groeniii/genomes/Af293_combined.fna

/mnt/LTR_userdata/auxie001/programs/gatk-4.2.6.1/gatk GenotypeGVCFs \
	-O global_vcf/global_raw.vcf.gz \
	-R /mnt/LTR_userdata/auxie001/groeniii/genomes/Af293_combined.fna \
	-V global_vcf/global_raw.gvcf.gz

#/mnt/LTR_userdata/auxie001/programs/gatk-4.2.6.1/gatk GenotypeGVCFs \
#	-O dutch_vcf/dutch_raw.vcf.gz \
#	-R /mnt/LTR_userdata/auxie001/groeniii/genomes/Af293_combined.fna \
#	-V dutch_vcf/dutch_raw.gvcf.gz

rm temp_names.list dutch_subset.list
