# KMAR population genomics	

###	Reference Genome							

We are adding the mitochondrial genome from the Tatarenkov 2015( GennBank: KT893707.1). 
I just concatenated the mitochondrial genome (in fasta) to the end of the genome file.
```sh
cat ../ LHSH01_cat.fa Kmar_KT893707_FDS08.fasta > LHSH01_cat_MITO.fa
```
I then replaced the “|” for “_” in the name of the contigs. This is to avoid problems with feeding the contig names in the downstream analyses.
Path to the reference in santacruz: /media/santacruz5/luana_projects/kmar_ref/REFERENCES/as2_wsu/mitochondria/LHSH01_cat_MITO.fa 


### Raw data						
Path to the raw data: 
```sh
/media/santacruz5/luana_projects/kmar_ref/1.1_fastq_rawData/1_ready_for_use_files
```

## Quality Control and Trimming 

### Running Fastqc
```sh
fastqc --outdir ../1.2_fastqc_output1 --nogroup --noextract *fastq.gz
```
Check the results of fastqc and decide if the data needs hard trimming or not

### Running TrimGalore 

Not all the files needed to be hard trimmed so I created a list of the files that needed to be hard trimmed and one for the ones that didn't.

Script for Hard Trimming: TrimGalore_Lua_cat_NEW_CLIP.sh
```sh
# this script loops through the file fastq_files_CLIP.txt, where the first column is the read 1, the second is the read 2 of a given sample, and the third column is the length of trimming
 
while read line || [ -n "$line" ];
do
        read1=$(echo $line | awk '{print $1}')
        read2=$(echo $line | awk '{print $2}')
		clip=$(echo $line | awk '{print $3}')
		
/home/repository_software/trim_galore_v0.4.4/trim_galore \
        	--fastqc_args "--noextract --nogroup --outdir /media/santacruz5/luana_projects/kmar_ref/2.2_fastqc_trimmed_28" \
        	--stringency 5 \
        	--gzip \
        	--quality 28\
        	--length 50 \
        	--output_dir /media/santacruz5/luana_projects/kmar_ref/2.1_trimgalore_28 \
        	--clip_R1 $clip \
        	--clip_R2 $clip \
        	--paired $read1 $read2 
        	
done < fastq_files_CLIP.txt
```

Script for Trimming without hard trimmining: TrimGalore_Lua_cat_NEW_noclip.sh
```sh
# this script loops through the file fastq_files_noclip.txt, where the first column is the read 1, the second is the read 2.

# Runs Trim Galore on the files specified in the "TrimGalore_fastq_Locations.txt".  Note: this is a tab delimited file.
# "awk is grabbing the 1st, 2nd, or 3rd column from the line above.
# For the while loop the "|| [ -n "$line" ]" Checks if the given string operand size is non-zero. If it is non-zero
# length then it returns true. We're using it to read the last line of the input file just in case it isn't
# terminated with a newline.
while read line || [ -n "$line" ];
do
        read1=$(echo $line | awk '{print $1}')
        read2=$(echo $line | awk '{print $2}')
        trim_galore \
        	--fastqc_args "--noextract --nogroup --outdir /media/santacruz5/luana_projects/kmar_ref/2.2_fastqc_trimmed_28" \
        	--stringency 5 \
        	--gzip \
        	--quality 28\
        	--length 50 \
        	--output_dir /media/santacruz5/luana_projects/kmar_ref/2.1_trimgalore_28 \
        	--paired $read1 $read2 
        	
done < fastq_files_noclip.txt
```

Getting trimming stats
```sh
cd /media/santacruz5/luana_projects/kmar_ref/2.1_trimgalore_28

grep "Quality-trimmed" *_trimming_report.txt
grep "Total basepairs processed" *_trimming_report.txt 
grep "Total written" *_trimming_report.txt 
grep "Total reads processed" *_trimming_report.txt 
grep "Reads with adapters" *_trimming_report.txt 
grep "passing filters" *_trimming_report.txt
```
Inspect Fasqc of trimmed files.


## Mapping with BWA 


Building a index for the new reference:
```sh
bwa index /media/santacruz5/luana_projects/kmar_ref/REFERENCES/as2_wsu/mitochondria/LHSH01_cat_MITO.fa
```
Dictionary for the new reference:
```sh
java -jar /data/kelley/projects/programs/picard-tools-1.141/picard.jar \
CreateSequenceDictionary \
          R= /data/kelley/projects/luana_projects/kmar_references/mitochondria_wsu/LHSH01_cat_MITO.fa \
          O= /data/kelley/projects/luana_projects/kmar_references/mitochondria_wsu/LHSH01_cat_MITO.dict
```
#### Creating a file with input names for bwa with the trimmed files

Directory:
```sh
/media/santacruz5/luana_projects/kmar_ref/1.1_fastq_rawData/1_ready_for_use_files
```
Substitute the suffixes of the files:
```sh
cat fastq_files_CLIP.txt fastq_files_noclip.txt > temp.txt
awk '{print $1, $2}' temp.txt > fastq_files.txt
sed 's/_1.q33.fastq.gz/_1.q33_val_1.fq.gz/g' fastq_files.txt | sed 's/_2.q33.fastq.gz/_2.q33_val_2.fq.gz/g' | sed 's/_1_cat.gz/_1_cat.gz_val_1.fq.gz/g' | sed 's/_2_cat.gz/_2_cat.gz_val_2.fq.gz/g' | sed 's/_1.cat.fastq.gz/_1_cat_val_1.fq.gz/g' > trimmed_files.txt
```
Creating an file to append as the extra column as the sample name:
```sh
ls *_1* | sed 's/_1.q33.fastq.gz//g' | sed 's/_1_cat.gz/ /g' | sed 's/_1.cat.fastq.gz/ /g' > sample_names.txt
paste trimmed_files.txt sample_names.txt > trimmed_files_names.txt
rm trimmed_files.txt
rm temp.txt
```

### Running bwa

Running script called run_bwa_NEW.sh
```sh
# Runs bwa mem on the files specified in the "trimmed_files_names.txt".  Note: this is a tab delimited file.
# "awk is grabbing the 1st, 2nd, or 3rd column from the line above.
# For the while loop the "|| [ -n "$line" ]" Checks if the given string operand size is non-zero. If it is non-zero
# length then it returns true. We're using it to read the last line of the input file just in case it isn't
# terminated with a newline.

mkdir /media/santacruz5/luana_projects/kmar_ref/3.1_bwa_mapping
cd /media/santacruz5/luana_projects/kmar_ref/3.1_bwa_mapping

while read line || [ -n "$line" ];
do

        read1=$(echo $line | awk '{print $1}')
        read2=$(echo $line | awk '{print $2}')
        output=$(echo $line | awk '{print $3}')

        
        bwa mem \
        /media/santacruz5/luana_projects/kmar_ref/REFERENCES/as2_wsu/mitochondria/LHSH01_cat_MITO.fa /media/santacruz5/luana_projects/kmar_ref/2.1_trimgalore/$read1 /media/santacruz5/luana_projects/kmar_ref/2.1_trimgalore/$read2 > $output
        
done < /media/santacruz5/luana_projects/kmar_ref/1.1_fastq_rawData/1_ready_for_use_files/trimmed_files_names.txt 
```

Bwa mem produces sam files.

#### Converting sam files to bam files
```sh
cd /media/santacruz5/luana_projects/kmar_ref/3.2_bamfiles
```
Running script called samtobam.sh:
```sh
# This shell script calls an R script called sam_loc.R
ls *.sam > input_sam.txt

# this R script creates a txt file called "sam_locations.txt" that will be used to run samtools
Rscript sam_loc.R 

#then run samtools to convert SAM to BAM files

while read line || [ -n "$line" ];
do
        input=$(echo $line | awk '{print $1}')
        output=$(echo $line | awk '{print $2}')

	samtools view -S -h -b -o ../3.2_bamfiles/$output $input

done < sam_locations.txt
```

## Calling Variants 

Variants calling was done in the HPC computer Kamiak.
Files were copied to the following directory:
```sh
/data/kelley/projects/luana_projects/kmar_lineages/3.6_readGroup_info
```

Creating a index for the reference genome:
```sh
samtools faidx /data/kelley/projects/luana_projects/kmar_references/mitochondria_wsu/LHSH01_cat_MITO.fa
```

### HaplotypeCaller
From the folder ```/data/kelley/projects/luana_projects/kmar_lineages/4.1_haplotypeCaller ``
Run the create_runHapCal.sh script, which takes the master_HapCall.txt file and create one script for each sample
this script also runs the R script hap_rg_loc.R, which has to be in the same folder - this script will create the file "samples_HapCal.txt"
```sh
./create_runHapCal.sh
```
Master script for the HaplotypeCaller

master_HapCall.txt:
```sh
module load java

java -jar /data/kelley/projects/programs/GenomeAnalysisTK-3.7/GenomeAnalysisTK.jar \
-T HaplotypeCaller \
-R /data/kelley/projects/luana_projects/kmar_references/mitochondria_wsu/LHSH01_cat_MITO.fa \
-I /data/kelley/projects/luana_projects/kmar_lineages/3.6_readGroup_info/MYSAMPLE \
--genotyping_mode DISCOVERY \
--emitRefConfidence BP_RESOLUTION \
-stand_call_conf 30 \
-o /data/kelley/projects/luana_projects/kmar_lineages/4.1_haplotypeCaller/MYOUTPUT.g.vcf
```

### GenotypeGVCFs

Uses the output of the HaplotypeCaller and makes the joint call based on the data of all samples.
The flag "-allSites" emits all sites.

genoGVCF_gitmo_kmar_kamiak_NEW.sh:
```sh
module load java/oracle_1.8.0_92

java -jar /data/kelley/projects/programs/GenomeAnalysisTK-3.7/GenomeAnalysisTK.jar \
   -T GenotypeGVCFs \
   -R /data/kelley/projects/luana_projects/kmar_references/mitochondria_wsu/LHSH01_cat_MITO.fa \
   --variant /data/kelley/projects/luana_projects/kmar_lineages/4.1_haplotypeCaller/BBSC_ACAGTG.g.vcf \
   --variant /data/kelley/projects/luana_projects/kmar_lineages/4.1_haplotypeCaller/BWN3_TCGATG.g.vcf \
   --variant /data/kelley/projects/luana_projects/kmar_lineages/4.1_haplotypeCaller/DAN2K_GGTCTA.g.vcf \
   --variant /data/kelley/projects/luana_projects/kmar_lineages/4.1_haplotypeCaller/FDS1-10_TGGATC.g.vcf \
   --variant /data/kelley/projects/luana_projects/kmar_lineages/4.1_haplotypeCaller/FW2-4_CCAGAT.g.vcf \
   --variant /data/kelley/projects/luana_projects/kmar_lineages/4.1_haplotypeCaller/HON.g.vcf \
   --variant /data/kelley/projects/luana_projects/kmar_lineages/4.1_haplotypeCaller/LION2_TGCTAG.g.vcf \
   --variant /data/kelley/projects/luana_projects/kmar_lineages/4.1_haplotypeCaller/Lion_TGACCA.g.vcf \
   --variant /data/kelley/projects/luana_projects/kmar_lineages/4.1_haplotypeCaller/LK1_GCCAAT.g.vcf \
   --variant /data/kelley/projects/luana_projects/kmar_lineages/4.1_haplotypeCaller/PR_TAAGCC-cat.g.vcf \
   --variant /data/kelley/projects/luana_projects/kmar_lineages/4.1_haplotypeCaller/R2.g.vcf \
   --variant /data/kelley/projects/luana_projects/kmar_lineages/4.1_haplotypeCaller/RHL_A80AP1ABXX_L6_B80A2CABXX_L4.g.vcf \
   --variant /data/kelley/projects/luana_projects/kmar_lineages/4.1_haplotypeCaller/SLC8E_GTCAGT.g.vcf \
   --variant /data/kelley/projects/luana_projects/kmar_lineages/4.1_haplotypeCaller/Vol_CGATGT.g.vcf \
   --variant /data/kelley/projects/luana_projects/kmar_lineages/4.1_haplotypeCaller/Gitmo_TTAGGC.g.vcf \
	-allSites
   -o /data/kelley/projects/luana_projects/kmar_lineages/4.2_genotypeGVCFs/gitmo_kmar/Gitmo_kmar_genotypeGVCF_allsites.vcf
```

## Processing Variants

### Extracting only SNPS from the vcf file

Directory with vcf:
```sh
/data/kelley/projects/luana_projects/kmar_lineages/5.1_raw_SNPS/all_jointGeno/all_JointGeno_raw_snps.vcf
```
Running script
```sh
sbatch selectVar_kamiak_NEW_all.sh

selectVar_kamiak_NEW_all.sh:

mkdir /data/kelley/projects/luana_projects/kmar_lineages/5.1_raw_SNPS/all_jointGeno

module load java/oracle_1.8.0_92

java -jar /data/kelley/projects/programs/GenomeAnalysisTK-3.7/GenomeAnalysisTK.jar \
-T SelectVariants \
-R /data/kelley/projects/luana_projects/kmar_references/mitochondria_wsu/LHSH01_cat_MITO.fa \
-V /data/kelley/projects/luana_projects/kmar_lineages/4.2_genotypeGVCFs/gitmo_kmar/Gitmo_kmar_genotypeGVCF.vcf \
-selectType SNP \
-o /data/kelley/projects/luana_projects/kmar_lineages/5.1_raw_SNPS/all_jointGeno/all_JointGeno_raw_snps.vcf
```

### Filtering Variants - Gatk

This script applies the gatk filters to the set of snps and the uses vcftools to extract the pass-only SNPs from that vcf.
```sh
sbatch variantFiltration_kamiak_all.sh
```
variantFiltration_kamiak_all.sh:
```sh
module load java/oracle_1.8.0_92

java -jar /data/kelley/projects/programs/GenomeAnalysisTK-3.7/GenomeAnalysisTK.jar \
-T VariantFiltration \
-R /data/kelley/projects/luana_projects/kmar_references/mitochondria_wsu/LHSH01_cat_MITO.fa \
--variant //data/kelley/projects/luana_projects/kmar_lineages/5.1_raw_SNPS/all_jointGeno/all_JointGeno_raw_snps.vcf \
--filterExpression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
--filterName "snps_filter" \
-o /data/kelley/projects/luana_projects/kmar_lineages/5.2_filtered_SNPS/all_JointGeno_filtered_snps.vcf


# use vcftools to create files with pass-only SNPs

module load gcc

/data/kelley/projects/programs/vcftools-vcftools-2543f81/bin/vcftools --vcf /data/kelley/projects/luana_projects/kmar_lineages/5.2_filtered_SNPS/all_JointGeno_filtered_snps.vcf --remove-filtered-all --recode --recode-INFO-all  --out /data/kelley/projects/luana_projects/kmar_lineages/5.2_filtered_SNPS/1.1_pass_only_snps/all_JointGeno_pass_only_snps

```
### Filtering variants for missingness and LD pruning

Filtering variants for missingness and LD prunning for PCA analysis.

This script removes the missingness of the whole set of snps (pass only, with all individuals) and then it thins it for LD in windowns of 5kb it creates the ped and bed files needed for the PCA analysis.
```sh
sbatch vcftools8_kamiak_NEW_all.sh
```
vcftools8_kamiak_NEW_all.sh:
```sh
## this script removes the missingness of the whole set of snps (pass only, with all individuals)
# then it thins it for LD in windowns of 5kb 
# it creates the ped and bed files needed for the PCA analysis

module load gcc

# mkdir /data/kelley/projects/luana_projects/kmar_lineages/5.2_filtered_SNPS/2.1_missingness_all

# exclude sites with missingness higher than 85%
vcftools --vcf /data/kelley/projects/luana_projects/kmar_lineages/5.2_filtered_SNPS/1.1_pass_only_snps/all_JointGeno_pass_only_snps.recode.vcf \
--max-missing 0.85 --recode \
--out /data/kelley/projects/luana_projects/kmar_lineages/5.2_filtered_SNPS/2.1_missingness_all/all_JointGeno.miss

# mkdir /data/kelley/projects/luana_projects/kmar_lineages/5.2_filtered_SNPS/3.1_thinned_snps_5kb_all

###
mkdir /data/kelley/projects/luana_projects/kmar_lineages/5.2_filtered_SNPS/3.1_thinned_snps_5kb_all/1_ped_file/all

# Thin sites so that no two sites are within the specified distance (5000) from one another
vcftools --vcf /data/kelley/projects/luana_projects/kmar_lineages/5.2_filtered_SNPS/2.1_missingness_all/all_JointGeno.miss.recode.vcf --out /data/kelley/projects/luana_projects/kmar_lineages/5.2_filtered_SNPS/3.1_thinned_snps_5kb_all/all_JointGeno.miss.thin5kb --thin 5000 --recode --recode-INFO-all

# mkdir /data/kelley/projects/luana_projects/kmar_lineages/5.2_filtered_SNPS/3.1_thinned_snps_5kb_all/1_ped_file

#make a ped file
vcftools --vcf /data/kelley/projects/luana_projects/kmar_lineages/5.2_filtered_SNPS/3.1_thinned_snps_5kb_all/all_JointGeno.miss.thin5kb.recode.vcf --plink --out /data/kelley/projects/luana_projects/kmar_lineages/5.2_filtered_SNPS/3.1_thinned_snps_5kb_all/1_ped_file/all/all_JointGeno.miss.thin5kbplink

# recode the file
/data/cornejo/projects/programs/plink-1.07-x86_64/plink --file /data/kelley/projects/luana_projects/kmar_lineages/5.2_filtered_SNPS/3.1_thinned_snps_5kb_all/1_ped_file/all/all_JointGeno.miss.thin5kbplink --recode12 --noweb --out /data/kelley/projects/luana_projects/kmar_lineages/5.2_filtered_SNPS/3.1_thinned_snps_5kb_all/1_ped_file/all/all_JointGeno.miss.thin5kbplinkrecode12

mkdir /data/kelley/projects/luana_projects/kmar_lineages/5.2_filtered_SNPS/3.1_thinned_snps_5kb_all/2_bed_file/all

# make a bed file
/data/cornejo/projects/programs/plink-1.07-x86_64/plink --file /data/kelley/projects/luana_projects/kmar_lineages/5.2_filtered_SNPS/3.1_thinned_snps_5kb_all/1_ped_file/all/all_JointGeno.miss.thin5kbplinkrecode12 --make-bed --out /data/kelley/projects/luana_projects/kmar_lineages/5.2_filtered_SNPS/3.1_thinned_snps_5kb_all/2_bed_file/all/all_JointGeno.miss.thin5kbplinkrecode12 --noweb
```
# Downstream analyses

## PCA 

PCA in plink with all samples (including Gitmo).
```sh
sbatch pca_plink_kamiak_NEW_all.sh
```
pca_plink_kamiak_NEW_all.sh:
```sh
mkdir /data/kelley/projects/luana_projects/kmar_lineages/7_PCA_all

plink --bfile /data/kelley/projects/luana_projects/kmar_lineages/5.2_filtered_SNPS/3.1_thinned_snps_5kb_all/2_bed_file/all/all_JointGeno.miss.thin5kbplinkrecode12 --pca --out /data/kelley/projects/luana_projects/kmar_lineages/7_PCA_all/all_JointGeno.miss.thin5kbplinkrecode12
```
PCA in plink without Gitmo
```sh
sbatch pca_plink_kamiak_NEW_noGit.sh
```
pca_plink_kamiak_NEW_noGit.sh:
```sh
mkdir /data/kelley/projects/luana_projects/kmar_lineages/7_PCA_all_nogit

plink --bfile /data/kelley/projects/luana_projects/kmar_lineages/5.2_filtered_SNPS/3.2_thinned_snps_5kb_nogitmo/2_bed_file/all_JointGeno_snps_noGit.miss.thin5kbplinkrecode12 --pca --out /data/kelley/projects/luana_projects/kmar_lineages/7_PCA_all/all_nogitmo/all_JointGeno_noGit.miss.thin5kbplinkrecode12
```

## Long Runs of Homozygosity


The analysis of the long runs of homozygosity in vcftools uses a vcf as input.
The first step to analyse the vcf for LROH is to remove the singletons of the file.

#### Removing Singletons

Remove sigletons from vcf file.
```sh
/data/kelley/projects/luana_projects/kmar_lineages/5.2_filtered_SNPS/2.2_missingness_noGitmo/all_JointGeno_snps_noGitmo.miss.recode.vcf

mkdir /data/kelley/projects/luana_projects/kmar_lineages/5.2_filtered_SNPS/4.1_no_singletons_noGit
```

Get a list of singletons using vcf tools
```sh
vcftools --vcf /data/kelley/projects/luana_projects/kmar_lineages/5.2_filtered_SNPS/2.2_missingness_noGitmo/all_JointGeno_snps_noGitmo.miss.recode.vcf --singletons --out /data/kelley/projects/luana_projects/kmar_lineages/5.2_filtered_SNPS/4.1_no_singletons_noGit/all_JointGeno_noGit.miss
```
Use the singletons to remove those sites from vcf.
Inspect the output of the singletons to check if it has chrom and position

```sh
cd /data/kelley/projects/luana_projects/kmar_lineages/5.2_filtered_SNPS/4.1_no_singletons_noGit
```
Excluding the private doubletons from the list of singletons
```sh
awk '{ if ($3 == "S") { print } }' all_JointGeno_noGit.miss.singletons > all_JointGeno_noGit.miss.singletons.noDOUB
```
Making a list with only the Chrom and position
```sh
awk '{print $1, $2}' all_JointGeno_noGit.miss.singletons.noDOUB > singletons.noDOUB_list.txt
```
Use vcftools to exclude the position in the list (singletons)
```sh
vcftools --vcf /data/kelley/projects/luana_projects/kmar_lineages/5.2_filtered_SNPS/2.2_missingness_noGitmo/all_JointGeno_snps_noGitmo.miss.recode.vcf --exclude-positions /data/kelley/projects/luana_projects/kmar_lineages/5.2_filtered_SNPS/4.1_no_singletons_noGit/singletons.noDOUB_list.txt --recode --out /data/kelley/projects/luana_projects/kmar_lineages/5.2_filtered_SNPS/4.1_no_singletons_noGit/all_JointGeno_noGit.nosigletons.noDOUB 
```
Count how many snps are in the file that I excluded the singletons
```sh
grep -v "#" all_JointGeno_noGit.nosigletons.noDOUB.recode.vcf |wc -l
$ 1978617
```
Count how many lines the original vcf had
```sh
grep -v "#" /data/kelley/projects/luana_projects/kmar_lineages/5.2_filtered_SNPS/2.1_missingness_all/all_JointGeno.miss.recode.vcf | wc -l
$ 2095404
```
Count how many lines the singleton file has
```sh
wc -l  all_JointGeno.miss.singletons.noDOUB
$ 154836
```sh
Create new directories
```sh
mkdir /data/kelley/projects/luana_projects/kmar_lineages/8_ROH/noGit_LROH_singletons_noDOUB
mkdir /data/kelley/projects/luana_projects/kmar_lineages/8_ROH/noGit_LROH_singletons_noDOUB/LROH_output
```
#### Creating list of contigs to loop through vcftools

Creating a file with only the contigs in the vcf
```sh
mkdir 8_ROH
cd 8_ROH
grep "##contig=<ID=" /data/kelley/projects/luana_projects/kmar_lineages/5.2_filtered_SNPS/all_JointGeno_filtered_snps.vcf > contigs.txt

# checking 
head contigs.txt 
tail contigs.txt 
wc -l contigs.txt 
$ 40758
```

Cleaning the file by deleting Deleting the “##contig=<ID”and “,length=” out of the list, which will leave only the name of the contig in this file
 
```sh
 sed 's/##contig=<ID=//' contigs.txt | sed 's/,length=.*//' > contigs_clean.txt
```
 
Testing the command in one contig
```sh
module load gcc

vcftools --vcf /data/kelley/projects/luana_projects/kmar_lineages/5.2_filtered_SNPS/all_JointGeno_filtered_snps.vcf --LROH --chr gi_1037232374_gb_LHSH01000001.1_ --out kmar
```
It works. But this file is empty because probably there are not run of homozigosity in this chromosome. 


To speed up, I am separating the data in 5 blocks:

total number of contigs: 40758
40758/5 = 40,758 (because we want an even number I will separate in 5 even groups)
4 x 8151 + 8154

Making the  5 different contigs:
```sh
head -n 8151 contigs_clean.txt > contigs_clean_1.txt
head -n 16302 contigs_clean.txt| tail -n 8151 > contigs_clean_2.txt
tail -n 24456 contigs_clean.txt| head -n 8151 > contigs_clean_3.txt
tail -n 24456 contigs_clean.txt| tail -n 16305 | head -n 8151 > contigs_clean_4.txt
tail -n 8154 contigs_clean.txt > contigs_clean_5.txt
```

Modify the main script to split it in to 5 scripts that take this file, to create one loop script for each of the contigs list
just change the file name that it takes at the end of the script, and rename it.

Run lroh
```sh
sbatch noGit_1loop_lroh_kamiak_singleton_noDOUB.sh
sbatch noGit_2loop_lroh_kamiak_singleton_noDOUB.sh
sbatch noGit_3loop_lroh_kamiak_singleton_noDOUB.sh
sbatch noGit_4loop_lroh_kamiak_singleton_noDOUB.sh
sbatch noGit_5loop_lroh_kamiak_singleton_noDOUB.sh
```

Example of one of the scripts, noGit_1loop_lroh_kamiak_singleton_noDOUB.sh:
```sh
Output_Directory="/data/kelley/projects/luana_projects/kmar_lineages/8_ROH/noGit_LROH_singletons_noDOUB/LROH_output"
module load gcc

while read line || [ -n "$line" ];
do

        Contig=$(echo $line | awk '{print $1}')

/data/kelley/projects/programs/vcftools-vcftools-2543f81/bin/vcftools --vcf /data/kelley/projects/luana_projects/kmar_lineages/5.2_filtered_SNPS/4.1_no_singletons_noGit/all_JointGeno_noGit.nosigletons.noDOUB.recode.vcf --LROH --chr $Contig --out $Output_Directory/kmar.$Contig

done < /data/kelley/projects/luana_projects/kmar_lineages/8_ROH/contigs_clean_1.txt
```

#### Wrangling the output of the LROH

Directory: 
```sh
/data/kelley/projects/luana_projects/kmar_lineages/8_ROH/noGit_LROH_singletons_noDOUB/LROH_output
```
Counting how many files are in each of the outputs from LROH in vcftools
```sh
wc -l kmar* > lines_in_files.txt
```
Count how many lines are the "lines_in_files.txt". The correct number should have the number of contigs +1, the last line is the total of lines counted.
```sh
wc -l lines_in_files.txt  
$ 40759 #exactly the same that the one with doubletons
```
Get the lines that have values greater than 1 in the first column
```sh
awk '$1 > 1' lines_in_files.txt > long_runs.txt # this file has a list of contigs that have long runs of homozygosity
```
Counting how many contigs have long runs
```sh
wc -l long_runs.txt 
$ 6515      #(6589 long_runs.txt (the one with the doubletons had 6869) (there are fewer long runs than it had with the sigletons in the vcf: 8414)
```
Move only contigs that have long runs to a new folder
```sh
sbatch /data/kelley/projects/luana_projects/kmar_lineages/SCRIPTS/move_contigs_while.sh
```
move_contigs_while.sh:
```sh
# moves contigs based on the long_runs.txt file
# creates a new folder and moves files that match the names in the list there

mkdir ../longruns_Contigs

while read line || [ -n "$line" ];
do
        filename=$(echo $line | awk '{print $2}')
mv $filename ../longruns_Contigs/
done < long_runs.txt
```
Create a file containing all the contigs in one:

Concatenate files excluding the first line (the first line contains the header)
Create a file with all the contigs by excluding the first line of all the files "tail -n +2"
```sh
cd /data/kelley/projects/luana_projects/kmar_lineages/8_ROH/noGit_LROH_singletons_noDOUB/longruns_Contigs
head -1 "kmar.gi_1037191428_gb_LHSH01040755.1_.LROH" > all_longrun_contigs.txt
tail -n +2 -q kmar* >> all_longrun_contigs.txt
```
Count how many lines are in the file all_longrun_contigs.txt
```sh
wc -l all_longrun_contigs.txt
$ 19001 all_longrun_contigs.txt 
```
Counting the number of lines in all the files:
```sh
wc -l kmar*
total lines :
$ 25514 total  (26500 (with doubletons 27551 total))
```
Counting the number of files:
```sh
ll kmar* | wc -l 
$ 6514
# Number of lines minus the number of files (each file has a header): 19,912 (same that all lines in the file all_longrun_contigs.txt minus one (the header line).
```

Making a file for each individual

I will used this subset of the data to run int o a R script to calculate the number of basepairs in each of the classes of the ROH.

```sh
cd /data/kelley/projects/luana_projects/kmar_lineages/SCRIPTS
sbatch subset_lroh.sh
```
subset_lroh.sh:
```sh
head -n 1 all_longrun_contigs.txt > R2_roh.txt; grep "R2" all_longrun_contigs.txt >> R2_roh.txt
head -n 1 all_longrun_contigs.txt > Vol_roh.txt; grep "Vol" all_longrun_contigs.txt >> Vol_roh.txt
head -n 1 all_longrun_contigs.txt > RHL_roh.txt; grep "RHL" all_longrun_contigs.txt >> RHL_roh.txt
head -n 1 all_longrun_contigs.txt > FW2-4_roh.txt; grep "FW2-4" all_longrun_contigs.txt >> FW2-4_roh.txt
head -n 1 all_longrun_contigs.txt > Lion_roh.txt; grep "Lion" all_longrun_contigs.txt >> Lion_roh.txt
head -n 1 all_longrun_contigs.txt > Gitmo_roh.txt; grep "Gitmo" all_longrun_contigs.txt >> Gitmo_roh.txt
head -n 1 all_longrun_contigs.txt > LK1_roh.txt; grep "LK1" all_longrun_contigs.txt >> LK1_roh.txt
head -n 1 all_longrun_contigs.txt > SLC8E_roh.txt; grep "SLC8E" all_longrun_contigs.txt >> SLC8E_roh.txt
head -n 1 all_longrun_contigs.txt > BWN3_roh.txt; grep "BWN3" all_longrun_contigs.txt >> BWN3_roh.txt
head -n 1 all_longrun_contigs.txt > DAN2K_roh.txt; grep "DAN2K" all_longrun_contigs.txt >> DAN2K_roh.txt
head -n 1 all_longrun_contigs.txt > FDS1-10_roh.txt; grep "FDS1-10" all_longrun_contigs.txt >> FDS1-10_roh.txt
head -n 1 all_longrun_contigs.txt > HON_roh.txt; grep "HON" all_longrun_contigs.txt >> HON_roh.txt
head -n 1 all_longrun_contigs.txt > LION2_roh.txt; grep "LION2" all_longrun_contigs.txt >> LION2_roh.txt
head -n 1 all_longrun_contigs.txt > PR_roh.txt; grep "PR" all_longrun_contigs.txt >> PR_roh.txt
```

Make new folder and move contigs kmar there
```sh
mkdir contigs
mv kmar.gi* contigs/
```
#### run R script

Copy the file all_longrun_contigs.txt to local computer and run R script.

R script:

plot_LROH_noGIT_singletons_noDOUB.R:
```r
library(ggplot2)

setwd("/Users/Luana_Lins/Documents/WSU/kmar_project/kmar_lineages/8_ROH/noGit_LROH_singleton_noDOUB/")
mydata<-read.table("all_longrun_contigs.txt", header=TRUE)
mydata$max_length<- mydata[, "MAX_END"] - mydata[, "MIN_START"]

write.table(mydata,"all_longrun_contigs_legth.txt", quote=FALSE, col.names = TRUE)
###### Ploting the klength distribution before subsetting ######

p<-ggplot(mydata, aes(max_length, fill=INDV)) +
  geom_histogram(bins=300) + scale_x_log10(name="Length of runs of homozigosity", labels = scales::comma) +
  labs(title = "All data") +
  theme(plot.title = element_text(hjust = 0.5))
p

pdf("length_dist_ROH_all.pdf", 7,5)
p
dev.off()
# width = 7, height = 5, units =in
png("length_dist_ROH_all.png", width = 7, height = 5, units = "in", res = 300)
p 
dev.off()

 # check which is the first quartile and use that value as a cut for all the data
summary(mydata)

#1 st Qu.:   536
#subsetting the data 
 mydata<-mydata[mydata$max_length>=527,]
 summary(mydata)
 nrow(mydata)
 
#### Plot trimmed data ########
p2<- ggplot(mydata, aes(max_length, fill=INDV)) +
   geom_histogram(bins=200) + scale_x_log10(name="Length of runs of homozigosity (bp)", labels = scales::comma) +
   labs(title = "Trimmed at 541 (first quartile)") +
   theme(plot.title = element_text(hjust = 0.5))
p2

pdf("length_dist_ROH_trimmed.pdf", 7,5)
p2
dev.off()
# save png
png("length_dist_ROH_trimmed.png", width = 7, height = 5, units = "in", res = 300)
p2
dev.off()

summary(mydata)
 
 ########## Estimating Clusters in the data ########

data_length<- as.data.frame(mydata["max_length"])
mycluster<-kmeans(data_length, 3, nstart =20)

mycluster
#make a dataframe with the cluster values (1,2 or 3)
clustering <- as.data.frame(mycluster$cluster)

# merge cluster adta with the whole dataset by the row names
merged<-merge(mydata,clustering, by="row.names",all.x=TRUE)

 # subset based on the cluster group
small_class   <-merged[merged$`mycluster$cluster`==1,]
medium_class  <-merged[merged$`mycluster$cluster`==2,]
large_class   <-merged[merged$`mycluster$cluster`==3,]


paste("small", min(small_class$max_length),"to", max(small_class$max_length))
paste("medium", min(medium_class$max_length),"to", max(medium_class$max_length))
paste("large", min(large_class$max_length),"to", max(large_class$max_length))

# small 527 to 64625
# medium 65477 to 242872
# large 244849 to 974281


#summary(small_class$max_length)
#summary(medium_class$max_length)
#summary(large_class$max_length)

#### this is still a work in progress to plot the density of each cluster together####
#sd(small_class$max_length)
#sd(medium_class$max_length)
#sd(large_class$max_length)
#plot(density(small_class$max_length), col="red")
#par(new=TRUE)
#plot(density(medium_class$max_length), col="steelblue")
#par(new=TRUE)
#plot(density(large_class$max_length), col="green")


###################

ref<-read.table("../reference/LHSH01_cat_MITO.fa.fai", header =FALSE)

#count the size of the genome
sum(ref$V2)
# 653977110

#subset file to only inclue contigs that are longer than the longest class of long runs of homozygosity (class 3, minimum length 258048 bp)
ref_class3 <- ref[ref$V2>=244849,]

#write contig names to a file
write.table(ref_class3$V1, file="contigs_class3.txt", quote = FALSE, sep="\t", col.names = FALSE, row.names = FALSE)

#count number of basepairs that are in contigs longer than class 3
sum(ref_class3$V2)

#read file with all runs of homozygosity

roh<- read.table("all_longrun_contigs.txt", header =TRUE)

mydata<-merge(ref_class3, roh, by.x="V1", by.y="CHROM", all.x=TRUE, all.y=FALSE)

#removing rows with NA on column1 that has the contig name
mydata<-mydata[complete.cases(mydata[ ,"AUTO_START"]),]

mydata$max_length<- mydata[, "MAX_END"] - mydata[, "MIN_START"]
# mydata$auto_length <- mydata [, "AUTO_END"] - mydata [, "AUTO_START"]

sample_names<-(unique(levels(mydata$INDV)))


# create empty dataframe to recieve tha calculated values for each sample
roh_length  <- data.frame()
roh_proportion <- data.frame()
genome <- sum(ref_class3$V2)
#roh_length<- data.frame(Sample=NA,Non_long_run=NA,Class_1=NA,Class_2=NA,Class_3=NA)
######################### start of the loop block ##############################
for (sample in sample_names) {
  
  subset<- mydata[mydata$INDV==sample,]
  
  # define classes  for subset data
  non_roh<- subset[subset$max_length<=526, ]
  class1 <- subset[subset$max_length>526 & subset$max_length<=64625, ]
  class2 <- subset[subset$max_length>=65477 & subset$max_length<=242872, ]
  class3 <- subset[subset$max_length>=244849, ]
  # calculate the lengths of each class in each sample (subset)
  # len_nonlroh_mydata <- as.numeric(sum(non_roh$max_length))
  len_c1 <- as.numeric(sum(class1$max_length))
  len_c2 <- as.numeric(sum(class2$max_length))
  len_c3 <- as.numeric(sum(class3$max_length))
  len_nonlroh_estimated <- genome - (len_c1+len_c2+len_c3)
  # calculate the proportion of the genome covered by each class  
  prop_nonlroh  <- as.numeric((len_nonlroh_estimated/genome)*100)
  prop_c1       <- as.numeric((len_c1/genome)*100)
  prop_c2       <- as.numeric((len_c2/genome)*100)
  prop_c3       <- as.numeric((len_c3/genome)*100)
  
  temp_roh_length <- c(sample,len_nonlroh_estimated,len_c1,len_c2,len_c3)
  temp_roh_prop <- c(sample,prop_nonlroh,prop_c1,prop_c2, prop_c3)
  roh_length<-rbind(roh_length, data.frame(t(temp_roh_length),row.names=NULL))
  roh_proportion<-rbind(roh_proportion, data.frame(t(temp_roh_prop),row.names=NULL))
}

colnames(roh_length)<-c("Sample","Non_long_run", "Class_1","Class_2","Class_3")
colnames(roh_proportion)<- c("Sample","% Non_long_run","% Class_1","% Class_2","% Class_3")

write.table(roh_proportion, file="prop_contigs_class3.txt", quote = FALSE, sep="\t", col.names = TRUE, row.names = FALSE)

######## Ploting ###################
library(ggplot2)
library(reshape)
library(scales)

### transforming the values to numeric #### 
roh_length$Non_long_run<- as.numeric(as.character(roh_length$Non_long_run))
roh_length$Class_1<- as.numeric(as.character(roh_length$Class_1))
roh_length$Class_2<- as.numeric(as.character(roh_length$Class_2))
roh_length$Class_3<- as.numeric(as.character(roh_length$Class_3))


##### melt data to plot absolute values ####
roh_melt<-melt(roh_length, id="Sample")
roh_melt
roh_melt_non<- roh_melt[roh_melt$variable!="Non_long_run",]
class3_only<- roh_melt_non[roh_melt_non$variable=="Class_3",]
# reoder the level (names of the samples) by the subset of only class 3
roh_melt_non$Sample <- factor(roh_melt_non$Sample, levels= class3_only$Sample[order(class3_only$value)])


# reorder #
g <- ggplot(roh_melt_non, aes(x=Sample, y=value)) + scale_fill_brewer(palette = "Spectral") + 
  geom_bar(aes(fill=variable), col="black", stat="identity") +
  theme(axis.text.x=element_text(angle=90, hjust = 1))+
  scale_y_continuous(labels=comma) +
  labs(title="Runs of Homozygosity", y="Number of basepairs in the genome",fill="Classes of ROH",
       x=NULL) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
g
pdf("roh_counts_class3_contigs_noGIT_singletons_noDOUB.pdf", 7, 5)
g
dev.off()

######### Plotting proportion of the genome that is in runs of homozygosity ####

roh_proportion$'% Non_long_run' <- as.numeric(as.character(roh_proportion$'% Non_long_run'))
roh_proportion$'% Class_1'      <- as.numeric(as.character(roh_proportion$'% Class_1'))
roh_proportion$'% Class_2'      <- as.numeric(as.character(roh_proportion$'% Class_2'))
roh_proportion$'% Class_3'      <- as.numeric(as.character(roh_proportion$'% Class_3'))


roh_prop_melt<- melt(roh_proportion, id="Sample")
roh_prop_melt
# subset to not include the non long runs, they take over the plot
roh_prop_melt_non<- roh_prop_melt[roh_prop_melt$variable!='% Non_long_run',]


# subset by coverage sample with coverage < 5 will not be plotted 
# samples to be excluded = Vol, LK1, Gitmo, Lion, BBSC
roh_prop_subset<- roh_prop_melt_non[roh_prop_melt_non$Sample!="Vol" & roh_prop_melt_non$Sample!="LK1" & roh_prop_melt_non$Sample!="Gitmo" & roh_prop_melt_non$Sample!="Lion" & roh_prop_melt_non$Sample!="BBSC",]
class3_prop_only<- roh_prop_subset[roh_prop_subset$variable=="% Class_3",]
roh_prop_subset$Sample <- factor(roh_prop_subset$Sample, levels= class3_prop_only$Sample[order(class3_prop_only$value)])


g3<- ggplot(roh_prop_subset, aes(x=Sample, y=(value/100))) + scale_fill_brewer(palette = "Spectral") +
  geom_bar(aes(fill=variable), col="black", stat="identity") +
  theme(axis.text.x=element_text(angle=90, hjust = 1))+
  scale_y_continuous(labels = percent) +
  labs(title="Runs of Homozygosity", y="Percentage of the K. marmoratus genome",fill="Classes of ROH",
       x=NULL) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
g3

pdf("roh_percentage_coverage_class3_contigs_noGIT_sigletons_NODOUB.pdf", 7, 5)
g3
dev.off()

#######################

# plot correlation of % of coverage and number of generations

# change the name of the lineages to match the generations file
prop<-read.table("prop_contigs_class3.txt", header =TRUE)
gen<- read.table("../../kmar_generations.txt", header =TRUE)

gen_prop <- merge(prop,gen, by.x="Sample", by.y="lineage", all.x=TRUE, all.y=TRUE)
# remove rows with NA
gen_prop<-gen_prop[complete.cases(gen_prop[ ,]),] 
cove<-read.table("../lineages_high_coverage.txt", header =TRUE)
kmar2<-merge(cove,gen_prop, by.x="sample", by.y="Sample", all.x=TRUE, all.y=FALSE)


g2<-ggplot(kmar2, aes(x=generations, y=X._Class_3, colour=location)) +
  geom_point(size =6) +  
  geom_smooth(method=lm,se=FALSE, colour="grey") +
  labs(title="R2=0.67, p=0.025", y="Proportion of the Genome Covered by ROH Class 3 (%)", x="Number of Generations in Laboratory") +
  scale_x_continuous(breaks=seq(0,11,1)) +
  # scale_y_continuous(breaks=seq(0,190000,10000)) +
  geom_text(aes(label=sample),hjust=0.3, vjust=-1, size=5, colour ="black") +
  theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA, size=1), panel.grid.major = element_blank(), axis.text.x=element_text(size=15), axis.text.y=element_text(size=15),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

g2

kmar2.lm <- lm(kmar2$X._Class_3 ~ kmar2$generations, data=kmar2)
summary(kmar2.lm)$r.squared 
summary(kmar2.lm)$coefficients

pdf("roh_gen_class3_prop_cove_noGIT.pdf", 6.5, 5)
g2
dev.off()
```



## Relatedness 

Calculating relatedness for the Florida lineages

vcftools7_kamiak_FLORIDA.sh:
```sh
module load gcc

#remove gitmo for the relatedness analysis from the snps only data filtered by gatk
vcftools --vcf /data/kelley/projects/luana_projects/kmar_lineages/5.2_filtered_SNPS/1.1_pass_only_snps/kmar_pass_only_snps.recode.vcf --keep Florida.txt --recode --out /data/kelley/projects/luana_projects/kmar_lineages/5.2_filtered_SNPS/1.2_pass_only_snps_noGitmo/Florida_snps

mkdir /data/kelley/projects/luana_projects/kmar_lineages/5.2_filtered_SNPS/2.3_missingness_Florida

# exclude sites with missingness higher than 85%
vcftools --vcf /data/kelley/projects/luana_projects/kmar_lineages/5.2_filtered_SNPS/1.2_pass_only_snps_noGitmo/Florida_snps.recode.vcf \
--max-missing 0.85 --recode \
--out /data/kelley/projects/luana_projects/kmar_lineages/5.2_filtered_SNPS/2.3_missingness_Florida/Florida.miss

mkdir /data/kelley/projects/luana_projects/kmar_lineages/6_relatedness_Florida

#run relatedness
vcftools --vcf /data/kelley/projects/luana_projects/kmar_lineages/5.2_filtered_SNPS/2.3_missingness_Florida/Florida.miss.recode.vcf --relatedness --out /data/kelley/projects/luana_projects/kmar_lineages/6_relatedness_Florida/Florida
vcftools --vcf /data/kelley/projects/luana_projects/kmar_lineages/5.2_filtered_SNPS/2.3_missingness_Florida/Florida.miss.recode.vcf --relatedness2 --out /data/kelley/projects/luana_projects/kmar_lineages/6_relatedness_Florida/Florida
```

#### Plotting Relatedness in R
```r
library(plotrix)
library(reshape)
library(corrplot)

rela<- read.table("Florida.relatedness2", header=TRUE)
rela<-rela[c("INDV1", "INDV2", "RELATEDNESS_PHI")]

mdat <- melt(rela)
mcor<-cast(mdat, INDV1 ~ INDV2)
# individuals names in column one become the row names(instead of the autmatic numbering)
rownames(mcor) <- mcor[,1]
#delete column one from the data, because now the data ther is in the row names
mcor[,1] <- NULL

#write the matrix as a text file
write.table(mcor, file="relatedness2_matrix.txt", quote = FALSE, sep="\t")

#plotting the matrix
#install.packages("corrplot")
# look here for more details on corroplot:
# https://cran.r-project.org/web/packages/corrplot/vignettes/corrplot-intro.html

mcor<- as.matrix(mcor)
pdf("relatedness2.pdf", 7,5)
corrplot(mcor,is.corr=FALSE, type="lower", order ="hclust", addrect=2, tl.srt=0, tl.col="gray20", tl.offset = 0.6, tl.cex = .8,
         addCoef.col = "black", number.digits = 2,number.cex = 0.5)
dev.off()
```
