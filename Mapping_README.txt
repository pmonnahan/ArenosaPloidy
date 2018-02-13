#### NGS GENOME SCAN PIPELINE YANT LAB 2016
###########################################

####### OVERVIEW ######
The aim of the pipeline is to produce a master VCF for population genetic analysis from raw fastq.gz files from sequencing.
The python3 scripts generate bash shell scripts for a batch of input files and send them to the NBI SLURM cluster for computation.

All steps require loading of python3
source python-3.3.3

We have several steps:

In processing order
P1 to P2	Preparation steps
M1 to M4	Mapping steps
V1 to V3	Variant calling steps
F1 to F5	Filtering steps



###### PREPARATION STEPS: P1, P2 ######

## P1_concat_fastqgz.py
Start with this concatenation script if the samples were sequenced on several flow cell lanes. Otherwise skip ahead to step P2.
Use the -print true option to verify the bash shell scripts before submitting them to the cluster.

synthax:
batch_concat_fastqgz.py -ind	<REQUIRED Relative or full path toist of individual fastq.gz files>
						-s		<REQUIRED Search: directory from which to start to look for folder containing the fastq.gz files>
						-o		<REQUIRED Output: directory in which to save the concatenated files>
						-lanes	<REQUIRED Number of lanes on which the individual was sequenced>
						-print	<If 'true' then bash shell are printed to the screen and NOT submitted to the cluster>
						-mem	<Requested memory>
							
IMPORTANT NOTE
The root directory for the search has to be below the output directory in the file tree.
-ind has to be a list of the complete name of the fastq.gz file (eg. <sample_index1-index2.Rx.fastq.gz>) with a single file per line.
The script can be run from a folder which is above the input containing folders in the hierarchy.


## P2_cutadapt_fastqgz.py
Start with this step if your samples have been sequenced on a single lane.

# Requires existing output directory

synthax:
python3 P2_cutadapt_fastqgz.py	-list	<REQUIRED: Relative or full path to text file with input fastq or fastq.gz files>
								-type	<REQUIRED: Type of adapters used in library preparation, either 'NextSeq' or 'TruSeq'>
								-mem	<Requested memorey [10 000]>
								-time	<Requested job time [0-2:00]>
								-print	<If 'true' bash shell scripts are printed to the screen and NOT submitted to the cluster>

IMPORTANT NOTE: -list has to be in the following format:
input_filename_R1 input_filename_R2 output_filename_R1 output_filename_R2

The columns have to be separated by space.
The output files have to be fastq and a NOT fastq.gz.



###### MAPPING STEPS: M1, M2, M3, M4 ######

## M1_bwa_mem.py
Mapping the short reads to a reference genome using bwa mem. Usually the input files are the 
output files of step P2_cutadapt_fastqgz.py. The output directory is created automatically, if is does not
already exist. Default values and directory paths are set for all options. Use the -print true option to verify
the bash shell scripts before submitting them to the cluster.
The tool has been found to run most efficiently on 2 cores with 8 GB RAM (parallelisation tested from 1 to 8 cores).

Run from directory above -fastqdir for easiest handling.

synthax <explanation [default value]>:
python3 M1_bwa_mem.py 	-fastqdir	<REQUIRED: Relative or full path to input file directory containing input fastq files [fastq_ca]>
						-R			<REQUIRED: Full path to reference fasta file [/nbi/Research-Groups/JIC/Levi-Yant/Lyrata_ref/alygenomes.fasta]>
						-o			<Relative or full path to output directory [aligned/]>
						-c			<Number of cores for parallelisation [2]>
						-mem		<Requested memory [8000]>
						-time		<Requested computation time [1-00:00]>
						-print		<If 'true' the script prints the bash shell scripts to screen and does NOT submit them to the cluster>
		
IMPORTANT NOTE
Be aware that the job is cancelled if the -time argument value has been reached. If in doubt, request more time.
Ignore slurmstepd errors in the output file if the job finishes with exitcode=0.0. These are warning messages created by
samtools which want to reserve 32GB RAM irrespective of actual need.


## M2_deduped.py
Remove duplicated reads (mainly from PCR steps). The input files are the sorted bam output files of step M1_bwa_mem.py. 
The output directory is created automatically, if is does not already exist. Default values and directory paths are set 
for all options. Use the -print true option to verifythe bash shell scripts before submitting them.

Run from directory above -aligndir for easiest handling.

synthax <explanation [default value]>:
python3 M2_deduped.py 	-aligndir	<REQUIRED: Relative or full path to input file directory containing input fastq files [aligned/]>
						-o			<Relative or full path to output directory [deduped/]>
						-c			<Number of cores for parallelisation [5]>
						-mem		<Requested memory [16000]>
						-time		<Requested computation time [0-1:00]>
						-print		<If 'true' the script prints the bash shell scripts to screen and does NOT submit them to the cluster>

IMPORTANT NOTE
Be aware that the job is cancelled if the -time argument value has been reached.


## M3_namefix.py
Add read groups, column headers and create indices of the bam files. The input files are the output bam files of step M2_deduped.py. 
The output directory is created automatically, if is does not already exist. Default values and directory paths are set for all options.
Use the -print true option to verify the bash shell scripts before submitting them.

Run from directory above -dedupdir for easiest handling.

synthax <explanation [default value]>:
python3 M3_namefix.py 	-dedupdir	<REQUIRED: Relative or full path to input file directory containing input dedup_bam files [deduped/]>
						-o			<Relative or full path to output directory [namefixed/]>
						-c			<Number of cores for parallelisation [5]>
						-time		<Requested computation time [0-1:00]>
						-print		<If 'true' the script prints the bash shell scripts to screen and does NOT submit them to the cluster>

IMPORTANT NOTE
Be aware that the job is cancelled if the -time argument value has been reached.


## M4_realign.py
Mapping indel in realation to the reference sequence and realignment of the bam files.
The output directory is created automatically, if is does not already exist. Default values and directory paths are set for all options.
Use the -print true option to verify the bash shell scripts before submitting them.

synthax <explanation [default value]>:
python3 M4_realidn.py 	-namefixdir	<REQUIRED: Relative or full path to input file directory containing input dedup_bam files [namefixed/]>
						-o			<Relative or full path to output directory [namefixed/]>
						-R			<REQUIRED: Full path to reference fasta file [/nbi/Research-Groups/JIC/Levi-Yant/Lyrata_ref/alygenomes.fasta]>
						-c			<Number of cores for parallelisation [12]>
						-mem		<Requested memory [24000]>
						-time		<Requested computation time [0-12:00]>
						-print		<If 'true' the script prints the bash shell scripts to screen and does NOT submit them to the cluster>

IMPORTANT NOTE
Be aware that the job is cancelled if the -time argument value has been reached.



###### VARIANT CALLING STEPS: V1, V2, V3 & V01, V02 ######

## V1_HC.PY
Call variants via HaplotypeCaller using analysis ready reads. Cohort1 and Cohort2 should be used to separate 
individuals of different ploidy. If all individuals have the same ploidy, enter only the cohort1 options only.
Input files are the output files of step M4_realign.py. The output directory is created automatically, if is does not 
already exist. Use the -print true option to verify the bash shell scripts before submitting them. The output files are bgzipped.
Default values are set for all options except -cohort1_ind.
The tool has been found to run most efficiently on 4 cores and 16 GB RAM total (parallelisation tested from 1 to 8 cores).

V1_HC.py creates haplotypes for scaffolds 1-8 within a single script by using the GATK -L option. If you want to split the 
analysis per scaffold to speed up, use V01_HC_scaffold.py. V01_HC_scaffold.py creates an output foler per individual 
in which the scaffold_x.g.vcf.gz files are saved.

Run from directory above -realigndir for easiest handling.

synthax <explanation [default value]>:
python3 V1_HC.py -realigndir		<REQUIRED: Relative or full path to input file directory containing input realigned bam files [realigned/]>
				 -R					<REQUIRED: Full path to reference fasta file [/nbi/Research-Groups/JIC/Levi-Yant/Lyrata_ref/alygenomes.fasta]>
				 -cohort1_ind		<REQUIRED: Relative or full path to txt file containing the indiviudal names, one name per line>
				 -cohort1_ploidy	<REQUIRED: Ploidy of individuals in cohort1>	
				 -cohort2_ind		<Relative or full path to txt file containing the individuals names of cohort2, one name per line>
				 -cohort2_ploidy	<Ploidy of individuals in cohort2
				 -minbaseq			<Minimum base quality to be used in analysis [25]>
				 -minmapq			<Minimum mapping quality to be used in analysis [25]>
				 -o					<Relative or full path to output directory [HC/]>
				 -c					<Number of cores for parallelisation [5]>
				 -mem				<Requested memory [16000]>
				 -time				<Requested computation time [6-00:00]>
				 -print				<If 'true' the script prints the bash shell scripts to screen and does NOT submit them to the cluster>

IMPORTANT NOTE
Be aware that the job is cancelled if the -time argument value has been reached. If in doubt, request more time.
Use versions V0x if you split the analysis per scaffold.


## V2_combine_GVCF.py
If you have more than ca. 200 input files, creating multisample GVCFs is necessary before genotyping. If you have < 200 files 
skip ahead to step V3_genotypeGVCF.py.

Run from directory above -HCdir for easiest handling.

synthax <explanation [default value]>:
python3 V2_combineGVCF.py	-HCdir	<REQUIRED: Relative or full path to directory with input g.vcf.gz files [HC/]>
							-scaf	<REQUIRED: Full path to file with scaffold list, one scaffold per line [/nbi/Research-Groups/JIC/Levi-Yant/Lyrata_ref/alygenomes.8scaffolds]>
							-R		<REQUIRED: Full path to reference fasta file [/nbi/Research-Groups/JIC/Levi-Yant/Lyrata_ref/alygenomes.8scaffolds]>
							-o		<Relative or full path to output directory [CombGVCF/]>
							-mem	<Total memory requested for each job (MB) [20000]>
							-time	<Time for job [2-00:00]>
							-print	<If 'true' the script prints the bash shell scripts to screen and does NOT submit dem to the cluster [false]>

IMPORTANT NOTE
Be aware that the job is cancelled if the -time argument value has been reached. If in doubt, request more time.
You will likely have to manually adjust the samples you want to combine. If it is too many (> 200) it will fail.

If you used version V01_HC_scaffold.py, use V02_combine_scaffolds.py. After V02 you can run V2 to create a multisample GVCF.


## V3_genotypeGVCF.PY
Genotype variants called by HaplotypeCaller. The input files are the output files of either V1, V2, or V02.
The output directory is created automatically, if is does not already exist. Use the -print true option to 
verify the bash shell scripts before submitting them. Default values have been set for all options except -i.
The tool has been found to run most efficiently on 5 cores and 20 GB RAM total (parallelisation tested from 1 to 8 cores).

Run from directory above -HCdir for easiest handling.

synthax <explanation [default value]>:
python3 V3_genotypeGVCF.py 	-HCdir	<REQUIRED: Relative or full path to input file directory containing input g.vcf.gz files [HC/]>
							-scaf	<REQUIRED: Relative or full path to list of scaffolds [/nbi/Research-Groups/JIC/Levi-Yant/Lyrata_ref/alygenomes.8scaffolds]>
							-R		<REQUIRED: Full path to reference fasta file [/nbi/Research-Groups/JIC/Levi-Yant/Lyrata_ref/alygenomes.fasta]>
							-o		<Relative or full path to output directory [Genotyped/]>
							-c		<Number of cores for parallelisation [5]>
							-mem	<Requested memory [20000]>
							-time	<Requested computation time [3-00:00]>
							-print	<If 'true' the script prints the bash shell scripts to screen and does NOT submit them to the cluster>
		

IMPORTANT NOTE
Be aware that the job is cancelled if the -time argument value has been reached. If in doubt, request more time.
The memory requests for this step can be huge (200GB). Be prepared to adjust this resource request a few times.



###### FILTERING STEPS: F1, F2, F3, F4, F5 ######

## F1_BI.py
In this step biallelic SNPs that are present in all samples are selected. Input files are the 
output of step V3_genotypeGVCF. The output files are given at suffix. The output folder is 
created automatically, if it does not exist. Use the -print true option to verify the bash 
shell scripts before submitting them. The -nChr gives the number of genotyped chromosomes 
present in the dataset, eg 10 tetraploid individuals have -nChr 40 (10 individuals * 4 copies of the genome).

synthax <explanation [default value]>:
python3 F1_BI.py	-i		<REQUIRED: Relative or full path to input file directory containing input vcf files>
					-o		<REQUIRED: Relative or full path to output directory>
					-R		<REQUIRED: Full path to reference fasta file [/nbi/Research-Groups/JIC/Levi-Yant/Lyrata_ref/alygenomes.fasta]>
					-nChr	<REQUIRED: Threshold number of sampled chromosomes>
					-c		<Number of cores for parallelisation [6]>
					-suffix	<Suffix appended to output filename [_chr_BI]>
					-mem	<Requested memory [8000]>
					-time	<Requested computation time [0-2:00]>
					-print	<If 'true' the script prints the bash shell scripts to screen and does NOT submit them to the cluster>


## F2_BestPractise.py
This scripts flags the variants with 'BP' that fail either of the set filters that follow GATK best practises
recommendations. Input files are the output of step F1_BI and the output files are given a suffix.
The output folder is created automatically, if it does not exist. Use the -print true option 
to verify the bash shell scripts before submitting them.

synthax <explanation [default value]>:
python3 F2_BestPractise.py	-i		<REQUIRED: Relative or full path to input file directory containing input vcf files>
							-o		<REQUIRED: Relative or full path to output directory>
							-R		<REQUIRED: Full path to reference fasta file [/nbi/Research-Groups/JIC/Levi-Yant/Lyrata_ref/alygenomes.fasta]>
							-suffix	<Suffix appended to output filename [_BP]>
							-mem	<Requested memory [6000]>
							-time	<Requested computation time [0-2:00]>
							-print	<If 'true' the script prints the bash shell scripts to screen and does NOT submit them to the cluster>


## F3_HetMask.py
This script flags sites in a vcf file that exhibit excess heterozygosity
The input files are the output from step F2_BP plus a list of sites exhibiting excess heterozygosity. The output folder is created automatically, 
if it does not exist. Use the -print true option to verify the bash shell scripts before submitting them.



## F5_excluedFiltered.py
This script removes every variant that does not have 'PASS' in the filter column.
Use the -print true option to verify the bash shell scripts before submitting them.

synthax <explanation [default value]>:
python3 F5_excludeFiltered.py	-i		<REQUIRED: Relative or full path to input file directory containing input vcf files>
								-o		<REQUIRED: Relative or full path to output directory>
								-R		<REQUIRED: Full path to reference fasta file [/nbi/Research-Groups/JIC/Levi-Yant/Lyrata_ref/alygenomes.fasta]>
								-suffix	<Suffix appended to output filename [_PASS]>
								-c		<Number of requested cores [6]>
								-mem	<Requested memory [12000]>
								-time	<Requested computation time [0-4:00]>
								-print	<If 'true' the script prints the bash shell scripts to screen and does NOT submit them to the cluster>

## ReadDepthFilter.py
This script removes sites that exhibit excess read depth.



