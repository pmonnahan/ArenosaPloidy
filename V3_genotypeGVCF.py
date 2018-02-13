#!/usr/bin/env python3
# written by Christian Sailer, 23 May 2016

import os, sys, argparse, subprocess

print()

#create variables that can be entered as arguments in command line
parser = argparse.ArgumentParser(description='This script uses the g.vcf files created by HaplotypeCaller (V1_HC.py) and assigns genotypes '+
                                 'via the GATK tool GenotypeGVCFs. A seperate SLURM shell script is created for each scaffold and sent to '+
                                 'the NBI SLURM cluster. If the output directory does not exist, it is created automatically.')

parser.add_argument('-HCdir', type=str, metavar='', default='CombGVCF/', help='REQUIRED: Full path to directory with input g.vcf.gz files [CombGVCF]')
parser.add_argument('-scaf', type=str, metavar='scaffold_path', default='/nbi/Research-Groups/JIC/Levi-Yant/Lyrata_ref/alygenomes.8scaffolds', help='REQUIRED: Full path to directory with scaffold list')
parser.add_argument('-R', type=str, metavar='reference_path', default='/nbi/Research-Groups/JIC/Levi-Yant/Lyrata_ref/alygenomes.fasta', help='REQUIRED: Full path to reference genome')
parser.add_argument('-o', type=str, metavar='genotype_dir', default='Genotyped/', help='Realtive path to output directory [Genotyped/]')
parser.add_argument('-c', type=int, metavar='num_cores', default='5', help='Number of requested cores for job [5]')
parser.add_argument('-mem', type=str, metavar='memory', default='20000', help='Total memory for each job (Mb) [20000]')
parser.add_argument('-time', type=str, metavar='time', default='1-00:00', help='Time for job [3-00:00]')
parser.add_argument('-print', type=str, metavar='print', default='false', help='If changed to true then shell files are printed to screen and not launched [false]')
args = parser.parse_args()

# print the number of scaffolds specified
scaffold_list = []
scaff = open(args.scaf, 'r')
for line in scaff:
    data = line
    data = data.replace('\n', '')
    scaffold_list.append(data)
for scaf in scaffold_list:
    print('\t' + scaf)
print('\nFound '+str(len(scaffold_list))+' scaffold names for reference genome\n')

# gather list of input bam files
in_gvcf_list = []
for file in os.listdir(args.HCdir):
    if file.endswith('g.vcf.gz'):
        in_gvcf_list.append(file)
in_gvcf_list.sort()

print('\nFound '+str(len(in_gvcf_list))+' input g.vcf.gz files')
for in_gvcf in in_gvcf_list:
    print('\t'+in_gvcf)
print('\n\nCreating shell files and sending jobs to NBI SLURM cluster\n\n')

#check if output directory exists, create it if necessary
if os.path.exists(args.o) == False:
    os.mkdir(args.o)

count = 0
for scaf in scaffold_list:
    sh_file = open(args.o+scaf+'.sh','w')
    scaf_short =  scaf.replace('fold', '')
    #write slurm shell file
    sh_file.write('#!/bin/bash -e\n'+
                  '#SBATCH -J GT.'+scaf_short+'\n'+
                  '#SBATCH -o '+args.o+scaf_short+'.out\n'+ 
                  '#SBATCH -e '+args.o+scaf_short+'.err\n'+
                  '#SBATCH -p nbi-medium\n'+
                  '#SBATCH -c '+str(args.c)+'\n'+
                  '#SBATCH -t '+args.time+'\n'
                  '#SBATCH --mem='+args.mem+'\n'+
                  #'source jre-1.8.0_45\n'+
                  'source GATK-nightly.2016.09.26\n'+
                  'srun java -XX:ParallelGCThreads=2 -Xmx32g -jar /nbi/software/testing/GATK/nightly.2016.09.26/x86_64/jars/GenomeAnalysisTK.jar -T GenotypeGVCFs -R '+
                  args.R+' -L '+scaf+' -o '+args.o+scaf_short+'.vcf --includeNonVariantSites -nt '+str(args.c)+' ')
    for in_gvcf in in_gvcf_list:
        sh_file.write('-V '+args.HCdir+in_gvcf+' ')
    sh_file.write('\n'+'printf "\\nFinished\\n\\n"\n')
         
    sh_file.close()

    #check if slurm shell file should be printed or sent to NBI SLURM 
    if args.print == 'false':
        #send slurm job to NBI SLURM cluster
        cmd = ('sbatch '+args.o+scaf+'.sh')
        p = subprocess.Popen(cmd, shell=True)
        sts = os.waitpid(p.pid, 0)[1]
    else:
        file = open(args.o+scaf+'.sh','r')
        data = file.read()
        print(data)
        
    count +=1
#if appropriate, report how many slurm shell files were sent to NBI SLURM 
if args.print == 'false':
    print('\nSent '+str(count)+' jobs to the NBI SLURM  cluster\n\nFinished!!\n\n')
