#!/usr/bin/env python3
# written by Christian Sailer, 31 May 2016

import os, sys, argparse, subprocess

print()

#create variables that can be entered as arguments in command line
parser = argparse.ArgumentParser(description='This script uses the g.vcf files created by HaplotypeCaller (V1_HC.py) and combines the '+
                                 'individual files into a single multisample GVCF. Adjust the number of files to combine manually. '+
                                 'A seperate SLURM shell script is created for each scaffold and sent to '+
                                 'the NBI SLURM cluster. If the output directory does not exist, it is created automatically.')

parser.add_argument('-HCdir', type=str, metavar='', default='HC/', help='REQUIRED: Full path to directory with input g.vcf.gz files')
parser.add_argument('-scaf', type=str, metavar='scaffold_path', default='/nbi/Research-Groups/JIC/Levi-Yant/Lyrata_ref/alygenomes.8scaffolds', help='REQUIRED: Full path to file with scaffold list')
parser.add_argument('-R', type=str, metavar='reference_path', default='/nbi/Research-Groups/JIC/Levi-Yant/Lyrata_ref/alygenomes.fasta', help='REQUIRED: Full path to reference genome')
parser.add_argument('-o', type=str, metavar='combGVCF_dir', default='CombGVCF/', help='Relative path to output directory [CombGVCF/]')
parser.add_argument('-mem', type=str, metavar='memory', default='20000', help='Total memory for each job (Mb) [20000]')
parser.add_argument('-time', type=str, metavar='time', default='2-00:00', help='Time for job [2-00:00]')
parser.add_argument('-print', type=str, metavar='print', default='true', help='If changed to true then shell files are printed to screen and not launched [true]')
args = parser.parse_args()

#check if output directory exists, create it if necessary
if os.path.exists(args.o) == False:
    os.mkdir(args.o)

in_gvcf_list = []
# gather list of input gvcf files by looping through the subdirList
count = 0
for file in os.listdir(args.HCdir):
    if file.endswith('.vcf'):
        in_gvcf_list.append(file)
in_gvcf_list.sort()

print('\nFound '+str(len(in_gvcf_list))+' input vcf files')
for in_gvcf in in_gvcf_list:
    print('\t'+in_gvcf)
print('\n\nCreating shell files and sending jobs to NBI SLURM cluster\n\n')

# combine all individual files
basename = "combined"
print(basename)
sh_file = open(args.o+basename+'.sh', 'w')
# write slurm shell file
sh_file.write('#!/bin/bash -e\n'+
              '#SBATCH -J c.'+basename+'\n'+
              '#SBATCH -o '+args.o+basename+'.out\n'+ 
              '#SBATCH -e '+args.o+basename+'.err\n'+
              '#SBATCH -p nbi-medium\n'+
              '#SBATCH -t '+args.time+'\n'
              '#SBATCH --mem='+args.mem+'\n'+
              'source jre-1.8.0_45\n'+
              'source GATK-3.6.0\n'+
              'srun java -XX:ParallelGCThreads=2 -Xmx20g -jar /nbi/software/testing/GATK/3.6.0/src/GenomeAnalysisTK.jar -T '+
              'CombineGVCFs -R '+args.R+' -o '+args.o+basename+'_vcf ')
for in_gvcf in in_gvcf_list:
    sh_file.write('-V '+args.HCdir+in_gvcf+' ')
sh_file.write('\nprintf "\\nFinished\\n\\n"\n')         
sh_file.close()

#check if slurm shell file should be printed or sent to NBI SLURM 
if args.print == 'false':
    #send slurm job to NBI SLURM cluster
    cmd = ('sbatch '+args.o+basename+'.sh')
    p = subprocess.Popen(cmd, shell=True)
    sts = os.waitpid(p.pid, 0)[1]
else:
    file = open(args.o+basename+'.sh','r')
    data = file.read()
    print(data)

count += 1

#if appropriate, report how many slurm shell files were sent to NBI SLURM 
if args.print == 'false':
    print('\nSent '+str(count)+' jobs to the NBI SLURM  cluster\n\nFinished!!\n\n')
