#!/usr/bin/env python3
# written in perl by Levi Yant, streamlined by Jeff DaCosta, translated to python3 and adjusted to lates GATK bests practises
# by Christian Sailer, 20 May 2016

import os, sys, argparse, subprocess

print()

#create variables that can be entered as arguments in command line
parser = argparse.ArgumentParser(description='This script uses the fastq files that have the adapters trimmed off (output of step P02). '+
                                 'The script creates one SLURM shell scripts to map the short reads to a reference using bwa mem per input'+
                                 ' file and sends it to the NBI SLUM cluster. If the output directory does not exist, it is created '+
                                 'automatically.')

parser.add_argument('-fastqdir', type=str, metavar='fastq_ca', default='fastq_ca/', help='REQUIRED: Full path to directory with input fastq files')
parser.add_argument('-R', type=str, metavar='reference_path', default='/nbi/Research-Groups/JIC/Levi-Yant/Lyrata_ref/alygenomes.fasta', help='REQUIRED: Full path to reference genome')
parser.add_argument('-o', type=str, metavar='aligned_dir', default='aligned/', help='Realtive path to output directory [aligned/]')
parser.add_argument('-c', type=int, metavar='num_cores', default='5', help='Number of requested cores for job [5]')
parser.add_argument('-mem', type=str, metavar='memory', default='16000', help='Total memory for each job (Mb) [16000]')
parser.add_argument('-time', type=str, metavar='time', default='2-00:00', help='Time for job [2-00:00]')
parser.add_argument('-print', type=str, metavar='print', default='false', help='If changed to true then shell files are printed to screen and not launched [false]')
args = parser.parse_args()


# gather list of input fastq files
in_fastq_list = []
for file in os.listdir(args.fastqdir):
    if file.endswith('.fastq'):
        in_fastq_list.append(file)
in_fastq_list.sort()

print('\nFound '+str(len(in_fastq_list))+' input fastq files')
for in_fastq in in_fastq_list:
    print('\t'+in_fastq)
print('\n\nCreating shell files and sending jobs to NBI SLURM cluster\n\n')

# check if output directory exists, create it if necessary
if os.path.exists(args.o) == False:
    os.mkdir(args.o)

# select only R1 to be used for the loop
R1_list = []
for file in os.listdir(args.fastqdir):
    if file.endswith('_R1.fastq'):
        R1_list.append(file)
R1_list.sort()

#loop through input fastq files
count = 0
for fq1 in R1_list:
    fastq_basename = fq1.replace('_R1.fastq', '')
    fq2 = fq1.replace('R1', 'R2')
    
    #write slurm shell file
    sh_file = open(args.o+fastq_basename+'.sh','w')
    sh_file.write('#!/bin/bash -e\n'+
                  '#SBATCH -J M1.'+fastq_basename+'\n'+
                  '#SBATCH -o '+args.o+fastq_basename+'.out\n'+ 
                  '#SBATCH -e '+args.o+fastq_basename+'.err\n'+
                  '#SBATCH -p nbi-long\n'+
                  '#SBATCH -c '+str(args.c)+'\n'+
                  '#SBATCH -t '+args.time+'\n'
                  '#SBATCH --mem='+args.mem+'\n'+
                  'source bwa-0.7.12\n'+
                  'source samtools-1.3\n'+
                  'srun bwa mem -t '+str(args.c)+' '+args.R+' '+args.fastqdir+fq1+' '+args.fastqdir+fq2+' > '+args.o+fastq_basename+'.sam'+'\n'+
                  'srun samtools view -bS '+args.o+fastq_basename+'.sam | samtools sort -T '+fastq_basename+' > '+
                  args.o+fastq_basename+'.sort.bam\n'+
                  'srun samtools index '+args.o+fastq_basename+'.sort.bam\n'+
                  'srun rm '+args.o+fastq_basename+'.sam\n'+
                  'srun samtools idxstats '+args.o+fastq_basename+'.sort.bam > '+args.o+fastq_basename+'.idxstats\n'+
                  'srun samtools flagstat '+args.o+fastq_basename+'.sort.bam > '+args.o+fastq_basename+'.flagstat\n'+
                  'printf "\\nFinished\\n\\n"\n')
    sh_file.close()

    #check if slurm shell file should be printed or sent to NBI SLURM 
    if args.print == 'false':
        #send slurm job to NBI SLURM cluster
        cmd = ('sbatch '+args.o+fastq_basename+'.sh')
        p = subprocess.Popen(cmd, shell=True)
        sts = os.waitpid(p.pid, 0)[1]
    else:
        file = open(args.o+fastq_basename+'.sh','r')
        data = file.read()
        print(data)
    count +=1
    
#if appropriate, report how many slurm shell files were sent to NBI SLURM 
if args.print == 'false':
    print('\nSent '+str(count)+' jobs to the NBI SLURM  cluster\n\nFinished!!\n\n')
