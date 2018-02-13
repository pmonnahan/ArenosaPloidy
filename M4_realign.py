#!/usr/bin/env python3
# written in perl by Levi Yant, streamlined by Jeff DaCosta, translated to python3 and adjusted to latest GATK best practises
# by Christian Sailer, 30 May 2016

import os, sys, argparse, subprocess

print()

#create variables that can be entered as arguments in command line
parser = argparse.ArgumentParser(description='This script uses the namefixed bam files (output of step M3) '+
                                 'and realignes the reads to the reference '+
                                 '(mapping indels). It creates one SLURM shell scripts per input file and '+
                                 'sends it to the NBI SLUM cluster. If the output '+
                                 'directory does not exist, it is created automatically.')

parser.add_argument('-namefixdir', type=str, metavar='namefix_dir', default='namefixed/', help='REQUIRED: Full path to directory with input sort.bam files [namefixed]')
parser.add_argument('-R', type=str, metavar='reference_path', default='/nbi/Research-Groups/JIC/Levi-Yant/Lyrata_ref/alygenomes.fasta', help='REQUIRED: Full path to reference genome')
parser.add_argument('-o', type=str, metavar='realigned_dir', default='realigned/', help='Realtive path to output directory [realigned/]')
parser.add_argument('-c', type=int, metavar='num_cores', default='5', help='Number of requested cores for job [12]')
parser.add_argument('-mem', type=str, metavar='memory', default='24000', help='Total memory for each job (Mb) [24000]')
parser.add_argument('-time', type=str, metavar='time', default='0-12:00', help='Time for job [0-12:00]')
parser.add_argument('-print', type=str, metavar='print', default='false', help='If changed to true then shell files are printed to screen and not launched [false]')
args = parser.parse_args()


# gather list of input sort.bam files
in_bam_list = []
for file in os.listdir(args.namefixdir):
    if file.endswith('_dedupNamed.bam'):
        in_bam_list.append(file)
in_bam_list.sort()

print('\nFound '+str(len(in_bam_list))+' input namefixed bam files')
for in_bam in in_bam_list:
    print('\t'+in_bam)
print('\n\nCreating shell files and sending jobs to NBI SLURM cluster\n\n')

# check if output directory exists, create it if necessary
if os.path.exists(args.o) == False:
    os.mkdir(args.o)

#loop through input sort.bam files
count = 0
for bam in in_bam_list:
    bam_basename = bam.replace('_dedupNamed.bam', '')
    
    #write slurm shell file
    sh_file = open(args.o+bam_basename+'.sh','w')
    sh_file.write('#!/bin/bash -e\n'+
                  '#SBATCH -J M4.'+bam_basename+'\n'+
                  '#SBATCH -o '+args.o+bam_basename+'.out\n'+ 
                  '#SBATCH -e '+args.o+bam_basename+'.err\n'+
                  '#SBATCH -p nbi-medium\n'+
                  '#SBATCH -c '+str(args.c)+'\n'+
                  '#SBATCH -t '+args.time+'\n'
                  '#SBATCH --mem='+args.mem+'\n'+
                  'source jre-1.8.0_45\n'+
                  'source GATK-3.6.0\n'+
                  'srun java -XX:ParallelGCThreads=2 -Xmx24g -jar /nbi/software/testing/GATK/3.6.0/src/GenomeAnalysisTK.jar -nt '+
                  str(args.c)+' -T RealignerTargetCreator -R '+args.R+' -o '+args.o+bam_basename+'.IndelRealigner.intervals '+
                  '-I '+args.namefixdir+bam_basename+'_dedupNamed.bam\n'+
                  'srun java -Xmx8g -jar /nbi/software/testing/GATK/3.6.0/src/GenomeAnalysisTK.jar -nt 1 -T IndelRealigner '+
                  '-targetIntervals '+args.o+bam_basename+'.IndelRealigner.intervals -I '+args.namefixdir+bam_basename+'_dedupNamed.bam '+
                  '-R '+args.R+' -o '+args.o+bam_basename+'_realigned.bam\n'+
                  'printf "\\nFinished\\n\\n"\n')
    sh_file.close()

    #check if slurm shell file should be printed or sent to NBI SLURM 
    if args.print == 'false':
        #send slurm job to NBI SLURM cluster
        cmd = ('sbatch '+args.o+bam_basename+'.sh')
        p = subprocess.Popen(cmd, shell=True)
        sts = os.waitpid(p.pid, 0)[1]
    else:
        file = open(args.o+bam_basename+'.sh','r')
        data = file.read()
        print(data)
    count +=1

#if appropriate, report how many slurm shell files were sent to NBI SLURM 
if args.print == 'false':
    print('\nSent '+str(count)+' jobs to the NBI SLURM  cluster\n\nFinished!!\n\n')
