#!/usr/bin/env python3
# written in perl by Levi Yant, streamlined by Jeff DaCosta, translated to python3 and adjusted to latest GATK best practises
# by Christian Sailer, 30 May 2016

import os, sys, argparse, subprocess

print()

#create variables that can be entered as arguments in command line
parser = argparse.ArgumentParser(description='This script uses the deduped bam files (output of step M2) '+
                                 'and adds columnheades and indices. '+
                                 'It creates one SLURM shell scripts per input file and '+
                                 'sends it to the NBI SLUM cluster. If the output '+
                                 'directory does not exist, it is created automatically.')

parser.add_argument('-dedupdir', type=str, metavar='dedup_dir', default='deduped/', help='REQUIRED: Full path to directory with input sort.bam files [deduped]')
parser.add_argument('-R', type=str, metavar='reference_path', default='/nbi/Research-Groups/JIC/Levi-Yant/Lyrata_ref/alygenomes.fasta', help='REQUIRED: Full path to reference genome')
parser.add_argument('-o', type=str, metavar='namefixed_dir', default='namefixed/', help='Realtive path to output directory [namefixed/]')
parser.add_argument('-c', type=int, metavar='num_cores', default='5', help='Number of requested cores for job [5]')
parser.add_argument('-mem', type=str, metavar='memory', default='5000', help='Total memory for each job (Mb) [5000]')
parser.add_argument('-time', type=str, metavar='time', default='0-1:00', help='Time for job [0-1:00]')
parser.add_argument('-print', type=str, metavar='print', default='false', help='If changed to true then shell files are printed to screen and not launched [false]')
args = parser.parse_args()


# gather list of input sort.bam files
in_bam_list = []
for file in os.listdir(args.dedupdir):
    if file.endswith('_dedup.bam'):
        in_bam_list.append(file)
in_bam_list.sort()

print('\nFound '+str(len(in_bam_list))+' input deduped bam files')
for in_bam in in_bam_list:
    print('\t'+in_bam)
print('\n\nCreating shell files and sending jobs to NBI SLURM cluster\n\n')

# check if output directory exists, create it if necessary
if os.path.exists(args.o) == False:
    os.mkdir(args.o)

#loop through input sort.bam files
count = 0
for bam in in_bam_list:
    bam_basename = bam.replace('_dedup.bam', '')
    
    #write slurm shell file
    sh_file = open(args.o+bam_basename+'.sh','w')
    sh_file.write('#!/bin/bash -e\n'+
                  '#SBATCH -J M3.'+bam_basename+'\n'+
                  '#SBATCH -o '+args.o+bam_basename+'.out\n'+ 
                  '#SBATCH -e '+args.o+bam_basename+'.err\n'+
                  '#SBATCH -p nbi-medium\n'+
                  '#SBATCH -c '+str(args.c)+'\n'+
                  '#SBATCH -t '+args.time+'\n'
                  '#SBATCH --mem='+args.mem+'\n'+
                  'source picard-1.134\n'+
                  'srun java -XX:ParallelGCThreads=2 -Xmx5g -jar /nbi/software/testing/picard/1.134/src/picard-tools-1.134/picard.jar '+
                  'AddOrReplaceReadGroups I='+args.dedupdir+bam_basename+'_dedup.bam O='+args.o+bam_basename+'_dedupNamed.bam ' +
                  'SORT_ORDER=coordinate CREATE_INDEX=True VALIDATION_STRINGENCY=LENIENT RGLB='+bam_basename+' RGPL=illumina RGPU=AAAAAA '+
                  'RGSM='+bam_basename+'\n' +
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
