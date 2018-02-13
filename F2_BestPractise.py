#!/usr/bin/env python3
# written by Jeff DaCosta, adjusted for NBI SLURM cluster by Christian Sailer, 7 April 2016
# updated 6 June 2016

import os, sys, argparse, subprocess

print()

#create variables that can be entered as arguments in command line
parser = argparse.ArgumentParser(description='This script uses the GATK VariantFiltration function to mask variants based on an expression '+
'that follows best practices guidelines. Variants that fail the expression get a "BP" in the FILTER column. '+
'A separate slurm shell script is created and sent to Odyssey for each input vcf file. Output vcf files '+
'are given a suffix and written to the specified output directory, which is automatically created if it does not exist.')

parser.add_argument('-i', type=str, metavar='inputdir_path', required=True, help='REQUIRED: Full path to directory with input vcf files')
parser.add_argument('-o', type=str, metavar='outputdir_path', required=True, help='REQUIRED: Full path to directory for output vcf files')
parser.add_argument('-R', type=str, metavar='reference_path', default='/nbi/Research-Groups/JIC/Levi-Yant/Lyrata_ref/alygenomes.fasta', help='REQUIRED: Full path to reference genome')
parser.add_argument('-prefix', type=str, metavar='suffix', default='BP_', help='Suffix to append to output filename [_BP]')
parser.add_argument('-mem', type=str, metavar='memory', default='6000', help='Total memory for each job (Mb) [6000]')
parser.add_argument('-time', type=str, metavar='time', default='2-00:00', help='Time for job [2-00:00]')
parser.add_argument('-print', type=str, metavar='print', default='false', help='If changed to true then shell files are printed to screen and not launched [false]')
args = parser.parse_args()

#gather list of input vcf files
in_vcf_list = []
for file in os.listdir(args.i):
    if file.endswith('.vcf.gz'):
        in_vcf_list.append(file)
in_vcf_list.sort()

print('\nFound '+str(len(in_vcf_list))+' input vcf files:')
for in_vcf in in_vcf_list:
    print('\t'+in_vcf)

print('\nCreating shell files and sending jobs to NBI SLURM cluster\n\n')

#check if output directory exists, create it if necessary
if os.path.exists(args.o) == False:
    os.mkdir(args.o)

#loop through input vcf files
count = 0
for in_vcf in in_vcf_list:
    basename = in_vcf.replace('.vcf.gz','')
    sh_file = open(args.o+args.prefix+basename+'.sh','w')


    #write slurm shell file
    sh_file.write('#!/bin/bash -e\n'+
                  '#SBATCH -J '+basename+'\n'+
                  '#SBATCH -o '+args.o+args.prefix+basename+'.out\n'+
                  '#SBATCH -e '+args.o+args.prefix+basename+'.err\n'+
                  '#SBATCH -p nbi-long\n'+
                  '#SBATCH -n 1\n'+
                  '#SBATCH -t '+args.time+'\n'
                  '#SBATCH --mem='+args.mem+'\n'+
                  'source  GATK-nightly.2016.09.26\n'+
                  'srun java -XX:ParallelGCThreads=2 -Xmx6g -jar /nbi/software/testing/GATK/nightly.2016.09.26/x86_64/jars/GenomeAnalysisTK.jar -T VariantFiltration -V '+
                  args.i+in_vcf+' -R '+args.R+' -filterName "QD" -filter "QD<2.0" -filterName "FS" -filter "FS>60.0" -filterName "MQ" -filter "MQ<40.0" -filterName "MQRS" -filter "MQRankSum < -12.5" '+
                  '-filterName "RPRS" -filter "ReadPosRankSum < -8.0" -filterName "HS" -filter "HaplotypeScore<13.0" -o '+args.o+args.prefix+basename+'.vcf.gz\n'+
                  'printf "\\nFinished\\n\\n"\n')
    sh_file.close()

    #check if slurm shell file should be printed or sent to NBI SLURM
    if args.print == 'false':
    	#send slurm job to NBI SLURM cluster
    	cmd = ('sbatch '+args.o+args.prefix+basename+'.sh')
    	p = subprocess.Popen(cmd, shell=True)
    	sts = os.waitpid(p.pid, 0)[1]
    else:
        file = open(args.o+args.prefix+basename+'.sh','r')
        data = file.read()
        print(data)

    count += 1

#if appropriate, report how many slurm shell files were sent to NBI SLURM
if args.print == 'false':
    print('\nSent '+str(count)+' jobs to the NBI SLURM cluster\n\nFinished!!\n\n')
