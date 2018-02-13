#!/usr/bin/env python3
# written by Jeff DaCosta, adjusted for NBI SLURM cluster and added the loop by Christian Sailer, 16 June 2016
# modified for excess heterozygostiy masking of the 300, 14 October 2016

import os, sys, argparse, subprocess

print()

#create variables that can be entered as arguments in command line
parser = argparse.ArgumentParser(description='This script uses the GATK VariantFiltration function to filter variants in a given vcf file. '+
'Masked variants are retained in the vcf file, but a mask name is added to the FILTER column. These variants can later be removed with the GATK_VariantFiltration_excludeFiltered.py script. '+
' This script needs to be run separately for each vcf input/mask file.')

parser.add_argument('-i', type=str, metavar='inputdir_path', required=True, help='REQUIRED: Full path to directory with input vcf file (step F2)')
parser.add_argument('-o', type=str, metavar='outputdir_path', required=True, help='REQUIRED: Full path to directory with output vcf file')
parser.add_argument('-R', type=str, metavar='reference_path', default='/nbi/Research-Groups/JIC/Levi-Yant/Lyrata_ref/alygenomes.fasta', help='REQUIRED: Full path to reference genome [/nbi/Research-Groups/JIC/Levi-Yant/Lyrata_ref/alygenomes.fasta]')
parser.add_argument('-mn', type=str, metavar='mask_name', default='ExHet', help='Name of mask to be entered in Filter column [ExHet]')
parser.add_argument('-m', type=str, metavar='mask_file', default='Excess.Het.List.Genes.With.5orMore.2PopsMin.AllHet.Sites.txt', help='Full or relative path to bedfile containing the excess het sites [Excess.Het.List.Genes.With.5orMore.2PopsMin.AllHet.Sites.txt]')
parser.add_argument('-prefix', type=str, metavar='prefix', default='HM_', help='Prefix to append to output filename [HM_]')
parser.add_argument('-mem', type=str, metavar='memory', default='12000', help='Total memory for each job (Mb) [12000]')
parser.add_argument('-time', type=str, metavar='time', default='0-2:00', help='Time for job [0-2:00]')
parser.add_argument('-print', type=str, metavar='print', default='true', help='If changed to false then shell files are launched [true]')
args = parser.parse_args()

#gather list of input vcf files
in_vcf_list = []
for file in os.listdir(args.i):
    if file.endswith('.vcf.gz'):
        in_vcf_list.append(file)
in_vcf_list.sort()

print('\nFound '+str(len(in_vcf_list))+' input vcf files')
for in_vcf in in_vcf_list:
    print('\t'+in_vcf)
print('Creating shell files and sending jobs to NBI SLURM cluster\n\n')

#check if output directory exists, create it if necessary
if os.path.exists(args.o) == False:
    os.mkdir(args.o)

#loop through input vcf files
count = 0
for in_vcf in in_vcf_list:
    basename = in_vcf.replace('.vcf.gz','')
    sh_file = open(args.o+basename+'.sh','w')
    sh_file.write('#!/bin/bash\n'+
                  '#SBATCH -J M.'+basename+'\n'+
                  '#SBATCH -o '+args.o+basename+'.out\n'+
                  '#SBATCH -e '+args.o+basename+'.err\n'+
                  '#SBATCH -p nbi-short\n'+
                  '#SBATCH -t '+args.time+'\n'
                  '#SBATCH --mem='+args.mem+'\n'+
                  'source GATK-nightly.2016.09.26\n'+
                  'srun java -XX:ParallelGCThreads=2 -Xmx12g -jar /nbi/software/testing/GATK/nightly.2016.09.26/x86_64/jars/GenomeAnalysisTK.jar -T '+
                  'VariantFiltration -V '+args.i+in_vcf+' -R '+args.R+
                  ' --mask '+args.m+' --maskName "'+args.mn+'" -o '+args.o+args.prefix+basename+'.vcf\n'+
                  'printf "\\nFinished\\n\\n"\n')
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
    count +=1

#if appropriate, report how many slurm shell files were sent to NBI SLURM
if args.print == 'false':
    print('\nSent '+str(count)+' jobs to the NBI SLURM cluster\n\nFinished!!\n\n')
