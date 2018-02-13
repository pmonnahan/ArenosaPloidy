#!/usr/bin/env python3
# written in perl by Levi Yant, streamlined by Jeff DaCosta,  by Christian Sailer, 2 May 2016

import os, sys, argparse, subprocess

print()

#create variables that can be entered as arguments in command line
parser = argparse.ArgumentParser(description='This script uses the analysis ready reads for variant calling via the GATK tool '+
                                 'HaplotypeCaller. The cohort1 and cohort2 option should be used to separate ploidies. A seperate '+
                                 'SLURM shell script is created for each input file and sent to the NBI SLURM cluster. If the output '+
                                 'directory does not exist, it is created automatically.')

parser.add_argument('-realigndir', type=str, metavar='', default='realigned/', help='REQUIRED: Full path to directory with input realigned bam files')
parser.add_argument('-R', type=str, metavar='reference_path', default='/nbi/Research-Groups/JIC/Levi-Yant/Lyrata_ref/alygenomes.fasta', help='REQUIRED: Full path to reference genome')
parser.add_argument('-cohort1_ind', type=str, metavar='cohort1_individuals', required=True, help='REQUIRED: Full path to text file containing the list of individuals of cohort 1. One individual per line.')
parser.add_argument('-cohort1_ploidy', type=int, metavar='cohort1_ploidy', required=True, help='REQUIRED: Ploidy of cohort 1 (integer)')                               
parser.add_argument('-cohort2_ind', type=str, metavar='cohort2_individuals', default='na', help='Full path to text file containing the list of individuals of cohort 2. One individual per line. [na]')
parser.add_argument('-cohort2_ploidy', type=int, metavar='cohort2_ploidy', help='Ploidy of cohort 2 (integer)')
parser.add_argument('-minbaseq', type=int, metavar='min_base_qual', default='25', help='Minimum base quality to be used in analysis [25]')
parser.add_argument('-minmapq', type=int, metavar='min_mapping_qual', default='25', help='Minimum mapping quality to be used in analysis [25]')
parser.add_argument('-suffix', type=str, metavar='', default="HC", help='Suffix added to end of output filename [HC]')
parser.add_argument('-o', type=str, metavar='HC_dir', default='HC/', help='Realtive path to output directory [HC/]')
parser.add_argument('-c', type=int, metavar='num_cores', default='4', help='Number of requested cores for job [4]')
parser.add_argument('-mem', type=str, metavar='memory', default='16000', help='Total memory for each job (Mb) [16000]')
parser.add_argument('-time', type=str, metavar='time', default='6-00:00', help='Time for job [6-00:00]')
parser.add_argument('-print', type=str, metavar='print', default='false', help='If changed to true then shell files are printed to screen and not launched [false]')
args = parser.parse_args()

# print what ploidies were entered
print('\ncohort1_ploidy = '+str(args.cohort1_ploidy)+'\n' +
      'cohort2_ploidy = '+str(args.cohort2_ploidy)+'\n')

count = 0

# test if two cohorts and 2 different ploidies are set
if args.cohort1_ploidy is args.cohort2_ploidy or args.cohort2_ind is 'na':
    if args.cohort1_ploidy is args.cohort2_ploidy:
        print('\nBoth cohorts have the same ploidy, so all samples will be genotyped together\n')
    else:
        print('\nOnly one cohort defined. All samples in cohort1 will be genotyped together\n')

    # gather list of input bam files
    in_bam_list = []
    cohort1 = open(args.cohort1_ind, 'r')
    for ind in cohort1:
        data = ind
        data = ind.replace('\n', '')
        in_bam_list.append(data)
        in_bam_list.sort()

    print('\nFound '+str(len(in_bam_list))+' input realigned bam files')
    for in_bam in in_bam_list:
        print('\t'+in_bam)
    print('\n\nCreating shell files and sending jobs to NBI SLURM cluster\n\n')

    #check if output directory exists, create it if necessary
    if os.path.exists(args.o) == False:
        os.mkdir(args.o)

    #loop through input bam files
    count = 0
    for in_bam in in_bam_list:
        bam_basename = in_bam

        sh_file = open(args.o+bam_basename+'.sh','w')
        #write slurm shell file
        sh_file.write('#!/bin/bash -e\n'+
                      '#SBATCH -J HC.'+bam_basename+'\n'+
                      '#SBATCH -o '+args.o+bam_basename+'_'+args.suffix+'.out\n'+ 
                      '#SBATCH -e '+args.o+bam_basename+'_'+args.suffix+'.err\n'+
                      '#SBATCH -p nbi-long\n'+
                      '#SBATCH -c '+str(args.c)+'\n'+
                      '#SBATCH -t '+args.time+'\n'
                      '#SBATCH --mem='+args.mem+'\n'+
                      #'source jre-1.8.0_45\n'+
                      'source GATK-3.5.0\n'+
                      'srun java XX:ParallelGCThreads=2 -Xmx16g -jar /nbi/software/testing/GATK/3.5.0/src/GenomeAnalysisTK.jar -T '+
                      'HaplotypeCaller -I '+args.realigndir+bam_basename+'_realigned.bam --min_base_quality_score '+str(args.minbaseq)+
                      ' --min_mapping_quality_score '+str(args.minmapq)+' -rf DuplicateRead -rf BadMate -rf BadCigar '+
                      '-ERC BP_Resolution -R '+args.R+' -L scaffold_1 -L scaffold_2 -L scaffold_3 -L scaffold_4 -L scaffold_5 -L '+
                      'scaffold_6 -L scaffold_7 -L scaffold_8 -o '+args.o+bam_basename+'_'+args.suffix+'_g.vcf.gz -ploidy '+
                      str(args.cohort1_ploidy)+' -stand_emit_conf 13 -stand_call_conf 25 --pcr_indel_model NONE -nct '+str(args.c)+'\n'+
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


else:
    print('\nCohorts have different ploidy. Samples will be analysed with different ploidy flags.\n')

    #check if output directory exists, create it if necessary
    if os.path.exists(args.o) == False:
        os.mkdir(args.o)

    # select samples for cohort 1
    cohort1_list = []
    cohort1 = open(args.cohort1_ind, 'r')
    for ind in cohort1:
        data = ind
        data = ind.replace('\n', '')
        cohort1_list.append(data)
        cohort1_list.sort()
        
    # select samples for cohort 2
    cohort2_list = []
    cohort2 = open(args.cohort2_ind, 'r')
    for ind in cohort2:
        data = ind
        data = ind.replace('\n', '')
        cohort2_list.append(data)
        cohort2_list.sort()
    print('\nFound '+str(len(cohort1_list))+' individuals in cohort '+str(args.cohort1_ind)+'.'+
          '\nFound '+str(len(cohort2_list))+' individuals in cohort '+str(args.cohort2_ind)+'.\n')

    # print individiduals of cohort1
    print('\nIndividuals of cohort 1')
    for ind in cohort1_list:
        print('\t'+ind)
    # print individiduals of cohort2
    print('\nIndividuals of cohort 2')
    for ind in cohort2_list:
        print('\t'+ind)

    count = 0

    # loop through bam files of cohort 1
    for ind in cohort1_list:
        bam_basename = ind
        sh_file = open(args.o+bam_basename+'.sh','w')
        #write slurm shell file
        sh_file.write('#!/bin/bash -e\n'+
                      '#SBATCH -J HC.'+bam_basename+'\n'+
                      '#SBATCH -o '+args.o+bam_basename+'_'+args.suffix+'.out\n'+ 
                      '#SBATCH -e '+args.o+bam_basename+'_'+args.suffix+'.err\n'+
                      '#SBATCH -p nbi-long\n'+
                      '#SBATCH -c '+str(args.c)+'\n'+
                      '#SBATCH -t '+args.time+'\n'
                      '#SBATCH --mem='+args.mem+'\n'+
                      'source jre-1.8.0_45\n'+
                      'source GATK-3.6.0\n'+
                      'srun java XX:ParallelGCThreads=2 -Xmx16g -jar /nbi/software/testing/GATK/3.6.0/src/GenomeAnalysisTK.jar -T '+
                      'HaplotypeCaller -I '+args.realigndir+bam_basename+'_realigned.bam --min_base_quality_score '+str(args.minbaseq)+
                      ' --min_mapping_quality_score '+str(args.minmapq)+' -rf DuplicateRead -rf BadMate -rf BadCigar '+
                      '-ERC BP_Resolution -R '+args.R+' -L scaffold_1 -L scaffold_2 -L scaffold_3 -L scaffold_4 -L scaffold_5 -L '+
                      'scaffold_6 -L scaffold_7 -L scaffold_8 -o '+args.o+bam_basename+'_'+args.suffix+'_g.vcf.gz -ploidy '+
                      str(args.cohort1_ploidy)+' -stand_emit_conf 13 -stand_call_conf 25 --pcr_indel_model NONE -nct '+str(args.c)+'\n'+
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

    # loop through bam files of cohort 2
    for ind in cohort2_list:
        bam_basename = ind

        sh_file = open(args.o+bam_basename+'.sh','w')
        #write slurm shell file
        sh_file.write('#!/bin/bash -e\n'+
                      '#SBATCH -J HC.'+bam_basename+'\n'+
                      '#SBATCH -o '+args.o+bam_basename+'_'+args.suffix+'.out\n'+ 
                      '#SBATCH -e '+args.o+bam_basename+'_'+args.suffix+'.err\n'+
                      '#SBATCH -p nbi-long\n'+
                      '#SBATCH -c '+str(args.c)+'\n'+
                      '#SBATCH -t '+args.time+'\n'
                      '#SBATCH --mem='+args.mem+'\n'+
                      'source jre-1.8.0_45\n'+
                      'source GATK-3.6.0\n'+
                      'srun java XX:ParallelGCThreads=2 -Xmx16g -jar /nbi/software/testing/GATK/3.6.0/src/GenomeAnalysisTK.jar -T '+
                      'HaplotypeCaller -I '+args.realigndir+bam_basename+'_realigned.bam --min_base_quality_score '+str(args.minbaseq)+
                      ' --min_mapping_quality_score '+str(args.minmapq)+' -rf DuplicateRead -rf BadMate -rf BadCigar '+
                      '-ERC BP_Resolution -R '+args.R+' -L scaffold_1 -L scaffold_2 -L scaffold_3 -L scaffold_4 -L scaffold_5 -L '+
                      'scaffold_6 -L scaffold_7 -L scaffold_8 -o '+args.o+bam_basename+'_'+args.suffix+'_g.vcf.gz -ploidy '+
                      str(args.cohort2_ploidy)+' -stand_emit_conf 13 -stand_call_conf 25 --pcr_indel_model NONE -nct '+str(args.c)+'\n'+
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
