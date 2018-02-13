#!/usr/bin/env python3

## For usage and help: ./batch_cutadapt1.8_PE.py -h
# updated July 29th by Christian Sailer

import os, sys, argparse, subprocess

#create variables that can be entered as arguments in command line
parser = argparse.ArgumentParser()
parser.add_argument('-list', type=str, metavar='', required=True, help='REQUIRED: Full path to text file with input fastq filename (column1) and output fastq filename (column2). Ensure that output directory exists.')
parser.add_argument('-type', type=str, metavar='', required=True, help='REQUIRED: Type of adapters: "Nextera" or "TruSeq"')
parser.add_argument('-mem', type=str, metavar='', default='10000', help='Total memory for each job (Mb) [10000]')
parser.add_argument('-time', type=str, metavar='', default='0-2:00', help='Time for job [0-2:00]')
parser.add_argument('-print', type=str, metavar='print', default='false', help='If changed to true then shell files are printed to screen and not launched [false]')
args = parser.parse_args()

#read input list file, populate lists of input and output file names
in1_list = []
in2_list = []                    
out1_list = []
out2_list = []
listfile = open(args.list,'r')
for line in listfile:
    data = line.split()
    in1_list.append(data[0])
    in2_list.append(data[1])
    out1_list.append(data[2])
    out2_list.append(data[3])
listfile.close()

print('\nFound '+str(len(in1_list))+' paired fastq files to process in cutadapter')
print('\nProcessing:')

count = 0

# adapter seqeunces
R1_Nex_adapt = 'CTGTCTCTTATACACATCTCCGAGCCCACGAGAC'
R2_Nex_adapt = 'CTGTCTCTTATACACATCTGACGCTGCCGACGA'
R1_TS_adapt = 'GATCGGAAGAGCACACGTCTGAACTCCAGTCAC'
R2_TS_adapt = 'AGATCGGAAGAGCGTCGTGTAGGGAAAGAG'

#create slurm shell file for each input file
for i in range(len(in1_list)):
    sh = open('ca.sh','w')
    print(in1_list[i]+' and '+in2_list[i])
    baseoutname = out1_list[i].replace('.R1.fastq','')
    sh.write('#!/bin/bash\n'+
             '#SBATCH -J CA.'+in1_list[i]+'\n'+
             '#SBATCH -e '+baseoutname+'.err\n'+
             '#SBATCH -o '+baseoutname+'.stats\n'+
             '#SBATCH -p nbi-short\n'+
             '#SBATCH -t '+args.time+'\n'+
             '#SBATCH --mem='+args.mem+'\n'+
             'source cutadapt-1.9.1\n'+
             'cutadapt -e 0.15 -O 4 -m 25 -a ')

    if args.type == 'Nextera':
        sh.write(R1_Nex_adapt)
    elif args.type == 'TruSeq':
        sh.write(R1_TS_adapt)
    else:
        print('ERROR: type must be either "Nextera" or "TruSeq"\n')
        quit()

    sh.write(' -A ')
    
    if args.type == 'Nextera':
        sh.write(R2_Nex_adapt)
    elif args.type == 'TruSeq':
        sh.write(R2_TS_adapt)
    else:
        print('ERROR: type must be either "Nextera" or "TruSeq"\n')
        quit()


    sh.write(' -o '+out1_list[i]+' -p '+out2_list[i]+' '+in1_list[i]+' '+in2_list[i])

    sh.close()

    #check if slurm shell file should be printed or sent to Odyssey
    if args.print == 'false':
    	#send slurm job to NBI cluster
    	cmd = ('sbatch ca.sh')
    	p = subprocess.Popen(cmd, shell=True)
    	sts = os.waitpid(p.pid, 0)[1]
    else:
        file = open('ca.sh','r')
        data = file.read()
        print(data)

    #send slurm job to NBI cluster
 #   cmd = ('sbatch ca.sh')
#    p = subprocess.Popen(cmd, shell=True)
#    sts = os.waitpid(p.pid, 0)[1]

    count +=1

os.remove('ca.sh')

print('\nSubmitted '+str(count)+' shell scripts to the cluster!\nFinished!!\n')

