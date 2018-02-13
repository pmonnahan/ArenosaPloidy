import argparse
import os
import subprocess

parser = argparse.ArgumentParser()
parser.add_argument('-v', type=str, metavar='vcf_directory', required=True, help='')
parser.add_argument('-dp', type=str, metavar='site_depth_cutoff', required=True, default='2704', help='')
parser.add_argument('-t', type=str, metavar='time', required=False, default='0-06:00', help='')
parser.add_argument('-mem', type=int, metavar='memory', required=False, default='8000', help='')
parser.add_argument('-p', type=bool, metavar='Print?', required=False, default=False, help='')

args = parser.parse_args()

mem1 = int(args.mem / 1000)
vcf_list = []
vcf_basenames = []

if args.v.endswith("/") is False:
    args.v += "/"

if os.path.exists(args.v + "RD/") is False:
    os.mkdir(args.v + "RD/")
outdir = args.v + "RD/"

for file in os.listdir(args.v):
    if file[-6:] == 'vcf.gz':
        vcf_list.append(file)
        vcf_basenames.append(file[:-7])
    elif file[-3:] == 'vcf':
        vcf_list.append(file)
        vcf_basenames.append(file[:-4])
for v, vcf in enumerate(vcf_list):
    # Select single population and biallelic SNPs for each scaffold and convert to variants table
    shfile1 = open("RD." + vcf + '.sh', 'w')
    shfile1.write('#!/bin/bash\n' +
                  '#SBATCH -J RD.' + vcf + '.sh' + '\n' +
                  '#SBATCH -e RD.' + vcf + '.gatk.err' + '\n' +
                  '#SBATCH -o RD.' + vcf + '.gatk.out' + '\n' +
                  '#SBATCH -p nbi-medium\n' +
                  '#SBATCH -n 1\n' +
                  '#SBATCH -t ' + args.t + '\n' +
                  '#SBATCH --mem=' + str(args.mem) + '\n' +
                  'source GATK-nightly.2016.09.26\n' +
                  'java -Xmx' + str(mem1) + 'g -jar /nbi/software/testing/GATK/nightly.2016.09.26/x86_64/jars/GenomeAnalysisTK.jar -T SelectVariants -R /nbi/Research-Groups/JIC/Levi-Yant/Lyrata_ref/alygenomes.fasta -V ' + args.v + vcf + ' -o ' + outdir + vcf_basenames[v] + '.RDF.vcf -select "DP < ' + args.dp + '"\n')

    shfile1.close()

    if args.p is False:  # send slurm job to NBI SLURM cluster
        cmd1 = ('sbatch RD.' + vcf + '.sh')
        p1 = subprocess.Popen(cmd1, shell=True)
        sts1 = os.waitpid(p1.pid, 0)[1]

    else:
        file1 = open("RD." + vcf + '.sh', 'r')
        data1 = file1.read()
        print(data1)

    os.remove("RD." + vcf + '.sh')
