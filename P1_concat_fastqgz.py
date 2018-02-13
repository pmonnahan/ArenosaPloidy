import os, sys, argparse, subprocess

# create variables that can be entered in the command line
parser = argparse.ArgumentParser()
parser.add_argument('-ind', type=str, metavar='', required=True, help='individual, text file containing the names of the individuals, one individual per row')
parser.add_argument('-s', type=str, metavar='', required=True, help='dir in which to start looking for fastq files, rootDir for os.walk')
parser.add_argument('-o', type=str, metavar='', required=True, help='text nameing the outfolder of concatenated files, starting after Levi-Yant')
parser.add_argument('-lanes', type=int, metavar='', required=True, help='number of lanes the samples specified in -ind were sequenced')
parser.add_argument('-mem', type=str, metavar='', default='5', help='total memory for each job in (Mb)')
parser.add_argument('-print', type=str, metavar='print', default='false', help='If changed to true then shell files are printed to screen and not launched [false]')
args = parser.parse_args()

# create objects to work with
name = []
folders = []
inlist = []
outlist = []
outfolder = [args.o]
baseoutname = []

# populate name list
ind = open(args.ind, 'r')
for x in ind:
    data = x
    name.append(data.replace('\n',''))
ind.close()

# create input list of files by searching the subdirectories of the rootDir (args.s)
for dirName, subdirList,fileList in os.walk(args.s, topdown = False):
    for fname in fileList:
        if fname in fileList:
            folders.append(dirName+ '/' +fname)
            # created a list with absolute paths to all individual files
# create a list of filepaths to all fastq.gz files specified in args.ind
for j in name:
    for i in range(len(folders)):
        if folders[i].endswith(j):
            infile = folders[i]
            inlist.append(infile)

# create output list
for j in name:
    outfile = '/nbi/Research-Groups/JIC/Levi-Yant/' +args.o + j
    outlist.append(outfile)

count = 0

# connect to cluster with an NBI SLURM cluster and
# run samtools merge on created list, args.lanes lines at a time
for x in range(len(name)):
   sh = open('concat'+str(x)+'.sh', 'w')
   baseoutname = outlist[x].replace('.fastqgz', '')
   sh.write('#!/bin/bash\n'+
              '#SBATCH -J CAT.'+name[x]+'\n' +
              '#SBATCH -e '+baseoutname+'.err\n' +
              '#SBATCH -o '+baseoutname+'.out\n' +
              '#SBATCH -p nbi-short\n' +
              #'#SBATCH -t '+args.time+'\n' +
              '#SBATCH --mem='+args.mem+'\n' +
              'cat ')
   for i in range(args.lanes):
       sh.write(inlist[x*args.lanes+i]+ ' ')
   sh.write(' > '+outlist[x]+' ')

   sh.close()

   # send slurm shell to NBI SLURM cluster 
   #cmd = ('sbatch concat'+str(x)+'.sh')
   #p = subprocess.Popen(cmd, shell=True)
   #sts = os.waitpid(p.pid, 0)[1]

   # check if slurm shell file should be printed or sent to NBI SLURM cluster  (from Jeff)
   if args.print == 'false':
        # send slurm shell to NBI SLURM cluster 
        cmd = ('sbatch concat'+str(x)+'.sh')
        p = subprocess.Popen(cmd, shell=True)
        sts = os.waitpid(p.pid, 0)[1]
   else:
       file = open('concat'+str(x)+'.sh','r')
       data = file.read()
       print(data)

   count +=1

os.remove('concat'+str(x)+'.sh')

print('\n'+str(count)+' Jobs submitted to NBI SLURM cluster !')


