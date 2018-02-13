###########################
# Author: Patrick Monnahan
# Purpose:  Thin sites from a set of VCF's to produce a pruned dataset of unlinked sites for use in clustering algorithms such as Structure.
###########################

import os
import sys
import subprocess
import argparse
import random
import transposer
import numpy
import csv
import scipy
import gzip

# create variables that can be entered in the command line
parser = argparse.ArgumentParser()
parser.add_argument('-v', type=str, metavar='vcf_path', required=True, help='path to vcfs')
parser.add_argument('-w', type=str, metavar='Window_Size', required=True, help='size of scaffold window')
parser.add_argument('-r', type=str, metavar='Number_Replications', required=True, help='Number of replicate data sets')
parser.add_argument('-d', type=str, metavar='Window_Distance', required=True, help='distance between any 2 windows')
parser.add_argument('-m', type=str, metavar='Missing_Data', required=True, help='amount of missing data to allow per site (between 0-1)')
parser.add_argument('-mf', type=str, metavar='Minimum_Frequency', required=True, help='minimum minor allele frequency')
parser.add_argument('-o', type=str, metavar='Output_Prefix', required=True, help='Vcfs retain original scaffold name but the concatenated Structure input file will be a text file with specified by output and within the LD_pruned directory')
parser.add_argument('-s', type=str, metavar='Subset?', required=False, default='false', help='if true, this will subsample polyploid data to create psuedo-diploid data')
parser.add_argument('-c', type=float, metavar='maximum_correlation', required=True, default='1.0', help='Maximum correlation between adjacent sites allowed for site to be used')
parser.add_argument('-gz', type=str, metavar='gzipped?', required=True, help='are vcfs gzipped (true) or not (false)')
# output population as 2nd column (must be integer for STRUCTURE*)
#

args = parser.parse_args()

if args.v.endswith("/") is False:
    args.v += "/"
if os.path.exists(args.v + 'LD_Pruned/') is False:  # Create folder for output if it doesn't already exist
    os.mkdir(args.v + 'LD_Pruned/')
vcf_list = []

for file in os.listdir(args.v):  # get names of vcf files in args.v directory
    if args.gz == 'true':
        if file[-3:] == '.gz':
            vcf_list.append(file)

    elif args.gz == 'false':
        if file[-3:] == 'vcf':
            vcf_list.append(file)

    else:
        print('error')

count = 0

for rep in range(int(args.r)):
    print("Started Rep ", rep)

    markernames = []
    tot_num_screened = 0
    tot_sites = 0

    structtempfile = open(args.v + 'LD_Pruned/' + args.o + ".rep" + str(rep) + ".LD_Pruned.TransposedStruct.txt", 'w')
    subtempfile = open(args.v + 'LD_Pruned/' + args.o + ".rep" + str(rep) + ".LD_Pruned.TransposedStructSubSample.txt", 'w')
    structfile = open(args.v + 'LD_Pruned/' + args.o + ".StructureInput.rep" + str(rep + 1) + ".LD_Pruned.txt", 'w')

    # structfile.write("Ind\tPopName\tPopFlag\t")
    structfile.write("\t\t")

    if args.s == 'true':  # Create files if subset is true.
        subfile = open(args.v + 'LD_Pruned/' + args.o + ".StructureInput.rep" + str(rep + 1) + ".LD_Pruned.Diploidized.txt", 'w')
        # subfile.write("Ind\tPopName\tPopFlag\t")
        subfile.write("\t\t")

    for iii, vcf in enumerate(vcf_list):
        if args.gz == 'true':
            newVCF = open(args.v + 'LD_Pruned/' + vcf[:-6] + "rep" + str(rep) + ".LD_Pruned.vcf", "w")  # new vcf if gzipped previously
            src = gzip.open(args.v + vcf)
        elif args.gz == 'false':
            newVCF = open(args.v + 'LD_Pruned/' + vcf[:-3] + "rep" + str(rep) + ".LD_Pruned.vcf", 'w')  # new vcf for storing info from the random draws
            src = open(args.v + vcf)

        line_num = 0
        first = 10000         # where to start looking in each scaffold
        current_window = []
        current_lines = []
        chosen_sites = []
        chosen_lines = []
        window_num = 0
        site_num = 0
        start = 0
        end = 0
        select = True
        names = []
        ploidy = []
        first_site = True
        r2_list = []
        vcf_sites = 0
        vcf_num_screened = 0
        # evaluate contents of each line of input file
        for line_idx, line in enumerate(src):  # Cycle over lines in the VCF file
            cols = line.replace('\n', '').split('\t')  # Split each line of vcf
            if len(cols) < 2:          # This should be info just before header
                newVCF.write(line)
            elif cols[0] == "#CHROM":  # This should be header

                newVCF.write(line)
                for j in cols[9:]:  # get names of individuals in vcf
                    names.append(j)

            else:
                genos = []
                dgenos = []
                scaff = cols[0]               # parse important info from each line
                position = int(cols[1])
                ref_base = cols[3]
                alt_base = cols[4]      # parsing all these things
                info = cols[7].split(";")
                AN = float(info[2].split("=")[1])
                AC = float(info[0].split("=")[1])
                newsite = []

                # Convert alleles to Structure input
                if ref_base == "A":
                    ref_base = 1
                elif ref_base == "T":
                    ref_base = 2
                elif ref_base == "G":
                    ref_base = 3
                elif ref_base == "C":
                    ref_base = 4
                else:
                    print('stop1')

                if len(alt_base) > 1:    # multiple bases at site (pass b/c don't want to deal with this type of site)
                    alt_base = -99
                elif alt_base == "A":
                    alt_base = 1
                elif alt_base == "T":
                    alt_base = 2
                elif alt_base == "G":
                    alt_base = 3
                elif alt_base == "C":
                    alt_base = 4

                if alt_base == -99:
                    pass

                elif first_site is True and iii == 0:  # Initial setup of info for output files.

                    for i, geno in enumerate(cols[9:]):  # Determine ploidy of each individual
                        geno = geno.split(":")[0]
                        geno = geno.split("/")
                        ploidy.append(len(geno))

                    for j, item in enumerate(names):  # Write individual name information for the temporary file that is to be transposed.  Names need to be repeated for each observed allele.
                        for jj in range(0, ploidy[j]):
                            structtempfile.write("%s\t" % item)
                    structtempfile.write("\n")

                    for j, item in enumerate(names):  # Write Population names for each individual.
                        for jj in range(0, ploidy[j]):
                            structtempfile.write("%s\t" % item[:3])
                    structtempfile.write("\n")

                    for j, item in enumerate(names):  # Write PopFlag for each individual which is a unique integer for each population.
                        if j == 0:
                                oldname = item[:3]
                                popcount = 1
                                for jj in range(0, ploidy[j]):
                                    structtempfile.write("%s\t" % str(popcount))
                        else:
                            if oldname != item[:3]:
                                popcount += 1
                                oldname = item[:3]
                                for jj in range(0, ploidy[j]):
                                    structtempfile.write("%s\t" % str(popcount))
                            else:
                                for jj in range(0, ploidy[j]):
                                    structtempfile.write("%s\t" % str(popcount))
                    structtempfile.write("\n")

                    if args.s == 'true':  # Same task as steps immediately above, but adjusted to accomodate subsetting.
                        for j, item in enumerate(names):
                            for jj in range(0, 2):
                                subtempfile.write("%s\t" % item)
                        subtempfile.write("\n")
                        for j, item in enumerate(names):
                            for jj in range(0, 2):
                                subtempfile.write("%s\t" % item[:3])
                        subtempfile.write("\n")
                        for j, item in enumerate(names):
                            if j == 0:
                                oldname = item[:3]
                                popcount = 1
                                for jj in range(0, 2):
                                    subtempfile.write("%s\t" % str(popcount))
                            else:
                                if oldname != item[:3]:
                                    popcount += 1
                                    oldname = item[:3]
                                    for jj in range(0, 2):
                                        subtempfile.write("%s\t" % str(popcount))
                                else:
                                    for jj in range(0, 2):
                                        subtempfile.write("%s\t" % str(popcount))
                        subtempfile.write("\n")

                    num_alleles = float(sum(ploidy))

                    first_site = False

                    # if first line of vcf passes filters, we go ahead and use it
                    if cols[6] == 'PASS' and 1.0 - (AN / (num_alleles)) < float(args.m) and AC / AN > float(args.mf) and AC / AN < 1.0 - float(args.mf):
                        current_window.append(cols)
                        current_lines.append(line)
                        start = position
                        end = position + int(args.w)
                        window_num += 1
                        line_num += 1
                        site_num = 1

                # if first line of vcf does not pass filter, then this statement will catch the first line that does
                elif position > first and site_num == 0:
                    if cols[6] == 'PASS' and 1.0 - (AN / (num_alleles)) < float(args.m) and AC / AN > float(args.mf) and AC / AN < 1.0 - float(args.mf):
                        current_window.append(cols)
                        current_lines.append(line)
                        start = position
                        end = position + int(args.w)
                        window_num += 1
                        line_num += 1
                        site_num = 1

                # All lines caught by this statement are within the current window
                elif position > start and position < end and site_num != 0:
                    if cols[6] == 'PASS' and 1.0 - (AN / (num_alleles)) < float(args.m) and AC / AN > float(args.mf) and AC / AN < 1.0 - float(args.mf):
                        current_window.append(cols)
                        current_lines.append(line)
                        site_num += 1

                # Here, we have moved past the current window and now need to select a site from all sites within the current window before resetting window bounds and moving on to next window.
                elif position > end and select is True and site_num != 0 and AC != 0 and AC != AN:
                    # select sites here
                    if len(current_window) != 0:

                        rx = random.randrange(0, len(current_window))
                        site = current_window[rx]
                        genox = []
                        for geno in site[9:]:
                            geno = geno.split(":")[0]
                            geno = geno.split("/")
                            ac = 0

                            for allele in geno:
                                if allele == '.':
                                    genos.append(-9)
                                    ac = -9
                                elif allele == '0':
                                    genos.append(ref_base)
                                elif allele == '1':
                                    ac += 1
                                    genos.append(alt_base)
                                else:
                                    print("allele not matched")

                            genox.append(ac)

                            for allele in numpy.random.choice(geno, 2, replace=False):
                                if allele == '.':
                                    dgenos.append(-9)
                                elif allele == '0':
                                    dgenos.append(ref_base)
                                elif allele == '1':
                                    dgenos.append(alt_base)
                                else:
                                    print("allele not matched")


                        if window_num == 1:
                            oldsite = list(genox)

                        else:
                            newsite = list(genox)

                            for kk, oo in enumerate(oldsite):
                                if oo == -9:

                                    oldsite.pop(kk)
                                    newsite.pop(kk)

                            for yy, nn in enumerate(newsite):
                                if nn == -9:
                                    oldsite.pop(yy)
                                    newsite.pop(yy)

                            # Do correlation calculations here
                            r2 = float(numpy.corrcoef(oldsite, newsite)[1][0])**2.0
                            r2_list.append(r2)

                            oldsite = list(genox)

                            if r2 < float(args.c):
                                vcf_sites += 1
                                tot_sites += 1
                                line1 = current_lines[rx]
                                chosen_sites.append(site)  # filling chosen_sites array with random choices.

                                markernames.append(str(site[0]) + "_" + str(site[1]))

                                for item in genos:
                                    structtempfile.write("%s\t" % item)
                                structtempfile.write("\n")
                                if args.s == 'true':
                                    for item in dgenos:
                                        subtempfile.write("%s\t" % item)
                                    subtempfile.write("\n")

                                newVCF.write(line1)
                            else:
                                tot_num_screened += 1
                                vcf_num_screened += 1

                        select = False
                    else:
                        start = position
                        end = position + int(args.w)
                        current_window = []
                        current_lines = []
                        site_num = 1

                # Here, a new site is found that is greater than user-specified distance between windows, so values are reset and new window is established.
                elif position > end + int(args.d) and site_num != 0:
                    if cols[6] == 'PASS' and 1.0 - (AN / (num_alleles)) < float(args.m) and AC / AN > float(args.mf) and AC / AN < 1.0 - float(args.mf):
                        select = True
                        site_num = 1
                        window_num += 1
                        start = position
                        end = position + int(args.w)
                        current_window = []
                        current_lines = []

                # Prints the line of vcf as we scroll through the vcf file
                if line_idx % 100000 == 0:
                    print(line_idx)

        print 'Finished vcf: ', vcf
        print 'mean r2: ', numpy.mean(r2_list)
        print 'Screened sites for vcf: ', vcf_num_screened
        print 'Tot sites for vcf: ', vcf_sites

    # This section transposes the temporary files created during this replicate and adds headers, so that they are formatted according to structure
    structtempfile.close()
    subtempfile.close()

    # Transposes file
    jj = transposer.transpose(i=args.v + 'LD_Pruned/' + args.o + ".rep" + str(rep) + ".LD_Pruned.TransposedStruct.txt", d="\t",)
    kk = transposer.transpose(i=args.v + 'LD_Pruned/' + args.o + ".rep" + str(rep) + ".LD_Pruned.TransposedStructSubSample.txt", d="\t",)

    # Write header names for each marker
    for marker in markernames:
        structfile.write("%s\t" % marker)
        if args.s == 'true':
            subfile.write("%s\t" % marker)
    structfile.write("\n")
    if args.s == 'true':
        subfile.write("\n")
    for item in jj:
        structfile.write("%s" % item)
        structfile.write("\n")
    if args.s == 'true':
        for item in kk:
            subfile.write("%s" % item)
            subfile.write("\n")

    # remove the temporary files that contained the info that needed to be transposed
    os.remove(args.v + 'LD_Pruned/' + args.o + ".rep" + str(rep) + ".LD_Pruned.TransposedStruct.txt")
    os.remove(args.v + 'LD_Pruned/' + args.o + ".rep" + str(rep) + ".LD_Pruned.TransposedStructSubSample.txt")


    print 'Finished Rep Number: ', rep
    print 'Number of screened, linked sites: ', tot_num_screened
    print 'Total number of non-screened sites: ', tot_sites
