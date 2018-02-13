###########################
# Author: Patrick Monnahan
# Purpose:  Repolarize reference/alternative alleles in A. arenosa mapped to north american A. lyrata.  We used genotype calls from mapping European A. lyrata, A. croatica, and A. haleri to the same N. American lyrata reference genome, and identified sites in which these samples showed an average frequency > 0.5.  Since all of these samples are more closely genetically related to A. arenosa, we repolarized the reference/alternative at these sites.
############################


import argparse
import gzip

parser = argparse.ArgumentParser(description = 'This program takes a merged vcf with highest-coverage individuals from specific pops as input, and outputs a lookup key from the merged vcf')
parser.add_argument('-v', type = str, metavar = 'vcf_path', required = True, help = 'path to vcf')
parser.add_argument('-gz', type = str, metavar = 'gzipped?', required = False, default = 'false', help = 'are vcfs gzipped (true) or not (false)')
parser.add_argument('-o', type = str, metavar = 'Output_Prefix', required = True, help = 'Vcfs retain original name but the concatenated lookup key will be a text file specified by output')
parser.add_argument('-mi', type = int, metavar = 'minimum amount individuals', required = True, help = 'minimum amount of individuals required to support alternative alleles')
parser.add_argument('-mp', type = float, metavar = 'min. proportion alt alleles', required = True, help = 'minimum proportion of alternative alleles to allow')
parser.add_argument('-ly', type = str, metavar = 'lyrata_only?', required = False, default = 'false', help = 'do you want to include lyrata only (true) or not (false)?')

args = parser.parse_args()


if args.gz == 'true' and args.v[-3:] == '.gz':
    gzip.gunzip(args.v)
    lookup_table_file = open(args.v+args.o+"repolarized.lookupKey.minAlleles_"+str(args.mi)+".txt", 'w')

lookup_table_file = open(args.o+"repolarized.lookupKey.minInd_"+str(args.mi)+".txt", 'w')

if args.ly == 'true':  
    args.mi = 2 # args.mi must = 2, since there are only two lyrata samples
    args.mp = 1.0

count = 0
type_counts = [0, 0, 0, 0]
count_file = open(args.o + "counts.txt",'w')

with open(args.v) as vcf:
    for line_idx, line in enumerate(vcf): # Cycle over lines in the VCF file
        cols = line.replace('\n', '').split('\t')  # Split each line of vcf
        if len(cols) < 2:               ## This should be info just before header
            pass
        elif cols[0] == "#CHROM": #This should be header
            
            names = [] # list to append pop names to (may not be necessary..)
            for j in cols[9:]: # get names of individuals in vcf
                names.append(j)

        else: 
            scaff = cols[0]               # parse important info from each line     
            position = int(cols[1])
            ref_base = cols[3]
            alt_base = cols[4]
            info = cols[7].split(";")
            AN = float(info[2].split("=")[1])
            AC = float(info[0].split("=")[1])
            newsite=[]
            num_ind = 0
            alt_ind = 0
            het_ind = 0
            min_ind = args.mi
            min_prop_alt = args.mp
            cros = 0
            hals = 0
            lyrs = 0
            cro_alt_alleles = 0
            lyr_alt_alleles = 0
            hal_alt_alleles = 0
            species = 0

            for j, ind in enumerate(cols[9:]):
                ind = ind.split(":")
                if len(ind) ==5:             
                    dp = ind[2]  
                    gt = ind[0]
                    gt = gt.split("/")

                    # below is the 'Lyrata Only' section.
                    if args.ly == 'true':
                        if j == 2 or j == 3:
                            if gt[0] != ".":
                                num_ind += 1
                                if sum([int(x) for x in gt]) == 2:
                                    alt_ind += 1
                                elif sum([int(x) for x in gt]) == 1:
                                    het_ind += 1


                    elif args.ly != 'true':
                        if j == 0 or j == 1:  # Croatica
                            if gt[0] != ".":
                                num_ind += 1
                                cros += 2
                                if sum([int(x) for x in gt]) == 2:
                                    alt_ind += 1
                                    cro_alt_alleles += 2
                                elif sum([int(x) for x in gt]) == 1:
                                    het_ind += 1
                                    cro_alt_alleles += 1
                        if j == 2 or j == 3:  # lyrata
                            if gt[0] != ".":
                                num_ind += 1
                                lyrs += 2
                                if sum([int(x) for x in gt]) == 2:
                                    alt_ind += 1
                                    lyr_alt_alleles += 2
                                elif sum([int(x) for x in gt]) == 1:
                                    het_ind += 1
                                    lyr_alt_alleles += 1
                        if j == 4 or j == 5:  # halleri
                            if gt[0] != ".":
                                num_ind += 1
                                hals += 2
                                if sum([int(x) for x in gt]) == 2:
                                    alt_ind += 1
                                    hal_alt_alleles += 2
                                elif sum([int(x) for x in gt]) == 1:
                                    het_ind += 1
                                    hal_alt_alleles += 1

            ## NOTE: for last 3 elif statements, notice the additional weighting of lyrata and haleri  over croatica, and haleri over lyrata.  This is because haleri > lyrata > croatica in terms of genetic similarity to arenosa 
            if all(k >= 2 for k in [cros, lyrs, hals]):  # 3 species genotyped for at least one indivdiual
                spec_freqs = [float(l) / float(m) for l, m in zip([cro_alt_alleles, lyr_alt_alleles, hal_alt_alleles], [cros, lyrs, hals])]
                if sum([1 for jj in spec_freqs if jj > 0.5]) >= 2: # Asking 
                    lookup_table_file.write(scaff + "\t" + str(position) + "\n")
                    type_counts[0] += 1
            elif all(k >= 2 for k in [cros, lyrs]):  # Genotypes for cro and lyr
                spec_freqs = [float(l) / float(m) for l, m in zip([cro_alt_alleles, lyr_alt_alleles, lyr_alt_alleles], [cros, lyrs, lyrs])] 
                if sum(spec_freqs) / 3.0 > 0.5:
                    lookup_table_file.write(scaff + "\t" + str(position) + "\n")
                    type_counts[1] += 1
            elif all(k >= 2 for k in [cros, hals]):  # Genotypes for cro and hal
                spec_freqs = [float(l) / float(m) for l, m in zip([cro_alt_alleles, hal_alt_alleles, hal_alt_alleles], [cros, hals, hals])]
                if sum(spec_freqs) / 3.0 > 0.5:
                    lookup_table_file.write(scaff + "\t" + str(position) + "\n")
                    type_counts[2] += 1
            elif all(k >= 2 for k in [lyrs, hals]):  # Genotypes for lyr and hal
                spec_freqs = [float(l) / float(m) for l, m in zip([lyr_alt_alleles, hal_alt_alleles, hal_alt_alleles], [lyrs, hals, hals])]
                if sum(spec_freqs) / 3.0 > 0.5:
                    lookup_table_file.write(scaff + "\t" + str(position) + "\n")
                    type_counts[3] += 1
    print(type_counts)

